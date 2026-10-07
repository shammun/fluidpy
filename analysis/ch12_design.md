# Chapter 12 — Turbulence: lesson design
(from `analysis/ch12_curation.md` (A 16 · B 210 · C 17 = 243 rows; CORE 16 · NOTE 217 · RECAP 8 · SKIP 2; 28 derivations D01–D28
(★ 10 · ★★ 14 · ★★★ 4: D07, D08, D10, D14); E1–E10 + backup B1; primers from **P280**) and `analysis/ch12.md` (§2 inventory with the
LaTeX the analyst read from the page images, §2b derivations a-D1…a-D46, §4 implementation rows 1–46, §5 the two new core modules,
§9 conventions C1–C5 and slips #1–#15). 2026-10-06, lesson-designer. **Audited and completed on 2026-10-07** after the first session was interrupted: see C.6, A.15 and the closing check after Part F.

**Equation provenance (honest statement).** Every equation placed here is taken from the analyst's page-image transcription in
`analysis/ch12.md` §2. For this design I re-read eight rendered pages against that transcription and found no difference:
p583 (term-by-term averages, the collected equation, $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-g[1-\alpha(\bar T-T_0)]\delta_{i3}+\frac1{\rho_0}\frac{\partial\bar\tau_{ij}}{\partial x_j}$ (12.30)); p589 (the four displays from $f(r)\equiv\overline{u_\parallel(\mathbf x+\mathbf r)u_\parallel(\mathbf x)}/\overline{u_\parallel^2}$ (12.38) to $R_{ij}=\overline{u^2}\big\{f(r)\delta_{ij}+\frac r2\frac{df}{dr}\big(\delta_{ij}-\frac{r_ir_j}{r^2}\big)\big\}$ (12.41), Λ_g = Λ_f/2, λ_g = λ_f/√2); p592 (the 1/Re ratio, $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47) with its labels — slip #2 seen: the label under the left side names Ē); p596 ($\frac{\lambda_T}L\propto\mathrm{Re}_L^{-1/2}$ (12.52), $\frac{S_{11}(k_1)}{u_K^2\eta}=\Phi(k_1\eta)$ (12.53), and (12.54), printed $S_{11}(k_1)=\mathit{const}\cdot\bar\varepsilon^{2/3}\cdot k_1^{5/3}$ — ⚠️ slip #1 seen: the exponent must be −5/3, $S_{11}=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$); p603 (V eliminated, the substituted equation, $\big\{\frac{\delta U'_{CL}}{U_{CL}}\big\}F^2-\big\{\frac{\delta U'_{CL}}{U_{CL}}+\delta'\big\}F'\int_0^\xi F\,d\xi=\big\{\frac\Psi{U_{CL}^2}\big\}G'$ (12.63), $\frac{\delta U'_{CL}}{U_{CL}}=C_1$, $\frac{\delta U'_{CL}}{U_{CL}}+\delta'=C_2$, $\frac\Psi{U_{CL}^2}=C_3$ (12.64), $J_s=\rho U_{CL}^2\delta\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi=\rho C_4^2x^{2\gamma+1}\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi$ (12.65)); p613 ($-\xi\frac{dF}{d\xi}=y^+\frac{df}{dy^+}$ (12.87), $U^+=\frac1\kappa\ln(y^+)+B$ (12.88), $F(\xi)=-\frac1\kappa\ln(\xi)+A$ (12.89));
p625 ($L_M\equiv-u_*^3/(\kappa\alpha g\overline{wT'})$ (12.110), $\mathrm{Rf}=z/L_M$ (12.111), the log-linear profile); p630 ($\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{X_\alpha u_\alpha}=2\int_0^t\overline{u_\alpha(t')u_\alpha(t)}\,dt'$ (12.116), r_α, $\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{u_\alpha^2}\int_0^tr_\alpha(\tau)\,d\tau$ (12.117), $\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\int_0^tdt'\int_0^{t'}r_\alpha(\tau)\,d\tau$ (12.118), the integration by parts, $\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\,t\int_0^t\big(1-\frac\tau t\big)r_\alpha(\tau)\,d\tau$ (12.119)).
**Inlining pass of 2026-10-07 (a second honest statement).** While every number-only mention of an equation in this file was completed with the equation itself, a further 34 rendered pages of ch12 were read against the transcription — p573, p574, p577–p582, p585–p588, p590, p591, p594, p595, p601, p602, p604, p605, p607, p610–p612, p617–p620, p622, p624, p627, p629, p631, p634 — again with no difference (so 42 of the chapter's 80 pages are now checked; the others are not, and the rule of the next sentence still holds for them), and nine pages of earlier chapters for the cross-references: ch01 p045 and p047, ch04 p142 and p163, ch07 p321, ch08 p343, ch09 p392, ch11 p532 and p546. One finding there: the book prints the pipe wall stress as $\tau_0=\frac a2\frac{dp}{dz}$ (8.8), with the sign of the pressure gradient, where recap R04 had written $-(a/2)\,dp/dx$; R04 now shows the printed form and says why Chapter 12's positive $\tau_0$ needs the minus sign.
The remaining pages were **not** re-read by the designer; builders render the page (`tools/render_pages.py ch12 --eq N.M`) before
setting any equation not in that list. Numbers marked *expect* were computed for this design with a scratch script (numpy/scipy
closed forms and quadrature; no fluidpy function existed yet); the builder compares an executed cell against them and reports
any difference instead of editing the expectation.

**Part C is written first and is the contract both the implementer and the builders keep.**

**Binding conventions for every builder (analysis §9, curation decisions 5–9, §8, §9).**
1. **Imports and aliases.** `from fluidpy import ch12_turbulence as ch12`; `from fluidpy.core import turbstats as TS,
   wall_turbulence as WT`. **`ch12` re-exports every public name of `TS` and `WT`** (the ch11 pattern), so explainer parity rows
   write `ch12.<name>` only. Reused: `from fluidpy.core import dimensional as DIM, stratification as STRAT, laminar as LAM,
   boundary_layer as BL, jets as JET, diffusion as DIFF`; `from fluidpy import ch11_instability as ch11`.
2. **Scalar-callable and parity-friendly** (the `shot.py` evaluator knows only `np`, `math`, `ch12`): every public function
   accepts Python floats and returns a float, a tuple or a `dict` of floats (arrays only when an array goes in). Parity `py:` rows
   use dict keys, integer indices and `np.…` only (never ndarray methods or builtins).
3. **Units.** SI. g = G0 = 9.80665 m/s² (`core.thermo.G0`) unless a cell says "g ≈ 10" for hand arithmetic. Heat flux `H` [W/m²]
   and `wT` [K m/s] positive upward. Spectra: two-sided in angular frequency [rad/s] or wavenumber [rad/m], 1/2π in the forward
   transform, ∫S = variance; `one_sided=True` doubles the density.
4. **Symbols (overloading resolved — the conventions block before C01 shows this table; each block repeats the row it needs).**

| Book symbol | Meanings in this chapter (and earlier) | Notebook / explainer symbol | Code name |
|---|---|---|---|
| over-bar, tilde, capital, lower case | $\overline{(\ )}$ ensemble average; $\tilde u$ total field; $U$, $\bar T$ mean; $u$, $T'$ **fluctuation** (so $u$ is not "the velocity" here) | same | `samples` (members on axis 0), `mean`, `fluct` |
| κ | von Kármán constant in $U^+=\frac1\kappa\ln(y^+)+B$ (12.88) and $L_M\equiv-u_*^3/(\kappa\alpha g\overline{wT'})$ (12.110); thermal diffusivity in $\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}+\frac{\partial}{\partial x_j}(\overline{u_jT'})=\kappa\frac{\partial^2\bar T}{\partial x_j^2}$ (12.31) and in the dissipation $\bar\varepsilon_T=\kappa\overline{(\partial T'/\partial x_j)^2}$ of the temperature-variance budget (N190 in C15); κ_m, κ_T, κ_mT nearby | $\kappa$ = von Kármán; $\kappa_{th}$ = molecular thermal diffusivity in our own lines (the book's κ quoted inside its equations with a ⚠️ line); $\kappa_T$ eddy diffusivity | `kappa`; `kappa_th`; `kappa_T` |
| k, K | thermal conductivity in $Q_j=-k\frac{\partial\bar T}{\partial x_j}+\rho_0C_p\overline{u_jT'}$ (12.32); wavenumber $k_1$; 3-D wavenumber $K$; the "k" of k–ε (the book's $\bar e$); kurtosis K (Ex. 12.1) | $k_{th}$ conductivity; $k_1$, $K$ wavenumbers; $\bar e$ turbulent kinetic energy ("k" only in the name k–ε) | `k_th`; `k1`, `K`; `e` |
| e, E, ε | $\bar e=\tfrac12\overline{u_i^2}$ (ch01: internal energy); $\bar E=\tfrac12U_i^2$ (ch11: disturbance energy); $\bar\varepsilon$ dissipation, $\bar\varepsilon_T$ its thermal twin | $\bar e$, $\bar E$, $\bar\varepsilon$, $\bar\varepsilon_T$ | `e`, `E_mean`, `eps`, `eps_T` |
| λ, Λ | Taylor microscales $\lambda_t,\lambda_f,\lambda_g,\lambda_T$ (earlier: wavelength, eigenvalue); integral scales $\Lambda_t,\Lambda_f,\Lambda_g$ (ch11: Λ = dissipation) | same, always with the subscript | `lambda_t`, `lambda_f`, `lambda_g`; `Lambda_t`, `Lambda_f`, `Lambda_g` |
| η | Kolmogorov length (earlier: similarity variable, surface elevation); $\eta_T$ Batchelor scale | $\eta$, $\eta_T$ | `eta`, `eta_T` |
| f, g, F, G | correlation coefficients $f(r)$, $g(r)$: $f(r)\equiv\overline{u_\parallel(\mathbf x+\mathbf r)u_\parallel(\mathbf x)}/\overline{u_\parallel^2}$ (12.38), g likewise with $u_\perp$; law of the wall $U^+=f(y^+)$ (12.80); Darcy $f$; $g$ gravity; tensor functions in $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40); jet profiles in $U=U_{CL}(x)F(y/\delta(x))$ (12.56) and $-\overline{uv}=\Psi(x)G(y/\delta(x))$ (12.57); defect function in $\frac{U_\infty-U}{u_*}=F(\xi)$ (12.84) | $f(r)$, $g(r)$; $f_w(y^+)$ in our own lines for the wall function; $f_D$ Darcy; $F_R$, $G_R$ in $R_{ij}=F_Rr_ir_j+G_R\delta_{ij}$ (12.40) in our own lines; $F$, $G$ jet; $F_d$ defect in our own lines | `f`, `g_corr`; `fD`; `F_R`, `G_R`; `F`, `G` |
| $R_{ij}$, r | correlation tensor (ch02–ch03: rotation tensor); r separation, radius; $r_{11}$, $r_\alpha$ correlation coefficients | $R_{ij}$ correlation; $r_{11}(\tau)$, $r_\alpha(\tau)$ | `R`, `r` |
| τ | time lag in $R_{11}(\tau)=\overline{u_1(t)u_1(t+\tau)}$ (12.17); decay time (Ex. 12.1); stress $\bar\tau$, $\tau_0$ | $\tau$ lag; $\tau_d$ decay time in our Example 12.1 lines; $\bar\tau$, $\tau_0$ stress; $\tau_c$ memory time of our test signals | `lag`, `tau`, `tau0`, `tau_c` |
| δ, θ, α | δ jet width, boundary-layer thickness, $\delta_{ij}$; θ momentum thickness (§12.9), potential temperature (ours, §12.11); α thermal expansion and the no-sum index of §12.12 | δ, $\delta_{ij}$; $\theta_m$ momentum thickness where both meet; θ potential temperature; α expansion; "component α (no sum)" said in words | `delta`; `theta_m`; `theta`; `alpha` |
| S, N, L, σ, Π, U₀ | spectra $S_e,S_{11},S_T,S(K)$ vs strain $\bar S_{ij},S'_{ij}$; N realizations vs buoyancy frequency; L outer scale, walk step, $L_M$; σ Gaussian width ($\sigma^2=2\nu t$ here, ch03's core radius had $4\nu t$), $\sigma_e,\sigma_\varepsilon$; Π wake strength (ch01: Π groups); $U_0$ probe speed (§12.4), nozzle speed (§12.8) | as in the book, each named at first use; $N_{bv}$ in our own lines when N members is on the same page | `n_members`, `N2`; `L`, `step`, `L_M`; `Pi` |
| h, d | h **full** channel height in $dP/dx=-2\tau_0/h$ (12.90); d slot width, nozzle or pipe diameter | h (full), δ = h/2 | `h`, `delta`, `d` |

5. **Normalisations** (analysis §9 C4): spectra as in 3; $\overline{u^2}$ in §12.6 is **one** component ($\bar e=\tfrac32\overline{u^2}$);
   the book's skewness and kurtosis are un-normalised central moments (`statistics(normalized=True)` gives the usual ones);
   Reynolds numbers always named ($\mathrm{Re}_L=\Delta UL/\nu$, $R_\lambda$ with which λ, $\mathrm{Re}_x$, $\mathrm{Re}_\tau=\delta^+$);
   y from the wall; jet variable ξ = y/x (δ = x by convention, virtual origin dropped).
6. **Assumptions that switch silently** (C5): Boussinesq with buoyancy in §12.5–12.7 → constant density in §12.8–12.10
   (ρ₀ → ρ, P = deviation from hydrostatic) → buoyancy back in §12.11; theory in ensemble averages, measurements in time averages
   (ergodicity, P282). A `> ⚠️` line marks each switch.
7. **Lapse rate (project rule).** Compute with Kundu's $\Gamma\equiv dT/dz$ ($\Gamma_a\approx-9.8$ K/km; stable when $dT/dz>\Gamma_a$) and
   always show the meteorological $\Gamma_{met}\equiv-dT/dz$ ($\Gamma_d\approx+9.8$ K/km; stable when $\Gamma_{met}<\Gamma_d$) beside it: a
   two-row table, both inequalities in every legend, badge and slider trace. $\bar T$, $T'$ in every buoyancy term are **potential**
   temperature; with the thermometer temperature $N^2=g\alpha(dT/dz-\Gamma_a)=g\alpha(\Gamma_d-\Gamma_{met})$.
8. **Book slips taught in corrected form**: `> ⚠️ **slip #k — the book prints** … **; the correct form is** …` where used; keys
   `ch12.book_slips()`.

| Slip | The book prints | Correct | Taught in | Planted wrong variant a test must fail |
|---|---|---|---|---|
| #1 | (12.54) $S_{11}=\mathit{const}\cdot\bar\varepsilon^{2/3}k_1^{5/3}$ | $S_{11}=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ | C08, D12, E4 | `inertial_spectrum_1d(printed=True)` |
| #2 | label under the left side of $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47): "change of Ē" | the turbulent $\bar e$ | C06, D10 | — |
| #3 | (12.70), printed $\dot M_s=\rho_s\int_{-\infty}^{+\infty}[U]_{y=0}\,dy$: the source integral taken at y = 0 | $\dot M_s=\rho_s\int_{-\infty}^{+\infty}[U]_{x=0}\,dy$, as in $J_s\equiv\rho_s\int_{-\infty}^{+\infty}[U^2]_{x=0}dy$ (12.62) | C09 (N118) | `scalar_flux_per_span` at y = 0 |
| #4 | after $-\overline{uv}=C_3U_{CL}^2G(y/x)=C_3C_5^2(J_s/\rho)x^{-1}G(y/x)$ (12.67): constants "C₃ and C₄"; "C₅, C₆ tabulated" | C₃ and C₅ ($C_5=C_4(\rho/J_s)^{1/2}$); the tabulated pair multiplies $U=C_5U_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}F(y/x)$ (12.72) and $\bar Y=C_7Y_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}H(y/x)$ (12.73) | C09 (N115), D15 | — |
| #5 | after $\frac{\delta U'_{CL}}{U_{CL}}=C_8\big(\frac{\delta U'_{CL}}{U_{CL}}+\delta'\big)=C_9\frac\Psi{U_{CL}^2}$ (12.74): the exponential family δ ~ e^{ax}, U_CL ~ e^{−ax} | middle coefficient ≡ 0 and $U_{CL}^2\delta$ not constant; only δ ~ x^m, U_CL ~ x^n with m + 2n = 0 keeps $J_s=\rho\int_{-\infty}^{+\infty}U^2dy=\mathit{const.}$ (12.62) | C09 (N125), E6 | `general_similarity_check` on the exponential family |
| #6 | (12.97), printed $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-\frac1\rho\frac{\partial P}{\partial x_j}+\frac{\partial}{\partial x_j}\big([\nu+\nu_T]\big(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\big)-\frac23\bar e\delta_{ij}\big)$: pressure gradient $\partial P/\partial x_j$ | $\partial P/\partial x_i$ (free index) | C12 (N165) | `rans_eddy_viscosity_residual(printed=True)` |
| #7 | p. 590 "Exercise 12.31" | Exercise 12.32 | C10 (R03), D17 | — |
| #8 | p. 595 $2\nu\overline{u_jS'_{ij}}$ | $2\nu\overline{u_iS'_{ij}}$, as in the exact budget written out in row #2 | C13 (N175) | — |
| #9 | p. 595 "five" constants, six entries | five: $C_\mu,C_{\varepsilon1},C_{\varepsilon2},\sigma_e,\sigma_\varepsilon$ | C13 (N178) | — |
| #10 | (12.75), printed $0=-U\frac{\partial\bar e}{\partial x}-V\frac{\partial\bar e}{\partial y}-\overline{uv}\frac{\partial U}{\partial y}-\frac{\partial}{\partial y}\big(\frac1{\rho_0}\overline{pv}+\frac12\overline{ev}\big)-\bar\varepsilon$, and (12.106), printed $\frac{\partial\bar e}{\partial t}+U\frac{\partial\bar e}{\partial x}=-\frac{\partial}{\partial z}\big(\frac1{\rho_0}\overline{pw}+\overline{ew}\big)-\overline{uw}\frac{\partial U}{\partial z}+g\alpha\overline{wT'}-\bar\varepsilon$: the triple correlation appears as $\tfrac12\overline{ev}$ and as $\overline{ew}$ | $\tfrac12\overline{u_i^2u_j}$ throughout | C06 (D10), C09 (N128), C14 (N182) | — |
| #11 | p. 604 "(11.119)" | $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\,t\int_0^t(1-\tau/t)r_\alpha\,d\tau$ (12.119) | C16 (N207), D27 | — |
| #12 | (12.129), printed $D_T\cong\overline{u_\alpha^2}\Lambda_t$ for $t\ll\Lambda_t$ | $D_T\cong\overline{u_\alpha^2}\Lambda_t$ for $t\gg\Lambda_t$ | C16 (N215), D28, E10 | `dispersion_regime` with swapped condition |
| #13 | Fig. 12.27 caption: width ∝ √x near, ∝ x far | ∝ x near the source, ∝ √x far | C16 (N212), E10 | `smoke_plume_width` slopes |
| #14 | Exercise 12.18a cites the scale definitions $\Lambda_f\equiv\int_0^\infty f(r)\,dr$, $\Lambda_g\equiv\int_0^\infty g(r)\,dr$, $\lambda_f^2\equiv-2/[d^2f/dr^2]_{r=0}$, $\lambda_g^2\equiv-2/[d^2g/dr^2]_{r=0}$ (12.39) | $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40) | C05 (N69), D07 | — |
| #15 | (12.112), printed $\frac{\partial}{\partial t}\big(\frac12\overline{T'^2}\big)+U\frac{\partial}{\partial x}\big(\frac12\overline{T'^2}\big)=-\overline{wT'}\frac{d\bar T}{dz}-\frac{\partial}{\partial z}\big(\frac12\overline{T'^2w}-\kappa\frac{\partial\overline{T'^2}}{\partial z}\big)-\bar\varepsilon_T$: molecular transport $\kappa\,\partial\overline{T'^2}/\partial z$ | $\kappa\,\partial(\tfrac12\overline{T'^2})/\partial z$ | C15 (N190) | `temperature_variance_sympy` |

9. **Colour code (curation §5; the same in figures, derivations, bars and explainers):** mean = purple (accent) · fluctuation and
   turbulent energy = teal · Reynolds stress / shear production = orange · viscous stress / dissipation = rose · buoyancy and
   temperature = blue · scalar = amber · laminar or reference ghosts = muted.
10. **Public-repo rule.** Table 12.1's constants and half-widths, the per-flow (κ, B) pairs, the wake strength, worked-example
    inputs and answers, the stirred-vessel and convection numbers, the fully-turbulent thresholds and all exercise data are
    **never** printed; examples use our own inputs and **labelled illustrative constants** (ξ½ = 0.10, κ = 0.41, B = 5.0 "a common
    textbook pair", Π passed explicitly) or cited public values (`WT.LOG_LAW_CONSTANTS`, `ch12.K_EPSILON_CONSTANTS`, Kolmogorov
    C ≈ 1.5 from Sreenivasan 1995). κ, B, Π and the free-shear constants are required keywords.
11. **Statistics are estimates.** Every sampled cell fixes its seed; assertions on sampled data use 5 standard errors (P281);
    synthetic fields are captioned "kinematic — prescribed spectrum, no cascade".

---

## Part C — functions the builders will call (the implementer's contract)

Names and signatures are those of `analysis/ch12.md` §4 (rows 1–46) and curation §8 (items 1–11). "Eq." = the book equation
implemented (written out in the docstring). Items marked **(+)** are not in analysis §4 (C.4 lists them). Return shapes are
part of the contract (convention 2).

### C.0 Reused (existing; called, not changed)
`DIM.pi_groups`, `DIM.solve_exponents` (ch01) · `STRAT.lapse_rate_stability`, `STRAT.brunt_vaisala_sq_from_lapse`, `STRAT.potential_temperature_gradient` (ch01; `lapse_rate_stability(dT_dz, Gamma_a=None, convention="kundu" or "meteorology")` returns `.verdict` and `.text`) · `ch11.gradient_richardson` (11.66) · `LAM.channel_flow`, `LAM.pipe_wall_stress`,
`LAM.pipe_friction_factor` (ch08) · `BL.blasius_profile`, `BL.blasius_skin_friction`, `BL.displacement_thickness`,
`BL.momentum_thickness` (ch09) · `JET.free_jet`, `JET.jet_momentum_flux` (ch09) · `DIFF.gaussian_spreading` (ch01) ·
`fluidpy.core.fd.stretched_grid` (ch10) · machinery: `setup_notebook`, `show_viz`, `animate`, `show_animation`, `slider_figure`,
`animate_figure`, `live`, `savefig`, `COLORS`.

### C.1 `fluidpy/core/turbstats.py` (NEW, `TS`; re-exported by `ch12`) — statistics, correlations, spectra

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `make_ensemble(n_members, t, mean_fn, sigma, tau_c, seed=0)` | array (N, nt): mean_fn(t) + Ornstein–Uhlenbeck fluctuation of std `sigma`, memory `tau_c` [s] (exact update) | N07; test bed for (12.1) |
| `ensemble_average(samples, axis=0)` | mean over members | (12.10) $\overline{u}=\frac1N\sum_nu(\mathbf x,t{:}n)$ |
| `moment(samples, m, axis=0)` | m-th moment | (12.1) |
| `central_moment(samples, m, axis=0)` | m-th central moment | (12.11) |
| `statistics(samples, normalized=False)` | dict mean, variance, skewness, kurtosis, std, rms | N22 (book: raw central moments) |
| `standard_error(samples)` | std/√N of the member mean | N20 (ours) |
| `standard_error_of_mean(sigma, N)` **(+)** | σ/√N (float) | N20 (parity form) |
| `time_average(t, u, window, m=1)` | centred sliding average (NaN within window/2 of the ends) | (12.2) |
| `volume_average(field, weights=None, m=1)` | float | (12.3) |
| `check_averaging_rules(samples_u, samples_v, t, A=2.0)` | dict of residuals: "sum", "constant", "d/dt", "integral dt", "d/dx", "integral dx", "average of average", "product" (the only non-zero one) | (12.4)–(12.9) |
| `product_average_split(samples_u, samples_v)` | (mean of product, product of means, covariance) | N18 $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$ |
| `reynolds_decompose(samples, axis=0)` | (mean, fluctuations) | (12.24)–(12.26) |
| `correlation(a, b, axis=0)` | $\overline{ab}$ of the fluctuations | (12.12), (12.13) |
| `correlation_coefficient(a, b, axis=0)` | r in [−1, 1] | (12.14), (12.15), (12.16) |
| `correlated_pair(n, r, sigma_u=1.0, sigma_v=1.0, seed=0)` | (u, v) Gaussian samples with correlation r | N25, N63 |
| `autocorrelation(u, dt, max_lag=None, method="fft", unbiased=True, demean=True)` | (lag [s], R [u²]) lag ≥ 0 | (12.17) |
| `cross_correlation(u, v, dt, max_lag=None)` | (lag of both signs, R) | (12.17), N30 |
| `integral_scale(lag, r, upto="first_zero")` | Λ [unit of lag] | (12.18), (12.39) |
| `correlation_time(lag, r)` | first zero t_c | N32 |
| `effective_samples(record_length, t_c)` | N ≈ Δt/t_c | N32 |
| `taylor_microscale(lag, r, fit_points=5)` | λ [unit of lag] | (12.19), (12.39) |
| `integral_scale_from_spectrum(S0, variance)` | Λ = πS(0)/variance | N37 |
| `spectrum_from_correlation(lag, R, omega)` | two-sided S [u² s/rad] | (12.20), (12.45) |
| `correlation_from_spectrum(omega, S, lag)` | R | (12.21) |
| `spectrum_variance(omega, S, two_sided=True)` | variance | (12.22), (12.55) |
| `periodogram(u, d, two_sided=True, segments=1, window="boxcar")` | (omega or k ≥ 0, S) with ∫S = variance | N38 (Ex. 12.8) |
| `taylor_frozen(t, u, U0)` | (x = U₀t, u) | N40 |
| `frequency_to_wavenumber_spectrum(omega, S, U0)` | ($k_1=\omega/U_0$, $S_{11}=U_0S_e$) | N40 |
| `synthetic_solenoidal_field(n, L, spectrum, seed=0, dim=2)` | (u, v[, w]) periodic, divergence-free, kinematic | N02 |
| `white_noise_field(n, seed=0)` | (u, v) | N02 |
| `spatial_correlation(a, b, dx, axis)` | (r, R) | (12.23) |
| `longitudinal_transverse_correlation(u, v, dx)` | (r, f, g) | (12.38) |
| `velocity_covariance(u_samples)` **(+)** | $\overline{u_iu_j}$ [m²/s²] | N49 |
| `reynolds_stress(u_samples, rho0=1.0)` | tensor $-\rho_0\overline{u_iu_j}$ [Pa] | N49 |
| `shell_spectrum(components, L, n_bins=None)` **(+)** | (K, E(K)) of a periodic field | N64 (check of the synthetic field) |
| `anisotropy_tensor(uu)` | $b_{ij}$ | N49 (ours) |
| `turbulent_kinetic_energy(u_samples)` | $\bar e$ [m²/s²] | N72 |
| `correlation_spectrum_pair(kind, sigma, tau_c, omega0=0.0)` **(+)** | dict: callables "r", "S"; floats "Lambda_t", "lambda_t" (NaN for "exponential"), "S0", "t_c" (inf if no zero); kind ∈ exponential, gaussian, damped_cosine | (12.17)–(12.22) closed forms |
| `smooth_signal(n, dt, spectrum, seed=0)` **(+)** | random-phase signal with a prescribed smooth spectrum | for (12.19) |

Closed forms `correlation_spectrum_pair` must return (σ² = variance): exponential $r=e^{-\lvert\tau\rvert/\tau_c}$, $S_e=\sigma^2\tau_c/[\pi(1+\omega^2\tau_c^2)]$,
$\Lambda_t=\tau_c$; gaussian $r=e^{-\tau^2/\tau_c^2}$, $S_e=\sigma^2\tau_c\,e^{-\omega^2\tau_c^2/4}/(2\sqrt\pi)$, $\Lambda_t=\tfrac{\sqrt\pi}2\tau_c$, $\lambda_t=\tau_c$;
damped cosine $r=e^{-\lvert\tau\rvert/\tau_c}\cos\omega_0\tau$, $S_e=\frac{\sigma^2\tau_c}{2\pi}\big[\frac1{1+(\omega-\omega_0)^2\tau_c^2}+\frac1{1+(\omega+\omega_0)^2\tau_c^2}\big]$,
$\Lambda_t=\tau_c/(1+\omega_0^2\tau_c^2)$, $t_c=\pi/(2\omega_0)$.

### C.2 `fluidpy/core/wall_turbulence.py` (NEW, `WT`; re-exported by `ch12`) — wall units, log law, composite profiles

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `friction_velocity(tau0, rho)` | $u_*$ [m/s] | (12.81) $u_*^2=\tau_0/\rho$ |
| `viscous_length(nu, u_star)` | $l_\nu=\nu/u_*$ [m] | (12.80) |
| `friction_reynolds_number(u_star, delta, nu)` | $\delta^+$ | N141 |
| `wall_units(y, U, tau0, rho, nu)` | (y⁺, U⁺, u_*, l_ν) | (12.80) |
| `from_wall_units(yplus, Uplus, u_star, nu)` | (y, U) | (12.80) inverse |
| `law_of_the_wall_groups()` | Π groups of (U, ρ, τ₀, ν, y) | (12.79) |
| `defect_law_groups()` | Π groups of (U_∞ − U, ρ, τ₀, δ, y) | (12.83) |
| `viscous_sublayer(yplus)` | U⁺ = y⁺ | (12.82) |
| `log_law(yplus, *, kappa, B)` | U⁺ | (12.88) |
| `log_law_defect(xi, *, kappa, A)` | $F_d(\xi)$ | (12.89) |
| `fit_log_law(yplus, Uplus, window=(30.0, 0.15), Re_tau=None)` | dict kappa, B, n_points, yplus_min, yplus_max, rms_residual | (12.88) fit |
| `log_law_indicator(yplus, Uplus)` | $y^+dU^+/dy^+$ | (12.87) |
| `layer_name(yplus, y_over_delta)` | "viscous sublayer" (y⁺ < 5) / "buffer layer" (5 ≤ y⁺ < 30) / "logarithmic layer" / "wake region" (y/δ > 0.15) | N139, N149 |
| `friction_law_from_overlap(Re_tau, *, kappa, A, B)` | $U_\infty^+=\frac1\kappa\ln\delta^++A+B$ | D19 (ours) |
| `log_law_crossing(*, kappa, B)` **(+)** | y⁺ where U⁺ = y⁺ meets the log law | E7 |
| `LOG_LAW_CONSTANTS` | dict flow → (kappa, B, citation), public values only | N156 |
| `spalding_yplus(Uplus, *, kappa, B)` | y⁺ | N153 |
| `spalding_uplus(yplus, *, kappa, B)` | U⁺ (brentq) | N153 |
| `coles_wake(xi, kind="cubic")` | W(ξ) | N154 |
| `composite_profile(y, delta, u_star, nu, *, kappa, B, Pi, wake="cubic")` | U [m/s] | N154 |
| `composite_profile_plus(yplus, Re_tau, *, kappa, B, Pi, wake="cubic")` **(+)** | U⁺ | N154 (wall-unit form; E7, F5) |
| `spalding_slope(Uplus, *, kappa, B)` **(+)** | dU⁺/dy⁺ | N153 |
| `velocity_defect(y, U, U_inf, u_star, delta)` | (ξ, (U_∞ − U)/u_*) | (12.84) |
| `total_stress(y, U, uv, mu, rho)` | $\bar\tau=\mu\,dU/dy-\rho\overline{uv}$ [Pa] | (12.76) |
| `channel_total_stress(y, h, tau0)` | $\tau_0(1-2y/h)$ [Pa] | D17 |
| `stress_partition(yplus, Re_tau, *, kappa, B)` | dict total, viscous, reynolds (fractions of τ₀; model partition) | N137 |
| `channel_momentum_residual(y, tau, dPdx)` | residual of $0=-dP/dx+d\bar\tau/dy$ | (12.76) |
| `channel_pressure_gradient(tau0, h)` | $-2\tau_0/h$ [Pa/m] | (12.90) |
| `pipe_pressure_gradient(tau0, d)` | $-4\tau_0/d$ | (12.91) |
| `wall_stress_from_pressure_gradient(dPdx, d)` | τ₀ | (12.91) |
| `boundary_layer_stress_from_profile(x, y, U, V, dPdx, rho)` | $\bar\tau(y)$ | (12.78) |
| `zpg_boundary_layer(x, U_inf, nu, *, kappa)` | dict θ, δ*, δ₉₉, H, C_f, u_*, δ⁺ | N150 |
| `skin_friction_zpg(Re_x, law="monkewitz", kappa=None)` | C_f | N150, N151 |
| `nagib_chauhan_kappa(B)`, `nagib_chauhan_B(kappa)` | float | (12.92) $\kappa B=1.6[\exp(0.1663B)-1]$ |
| `rough_wall_log_law(y, u_star, y0, *, kappa)` | U [m/s] (NaN for y ≤ y₀) | (12.93) |
| `friction_velocity_from_wind(U_ref, z_ref, z0, *, kappa)` | u_* | (12.93) inverted |
| `drag_coefficient_neutral(z_ref, z0, *, kappa)` | $C_D=[\kappa/\ln(z/z_0)]^2$ | D20 (ours) |
| `pipe_bulk_velocity_loglaw(a, u_star, nu, *, kappa, B)` | U_av | N160 |
| `pipe_friction_factor_turbulent(Re_d, kappa=None, B=None)` | Darcy f (Prandtl's law when both None) | N160 |
| `dimensionless_shear(zeta, beta=5.0, unstable="log_linear")` **(+)** | $\phi_m(\zeta)$ | N188 (ours) |
| `surface_layer_wind(z, u_star, z0, L_M=np.inf, *, kappa, beta=5.0, unstable="log_linear")` | U(z) [m/s]; `unstable="businger_dyer"` uses the integrated Businger–Dyer form for L_M < 0 | N188 |

### C.3 `fluidpy/ch12_turbulence.py` (chapter module; re-exports C.1 and C.2)

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `conventions()` **(+)** | pandas frame: the symbol table of convention 4 | decision 6 |
| `book_slips()` **(+)** | list of dicts id, where, printed, corrected, test | §9 slips |
| `time_average_exp_cos(t, A, B, tau, omega, window)` | (average, mean factor sinh(x)/x, wave factor sin(x)/x) | Example 12.1 |
| `divergence_rms(u, v, dx)` | rms of ∇·u | N02 |
| `mean_divergence(samples_u, samples_v, dx, dy)` | field | (12.27), (12.28) |
| `frozen_field_probe(field, U0, u_rms, dt)` | dict: record, reconstruction error | N40 |
| `rans_sympy()` | dict of averaged equations + printed steps | (12.27)–(12.34) |
| `reynolds_stress_budget_sympy()` | dict: residual of one component; trace/2 − (12.47) | (12.35) |
| `tke_budget_sympy()` | dict: each step of D10 and the final residual | (12.47) |
| `mean_energy_budget_sympy()` | dict: residual | (12.46) |
| `temperature_variance_sympy()` | dict: residual; factor of the molecular transport | (12.112) |
| `overlap_matching_sympy()` | dict: (12.85)–(12.89) | D19 |
| `plane_jet_similarity_sympy()` | dict: the three coefficients; residual for a Gaussian F | (12.63) |
| `gradient_moments_isotropic_sympy()` | (2, 4, −1) and 30 | (12.43) |
| `isotropic_tensor_divergence_sympy(f=None)` | 0 | (12.41) |
| `mean_stress_tensor(P, gradU, mu, rho0, uu)` | $\bar\tau_{ij}$ [Pa] | (12.30) |
| `rans_momentum_residual(fields, …, g=0.0)` | residual of (12.30) | (12.30) |
| `rans_2d_residual(U, V, P, uu, uv, vv, x, y, nu, rho)` | (x-, y-residual) | (12.59), (12.60) |
| `rans_eddy_viscosity_residual(…, printed=False)` | residual | (12.97) |
| `mean_heat_flux(gradT, uT, k_th, rho0, cp)` | $Q_j$ [W/m²] | (12.32) |
| `mean_scalar_flux(gradY, uY, kappa_m)` | flux | (12.34) |
| `turbulent_heat_flux(w_rms, T_rms, r_wT, rho, cp)` | H [W/m²] | N54 |
| `mixture_density(v, rho_s, rho)`, `mass_fraction_from_volume_fraction(v, rho_s, rho)` | ρ_m; Y | N55 |
| `displaced_parcel_uv(dUdy, l_rms, v_rms, n=100000, seed=0, correlation=1.0)` | dict: uv, estimate $-\overline{v\ell}\,dU/dy$, r_uv, stderr, samples (u, v, ℓ) | N50 |
| `parcel_uv_expected(dUdy, l_rms, v_rms, correlation=1.0)` **(+)** | −correlation · v_rms · l_rms · dU/dy [m²/s²] | N50 |
| `closure_count(level=3)` | table of unknowns vs equations | N60 |
| `isotropy_report(u, v, w=None, dx=1.0)` | dict of ratios | (12.36), (12.37) |
| `isotropic_correlation_tensor(rvec, f, g=None, u2=1.0, incompressible=True)` | 3 × 3 $R_{ij}$ | (12.40), (12.41) |
| `transverse_from_longitudinal(r, f)` | $g=f+\tfrac r2f'$ | N70 |
| `isotropic_scales(r, f)` | dict Λ_f, Λ_g, λ_f, λ_g and the two ratios | (12.39), N70 |
| `dissipation_rate(grad_u_samples, nu)` | $\bar\varepsilon$ [m²/s³] | (12.42) |
| `dissipation_isotropic(nu, u2, lambda_f=None, lambda_g=None, dudx_sq=None)` | $\bar\varepsilon$ (exactly one of the three) | (12.43) |
| `taylor_reynolds_number(u2, lam, nu)` | $R_\lambda$ | (12.44) |
| `mean_energy_budget(y, U, uv, nu, dPdx=None, rho=1.0)` | dict of terms + residual | (12.46) |
| `tke_budget(y, U, uv, eps, transport=None, buoyancy=None)` | dict production, dissipation, transport, residual | (12.47) |
| `shear_production(uu, gradU)` | $-\overline{u_iu_j}\,\partial U_i/\partial x_j$ [m²/s³] | N81 |
| `reynolds_stress_production(uu, gradU)` | tensor | (12.35) |
| `buoyant_production(wT, alpha, g=G0)` | $g\alpha\overline{wT'}$ | N81 |
| `mean_to_turbulent_dissipation_ratio(Re, urms_over_U=1.0)` | ~1/Re | N80 |
| `mixing_potential_energy_change(z, T_initial, alpha, rho0, g=G0)` | ΔPE [J/m²] | N82 |
| `channel_energy_budget(Re_tau, kappa, A_plus, n=400)` **(+)** | dict of arrays (wall units): yplus, Uplus, uv_plus, pressure_work, mean_dissipation, production, turb_sink (residual, labelled model), integrals | (12.46), (12.47) for U(y) |
| `channel_energy_budget_at(yplus, Re_tau, kappa, A_plus)` **(+)** | dict of floats (same keys, one height) | E5 parity |
| `dissipation_outer_scaling(dU, L, c=1.0)` | $c(\Delta U)^3/L$ | (12.48), (12.49) |
| `kolmogorov_scales(nu, eps)` | (η [m], u_K [m/s], τ_η [s]) | (12.50) |
| `scale_separation(Re_L, c_eps=1.0)` | dict eta_over_L, lambdaT_over_L, uK_over_dU, tau_eta_over_T | (12.51), (12.52) |
| `dns_grid_points(Re_L)` | $\mathrm{Re}_L^{9/4}$ | N89 (ours) |
| `scale_ordering(Re_L, urms_over_dU=1.0)` | dict + bool | N92 |
| `cascade_tiers(L, dU, nu, ratio=2.0)` | dict of arrays size, velocity, turnover, Re | N87 |
| `richardson_diffusivity(l, eps, c=1.0)` | $c\,\bar\varepsilon^{1/3}l^{4/3}$ [m²/s] | N06 |
| `inertial_range_decades(Re_L)` **(+)** | $\log_{10}(L/\eta)=\tfrac34\log_{10}\mathrm{Re}_L$ | (12.51) |
| `scale_table(case)` | dict for our cases "kitchen_mixer", "wind_tunnel", "atmospheric_boundary_layer", "ocean_thermocline" | N90 |
| `kolmogorov_normalize_spectrum(k1, S11, nu, eps)` | ($k_1\eta$, $S_{11}/(u_K^2\eta)$) | (12.53) |
| `inertial_spectrum_1d(k1, eps, C1=None, two_sided=True, printed=False)` | $S_{11}$ [m³/s²] | (12.54), (12.55) |
| `inertial_spectrum_3d(K, eps, C=1.5)` | S(K) | N95 |
| `fit_inertial_range(k, S, band)` | (slope, constant) | (12.54) |
| `kolmogorov_constants(C=1.5)` | dict C, C1_one_sided = 18C/55, C1_two_sided | N95 |
| `model_spectrum(K, eps, nu, L=None, kind="pao")` | E(K) [m³/s²] | N97 |
| `one_dimensional_from_3d(k1, E_fn)` | one-sided $E_{11}(k_1)$ | N95 |
| `general_similarity_check(delta_fn, UCL_fn, Psi_fn, x)` | dict: the three coefficients, ratios, $U_{CL}^2\delta$ | (12.64), (12.74) |
| `free_shear_exponents(flow, return_equations=False)` | dict of `Fraction`: width, velocity, scalar, reynolds (+ the two equations as text) | Table 12.1 exponents (D16) |
| `local_reynolds_number_exponent(flow)` | Fraction | D16 |
| `gaussian_profile(xi, xi_half)` | $\exp\{-\ln2\,\xi^2/\xi_{1/2}^2\}$ | N122 |
| `profile_integrals(xi_half_U, xi_half_Y=None)` | dict I1 = ∫F, I2 = ∫F², IHF | N122 |
| `plane_jet_mean_velocity(x, y, Js, rho, *, C5, xi_half, x0=0.0)` | U [m/s]; `C5="from_invariant"` sets $C_5=I_2^{-1/2}$ | (12.66) |
| `plane_jet_cross_velocity(x, y, Js, rho, *, C5, xi_half, x0=0.0)` | V | (12.58) |
| `plane_jet_stress_profile(xi, F, C3=1.0)` | G(ξ) with $C_3G=-\tfrac12F\int_0^\xi F\,d\xi$ (sign: see C.4) | (12.63) integrated |
| `plane_jet_reynolds_stress(x, y, Js, rho, *, C5, xi_half, C3=1.0)` | $-\overline{uv}$ | (12.67) |
| `plane_jet_volume_flux(x, Js, rho, *, C5, xi_half)` | $\dot V$ [m²/s] | (12.68) |
| `plane_jet_entrainment_velocity(x, Js, rho, *, C5, xi_half)` | $\tfrac12d\dot V/dx$ | D15 |
| `plane_jet_mass_fraction(x, y, Ms, Js, rho, *, C6, xi_half_Y)` | $\bar Y$ | (12.71) |
| `jet_momentum_flux_per_span(y, U, rho)` | $J_s$ [N/m] | (12.62) |
| `scalar_flux_per_span(y, U, Y, rho)` | $\dot M_s$ | (12.70) |
| `slot_momentum_flux(rho_s, U0, d)`, `slot_mass_flux(rho_s, U0, d)` | $\rho_sU_0^2d$; $\rho_sU_0d$ | (12.72), (12.73) |
| `virtual_origin_fit(x, half_width)` | (slope, x₀) | N113 |
| `thin_shear_layer_terms(x, y, Js, rho, nu, *, C5, xi_half)` | dict of term sizes | (12.61) |
| `plane_jet_eddy_viscosity_profile(xi, xi_half)` | sech² shape | N130 |
| `free_shear_flow(flow, x, cross, *, constants, **params)` | dict U, Y, width, centreline | Table 12.1 forms |
| `free_shear_centerline(flow, x, *, constants, **params)` | dict | (12.72), (12.73) |
| `free_shear_profile(flow, xi, *, xi_half)` **(+)** | Gaussian in the flow's similarity variable | N122 |
| `wrong_exponent_fluxes(x, n, m)` **(+)** | (momentum flux, volume flux) relative to x = 1 for $U_{CL}\propto x^n$, δ ∝ x^m: ($x^{2n+m}$, $x^{n+m}$) | (12.62), (12.68) |
| `FREE_SHEAR_CONSTANTS` | dict (no default set until a public table is confirmed) | N123 |
| `stoichiometric_mass_fraction(fuel_M, oxidiser_M, moles_O2_per_fuel, x_O2)` | Y | N126 |
| `round_jet_distance_for_mass_fraction(Y_target, d, rho_s, rho, Y0, C_Y)` | x [m] | N126 |
| `jet_tke_budget(xi, xi_half, C3, model="eddy_viscosity")` | dict of terms (qualitative) | (12.75) |
| `eddy_viscosity_stress(gradU, nu_T, e)` | $\overline{u_iu_j}$ | (12.94) |
| `eddy_viscosity_from_data(uv, dUdy)` | ν_T | (12.94) inverted |
| `gradient_diffusion_flux(grad, K)` | −K grad | (12.95), (12.96) |
| `eddy_diffusivity_estimate(l_T, u_T, c=1.0)` | $c\,l_Tu_T$ | (12.98) |
| `mixing_length_stress(dUdy, l_T)` | $-\overline{uv}=l_T^2\lvert dU/dy\rvert\,dU/dy$ | N169 |
| `mixing_length_eddy_viscosity(dUdy, l_T)` | $l_T^2\lvert dU/dy\rvert$ | N169 |
| `mixing_length_wall_profile(yplus, kappa, damping=None, A_plus=26.0)` | dict Uplus, slope, lT_plus, nuT_over_nu, uv_plus (floats for a float y⁺) | (12.100), (12.101) |
| `mixing_length_intercept(kappa, A_plus=None)` **(+)** | B; without damping exactly $[\ln(4\kappa)-1]/\kappa$ | D21 (ours) |
| `shear_flow_eddy_viscosity_solve(y, nu_T_fn, dPdx, rho, nu, bc)` | U(y) | (12.99) |
| `channel_mixing_length(Re_tau, kappa, A_plus=26.0, n=400, core_cap=0.09)` | dict yplus, Uplus, uv_plus, production, Re_bulk, Cf | N181 |
| `convective_velocity_scale(L, dT, T, g=G0)`, `convective_eddy_diffusivity(L, dT, T, g=G0)` | w; κ_T | (12.102) |
| `one_equation_closure(e, l_T, c, C_eps, sigma_e)` | dict nu_T, eps | N175 |
| `k_epsilon_eddy_viscosity(e, eps, C_mu=0.09)` | $C_\mu\bar e^2/\bar\varepsilon$ | (12.104) |
| `k_epsilon_length_scale(e, eps)` | $\bar e^{3/2}/\bar\varepsilon$ | (12.104) |
| `k_epsilon_rhs(e, eps, production, constants=K_EPSILON_CONSTANTS)` | (dē/dt, dε̄/dt) sources | (12.103), (12.105) |
| `K_EPSILON_CONSTANTS` | dict (Launder & Sharma 1974) | N178 |
| `k_epsilon_decay(e0, eps0, t, C_eps2=1.92)` | (ē(t), ε̄(t), n, t₀) | D23 |
| `k_epsilon_loglayer_kappa(C_mu, C_eps1, C_eps2, sigma_eps)` | κ | D23 |
| `stratified_tke_budget(z, U, uw, wT, eps, alpha, g=G0)` | dict terms + residual | (12.106) |
| `flux_richardson(wT, uw, dUdz, alpha, g=G0)` | Rf (NaN at zero shear production) | (12.107) |
| `turbulence_regime(Rf, Rf_cr=0.25)` | "convective" (Rf < 0) / "shear-driven" (0 ≤ Rf < Rf_cr) / "decaying" | N183 |
| `gradient_richardson_thermal(dTdz, dUdz, alpha, g=G0, *, Gamma_a, convention="kundu")` — `Gamma_a` **required** (review M1): −g/C_p ≈ −9.76e-3 K/m for an in-situ gradient, `0.0` only for a potential-temperature gradient | dict Ri, N2, dthetadz, verdict_kundu, verdict_met | (12.108) |
| `turbulent_prandtl(nu_T, kappa_T)`, `flux_from_gradient_richardson(Ri, Pr_T)` | float | (12.109) |
| `monin_obukhov_length(u_star, wT, alpha=None, T=None, *, kappa, g=G0)` | L_M [m] (±inf when wT = 0) | (12.110) |
| `monin_obukhov_from_fluxes(tau, H, rho, cp, T, *, kappa, g=G0)` | L_M | (12.110) |
| `flux_richardson_surface_layer(z, L_M)` | z/L_M | (12.111) |
| `gradient_richardson_surface_layer(z, L_M, Pr_T=1.0, beta=5.0)` **(+)** | Ri(z) | (12.109) + (12.111) |
| `surface_layer_regime(z, L_M)` | "forced convection" / "free convection" / "stable" / "neutral" | N187 |
| `surface_layer_state(u_star, H, T, z0, z, rho, cp, kappa, beta=5.0, Pr_T=1.0, convention="kundu")` **(+)** | dict wT, L_M, Rf, Ri, phi_m, U, U_neutral, regime, z_crit, verdict_kundu, verdict_met | (12.107)–(12.111) |
| `temperature_variance_budget(z, T_mean, wT, eps_T)` | dict | (12.112) |
| `scalar_spectrum(K, eps, eps_T, nu, kappa_th, C_T=1.0)` | $S_T$ | (12.113), (12.114) |
| `batchelor_scale(nu, kappa_th, eps)` | η_T | N193 |
| `langevin_particles(n, t, u_rms, Lambda_t, seed=0, dim=1)` | (X, u) arrays (n, nt) | N196 |
| `dispersion_rate_from_particles(t, X, u)` | (d⟨X²⟩/dt, 2⟨Xu⟩) | (12.115), (12.116) |
| `taylor_dispersion_rate(t, r_fn, u2)` | $2\overline{u^2}\int_0^tr\,d\tau$ | (12.117) |
| `taylor_dispersion(t, r_fn, u2, form="single")` | $\overline{X^2}$ ("double" = (12.118)) | (12.118), (12.119) |
| `taylor_dispersion_exponential(t, u2, Lambda_t)` | $2\overline{u^2}\Lambda_t^2[t/\Lambda_t-1+e^{-t/\Lambda_t}]$ | N208 |
| `taylor_dispersion_gaussian(t, u2, t_c)` | $2\overline{u^2}\big[\tfrac{\sqrt\pi}2t_c\,t\,\mathrm{erf}(t/t_c)-\tfrac{t_c^2}2(1-e^{-t^2/t_c^2})\big]$ | N208 |
| `dispersion_regime(t, Lambda_t)` | "ballistic" (t < 0.3Λ) / "transition" / "diffusive" (t > 3Λ) | (12.121), (12.123) |
| `dispersion_local_slope(t, Lambda_t)` **(+)** | d ln⟨X²⟩/d ln t (2 → 1) | E10 |
| `eddy_diffusivity_taylor(t, r_fn, u2)` | D_T [m²/s] | (12.127) |
| `eddy_diffusivity_exponential(t, u2, Lambda_t)` **(+)** | $\overline{u^2}\Lambda_t(1-e^{-t/\Lambda_t})$ | (12.127) |
| `random_walk(n_steps, n_walkers, L=1.0, dim=2, persistence=0.0, seed=0)` | positions (walkers, steps + 1, dim) | (12.124), (12.125) |
| `smoke_plume_width(x, U, w_rms, Lambda_t)` | Z_rms(x) with t = x/U | N212 |
| `plume_concentration(x, z, Q, U, sigma_z)` | Gaussian cross-section | N212 (ours) |
| `diffusivity_from_variance(t, variance)` | ½ dσ²/dt | (12.126) |

`k_epsilon_channel` (analysis §4 row 39) is **not called** by the notebook or any explainer (optional; not part of this contract).

### C.3b `scripts/ch12_*.py` (runnable demos; `--no-show`)
As analysis §4: `ch12_averaging.py`, `ch12_correlation_spectrum.py`, `ch12_reynolds_stress.py`, `ch12_isotropic.py`,
`ch12_cascade_spectrum.py`, `ch12_plane_jet.py`, `ch12_free_shear_table.py`, `ch12_wall_layers.py`, `ch12_mixing_length.py`,
`ch12_k_epsilon.py`, `ch12_stratified.py`, `ch12_dispersion.py`; plus **(+)** `scripts/ch12_drawings.py` (pure drawing helpers the
notebook imports: `scales_sketch`, `parcel_sketch`, `stress_element_sketch`, `fg_geometry_sketch`, `free_shear_sketch`,
`wall_layers_sketch`, `surface_layer_sketch`, `plume_sketch` — no physics) and **(+)** `scripts/ch12_tables.py` (writes
`reference/ch12/explainer_tables.json`: the `channel_energy_budget` table for E5 at Re_τ = 180, 550, 1000, 5200 and the
van Driest U⁺(y⁺) table for E8, "ours, computed by fluidpy").

### C.4 Not in `analysis/ch12.md` §4 (flag list for the implementer)
1. From curation §8 (already known to the implementer): `TS.correlation_spectrum_pair`, `TS.smooth_signal`,
   `ch12.parcel_uv_expected`, `ch12.channel_energy_budget`, `ch12.inertial_range_decades`, `ch12.scale_table` cases,
   `ch12.wrong_exponent_fluxes`, `ch12.free_shear_profile`, `ch12.mixing_length_intercept`, `WT.log_law_crossing`,
   `ch12.surface_layer_state`, `ch12.gradient_richardson_surface_layer`, `ch12.eddy_diffusivity_exponential`,
   `ch12.dispersion_local_slope`, `ch12.book_slips`, `ch12.conventions`, `free_shear_exponents(return_equations=True)`.
2. **New in this design (one function, two scripts):** `ch12.channel_energy_budget_at(yplus, Re_tau, kappa, A_plus)` (E5 parity and inspector), `scripts/ch12_drawings.py`, `scripts/ch12_tables.py`. **Adopted from the implementer's files as they stood while this design was written** (`core/turbstats.py`, `core/wall_turbulence.py`; not in analysis §4 either): `TS.standard_error_of_mean`, `TS.velocity_covariance`, `TS.shell_spectrum`, `WT.composite_profile_plus`, `WT.spalding_slope`, `WT.dimensionless_shear`; `fluidpy/ch12_turbulence.py` did not exist yet, so every C.3 signature is this design's reading of analysis §4 and curation §8 and must be reconciled with the implementer's file in the verify phase. The re-export of `TS` and `WT` from `ch12` (convention 1) is also a design requirement.
3. **Return shapes fixed here** (analysis §4 leaves them open): `mixing_length_wall_profile` → dict; `kolmogorov_scales` →
   3-tuple; `gradient_richardson_thermal` → dict with both verdict strings; `correlation_spectrum_pair` → dict;
   `displaced_parcel_uv` → dict; `free_shear_exponents` → dict of `Fraction` (parity rows multiply by 1.0);
   `time_average_exp_cos` → 3-tuple; `layer_name`, `turbulence_regime`, `surface_layer_regime`, `dispersion_regime` → the
   exact strings in the tables above (explainers pin them with exact-text rows).
4. **A sign to settle with the analysis (row #121).** Integrating (12.63),
   $\{\delta U'_{CL}/U_{CL}\}F^2-\{\delta U'_{CL}/U_{CL}+\delta'\}F'\int_0^\xi F\,d\xi=\{\Psi/U_{CL}^2\}G'$, with $C_1=-\tfrac12$, $C_2=+\tfrac12$
   gives $C_3G=-\tfrac12F\int_0^\xi F\,d\xi$ (D15 step 8): G < 0 for ξ > 0, so $-\overline{uv}=\Psi G<0$ where $dU/dy<0$ — the sign an
   eddy viscosity gives. Analysis row #121 writes the right side without the minus sign; the design uses the derived sign and
   asks the verifier to confirm it with `plane_jet_similarity_sympy`.
5. `surface_layer_wind(unstable="log_linear")` returns a non-physical (even negative) wind for strongly unstable z/∣L_M∣ ≳ 1
   (*expect* U(10 m) = −0.16 m/s for u_* = 0.3, H = +300 W/m², z₀ = 0.03 m); E9 and F7 therefore use
   `unstable="businger_dyer"` on the unstable side and show the book's log-linear form as a dashed line labelled "outside its
   range".

### C.5 Parity conventions (the explainers reproduce these exactly)
1. Stochastic stages never use a sampled value for parity: E1 the window factors and σ/√n; E3 `parcel_uv_expected`; E10 the
   exponential closed forms.
2. Spalding inverse: Newton from the log-law value, 30 iterations, agrees with `spalding_uplus` (brentq) to 1e-9.
3. van Driest profile (E5, E8): JS integrates $dU^+/dy^+=2/(1+\sqrt{1+4l^{+2}})$, $l^+=\kappa y^+(1-e^{-y^+/A^+})$, by Simpson on a
   log grid (y⁺ from 1e-4, 800 intervals per decade ⇒ rtol 1e-5 against `mixing_length_wall_profile`); the intercept is read at
   y⁺ = 2 × 10⁴.
4. Businger–Dyer (E9): $\phi_m=(1-16\zeta)^{-1/4}$ for ζ < 0, $\psi_m=2\ln\frac{1+x}2+\ln\frac{1+x^2}2-2\arctan x+\frac\pi2$,
   $x=(1-16\zeta)^{1/4}$, $U=\frac{u_*}\kappa[\ln\frac z{z_0}-\psi_m(\zeta)]$ (AMS Glossary form; cited, not from the book).
5. Verdict strings (E9) come from `surface_layer_state(...)["verdict_kundu"]` / `["verdict_met"]` and are pinned by exact-text rows.

### C.6 Audit of the contract (2026-10-07) — who calls what, and what must be reconciled
1. **Every name in C.1–C.3 is now called or named by Part A, B or F.** Forty-four names of the contract were in no storyboard
   row although the curation's chapter map names most of them as the way a NOTE is stated; A.15a attaches forty-three of
   them to their notes with arguments and an *expect*. One name stays **contract-only on purpose**: `ch12.rans_momentum_residual` (analysis §4;
   evidence for the Reynolds-averaged momentum equation by a manufactured solution in `tests/test_ch12.py` — its dict of
   callables is too heavy for a teaching cell; the notebook shows `rans_sympy` instead). `k_epsilon_channel` is outside the
   contract (line above C.3b).
2. **No function is added, renamed or re-signed by this audit.** Part C is unchanged apart from this subsection.
3. **Contract vs the implementer's files** (checked by importing `fluidpy.ch12_turbulence` and comparing
   `inspect.signature` with the 212 signatures above that carry no "…"). At 09:00 on 2026-10-07 six differences stood
   (`channel_energy_budget_at` missing; `plane_jet_mass_fraction`, `plane_jet_reynolds_stress`,
   `rans_eddy_viscosity_residual(printed=…)`, `mean_heat_flux` signatures; `batchelor_scale` raising a `TypeError`); the
   implementer closed all six while this audit ran. **State at 09:06: every C.1–C.3 name exists on `ch12` and no signature
   contradicts the contract.** Three things a builder should still know:
   - `WT.total_stress` returns a dict (total, viscous, reynolds) where C.2 says "τ̄ [Pa]"; builders use `["total"]`.
   - `ch12.free_shear_flow` and `ch12.free_shear_centerline` spell the contract's `**params` out as keywords (`d`, `U0`,
     `rho_s`, `rho`, `Y0`, `g`, `theta`, `U_inf`, …) and ask for their own key names inside `constants`: read the
     docstring; the values stay labelled illustrative (convention 10).
   - `ch12.surface_layer_regime(z, L_M)` returns the four contract strings, with the book's meaning of "forced
     convection": any height well below ∣L_M∣, **stable or unstable** (so z = 10 m, L_M = +100 m is "forced convection";
     "stable" is returned for z ≳ L_M > 0). E9's status must not read "forced convection" as "unstable".

---

## Part A — notebook storyboard (`notebooks/build_ch12.py` → `notebooks/ch12_turbulence.ipynb`)

**One line per book section** (cell numbers are estimates; ≈ 640 cells in all):
- §12.1 → N01 N02 N03 N04 (open the C01 story; P286) — cells ≈ 10–22
- §12.2 → N05 (timeline), N06 named (stated at the end of C16) — cells ≈ 23–26
- §12.3 → C01 (N07–N23 · D01 · P280–P285 · A1 · F1 · figs · **E1**) — cells ≈ 27–85
- §12.4 → C02 (N24–N34 · D02 D03 D04 · P287 P288 P289 · A2 · figs), C03 (N35–N40 · D05 · P290–P293 · F2 · figs · **E2**) — cells ≈ 86–165
- §12.5 → R01 R02, C04 (N41–N61 · D06 · P294 P295 P296 · A3 · figs · **E3**) — cells ≈ 166–215
- §12.6 → C05 (N62–N77 · D07 D08 · P297 · figs) — cells ≈ 216–265
- §12.7 → C06 (N78–N83 · D09 D10 · fig · **E5**), C07 (N84–N92 · D11 · A4 · fig), C08 (N93–N99 · D12 · P298 · F3 · live · **E4**) — cells ≈ 266–345
- §12.8 → C09 (N100–N130 · D13 D14 D15 D16 · F4 · figs · **E6**) — cells ≈ 346–410
- §12.9 → R03 R04, C10 (N131–N142 · D17 D18 · P299 P300 · figs), C11 (N143–N161 · D19 D20 · P301 · F5 · figs · **E7**) — cells ≈ 411–480
- §12.10 → C12 (N162–N174 N181 · D21 · P302 · F6 · figs · **E8**), C13 (N175–N180 · D22 D23 · P303 · fig) — cells ≈ 481–535
- §12.11 → R05 R06, C14 (N182–N185 · D24 · fig), C15 (N186–N195 · D25 · P304 · F7 · live · figs · **E9**) — cells ≈ 536–585
- §12.12 → R07 R08, C16 (N196–N216 N06 · D26 D27 D28 · P305 P306 · A5 · F8 · figs · **E10**) — cells ≈ 586–632
- §12.13 → **C16 (continued)**: N217; S01, S02; summary — cells ≈ 633–640

Every CORE block follows: problem in plain words → idea → primers → maths (notes and derivations, Part F) → tiny example →
code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would change if…". Code
drafts give intent + exact calls; the builder comments every line (novice grade, units, the equation written out next to its
number). *expect* = numbers the executed cell must print. *see / read / change* = the three figure notes. Note ids open every
note in bold (**N09 [B]**). Wherever a draft below names an equation, the equation is written out with its number, and the
builder keeps it so. **A.15 is part of every block:** before building a CORE block, read its rows in A.15a (function calls the
curation attaches to that block's notes) and apply A.15b (repeat mentions of an equation by number).

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…16…], prerequisites=[…])`. **Big idea (draft):** "Stir milk into coffee and in two seconds
   the cup is uniform; without stirring, molecular diffusion would need hours. The difference is turbulence — a flow so
   irregular that no one can predict where a given drop will be, yet so regular *on average* that engineers design aircraft
   with it and climate models run on it. This chapter is about that 'on average'. We split every field into a mean and a
   fluctuation, average the equations, and find that the fluctuations push back on the mean through one new term — a
   correlation, the Reynolds stress — for which there is no equation. Everything else is a way of living with that gap:
   measuring the fluctuations (correlations, spectra), following their energy from the big eddies that take it from the mean
   flow down to millimetre eddies that turn it into heat (the cascade and Kolmogorov's −5/3 law), using symmetry and
   conservation where they are enough (jets, the logarithmic law near a wall), guessing the missing term where they are not
   (eddy viscosity, mixing length, k–ε), adding buoyancy (Richardson numbers, the Monin–Obukhov length — the vocabulary of
   every surface-flux scheme in a climate model), and finally asking how turbulence spreads things (Taylor's dispersion)."
   **Roadmap (one line per CORE):** C01 averages and the one rule that fails · C02 the autocorrelation and the memory time ·
   C03 the spectrum is the same information · C04 Reynolds-averaged equations and the Reynolds stress · C05 isotropic
   turbulence: one function, one gradient · C06 the two energy budgets and the term they share · C07 Kolmogorov scales and the
   cascade · C08 the −5/3 law · C09 the self-similar jet · C10 wall units and the law of the wall · C11 the logarithmic law ·
   C12 eddy viscosity and the mixing length · C13 the k–ε model · C14 flux and gradient Richardson numbers · C15 the
   Monin–Obukhov length · C16 Taylor's dispersion. **Prerequisites:** index notation, δ_ij, trace, symmetric contraction
   (Ch. 2), Lagrangian vs Eulerian (Ch. 3), the Boussinesq equations and the dissipation $\varepsilon=2\nu S_{ij}S_{ij}$ (4.58)
   (Ch. 4), vortex stretching (Ch. 5), similarity solutions and the laminar jet (Ch. 8–9), thin-layer scaling (Ch. 9), the Π
   theorem, potential temperature and the lapse-rate conventions (Ch. 1), the disturbance-energy budget and Ri (Ch. 11).
2. `nb.explainer_index([...])` — 10 rows: ("reynolds_averaging_window", "How long must you average to get 'the mean'?", "a mean
   is only as good as its number of independent samples: N members, or Δt/Λ_t memory times") · ("correlation_and_spectrum",
   "Why are a correlation and a spectrum the same information?", "stretch the memory and the spectrum squeezes; its height at
   zero is the integral scale, its area the variance") · ("reynolds_stress_parcels", "How can fluctuations that average to zero
   push the mean flow?", "going up correlates with being slow: ⟨uv⟩ < 0 in a positive shear") · ("energy_cascade_spectrum",
   "What does a higher Reynolds number change?", "the supply ΔU³/L is fixed by the big eddies; viscosity only sets how far down
   the ladder the energy must go") · ("turbulent_energy_budget", "Where does turbulent energy come from and go?", "one term,
   two signs: the mean flow's loss is the turbulence's income") · ("turbulent_jet_similarity", "How does a jet spread without a
   turbulence model?", "one invariant + shape-preservation fix both exponents") · ("law_of_the_wall", "Where does the logarithm
   come from?", "it is what is left when the answer may depend on neither the viscous nor the outer length") ·
   ("mixing_length_closure", "What does a closure constant do?", "κ sets the slope, the wall damping sets the intercept") ·
   ("stratified_surface_layer", "When does stratification kill turbulence?", "Rf = z/L_M: below ∣L_M∣ shear rules, above it
   buoyancy") · ("taylor_dispersion", "Why t first and √t later?", "a random walk whose step is the particle's memory").
3. `nb.setup()`.
4. `nb.code` — **chapter imports**: `import numpy as np`, `import sympy as sp`, `import matplotlib.pyplot as plt`, `import pandas
   as pd`, `from fractions import Fraction`, `from scipy import signal`, `from scipy.integrate import quad, solve_ivp,
   cumulative_trapezoid`, `from scipy.optimize import brentq`, `from fluidpy import ch12_turbulence as ch12`, `from fluidpy.core
   import turbstats as TS, wall_turbulence as WT`, `from fluidpy.core import dimensional as DIM, stratification as STRAT, laminar
   as LAM, boundary_layer as BL, jets as JET, diffusion as DIFF`, `from fluidpy import ch11_instability as ch11`, `from
   fluidpy.core.interact import slider_figure, animate_figure, live`, `from fluidpy.core.anim import animate`, `from
   fluidpy.core.style import COLORS, savefig`, `from fluidpy.core.thermo import G0`, `import sys; sys.path.insert(0, "scripts");
   from ch12_drawings import *`. *explain:* one line per import ("`TS` and `WT` are the two new core modules; `ch12` re-exports
   both, so `ch12.log_law` and `WT.log_law` are the same function").
5. `nb.md` — **⚠️ Conventions in this chapter** (the conventions block of curation decision 6, **before C01**): (a)
   `ch12.conventions()` shown as the symbol table of convention 4 (κ, k / K, e / E / ε, λ / Λ, η, f / g / F / G, R_ij, τ, δ, θ, α,
   S, N, L, σ, Π, U₀, h / d — book meanings · what we write · code name); (b) "**u is a fluctuation**: tilde = total, capital
   or over-bar = mean, lower case or prime = fluctuation: $\tilde u_i=U_i+u_i$ (12.24)"; (c) the normalisations of convention 5
   as a five-row table; (d) the three assumption switches of convention 6; (e) the two lapse-rate conventions as a two-row
   table with a worked conversion (an atmosphere cooling 6.5 K/km with height: Kundu $dT/dz=-6.5>\Gamma_a=-9.8$ K/km ⇔
   meteorology $\Gamma_{met}=6.5<\Gamma_d=9.8$ K/km ⇒ **stable** — negate the number, flip the inequality, P48); (f)
   `pd.DataFrame(ch12.book_slips())` as "slip #k · the book prints · correct · where we fix it".
6. `nb.md` — **🔁 Tools from earlier chapters used in this one** (one line each): seeded generators (P10), Gaussian (P11),
   1/√N scatter (P12), log–log slopes and `np.polyfit` (P13), `assert np.allclose` (P15), animate (P16), `slider_figure` (P17),
   `show_viz` (P18), `np.gradient` (P22), ∂ (P25), Taylor series (P26/P98), `solve_ivp` (P31/P94), trapezoid/Simpson (P37/P203),
   product rule (P38), sympy (P40/P117), Euler's formula (P45), inequalities under a sign change (P48), chain rule (P49),
   `Fraction` (P60), `np.einsum` (P62), iterated integrals (P83), fundamental theorem of calculus (P84), `quad` (P87),
   substitution (P106), `brentq` (P108), differentiation under the integral (P109), isotropic tensors (P119/P120),
   order-of-magnitude scaling (P130), scaled variables (P133), Fourier modes and the FFT (P142), improper integrals (P144),
   quadratic formula (P159), separation argument (P167), sinh/cosh (P168), sum-to-product (P171), Fourier integral (P172),
   `np.fft.rfft` (P181), two-length scaling (P188), Leibniz rule (P189), `cumulative_trapezoid` (P190), Picard iteration (P192),
   Gaussian integral (P194), erf (P195), exponent matching (P197), chain rule with a moving similarity variable (P206),
   integration by parts (P218a), caching (P252), necessary vs sufficient (P255), even and odd functions (P261); Π theorem
   (ch01 C69), potential temperature and N² (ch01), summation convention, δ_ij, trace, symmetric contraction, principal axes
   (ch02), Lagrangian description (ch03), flux form and Boussinesq (ch04), thin-layer scaling and thicknesses (ch09).

---

### A.1 §12.1 Introduction — N01 N02 N03 N04 (open the C01 story)
1. `nb.section("12.1", "Introduction", intro="**What is this section about?** What 'turbulent' means: not just 'irregular',
   but five properties that come together. Then the plan — since no one can predict a single turbulent flow in detail, we
   predict its statistics.")`
2. `nb.note` — **N01 [B]** "**Five marks of turbulence**" as a five-row table (mark · what it means · which block quantifies
   it): fluctuations (irregular in space and time → C01–C03) · nonlinearity (appears above a critical Re, Ra or 1/Ri, Ch. 11;
   the nonlinear term makes the Reynolds stress → C04) · three-dimensional fluctuating vorticity, eddies of many sizes
   (vortex stretching, $D\boldsymbol\omega/Dt=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega$ (5.13) → C07) ·
   dissipation (it dies without an energy supply → C06) · diffusivity (rapid mixing → C12, C16). "Random waves on a pond are
   irregular but not dissipative or mixing: not turbulence."
3. `nb.primer("random-phase synthetic fields (np.fft.ifftn)", "A field that looks turbulent can be built by giving every
   Fourier mode a chosen amplitude and a random phase, then transforming back. Building the velocity from a stream function
   makes it divergence-free. It has the spectrum we chose but no dynamics — no cascade — so every caption says 'kinematic'.",
   code="rng = np.random.default_rng(0)                 # seeded: same field every run\nphase = np.exp(2j*np.pi*rng.random((8, 8)))    # unit amplitude, random phase per mode\nfield = np.fft.ifft2(phase).real               # back to physical space\nprint(field.shape, round(field.std(), 3))")` (**P286**)
4. `nb.note` — **N02 [B]** "A turbulent velocity field still obeys mass conservation, $\nabla\cdot\mathbf u=0$; white noise does
   not." + `nb.figure`: two panels, `TS.synthetic_solenoidal_field(128, 1.0, spectrum=lambda K: K**4*np.exp(-(K/20)**2),
   seed=1)` (quiver over speed colour) vs `TS.white_noise_field(128, seed=1)`; titles carry `ch12.divergence_rms`. *expect:*
   divergence rms ≈ 1e-13 (round-off) vs O(100). *see:* swirls vs salt-and-pepper; *read:* eddies are what continuity looks
   like; *change:* "…we move the spectrum's peak to higher K: smaller eddies, still divergence-free".
5. `nb.note` — **N03 [C]** "In a turbulent boundary layer the largest eddies are as big as the layer, $l\sim\delta$, and the
   edge between turbulent and irrotational fluid is ragged (sketch: the synthetic field masked by a wavy line). The outer
   scale L returns in C07, entrainment across that edge in C09, δ in C10."
6. `nb.note` — **N04 [C]** "**Scope.** Incompressible, no Coriolis force, three-dimensional fluctuations. Rotating stratified
   flows have nearly two-dimensional (geostrophic) turbulence, which sends energy to *larger* scales — the opposite of C07 —
   and waits for Ch. 13."

### A.2 §12.2 Historical Notes — N05, N06 named
1. `nb.section("12.2", "Historical Notes", intro="**What is this section about?** Who found what — as a map of the chapter.")`
2. `nb.note` — **N05 [C]** timeline table (year · name · idea · block): 1883 Reynolds, pipe experiment and averaging → C04 ·
   1915–1921 Taylor, eddy transport and dispersion → C16 · 1925 Prandtl, mixing length; 1930 von Kármán, logarithmic profile
   → C12, C11 · 1922–1926 Richardson, cascade and the 4/3 law → C07, C16 · 1935 Taylor, isotropic turbulence → C05 · 1941
   Kolmogorov and Obukhov, universal small scales → C07, C08 · 1954 Monin and Obukhov → C15.
3. `nb.md` — "**N06** (Richardson's four-thirds law, $K\sim\bar\varepsilon^{1/3}l^{4/3}$) is named here and stated with a number at
   the end of C16, where both of its ingredients exist."

---

### A.3 §12.3 Nomenclature and Statistics for Turbulent Flow — C01
#### C01 — Ensemble averages (12.1), the time average (12.2) and the commutation rules (12.4)–(12.9)
1. `nb.section("12.3", "Nomenclature and Statistics for Turbulent Flow", intro="**What is this section about?** What 'the
   mean' of a turbulent signal is, how to compute it from many runs or from one long record, which operations an average
   passes through — and the one it does not, which is where the rest of the chapter comes from.")`
2. `nb.core("C01", "Ensemble averages $\\langle u^m\\rangle=\\lim_{N\\to\\infty}\\frac1N\\sum_{n=1}^N(u(\\mathbf x,t{:}n))^m$ (12.1) and
   the rules of averaging", question="What is 'the mean' of a signal that never repeats, and what may an average be moved
   through?")`
3. `nb.md` — **The problem in plain words:** "A wind gauge on a mast reads 4.1, 6.3, 3.8, 5.5 m/s within four seconds. What
   is 'the wind'? A climate normal is a 30-year average for the same reason: the weather never repeats, so we describe it by
   its mean and by how it scatters about the mean. In a laboratory you can repeat an experiment a hundred times and average
   the hundred runs at the same instant; outdoors you have one run and can only average over time. When do the two agree?"
4. `nb.md` — **The idea** (ASCII):
   ```
   realization 1:  ~~~/\~~\/~~~      average ACROSS runs at one time   → ensemble average (theory uses this)
   realization 2:  ~\/~~~/\~~~~      average ALONG one run over Δt     → time average (measurements use this)
   realization N:  ~~~\/~/\~~~~      agree when: stationary, and Δt ≫ memory, Δt ≪ drift
   an average is a SUM (÷N):  it passes through +, ×constant, ∂/∂t, ∂/∂x, ∫   — all linear
                              it does NOT pass through a product:  mean(ũṽ) = Ū V̄ + mean(uv)
   ```
5. `nb.primer("random variable, probability density and histogram", …)` (**P280**) — "A random variable takes a different
   value in every realization; its probability density says how often each value occurs; a histogram (`ax.hist`,
   `np.histogram`) estimates it. A Gaussian has skewness 0 and (normalised) kurtosis 3." code: 4 lines sampling 10⁴ Gaussian
   numbers, printing mean, std and `np.histogram(x, bins=5)[0]`.
6. `nb.primer("Ornstein–Uhlenbeck signal (a Langevin equation)", …)` (**P283**) — "Our test signal: a random signal with a
   chosen standard deviation σ and memory time τ_c. Each step keeps a fraction $e^{-\Delta t/\tau_c}$ of the old value and adds
   fresh noise: $u_{n+1}=u_ne^{-\Delta t/\tau_c}+\sigma\sqrt{1-e^{-2\Delta t/\tau_c}}\,\xi_n$. The update is exact for any Δt." code: 4-line loop,
   σ = 1, τ_c = 0.5 s, printing the std (≈ 1).
7. `nb.note` — **N07 [B]** "**Realization, ensemble, ensemble average.** One run is a realization; the set of all runs under
   the same conditions is the ensemble; the over-bar is the average over N of them and ⟨ ⟩ its limit N → ∞. Chapter 11's Lorenz
   system showed why we must: two runs that start 10⁻⁸ apart end up different, so only statistics are reproducible." + `nb.code`:
   `ens = TS.make_ensemble(64, t, mean_fn=lambda t: np.exp(-t/10), sigma=0.3, tau_c=0.5, seed=0)`; `TS.ensemble_average(ens)`.
8. `nb.md` — the equations: $\langle u^m(\mathbf x,t)\rangle=\lim_{N\to\infty}\frac1N\sum_{n=1}^N\big(u(\mathbf x,t{:}n)\big)^m$ (12.1), with "t:n" read
   "time t in realization n"; **N19 [B]** the mean is the m = 1 case, $\overline{u(\mathbf x,t)}\equiv\frac1N\sum_{n=1}^Nu(\mathbf x,t{:}n)$ (12.10).
9. `nb.note` — **N08 [B]** stationary = statistics independent of the time origin; then
   $\overline{u^m(\mathbf x)}=\frac1{\Delta t}\int_{t-\Delta t/2}^{t+\Delta t/2}u^m(\mathbf x,t)\,dt$ (12.2). + `nb.primer("sliding-window average with np.cumsum", …)`
   (**P284**): "the sum over a window = difference of two cumulative sums; within half a window of each end there is no
   full window (NaN)"; code: 4 lines on `np.arange(10.)` with window 3.
10. `nb.note` — **N09 [B]** homogeneous = stationary in space; $\overline{u^m(t)}=\frac1V\int_Vu^m(\mathbf x,t)\,dV$ (12.3); one number:
    `TS.volume_average(u_field, m=2)` of the N02 field.
11. `nb.note` — **N10 [B]** + `nb.figure` (our Fig. 12.2): two members, constant mean vs decaying mean, the mean dashed purple,
    fluctuations teal. *see / read / change* ("…the decay time were shorter than the memory: no window could separate them").
12. `nb.primer("ergodicity", …)` (**P282**): "For a stationary signal, a long time average of one run equals the ensemble
    average. The book assumes it silently every time a measurement is compared with theory." code: one OU run of 2000 s, its
    time mean (≈ 0 ± 0.02) vs the ensemble mean of 2000 short runs at one time.
13. `nb.note` — **N11 [B]** "Outdoors we average over time (or space): the window must be long against the memory and short
    against the drift, $t_c\ll\Delta t\ll$ drift time. **Climate hook:** a 30-year normal is such a window; an 'eddy flux'
    $\overline{v'T'}$ in the general circulation is a covariance of deviations from a time or zonal mean (Ch. 13)."
14. `nb.derivation("D01", …)` — Part F D01 (8 steps), ref "12.4". Carries **N12** $\overline{u^m+v^m}=\overline{u^m}+\overline{v^m}$ (12.4), **N13** $\overline{Au^m}=A\overline{u^m}$ (12.5), **N14** $\overline{\partial u^m/\partial t}=\frac{\partial}{\partial t}\overline{u^m}$ (12.6), **N15**
    $\overline{\int_a^bu^m\,dt}=\int_a^b\overline{u^m}\,dt$ (12.7), **N16** $\overline{\partial u^m/\partial x_j}=\frac{\partial}{\partial x_j}\overline{u^m}$ (12.8), **N17** $\overline{\int u^m\,d\mathbf x}=\int\overline{u^m}\,d\mathbf x$ (12.9), **N18** (the product:
    $\overline{u^m}\neq\bar u^m$ and $\overline{uv}\neq\bar u\,\bar v$ in general). Each note is one line under the derivation: e.g. **N14 [B]**
    "$\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t$ (12.6): ⚠️ exact for ensemble averages; for a finite time window it holds only when the
    window is short against the drift. Used again in C16."
15. `nb.code` — `res = TS.check_averaging_rules(ens_u, ens_v, t, A=2.0)`; print the dict. *expect:* every entry ≲ 1e-13 except
    `"product"` (≈ the covariance, ≠ 0). `TS.product_average_split(ens_u, ens_v)` → (mean product, product of means,
    covariance) with first = second + third to round-off. *explain:* 1. the same discrete operator is applied to each member
    and to the mean; 2. only the product leaves a remainder.
16. `nb.worked_example("two realizations are enough to break the product rule", "Run 1: ũ = 1, ṽ = 2. Run 2: ũ = 3, ṽ = 6. 1.
    Means: Ū = 2, V̄ = 4, so Ū V̄ = 8. 2. Mean of the product: (1·2 + 3·6)/2 = 10. 3. Fluctuations: u = −1, +1; v = −2, +2; mean of
    uv = (2 + 2)/2 = 2. 4. Check $\\overline{\\tilde u\\tilde v}=\\bar u\\bar v+\\overline{uv}$: 10 = 8 + 2 ✓. The 2 that is left over is the kind of term the
    Reynolds stress is made of.")` *expect:* 10, 8, 2.
17. `nb.check_agree` — **from scratch (curation §7):** a loop over members accumulating Σu and Σu², then mean and variance; the
    product split by hand; `assert np.allclose` against `TS.ensemble_average`, `TS.central_moment(ens, 2)`,
    `TS.product_average_split`.
18. `nb.primer("standard error of a mean", …)` (**P281**): "The mean of N independent samples scatters about the true mean
    with standard deviation σ/√N. 'Within 5 standard errors' is our test for every sampled number in this chapter." code:
    means of 2000 batches of 64 Gaussian samples: std ≈ 0.125 = 1/√64.
19. `nb.note` — **N20 [B]** + `nb.animation` **A1** (`player="frames"`): members of the decaying ensemble appear one by one
    (grey), the running ensemble mean (purple) settles on $e^{-t/10}$ (ghost), a sliding time average of member 1 (teal)
    beside it; stops at N = 2, 4, 8, 64; 40 frames (FAST 24). Then a log–log panel: `TS.standard_error` vs N with the −½
    line. *expect:* fitted slope −0.50 ± 0.05; σ/√64 = 0.3/8 = 0.0375. *see / read / change* ("…members were correlated
    with each other: the scatter would fall more slowly than N^{−1/2}").
20. `nb.note` — **N21 [B]** $\overline{(u-\langle u\rangle)^m}\equiv\frac1N\sum_n\big(u(\mathbf x,t{:}n)-\langle u(\mathbf x,t)\rangle\big)^m$ (12.11); **N22 [B]** variance,
    skewness, kurtosis, standard deviation, rms as a two-row table "book (raw central moments, units u², u³, u⁴) vs usual
    (divided by σ³, σ⁴; Gaussian 0 and 3)"; `TS.statistics(x, normalized=False)` and `(…, normalized=True)`; the identities
    $\overline{(u-\bar u)^2}=\overline{u^2}-\bar u^2$. *expect:* for 10⁵ Gaussian samples skewness 0 ± 0.02, kurtosis 3 ± 0.05.
21. `nb.primer("the sinc function and np.sinc", …)` (**P285**): "sin(x)/x is 1 at x = 0 and zero at x = π, 2π, …; numpy's
    `np.sinc(x)` is sin(πx)/(πx), so sin(a)/a = `np.sinc(a/np.pi)`." code: 3 lines.
22. `nb.note` — **N23 [B]** (Example 12.1, our numbers): for $u=Ae^{-t/\tau_d}+B\cos\omega t$ the window average is
    $\overline{u(t)}=\Big[\frac{\sinh(\Delta t/2\tau_d)}{\Delta t/2\tau_d}\Big]Ae^{-t/\tau_d}+\Big[\frac{\sin(\omega\Delta t/2)}{\omega\Delta t/2}\Big]B\cos\omega t$; a good window
    has $1\ll\omega\Delta t\ll\omega\tau_d$. + `nb.code`: sympy does the window integral (6 commented lines) and
    `ch12.time_average_exp_cos(5.0, 1.0, 0.5, 10.0, 2*np.pi, 0.4)`. *expect:* factors 1.0000667 and 0.7568; average 0.98498.
    `nb.worked_example("a window of 0.4 periods", "A = 1, B = 0.5, τ_d = 10 s, period 1 s (ω = 2π rad/s), Δt = 0.4 s. 1. Mean
    factor: x = Δt/2τ_d = 0.02, sinh(x)/x ≈ 1 + x²/6 = 1.00007 — the mean passes untouched. 2. Wave factor: x = ωΔt/2 = 0.4π =
    1.2566, sin(x)/x = 0.9511/1.2566 = 0.757 — 76 % of the wave leaks into the 'mean'. 3. With Δt = 1 s (one whole period)
    the wave factor is sin(π)/π = 0: the wave is removed exactly. 4. With Δt = 20 s the mean factor is sinh(1)/1 = 1.175: now the
    window is too long and distorts the decaying mean by 17.5 %.")`
23. `nb.plotly` — **F1** `slider_figure` over Δt ∈ 30 values 0.1…20 s (log-spaced): traces "signal", "sliding average
    (`TS.time_average`)", "true mean $Ae^{-t/\tau_d}$", "formula (`time_average_exp_cos`)"; title carries both factors. *see / read /
    change*.
24. `nb.explainer("reynolds_averaging_window", heading="How long must you average to get 'the mean'?", why="Dragging the window
    and the number of members shows the two estimates of the mean fail in different ways — scatter ∝ N^{−1/2}, and the two
    window factors — and the product bars show ⟨ũṽ⟩ − ŪV̄ refusing to vanish.", tries=["Set the window to exactly one wave period:
    the teal curve loses the wave completely.", "Halve the window: 76 % of the wave is back.", "Raise N from 8 to 64: the scatter
    of the purple curve falls by √8 ≈ 2.8.", "Open the product bars and set the correlation to zero: the leftover vanishes."])`
25. `nb.md` — **What would change if…** "…we multiplied two signals before averaging — as the term $\tilde u_j\tilde u_i$ in the
    momentum equation does? The leftover covariance appears, and C04 shows it acts on the mean flow like a stress. First we need
    to measure how a signal is correlated with itself (C02)."

---

### A.4 §12.4 Correlations and Spectra — C02, C03
#### C02 — The autocorrelation in lag form (12.17) and the integral time scale (12.18)
1. `nb.section("12.4", "Correlations and Spectra", intro="**What is this section about?** Two ways to describe how a
   fluctuating signal is organised: how long it remembers itself (the correlation) and which frequencies carry its variance
   (the spectrum). They are the same information in two forms.")`
2. `nb.core("C02", "The autocorrelation $R_{11}(\\tau)=\\overline{u_1(t)u_1(t+\\tau)}=R_{11}(-\\tau)$ (12.17) and the integral time scale
   $\\Lambda_t=\\int_0^\\infty r_{11}(\\tau)\\,d\\tau$ (12.18)", question="How long does a turbulent signal remember itself?")`
3. `nb.md` — **Plain words:** "If the wind is strong now, it is probably still strong half a second from now, and tells you
   nothing about the wind ten minutes from now. Somewhere in between the signal 'forgets'. That forgetting time decides how
   long a record you need (C01), how wide an eddy is, and — in C16 — how smoke spreads."
4. `nb.md` — **The idea:** "Slide a copy of the signal along itself by a lag τ, multiply point by point, average. Aligned
   (τ = 0) every product is a square: the average is the variance. Far apart the signs are unrelated: the average is 0. The
   area under the normalised curve is the memory time Λ_t." (ASCII: two wiggles shifted by τ, product shaded.)
5. `nb.note` — **N24 [B]** $R_{ij}(\mathbf x_1,t_1,\mathbf x_2,t_2)\equiv\overline{u_i(\mathbf x_1,t_1)u_j(\mathbf x_2,t_2)}$ (12.12). ⚠️ "$R_{ij}$ was the rotation tensor in
   Ch. 2–3; here it is a correlation." `TS.correlation(a, b)`. **N26 [B]** the autocorrelation is the i = j = 1 case,
   $R_{11}\equiv\overline{u_1(\mathbf x_1,t_1)u_1(\mathbf x_2,t_2)}$ (12.13).
6. `nb.primer("covariance, the covariance matrix and its ellipse", …)` (**P287**): "The covariance of two zero-mean variables
   is the mean of their product; the 2 × 2 covariance matrix has the variances on the diagonal; its eigenvectors (Ch. 2's
   principal axes) point along the scatter cloud. `TS.correlated_pair(n, r)` draws such a pair." code: 4 lines, r = −0.6,
   printing `np.cov`.
7. `nb.note` — **N25 [B]** four-row table with scatter thumbnails (`TS.correlated_pair`, r = 0, 0.3, 0.9, −0.9): uncorrelated,
   weakly, strongly, anticorrelated. **N27 [B]** $r_{12}\equiv R_{12}/\sqrt{R_{11}R_{22}}=\overline{u_1u_2}\big/\big(\sqrt{\overline{u_1^2}}\sqrt{\overline{u_2^2}}\big)$ (12.14) (arguments as in $R_{ij}(\mathbf x_1,t_1,\mathbf x_2,t_2)\equiv\overline{u_i(\mathbf x_1,t_1)u_j(\mathbf x_2,t_2)}$ (12.12)); `TS.correlation_coefficient` vs `np.corrcoef`. **N28 [B]** $r_{11}\equiv R_{11}(\mathbf x_1,t_1,\mathbf x_2,t_2)/\sqrt{R_{11}(1,1)R_{11}(2,2)}$ (12.15),
   equal to 1 at zero separation.
8. `nb.primer("a quadratic that is never negative has discriminant ≤ 0", …)` (**P288**): "If $a\lambda^2+b\lambda+c\ge0$ for every
   real λ (a > 0), the parabola never dips below the axis, so it has at most one real root: $b^2-4ac\le0$." code: 3 lines
   checking a = 1, b = 2, c = 2 (disc −4) vs c = 0.5 (disc +2, dips below).
9. `nb.derivation("D02", …)` — Part F D02 (6 steps), ref "12.16": **N29 [B]**
   $\lvert\overline{uv}\rvert\le\sqrt{\overline{u^2}}\sqrt{\overline{v^2}}$ (12.16). ⚠️ trap: the page prints no absolute value.
10. `nb.derivation("D03", …)` — Part F D03 (5 steps), ref "12.17".
11. `nb.primer("correlation by FFT with zero padding", …)` (**P289**): "The direct sum over lags costs N × lags; the FFT gives
    all lags at once (transform, multiply by the conjugate, transform back), provided the record is zero-padded so that it
    does not wrap around. 'Unbiased' divides lag m by N − m, the number of products actually available; far lags are noisy."
    code: 5 lines comparing `np.correlate` and the FFT route on 16 numbers.
12. `nb.note` — **N31 [B]** $\Lambda_t\equiv\int_0^\infty r_{11}(\tau)\,d\tau=\frac1{R_{11}(0)}\int_0^\infty R_{11}(\tau)\,d\tau$ (12.18): the width of the rectangle of
    height 1 with the same area. `TS.integral_scale(lag, r)` to the first zero. **N32 [B]** correlation time $t_c$ = first
    zero of $r_{11}$; a record of length Δt holds about $N\approx\Delta t/t_c$ independent samples (`TS.effective_samples`) — "the
    N of C01's N^{−1/2}". **N34 [B]** is the figure of row 16.
13. `nb.worked_example("memory of an exponential correlation", "Take $r(\\tau)=e^{-\\tau/\\tau_c}$ with τ_c = 0.5 s. 1.
    $\\Lambda_t=\\int_0^\\infty e^{-\\tau/0.5}d\\tau=0.5$ s. 2. A 100 s record holds about 100/0.5 = 200 memory times. 3. So the standard error of
    its time mean is about σ/√(200/2) = 0.1σ (two memory times per independent sample for this shape). 4. For a Gaussian
    shape $r=e^{-\\tau^2/\\tau_c^2}$: $\\Lambda_t=\\tfrac{\\sqrt\\pi}2\\tau_c=0.443$ s.")` *expect:* 0.5; 200; 0.443.
14. `nb.code` — `u = TS.make_ensemble(1, t, lambda t: 0*t, 1.0, 0.5, seed=3)[0]` (T = 2000 s, dt = 0.01); `lag, R =
    TS.autocorrelation(u, dt, max_lag=500)`; `r = R/R[0]`; `TS.integral_scale(lag, r)`, `TS.correlation_time(lag, r)`,
    `TS.effective_samples(t[-1], …)`. *expect:* Λ_t = 0.50 ± 0.05 s (seeded). *explain:* 3 numbered points.
15. `nb.check_agree` — **from scratch:** the direct-sum autocorrelation (loop over 200 lags) and Λ_t by `np.trapezoid` to the
    first zero; `assert np.allclose(R_loop, R[:200])`, `assert np.isclose(Lam_mine, TS.integral_scale(lag, r))`.
16. `nb.figure` — **N34** (our Fig. 12.5): r(τ) teal, the exact $e^{-\tau/\tau_c}$ muted dashed, the equal-area rectangle (height
    1, width Λ_t) purple, the first zero marked. *see / read* ("the rectangle's width is the memory") */ change* ("…the record
    were 20 s instead of 2000 s: the tail would wander about zero and Λ_t would depend on where we stop integrating").
17. `nb.note` — **N30 [B]** + `nb.animation` **A2** (video, 60 frames, FAST 30): a copy of the signal delayed by t_o = 0.8 s
    slides past the original; the product is shaded; `TS.cross_correlation` is drawn point by point and peaks at τ = t_o.
    "The cross-correlation is not even: $R_{uv}(\tau)=R_{vu}(-\tau)$ (D03 step 5)." *see / read / change*.
18. `nb.derivation("D04", …)` — Part F D04 (7 steps), ref "12.19": **N33 [B]** $\lambda_t^2\equiv-2\big/\big[d^2r_{11}/d\tau^2\big]_{\tau=0}$ (12.19). +
    `nb.code`: `u_s = TS.smooth_signal(2**16, 0.01, spectrum=lambda w: np.exp(-(w*0.5)**2/4), seed=4)`; `TS.taylor_microscale(lag, r)`
    *expect:* λ_t ≈ 0.50 s (Gaussian pair, λ_t = τ_c); the osculating parabola added to the figure of row 16 (rose).
    ⚠️ "An Ornstein–Uhlenbeck signal has a cusp, $r\approx1-\lvert\tau\rvert/\tau_c$: no parabola fits, its microscale does not exist."
19. `nb.md` — **What would change if…** "…we asked which *frequencies* carry the variance instead of how long the memory is?
    That is the Fourier transform of this same curve (C03)."

#### C03 — The energy spectrum (12.20)–(12.22)
1. `nb.core("C03", "The energy spectrum $S_e(\\omega)=\\frac1{2\\pi}\\int_{-\\infty}^{+\\infty}R_{11}(\\tau)e^{-i\\omega\\tau}d\\tau$ (12.20) and the variance it
   distributes, $\\overline{u_1^2}=\\int_{-\\infty}^{+\\infty}S_e(\\omega)\\,d\\omega$ (12.22)", question="Which frequencies carry the variance, and why is that
   the same information as the memory?")`
2. `nb.md` — **Plain words:** "A graphic equaliser shows how loud each pitch is. A spectrum does the same for a wind record: how
   much of the variance comes from slow gusts, how much from fast flutter. Ocean and atmosphere records are almost always
   shown this way, and the −5/3 law of C08 is a statement about a spectrum."
3. `nb.md` — **The idea:** table | long memory | slow wiggles | spectrum squeezed near ω = 0 | vs | short memory | fast wiggles
   | spectrum flat out to ω ≈ 1/Λ_t |; "area under S = variance, always; height at zero = variance × Λ_t/π."
4. `nb.primer("Fourier-transform pair: where the 2π sits", …)` (**P290**): "A transform pair needs one factor 1/2π in total.
   The book puts it in the forward transform and uses angular frequency ω [rad/s]. For a real, even function the complex
   exponential reduces to a cosine: $\frac1{2\pi}\int R\,e^{-i\omega\tau}d\tau=\frac1\pi\int_0^\infty R\cos\omega\tau\,d\tau$." code: 4 lines doing the cosine
   integral of $e^{-\lvert\tau\rvert/0.5}$ at ω = 2 with `quad`: 0.0796 = τ_c/[π(1 + ω²τ_c²)].
5. `nb.derivation("D05", …)` — Part F D05 (9 steps), ref "12.20". Carries **N35 [B]** $R_{11}(\tau)\equiv\int_{-\infty}^{+\infty}S_e(\omega)e^{+i\omega\tau}d\omega$ (12.21),
   **N36 [B]** $\overline{u_1^2}\equiv\int_{-\infty}^{+\infty}S_e(\omega)\,d\omega$ (12.22), **N37 [B]** $S_e(0)=\overline{u_1^2}\Lambda_t/\pi$.
6. `nb.worked_example("the pair for an exponential memory", "σ² = 1 m²/s², τ_c = 0.5 s. 1. $S_e(\\omega)=\\sigma^2\\tau_c/[\\pi(1+\\omega^2\\tau_c^2)]$.
   2. At ω = 0: 0.5/π = 0.159 m²/s — and $\\overline{u^2}\\Lambda_t/\\pi$ = 1 × 0.5/π ✓. 3. At ω = 2 rad/s (= 1/τ_c): 0.5/(2π) = 0.0796, half the
   peak: the spectrum is flat to about 1/τ_c. 4. Area: $\\int_{-\\infty}^\\infty S_e\\,d\\omega=\\frac{\\sigma^2}\\pi[\\arctan]_{-\\infty}^{\\infty}=\\sigma^2$ ✓.")`
7. `nb.code` — `pair = TS.correlation_spectrum_pair("exponential", 1.0, 0.5)`; `S = TS.spectrum_from_correlation(lag, R, omega)`;
   `TS.spectrum_variance(omega, S)`; `TS.integral_scale_from_spectrum(S[0], R[0])`; round trip
   `TS.correlation_from_spectrum(omega, S, lag)`. *expect:* S(0) ≈ 0.159, variance ≈ 1.0 (± 3 %), Λ_t ≈ 0.5.
8. `nb.primer("the discrete Fourier transform as a Riemann sum; Parseval", …)` (**P291**): "`np.fft.rfft(u)*dt` approximates
   $\int u\,e^{-i\omega t}dt$ at $\omega_k=2\pi k/(N\,dt)$; the two-sided density in the book's normalisation is $\lvert\hat u\rvert^2/(2\pi T)$; Parseval:
   its sum times Δω (both signs of ω) is the variance. One-sided = 2 × two-sided." code: 5 lines; printed ratio 1.000.
9. `nb.primer("segment averaging (Welch) and leakage", …)` (**P292**): "One raw periodogram is as noisy as its own mean
   however long the record; cutting the record into K segments and averaging their periodograms reduces the scatter by
   √K at the price of resolution. `scipy.signal.welch` does this (per Hz, one-sided: convert by 1/(4π))." code: 3 lines.
10. `nb.note` — **N38 [B]** "(Exercise 12.8, stated) the same spectrum from a finite record:
    $S_e(\omega)=\lim_{T\to\infty}\frac1{2\pi T}\big\lvert\int_{-T/2}^{+T/2}u(t)e^{-i\omega t}dt\big\rvert^2$ — the name of this equivalence is the Wiener–Khinchin theorem; the
    proof expands the square as a double integral whose triangular weight $(1-\lvert\tau\rvert/T)\to1$." + `nb.figure` **three routes**:
    (1) correlation → `spectrum_from_correlation` (bold teal), (2) `TS.periodogram(u, dt, segments=32)` (dots), (3)
    `signal.welch` converted (thin orange), exact $S_e$ muted; log–log. *see / read / change*.
11. `nb.check_agree` — **from scratch:** the cosine transform of R(τ) by `np.trapezoid` at 50 frequencies; a periodogram from
    `np.fft.rfft` scaled by hand; `assert np.allclose(S_mine, S, rtol=1e-10)`; `assert np.isclose(np.sum(P_mine)*dw*2, u.var(),
    rtol=1e-12)` (Parseval).
12. `nb.note` — **N39 [B]** homogeneous field: $R_{ij}(\mathbf r)\equiv\overline{u_i(\mathbf x)u_j(\mathbf x+\mathbf r)}$ (12.23); `TS.spatial_correlation` on the N02 field.
13. `nb.primer("change of variable in a density; units of a spectral density", …)` (**P293**): "A density is 'variance per
    unit ω'. Changing the axis to $k=\omega/U_0$ must keep the area: $S(k)\,dk=S(\omega)\,d\omega$ ⇒ $S(k)=U_0S(\omega)$. Units: $S_e$ [m²/s² per
    rad/s] = m²/s; $S_{11}$ [m²/s² per rad/m] = m³/s²." code: 3 lines checking the two areas agree.
14. `nb.note` — **N40 [B]** Taylor's frozen-turbulence hypothesis: a probe moving at $U_0$ (or a mean wind sweeping the eddies
    past a mast) turns time into distance, $x=U_0t$, $\partial u_1/\partial x_1\approx-(1/U_1)\,\partial u_1/\partial t$, so $k_1=\omega/U_0$ and
    $S_{11}(k_1)=U_0S_e(\omega)$. `TS.taylor_frozen`, `TS.frequency_to_wavenumber_spectrum`, `ch12.frozen_field_probe` (error table for
    u_rms/U₀ = 0.05, 0.2, 0.5). **Climate hook:** tower and aircraft spectra are all converted this way.
15. `nb.plotly` — **F2** `slider_figure` over τ_c (25 values 0.1…5 s): left r(τ) with the Λ_t rectangle, right $S_e(\omega)$ log–log
    with S(0) marked (`correlation_spectrum_pair`). *see / read / change*.
16. `nb.explainer("correlation_and_spectrum", heading="Why are a correlation and a spectrum the same information?", why="Stretching
    the memory narrows the spectrum in front of you while its area stays the variance and its height at zero follows Λ_t;
    sliding the lag shows the product being averaged.", tries=["Double the memory time: S(0) doubles, the corner frequency
    halves, the shaded area does not move.", "Switch to 'damped cosine' and raise ω₀: a peak leaves zero frequency and Λ_t
    collapses — an oscillating correlation has little net area.", "Click the spectrum: the inspector shows the cosine-transform
    arithmetic at that ω.", "Switch the axis to wavenumber with U₀ = 10 m/s."])`
17. `nb.md` — **What would change if…** "…we now put the decomposition mean + fluctuation into the equations of motion and
    averaged? Every term passes through the average except the product — C04."

---

### A.5 §12.5 Averaged Equations of Motion — R01 R02, C04
#### C04 — The Reynolds-averaged momentum equation (12.30) and the Reynolds stress
1. `nb.section("12.5", "Averaged Equations of Motion", intro="**What is this section about?** We average the equations of
   motion. The mean flow obeys almost the same equations as before — plus one new term made of fluctuations, which acts like
   a stress and for which the averaging gives no equation.")`
2. `nb.core("C04", "The Reynolds-averaged momentum equation $\\frac{\\partial U_i}{\\partial t}+U_j\\frac{\\partial U_i}{\\partial x_j}=-g[1-\\alpha(\\bar T-T_0)]\\delta_{i3}+\\frac1{\\rho_0}\\frac{\\partial\\bar\\tau_{ij}}{\\partial x_j}$,
   $\\bar\\tau_{ij}=-P\\delta_{ij}+2\\mu\\bar S_{ij}-\\rho_0\\overline{u_iu_j}$ (12.30)", question="How can fluctuations that average to zero push on the mean
   flow?")`
3. `nb.md` — **Plain words:** "An aircraft designer wants the mean drag, a climate modeller the mean wind; neither wants every
   gust (**N41 [C]**: the flight lasts hours, the gusts milliseconds, and the range of scales is too wide to compute — C07
   puts a number on it). So we ask what equation the *mean* obeys. The surprise: the gusts do not average out of it."
4. `nb.md` — **The idea** (ASCII): a shear U(y) with two parcels swapping levels:
   ```
   fast  ────────────►   y + ℓ     a parcel coming DOWN  (v < 0) arrives too fast   (u > 0)   → uv < 0
   slow  ─────►          y         a parcel going  UP    (v > 0) arrives too slow   (u < 0)   → uv < 0
   every exchange carries x-momentum downward: mean(uv) < 0  ⇒  a stress −ρ₀·mean(uv) > 0 on the mean flow
   ```
5. `nb.recap("R01", "The Boussinesq equations in flux form", "From Ch. 4: continuity $\\partial\\tilde u_i/\\partial x_i=0$ (4.10), momentum $\\frac{\\partial\\tilde u_i}{\\partial t}+\\tilde u_j\\frac{\\partial\\tilde u_i}{\\partial x_j}=-\\frac1{\\rho_0}\\frac{\\partial\\tilde p}{\\partial x_i}-g[1-\\alpha(\\tilde T-T_0)]\\delta_{i3}+\\nu\\frac{\\partial^2\\tilde u_i}{\\partial x_j^2}$ (4.86) and heat $\\frac{\\partial\\tilde T}{\\partial t}+\\tilde u_j\\frac{\\partial\\tilde T}{\\partial x_j}=\\kappa\\frac{\\partial^2\\tilde T}{\\partial x_j^2}$ (4.89). Adding $\\tilde u_i\\,\\partial\\tilde u_j/\\partial x_j=0$ turns $\\tilde u_j\\,\\partial\\tilde u_i/\\partial x_j$ into $\\partial(\\tilde u_j\\tilde u_i)/\\partial x_j$: $\\frac{\\partial\\tilde u_i}{\\partial t}+\\frac{\\partial}{\\partial x_j}(\\tilde u_j\\tilde u_i)=-\\frac1{\\rho_0}\\frac{\\partial\\tilde p}{\\partial x_i}-g[1-\\alpha(\\tilde T-T_0)]\\delta_{i3}+\\nu\\frac{\\partial^2\\tilde u_i}{\\partial x_j^2}$.
   This is step 1 of D06.", where="Ch. 4 §4.9, C13 (D27–D28)")`
6. `nb.note` — **N42 [B]** $\tilde u_i=U_i+u_i$, $\tilde p=P+p$, $\tilde\rho=\bar\rho+\rho'$, $\tilde T=\bar T+T'$ (12.24); **N43 [B]**
   $\overline{\tilde u_i}=U_i$, $\overline{\tilde p}=P$, $\overline{\tilde\rho}=\bar\rho$, $\overline{\tilde T}=\bar T$ (12.25); **N44 [B]** $\overline{u_i}=0$, $\bar p=0$, $\overline{\rho'}=0$, $\overline{T'}=0$ (12.26).
   `mean, fluct = TS.reynolds_decompose(ens)`; `np.abs(fluct.mean(axis=0)).max()` *expect:* ≲ 1e-15.
7. `nb.primer("a sympy averaging operator", …)` (**P294**): "We teach sympy the rules of D01: the average is linear, leaves a
   mean alone, kills a single fluctuation, and keeps a product of fluctuations as a new symbol." code: 5 lines expanding
   (U + u)(V + v) and averaging to `U*V + uv_bar`.
8. `nb.derivation("D06", …)` — Part F D06 (13 steps), ref "12.30". Carries **N45 [B]** $\partial U_i/\partial x_i=0$ (12.27), **N46 [B]**
   $\partial u_i/\partial x_i=0$ (12.28), **N47 [B]** $\frac{\partial(U_i+u_i)}{\partial t}+\frac{\partial}{\partial x_j}\big((U_j+u_j)(U_i+u_i)\big)=-\frac1{\rho_0}\frac{\partial(P+p)}{\partial x_i}-g[1-\alpha(\bar T+T'-T_0)]\delta_{i3}+\nu\frac{\partial^2(U_i+u_i)}{\partial x_j^2}$ (12.29), **N48 [B]** the term-by-term averages.
9. `nb.code` — `eqs = ch12.rans_sympy()` printing the steps (cached, P252); *expect:* the collected equation has exactly one
   term with no counterpart in the instantaneous equation, $\partial\overline{u_iu_j}/\partial x_j$.
10. `nb.primer("counting the independent components of a symmetric tensor", …)` (**P295**): "A symmetric 3 × 3 tensor has 6
    independent components (3 diagonal + 3 above it); a fully symmetric triple product $\overline{u_iu_ju_k}$ has 10." code:
    `len({tuple(sorted(c)) for c in itertools.product(range(3), repeat=3)})` → 10.
11. `nb.note` — **N49 [B]** "**The Reynolds stress tensor** $-\rho_0\overline{u_iu_j}$: symmetric, six components; the diagonal ones
    are normal stresses that add to the mean pressure, the off-diagonal ones shear stresses. It is a covariance matrix
    (P287). Except within a hair of a wall it is far larger than the viscous stress." `TS.reynolds_stress(samples, rho0=1.2)`,
    `TS.anisotropy_tensor`.
12. `nb.recap("R02", "Sign convention of a stress on a face", "As for the stress cube of Ch. 2 (Fig. 2.4): a positive
    $-\\rho_0\\overline{uv}$ on the face whose outward normal is +y points along +x.", where="Ch. 2 §2.4")` — drawn inside the figure of row 15
    (`stress_element_sketch`).
13. `nb.note` — **N50 [B]** the parcel argument: a parcel displaced by ℓ keeps its old mean speed for a while, so where it
    arrives $u\approx-\ell\,dU/dy$; multiplying by v and averaging, $\overline{uv}\approx-\overline{v\ell}\,dU/dy<0$ when $dU/dy>0$ (v and ℓ have the same
    sign). **N51 [B]** "Reynolds stress = momentum flux: $\rho_0\overline{(U+u)v}=\rho_0U\bar v+\rho_0\overline{uv}=\rho_0\overline{uv}$ is the mean flux of x-momentum
    in the y-direction carried by the fluctuations."
14. `nb.worked_example("the stress carried by parcels", "dU/dy = 2 s⁻¹, parcels move ℓ_rms = 0.1 m with v_rms = 0.5 m/s and
    keep all their momentum. 1. u ≈ −ℓ dU/dy has rms 0.2 m/s. 2. $\\overline{uv}=-v_{rms}\\ell_{rms}\\,dU/dy$ = −0.5 × 0.1 × 2 = −0.1 m²/s². 3. In
    air (ρ₀ = 1.2 kg/m³): Reynolds shear stress $-\\rho_0\\overline{uv}$ = 0.12 Pa. 4. Viscous stress μ dU/dy = 1.8 × 10⁻⁵ × 2 = 3.6 × 10⁻⁵ Pa —
    3000 times smaller. 5. The viscosity that would do the same job: $\\nu_T=-\\overline{uv}/(dU/dy)$ = 0.05 m²/s, against ν = 1.5 × 10⁻⁵.")`
15. `nb.code` + `nb.figure` — `out = ch12.displaced_parcel_uv(2.0, 0.1, 0.5, n=100000, seed=0)`; `ch12.parcel_uv_expected(2.0,
    0.1, 0.5)`; assert within 5 standard errors (P281). *expect:* uv = −0.100 ± 0.001; r_uv ≈ −1 (correlation 1) and ≈ −0.5
    with `correlation=0.5`. Figure (our Figs. 12.6–12.8): left `parcel_sketch` + stress element; right the (u, v) scatter
    with its covariance ellipse, quadrants 2 and 4 shaded orange. *see / read* ("the tilt of the cloud *is* the stress") */
    change* ("…dU/dy < 0: the cloud tilts the other way, $\overline{uv}>0$; …no shear: a round cloud, no stress — the isotropic case of C05").
16. `nb.animation` — **A3** (video, 80 frames, FAST 40): parcels exchanged across the shear, each dropping a dot in the (u, v)
    plane; the running $\overline{uv}$ converges to −0.1 with its ±5 standard-error band. *see / read / change*.
17. `nb.check_agree` — **from scratch:** the 2 × 2 covariance with explicit sums; `assert np.allclose(mine, TS.velocity_covariance(…))`, and `-1.2*mine` against `TS.reynolds_stress(…, rho0=1.2)`;
    `assert abs(out["uv"] - (-np.mean(v*l)*2.0)) < 1e-12`.
18. `nb.note` — heat and scalar (each 2–4 lines with its equation): **N52 [B]**
    $\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}+\frac{\partial}{\partial x_j}(\overline{u_jT'})=\kappa\frac{\partial^2\bar T}{\partial x_j^2}$ (12.31) (⚠️ this κ is the thermal diffusivity, `kappa_th`); **N53 [B]**
    $\rho_0C_p\big(\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}\big)=-\frac{\partial Q_j}{\partial x_j}$, $Q_j=-k\frac{\partial\bar T}{\partial x_j}+\rho_0C_p\overline{u_jT'}$ (12.32) (`mean_heat_flux`; Fourier's law
    $\mathbf q=-k\nabla T$ (1.2) recalled); **N54 [B]** turbulent heat flux $\rho_0C_p\overline{u_jT'}$, upward over a heated surface.
19. `nb.primer("kinematic vs dynamic fluxes", …)` (**P296**): "Observers report H [W/m²] and τ₀ [Pa]; the equations carry
    $\overline{wT'}$ [K m/s] and $u_*^2$ [m²/s²]. $H=\rho c_p\overline{wT'}$, $\tau_0=\rho u_*^2$; for a perfect gas α = 1/T." code: 3 lines: wT = 0.1 K m/s in air ⇒
    H = 1.2 × 1005 × 0.1 = 120.6 W/m².
20. `nb.note` — **N55 [C]** passive scalar and mixture density $\rho_m=\tilde\upsilon\rho_s+(1-\tilde\upsilon)\rho$ (`mixture_density`; used for jet dilution in C09);
    **N56 [B]** $\frac{\partial}{\partial t}(\rho_m\tilde Y)+\frac{\partial}{\partial x_j}(\rho_m\tilde Y\tilde u_j)=\frac{\partial}{\partial x_j}\big(\rho_m\kappa_m\frac{\partial\tilde Y}{\partial x_j}\big)$ (12.33); **N57 [B]**
    $\frac{\partial\bar Y}{\partial t}+U_j\frac{\partial\bar Y}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\kappa_m\frac{\partial\bar Y}{\partial x_j}-\overline{u_jY'}\big)$ (12.34) — "the D06 moves a third time, at constant $\rho_m$".
21. `nb.note` — **N58 [B]** "**RANS** = the mean continuity $\partial U_i/\partial x_i=0$ (12.27), the mean momentum equation $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-g[1-\alpha(\bar T-T_0)]\delta_{i3}+\frac1{\rho_0}\frac{\partial\bar\tau_{ij}}{\partial x_j}$ (12.30), the mean heat equation $\rho_0C_p\big(\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}\big)=-\frac{\partial Q_j}{\partial x_j}$ (12.32) and the mean scalar equation $\frac{\partial\bar Y}{\partial t}+U_j\frac{\partial\bar Y}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\kappa_m\frac{\partial\bar Y}{\partial x_j}-\overline{u_jY'}\big)$ (12.34). ⚠️ From §12.8 to
    §12.10 the book silently drops gravity and the subscript 0: $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-\frac1\rho\frac{\partial P}{\partial x_i}+\frac{\partial}{\partial x_j}\big(2\nu\bar S_{ij}-\overline{u_iu_j}\big)$, P = deviation from
    hydrostatic."
22. `nb.note` — **N59 [B]** the transport equation for the Reynolds stress, stated term by term in a table (term · name ·
    colour):
    $\frac{\partial\overline{u_iu_j}}{\partial t}+U_k\frac{\partial\overline{u_iu_j}}{\partial x_k}+\frac{\partial\overline{u_iu_ju_k}}{\partial x_k}=-\overline{u_iu_k}\frac{\partial U_j}{\partial x_k}-\overline{u_ju_k}\frac{\partial U_i}{\partial x_k}-\frac1\rho\Big(\overline{u_i\frac{\partial p}{\partial x_j}}+\overline{u_j\frac{\partial p}{\partial x_i}}\Big)-2\nu\overline{\frac{\partial u_i}{\partial x_k}\frac{\partial u_j}{\partial x_k}}+\nu\frac{\partial^2\overline{u_iu_j}}{\partial x_k^2}+g\alpha\big(\overline{u_jT'}\delta_{i3}+\overline{u_iT'}\delta_{j3}\big)$ (12.35).
    `ch12.reynolds_stress_budget_sympy()` *expect:* residual 0 for the 12-component; half the trace reproduces the turbulent-energy budget of C06, $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47) (D10 derives that half-trace in full).
23. `nb.note` — **N60 [B]** "**The closure problem.** `ch12.closure_count(3)`: the mean equations are 4 (continuity + 3
    momentum) for 4 + 6 = 10 unknowns (U_i, P, six $\overline{u_iu_j}$); the Reynolds-stress equation of N59 above gives 6 more equations but brings 10 triple
    correlations and the pressure and dissipation terms. Averaging always loses. Three responses: model the unknown
    correlations (RANS models: C12, C13), compute every scale (DNS: its cost is in C07), or compute the large eddies and
    model the small (LES, Ch. 10 pointer)." **N61 [C]** "With Reynolds stresses the mean momentum equation has no Bernoulli
    integral, even for steady 'inviscid' mean flow (Ch. 4 `which_bernoulli`)."
24. `nb.explainer("reynolds_stress_parcels", heading="How can fluctuations that average to zero push the mean flow?", why="You
    release parcels in a shear and watch each land in the (u, v) scatter; turning the shear down, to zero and negative tilts
    the cloud and the stress bar follows — the sign argument becomes a count.", tries=["Set the shear to zero: the cloud is
    round and the orange bar vanishes.", "Reverse the shear: the cloud tilts the other way.", "Lower 'keeps its momentum' to
    0.3: the cloud fattens, r_uv and the stress fall together.", "Click a parcel: its ℓ, u, v and uv."])`
25. `nb.md` — **What would change if…** "…the turbulence had no preferred direction at all? Then the cloud is round, the shear
    stresses vanish, and the whole two-point tensor collapses to one function — C05."

---

### A.6 §12.6 Homogeneous Isotropic Turbulence — C05
#### C05 — Isotropic correlations (12.38)–(12.41) and the isotropic dissipation (12.43)
1. `nb.section("12.6", "Homogeneous Isotropic Turbulence", intro="**What is this section about?** The simplest turbulence:
   statistically the same at every point and in every direction. Symmetry then does most of the work — nine correlation
   functions become one, and the dissipation rate becomes one measurable gradient.")`
2. `nb.core("C05", "Isotropic dissipation $\\bar\\varepsilon=-15\\nu\\overline{u^2}[d^2f/dr^2]_{r=0}=30\\nu\\overline{u^2}/\\lambda_f^2=15\\nu\\overline{u^2}/\\lambda_g^2$ (12.43) from the correlation coefficients f and g", question="How much can symmetry alone tell us about turbulence?")`
3. `nb.md` — **Plain words:** "Behind a grid in a wind tunnel the turbulence has no mean shear to remember and soon looks
   the same in every direction (**N62 [B]**: homogeneous = no preferred place; isotropic = no preferred direction; grid
   turbulence is the laboratory approximation). The small eddies of every flow look like this too (local isotropy), which
   is why results of this section are used far beyond wind tunnels — for instance to measure the dissipation rate in the
   ocean from one velocity gradient."
4. `nb.md` — **The idea:** "Two velocity vectors a distance r apart. Only two arrangements are different: both components
   *along* the line joining the points (f), or both *across* it (g). Incompressibility ties g to f. So one curve f(r) holds
   everything." + `fg_geometry_sketch` (**N71 [B]**, our Fig. 12.9).
5. `nb.note` — **N63 [B]** scatter of (u, v) (`TS.correlated_pair`): round cloud ⇒ $\overline{uv}=0$; cloud along v = −u ⇒ $\overline{uv}<0$.
   **N64 [B]** $\frac{\partial}{\partial x_i}\overline{u_j^n}=0$, $\overline{u_1^2}=\overline{u_2^2}=\overline{u_3^2}$, $\overline{(\partial u_1/\partial x_1)^n}=\overline{(\partial u_2/\partial x_2)^n}=\overline{(\partial u_3/\partial x_3)^n}$ (12.36); **N65 [B]**
   the six cross-gradient moments are equal, $\overline{(\partial u_1/\partial x_2)^n}=\overline{(\partial u_1/\partial x_3)^n}=\dots=\overline{(\partial u_3/\partial x_2)^n}$ (12.37). `ch12.isotropy_report(u, v, w, dx)` on a 3-D
   synthetic field (48³ FAST / 64³, cached). *expect:* normal-stress ratios 1.00 ± 0.05; shear correlations ≲ 0.03.
6. `nb.note` — **N66 [B]** $f(r)\equiv\overline{u_\parallel(\mathbf x+\mathbf r)u_\parallel(\mathbf x)}\big/\overline{u_\parallel^2}$ and $g(r)\equiv\overline{u_\perp(\mathbf x+\mathbf r)u_\perp(\mathbf x)}\big/\overline{u_\perp^2}$ (12.38), with
   $\overline{u_\parallel^2}=\overline{u_\perp^2}=\overline{u^2}$, f(0) = g(0) = 1; **N67 [B]** $\Lambda_f\equiv\int_0^\infty f\,dr$, $\Lambda_g\equiv\int_0^\infty g\,dr$, $\lambda_f^2\equiv-2/[d^2f/dr^2]_{r=0}$,
   $\lambda_g^2\equiv-2/[d^2g/dr^2]_{r=0}$ (12.39) — "C02's functions with r for τ". **N72 [B]** $R_{ii}(0)=\overline{u_iu_i}=2\bar e$; ⚠️ "$\overline{u^2}$ is **one**
   component: $\bar e=\tfrac32\overline{u^2}$; and $\bar e$ is not Ch. 1's internal energy e."
7. `nb.primer("derivatives of functions of r = ∣r∣", …)` (**P297**): "$\partial r/\partial r_j=r_j/r$ (from $r^2=r_kr_k$), $\partial r_i/\partial r_j=\delta_{ij}$,
   $\delta_{jj}=3$, $r_jr_j=r^2$; so $\partial F(r)/\partial r_j=F'(r)\,r_j/r$." code: 4 sympy lines verifying the first and the divergence of $r_i$ (= 3).
8. `nb.derivation("D07", …)` — Part F D07 (14 steps, ★★★, `check_src`), ref "12.41". Carries **N68 [B]**
   $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40), **N69 [B]** $R_{ij}=\overline{u^2}\big\{f\delta_{ij}+\frac r2\frac{df}{dr}\big(\delta_{ij}-\frac{r_ir_j}{r^2}\big)\big\}$ (12.41) (slip #14 callout), **N70 [B]**
   $g=f+\frac r2f'$, $\Lambda_g=\Lambda_f/2$, $\lambda_g=\lambda_f/\sqrt2$ (⚠️ three-dimensional results).
9. `nb.figure` — f and g from the 3-D synthetic field (`TS.longitudinal_transverse_correlation`) with `transverse_from_longitudinal(r, f)` dashed on top of g. *see* (g dips below zero, f does not) */ read* ("fluid crossing the line must
   come back: the transverse correlation has a negative lobe so that its integral is half of f's") */ change*.
10. `nb.note` — **N73 [B]** $\bar\varepsilon=\frac\nu2\overline{(\partial u_i/\partial x_j+\partial u_j/\partial x_i)^2}$ (12.42) — the average of Ch. 4's $\varepsilon=2\nu S_{ij}S_{ij}$ (4.58) for the
    fluctuations; `ch12.dissipation_rate`.
11. `nb.derivation("D08", …)` — Part F D08 (13 steps, ★★★, `check_src`), ref "12.43". Carries **N74 [B]**
    $\overline{\frac{\partial u_i}{\partial x_k}\frac{\partial u_j}{\partial x_l}}=-\big(\frac{\partial^2R_{ij}}{\partial r_k\partial r_l}\big)_{r=0}$ and the moments 2 : 4 : −1.
12. `nb.worked_example("dissipation from a Taylor microscale", "Air, ν = 1.5 × 10⁻⁵ m²/s; $\\overline{u^2}$ = 1 m²/s² (one component); λ_f = 1 cm.
    1. $\\bar\\varepsilon=30\\nu\\overline{u^2}/\\lambda_f^2$ = 30 × 1.5 × 10⁻⁵ × 1/10⁻⁴ = 4.5 m²/s³. 2. λ_g = λ_f/√2 = 7.07 mm; 15 × 1.5 × 10⁻⁵/(5 × 10⁻⁵) = 4.5 ✓.
    3. One gradient: $\\overline{(\\partial u_1/\\partial x_1)^2}=2\\overline{u^2}/\\lambda_f^2$ = 2 × 10⁴ s⁻²; 15ν × 2 × 10⁴ = 4.5 ✓. 4. $R_\\lambda=\\lambda_g\\sqrt{\\overline{u^2}}/\\nu$ = 471. 5. Kolmogorov
    scale (C07): η = (ν³/ε̄)^{1/4} = 0.165 mm — 43 times smaller than λ_g.")`
13. `nb.code` — `ch12.dissipation_isotropic(1.5e-5, 1.0, lambda_f=0.01)`, `(…, lambda_g=0.01/np.sqrt(2))`, `(…, dudx_sq=2e4)`;
    `ch12.taylor_reynolds_number(1.0, 0.01/np.sqrt(2), 1.5e-5)`; `ch12.gradient_moments_isotropic_sympy()`. *expect:* 4.5 three
    times; 471.4; (2, 4, −1), 30.
14. `nb.check_agree` — **from scratch:** a Gaussian $f=e^{-r^2/\lambda_f^2}$ on a grid: $f''(0)$ by a centred difference, g by
    `np.gradient`, λ_g from g; ε̄ three ways; `assert np.allclose` with `dissipation_isotropic` and `isotropic_scales`
    (Λ_g/Λ_f = 0.5, λ_g/λ_f = 0.7071).
15. `nb.note` — **N75 [B]** $R_\lambda\equiv\lambda_{(g\text{ or }f)}\sqrt{\overline{u^2}}/\nu$ (12.44) — "say which λ: a factor √2". **N76 [B]**
    $S_{11}(k_1)=\frac1{2\pi}\int_{-\infty}^{+\infty}R_{11}(r_1)e^{-ik_1r_1}dr_1=\frac{\overline{u_1^2}}{2\pi}\int_{-\infty}^{+\infty}f(r_1)e^{-ik_1r_1}dr_1$ (12.45): C03's pair with $(r_1,k_1)$ for (τ, ω);
    units m³/s² (P293). **N77 [B]** the measurement recipe as a four-step pipeline cell: `TS.taylor_frozen` →
    `TS.autocorrelation` → `TS.spectrum_from_correlation`, against `TS.periodogram`; a figure of $S_{11}(k_1)$ both ways.
16. `nb.md` — **What would change if…** "…there *is* a mean shear? Then the cloud of C04 tilts, the turbulence can take energy
    from the mean flow — and we need a budget for that energy (C06)."

---

### A.7 §12.7 Turbulent Energy Cascade and Spectrum — C06, C07, C08
#### C06 — The turbulent kinetic-energy budget (12.47) and the mean-flow budget (12.46)
1. `nb.section("12.7", "Turbulent Energy Cascade and Spectrum", intro="**What is this section about?** Energy. The mean flow
   loses it to the large eddies, they hand it to smaller ones, and the smallest turn it into heat. We write the two budgets,
   find the size of the smallest eddies, and derive the most famous formula of the subject, the −5/3 spectrum.")`
2. `nb.core("C06", "The turbulent kinetic-energy budget $\\frac{\\partial\\bar e}{\\partial t}+U_j\\frac{\\partial\\bar e}{\\partial x_j}=\\frac{\\partial}{\\partial x_j}\\big(-\\frac1{\\rho_0}\\overline{pu_j}+2\\nu\\overline{u_iS'_{ij}}-\\frac12\\overline{u_i^2u_j}\\big)-2\\nu\\overline{S'_{ij}S'_{ij}}-\\overline{u_iu_j}\\frac{\\partial U_i}{\\partial x_j}+g\\alpha\\overline{u_3T'}$
   (12.47)", question="Where does the energy of the turbulence come from, and where does it go?")`
3. `nb.md` — **Plain words:** "Switch off the fan and the turbulence in a room dies within seconds: it needs feeding. In a
   river the feed is the mean current rubbing against the bed; over a hot road it is buoyancy. An ocean-mixing or
   boundary-layer scheme in a climate model is, at heart, a bookkeeping of this energy."
4. `nb.md` — **The idea** (ASCII, colours as convention 9):
   ```
   mean flow  Ē = ½U_i²   ──[ shear production  −mean(u_i u_j) ∂U_i/∂x_j ]──►  turbulence  ē = ½ mean(u_i²)  ──[ ε̄ ]──► heat
        │  (tiny direct viscous loss 2ν S̄S̄ ~ 1/Re of production)                 ▲  buoyancy gα·mean(wT′): + heated below, − stable
   the SAME term appears in both budgets with OPPOSITE signs
   ```
5. `nb.derivation("D09", …)` — Part F D09 (9 steps), ref "12.46": **N78 [B]**
   $\frac{\partial\bar E}{\partial t}+U_j\frac{\partial\bar E}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac{U_jP}{\rho_0}+2\nu U_i\bar S_{ij}-\overline{u_iu_j}U_i\big)-2\nu\bar S_{ij}\bar S_{ij}+\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}-\frac g{\rho_0}\bar\rho U_3$ (12.46).
6. `nb.note` — **N79 [B]** "Divergence terms only move energy from place to place (Gauss's theorem, Ch. 2: their volume integral
   is a surface flux). The exchange term sees only the symmetric part of the mean gradient:
   $\overline{u_iu_j}\,\partial U_i/\partial x_j=\overline{u_iu_j}\bar S_{ij}$ (Ch. 2 §2.10); for U(y) it is $\overline{uv}\,dU/dy$ — negative, a loss, because $\overline{uv}<0$ when $dU/dy>0$
   (C04)." **N80 [B]** $\frac{2\nu\bar S_{ij}\bar S_{ij}}{\overline{u_iu_j}(\partial U_i/\partial x_j)}\sim\frac{\nu(U/L)^2}{u_{rms}^2(U/L)}\sim\frac\nu{UL}=\frac1{\mathrm{Re}}\ll1$ (`mean_to_turbulent_dissipation_ratio(1e6)` → 1e-6).
7. `nb.derivation("D10", …)` — Part F D10 (15 steps, ★★★, `check_src`), ref "12.47". ⚠️ slip #2 and slip #10 callouts.
8. `nb.note` — **N81 [B]** table (term · sign · name · colour): $-\overline{u_iu_j}\,\partial U_i/\partial x_j$ shear production (orange; "+" here, "−" in the mean-flow budget of D09 above) · $\bar\varepsilon=2\nu\overline{S'_{ij}S'_{ij}}>0$ dissipation (rose; **not** small) · $g\alpha\overline{u_3T'}$ buoyant production (> 0, upward heat
   flux) or destruction (< 0) (blue). "This is the turbulent twin of Ch. 11's disturbance-energy equation $\frac{d}{dt}\int\frac12u_i^2\,dV=-\int u_iu_j\frac{\partial U_i}{\partial x_j}\,dV-\Lambda$ (11.88), where the
   same $-\overline{uv}\,dU/dy$ fed a growing wave." `ch12.shear_production`, `ch12.buoyant_production`. **N83 [B]** "Isotropic turbulence has no
   shear production: $\overline{u_iu_j}=\overline{u_1^2}\delta_{ij}$ contracts with $\partial U_i/\partial x_j$ to $\overline{u_1^2}\,\partial U_i/\partial x_i=0$." `assert` in code.
9. `nb.worked_example("production in a log layer", "u_* = 0.3 m/s, z = 10 m, κ = 0.4: $-\\overline{uw}=u_*^2$ = 0.09 m²/s², dU/dz = u_*/(κz) =
    0.075 s⁻¹. 1. Production = 0.09 × 0.075 = 6.75 × 10⁻³ m²/s³ (= u_*³/κz). 2. Night, heat flux $\\overline{wT'}$ = −0.02 K m/s, α = 1/300 K⁻¹,
    g ≈ 9.81: buoyancy term gα$\\overline{wT'}$ = −6.5 × 10⁻⁴ m²/s³ — it removes about 10 % of the production. 3. Steady, no transport:
    ε̄ = 6.75 × 10⁻³ − 0.65 × 10⁻³ = 6.1 × 10⁻³ m²/s³.")`
10. `nb.code` — `ch = ch12.channel_mixing_length(1000, kappa=0.41, A_plus=26.0)`; `bud = ch12.channel_energy_budget(1000, 0.41,
    26.0)`; `ch12.tke_budget(...)`, `ch12.mean_energy_budget(...)`. *expect:* production peaks at y⁺ ≈ 10.4 with P⁺ = 0.250
    (where viscous and Reynolds stress are equal — model; DNS puts it near y⁺ ≈ 12); integral identity pressure work = mean
    dissipation + production to 1e-6. *explain:* numbered.
11. `nb.check_agree` — **from scratch:** production $-\overline{uv}\,dU/dy$ and its integral by `np.gradient` and `np.trapezoid` on the model
    channel; `assert np.allclose(P_mine, bud["production"], rtol=1e-3)`.
12. `nb.figure` — the two budgets across the channel (wall units, semilog y⁺): top the mean-flow terms (pressure work
    purple, direct dissipation rose, loss to turbulence orange, negative), bottom the turbulence terms (production orange,
    positive — the mirror image; modelled sink rose). *see / read* ("the orange curves are mirror images: one term, two
    budgets") */ change* ("…Re_τ = 5200: the peak stays at y⁺ ≈ 10 but the direct dissipation's share of the total falls").
13. `nb.note` — **N82 [B]** (our Fig. 12.10) `ch12.mixing_potential_energy_change`: an unstably stratified layer mixed to a
    uniform state loses potential energy, which appears as turbulent energy; ⚠️ T is **potential** temperature; caption
    states the profile in both conventions ("$dT/dz=-15<\Gamma_a=-9.8$ K/km ⇔ $\Gamma_{met}=15>\Gamma_d=9.8$ K/km: unstable").
14. `nb.explainer("turbulent_energy_budget", heading="Where does turbulent energy come from and go?", why="Clicking a height
    shows the two budgets' bars side by side with the exchange term mirrored; sliding Re_τ shrinks the mean flow's direct
    viscous loss — 'same term, opposite sign' is seen, not read.", tries=["Click y⁺ ≈ 10: production is at its maximum, ¼ in
    wall units.", "Go to the centreline: production vanishes with the shear.", "Raise Re_τ from 180 to 5200 and watch the
    integrated direct dissipation's share fall.", "Turn the wall damping off and see where the peak moves."])`
15. `nb.md` — **What would change if…** "…we ask how *much* energy flows, and how small the eddies that finally dissipate it
    are? C07."

#### C07 — Kolmogorov scales (12.50) and the scale separation (12.51)
1. `nb.core("C07", "Kolmogorov scales $\\eta=(\\nu^3/\\bar\\varepsilon)^{1/4}$, $u_K=(\\nu\\bar\\varepsilon)^{1/4}$ (12.50) with $\\bar\\varepsilon\\sim(\\Delta U)^3/L$ (12.49) and
   $\\eta/L\\sim\\mathrm{Re}_L^{-3/4}$ (12.51)", question="If viscosity does the dissipating, why does the dissipation rate not depend on
   viscosity?")`
2. `nb.md` — **Plain words:** "Stir a cup of water and a cup of honey with the same spoon at the same speed until both are
   turbulent (you would need a much bigger spoon for honey). The power you put in is set by the spoon and the speed, not by
   the fluid. All of it ends as heat. So the fluid must arrange for viscosity to remove exactly that power — by making
   eddies small enough. **Climate hook:** in the atmosphere those eddies are about a millimetre across while a model grid
   box is kilometres wide: nine decades that can only be parameterised."
3. `nb.md` — **The idea** (the cascade ladder, **N87 [B]**):
   ```
   L (big eddies, speed ΔU) → L/2 → L/4 → …  each tier strained by the one above, no viscosity yet  … → η (viscous: Re = 1)
   energy passes down at the SAME rate ε̄ ~ ΔU³/L through every tier:   u′(l′) ~ (ε̄ l′)^{1/3},   Re(l′) = u′l′/ν ↓
   ```
4. `nb.note` — **N84 [B]** (`scales_sketch`, our Fig. 12.11) L = cross-stream extent of the velocity difference ΔU, so mean
   gradients are ~ ΔU/L.
5. `nb.derivation("D11", …)` — Part F D11 (9 steps), ref "12.50". Carries **N85 [B]**
   $\dot W\sim\overline{u_iu_j}(\partial U_i/\partial x_j)\sim(\Delta U)^2[\Delta U/L]=(\Delta U)^3/L$ (12.48), **N86 [B]** $\dot W=\bar\varepsilon$ so $\bar\varepsilon\sim(\Delta U)^3/L$ (12.49), **N88 [B]**
   $\eta u_K/\nu=1$, **N89 [B]** $\eta/L\sim\mathrm{Re}_L^{-3/4}$ with $\mathrm{Re}_L=\Delta UL/\nu$ (12.51).
6. `nb.worked_example("a stirred tank of water", "ΔU = 1 m/s, L = 1 m, ν = 10⁻⁶ m²/s. 1. Re_L = ΔU L/ν = 10⁶. 2. ε̄ ~ ΔU³/L = 1 m²/s³
   (1 W per kg). 3. η = (ν³/ε̄)^{1/4} = (10⁻¹⁸)^{1/4} = 3.2 × 10⁻⁵ m = 32 µm. 4. u_K = (νε̄)^{1/4} = 3.2 cm/s; τ_η = (ν/ε̄)^{1/2} = 1 ms. 5. Check:
   η u_K/ν = 3.2 × 10⁻⁵ × 3.2 × 10⁻²/10⁻⁶ = 1 ✓; η/L = 3.2 × 10⁻⁵ = (10⁶)^{−3/4} ✓ — 4.5 decades. 6. A simulation resolving it all needs
   ~ Re^{9/4} = 3 × 10¹³ grid points.")`
7. `nb.code` — `eps = ch12.dissipation_outer_scaling(1.0, 1.0)`; `eta, uK, tau_eta = ch12.kolmogorov_scales(1e-6, eps)`;
   `ch12.scale_separation(1e6)`; `ch12.dns_grid_points(1e6)`; `ch12.cascade_tiers(1.0, 1.0, 1e-6)`; `pd.DataFrame([ch12.scale_table(c)
   for c in cases])`. *expect:* 3.162e-5 m, 0.03162 m/s, 1e-3 s; η/L = 3.162e-5; 3.16e13; our table — atmospheric boundary
   layer (ΔU = 5 m/s, L = 1 km, ν = 1.5 × 10⁻⁵): ε̄ = 0.125 m²/s³, η = 0.41 mm, Re_L = 3.3 × 10⁸, 6.4 decades; ocean thermocline
   (0.1 m/s, 10 m, 10⁻⁶): ε̄ = 10⁻⁴, η = 0.32 mm; kitchen mixer (1 m/s, 5 cm): η = 15 µm. (**N90 [B]**; the book's own example
   stays private.)
8. `nb.check_agree` — **from scratch:** the 2 × 2 exponent solve for $\eta=\nu^a\bar\varepsilon^b$ with `np.linalg.solve([[2, 2], [-1, -3]], [1, 0])` →
   (0.75, −0.25); `assert np.allclose` with `DIM.solve_exponents` and with `kolmogorov_scales` for the atmospheric case.
9. `nb.animation` — **A4** (`player="frames"`, one frame per half-decade of Re_L from 10³ to 10⁸): the ladder on a log length
   axis (L fixed at the left, η sliding right, λ_T between), tiers added at the small end; beside it the model spectrum with
   its −5/3 range widening. *see / read* ("the left end never moves: a higher Re makes the cascade longer, not stronger") */
   change*.
10. `nb.note` — **N91 [B]** eliminating $\bar\varepsilon$ between $\bar\varepsilon=15\nu\overline{u^2}/\lambda_g^2$ (12.43) and $\bar\varepsilon\sim(\Delta U)^3/L$ (12.49):
    $\frac{(\Delta U)^3}L\propto\frac{\nu\overline{u^2}}{\lambda_T^2}\Rightarrow\frac{\lambda_T^2}{L^2}\propto\frac{\overline{u^2}}{(\Delta U)^2}\Big(\frac\nu{\Delta UL}\Big)\propto\frac1{\mathrm{Re}_L}$, or $\frac{\lambda_T}L\propto\mathrm{Re}_L^{-1/2}$ (12.52) — "λ_T
    is not the size of the dissipating eddies: $\bar\varepsilon=15\nu\overline{u^2}/\lambda_g^2$ (12.43) pairs it with the large-eddy velocity, not with $u_K$; $\lambda_T/\eta\sim\mathrm{Re}_L^{1/4}$".
    **N92 [B]** $R_\lambda\sim\mathrm{Re}_L^{1/2}$ and the ordering $\eta<\lambda_T<\Lambda_{(f\text{ or }g)}<L$ (`scale_ordering`; the book's thresholds private).
11. `nb.md` — **What would change if…** "…we ask how the energy is shared among the tiers in between? That is a spectrum, and
    dimensional analysis gives its shape (C08)."

#### C08 — Kolmogorov's −5/3 law (12.54) and the universal spectrum (12.53)
1. `nb.core("C08", "Kolmogorov's inertial-subrange law $S_{11}(k_1)=C_1\\bar\\varepsilon^{2/3}k_1^{-5/3}$ for $2\\pi/L\\ll k_1\\ll2\\pi/\\eta$ (12.54, slip #1
   corrected)", question="Why do spectra from a tidal channel, a jet and the atmosphere all fall with the same slope?")`
2. `nb.md` — **Plain words:** "Plot the spectrum of the wind on a mast, of the current in a tidal channel, of the flow behind
   a grid: on log–log paper, over the middle range of scales, all three are straight lines with the same slope. Something
   universal is going on. Eddies in that range are too small to know about the geometry and too big to feel viscosity; all
   they know is the rate at which energy is handed through them."
3. `nb.md` — **The idea:** "Only $\bar\varepsilon$ [m²/s³] and $k_1$ [1/m] are available; $S_{11}$ has units m³/s². There is exactly one way to
   build m³/s² from them." Table of the three ranges: energy-containing ($k_1\sim2\pi/L$, depends on the flow) · inertial
   (−5/3, universal) · dissipation ($k_1\eta\sim1$, universal in Kolmogorov units).
4. `nb.derivation("D12", …)` — Part F D12 (9 steps), ref "12.54". Carries **N93 [B]**
   $\frac{S_{11}(k_1)}{\nu^{5/4}\bar\varepsilon^{1/4}}=\Phi\big(\frac{k_1\nu^{3/4}}{\bar\varepsilon^{1/4}}\big)$, or $\frac{S_{11}(k_1)}{u_K^2\eta}=\Phi(k_1\eta)$ for $k_1\gg2\pi/L$ (12.53), and **N94 [B]**
   $\int_{-\infty}^{+\infty}S_{11}(k_1)\,dk_1=\overline{u_1^2}=\int_0^{+\infty}2S_{11}(k_1)\,dk_1$ (12.55).
5. `> ⚠️ **slip #1 — the book prints** $S_{11}(k_1)=\mathit{const}\cdot\bar\varepsilon^{2/3}\cdot k_1^{5/3}$ **; the correct form is** $S_{11}(k_1)=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$.` "Three
   checks: units (D12 step 5), the sentence after it calls it the $k^{-5/3}$ law, and a spectrum that *rose* with $k_1$ would
   have infinite variance." + `nb.code`: `inertial_spectrum_1d(k, 1.0, printed=True)` vs corrected, units printed with `pint`.
6. `nb.note` — **N95 [B]** three-dimensional spectrum $\bar e=\int_0^\infty S(K)\,dK$, $S(K)=C\bar\varepsilon^{2/3}K^{-5/3}$ with C ≈ 1.5 (Sreenivasan 1995);
   isotropy gives the one-sided one-dimensional constant $C_1=\tfrac{18}{55}C$ (`kolmogorov_constants(1.5)` → 0.491; two-sided
   0.245). **N98 [C]** the law was first confirmed in a tidal channel at very large Reynolds number (our analogue: the sea-water
   spectra of C15). **N99 [C]** universal small scales are why closure models and LES can work (C12, C13; Ch. 10).
7. `nb.primer("quadrature on a logarithmic grid", …)` (**P298**): "A spectrum spans decades; with $K=e^s$, $dK=K\,ds$, so
   $\int E\,dK=\int E\,K\,ds$ on an evenly spaced s." code: 4 lines integrating $K^{-5/3}$ from 1 to 10⁶ (exact 1.5(1 − 10⁻⁴)).
8. `nb.note` — **N97 [B]** a model spectrum for the whole range: Pao's $E(K)=C\bar\varepsilon^{2/3}K^{-5/3}\exp\{-\tfrac32C(K\eta)^{4/3}\}$ (Pao 1965), with
   Pope's energy-range factor when L is given (`model_spectrum`); two integral checks $\int E\,dK=\bar e$ and $2\nu\int K^2E\,dK=\bar\varepsilon$.
   *expect:* dissipation integral = ε̄ within 1 % (Pao's form satisfies it exactly).
9. `nb.worked_example("one point on the −5/3 line", "ε̄ = 1 m²/s³, k₁ = 100 rad/m (a 6 cm eddy), two-sided constant C₁ = (9/55) × 1.5
   = 0.245. 1. ε̄^{2/3} = 1. 2. k₁^{−5/3} = 100^{−5/3} = 4.64 × 10⁻⁴. 3. S₁₁ = 0.245 × 4.64 × 10⁻⁴ = 1.14 × 10⁻⁴ m³/s². 4. Ten times the
   wavenumber: 10^{−5/3} = 1/46.4 — the spectrum falls 46-fold per decade. 5. Units: (m²/s³)^{2/3}(1/m)^{−5/3} = m^{4/3+5/3}s⁻² = m³/s² ✓.")`
10. `nb.code` — `K = np.logspace(0, 6, 400)`; `E = ch12.model_spectrum(K, 1.0, 1e-6, L=1.0, kind="pope")`; `ch12.fit_inertial_range(K,
    E, band=(50, 2000))`; `ch12.kolmogorov_normalize_spectrum`. *expect:* slope −1.667 ± 0.02 in the band.
11. `nb.check_agree` — **from scratch:** slope by `np.polyfit(np.log(K[m]), np.log(E[m]), 1)`; `assert abs(slope + 5/3) < 0.02`;
    `assert np.isclose(slope, ch12.fit_inertial_range(K, E, (50, 2000))[0])`.
12. `nb.plotly` — **F3** (**N96 [B]**, our Fig. 12.12) `slider_figure` over Re_L (10³…10⁸, 21 steps): model spectra in
    Kolmogorov scaling $S/(u_K^2\eta)$ vs $K\eta$ with the −5/3 line and the inertial range shaded; + `nb.live`: free ΔU, L, ν.
    *see / read* ("the right-hand end is the same curve for every Re; the plateau on the left moves out") */ change*.
13. `nb.explainer("energy_cascade_spectrum", heading="What does a higher Reynolds number change?", why="Dragging Re_L slides η
    away from L, adds tiers and widens the −5/3 range while its level stays put; changing ν at fixed ΔU and L moves the
    right end only.", tries=["Switch air ↔ water at the same ΔU and L: ε̄ does not move, η does.", "Load 'atmospheric
    boundary layer': 6.4 decades between L and η.", "Toggle the printed +5/3 ghost and read the units inspector.", "Double
    ΔU: ε̄ × 8, η × 2^{−3/4} = 0.59."])`
14. `nb.md` — **What would change if…** "…the turbulence lives in a jet or a wake, where the big eddies and the mean flow
    change downstream? Conservation and symmetry still give the answer (C09)."

---

### A.8 §12.8 Free Turbulent Shear Flows — C09
#### C09 — The self-similar plane jet (12.66)
1. `nb.section("12.8", "Free Turbulent Shear Flows", intro="**What is this section about?** Jets, wakes, plumes and mixing
   layers — turbulence with no wall nearby. Far downstream they forget their source and keep one shape that only stretches.
   One conserved quantity plus that shape-keeping gives how fast they widen and slow down, without any turbulence model.")`
   ⚠️ line: "From here to §12.10 the density is constant and P is the deviation from hydrostatic (convention 6)."
2. `nb.core("C09", "The far field of the plane turbulent jet $U(x,y)=C_5(J_s/\\rho)^{1/2}x^{-1/2}F(y/x)$ (12.66)", question="How can
   anyone predict how a turbulent jet spreads without solving for the turbulence?")`
3. `nb.md` — **Plain words:** "Smoke from a chimney, the exhaust of an engine, a river entering the sea, the wake of an island
   in the trade winds: each widens downstream by swallowing the still fluid around it (**N101 [B]** *entrainment*), and each
   looks the same at every distance once you rescale it (*self-preservation*). **N100 [B]**: jet, wake, mixing layer; a plume
   is a jet driven by buoyancy (figure `free_shear_sketch` from `ch12.free_shear_flow` with labelled illustrative constants).
   A *free* shear flow has no wall; §12.9 adds one."
4. `nb.md` — **The idea:**
   ```
   conservation:   momentum flux  J_s = ρ ∫U² dy   is the same at every x        →  U_CL² · δ = const
   shape-keeping:  U = U_CL(x) · F(y/δ(x))  must satisfy the equations at every x →  dδ/dx = const
   together:       δ ∝ x,   U_CL ∝ x^{-1/2},   volume flux ∝ U_CL·δ ∝ x^{+1/2}  (it grows: entrainment)
   ```
5. `nb.note` — **N102 [B]** $U(x,y)=U_{CL}(x)F(y/\delta(x))$ (12.56), F even, F(0) = 1; **N103 [B]** $-\overline{uv}=\Psi(x)G(y/\delta(x))$ (12.57), G odd —
   ⚠️ "these F, G are jet profiles, not the tensor functions of $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40)"; **N104 [B]** $\partial U/\partial x+\partial V/\partial y=0$ (12.58); **N105 [B]**
   $U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}=-\frac1\rho\frac{\partial P}{\partial x}+\nu\big(\frac{\partial^2U}{\partial x^2}+\frac{\partial^2U}{\partial y^2}\big)-\frac{\partial\overline{u^2}}{\partial x}-\frac{\partial\overline{uv}}{\partial y}$ (12.59); **N106 [B]**
   $U\frac{\partial V}{\partial x}+V\frac{\partial V}{\partial y}=-\frac1\rho\frac{\partial P}{\partial y}+\nu\big(\frac{\partial^2V}{\partial x^2}+\frac{\partial^2V}{\partial y^2}\big)-\frac{\partial\overline{uv}}{\partial x}-\frac{\partial\overline{v^2}}{\partial y}$ (12.60) (`rans_2d_residual`).
6. `nb.derivation("D13", …)` — Part F D13 (10 steps), ref "12.62". Carries **N107 [B]**
   $U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}\cong-\frac{\partial\overline{uv}}{\partial y}$, $0\cong-\frac1\rho\frac{\partial}{\partial y}(P+\rho\overline{v^2})$ (12.61), **N108 [B]** the conservative form, **N109 [B]**
   $J_s\equiv\rho_s\int_{-\infty}^{+\infty}[U^2]_{x=0}dy\cong\rho\int_{-\infty}^{+\infty}U^2dy=\mathit{const.}$ (12.62). `thin_shear_layer_terms` prints the sizes of kept and dropped terms.
7. `nb.derivation("D14", …)` — Part F D14 (14 steps, ★★★, `check_src`), ref "12.63". Carries **N110 [B]** (V eliminated) and
   **N111 [B]** $\big\{\frac{\delta U'_{CL}}{U_{CL}}\big\}F^2-\big\{\frac{\delta U'_{CL}}{U_{CL}}+\delta'\big\}F'\int_0^\xi F\,d\xi=\big\{\frac\Psi{U_{CL}^2}\big\}G'$ (12.63).
8. `nb.derivation("D15", …)` — Part F D15 (10 steps), ref "12.66". Carries **N112 [B]** $\frac{\delta U'_{CL}}{U_{CL}}=C_1$,
   $\frac{\delta U'_{CL}}{U_{CL}}+\delta'=C_2$, $\frac\Psi{U_{CL}^2}=C_3$ (12.64), **N113 [B]** $\delta=(C_2-C_1)(x-x_o)$ (virtual origin $x_o$, of order the slot width;
   `virtual_origin_fit`), **N114 [B]** $J_s=\rho\int_{-\infty}^{+\infty}U^2dy=\rho U_{CL}^2\delta\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi=\rho C_4^2x^{2\gamma+1}\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi$ (12.65) ⇒ $2\gamma+1=0$, **N115 [B]** $-\overline{uv}=C_3U_{CL}^2G(y/x)=C_3C_5^2(J_s/\rho)x^{-1}G(y/x)$ (12.67) with the
   slip #4 callout, **N116 [B]** $\dot V(x)=\int_{-\infty}^{+\infty}U\,dy=C_5(J_s/\rho)^{1/2}x^{+1/2}\int_{-\infty}^{+\infty}F(\xi)\,d\xi$ (12.68).
9. `nb.note` — **N122 [B]** Gaussian fit $F(y/x)=\exp\{-\ln(2)(y/x)^2/(\xi_{1/2})_U^2\}$ (`gaussian_profile`, `profile_integrals`; Gaussian
   integral P194): $\int F\,d\xi=\xi_{1/2}\sqrt{\pi/\ln2}$, $\int F^2d\xi=\xi_{1/2}\sqrt{\pi/(2\ln2)}$; with a Gaussian F the invariant fixes
   $C_5=\big(\int F^2d\xi\big)^{-1/2}$.
10. `nb.worked_example("a plane jet with an illustrative half-width", "Take J_s/ρ = 1 m³/s² and an **illustrative** half-width
    ξ½ = 0.10 (not the book's value). 1. ∫F²dξ = 0.10 × √(π/(2 ln 2)) = 0.1505, so C₅ = 0.1505^{−1/2} = 2.577. 2. Centreline
    speed at x = 1 m: U_CL = 2.577 m/s; at x = 4 m: 2.577/√4 = 1.289 m/s. 3. Half-width: 0.10 m at x = 1 m, 0.40 m at x = 4 m.
    4. Momentum flux at both stations: ρ U_CL² x ∫F²dξ = ρ × 1.000 ✓ (the same). 5. Volume flux: ∫F dξ = 0.2129; at x = 1 m
    V̇ = 2.577 × 0.2129 = 0.549 m²/s; at x = 4 m it has doubled to 1.097 m²/s — the extra fluid was entrained. 6. Entrainment
    velocity at x = 1 m: ½ dV̇/dx = V̇/(4x) = 0.137 m/s on each side.")`
11. `nb.code` — `U = ch12.plane_jet_mean_velocity(x, y, 1.0, 1.0, C5="from_invariant", xi_half=0.10)` at x = 1, 2, 4 m;
    `ch12.jet_momentum_flux_per_span(y, U, 1.0)`; `ch12.plane_jet_volume_flux`; `ch12.plane_jet_entrainment_velocity`. *expect:* J_s
    = 1.000 at all three; V̇ = 0.549, 0.776, 1.097; v_e(1 m) = 0.137. *explain:* numbered.
12. `nb.check_agree` — **from scratch:** $\rho\int U^2dy$ and $\int U\,dy$ by `np.trapezoid` at the three stations; `assert np.allclose(J,
    1.0, rtol=1e-6)`; `assert np.allclose(V[2]/V[0], 2.0, rtol=1e-6)`.
13. `nb.plotly` — **F4** `slider_figure` over x (0.5…8 m, 24 steps): left the raw profile U(x, y) (turbulent, purple; the
    laminar jet of Ch. 9, $\delta\propto x^{2/3}$, $U_{CL}\propto x^{-1/3}$, muted, `JET.free_jet`), right $U/U_{CL}$ vs y/x — one curve. *see / read /
    change* ("…the jet were laminar: it widens more slowly at first and its exponents are 2/3 and −1/3, because viscosity
    — not the eddies — sets the spreading").
14. `nb.note` — the scalar field: **N117 [B]** $\bar Y(x,y)=Y_{CL}(x)H(y/x)$ (12.69); **N118 [B]**
    $\dot M_s=\rho_s\int_{-\infty}^{+\infty}[U]_{x=0}dy\cong\rho\int_{-\infty}^{+\infty}\bar YU\,dy=\rho Y_{CL}C_5(J_s/\rho)^{1/2}x^{+1/2}\int_{-\infty}^{+\infty}H(\xi)F(\xi)\,d\xi$ (12.70) with the slip #3 callout;
    **N119 [B]** $\bar Y=C_6\big(\dot M_s/\sqrt{\rho J_s}\big)x^{-1/2}H(y/x)$ (12.71); **N120 [B]** $U=C_5U_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}F(y/x)$ (12.72); **N121 [B]**
    $\bar Y=C_7Y_0(\rho_s/\rho)^{1/2}(x/d)^{-1/2}H(y/x)$ (12.73). `scalar_flux_per_span` flat in x. *expect:* flat to 1e-6.
15. `nb.derivation("D16", …)` — Part F D16 (12 steps), ref "Table 12.1 (exponents)". **N124 [B]**. Then **N123 [B]**: our own
    table from `pd.DataFrame({f: ch12.free_shear_exponents(f) for f in flows})` (exact fractions: width, velocity, scalar,
    local Reynolds number) — "the exponents are mathematics and public; the book's constants and half-widths are data and
    are not reproduced; our examples pass labelled illustrative constants." *expect:* plane jet (1, −1/2, −1/2, 1/2); round
    jet (1, −1, −1, 0); plane wake (1/2, −1/2, −1/2, 0); round wake (1/3, −2/3, −2/3, −1/3); plane plume (1, 0, −1, 1); round
    plume (1, −1/3, −5/3, 2/3); shear layer (1, 0, –, 1).
16. `nb.note` — **N125 [B]** more general similarity: $\frac{\delta U'_{CL}}{U_{CL}}=C_8\big(\frac{\delta U'_{CL}}{U_{CL}}+\delta'\big)=C_9\frac\Psi{U_{CL}^2}$ (12.74), e.g. $\delta\sim x^m$,
    $U_{CL}\sim x^n$, $\Psi\sim x^{2n+m-1}$ — so the constants need not be universal. `> ⚠️ **slip #5**` + `general_similarity_check` on
    both families. *expect:* power family coefficients (n, m + n, 1) × $x^{m-1}$; exponential family middle coefficient 0 and
    $U_{CL}^2\delta\propto e^{-ax}$ (not constant).
17. `nb.note` — **N126 [B]** (Example 12.2 with our own numbers): a round methane jet (d = 5 mm, U₀ = 30 m/s) into air:
    `stoichiometric_mass_fraction(16.04, 28.97, 2, 0.21)` → 0.055; `round_jet_distance_for_mass_fraction(0.055, 0.005, 0.656,
    1.184, 1.0, C_Y=4.0)` with C_Y an **illustrative** constant (not a tabulated one; recompute this row's numbers in the cell); mole vs mass fraction glossed. The book's inputs and
    answers stay private.
18. `nb.note` — **N127 [C]** profiles of $\bar e$ and the stresses across the jet (qualitative panel; $\overline{uv}=0$ on the axis by
    symmetry, largest near the steepest mean shear). **N128 [B]**
    $0=-U\frac{\partial\bar e}{\partial x}-V\frac{\partial\bar e}{\partial y}-\overline{uv}\frac{\partial U}{\partial y}-\frac{\partial}{\partial y}\big(\frac1{\rho_0}\overline{pv}+\frac12\overline{u_i^2v}\big)-\bar\varepsilon$ (12.75, written with the triple correlation $\tfrac12\overline{u_i^2u_j}$ of the exact budget of C06; ⚠️ slip #10 callout: the page prints $\tfrac12\overline{ev}$). **N129 [C]** three balances from `jet_tke_budget` (labelled model): axis — advection ≈ dissipation;
    mid-layer — production ≈ dissipation; edge — transport ≈ advection.
19. `nb.note` — **N130 [B]** (ours) "Close the jet with a constant eddy viscosity $\nu_T\propto U_{CL}\delta$: the equation is then the
    laminar jet's with ν → ν_T, so $F=\mathrm{sech}^2(a\xi)$ — the laminar shape with turbulent exponents"
    (`plane_jet_eddy_viscosity_profile` vs the Gaussian: they differ by < 5 % of $U_{CL}$ inside the half-width). Forward
    pointer: C12.
20. `nb.explainer("turbulent_jet_similarity", heading="How does a jet spread without a turbulence model?", why="Toggling raw ↔
    rescaled makes five profiles collapse onto one F(ξ); a wrong-exponent slider breaks either the collapse or the invariant
    bar; switching the flow changes the invariant and the exponents follow.", tries=["Toggle 'rescaled': five curves
    become one.", "Set the trial decay exponent to −0.4: the momentum bar grows as x^{0.2} — not allowed.", "Switch to 'round
    jet': the velocity now falls as 1/x and the local Reynolds number stays constant.", "Switch to 'plane wake' and watch
    the width grow as √x."])`
21. `nb.md` — **What would change if…** "…there is a wall? Then a second length — a tiny viscous one — enters, the flow is
    never Reynolds-number independent, and the argument must be made twice (C10, C11)."

---

### A.9 §12.9 Wall-Bounded Turbulent Shear Flows — R03 R04, C10, C11
#### C10 — The law of the wall (12.80) in wall units
1. `nb.section("12.9", "Wall-Bounded Turbulent Shear Flows", intro="**What is this section about?** Turbulence next to a wall:
   pipes, channels, the boundary layer on a wing, the wind over the ground. Two lengths matter — a viscous one of microns
   to millimetres at the wall and the thickness of the whole layer — and where the two descriptions overlap the mean
   velocity is a logarithm.")`
2. `nb.core("C10", "The law of the wall $U^+\\equiv U/u_*=f(yu_*/\\nu)=f(y^+)$ (12.80) with $u_*^2\\equiv\\tau_0/\\rho$ (12.81)", question="Why do
   profiles from different flows and Reynolds numbers fall on one curve near a wall?")`
3. `nb.md` — **Plain words:** "Measure the wind 1 mm above a smooth plate in a wind tunnel, the water speed 0.1 mm from a
   pipe wall, the flow just above a ship's hull. Plotted in m/s against mm they have nothing in common. Divide each speed
   by a 'friction velocity' and each distance by a 'viscous length', both built from the wall stress alone, and they fall
   on one curve. **N131 [B]**: a wall brings two length scales, $l_\nu\ll\delta$, and — unlike a jet — the flow never becomes
   independent of the Reynolds number on a smooth wall."
4. `nb.md` — **The idea:** `wall_layers_sketch` + table | layer | what matters | length | velocity law |: viscous sublayer
   (ν, τ₀; y⁺ < 5; U⁺ = y⁺) · buffer (both stresses; 5–30) · logarithmic / overlap (neither length; C11) · wake / outer (δ; C11).
5. `nb.note` — **N132 [B]** + `nb.figure` (our Fig. 12.16): mean turbulent channel profile (`WT.composite_profile`, κ = 0.41,
   B = 5.0, Π = 0 labelled illustrative) vs the laminar parabola (`LAM.channel_flow`) at the same bulk speed: blunter in the
   core, far steeper at the wall. *see / read / change*.
6. `nb.derivation("D17", …)` — Part F D17 (9 steps), ref "12.76". Carries **N133 [B]** $0=-\frac{\partial P}{\partial x}+\frac{\partial\bar\tau}{\partial y}$,
   $0=-\frac{\partial}{\partial y}(P+\rho\overline{v^2})$, $\bar\tau=\mu\frac{\partial U}{\partial y}-\rho_0\overline{uv}$ (12.76); **N134 [B]** $P(x,y)-P(x,0)=-\rho\overline{v^2}$; **N135 [B]**
   $\frac{\partial}{\partial x}P(x,y)-\frac{d}{dx}P(x,0)=-\rho\frac{\partial}{\partial x}\overline{v^2}=0$ (12.77); **N136 [B]** $\bar\tau(y)=\tau_0(1-2y/h)$ (h = **full** height).
7. `nb.recap("R03", "Channel: pressure gradient fixed by the wall stress", "Ch. 8's force balance on a slug of fluid (D07)
   holds for the turbulent mean too: $dP/dx=-2\\tau_o/h$ (12.90) — the last step of D17. (⚠️ slip #7: the book cites Exercise
   12.31 for the proof; it is Exercise 12.32.)", where="Ch. 8 §8.2, D07")` + `WT.channel_pressure_gradient(0.3, 0.1)` → −6 Pa/m.
8. `nb.recap("R04", "Pipe: the same balance", "$dP/dx=-4\\tau_o/d$ (12.91), Ch. 8's $\\tau_0=\\frac a2\\frac{dp}{dz}$ (8.8) with a = d/2 and z → x. (There $\\tau_0$ carries the sign of the pressure gradient; here $\\tau_0>0$ is its magnitude, hence the minus sign.)", where="Ch. 8 §8.2")` + parity `assert np.isclose(WT.wall_stress_from_pressure_gradient(-12.0, 0.1),
   LAM.pipe_wall_stress(...))`.
9. `nb.note` — **N137 [B]** + `nb.figure` (our Fig. 12.17): `WT.stress_partition` across a channel at Re_τ = 1000: total
   (straight line, black), viscous part (rose, confined to y⁺ ≲ 30), Reynolds part (orange); inset: a zero-pressure-gradient
   boundary layer from `boundary_layer_stress_from_profile` — nearly constant stress near the wall. **N138 [B]**
   $U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}=-\frac1\rho\frac{\partial P}{\partial x}+\frac1\rho\frac{\partial\bar\tau}{\partial y}$ (12.78) (laminar form: Ch. 9's $u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2u}{\partial y^2}$ (9.9)). *see / read / change*.
10. `nb.primer("inner, outer and overlap: two descriptions that must agree where both hold", …)` (**P300**): "Close to the
    wall one set of variables works (inner), far away another (outer). If there is a region where both are valid, the two
    formulas must give the same answer there — *matching*. It is the separation argument (P167) applied to two
    approximations." code: 4 lines — f_in(y) = y/(1 + y) rescaled both ways and compared in the middle.
11. `nb.note` — **N139 [B]** inner layer (viscosity matters; length y or $l_\nu$), outer layer (it does not; length δ), overlap
    (`WT.layer_name`). **N140 [B]** the inner list $U=U(\rho,\tau_0,\nu,y)$ (12.79) — δ and $U_\infty$ deliberately left out.
12. `nb.derivation("D18", …)` — Part F D18 (7 steps), ref "12.80". Carries **N141 [B]** $u_*^2\equiv\tau_0/\rho$ (12.81), $l_\nu=\nu/u_*$, $\delta^+=\delta u_*/\nu$ and
    **N142 [B]** $\mu(dU/dy)=\tau_0\Rightarrow U=\tau_0y/\mu$, or $U^+=y^+$ (12.82).
13. `nb.worked_example("wall units for air over a plate", "τ₀ = 0.3 Pa, ρ = 1.2 kg/m³, ν = 1.5 × 10⁻⁵ m²/s. 1. u_* = √(0.3/1.2) =
    0.5 m/s. 2. l_ν = ν/u_* = 3 × 10⁻⁵ m = 30 µm. 3. The sublayer (y⁺ < 5) is 0.15 mm thick. 4. At y = 0.09 mm: y⁺ = 3, U⁺ = 3, so
    U = 1.5 m/s. 5. A 3 cm boundary layer has δ⁺ = 0.03/3 × 10⁻⁵ = 1000.")`
14. `nb.code` — `yp, Up, us, lnu = WT.wall_units(y, U, 0.3, 1.2, 1.5e-5)`; `WT.friction_reynolds_number(0.5, 0.03, 1.5e-5)`;
    `WT.law_of_the_wall_groups()`; `WT.viscous_sublayer(3.0)`. *expect:* 0.5, 3e-5, 1000; two Π groups $U/\sqrt{\tau_0/\rho}$ and
    $y\sqrt{\tau_0/\rho}/\nu$; 3.0.
15. `nb.check_agree` — **from scratch:** four lines (u_*, l_ν, y⁺, U⁺) by hand; `assert np.allclose`.
16. `nb.primer("semi-log axes: a logarithm is a straight line", …)` (**P299**): "On `ax.semilogx`, $U^+=\frac1\kappa\ln y^++B$ is a
    straight line; it rises by $\ln(10)/\kappa=2.303/\kappa$ per decade (5.6 for κ = 0.41). The linear law U⁺ = y⁺ looks curved there."
    code: 3 lines printing U⁺ at y⁺ = 10, 100, 1000 and the differences (5.62).
17. `nb.figure` — the law of the wall on semi-log axes: Spalding's curve (`WT.spalding_uplus`, κ = 0.41, B = 5.0) bold, U⁺ = y⁺
    (rose dashed), the log line (teal dashed, anticipating C11), bands for the layers. *see / read* ("the sublayer line is
    accurate to 5 % up to y⁺ ≈ 5–7") */ change* ("…the wall were rough: the sublayer disappears and viscosity with it (C11)").
18. `nb.md` — **What would change if…** "…we go far from the wall, where viscosity cannot matter? A second law holds there,
    and making the two agree gives the logarithm (C11)."

#### C11 — The logarithmic law (12.88) by overlap matching
1. `nb.core("C11", "The logarithmic law $U^+=\\frac1\\kappa\\ln(y^+)+B$ (12.88); rough wall $U^+=\\frac1\\kappa\\ln(y/y_0)$ (12.93)", question="Where does
   the logarithm come from?")`
2. `nb.md` — **Plain words:** "A wind profile over open country, plotted against the logarithm of height, is a straight
   line — the fact behind every wind-turbine siting study, every 'wind at 10 m' in a weather report and every bulk formula
   for surface drag in a climate model. No turbulence model is needed to get it: only the statement that two descriptions
   of the same flow must agree."
3. `nb.md` — **The idea:** "Near the wall: velocity depends on y and the viscous length, not on δ. Far away: on y and δ, not
   on viscosity. In between, the *gradient* can depend on neither length. The only thing left with the right units is
   $dU/dy=u_*/(\kappa y)$ — and the integral of 1/y is a logarithm."
4. `nb.derivation("D19", …)` — Part F D19 (12 steps), ref "12.88". Carries **N143 [B]** $U=U(\rho,\tau_0,\delta,y)$ (12.83), **N144 [B]**
   $\frac{U_\infty-U}{u_*}=F(y/\delta)=F(\xi)$ (12.84), **N145 [B]** $\frac{dU}{dy}=\frac{u_*^2}\nu\frac{df}{dy^+}$ (12.85), **N146 [B]** $-\frac{dU}{dy}=\frac{u_*}\delta\frac{dF}{d\xi}$ (12.86), **N147
   [B]** $-\xi\frac{dF}{d\xi}=y^+\frac{df}{dy^+}$ (12.87), **N148 [B]** $F(\xi)=-\frac1\kappa\ln(\xi)+A$ (12.89). ⚠️ "κ here is the von Kármán constant, not the thermal diffusivity of $\frac{\partial\bar T}{\partial t}+U_j\frac{\partial\bar T}{\partial x_j}+\frac{\partial}{\partial x_j}(\overline{u_jT'})=\kappa\frac{\partial^2\bar T}{\partial x_j^2}$ (12.31), our $\kappa_{th}$."
5. `nb.worked_example("reading a log law", "Illustrative pair κ = 0.41, B = 5.0. 1. At y⁺ = 100: U⁺ = ln(100)/0.41 + 5.0 = 11.23 +
   5.0 = 16.23. 2. At y⁺ = 1000: 21.85 — one decade adds 2.303/0.41 = 5.62. 3. Where do U⁺ = y⁺ and the log law cross? Solve
   y⁺ = ln(y⁺)/0.41 + 5.0: y⁺ ≈ 10.8 — the middle of the buffer layer. 4. With u_* = 0.5 m/s, l_ν = 30 µm: at y = 3 mm (y⁺ = 100) the
   mean speed is 16.23 × 0.5 = 8.1 m/s.")`
6. `nb.code` — `WT.log_law(100.0, kappa=0.41, B=5.0)`; `WT.log_law_crossing(kappa=0.41, B=5.0)`; `WT.friction_law_from_overlap(1000,
   kappa=0.41, A=1.0, B=5.0)`; `ch12.overlap_matching_sympy()`. *expect:* 16.232; 10.80; 22.85 (= ln(1000)/0.41 + 6.0).
7. `nb.primer("composite profile: inner law + outer correction", …)` (**P301**): "Add to an inner formula a correction that
   vanishes near the wall and grows to its full size at the edge: the sum is usable across the whole layer." code: 3 lines
   evaluating `WT.coles_wake(np.array([0, 0.5, 1.0]))` → 0, 0.5, 1.
8. `nb.note` — **N153 [B]** Spalding's single formula for the whole inner layer,
   $y^+=U^++e^{-\kappa B}\big[e^{\kappa U^+}-1-\kappa U^+-\tfrac12(\kappa U^+)^2-\tfrac16(\kappa U^+)^3\big]$ (`spalding_uplus` by `brentq`, P108; *expect* U⁺(1, 5, 12, 30, 100) = 1.000, 4.866,
   9.157, 12.634, 16.077 for κ = 0.41, B = 5.0). **N154 [B]** Coles' outer profile $U^+=\frac1\kappa\ln(y^+)+B+\frac{2\Pi}\kappa W(y/\delta)$,
   $W=3\xi^2-2\xi^3$ or $\sin^2(\pi\xi/2)$ (`composite_profile`; Π passed explicitly — the book's value is private). **N155 [C]** open
   questions (a power-law overlap; stress-gradient layers); the sublayer is universal, the wake is not — hence slightly
   different constants per flow. **N156 [B]** `pd.DataFrame(WT.LOG_LAW_CONSTANTS)` — cited public values only; κ and B are
   required keywords everywhere. **N157 [B]** an empirical link across flows and pressure gradients,
   $\kappa B=1.6[\exp(0.1663B)-1]$ (12.92) (`nagib_chauhan_kappa(5.0)`).
9. `nb.plotly` — **F5** (**N149 [B]**, our Fig. 12.18) `slider_figure` over Re_τ (180…10⁵, 20 steps): `composite_profile` in wall
   units on a semi-log axis with the layer bands; dots: public channel DNS (Lee & Moser 2015, `reference/ch12/`) at Re_τ =
   180, 1000, 5200. *see / read* ("the left part never moves; the wake peels off at y⁺ ≈ 0.15 δ⁺") */ change*.
10. `nb.check_agree` — **from scratch:** κ and B from a straight-line fit of U⁺ against ln y⁺ (`np.polyfit`) on a composite
    profile at δ⁺ = 5000 inside 30 < y⁺ < 0.15 δ⁺; `fit = WT.fit_log_law(yp, Up, Re_tau=5000)`; `assert np.allclose((1/slope, intercept), (fit["kappa"], fit["B"]), rtol=1e-8)`.
    *expect:* **pending verify ruling (`reports/ch12_verification.md`, flagged item 1).** The first draft of this row expected
    κ = 0.41 ± 0.01, B = 5.0 ± 0.2; the implemented function returns κ = 0.383, B = 3.99 on this window, because Spalding's
    composite profile is still about 0.66 below the log law at y⁺ = 30, and 0.399, 4.56 when the window starts at y⁺ = 100.
    The builder takes the window and the expected pair from the verifier's ruling and never edits them to make the cell
    pass; the assertion of this row (hand fit = library fit) holds whatever the ruling. **Teaching point** (one sentence in
    "What does this show?", with a two-row table of the fitted pair for a window starting at y⁺ = 30 and at y⁺ = 100): "The
    constants you fit depend on where the window starts. At y⁺ = 30 the buffer layer has not quite joined the logarithm, so
    a window that starts there returns a smaller κ and B than the pair the profile was built with — one reason published
    values differ from one experiment to the next (N155, N156). The wake bends the top of the window in the same way."
11. `nb.figure` — the indicator $y^+dU^+/dy^+$ (`WT.log_law_indicator`) at δ⁺ = 180, 1000, 10⁴: a plateau at 1/κ = 2.44 appears
    only when δ⁺ is large. *see / read / change*.
12. `nb.derivation("D20", …)` — Part F D20 (5 steps), ref "12.93": **N159 [B]**. **N158 [B]** figure (our Fig. 12.19): `WT.log_law`
    and `WT.rough_wall_log_law` on one semi-log plot — smooth: the logarithm sits on a sublayer; rough: the elements poke
    through it and the logarithm extrapolates to U = 0 at $y_0$ ("hydrodynamically smooth / rough").
13. `nb.worked_example("wind over grass", "z₀ = 0.03 m, u_* = 0.4 m/s, κ = 0.41. 1. U(10 m) = (0.4/0.41) ln(10/0.03) = 0.976 × 5.81
    = 5.67 m/s. 2. Neutral drag coefficient at 10 m: C_D = [κ/ln(z/z₀)]² = (0.41/5.81)² = 4.98 × 10⁻³. 3. Check: τ₀ = ρu_*² =
    1.2 × 0.16 = 0.192 Pa = ρ C_D U² = 1.2 × 4.98 × 10⁻³ × 32.1 ✓. **Climate hook:** this C_D is the neutral value every bulk
    surface-flux formula starts from; C15 corrects it for stability.")` + `WT.rough_wall_log_law`, `WT.friction_velocity_from_wind`,
    `WT.drag_coefficient_neutral`.
14. `nb.note` — **N150 [B]** zero-pressure-gradient boundary-layer fits (displayed equations of the book; orchestrator's open
    decision (a) of the curation): $\theta\approx0.016\,x\,\mathrm{Re}_x^{-0.15}$, $\delta^*\approx\theta\exp\{7.11\kappa/\ln(\mathrm{Re}_\theta)\}$,
    $\delta_{99}=0.2\,\delta^*[\kappa^{-1}\ln(\mathrm{Re}_{\delta^*})+3.30]$, $C_f\cong2.0/[\kappa^{-1}\ln(\mathrm{Re}_{\delta^*})+3.30]^2$ (`zpg_boundary_layer`; ours: $C_f=2/(U_\infty^+)^2$). **N151 [B]**
    $C_f\cong0.370(\log_{10}\mathrm{Re}_x)^{-2.584}$ and $C_f\cong0.455/[\ln(0.06\,\mathrm{Re}_x)]^2$ (`skin_friction_zpg(law=)`; log₁₀ vs ln flagged). **N152 [B]**
    (Example 12.3, our own numbers): x = 3 m, $U_\infty$ = 60 m/s, ν = 1.5 × 10⁻⁵ ⇒ Re_x = 1.2 × 10⁷; `zpg_boundary_layer(3.0, 60.0,
    1.5e-5, kappa=0.384)` beside the Blasius values. *expect:* θ ≈ 4.2 mm (the cell prints δ*, δ₉₉, C_f; the builder records
    them). **N161 [B]** one log–log figure: `BL.blasius_skin_friction` vs `skin_friction_zpg` over Re_x = 10⁴…10¹⁰.
15. `nb.note` — **N160 [B]** (Exercise 12.34, stated): integrating the log law over a pipe's cross-section (one integration
    by parts, P218a) gives $U_{av}\cong u_*[(1/\kappa)\ln(au_*/\nu)+B-3/(2\kappa)]$ and, with Darcy's $f_D=8(u_*/U_{av})^2$, Prandtl's law
    $f_D^{-1/2}=2.0\log_{10}(\mathrm{Re}_df_D^{1/2})-0.8$; `pipe_bulk_velocity_loglaw` against quadrature; `pipe_friction_factor_turbulent(1e5)`
    vs laminar 64/Re (`LAM.pipe_friction_factor`). *expect:* derived slope 1.99 for κ = 0.41 (vs 2.0); f_D(10⁵) ≈ 0.018.
16. `nb.explainer("law_of_the_wall", heading="Where does the logarithm come from?", why="Sliding Re_τ stretches the log region
    while the inner curve stays fixed in wall units; switching to outer units makes the opposite collapse; the indicator
    shows the plateau 1/κ appear only when δ⁺ is large.", tries=["Raise Re_τ from 180 to 5200: the straight part grows
    from nothing to more than a decade.", "Switch to outer scaling: now the wakes collapse and the wall region fans out.",
    "Click y⁺ = 12: buffer layer — neither law holds.", "Switch to a rough wall and change y₀⁺."])`
17. `nb.md` — **What would change if…** "…we tried to *compute* this profile from the averaged equations? We must guess the
    Reynolds stress. The simplest guess that returns the logarithm is next (C12)."

---

### A.10 §12.10 Turbulence Modeling — C12, C13
#### C12 — Eddy viscosity (12.94) and the mixing length (12.99)–(12.101)
1. `nb.section("12.10", "Turbulence Modeling", intro="**What is this section about?** Closing the gap: replacing the unknown
   Reynolds stress by a formula in terms of the mean flow. We build the simplest such model, see exactly what its one
   constant does, then meet the two-equation model used in most engineering codes.")`
2. `nb.core("C12", "The turbulent-viscosity hypothesis $\\overline{u_iu_j}=\\frac23\\bar e\\delta_{ij}-\\nu_T\\big(\\frac{\\partial U_i}{\\partial x_j}+\\frac{\\partial U_j}{\\partial x_i}\\big)$ (12.94) and the
   mixing length", question="What does 'modelling the Reynolds stress' mean in practice, and how can one constant give the
   whole profile?")`
3. `nb.md` — **Plain words:** "Every ocean and atmosphere model has a number called the eddy viscosity (Ch. 13's Ekman layer
   cannot be written down without it). It is not a property of water or air — honey has a viscosity, a storm does not — but
   of the *flow*: big energetic eddies mix momentum fast. **N162 [B]**: turbulent-viscosity hypothesis for momentum,
   gradient-diffusion hypothesis for heat and scalars, by analogy with Newton's, Fourier's and Fick's laws of Ch. 1 and 4."
4. `nb.md` — **The idea:** "ν ~ (speed of the carriers) × (distance they travel before mixing). For molecules: sound speed ×
   mean free path. For turbulence: eddy speed $u_T$ × eddy size $l_T$. Near a wall an eddy cannot be bigger than its distance
   from the wall: $l_T=\kappa y$. One guessed length — and the log law comes out."
5. `nb.note` — **N163 [B]** $\overline{u_iT'}=-\kappa_T\,\partial\bar T/\partial x_i$ (12.95); **N164 [B]** $\overline{u_iY'}=-\kappa_{mT}\,\partial\bar Y/\partial x_i$ (12.96)
   (`gradient_diffusion_flux`); **N165 [B]** with $\overline{u_iu_j}=\frac23\bar e\delta_{ij}-\nu_T\big(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\big)$ (12.94) in the constant-density RANS equation,
   $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-\frac1\rho\frac{\partial P}{\partial x_i}+\frac{\partial}{\partial x_j}\Big([\nu+\nu_T]\Big(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\Big)-\frac23\bar e\delta_{ij}\Big)$ (12.97) — `> ⚠️ **slip #6 — the book prints**
   $\partial P/\partial x_j$ **; the correct form is** $\partial P/\partial x_i$` (i is the free index); the $\tfrac23\bar e$ term acts as an extra pressure. **N166 [B]**
   "Why the analogy is imperfect: molecules travel a mean free path far smaller than the flow (Kn ≪ 1, Ch. 1); eddies are as
   big as the shear layer itself, $l_T/L=O(1)$ — there is no scale separation, so ν_T cannot be a universal constant."
   **N167 [B]** $\nu_T,\ \kappa_T,\ \text{or}\ \kappa_{mT}\sim l_Tu_T$ (12.98) (`eddy_diffusivity_estimate`); "C16 will derive $D_T=\overline{u^2}\Lambda_t$, the
   same product".
6. `nb.derivation("D21", …)` — Part F D21 (11 steps), ref "12.101". Carries **N168 [B]**
   $0=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big(\nu\frac{\partial U}{\partial y}-\overline{uv}\big)=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big([\nu+\nu_T]\frac{\partial U}{\partial y}\big)$ (12.99), **N169 [B]** $-\overline{uv}=l_T^2(dU/dy)^2$, **N170 [B]**
   $0=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big(\nu\frac{dU}{dy}+\kappa^2y^2\big(\frac{dU}{dy}\big)^2\big)$ (12.100), **N171 [B]** the first integral and the exact slope, **N172 [B]**
   $\frac{dU}{dy}\cong\sqrt{\frac{\tau_0}\rho}\frac1{\kappa y}$, or $\frac U{u_*}\cong\frac1\kappa\ln y+\mathit{const.}$ (12.101).
7. `nb.worked_example("the model at y⁺ = 10", "κ = 0.41. 1. The first integral in wall units: s + κ²y⁺²s² = 1 with s = dU⁺/dy⁺.
   2. κ²y⁺² = 0.1681 × 100 = 16.81. 3. Positive root: s = 2/(1 + √(1 + 4 × 16.81)) = 2/(1 + 8.261) = 0.216. 4. So the viscous
   stress carries 21.6 % and the Reynolds stress 78.4 % of τ₀ here. 5. Compare the pure log slope 1/(κy⁺) = 0.244: already
   close. 6. Eddy viscosity ν_T/ν = κ²y⁺²s = 3.63.")`
8. `nb.code` — `ch12.mixing_length_wall_profile(10.0, 0.41)` → dict; the same with `damping="van_driest"`;
   `ch12.mixing_length_intercept(0.41)`, `ch12.mixing_length_intercept(0.41, A_plus=26.0)`. *expect:* slope 0.2160; U⁺(10) = 4.67
   (no damping) and 8.42 (A⁺ = 26); **B = −1.23 without damping (exactly [ln(4κ) − 1]/κ) and 5.28 with A⁺ = 26**. *explain:* "the
   plain model has the right slope 1/κ but its line sits 6 units too low: with $l_T=\kappa y$ the eddies mix right down to the
   wall and the sublayer is far too thin. Van Driest's damping, $l_T=\kappa y[1-e^{-y^+/A^+}]$ (glossed: an empirical factor that
   switches the eddies off within ~ A⁺ wall units), fixes it."
9. `nb.check_agree` — **from scratch:** the positive root of $s+\kappa^2y^{+2}s^2=1$ point by point on a log grid, then
   `cumulative_trapezoid`; `assert np.allclose(U_mine, [ch12.mixing_length_wall_profile(y, 0.41)["Uplus"] for y in grid],
   rtol=1e-4)`.
10. `nb.plotly` — **F6** `slider_figure` over A⁺ (0…40, 21 steps) with κ ∈ {0.384, 0.40, 0.41} as a dropdown: the model's
    U⁺(y⁺) on a semi-log axis, the reference log line (κ = 0.41, B = 5.0) and the implied B in the title. *expect:* B(0.41, 20)
    = 4.01, B(0.41, 26) = 5.28, B(0.41, 30) = 6.07; B(0.384, 26) = 5.10. *see / read* ("κ turns the line, A⁺ slides it") */ change*.
11. `nb.primer("a nonlinear diffusion problem by Picard iteration on the eddy viscosity", …)` (**P302**): "ν_T depends on the
    answer (dU/dy). Freeze ν_T from the last guess, solve the *linear* problem $\frac{d}{dy}([\nu+\nu_T]\frac{dU}{dy})=\frac1\rho\frac{dP}{dx}$, update ν_T, repeat
    until nothing changes (P192's idea applied to a boundary-value problem)." code: 6 lines on a 20-point grid, printing
    the change per sweep (falls by ~ 2 each time).
12. `nb.note` — **N181 [B]** (ours) the whole channel: the linear total stress of D17 plus the damped mixing length,
    $(1-y/\delta)u_*^2=\nu\frac{dU}{dy}+l_T^2\big(\frac{dU}{dy}\big)^2$, $l_T=\kappa y[1-e^{-y^+/A^+}]$ capped in the core. `ch12.channel_mixing_length(1000,
    kappa=0.41, A_plus=26.0)` + `nb.figure`: model U⁺ (purple) against public DNS dots (Lee & Moser 2015) and the laminar
    parabola at the same pressure gradient (muted; its centreline U⁺ = Re_τ/2 = 500 — "the turbulent flow is 20 times
    slower for the same push"); labelled **approximate**. *see / read / change*.
13. `nb.note` — **N173 [B]** a buoyant fluctuation accelerates as $Dw/Dt\sim g\alpha T'\sim g\Delta T/T$ (12.102) (α = 1/T for a perfect
    gas); **N174 [B]** free fall over the layer depth: $w^2/L\sim g\Delta T/T\Rightarrow w\sim\sqrt{gL\Delta T/T}$, $\kappa_T\sim wL$. Our layer: L = 1 km, ΔT = 2 K,
    T = 300 K ⇒ w ≈ 8.1 m/s, κ_T ≈ 8 × 10³ m²/s — 4 × 10⁸ times the molecular 2 × 10⁻⁵ m²/s (`convective_velocity_scale`,
    `convective_eddy_diffusivity`).
14. `nb.explainer("mixing_length_closure", heading="What does a closure constant do?", why="Sliding κ rotates the log line and
    sliding A⁺ shifts it, with the implied intercept B read off live; switching the model off returns the laminar parabola
    at the same pressure gradient.", tries=["Set A⁺ = 0: the line keeps its slope but B drops to −1.2.", "Set A⁺ = 26:
    B = 5.3.", "Change κ from 0.41 to 0.38 and watch the slope steepen.", "Switch the model off: the laminar parabola is 20
    times faster at Re_τ = 1000."])`
15. `nb.md` — **What would change if…** "…no single length like κy is available — a separated flow, a jet in a cross-wind? Then
    the model must carry its own velocity and length scales: C13."

#### C13 — The k–ε model (12.103)–(12.105)
1. `nb.core("C13", "The modelled energy equation $\\frac{\\partial\\bar e}{\\partial t}+U_j\\frac{\\partial\\bar e}{\\partial x_j}=\\frac{\\partial}{\\partial x_j}\\big(\\frac{\\nu_T}{\\sigma_e}\\frac{\\partial\\bar e}{\\partial x_j}\\big)-\\bar\\varepsilon-\\overline{u_iu_j}\\frac{\\partial U_i}{\\partial x_j}$ (12.103) —
   the k–ε model", question="Where do the constants of a turbulence model come from?")`
2. `nb.md` — **Plain words:** "Most engineering flow software, and the turbulence closures inside ocean and boundary-layer
   schemes, carry two extra fields: how much turbulent energy there is and how fast it is being dissipated. From those two
   they build an eddy viscosity everywhere. The price is five constants. Are they fudge factors? Two of them are pinned by
   experiments you can picture. (⚠️ the 'k' of k–ε is the book's $\bar e$; the book's k is a wavenumber or a conductivity.)"
3. `nb.md` — **The idea:** "velocity scale $u_T=\sqrt{\bar e}$; length scale = how far an eddy of that speed goes in its lifetime
   $\bar e/\bar\varepsilon$: $l_T=\bar e^{3/2}/\bar\varepsilon$; so $\nu_T\propto u_Tl_T=\bar e^2/\bar\varepsilon$."
4. `nb.note` — **N175 [B]** one-equation ingredients: $u_T=c\sqrt{\bar e}$, $\bar\varepsilon=C_\varepsilon\bar e^{3/2}/l_T$, and the three transport terms of the exact budget (C06) lumped into one gradient flux, $-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}=\frac{\nu_T}{\sigma_e}\frac{\partial\bar e}{\partial x_j}$ (slip #8 callout: the book prints $u_j$ in the
   viscous term) (`one_equation_closure`).
5. `nb.derivation("D22", …)` — Part F D22 (6 steps), ref "12.103": **N176 [B]** $\nu_T=C_\mu[\bar e^{3/2}/\bar\varepsilon]\sqrt{\bar e}=C_\mu\bar e^2/\bar\varepsilon$ (12.104).
6. `nb.note` — **N177 [B]** the modelled dissipation equation (built by analogy, not derived), term by term:
   $\frac{\partial\bar\varepsilon}{\partial t}+U_j\frac{\partial\bar\varepsilon}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\frac{\nu_T}{\sigma_\varepsilon}\frac{\partial\bar\varepsilon}{\partial x_j}\big)-C_{\varepsilon1}\big(\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}\big)\frac{\bar\varepsilon}{\bar e}-C_{\varepsilon2}\frac{\bar\varepsilon^2}{\bar e}$ (12.105) ("production of ε̄ ∝ production of ē ÷
   the time scale ē/ε̄; destruction ∝ ε̄ ÷ the same time scale"). **N178 [B]** `ch12.K_EPSILON_CONSTANTS` (Launder & Sharma 1974):
   $C_\mu=0.09$, $C_{\varepsilon1}=1.44$, $C_{\varepsilon2}=1.92$, $\sigma_e=1.0$, $\sigma_\varepsilon=1.3$ (slip #9 callout). **N179 [C]** the closed set = mean continuity $\partial U_i/\partial x_i=0$ (12.27), the mean momentum equation $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=\frac1\rho\frac{\partial\bar\tau_{ij}}{\partial x_j}$ (12.30, constant density), the closure $\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ (12.94), the ē equation of D22 with $\nu_T=C_\mu\bar e^2/\bar\varepsilon$ (12.104), and the ε̄ equation of N177 above; at a wall the log law of C11 is imposed as a boundary condition ("wall function");
   results are sensitive to inlet values; Reynolds-stress closures model the transport equation for $\overline{u_iu_j}$ (N59 in C04) instead (Ch. 10 pointer).
7. `nb.primer("dividing two ODEs to eliminate time", …)` (**P303**): "If $dy/dt=f$ and $dx/dt=g$, then $dy/dx=f/g$: time
   disappears and y is found as a function of x. If the result is $dy/dx=n\,y/x$, the solution is the power law $y\propto x^n$."
   code: 4 sympy lines: `dsolve(Eq(y(x).diff(x), 2*y(x)/x))` → C₁x².
8. `nb.derivation("D23", …)` — Part F D23 (12 steps), ref "12.105": **N180 [B]**.
9. `nb.worked_example("what the standard constants imply", "1. Decay exponent n = 1/(C_ε2 − 1) = 1/0.92 = 1.087 — grid-turbulence
   experiments give n ≈ 1.1–1.3: C_ε2 was chosen to match. 2. Log layer: κ² = √C_μ (C_ε2 − C_ε1) σ_ε = 0.3 × 0.48 × 1.3 = 0.187, κ =
   0.433 — a little above the measured 0.38–0.41: the constants are a compromise. 3. In the log layer ē/u_*² = 1/√C_μ = 3.33.
   4. Decay from ē₀ = 1 m²/s², ε̄₀ = 1 m²/s³: t₀ = n ē₀/ε̄₀ = 1.087 s; at t = 10 s, ē = (1 + 10/1.087)^{−1.087} = 0.080 m²/s².")`
10. `nb.code` — `ch12.k_epsilon_decay(1.0, 1.0, t)`; `ch12.k_epsilon_loglayer_kappa(0.09, 1.44, 1.92, 1.3)`;
    `ch12.k_epsilon_eddy_viscosity(1.0, 1.0)`. *expect:* n = 1.0870, t₀ = 1.0870; ē(1) = 0.4921, ε̄(1) = 0.2563; ē(10) = 0.0801,
    ε̄(10) = 0.00785; κ = 0.4327; ν_T(0) = 0.09, ν_T(10) = 0.0735 m²/s.
11. `nb.check_agree` — **from scratch:** an RK4 loop (P95) for $d\bar e/dt=-\bar\varepsilon$, $d\bar\varepsilon/dt=-C_{\varepsilon2}\bar\varepsilon^2/\bar e$ with Δt = 0.01 to t = 10;
    `assert np.allclose(e_rk4, e_closed, rtol=1e-6)`.
12. `nb.figure` — ē(t) and ε̄(t) on log–log axes: closed form (lines), `solve_ivp` (dots), slopes −n and −(n + 1) marked;
    second panel: ν_T(t) and $l_T=\bar e^{3/2}/\bar\varepsilon$ (grows slowly: the small eddies die first). *see / read* ("a straight line of
    slope −1.09: one constant, one measurable exponent") */ change* ("…C_ε2 = 2.0: n = 1.0 exactly — energy ∝ 1/t").
13. `nb.md` — **What would change if…** "…gravity acts on density fluctuations? The budget gains a buoyancy term that can
    feed the turbulence or starve it (C14)." (Backup explainer B1 `k_epsilon_decay` is storyboarded in Part B; it is built
    only if one of E1–E10 fails review. C13's visual is the figure of row 12.)

---

### A.11 §12.11 Turbulence in a Stratified Medium — R05 R06, C14, C15
#### C14 — The flux and gradient Richardson numbers (12.107)–(12.109)
1. `nb.section("12.11", "Turbulence in a Stratified Medium", intro="**What is this section about?** Turbulence where density
   changes with height — the atmosphere and the ocean. Buoyancy can feed turbulence (heated ground at noon) or drain it (a
   clear calm night). Two numbers measure the contest, and one length says at what height buoyancy takes over.")` ⚠️
   "Buoyancy is back (convention 6)."
2. `nb.recap("R05", "Stratification, potential temperature and the two lapse-rate conventions", "Stability is set by how the
   temperature gradient compares with the adiabatic one (Ch. 1 §1.10), so 'temperature' in this section means **potential**
   temperature θ: stable when dθ/dz > 0. In terms of the thermometer temperature, with Kundu's $\\Gamma\\equiv dT/dz$:
   $N^2=g\\alpha(dT/dz-\\Gamma_a)$, $\\Gamma_a=-g/C_p\\approx-9.8$ K/km, stable when $dT/dz>\\Gamma_a$; in meteorology's $\\Gamma_{met}\\equiv-dT/dz$:
   $N^2=g\\alpha(\\Gamma_d-\\Gamma_{met})$, $\\Gamma_d\\approx+9.8$ K/km, stable when $\\Gamma_{met}<\\Gamma_d$. Negate the number, flip the inequality (P48). Where these come from: Ch. 1's $N^2=-\\frac g{\\rho(z_o)}\\big(\\frac{d\\rho}{dz}-\\frac{d\\rho_a}{dz}\\big)$ (1.29) and $\\frac T\\theta\\frac{d\\theta}{dz}=\\frac{dT}{dz}+\\frac g{C_p}=\\Gamma-\\Gamma_a$ (1.32); Ch. 7's incompressible form $N^2\\equiv-\\frac g{\\rho_0}\\frac{d\\bar\\rho}{dz}$ (7.127).", where="Ch. 1 §1.10; Ch. 7; Ch. 11 §11.7")` + the two-row table and `STRAT.lapse_rate_stability` badge
   for three profiles (−6.5 K/km standard atmosphere: stable in both; −9.8: neutral; −12: unstable).
3. `nb.core("C14", "The flux Richardson number $\\mathrm{Rf}=\\frac{-g\\alpha\\overline{wT'}}{-\\overline{uw}(dU/dz)}$ (12.107), $\\mathrm{Ri}\\equiv\\frac{N^2}{(dU/dz)^2}$ (12.108) and
   $\\mathrm{Ri}=(\\nu_T/\\kappa_T)\\mathrm{Rf}$ (12.109)", question="When does stratification switch turbulence off?")`
4. `nb.md` — **Plain words:** "On a clear night the ground cools, the air near it becomes heavy, and the wind you felt at
   dusk dies near the surface while it still blows at tree-top height: the turbulence that carried momentum down has been
   suppressed. In the ocean the same contest decides how much heat is mixed below the thermocline. A mixing scheme must
   decide, level by level, whether the shear can keep turbulence alive against the stratification."
5. `nb.md` — **The idea:** "Lifting heavy fluid and pushing light fluid down costs energy; the turbulence pays it out of what
   the shear supplies. Rf = cost / income. If the cost takes more than about a quarter of the income, dissipation takes the
   rest and more — the turbulence starves." Table: | Rf < 0 | heat flux upward: buoyancy *adds* energy (convective) | 0 < Rf
   < ≈ ¼ | shear-driven, weakened | Rf ≳ ¼ | decaying |.
6. `nb.derivation("D24", …)` — Part F D24 (8 steps), ref "12.107". Carries **N182 [B]**
   $\frac{\partial\bar e}{\partial t}+U\frac{\partial\bar e}{\partial x}=-\frac{\partial}{\partial z}\big(\frac1{\rho_0}\overline{pw}+\frac12\overline{u_i^2w}\big)-\overline{uw}\frac{\partial U}{\partial z}+g\alpha\overline{wT'}-\bar\varepsilon$ (12.106, slip #10 callout) and **N184 [B]** $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109).
7. `nb.recap("R06", "The gradient Richardson number", "Ch. 11's $\\mathrm{Ri}\\equiv N^2/(dU/dz)^2$ (11.66), with linear stability guaranteed if, everywhere in the flow, $\\mathrm{Ri}>\\tfrac14$ (11.67) — here
   $\\mathrm{Ri}=\\alpha g(d\\bar T/dz)/(dU/dz)^2$ (12.108) with $\\bar T$ the potential temperature. New: `gradient_richardson_thermal` takes the
   thermometer gradient and $\\Gamma_a$ (Kundu sign) and returns Ri with the verdict in both conventions.", where="Ch. 11 §11.7")`
8. `nb.note` — **N183 [B]** "Turbulence stops being self-supporting near $\mathrm{Rf}_{cr}\approx0.25$ (an **observation**); large −Rf means
   convection dominates (`turbulence_regime`). `> ⚠️ Common confusion:` two different 'one quarter' statements — Ch. 11's
   theorem (Ri > ¼ everywhere ⇒ a laminar stratified shear flow is linearly stable: sufficient, exact) and this observed
   $\mathrm{Rf}_{cr}$ for existing turbulence, which by $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109) corresponds to $\mathrm{Ri}=\mathrm{Pr}_T/4$." **N185 [B]** turbulent Prandtl number
   $\nu_T/\kappa_T$: > 1 when stable (internal waves, Ch. 7, carry momentum but not heat), small when unstable, ≈ 1 when neutral
   (the *Reynolds analogy*). ⚠️ "ν is set like an italic v in the book's $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109)."
9. `nb.worked_example("a stable night", "$\\overline{wT'}$ = −0.02 K m/s (downward), $\\overline{uw}$ = −0.09 m²/s² (u_* = 0.3 m/s), dU/dz = 0.1 s⁻¹, α =
   1/300 K⁻¹, g ≈ 9.81. 1. Buoyant destruction −gα$\\overline{wT'}$ = 9.81 × 0.02/300 = 6.54 × 10⁻⁴ m²/s³. 2. Shear production −$\\overline{uw}$ dU/dz =
   0.09 × 0.1 = 9.0 × 10⁻³. 3. Rf = 0.073 — shear-driven. 4. Thermometer gradient dT/dz = +10 K/km (an inversion). Kundu:
   dθ/dz = dT/dz − Γ_a = 0.010 − (−0.0098) = 0.0198 K/m. Meteorology: Γ_met = −10 K/km, Γ_d − Γ_met = 9.8 + 10 = 19.8 K/km — the
   same number. 5. N² = gα dθ/dz = 6.47 × 10⁻⁴ s⁻²; Ri = N²/(dU/dz)² = 0.065. 6. ν_T = 0.09/0.1 = 0.9 m²/s; κ_T = 0.02/0.0198 = 1.01 m²/s;
   Pr_T = 0.89; and Ri/Rf = 0.065/0.073 = 0.89 ✓ — $\\mathrm{Ri}=(\\nu_T/\\kappa_T)\\mathrm{Rf}$ (12.109).")`
10. `nb.code` — `ch12.flux_richardson(-0.02, -0.09, 0.1, 1/300)`; `ch12.gradient_richardson_thermal(0.010, 0.1, 1/300,
    Gamma_a=-0.0098)`; `ch12.turbulent_prandtl(0.9, 1.0101)`; `ch12.turbulence_regime(0.0727)`. *expect:* 0.0727; Ri = 0.0647,
    verdict_kundu "stable ⇔ dT/dz > Γa: 10.0 > −9.8 K/km", verdict_met "stable ⇔ Γ < Γa: −10.0 < 9.8 K/km" (the `.text` of the existing `STRAT.lapse_rate_stability(0.010, Gamma_a=-0.0098, convention="kundu" / "meteorology")`, run for this design; in the meteorological line the function's "Γa" is the positive dry-adiabatic rate we call Γ_d); 0.891; "shear-driven".
11. `nb.check_agree` — **from scratch:** Rf and Ri by hand (Kundu sign), then Ri again from dθ/dz directly; `assert np.isclose`
    for both routes and with the functions; `assert np.isclose(Ri/Rf, 0.891, rtol=1e-3)`.
12. `nb.figure` — budget bars for three cases (unstable noon H = +150 W/m²; neutral; stable night H = −30 W/m²) at z = 10 m:
    shear production (orange), buoyancy (blue, ± ), dissipation as the residual (rose); Rf printed on each; legend lines in
    both conventions. *see / read / change* ("…the wind dropped by half at night: production ∝ u_*³ falls 8-fold, Rf rises
    8-fold to 0.96 — the turbulence collapses, which is why calm clear nights are so still").
13. `nb.md` — **What would change if…** "…we ask at what *height* buoyancy catches up with shear? In the log layer the answer
    is a single length (C15)."

#### C15 — The Monin–Obukhov length (12.110) and the stratified surface layer
1. `nb.core("C15", "The Monin–Obukhov length $L_M\\equiv-u_*^3/(\\kappa\\alpha g\\overline{wT'})$ (12.110), with $\\mathrm{Rf}=z/L_M$ (12.111)",
   question="What does the Monin–Obukhov length measure?")`
2. `nb.md` — **Plain words:** "A weather or climate model knows the wind and temperature at its lowest level, some tens of
   metres up, and needs the stress and heat flux at the ground. The bridge is a wind profile that is logarithmic near the
   ground and bends with stability — and the one length that says how fast it bends is $L_M$. Every bulk surface-flux
   formula is built on it."
3. `nb.md` — **The idea:** shear production in the log layer falls with height as $u_*^3/(\kappa z)$; the buoyancy term
   $g\alpha\overline{wT'}$ does not change with height in the surface layer. They are equal in size at $z=\lvert L_M\rvert$. `surface_layer_sketch`
   (**N187 [B]**, our Fig. 12.21): unstable day — forced convection below $\lvert L_M\rvert$, free convection with thermal plumes above.
4. `nb.derivation("D25", …)` — Part F D25 (7 steps), ref "12.111". Carries **N186 [B]** $\mathrm{Rf}=z/L_M$ (12.111) and **N188 [B]**
   $U=\frac{u_*}\kappa\big[\ln\frac z{z_o}+5\frac z{L_M}\big]$ (the book's coefficient 5; Businger–Dyer's stable-side value is 4.7 — glossed, AMS Glossary).
5. `nb.primer("the stability parameter ζ = z/L and integrating a flux–profile relation", …)` (**P304**) — placed **before**
   D25: "Surface-layer schemes write the dimensionless shear $\phi_m=(\kappa z/u_*)\,dU/dz$ as a function of ζ = z/L alone: 1 when
   neutral, > 1 stable, < 1 unstable. Integrating $dU/dz=(u_*/\kappa z)\phi_m$ from $z_0$ gives the wind profile." code: 4 lines integrating
   φ_m = 1 + 5ζ numerically from z₀ = 0.03 m to 10 m with L = 83 m and comparing with ln(z/z₀) + 5(z − z₀)/L.
6. `nb.worked_example("a clear night over grass", "u_* = 0.3 m/s, H = −30 W/m² (downward), ρ = 1.2 kg/m³, c_p = 1005 J/(kg K), T = 300 K,
   κ = 0.4, z₀ = 0.03 m. 1. $\\overline{wT'}$ = H/(ρc_p) = −30/1206 = −0.0249 K m/s (P296). 2. α = 1/T = 3.33 × 10⁻³ K⁻¹. 3. L_M = −u_*³/(κ α g
   $\\overline{wT'}$) = −0.027/(0.4 × 3.33 × 10⁻³ × 9.81 × (−0.0249)) = +83 m (positive: stable). 4. Rf at 10 m = z/L_M = 0.12. 5. Rf reaches ¼ at
   z = L_M/4 = 21 m: above that the turbulence cannot sustain itself. 6. Wind at 10 m: (0.3/0.4)[ln(333) + 5 × 10/83] = 0.75 ×
   (5.81 + 0.60) = 4.81 m/s, against 4.36 m/s for a neutral profile — more shear for the same stress.")`
7. `nb.code` — `L = ch12.monin_obukhov_from_fluxes(tau=1.2*0.3**2, H=-30.0, rho=1.2, cp=1005.0, T=300.0, kappa=0.4)`;
   `ch12.flux_richardson_surface_layer(10.0, L)`; `WT.surface_layer_wind(10.0, 0.3, 0.03, L, kappa=0.4)`; `st =
   ch12.surface_layer_state(0.3, -30.0, 300.0, 0.03, 10.0, 1.2, 1005.0, 0.4)`. *expect:* L = 83.0 m (83.01 with G0, 82.98 with the hand value g = 9.81); Rf = 0.120; U = 4.81 (neutral 4.357); z_crit = 20.7 m; for H = +150: L = −16.6 m, Rf(10 m) = −0.60 (free convection above ~ 17 m).
8. `nb.check_agree` — **from scratch:** L_M from u_* and H in one line; the wind at 10 m from the log-linear formula;
   `assert np.isclose` with `monin_obukhov_from_fluxes` and `surface_layer_wind`.
9. `nb.plotly` — **F7** (**N189 [B]**, our Fig. 12.22) `slider_figure` over 1/L_M (from −0.1 to +0.05 m⁻¹, 31 steps through 0):
   ln z against U (stable: log-linear; unstable: Businger–Dyer; the book's log-linear form dashed where it leaves its
   range), neutral straight line as ghost, Rf(z) on a second panel with the ¼ line; trace names carry the regime and the
   verdict in both conventions. + `nb.live` (free u_*, H, z₀). *see / read* ("stable: the profile bends toward larger U —
   more shear; unstable: less shear, the air is well mixed") */ change*.
10. `nb.note` — the temperature field: **N190 [B]**
    $\frac{\partial}{\partial t}\big(\frac12\overline{T'^2}\big)+U\frac{\partial}{\partial x}\big(\frac12\overline{T'^2}\big)=-\overline{wT'}\frac{d\bar T}{dz}-\frac{\partial}{\partial z}\big(\frac12\overline{T'^2w}-\kappa\frac{\partial}{\partial z}\big(\frac12\overline{T'^2}\big)\big)-\bar\varepsilon_T$, $\bar\varepsilon_T=\kappa\overline{(\partial T'/\partial x_j)^2}$ (12.112, slip #15
    corrected: the page prints $\kappa\,\partial\overline{T'^2}/\partial z$; here κ is the thermal diffusivity) — "the recipe of D10 with T′ for $u_i$"
    (`temperature_variance_sympy` shows the ½). **N191 [B]** $\overline{T'^2}\equiv\int_0^\infty S_T(K)\,dK$. **N192 [B]** Obukhov–Corrsin:
    $S_T\propto\bar\varepsilon_T\bar\varepsilon^{-1/3}K^{-5/3}$ for $2\pi/L\ll K\ll2\pi/\eta$ (12.113) — D12's argument with temperature as an extra dimension
    (`DIM.solve_exponents` in the cell). **N193 [B]** Batchelor scale $\eta_T=\eta(\kappa/\nu)^{1/2}$ for ν/κ ≫ 1 (sea water, Pr ≈ 7, η = 1 mm:
    η_T = 0.38 mm). **N194 [B]** $S_T\propto K^{-1}$ for $2\pi/\eta\ll K\ll2\pi/\eta_T$ (12.114). **N195 [B]** + `nb.figure` (our Fig. 12.23):
    `model_spectrum` and `scalar_spectrum` for sea water — both −5/3, then temperature flattens to −1 while velocity rolls
    off. *see / read / change*.
11. `nb.explainer("stratified_surface_layer", heading="When does stratification kill turbulence?", why="Dragging the surface
    heat flux through zero flips the sign of L_M, bends the wind profile to the other side of the neutral logarithm and
    moves the Rf = ¼ height up and down, while the badge states the regime in both lapse-rate conventions.", tries=["Load
    'clear calm night' and lower u_*: the Rf = ¼ height drops toward the ground.", "Drag H through zero: L_M jumps from
    +∞ to −∞ and the profile crosses the neutral line.", "Toggle the convention: the same layer, the other inequality.",
    "Switch z₀ from sea to forest at the same 10 m wind."])`
12. `nb.md` — **What would change if…** "…instead of momentum and heat we follow a puff of smoke released into the turbulence?
    C16."

---

### A.12 §12.12 Taylor's Theory of Turbulent Dispersion — R07 R08, C16
#### C16 — Taylor's dispersion formula (12.119)
1. `nb.section("12.12", "Taylor's Theory of Turbulent Dispersion", intro="**What is this section about?** How far turbulence
   carries a marked particle. One exact formula shows the cloud spreads in proportion to time at first and to the square
   root of time later — and that an 'eddy diffusivity' is not a constant.")`
2. `nb.core("C16", "Taylor's formula $\\overline{X_\\alpha^2}(t)=2\\overline{u_\\alpha^2}\\,t\\int_0^t\\big(1-\\frac\\tau t\\big)r_\\alpha(\\tau)\\,d\\tau$ (12.119)", question="Why does a puff
   spread like t at first and like √t later?")`
3. `nb.md` — **Plain words:** "Smoke leaves a chimney as a narrow cone that opens like a wedge, then, far downwind, widens
   ever more slowly. An oil slick, volcanic ash, a cloud of drifting buoys in the ocean behave the same way. The switch
   happens when each smoke particle has forgotten the velocity it started with."
4. `nb.md` — **The idea:**
   ```
   t ≪ Λ_t : the particle still has its first velocity  → flies straight:  X ≈ u t          → X_rms = u_rms · t
   t ≫ Λ_t : many independent "steps" of duration ~2Λ_t  → random walk:    X_rms = step·√n  → X_rms = u_rms √(2Λ_t t)
   ```
5. `nb.note` — **N196 [B]** Lagrangian position $\mathbf X(\mathbf a,t)$ of the particle that started at **a** (Ch. 3); stationary homogeneous
   turbulence, zero mean velocity; `ch12.langevin_particles` (the OU velocity of P283 integrated in time) — figure of
   particle paths (our Fig. 12.24). **N199 [B]** $r_\alpha(\tau)\equiv\overline{u_\alpha(t)u_\alpha(t+\tau)}/\overline{u_\alpha^2}$ — ⚠️ "α is a component label (no sum), not the
   expansion coefficient"; here $\Lambda_t$ is the **Lagrangian** integral scale.
6. `nb.primer("a double integral over a triangle", …)` (**P305**): "$\int_0^tdt'\int_0^{t'}g(\tau)\,d\tau$ covers the triangle 0 < τ < t′ < t. Each
   value of τ is counted for every t′ between τ and t, a strip of length t − τ, so the double integral equals
   $\int_0^t(t-\tau)g(\tau)\,d\tau$." code: 4 lines with g = e^{−τ}, t = 2, both forms by `quad`/`dblquad` (1.1353).
7. `nb.derivation("D26", …)` — Part F D26 (11 steps), ref "12.119". Carries **N197 [B]** $\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{X_\alpha\frac{dX_\alpha}{dt}}$ (12.115),
   **N198 [B]** $\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{X_\alpha u_\alpha}=2\int_0^t\overline{u_\alpha(t')u_\alpha(t)}\,dt'$ (12.116), **N200 [B]** $\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{u_\alpha^2}\int_0^tr_\alpha(\tau)\,d\tau$ (12.117), **N201 [B]** $\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\int_0^tdt'\int_0^{t'}r_\alpha(\tau)\,d\tau$ (12.118), **N202 [B]** the integration
   by parts.
8. `nb.derivation("D27", …)` — Part F D27 (8 steps), ref "12.121". Carries **N204 [B]** $\overline{X_\alpha^2}\simeq\overline{u_\alpha^2}t^2$ (12.120), **N205 [B]**
   $(X_\alpha)_{rms}=(u_\alpha)_{rms}t$ for $t\ll\Lambda_t$ (12.121), **N206 [B]** $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\Lambda_tt$ (12.122), **N207 [B]** $(X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}$ for
   $t\gg\Lambda_t$ (12.123) (slip #11 callout), **N208 [B]** the closed form for $r=e^{-\tau/\Lambda_t}$; the Gaussian case
   ($r=e^{-\tau^2/t_c^2}$, $\Lambda_t=\tfrac{\sqrt\pi}2t_c$, Exercise 12.38) stated with `taylor_dispersion_gaussian`. **N203 [B]** figure (our Fig.
   12.25): C02's correlation curve with "small t" and "large t" marked.
9. `nb.worked_example("a puff in turbulence with a 10 s memory", "u_rms = 1 m/s, Λ_t = 10 s, exponential correlation: $\\overline{X^2}=2\\overline{u^2}\\Lambda_t^2[t/\\Lambda_t-1+e^{-t/\\Lambda_t}]$.
   1. t = 1 s: 200 × (0.1 − 1 + 0.90484) = 0.967 m² — ballistic value u²t² = 1.000 (3 % high): X_rms ≈ 0.98 m. 2. t = 10 s: 200 ×
   0.3679 = 73.6 m² (X_rms = 8.6 m) — neither limit works (ballistic 100, diffusive 200). 3. t = 100 s: 200 × 9.000 = 1800 m²
   (X_rms = 42 m); diffusive formula 2u²Λ_t t = 2000 — 11 % high; with the constant offset −2u²Λ_t² = −200 it is exact. 4. Eddy
   diffusivity at t = 10 s: D_T = u²Λ_t(1 − e⁻¹) = 6.3 m²/s; final value 10 m²/s.")`
10. `nb.code` — `X, u = ch12.langevin_particles(10000 if not FAST else 2000, t, 1.0, 10.0, seed=0)`; `X2 = (X**2).mean(axis=0)`;
    `ch12.taylor_dispersion_exponential(t, 1.0, 10.0)`; `ch12.taylor_dispersion(t, lambda s: np.exp(-s/10), 1.0, form="double")`
    and `form="single"`; `ch12.dispersion_rate_from_particles(t, X, u)`; `ch12.dispersion_regime(1.0, 10.0)`. *expect:* closed form
    0.9675, 73.58, 1800.0 at t = 1, 10, 100 s; particles within 5 standard errors; both integral forms equal to 1e-8; both
    sides of $\frac{d}{dt}(\overline{X_\alpha^2})=2\overline{X_\alpha u_\alpha}$ (12.116) agree within sampling error.
11. `nb.check_agree` — **from scratch:** Taylor's integral $2\overline{u^2}t\int_0^t(1-\tau/t)e^{-\tau/\Lambda_t}d\tau$ by `np.trapezoid`; `assert np.allclose(mine,
    ch12.taylor_dispersion_exponential(t, 1.0, 10.0), rtol=1e-5)`; a random-walk loop (cumulative sum of unit steps in random
    directions) giving $\overline{R_n^2}=nL^2$ within 5 standard errors of `ch12.random_walk`.
12. `nb.animation` — **A5** (video, 90 frames, FAST 40): 400 Langevin particles leaving a point source with the ±X_rms
    envelope and the two asymptotes (dashed); beside it $\overline{X^2}(t)$ on log–log axes (dots: particles; line: closed form;
    slopes 2 and 1). *see / read / change* ("…Λ_t = 100 s: the cloud would still be a wedge at the end of the movie").
13. `nb.primer("proof by induction", …)` (**P306**): "Show a statement for n = 1; show that if it holds for n − 1 it holds for
    n; then it holds for every n." code: 3 lines checking $\sum_{k=1}^nk=n(n+1)/2$ for n = 1…5.
14. `nb.note` — the random walk: **N209 [B]** from $\mathbf R_n=\mathbf R_{n-1}+\mathbf L$: $\overline{R_n^2}=\overline{R_{n-1}^2}+L^2+2\overline{\mathbf R_{n-1}\cdot\mathbf L}$ (12.124) — the cross term
    averages to zero when the new step is uncorrelated with the past (the product rule of C01 again); **N210 [B]** by
    induction $\overline{R_n^2}=\overline{R_{n-1}^2}+L^2=\dots=nL^2$ (figure: three walker paths and the ensemble rms, our Fig. 12.26); **N211 [B]**
    $(R_n)_{rms}=L\sqrt n$ (12.125) — "with $n=t/\Delta t$, $L=u_{rms}\Delta t$ and $\Delta t=2\Lambda_t$ this is $(X_\alpha)_{rms}=u_{rms}\sqrt{2\Lambda_tt}$ (12.123): a particle makes
    one independent step every two memory times". *expect:* for 4000 walkers, n = 100: rms = 10.0 ± 0.2.
15. `nb.note` — **N212 [B]** + `nb.figure` (our Fig. 12.27): a time-averaged plume in a wind U = 5 m/s with w_rms = 0.5 m/s, Λ_t =
    20 s: `ch12.smoke_plume_width` (t = x/U) and `plume_concentration` as colour; *expect:* Z_rms = 0.98, 8.6, 42.4, 99.0 m at
    x = 10, 100, 1000, 5000 m (linear near, √x far). `> ⚠️ **slip #13 — the book's caption prints** width ∝ x^{1/2} near the
    source and ∝ x far away **; the correct statement is** ∝ x near ($(X_\alpha)_{rms}=(u_\alpha)_{rms}t$ (12.121)) and ∝ x^{1/2} far
    ($(X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}$ (12.123))`. *see / read / change*.
16. `nb.recap("R07", "The molecular yardstick", "A line vortex switched on at t = 0 spreads as $u_\\theta=(\\Gamma/2\\pi r)\\exp(-r^2/4\\nu t)$,
    a Gaussian of standard deviation $\\sigma=\\sqrt{2\\nu t}$ per coordinate. ⚠️ Ch. 3's core radius used $\\sigma^2=4\\nu t$.", where="Ch. 8
    (Exercise 8.26 spin-up); Ch. 3, Ch. 5 Lamb–Oseen")`
17. `nb.recap("R08", "A diffusivity from the growth of a variance", "Ch. 1's $\\sigma^2=2Dt$ read backwards: $\\nu=\\frac12\\frac{d\\sigma^2}{dt}$ (12.126)
    (`diffusivity_from_variance`). This is the starting line of D28.", where="Ch. 1 §1.5; `core.diffusion.gaussian_spreading`")`
18. `nb.derivation("D28", …)` — Part F D28 (7 steps), ref "12.127". Carries **N213 [B]**
    $D_T\equiv\frac12\frac{d}{dt}(\overline{X_\alpha^2})=\overline{u_\alpha^2}\int_0^tr_\alpha(\tau)\,d\tau$ (12.127), **N214 [B]** $D_T\cong\overline{u_\alpha^2}t$ for $t\ll\Lambda_t$ (12.128), **N215 [B]**
    $D_T\cong\overline{u_\alpha^2}\Lambda_t$ for $t\gg\Lambda_t$ (12.129) with `> ⚠️ **slip #12 — the book prints** the condition $t\ll\Lambda_t$ on $D_T\cong\overline{u_\alpha^2}\Lambda_t$ (12.129) **; the
    correct condition is** $t\gg\Lambda_t$`.
19. `nb.plotly` — **F8** `slider_figure` over Λ_t (1…100 s, 25 steps): X_rms(t) on log–log axes with both limits and a
    constant-D Gaussian ghost of the same late-time width; D_T(t) on a second panel rising to $\overline{u^2}\Lambda_t$. *see / read / change*.
20. `nb.note` — **N216 [B]** "Turbulent diffusion is not molecular diffusion with a bigger constant: near the source a
    constant-D model spreads the puff as √t — far too wide a start, then too slow; and a growing patch is stirred by ever
    larger eddies (Ch. 10's distinction between a numerical and an eddy diffusivity recalled)." **N06 [B]** Richardson's
    four-thirds law, stated here where both ingredients exist: an effective diffusivity for a patch of size l can depend
    only on l and ε̄ in the inertial range (C08), so $K\sim\bar\varepsilon^{1/3}l^{4/3}$ (`DIM.solve_exponents` → (1/3, 4/3);
    `ch12.richardson_diffusivity(1000.0, 1e-8)` ≈ 22 m²/s for a 1 km ocean patch with ε̄ = 10⁻⁸ m²/s³; 10³ m²/s for the
    atmosphere with ε̄ = 10⁻³). **Climate hook:** why tracer diffusivities in ocean models depend on the grid size.
21. `nb.explainer("taylor_dispersion", heading="Why t first and √t later?", why="The particles spread on one clock with the
    ⟨X²⟩ curve changing slope from 2 to 1 as t passes Λ_t; dragging Λ_t moves the bend; a constant-diffusivity ghost shows
    how wrong Fickian diffusion is near the source; the plume mode maps t to x/U.", tries=["Press ▶ and watch the local slope
    readout fall from 2 to 1.", "Drag Λ_t up tenfold: the bend moves tenfold later and the final D_T is ten times
    larger.", "Switch to 'plume in a wind': the same curve drawn in space.", "Compare with the constant-D ghost at t = 0.1 Λ_t."])`
22. `nb.md` — **What would change if…** "…the planet rotates and the fluid is stratified on the large scale? The eddies become
    flat, nearly two-dimensional, and their energy moves to *larger* scales instead of smaller ones — C07's cascade turned
    round. The tools of this chapter (averages, Reynolds stresses, eddy viscosity, Ri and $L_M$, dispersion) all carry over;
    the next chapter adds the Coriolis force to them."

### A.13 §12.13 Concluding Remarks — C16 (continued): N217; S01, S02; summary
1. `nb.section("12.13", "Concluding Remarks", intro="**What is this section about?** Where the subject goes from here.")`
2. `nb.note` — **N217 [C]** "Turbulence remains an open research field. What this chapter left out comes next: Ch. 13 adds
   rotation and large-scale stratification — eddy-viscosity Ekman layers (C12), surface layers (C11, C15), and geostrophic
   turbulence, whose energy flows toward *larger* scales with a −3 enstrophy range instead of C08's −5/3. Direct and
   large-eddy simulation are named in Ch. 10."
3. `nb.pointer("**S01** — Exercises 12.1–12.38 are not reproduced. Results this lesson derived or stated: the moment identities (N22), the real, even spectrum (D05), the periodogram (N38), the Taylor microscale (D04), frozen turbulence (N40), the mean scalar equation (N57), the mean-flow energy budget (D09), the Reynolds-stress equation (N59), the isotropic tensor and dissipation (D07, D08), the exponents of Table 12.1 (D16), laminar vs turbulent skin friction (N161), the channel force balance (D17), the pipe friction law (N160), dispersion for a Gaussian correlation (N208).")`
4. `nb.pointer("**S02** — Literature: public sources are cited where used (Sreenivasan 1995; Lee & Moser 2015; Launder & Sharma 1974; McKeon et al. 2005 for Prandtl's law; AMS Glossary for Businger–Dyer; Pao 1965; Spalding 1961; Coles 1956; Taylor 1921). For more: Pope, *Turbulent Flows*; Tennekes & Lumley, *A First Course in Turbulence*.")`
5. `nb.summary(clicked=[…16…], feeds_forward=[…], left_out=[…])` — **clicked** (one per CORE): C01 "An average passes through
   every linear operation and fails on a product: $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$." · C02 "The area under the autocorrelation is the
   memory time; a record is worth Δt/Λ_t independent samples." · C03 "The spectrum is the Fourier transform of the
   autocorrelation: long memory, narrow spectrum; area = variance." · C04 "Averaging the momentum equation leaves one new
   term, the Reynolds stress $-\rho_0\overline{u_iu_j}$ — a covariance, negative $\overline{uv}$ in a positive shear, with no equation of its own." ·
   C05 "With no preferred direction, one function f(r) fixes every two-point correlation and $\bar\varepsilon=15\nu\overline{u^2}/\lambda_g^2$." · C06
   "Shear production leaves the mean-flow budget and enters the turbulent one; dissipation ends it." · C07 "The big eddies
   set ε̄ ~ ΔU³/L; viscosity only sets η = (ν³/ε̄)^{1/4}, and η/L ~ Re^{−3/4}." · C08 "Between L and η only ε̄ and k matter, so
   $S_{11}\propto\bar\varepsilon^{2/3}k_1^{-5/3}$." · C09 "A conserved momentum flux plus shape-keeping give δ ∝ x and U_CL ∝ x^{−1/2}." · C10 "Near a
   smooth wall τ₀ and ν alone set the scales: U⁺ = f(y⁺), and U⁺ = y⁺ in the sublayer." · C11 "Where the answer may depend on
   neither the viscous nor the outer length, y dU/dy = u_*/κ: a logarithm." · C12 "An eddy viscosity is a length times a
   velocity of the flow; l_T = κy returns the log law, wall damping sets its intercept." · C13 "k–ε builds ν_T = C_μē²/ε̄ from
   two transport equations; its constants are tied to the decay exponent and to κ." · C14 "Rf = buoyant destruction / shear
   production; turbulence starves near Rf ≈ ¼; Ri = Pr_T Rf." · C15 "L_M is the height where buoyancy matches shear: Rf = z/L_M."
   · C16 "A cloud spreads as t while particles remember their velocity and as √t afterwards; D_T grows to $\overline{u^2}\Lambda_t$."
   **feeds_forward:** Ch. 13 (eddy viscosity in Ekman layers; surface layer and bulk drag; Ri-dependent mixing; geostrophic
   turbulence spectra; `TS.periodogram` for wave and eddy spectra) · Ch. 14 (turbulent skin friction on airfoils) · Ch. 15
   (friction in ducts). **left_out:** DNS and LES (Ch. 10 pointer; Pope) · intermittency and structure functions (Frisch) ·
   Reynolds-stress closures (Pope ch. 11) · two-dimensional turbulence (Ch. 13).

### A.14 Placement check (every curation id has exactly one home)
- **CORE (16):** C01 A.3 · C02, C03 A.4 · C04 A.5 · C05 A.6 · C06, C07, C08 A.7 · C09 A.8 · C10, C11 A.9 · C12, C13 A.10 · C14,
  C15 A.11 · C16 A.12.
- **RECAP (8):** R01, R02 in C04 · R03, R04 in C10 · R05 before C14 (its A parent; placed as the section's opening recap), R06
  in C14 · R07, R08 in C16.
- **SKIP (2):** S01, S02 in A.13.
- **NOTE (217):** N01–N04 A.1 · N05 A.2 · N06 named in A.2 and stated in C16 row 20 · N07–N23 C01 · N24–N34 C02 · N35–N40 C03 ·
  N41–N61 C04 · N62–N77 C05 · N78–N83 C06 · N84–N92 C07 · N93–N99 C08 · N100–N130 C09 · N131–N142 C10 · N143–N161 C11 ·
  N162–N174 and N181 C12 · N175–N180 C13 · N182–N185 C14 · N186–N195 C15 · N196–N216 C16 · N217 A.13.
- **Derivations (28):** D01 C01 · D02 D03 D04 C02 · D05 C03 · D06 C04 · D07 D08 C05 · D09 D10 C06 · D11 C07 · D12 C08 · D13 D14
  D15 D16 C09 · D17 D18 C10 · D19 D20 C11 · D21 C12 · D22 D23 C13 · D24 C14 · D25 C15 · D26 D27 D28 C16.
- **Primers (27):** P280–P285 C01 (P286 in A.1, which opens C01) · P287–P289 C02 · P290–P293 C03 · P294–P296 C04 · P297 C05 ·
  P298 C08 · P299, P300 C10 · P301 C11 · P302 C12 · P303 C13 · P304 C15 · P305, P306 C16.
- **Explainers (10, each embedded once):** E1 C01 · E2 C03 · E3 C04 · E5 C06 · E4 C08 · E6 C09 · E7 C11 · E8 C12 · E9 C15 · E10 C16.
  C05 and C13 have no explainer (static figures; B1 is the backup for C13). C02's visuals are its figure and A2; C07's is A4;
  C10's its two figures; C14's the budget-bar figure.
- **Animations:** A1 C01 · A2 C02 · A3 C04 · A4 C07 · A5 C16. **Slider figures:** F1 C01 · F2 C03 · F3 C08 · F4 C09 · F5 C11 · F6 C12 ·
  F7 C15 · F8 C16. **Live cells:** with F3 and F7.
- **Coded-note calls:** the 39 rows of A.15a place every remaining Part C function in the note the curation names for it.

### A.15 Addendum to the blocks above (binding; added by the audit of 2026-10-07)

**A.15a — function calls the curation attaches to notes.** The curation's chapter map states these notes "with" a named
function; the rows above give the note's text and equation but not the call. The builder adds each call to the cell of
that note (one commented line, then the printed value in the note's sentence). *expect* values were computed on 2026-10-07
with the implementer's module as it stood (closed forms checked by hand where marked ✎); where the cell needs arrays the
intent is given and the builder reads the docstring for the layout. Constants are our own or labelled illustrative
(convention 10).

| Block · note (row) | Call (exact) | *expect* and the sentence it supports |
|---|---|---|
| C04 · N45 (row 8, under D06) | `ch12.mean_divergence(ens_u, ens_v, dx, dy)` on 16 members of `TS.synthetic_solenoidal_field` | the divergence of the ensemble mean is round-off (`np.abs(...).max()` ≲ 1e-10 of u_rms/dx): "every member is divergence-free, so the mean is, $\partial U_i/\partial x_i=0$ (12.27), and so is each fluctuation" |
| C04 · N54 (row 18) | `ch12.turbulent_heat_flux(0.5, 0.4, 0.5, 1.2, 1005.0)` | 120.6 W/m² ✎ (= ρ c_p r w_rms T_rms = 1.2 × 1005 × 0.5 × 0.5 × 0.4) — the same flux as primer P296's 0.1 K m/s |
| C04 · N55 (row 20) | `ch12.mixture_density(0.1, 1.8, 1.2)`; `ch12.mass_fraction_from_volume_fraction(0.1, 1.8, 1.2)` | 1.26 kg/m³; 0.1429 ✎ (= 0.18/1.26): volume and mass fractions differ when the densities do |
| C04 · N57 (row 20) | `ch12.mean_scalar_flux(-0.02, 0.001, 2e-5)` | 1.0004 × 10⁻³ m/s ✎: molecular part 4 × 10⁻⁷, turbulent part 10⁻³ — the turbulent flux is 2500 times larger |
| C04 · N59 (the row of `reynolds_stress_budget_sympy`) | `ch12.reynolds_stress_production(uu, gradU)` with $\overline{u_iu_j}$ = [[0.5, −0.1, 0], [−0.1, 0.2, 0], [0, 0, 0.3]] m²/s² and dU/dy = 2 s⁻¹ | P₁₁ = 0.4, P₁₂ = −0.4 m²/s³, all else 0 ✎; half its trace, 0.2, equals `ch12.shear_production(uu, gradU)` — shear feeds only the streamwise normal stress |
| C05 · N64 (row 5) | `K, E = TS.shell_spectrum((u, v, w), L)` on the cached 3-D field | `np.sum(E)*(K[1]-K[0])` equals `0.5*(u**2+v**2+w**2).mean()` to round-off: the field has the spectrum we prescribed |
| C05 · N72 (row 6) | `TS.turbulent_kinetic_energy(u_samples)` for 10⁵ samples of three independent unit-variance components (seed 0) | 1.50 ± 0.01 = $\tfrac32\overline{u^2}$ |
| C05 · N68 (row 8, under D07) | `ch12.isotropic_correlation_tensor(np.array([0.5, 0.0, 0.0]), lambda r: np.exp(-r**2))` | diag(0.7788, 0.5841, 0.5841) ✎: R₁₁ = f(0.5) = e^{−0.25}; R₂₂ = R₃₃ = g(0.5) = (1 − r²)e^{−r²}; off-diagonals 0 |
| C05 · N69 (row 8, under D07) | `out = ch12.isotropic_tensor_divergence_sympy()` | `out["divergence"] == [0, 0, 0]` and `out["check"] is True`: $R_{ij}=\overline{u^2}\{f\delta_{ij}+\frac r2f'(\delta_{ij}-r_ir_j/r^2)\}$ (12.41) is divergence-free for any f |
| C08 · N95 (row 6) | `ch12.one_dimensional_from_3d(10.0, lambda K: 1.5*K**(-5/3)) / (1.5*10.0**(-5/3))`; `ch12.inertial_spectrum_3d(10.0, 1.0)` | 0.32727 = 18/55 ✎; 0.03232 m³/s² ✎ (= 1.5 × 10^{−5/3}) |
| C09 · N103 (row 5) | `ch12.plane_jet_reynolds_stress(1.0, 0.1, 1.0, 1.0, C5="from_invariant", xi_half=0.10)`; the same at (4.0, 0.4) and at (1.0, 0.0) | −0.1345 m²/s² ($-\overline{uv}<0$ where $dU/dy<0$: the sign an eddy viscosity gives); −0.0336 = one quarter (∝ x⁻¹, same ξ); 0 on the axis |
| C09 · N104 (row 5) | `ch12.plane_jet_cross_velocity(1.0, 0.1, 1.0, 1.0, C5="from_invariant", xi_half=0.10)`; the same at y = 1.0 | +0.0245 m/s (outward inside the jet); −0.1372 m/s at the edge = −`plane_jet_entrainment_velocity(1.0, 1.0, 1.0, …)`: the surroundings flow **in** |
| C09 · N111 (row with D14) | `out = ch12.plane_jet_similarity_sympy()`; `ch12.plane_jet_stress_profile(0.1, xi_half=0.10)` | `out["check"] is True`; `out["coefficients"]` are the three brackets $\{\delta U'_{CL}/U_{CL}\}$, $\{\delta U'_{CL}/U_{CL}+\delta'\}$, $\{\Psi/U_{CL}^2\}$ of the similarity equation, as written in D14's Result; `out["exponential_family"]` has a zero middle coefficient (slip #5); G(0.1) = −0.02025 |
| C09 · N117, N119 (row 14) | `ch12.plane_jet_mass_fraction(1.0, 0.0, 1.0, 1.0, 1.0, C6=1.0, xi_half_Y=0.15)`; the same at x = 4.0; then with `C6="from_invariant", C5="from_invariant", xi_half=0.10` | 1.0 and 0.5 ✎ (∝ x^{−1/2}, C₆ = 1 labelled illustrative); 2.240 when C₆ is fixed by the flux invariant $\dot M_s\cong\rho\int\bar YU\,dy$ (12.70) for our illustrative widths |
| C09 · N120, N121 (row 14) | `ch12.slot_momentum_flux(1.2, 10.0, 0.01)`; `ch12.slot_mass_flux(1.2, 10.0, 0.01)`; `ch12.free_shear_centerline("plane_jet", x, constants={…}, d=0.01, U0=10.0, rho_s=1.2, rho=1.2)` | 1.2 N/m ✎; 0.12 kg/(m s) ✎; centreline speed at x = 0.4 m is twice that at 1.6 m (x^{−1/2}), whatever the illustrative constants |
| C09 · N123 (row 15) | `ch12.FREE_SHEAR_CONSTANTS` printed | `{}` — "no public table adopted; every constant in this notebook is passed explicitly and labelled" |
| C09 · N124 (row 15, under D16) | `[ch12.local_reynolds_number_exponent(f) for f in flows]` | ½, 0, 0, −⅓, 1, ⅔, 1 for plane jet, round jet, plane wake, round wake, plane plume, round plume, shear layer (= D16's n + m) |
| C09 · N127 (row 18) | `ch12.plane_jet_stress_profile(xi, xi_half=0.10)` as the shear-stress curve of the qualitative panel | odd in ξ, zero on the axis, extreme near ξ ≈ ± ξ½ |
| C10 · N133 (row 6, under D17) | `st = WT.total_stress(y, U, uv, mu, rho)` on the model channel | `st["total"]` is a straight line through 0 at mid-height; `st["viscous"] + st["reynolds"]` equals it |
| C10 · N136 (row 6) | `WT.channel_total_stress(0.025, 0.1, 0.3)` | 0.15 Pa ✎ (a quarter of the way up a channel of **full** height h = 0.1 m: half the wall stress) |
| C10 · R04 (row 8) | `WT.pipe_pressure_gradient(0.3, 0.1)` | −12 Pa/m ✎ (= −4τ₀/d) — the inverse of the `wall_stress_from_pressure_gradient(-12.0, 0.1)` already in the row |
| C10 · N142 (row 13, worked example) | `WT.from_wall_units(3.0, 3.0, 0.5, 1.5e-5)` | (9.0 × 10⁻⁵ m, 1.5 m/s) ✎ — step 4 of the worked example backwards |
| C11 · N143 (start of D19) | `WT.defect_law_groups()` | two groups: $U\sqrt{\rho/\tau_0}$ (= U/u_*) and y/δ — the outer twin of `law_of_the_wall_groups()` |
| C11 · N144 (under D19) | `WT.velocity_defect(np.array([0.01, 0.05]), np.array([8.0, 9.0]), 10.0, 0.5, 0.1)` | ξ = (0.1, 0.5), defect = (4.0, 2.0) ✎ |
| C11 · N148 (under D19) | `WT.log_law_defect(0.1, kappa=0.41, A=1.0)` | 6.616 ✎ (= −ln(0.1)/0.41 + 1; A = 1.0 illustrative) |
| C11 · N153 (row with Spalding) | `WT.spalding_yplus(10.0, kappa=0.41, B=5.0)`; `WT.spalding_uplus(100.0, kappa=0.41, B=5.0)`; `WT.spalding_slope(10.0, kappa=0.41, B=5.0)` | 14.55; 16.08; 0.288 — one formula from the sublayer to the log layer (κ, B the illustrative pair) |
| C11 · N157 (row 8) | `B = WT.nagib_chauhan_B(0.384)`; `WT.nagib_chauhan_kappa(B)` | round trip returns 0.384 (print B from the cell, do not type it; review M2) — the two directions of $\kappa B=1.6[\exp(0.1663B)-1]$ (12.92) |
| C12 · N162 (row with the eddy-viscosity hypothesis) | `ch12.eddy_viscosity_stress(gradU, 0.01, 0.3)` with dU/dy = 2 s⁻¹ | $\overline{uv}$ = −0.02 m²/s² ✎ (= −ν_T dU/dy); each normal stress ⅔ē = 0.2 |
| C12 · N165 (row 5) | `ch12.rans_eddy_viscosity_residual(U, P, nu_T, e, x, nu=…, rho=…)` on a manufactured field (ν_T constant, P chosen to balance), then with `printed=True` | residual ≲ 1e-6; the printed form (pressure gradient with index j, slip #6) leaves an O(1) residual, the same number in every component |
| C12 · N168 (start of D21) | `sol = ch12.shear_flow_eddy_viscosity_solve(y, lambda yf, s: 0.0*yf, dPdx, rho, nu)` then with the mixing-length rule | ν_T ≡ 0 returns the laminar parabola to 1e-10 (Ch. 8); with $\nu_T=l_T^2\lvert dU/dy\rvert$ the profile flattens (`sol["converged"]` True) — the Picard primer P302 at work |
| C13 · N177 (row with the ε̄ equation) | `ch12.k_epsilon_rhs(1.0, 1.0, 0.5)`; `ch12.k_epsilon_rhs(1.0, 1.0, 0.0)` | (−0.5, −1.2) ✎: dē/dt = P − ε̄, dε̄/dt = (1.44 × 0.5 − 1.92) ε̄/ē; and (−1.0, −1.92) with no production — the pair D23 integrates |
| C14 · N182 (under D24) | `b = ch12.stratified_tke_budget(z, U, uw, wT, eps, alpha)` with the worked example's numbers | `b["Rf"]` = 0.0727 at the worked height, `b["regime"]` = "shear-driven"; buoyancy removes 7 % of the shear production |
| C14 · N184 (row 6) | `ch12.flux_from_gradient_richardson(0.0647, 0.891)` | 0.0726 ✎ (= Ri/Pr_T): the worked example's Rf recovered from Ri |
| C15 · the definition of $L_M$ (code row) | `ch12.monin_obukhov_length(0.3, 0.1, T=300.0, kappa=0.41)`; with wT = −0.02; with wT = 0.0 | −20.1 m (unstable); +100.7 m (stable); inf (neutral) — and the first equals `monin_obukhov_from_fluxes(0.108, 120.6, 1.2, 1005.0, 300.0, kappa=0.41)` |
| C15 · N187 (row 3) | `ch12.surface_layer_regime(z, L_M)` for (2 m, −20 m), (100 m, −20 m), (10 m, +100 m), (300 m, +100 m), (10 m, inf) | "forced convection", "free convection", "forced convection", "stable", "neutral" — below ∣L_M∣ the layer is shear-driven whatever the sign of the heat flux |
| C15 · N190 (row with the temperature-variance budget) | `ch12.temperature_variance_budget(z, T_mean, wT, eps_T)` | production $-\overline{wT'}\,d\bar T/dz>0$ whenever the flux runs down the gradient; `["dissipation"]` = −ε̄_T |
| C15 · N193 (row with the scalar spectrum) | `ch12.batchelor_scale(1e-6, 1.4e-7, 1e-6)` | 3.74 × 10⁻⁴ m ✎ (η = 1 mm for sea water with ε̄ = 10⁻⁶ m²/s³; × (κ_th/ν)^{1/2} = 0.374): temperature has finer structure than velocity |
| C16 · N200 (under D26) | `ch12.taylor_dispersion_rate(t, lambda s: np.exp(-s/10), 1.0)` at t = 1, 10, 100 s | 1.903, 12.64, 20.00 m²/s — rises from 2u²t to 2u²Λ_t |
| C16 · N213 (under D28) | `ch12.eddy_diffusivity_taylor(t, lambda s: np.exp(-s/10), 1.0)` at the same times | 0.952, 6.321, 10.00 m²/s = half the line above = `eddy_diffusivity_exponential(t, 1.0, 10.0)` (the worked example's 6.3 and 10) |

**A.15b — repeat mentions of an equation by number: completed (inlining pass of 2026-10-07).** The rule "the equation is
written beside its number" also holds for a second mention inside the same block, for a recap's `where=`, for a table cell
and for a reference to an earlier chapter. The audit had left about 14 such mentions in Part A for the builders to finish;
they are now finished in the rows themselves, each in one of two ways:

- **the compact equation stands beside the number** — copied from the place in this file where the block writes it out,
  compared with `analysis/ch12.md` §2 and with the rendered page (header, second provenance statement);
- **or, where the equation is long and stands in full a few rows above, the number is gone and the row points to it in
  words** ("the mean-flow budget of D09 above", "the Reynolds-stress equation of N59 above", "the exact budget (C06)").

Builders copy the rows as they stand; nothing is left for them to inline. Rows touched: C02 N27 · C04 N58, N59, N60 · the
`nb.core` title of C05 (the range of numbers is gone; rows 3–8 of that block write each of the four equations out) · C06 N81
(twice) · C07 N91 (twice) · C09 N128 · C10 R04 and N138 · C11 row 4 · C12 N165 · C13 N175 and N179 · C14 R05 and R06 (their
`where=` now names chapter and section only; the earlier-chapter equations are in the recap text) · C16 N211 and N215 · five
locator cells of A.15a. One correction came out of reading the earlier chapter's page: R04 had quoted Chapter 8's pipe wall
stress with a minus sign; the book prints $\tau_0=\frac a2\frac{dp}{dz}$ (8.8), where $\tau_0$ carries the sign of the pressure
gradient, and R04 now shows that form and says why Chapter 12's positive $\tau_0$ needs the minus sign.

---

## Part B — explainer storyboards

Common to all ten (and the backup): created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch12">`; tabs
Walkthrough · Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function
that mirrors a `ch12` callable and is proved by `selftest()` parity rows (`py:` expressions use only `ch12.…`, `np.…`, numbers,
strings, dict keys and integer indices — Part C convention 2 and C.5). Explain is "Explanation & interpretation" in numbered
sections built with `Viz.work.step / line / box / table / hint / interpret`, modelled on `forced_damped_vibrations.html` (the
regime-dependent reading) and `amplitude_phase_second_order_II_3.html` (a numbered derivation with live numbers, one section per
optional view): **0** what the views show and what each colour means · **1…n** every displayed quantity from the controls
("formula = substituted = result — why", results boxed) · a section or hint per view hidden on phones · the values at the
current time (live) · **Reading the current setting** (regime-dependent). Derivation steps are copied from Part F (same `did`
titles, same step count, same order; phones shorten *why* to its first sentence; plain-text *why* and *watch* never contain raw
TeX). Every tour, Explain, Derivation, notes, status, equation and quiz text that names a book equation **writes it out** next
to its number. Drafts below write equations in Unicode for readability; builders set each in TeX (backslashes doubled in JS
strings; never through a shell heredoc). Colours as convention 9 (mean purple, fluctuation teal, Reynolds stress / production
orange, viscous / dissipation rose, buoyancy and temperature blue, scalar amber, ghosts muted). Walkthrough texts ≤ 45 words,
step 1 ≤ 24 words, ≤ 2 extras per step, `play: false` on steps that quote numbers. A view hidden on portrait phones never
carries a step's key number. Random signals and particles use `Viz.rng(seed)` with the exact Ornstein–Uhlenbeck update; parity
never uses a sampled value (C.5). Parallel builders use private scratch subfolders (`<scratchpad>/<slug>/`).

### E1 · reynolds_averaging_window
- **Title:** "How long must you average to get 'the mean'?" · **Summary:** "Many runs or one long record: drag the number of
  members and the window and watch each estimate of the mean converge — or fail." · **CORE:** C01, C02 (also N07, N08, N11, N18,
  N20, N23, N32) · **Reference:** `forced_damped_vibrations.html` (signal + toggled curves on one time slider; the Explain panel).
- **meta:** `viz:order 1` · `viz:sections 12.3 12.4` · `viz:equations 12.1 12.2 12.6 12.17` · `viz:fluidpy ch12.time_average_exp_cos
  ch12.standard_error_of_mean ch12.effective_samples ch12.make_ensemble ch12.time_average ch12.product_average_split` ·
  `viz:derivations D01 D03`.
- **Physics (JS ↔ Python):** `windowFactors(dt, taud, omega)` → [sinh(x)/x with x = Δt/2τ_d, sin(x)/x with x = ωΔt/2] ↔
  `ch12.time_average_exp_cos(t, A, B, tau, omega, window)[1]`, `[2]`; `avgExpCos(t, A, B, taud, omega, dt)` ↔ `[0]`; `stdErr(sigma,
  n)` ↔ `ch12.standard_error_of_mean(sigma, N)`; `nEff(T, tc)` ↔ `ch12.effective_samples(T, tc)`; `ouMember(seed, …)` (exact
  update $u_{n+1}=u_ne^{-\Delta t/\tau_c}+\sigma\sqrt{1-e^{-2\Delta t/\tau_c}}\,\xi_n$) ↔ `TS.make_ensemble` (same law, not the same numbers — no parity);
  `slidingMean(u, w)` (cumulative sums) ↔ `TS.time_average`; `productSplit(us, vs)` ↔ `TS.product_average_split`.
- **Views** (rows [1.2, 1]): 1. `record` "The record" (row 0, full width): t 0…20 s; members grey (≤ 8 drawn), the ensemble
  mean of N members bold purple, the sliding time average of member 1 bold teal with its window drawn as a teal band
  centred on the cursor, the true mean as a muted dashed ghost; bold only up to the transport time ("so far"), faint beyond.
  Pointer: drag the cursor. 2. `err` "Error of each estimate" (row 1, flex 1): log–log; purple dots = rms error of the
  ensemble mean vs N with the σ/√N line; a teal curve = error of the time average vs Δt (the leak ∣sin x/x∣·B and the drift
  ∣sinh x/x − 1∣·A drawn as two thin lines whose sum is the curve), current N and Δt marked ◆. 3. `prod` "Mean of a product"
  (row 1, flex 0.7, `hidePortrait`): three bars $\bar u\bar v$ (purple), $\overline{uv}$ (orange), $\overline{\tilde u\tilde v}$ (outline) with the sum check.
- **Controls:** `N` "Members $N$" (2, 4, 8, 16, 32, 64; default 8) · `dt` "Window $\Delta t$" (0.05…20 s log, default 0.4, s) ·
  `tauc` "Memory $\tau_c$" (0.05…2 s, default 0.2) · `sig` chips stationary / decaying mean / mean + wave (default mean + wave: A = 1,
  B = 0.5, τ_d = 10 s, period 1 s) · `ruv` "Correlation of $u$ and $v$" (0…1, default 0.6; optional) · transport `t`.
- **Transport:** `t` 0…20 s, rate 2, `end: 'hold'`; end card: "N = 8: scatter 0.106 ≈ σ/√8 · window 0.4 s: 76 % of the wave
  leaked".
- **Presets:** "too short a window" {dt: 0.4} · "one whole period" {dt: 1.0} · "good window" {sig: decaying, dt: 3, tauc: 0.2}
  · "too long" {dt: 20} · "N = 2" · "N = 64".
- **Status:** wave signal: "⚠️ window 0.4 periods: 76 % of the wave leaks into the mean" (∣sin x/x∣ > 0.1) / "✅ good window:
  1 ≪ ωΔt ≪ ωτ_d (leak 0 %, drift +0.4 %)" / "⚠️ window too long: the decaying mean is distorted by +17.5 %"; stationary
  signal: "window holds Δt/Λ_t = 15 memory times → error ≈ σ/√(15/2)".
- **Readouts:** "Leak sin x/x" · "Drift sinh x/x" · "σ/√N" · "Δt/Λ_t" · "⟨ũṽ⟩ − ŪV̄".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **terms**
  ($\bar u\bar v$ + $\overline{uv}$ = $\overline{\tilde u\tilde v}$; click a bar to isolate it).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The record:** grey = single runs; <b class=c-accent>purple</b> = average of your N runs at each
     instant (the ensemble average ⟨u^m⟩ = lim (1/N) Σ u(x, t:n)^m (12.1) with m = 1 and finite N); <b class=c-teal>teal</b> =
     average of run 1 over a window Δt (ū = (1/Δt)∫u dt (12.2)); dashed = the true mean. **Error:** how far each estimate is
     from the truth. **Product** (hidden on phones): the mean of a product split into its two parts."
  1. *The ensemble estimate* — "scatter of the mean of N independent runs = σ/√N = 0.30/√8 = **0.106** — halving it costs four
     times the runs."
  2. *The window and the wave* — "x = ωΔt/2 = 2π × 0.4/2 = 1.257; sin x/x = 0.951/1.257 = **0.757**: this fraction of the wave
     survives the averaging. It is zero whenever the window holds a whole number of periods."
  3. *The window and the drifting mean* — "x = Δt/2τ_d = 0.4/20 = 0.02; sinh x/x = **1.00007**: the mean passes untouched while
     Δt ≪ τ_d."
  4. *How many independent samples one record holds* — "Λ_t = τ_c = 0.2 s; Δt/Λ_t = **2** memory times in the window; for a
     stationary signal the time average scatters by about σ/√(Δt/2Λ_t)."
  5. *The product* (or the hint "turn the phone sideways for the bars") — "Ū V̄ = … ; ⟨uv⟩ = r σ_u σ_v = 0.6 × 0.3 × 0.3 =
     **0.054**; ⟨ũṽ⟩ = Ū V̄ + ⟨uv⟩ — the average passed through everything except this product."
  6. *At the current time* — "t = `Viz.live('t')` s: true mean `Viz.live('true')`, ensemble estimate `Viz.live('ens')`, window
     estimate `Viz.live('win')`."
  7. *Reading the current setting* — short window: "The window is shorter than the wave's period, so the 'mean' still wobbles
     with it: you have not averaged, only smoothed." · whole period: "Exactly one period fits: the wave cancels itself. Any
     whole number of periods does this — which is why tidal analyses use windows of whole tidal cycles." · good: "Long against
     the fluctuation, short against the drift: the window sees a constant mean plus noise. This is the 1 ≪ ωΔt ≪ ωτ_d of
     Example 12.1, and the reason a 30-year climate normal is neither 3 years nor 300." · too long: "The window now spans a
     good part of the decay and overestimates the mean: sinh x/x > 1." · stationary: "Nothing drifts, so longer is always
     better: the error falls as (Δt/Λ_t)^{−1/2}."
- **Derivation tab:** **D01** (8 steps) `view: 'prod'` on wide screens, `'record'` on phones; goal `set` {N: 2, sig: stationary};
  step 2 `live` "(1/2)(u₁ + u₂) with your two members = …"; step 7 `set` {ruv: 0.6, N: 64}, `live` "⟨ũṽ⟩ = … = Ū V̄ + ⟨uv⟩ = … +
  …"; step 8 `watch` "set the correlation to 0: the orange bar vanishes and the product rule 'works' — by accident".
  **D03** (5 steps) `view: 'record'`, `set` {sig: stationary}; step 2 `watch` "slide the time cursor: the statistics of the grey
  band do not change — that is what stationary means". Interpret: `s => "With N = " + s.N + " the scatter is " + …`.
- **Code:**
  ```python
  ens = TS.make_ensemble({{N}}, t, mean_fn, sigma=0.3, tau_c={{tauc}}, seed=0)
  U = TS.ensemble_average(ens)              # (12.10): scatter sigma/sqrt(N) = {{se}}
  ubar = TS.time_average(t, ens[0], {{dt}})  # (12.2): window of {{dt}} s
  avg, f_mean, f_wave = ch12.time_average_exp_cos(t, 1.0, 0.5, 10.0, 2*np.pi, {{dt}})
  # f_mean = sinh(x)/x = {{fm}}   f_wave = sin(x)/x = {{fw}}
  mp, pm, cov = TS.product_average_split(ens_u, ens_v)   # {{mp}} = {{pm}} + {{cov}}
  ```
- **Walkthrough (6 steps):** 1. "What is 'the mean'?" — "The signal never repeats. Press ▶ and watch two ways of averaging
  it." `play: true` · 2. "Across runs" — "Purple averages N runs at each instant. Its scatter is σ/√N: try 2, 8, 64."
  `controls: ['N']`, `readouts: ['se']` · 3. "Along one run" — "Teal averages one run over a window. At 0.4 periods, 76 % of
  the wave is still there." `set` {dt: 0.4}, `eq: 'win'` · 4. "One whole period" — "Set Δt to one period: sin(π)/π = 0. The
  wave is gone." `set` {dt: 1}, `code: {id: 'avg', lines: [4, 5]}` · 5. "The rule that fails" — "An average passes through
  sums and derivatives, not products: ⟨ũṽ⟩ = Ū V̄ + ⟨uv⟩." `terms: true`, `derive: {id: 'D01', step: 7}` · 6. "Your turn" —
  "Predict the window that distorts the decaying mean by 5 %. Then find it." `controls: ['dt', 'sig']`.
- **Equations:** `ens` "Ensemble average" ref 'Eq. (12.1)' ⟨u^m⟩ = lim_{N→∞} (1/N) Σ_n (u(x, t:n))^m, live "N = 8" · `win` "Time
  average" ref 'Eq. (12.2)' ū^m = (1/Δt)∫_{t−Δt/2}^{t+Δt/2} u^m dt, live with the two factors · `ddt` "Averaging and ∂/∂t commute"
  ref 'Eq. (12.6)' $\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t$, note "exact for ensembles; approximate for a window" · `prod` "The product"
  ref 'D01 (no book number)' $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$ live · `lag` "Stationary correlation" ref 'Eq. (12.17)' R₁₁(τ) = ⟨u₁(t)u₁(t + τ)⟩.
- **Check yourself:** (1) "How many members halve the scatter of N = 16?" — "64: scatter ∝ N^{−1/2}." `set {N: 16}` · (2) "Which
  windows remove the wave completely?" — "Whole numbers of periods: sin(ωΔt/2) = 0." `set {dt: 2}` · (3) "With correlation 0,
  does ⟨ũṽ⟩ = Ū V̄?" — "Yes — but only then; the rule fails in general." `set {ruv: 0}` · (4) "Why is a 20 s window bad for the
  decaying mean?" — "Δt/2τ_d = 1: sinh(1)/1 = 1.175, a 17.5 % overestimate."
- **Selftest parity rows:** `{name: 'wave factor', js: windowFactors(0.4, 10, 2*Math.PI)[1], py: 'ch12.time_average_exp_cos(5.0, 1.0,
  0.5, 10.0, 2*np.pi, 0.4)[2]', rtol: 1e-12}` (0.75683) · `{name: 'mean factor', js: windowFactors(20, 10, 2*Math.PI)[0], py:
  'ch12.time_average_exp_cos(5.0, 1.0, 0.5, 10.0, 2*np.pi, 20.0)[1]', rtol: 1e-12}` (1.17520) · `{name: 'average', js: avgExpCos(5, 1,
  0.5, 10, 2*Math.PI, 0.4), py: 'ch12.time_average_exp_cos(5.0, 1.0, 0.5, 10.0, 2*np.pi, 0.4)[0]', rtol: 1e-12}` (0.98498) · `{name:
  'std error', js: stdErr(0.3, 8), py: 'ch12.standard_error_of_mean(0.3, 8)', rtol: 1e-12}` · invariant `{name: 'OU variance', js:
  ouVariance(1.0, 0.2, 40000), expect: 1.0, rtol: 0.05}`.
- **Fit plan:** 360×640: status (1 line) · `record` (60 %) over `err` (40 %); `prod` hidden (its number is the readout "⟨ũṽ⟩ −
  ŪV̄"); transport without step buttons; walkthrough card paged. 844×390: `record` | `err`. Desktop / 1000×700: rows [1.2, 1].

### E2 · correlation_and_spectrum
- **Title:** "Why are a correlation and a spectrum the same thing?" · **Summary:** "Slide a signal past itself to build its
  autocorrelation; stretch its memory and watch the spectrum squeeze while its area stays the variance." · **CORE:** C02, C03
  (also N30, N31, N33, N37, N38, N40) · **Reference:** `amplitude_phase_second_order_II_3.html` (linked windows, crosshair
  readouts, numbered live derivation).
- **meta:** `viz:order 2` · `viz:sections 12.4` · `viz:equations 12.17 12.18 12.19 12.20 12.21 12.22` · `viz:fluidpy
  ch12.correlation_spectrum_pair ch12.integral_scale ch12.taylor_microscale ch12.spectrum_from_correlation ch12.spectrum_variance
  ch12.frequency_to_wavenumber_spectrum` · `viz:derivations D03 D04 D05`.
- **Physics (JS ↔ Python):** `pair(kind, sigma, tauc, w0)` → {r(τ), S(ω), Lambda, lambda, S0, tc} with the closed forms of
  Part C.1 ↔ `ch12.correlation_spectrum_pair(kind, sigma, tau_c, omega0)`; `cosTransform(R, w)` (Simpson to 40 τ_c) ↔
  `TS.spectrum_from_correlation`; `area(S)` ↔ `TS.spectrum_variance`; `toWavenumber(w, S, U0)` ↔
  `TS.frequency_to_wavenumber_spectrum`; a record of 2¹² samples is synthesised for the faint periodogram (random phases with
  amplitude √(2SΔω), `Viz.rng`).
- **Views** (rows [0.8, 1.2]): 1. `sig` "The signal and its shifted copy" (row 0): u(t) teal, u(t + τ) orange, their product
  shaded (green where positive, rose where negative); title carries the running mean of the product = R(τ). Pointer: drag
  the copy sideways (sets τ). 2. `corr` "r(τ)" (row 1, flex 1): the curve teal, a dot at the current lag, the equal-area
  rectangle (height 1, width Λ_t) purple, the first zero t_c ▲, the osculating parabola 1 − τ²/λ_t² rose dashed (Gaussian
  mode only; "no parabola: cusp" label otherwise). 3. `spec` "S_e(ω)" (row 1, flex 1; log–log): exact transform bold blue,
  periodogram faint, the level S(0) dashed with its marker, corner frequency 1/Λ_t ▼, area shaded; pointer: click → inspector.
- **Controls:** `kind` chips exponential / Gaussian / damped cosine (mode) · `tauc` "Memory $\tau_c$" (0.05…5 s log, default 0.5)
  · `w0` "Oscillation $\omega_0$" (0…20 rad/s, default 6; damped cosine only) · `lag` "Lag $\tau$" (0…5 τ_c, default 0.3 s) · `axis`
  toggle frequency / wavenumber + `U0` "Probe speed $U_0$" (1…50 m/s, default 10; optional).
- **Transport:** `lag` sweeps 0 → 5 τ_c (rate 0.5 τ_c/s, `end: 'hold'`): the correlation curve is drawn point by point (bold "so
  far", faint ahead); end card "area under r = Λ_t = 0.50 s; S(0) = σ²Λ_t/π = 0.159".
- **Presets:** "short memory" {kind: exponential, tauc: 0.1} · "long memory" {tauc: 2} · "a hidden wave" {kind: damped cosine,
  tauc: 2, w0: 6} · "Gaussian pair" {kind: gaussian, tauc: 0.5}.
- **Status:** "memory Λ_t = 0.50 s ⇒ spectrum flat to ω ≈ 1/Λ_t = 2.0 rad/s" · damped cosine with ω₀τ_c > 1: "🌊 hidden wave: peak
  at ω₀ = 6 rad/s; Λ_t only 0.05 s — the lobes cancel".
- **Readouts:** "Λ_t" [s] · "t_c" · "λ_t" · "S(0)" · "∫S dω".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (three shapes) · **presets** · **status** ·
  **transport** · **inspector** (click the spectrum at ω = 2: "S_e(2) = (1/π)∫₀^∞ R cos(2τ) dτ = σ²τ_c/[π(1 + ω²τ_c²)] = 0.5/(π × 2) =
  **0.0796** = ½ S(0)").
- **Explain:**
  0. *What the views show* — "Top: the signal (<b class=c-teal>teal</b>) and a copy delayed by τ (<b class=c-orange>orange</b>);
     the shaded product, averaged, is one point of the correlation curve. Bottom left: that curve for every lag. Bottom right:
     its Fourier transform — how the variance is shared among frequencies."
  1. *The correlation at your lag* — "R₁₁(τ) = ⟨u₁(t)u₁(t + τ)⟩ (12.17); r(0.3) = e^{−0.3/0.5} = **0.549**."
  2. *The memory time* — "Λ_t = ∫₀^∞ r₁₁ dτ (12.18) = τ_c = **0.50 s** (exponential); = (√π/2)τ_c (Gaussian); = τ_c/(1 + ω₀²τ_c²)
     (damped cosine)." boxed; "first zero t_c = π/2ω₀ (damped cosine) or none."
  3. *The Taylor microscale* — "λ_t² = −2/r″(0) (12.19): Gaussian r ≈ 1 − τ²/τ_c² ⇒ λ_t = τ_c = **0.50 s**. For the other two shapes
     r has a corner at τ = 0 and no parabola fits."
  4. *The spectrum* — "S_e(ω) = (1/2π)∫R₁₁e^{−iωτ}dτ (12.20) = σ²τ_c/[π(1 + ω²τ_c²)]; S(0) = σ²Λ_t/π = 1 × 0.5/π = **0.159**; it has
     fallen to half at ω = 1/τ_c = 2 rad/s."
  5. *Its area* — "∫S_e dω = ⟨u₁²⟩ (12.22) = **1.00** — whatever τ_c is. Halve the memory: S(0) halves and the corner doubles."
  6. *Wavenumber axis* (optional view) — "k₁ = ω/U₀, S₁₁(k₁) = U₀S_e(ω): with U₀ = 10 m/s the corner sits at k₁ = 0.2 rad/m (a
     31 m eddy)."
  7. *Reading the current setting* — short memory: "The signal forgets quickly, so it wiggles fast: variance is spread over a
     wide band of frequencies (nearly 'white')." · long memory: "Slow, persistent swings: the variance sits near zero
     frequency ('red' spectrum, like most climate records)." · hidden wave: "The correlation oscillates and the spectrum has
     a peak at ω₀; the net area under r — the integral scale — is small because positive and negative lobes cancel. A small
     Λ_t does not always mean a short memory." · Gaussian: "A smooth signal: its correlation is flat at the top (a parabola),
     and its spectrum falls off faster than any power."
- **Derivation tab:** **D03** (5 steps) `view: 'sig'`; step 3 `watch` "drag the copy to −τ: the shaded pattern is the mirror
  image, the average is the same". **D04** (7 steps) `view: 'corr'`, goal `set` {kind: gaussian}; step 4 `live` "r ≈ 1 − τ²/λ_t² with
  λ_t = 0.50 s"; step 7 `watch` "switch to 'exponential': the top is a corner, the parabola disappears". **D05** (9 steps)
  `view: 'spec'`; step 5 `live` "S_e(ω) = (1/π)∫₀^∞ R cos ωτ dτ = …"; step 9 `live` "S_e(0) = ⟨u²⟩Λ_t/π = …". Interpret from section 7.
- **Code:**
  ```python
  pair = TS.correlation_spectrum_pair("{{kind}}", 1.0, {{tauc}}, {{w0}})
  lag, R = TS.autocorrelation(u, dt)                # (12.17)
  Lam = TS.integral_scale(lag, R/R[0])              # (12.18): {{Lam}} s
  S = TS.spectrum_from_correlation(lag, R, omega)   # (12.20)
  var = TS.spectrum_variance(omega, S)              # (12.22): {{var}}
  S0 = pair["S0"]                                   # variance*Lam/pi = {{S0}}
  ```
- **Walkthrough (7 steps):** 1. "How long does it remember?" — "Slide the orange copy along the teal signal. Aligned, the
  product is all positive." `controls: ['lag']` · 2. "Build the curve" — "Each lag gives one average. Press ▶ to draw
  r(τ)." `play: true` · 3. "The memory time" — "The purple rectangle has the same area as the curve: its width is Λ_t =
  ∫₀^∞ r₁₁ dτ." `highlight: ['readout:Lam']`, `eq: 'int'` · 4. "Same information" — "Fourier-transform r: the spectrum. Its
  area is the variance, its height at zero σ²Λ_t/π." `derive: {id: 'D05', step: 9}` · 5. "Stretch the memory" — "Drag τ_c up:
  the spectrum squeezes left and grows taller; the shaded area never changes." `controls: ['tauc']` · 6. "A hidden wave" —
  "An oscillating correlation is a peak away from zero frequency." `set` {kind: 'damped_cosine', tauc: 2, w0: 6} · 7. "Your
  turn" — "Predict S(0) for τ_c = 1 s, then check. Click the spectrum for the arithmetic." `inspect: true`.
- **Equations:** `lag` ref 'Eq. (12.17)' R₁₁(τ) = ⟨u₁(t)u₁(t + τ)⟩ = R₁₁(−τ) · `int` ref 'Eq. (12.18)' Λ_t ≡ ∫₀^∞ r₁₁(τ)dτ, live · `mic`
  ref 'Eq. (12.19)' λ_t² ≡ −2/[d²r₁₁/dτ²]_{τ=0} · `spec` ref 'Eq. (12.20)' S_e(ω) ≡ (1/2π)∫R₁₁(τ)e^{−iωτ}dτ, live · `var` ref 'Eq. (12.22)'
  ⟨u₁²⟩ ≡ ∫S_e(ω)dω, live "= 1.00". (The inverse R₁₁(τ) = ∫S_e(ω)e^{+iωτ}dω (12.21) is quoted in D05.)
- **Check yourself:** (1) "Double τ_c. What happens to S(0) and to the area?" — "S(0) doubles (∝ Λ_t); the area stays σ²." · (2)
  "For which shape does the Taylor microscale exist?" — "Only the Gaussian: the others have a corner at τ = 0." `set {kind:
  gaussian}` · (3) "Damped cosine, τ_c = 2 s, ω₀ = 6: why is Λ_t only 0.014 s?" — "Λ_t = τ_c/(1 + ω₀²τ_c²) = 2/145: the lobes
  cancel." · (4) "At what frequency has an exponential-memory spectrum fallen to half?" — "ω = 1/τ_c."
- **Selftest parity rows:** `{name: 'S exp at 2', js: pair('exponential', 1, 0.5, 0).S(2), py: 'ch12.correlation_spectrum_pair
  ("exponential", 1.0, 0.5)["S"](2.0)', rtol: 1e-12}` (0.079577) · `{name: 'Lambda gaussian', js: pair('gaussian', 1, 0.5, 0).Lambda, py:
  'ch12.correlation_spectrum_pair("gaussian", 1.0, 0.5)["Lambda_t"]', rtol: 1e-12}` (0.443113) · `{name: 'S damped at w0', js:
  pair('damped_cosine', 1, 0.5, 6).S(6), py: 'ch12.correlation_spectrum_pair("damped_cosine", 1.0, 0.5, 6.0)["S"](6.0)', rtol: 1e-12}`
  (0.081728) · `{name: 'Lambda from S0', js: pair('exponential', 1, 0.5, 0).Lambda, py: 'ch12.integral_scale_from_spectrum(0.15915494309, 1.0)', rtol: 1e-9}` (0.5) · invariant `{name: 'area = variance', js: area('gaussian', 1, 0.5), expect: 1, rtol: 1e-6}`.
- **Fit plan:** 360×640: `sig` (30 %) over `corr` | `spec` side by side is too narrow → `corr` (35 %) and `spec` (35 %) stacked;
  status one line; the periodogram hidden on phones. 844×390: `corr` | `spec`, `sig` as a thin strip. Desktop: rows [0.8, 1.2].

### E3 · reynolds_stress_parcels
- **Title:** "How can fluctuations push the mean flow?" · **Summary:** "Release parcels in a shear: each lands in the (u, v)
  scatter, the cloud tilts, and the tilt is the Reynolds stress." · **CORE:** C04 (also N42, N48–N51, N60, N63) · **Reference:**
  `forward_noising_lab.html` (click a sample to see its arithmetic) + `forced_damped_vibrations.html` (Explain).
- **meta:** `viz:order 3` · `viz:sections 12.5 12.6` · `viz:equations 12.14 12.24 12.30` · `viz:fluidpy ch12.displaced_parcel_uv
  ch12.parcel_uv_expected ch12.reynolds_stress ch12.correlation_coefficient ch12.mean_stress_tensor` · `viz:derivations D06`.
- **Physics (JS ↔ Python):** each parcel draws v ~ N(0, v_rms²) and ℓ = l_rms[c·v/v_rms + √(1 − c²)·ξ] (c = "keeps its
  momentum"), then u = −ℓ dU/dy; `uvExpected(dUdy, l, v, c)` = −c·v·l·dU/dy ↔ `ch12.parcel_uv_expected(dUdy, l_rms, v_rms,
  correlation)`; `ruvExpected(c, dUdy)` = −c·sign(dU/dy); `stresses(rho, mu, dUdy, uv)` → [μ dU/dy, −ρ⟨uv⟩, total] ↔
  `ch12.mean_stress_tensor`; `nuT` = −⟨uv⟩/(dU/dy) ↔ `ch12.eddy_viscosity_from_data`; sample statistics ↔ `TS.reynolds_stress`,
  `TS.correlation_coefficient` (sampled: shown, never parity).
- **Views** (rows [1, 1]): 1. `flow` "Parcels in a shear" (row 0, flex 1.3): the mean profile U(y) as purple arrows; parcels
  hop between levels, coloured orange when their uv < 0 and blue when > 0; the stress element of Fig. 12.7 in a corner.
  Pointer: click a parcel → inspector. 2. `scatter` "(u, v)" (row 0, flex 1): dots accumulate; covariance ellipse and
  principal axes; quadrants 2 and 4 tinted orange; quadrant counts in the corners. 3. `bars` "Stress on the mean flow" (row 1,
  `hidePortrait`): viscous μ dU/dy (rose), Reynolds −ρ⟨uv⟩ (orange), total (outline), log scale with the ratio printed.
- **Controls:** `shear` "Shear $dU/dy$" (−4…4 s⁻¹, default 2) · `lrms` "Displacement $\ell_{rms}$" (0.01…0.5 m, default 0.1) ·
  `vrms` "$v_{rms}$" (0.05…2 m/s, default 0.5) · `keep` "Keeps its momentum" (0…1, default 0.8) · `fluid` chips air / water
  (optional) · transport `n` (parcels released).
- **Transport:** `n` 0…2000 parcels, rate 100/s, `end: 'hold'`; end card: "2000 parcels: ⟨uv⟩ = −0.079 ± 0.014 m²/s² (5 s.e.; expected −0.080); 80 % landed in quadrants 2 and 4" (the fraction is ½ + arcsin(0.8)/π for a Gaussian pair).
- **Presets:** "strong shear" {shear: 4} · "no shear — isotropic cloud" {shear: 0} · "reversed shear" {shear: −2} · "parcels that
  forget" {keep: 0.3}.
- **Status:** "dU/dy > 0 ⇒ ⟨uv⟩ = −0.080 m²/s² < 0: momentum flux toward the slow side" · "no shear: round cloud, no shear
  stress" · "keep = 0.3: r_uv = −0.30 — weakly correlated".
- **Readouts:** "⟨uv⟩ expected" · "⟨uv⟩ sampled ± 5 s.e." · "r_uv" · "−ρ⟨uv⟩" [Pa] · "ν_T implied".
- **Depth features:** Explain + Code + Derivation · **linked views** · **transport** (end card) · **presets** · **status** ·
  **terms** (viscous + Reynolds = total) · **inspector** (click a parcel: "ℓ = +0.12 m, v = +0.61 m/s ⇒ u = −ℓ dU/dy = −0.12 × 2 =
  −0.24 m/s; uv = −0.146 m²/s² (quadrant 2)").
- **Explain:**
  0. *What the views show* — "Left: parcels swapping levels in a mean shear (<b class=c-accent>purple</b> arrows). Each keeps
     its old mean speed for a while. Right: each parcel's (u, v). <b class=c-orange>Orange</b> = events with uv < 0. Bottom
     (hidden on phones): the stress those events exert, against the viscous stress."
  1. *One parcel* — "displaced by ℓ it arrives with u ≈ −ℓ dU/dy: too slow if it came from below (v > 0), too fast if from above
     (v < 0)."
  2. *The average* — "⟨uv⟩ ≈ −⟨vℓ⟩ dU/dy = −c·v_rms·ℓ_rms·dU/dy = −0.8 × 0.5 × 0.1 × 2 = **−0.080 m²/s²**." boxed; "sampled: `Viz.live
     ('uv')` ± `Viz.live('se')` (5 standard errors)."
  3. *The correlation coefficient* — "r₁₂ = ⟨u₁u₂⟩/(√⟨u₁²⟩√⟨u₂²⟩) (12.14) = −0.080/(0.2 × 0.5) = **−0.80**; in real shear flows it is about −0.4."
  4. *The stress* — "τ̄_ij = −Pδ_ij + 2μS̄_ij − ρ₀⟨u_iu_j⟩ (12.30): Reynolds part −ρ₀⟨uv⟩ = 1.2 × 0.080 = **0.096 Pa**; viscous part μ dU/dy = 1.8 × 10⁻⁵ × 2 = 3.6 × 10⁻⁵ Pa — **2700 times smaller**. (Water: 80 Pa and 2 × 10⁻³ Pa.)"
  5. *The viscosity it mimics* — "ν_T = −⟨uv⟩/(dU/dy) = c·v_rms ℓ_rms = **0.04 m²/s** against ν = 1.5 × 10⁻⁵: a property of the
     motion, not of the air."
  6. *Reading the current setting* — positive shear: "Going up goes with being slow. No single parcel pushes the flow; the
     push is the *correlation*. Momentum is carried from the fast side to the slow side, as friction would — only thousands
     of times more effectively." · zero shear: "With no mean gradient to remember, u and v are unrelated: the cloud is round
     and the shear stress is zero — the isotropic case." · negative shear: "The cloud tilts the other way; ⟨uv⟩ > 0. The
     stress always opposes the shear." · low keep: "Parcels that mix quickly with their surroundings carry less momentum:
     the correlation and the stress fall together."
- **Derivation tab:** **D06** (13 steps) `view: 'scatter'`; step 6 (expand the product) `live` "(U + u)(V + v) with the clicked
  parcel's numbers"; step 8 (the cross terms vanish) `watch` "the dots are centred: ⟨u⟩ = ⟨v⟩ = 0, but the cloud is tilted";
  step 9 `set` {shear: 2, keep: 1}, `live` "⟨uv⟩ = −0.100"; step 13 `live` "−ρ₀⟨uv⟩ = 0.12 Pa". Interpret: section 6.
- **Code:**
  ```python
  out = ch12.displaced_parcel_uv({{shear}}, {{lrms}}, {{vrms}}, n=100000, seed=0,
                                 correlation={{keep}})
  uv = out["uv"]                               # sampled <uv> = {{uvs}} m2/s2
  uv_exp = ch12.parcel_uv_expected({{shear}}, {{lrms}}, {{vrms}}, {{keep}})  # {{uve}}
  tau_R = -1.2 * uv_exp                        # Reynolds stress in air: {{tauR}} Pa
  tau_v = 1.8e-5 * {{shear}}                   # viscous stress: {{tauv}} Pa
  nu_T = -uv_exp / {{shear}}                   # implied eddy viscosity: {{nuT}} m2/s
  ```
- **Walkthrough (6 steps):** 1. "Zero on average, yet a push?" — "Fluctuations average to zero. Press ▶ and watch where the
  parcels land." `play: true` · 2. "Up means slow" — "A parcel from below brings its lower speed: u < 0 with v > 0. Click
  one." `inspect: true` · 3. "The cloud tilts" — "Both kinds of event give uv < 0: ⟨uv⟩ = −0.080 m²/s²." `highlight:
  ['view:scatter']`, `eq: 'parcel'` · 4. "It is a stress" — "−ρ₀⟨uv⟩ = 0.096 Pa: 2700 times the viscous stress here." `terms:
  true`, `derive: {id: 'D06', step: 13}` · 5. "Take the shear away" — "No shear, round cloud, no stress. Reverse it and the
  tilt flips." `set` {shear: 0}, `controls: ['shear']` · 6. "Your turn" — "Predict ⟨uv⟩ when parcels keep only half their
  momentum; then set it." `controls: ['keep']`.
- **Equations:** `dec` ref 'Eq. (12.24)' ũ_i = U_i + u_i · `rans` ref 'Eq. (12.30)' ∂U_i/∂t + U_j∂U_i/∂x_j = −g[1 − α(T̄ − T₀)]δ_i3 +
  (1/ρ₀)∂τ̄_ij/∂x_j, τ̄_ij = −Pδ_ij + 2μS̄_ij − ρ₀⟨u_iu_j⟩ (aligned over three lines), live "−ρ₀⟨uv⟩ = …" · `flux` ref 'p. 557 (no
  number)' ρ₀⟨(U + u)v⟩ = ρ₀⟨uv⟩ · `parcel` ref 'Fig. 12.6 argument' u ≈ −ℓ dU/dy, ⟨uv⟩ ≈ −⟨vℓ⟩ dU/dy, live · `r12` ref 'Eq. (12.14)'
  r₁₂ = ⟨u₁u₂⟩/(√⟨u₁²⟩√⟨u₂²⟩), live.
- **Check yourself:** (1) "Double the shear. What happens to ⟨uv⟩ and to r_uv?" — "⟨uv⟩ doubles; r_uv does not change (u
  doubles too)." `set {shear: 4}` · (2) "Can ⟨uv⟩ be non-zero when ⟨u⟩ = ⟨v⟩ = 0?" — "Yes: it is a covariance; the dots are
  centred but tilted." · (3) "In which quadrants do most parcels land when dU/dy < 0?" — "1 and 3: ⟨uv⟩ > 0." `set {shear: -2}` ·
  (4) "Why is ν_T here 0.04 m²/s whatever the shear?" — "ν_T = c·v_rms ℓ_rms: a velocity times a length of the motion (12.98's
  ν_T ~ l_T u_T)."
- **Selftest parity rows:** `{name: 'uv expected', js: uvExpected(2, 0.1, 0.5, 1), py: 'ch12.parcel_uv_expected(2.0, 0.1, 0.5, 1.0)',
  rtol: 1e-12}` (−0.1) · `{name: 'uv expected, forgetful', js: uvExpected(-3, 0.2, 0.4, 0.3), py: 'ch12.parcel_uv_expected(-3.0, 0.2,
  0.4, 0.3)', rtol: 1e-12}` (+0.072) · `{name: 'nu_T', js: nuT(-0.1, 2), py: 'ch12.eddy_viscosity_from_data(-0.1, 2.0)', rtol: 1e-12}` ·
  invariant `{name: 'sampled within 5 s.e.', js: sampledUvZ(2, 0.1, 0.5, 1, 20000), expect: 0, atol: 5}`.
- **Fit plan:** 360×640: `flow` (50 %) over `scatter` (50 %); `bars` hidden — the stress and the ratio are readouts; end card
  wraps (`Viz.card`). 844×390: `flow` | `scatter`. Desktop: bars under the scatter.

### E4 · energy_cascade_spectrum
- **Title:** "What does a higher Reynolds number change?" · **Summary:** "The big eddies fix the energy supply; viscosity only
  decides how far down the ladder it must go. Drag ΔU, L and ν and watch the cascade and the −5/3 range." · **CORE:** C07,
  C08 (also N85–N92, N94–N97, slip #1) · **Reference:** `overfitting_curves.html` (a marker and a verdict on a two-slider
  figure) + `angular_frequency_explorer_1.html` (table with the current row, presets).
- **meta:** `viz:order 4` · `viz:sections 12.7` · `viz:equations 12.49 12.50 12.51 12.52 12.53 12.54` · `viz:fluidpy
  ch12.dissipation_outer_scaling ch12.kolmogorov_scales ch12.scale_separation ch12.cascade_tiers ch12.model_spectrum
  ch12.inertial_spectrum_1d ch12.dns_grid_points ch12.inertial_range_decades` · `viz:derivations D11 D12`.
- **Physics (JS ↔ Python):** `epsOuter(dU, L)` ↔ `ch12.dissipation_outer_scaling`; `kolmogorov(nu, eps)` → [η, u_K, τ_η] ↔
  `ch12.kolmogorov_scales`; `separation(Re)` ↔ `ch12.scale_separation(Re)["eta_over_L"]`, `["lambdaT_over_L"]`; `tiers(L, dU, nu)`
  (sizes L/2ⁿ, u′ = (ε̄l′)^{1/3}, turnover l′/u′, Re(l′)) ↔ `ch12.cascade_tiers`; `Epao(K, eps, nu, L)` ↔ `ch12.model_spectrum(K, eps,
  nu, L, kind="pope")`; `S11inertial(k, eps, printed)` ↔ `ch12.inertial_spectrum_1d(k1, eps, two_sided=True, printed=…)`;
  `gridPoints(Re)` ↔ `ch12.dns_grid_points`; `decades(Re)` ↔ `ch12.inertial_range_decades`.
- **Views** (rows [0.8, 1.2]): 1. `ladder` "The cascade" (row 0): a log length axis from 10 L down to 10⁻⁷ m; tiers as circles
  shrinking by 2, each labelled (on wide screens) with u′ and turnover time; L (purple ▼), λ_T (teal ▼), η (rose ▼); an
  orange arrow "ε̄ = ΔU³/L" running through every tier; a grey bar "model grid 1 km" when the notes are on. 2. `spec`
  "Spectrum" (row 1, flex 1.4; log–log): the model spectrum bold teal, the −5/3 line orange dashed, the inertial range
  shaded, a ghost of the printed +5/3 (muted, toggle), 2π/L and 2π/η marked; mode toggle physical units ↔ Kolmogorov-scaled
  (then all Re fall on one curve and only the left end moves). Pointer: click → inspector. 3. `cases` "Real cases" (row 1,
  flex 0.8, `hidePortrait`): our table (kitchen mixer, wind tunnel, atmospheric boundary layer, ocean thermocline; columns
  Re_L, ε̄, η, decades) with the current row highlighted.
- **Controls:** `dU` "$\Delta U$" (0.01…50 m/s log, default 1) · `L` "$L$" (0.01…10⁴ m log, default 1) · `nu` chips air 1.5 × 10⁻⁵ /
  water 10⁻⁶ m²/s · `scaled` toggle physical / Kolmogorov-scaled · `printed` toggle "show the printed +5/3" (optional).
- **Presets:** "kitchen mixer" {dU: 1, L: 0.05, nu: water} · "wind tunnel" {dU: 2, L: 0.2, nu: air} · "atmospheric boundary layer"
  {dU: 5, L: 1000, nu: air} · "ocean thermocline" {dU: 0.1, L: 10, nu: water}.
- **Status:** "Re_L = 1.0 × 10⁶: η/L = 3.2 × 10⁻⁵ — 4.5 decades below L" · Re_L < 10³: "⚠️ Re_L too low for an inertial range:
  L and η are less than 2.3 decades apart".
- **Readouts:** "Re_L" · "ε̄" [m²/s³] · "η" [m] · "u_K" · "τ_η" · "λ_T/L" · "DNS points".
- **Depth features:** Explain + Code + Derivation · **linked views** · **presets** · **status** · **inspector** (click the
  spectrum at k₁ = 100: "ε̄^{2/3} k₁^{−5/3}: (m²/s³)^{2/3} × (1/m)^{−5/3} = m^{4/3 + 5/3} s⁻² = m³/s² ✓; value 0.245 × 1 × 4.64 × 10⁻⁴ =
  1.14 × 10⁻⁴ m³/s². The printed k₁^{+5/3} would have units m^{−1/3}s⁻²") · **notes** ("A climate model with a 1 km grid sits
  `…` decades above η: everything to the right of the grey bar must be parameterised") · **modes** (physical / scaled).
- **Explain:**
  0. *What the views show* — "Top: eddies of every size between the big ones (<b class=c-accent>L</b>) and the smallest
     (<b class=c-rose>η</b>); the <b class=c-orange>orange</b> arrow is the energy handed down, the same through every tier.
     Bottom: how the energy is shared among sizes; the straight orange piece is the −5/3 law."
  1. *The supply* — "ε̄ ~ (ΔU)³/L (12.49) = 1³/1 = **1.0 m²/s³** — no viscosity in it."
  2. *The smallest eddies* — "η = (ν³/ε̄)^{1/4} (12.50) = (10⁻¹⁸/1)^{1/4} = **3.2 × 10⁻⁵ m**; u_K = (νε̄)^{1/4} = **0.032 m/s**; τ_η =
     (ν/ε̄)^{1/2} = 1.0 ms; check ηu_K/ν = **1**."
  3. *The separation* — "Re_L = ΔU L/ν = 10⁶; η/L ~ Re_L^{−3/4} (12.51) = **3.2 × 10⁻⁵**: ¾ log₁₀Re = 4.5 decades; λ_T/L ∝ Re_L^{−1/2}
     (12.52) = √(15/Re) = 3.9 × 10⁻³."
  4. *The spectrum* — "between 2π/L and 2π/η: S₁₁ = C₁ε̄^{2/3}k₁^{−5/3} (12.54, exponent corrected); in Kolmogorov units
     S₁₁/(u_K²η) = Φ(k₁η) (12.53) — one curve for every flow."
  5. *The cost of resolving it all* — "(L/η)³ = Re_L^{9/4} = **3 × 10¹³** grid points."
  6. *Real cases* (or the hint "turn the phone sideways for the table").
  7. *Reading the current setting* — raise ΔU or L: "More supply or bigger eddies: Re_L rises and η retreats further from L.
     The cascade gets longer; the −5/3 line gets longer; its level is set by ε̄ alone." · change ν only: "The supply does not
     notice. Only the right-hand end moves: thicker fluid stops the cascade sooner. This is why the dissipation *rate* is
     independent of viscosity even though viscosity does the dissipating." · low Re: "L and η are too close: there is no
     range that is free of both the geometry and viscosity, and no −5/3 line." · atmosphere preset: "6.4 decades: millimetre
     eddies under kilometre gusts. No computer resolves this; sub-grid mixing is modelled (C12)."
- **Derivation tab:** **D11** (9 steps) `view: 'ladder'`; step 3 `live` "ε̄ = ΔU³/L = …"; step 6 `live` "η = (ν³/ε̄)^{1/4} = …";
  step 9 `live` "η/L = Re^{−3/4} = …", `watch` "switch air ↔ water: only the rose marker moves". **D12** (9 steps) `view:
  'spec'`; step 5 `live` "units: m³/s² = (m²/s³)^a (1/m)^b ⇒ a = 2/3, b = −5/3"; step 6 `set` {printed: true}, `watch` "the muted
  +5/3 ghost climbs forever: infinite variance"; step 8 `set` {scaled: true}.
- **Code:**
  ```python
  eps = ch12.dissipation_outer_scaling({{dU}}, {{L}})    # (12.49): {{eps}} m2/s3
  eta, uK, tau_eta = ch12.kolmogorov_scales({{nu}}, eps) # (12.50): {{eta}} m
  Re = {{dU}} * {{L}} / {{nu}}                            # Re_L = {{Re}}
  sep = ch12.scale_separation(Re)["eta_over_L"]          # (12.51): {{sep}}
  S = ch12.inertial_spectrum_1d(k1, eps)                 # (12.54): C1 eps^(2/3) k1^(-5/3)
  n_dns = ch12.dns_grid_points(Re)                       # Re^(9/4) = {{ndns}}
  ```
- **Walkthrough (7 steps):** 1. "Who sets the dissipation?" — "Viscosity dissipates the energy. So why does the rate not
  depend on viscosity?" · 2. "The supply" — "The big eddies take ε̄ ~ ΔU³/L from the mean flow. Change ν: this number does
  not move." `controls: ['nu']`, `readouts: ['eps']` · 3. "The ladder" — "Each tier strains the next. The same ε̄ passes
  through all of them." `highlight: ['view:ladder']` · 4. "Where it ends" — "Only ν and ε̄ matter down there: η =
  (ν³/ε̄)^{1/4}." `derive: {id: 'D11', step: 6}` · 5. "More Reynolds number" — "Raise ΔU: η/L ~ Re^{−3/4}. The cascade gets
  longer." `controls: ['dU']`, `eq: 'sep'` · 6. "The −5/3 law" — "Between L and η only ε̄ and k₁ matter: S₁₁ = C₁ε̄^{2/3}k₁^{−5/3}.
  Click the line." `inspect: true`, `derive: {id: 'D12', step: 5}` · 7. "Your turn" — "Load the atmosphere. Predict the
  decades between L and η, then read them." `set` {dU: 5, L: 1000, nu: 1.5e-5}, `notes: true`.
- **Equations:** `sup` ref 'Eq. (12.49)' ε̄ ~ (ΔU)³/L, live · `kol` ref 'Eq. (12.50)' η = (ν³/ε̄)^{1/4}, u_K = (νε̄)^{1/4}, live · `sep` ref
  'Eq. (12.51)' η/L ~ Re_L^{−3/4}, live · `tay` ref 'Eq. (12.52)' λ_T/L ∝ Re_L^{−1/2} · `law` ref 'Eq. (12.54), exponent corrected' S₁₁(k₁) =
  C₁ε̄^{2/3}k₁^{−5/3} for 2π/L ≪ k₁ ≪ 2π/η, note "the page prints +5/3 (slip #1)" · `univ` ref 'Eq. (12.53)' S₁₁(k₁)/(u_K²η) = Φ(k₁η).
- **Check yourself:** (1) "Switch air → water at the same ΔU and L. Which of ε̄, η, Re_L change?" — "η (smaller) and Re_L
  (larger); ε̄ does not." · (2) "Double ΔU: by what factor does η change?" — "ε̄ × 8 ⇒ η × 8^{−1/4} = 0.59." · (3) "How many decades
  of inertial range at Re_L = 10⁸?" — "¾ × 8 = 6 between L and η." · (4) "Why can the printed k₁^{+5/3} not be right?" — "Units
  (m^{−1/3}s⁻² instead of m³/s²) and an infinite variance."
- **Selftest parity rows:** `{name: 'eta', js: kolmogorov(1e-6, 1)[0], py: 'ch12.kolmogorov_scales(1e-6, 1.0)[0]', rtol: 1e-12}`
  (3.1623e-5) · `{name: 'u_K', js: kolmogorov(1.5e-5, 0.125)[1], py: 'ch12.kolmogorov_scales(1.5e-5, 0.125)[1]', rtol: 1e-12}` (0.037004) ·
  `{name: 'eta/L', js: separation(1e7), py: 'ch12.scale_separation(1e7)["eta_over_L"]', rtol: 1e-12}` (5.6234e-6) · `{name: 'S11
  inertial', js: S11inertial(100, 1, false), py: 'ch12.inertial_spectrum_1d(100.0, 1.0)', rtol: 1e-10}` (1.1393e-4) · `{name: 'model
  spectrum', js: Epao(50, 1, 1e-6, 1), py: 'ch12.model_spectrum(50.0, 1.0, 1e-6, L=1.0, kind="pope")', rtol: 1e-8}` · invariant `{name:
  'eta uK / nu', js: kolmogorov(1e-6, 3)[0]*kolmogorov(1e-6, 3)[1]/1e-6, expect: 1, rtol: 1e-12}`.
- **Fit plan:** 360×640: status (2 lines → shortened to "Re 1e6 · η/L 3.2e-5 · 4.5 dec") · `ladder` (35 %) over `spec` (65 %);
  tier labels hidden, only L, λ_T, η markers; `cases` hidden (the current case name in the status). Desktop: rows [0.8, 1.2].

### E5 · turbulent_energy_budget
- **Title:** "Where does turbulent energy come from and go?" · **Summary:** "Click a height in a channel: the mean flow's
  loss and the turbulence's gain are the same bar with opposite signs." · **CORE:** C06 (also N78, N80, N81, N136, N149, N181) ·
  **Reference:** `fid_formula_lab.html` (clickable terms and bars that add up) + `forced_damped_vibrations.html` (Explain).
- **meta:** `viz:order 5` · `viz:sections 12.7 12.9 12.10` · `viz:equations 12.46 12.47` · `viz:fluidpy ch12.channel_energy_budget
  ch12.channel_energy_budget_at ch12.mean_energy_budget ch12.tke_budget ch12.shear_production ch12.channel_mixing_length
  ch12.mixing_length_wall_profile` · `viz:derivations D09 D10`.
- **Physics (JS ↔ Python):** inner-layer model in wall units: l⁺ = κy⁺(1 − e^{−y⁺/A⁺}) (capped at 0.09 Re_τ in the core), total
  stress τ⁺ = 1 − y⁺/Re_τ, slope s = dU⁺/dy⁺ = 2τ⁺/(1 + √(1 + 4l⁺²τ⁺)), −⟨uv⟩⁺ = τ⁺ − s; terms at a height (all in units u_*⁴/ν):
  pressure work W = U⁺/Re_τ; mean viscous dissipation s²; production P = (τ⁺ − s)s (the loss of the mean flow and the gain of
  the turbulence); viscous transport of mean energy d(U⁺s)/dy⁺ and Reynolds-stress transport as the remainder;
  `budgetAt(yp, Ret, kappa, Ap)` ↔ `ch12.channel_energy_budget_at(yplus, Re_tau, kappa, A_plus)`; the profile table for the
  curves from `reference/ch12/explainer_tables.json` (Re_τ = 180, 550, 1000, 5200; "ours, computed by fluidpy") ↔
  `ch12.channel_energy_budget`; integrals ∫W dy⁺ = ∫s² dy⁺ + ∫P dy⁺.
- **Views** (rows [1, 1]): 1. `prof` "The channel" (row 0, flex 1.2; semi-log y⁺): U⁺ (purple), −⟨uv⟩⁺ (orange), viscous stress s
  (rose), production P ×4 (orange filled), the cursor line; the height where s = −⟨uv⟩⁺ marked ◆ ("production peak").
  Pointer: click or drag → cursor. 2. `bars` "Two budgets at this height" (row 0, flex 1): left group **mean flow**
  (pressure work + purple; transport in/out grey; direct dissipation − rose; loss to turbulence − orange), right group
  **turbulence** (production + orange — same length as the bar to its left, joined by a bracket; dissipation and transport,
  modelled, − rose hatched, labelled "model"). 3. `whole` "Across the whole channel" (row 1, `hidePortrait`): a stacked bar
  work by the pressure gradient = direct dissipation + production, with percentages.
- **Controls:** `Ret` "$\mathrm{Re}_\tau$" (180, 550, 1000, 5200 chips) · `yp` "Height $y^+$" (0.5…Re_τ log, default 12) · `damp`
  toggle wall damping (A⁺ = 26 / none) · `axis` toggle wall units / outer units (optional).
- **Presets:** "sublayer" {yp: 2} · "buffer-layer peak" {yp: 11} · "log layer" {yp: 100} · "centreline" {yp: Re_τ}.
- **Status:** "y⁺ ≈ 10.4: production at its maximum, ¼ in wall units — viscous and Reynolds stress are equal" · "y⁺ = 2: sublayer
  — the mean flow loses energy directly to viscosity; almost no production" · "centreline: no shear, no production; turbulence
  here lives on transport".
- **Readouts:** "Production P⁺" · "Direct dissipation s²" · "−⟨uv⟩⁺" · "dU⁺/dy⁺" · "P / total loss".
- **Depth features:** Explain + Code + Derivation · **linked views** · **terms** (two bar groups that share one mirrored term;
  click a term to light its curve) · **inspector** (click a height: "y⁺ = 100, Re_τ = 1000: τ⁺ = 0.90; l⁺ = 40.1; s = 2τ⁺/(1 + √(1 + 4l⁺²τ⁺)) = 0.0233; −⟨uv⟩⁺ = 0.877; P = 0.877 × 0.0233 = 0.0205; direct dissipation s² = 5.4 × 10⁻⁴ — 2.7 % of P") · **presets** · **status**.
- **Explain:**
  0. *What the views show* — "Left: the mean velocity (<b class=c-accent>purple</b>), the two parts of the shear stress —
     <b class=c-rose>viscous</b> and <b class=c-orange>Reynolds</b> — and the production. Right: at the cursor height, what the
     mean flow gains and loses, and what the turbulence gains and loses. **Exact:** the stresses and the production of this
     model profile. **Modelled:** how the turbulence disposes of it (hatched)."
  1. *The stresses here* — "total τ̄⁺ = 1 − y⁺/Re_τ = …; viscous dU⁺/dy⁺ = …; Reynolds −⟨uv⟩⁺ = the rest = …".
  2. *The shared term* — "for U(y): production = −⟨uv⟩ dU/dy = **…** (in u_*⁴/ν). It is +⟨u_iu_j⟩∂U_i/∂x_j — a loss — in the
     mean-flow budget ∂Ē/∂t + U_j∂Ē/∂x_j = ∂/∂x_j(…) − 2νS̄_ijS̄_ij + ⟨u_iu_j⟩∂U_i/∂x_j − (g/ρ₀)ρ̄U₃ (12.46), and −⟨u_iu_j⟩∂U_i/∂x_j — a gain — in the turbulent budget ∂ē/∂t + U_j∂ē/∂x_j = ∂/∂x_j(…) − 2ν⟨S′_ijS′_ij⟩ − ⟨u_iu_j⟩∂U_i/∂x_j + gα⟨u₃T′⟩ (12.47)."
  3. *Where it peaks* — "P = (τ⁺ − s)s is largest where s = τ⁺/2: viscous stress = Reynolds stress, P_max = τ⁺²/4 ≈ **0.25** at
     y⁺ ≈ 10.4 (this model; DNS: y⁺ ≈ 12)."
  4. *The mean flow's direct loss* — "2νS̄_ijS̄_ij → (dU⁺/dy⁺)² = …; ratio to production … — of order 1/Re away from the wall
     (N80), but it wins inside the sublayer."
  5. *The whole channel* (or hint) — "work by the pressure gradient = direct dissipation + production: at Re_τ = 1000 …% + …%."
  6. *Reading the current setting* — sublayer: "The wall has killed the fluctuations: the mean flow pays viscosity directly."
     · buffer: "The busiest place in the flow: the mean shear is still large and the Reynolds stress has just grown to match
     the viscous stress. A large share of all production happens in this thin layer." · log layer: "Production ≈ u_*³/(κy) falls with
     height; nearly all of it is dissipated locally — the equilibrium the k–ε constants are tuned to." · centreline: "No mean
     shear, no production: the turbulence here was made elsewhere and carried in by the transport terms."
- **Derivation tab:** **D09** (9 steps) `view: 'bars'`; step 5 `watch` "the grey transport bars: they move energy in y but
  integrate to zero across the channel"; step 8 `live` "+⟨uv⟩dU/dy = −P = …". **D10** (15 steps, all shown) `view: 'bars'`; step 8
  `live` "−⟨uv⟩dU/dy = +…" with `highlight: ['term:production']`; step 13 `watch` "the hatched bar: the positive 2ν⟨S′S′⟩ is
  what the turbulence loses"; step 15 `set` {yp: 11}. Interpret: section 6.
- **Code:**
  ```python
  b = ch12.channel_energy_budget_at({{yp}}, {{Ret}}, 0.41, {{Ap}})
  P = b["production"]            # -<uv> dU/dy = {{P}}  (wall units)
  D = b["mean_dissipation"]      # (dU+/dy+)^2 = {{D}}
  W = b["pressure_work"]         # U+/Re_tau = {{W}}
  # mean flow loses P (12.46); turbulence gains P (12.47)
  whole = ch12.channel_energy_budget({{Ret}}, 0.41, {{Ap}})["integrals"]
  ```
- **Walkthrough (6 steps):** 1. "Who feeds the turbulence?" — "Turbulence dies without a supply. Click a height and look at
  the bars." `inspect: true` · 2. "One term, two signs" — "The orange bar is a loss on the left and the same gain on the
  right: −⟨uv⟩ dU/dy." `terms: true`, `eq: 'tke'` · 3. "Where it peaks" — "Move to y⁺ ≈ 10: viscous and Reynolds stress are
  equal and production is ¼." `set` {yp: 10.4} · 4. "The sublayer" — "At y⁺ = 2 the mean flow pays viscosity directly."
  `set` {yp: 2} · 5. "Why the term appears twice" — "Multiply the fluctuation equation by u_i and average: the same product
  comes out with a minus sign." `derive: {id: 'D10', step: 8}` · 6. "Your turn" — "Predict what the centreline bars look
  like; then go there." `controls: ['yp', 'Ret']`.
- **Equations:** `mean` ref 'Eq. (12.46)' (aligned, three lines) ∂Ē/∂t + U_j∂Ē/∂x_j = ∂/∂x_j(−U_jP/ρ₀ + 2νU_iS̄_ij − ⟨u_iu_j⟩U_i) − 2νS̄_ijS̄_ij
  + ⟨u_iu_j⟩∂U_i/∂x_j − (g/ρ₀)ρ̄U₃ · `tke` ref 'Eq. (12.47)' (aligned) ∂ē/∂t + U_j∂ē/∂x_j = ∂/∂x_j(−⟨pu_j⟩/ρ₀ + 2ν⟨u_iS′_ij⟩ − ½⟨u_i²u_j⟩) −
  2ν⟨S′_ijS′_ij⟩ − ⟨u_iu_j⟩∂U_i/∂x_j + gα⟨u₃T′⟩ · `prod` ref 'special case U(y)' P = −⟨uv⟩ dU/dy, live · `ratio` ref 'p. 565 (no number)'
  2νS̄_ijS̄_ij/[⟨u_iu_j⟩∂U_i/∂x_j] ~ 1/Re, live.
- **Check yourself:** (1) "At what height are the viscous and Reynolds stresses equal, and what is special there?" — "y⁺ ≈ 10:
  production is at its maximum." · (2) "Is production ever negative in this channel?" — "No: −⟨uv⟩ and dU/dy have the same
  sign everywhere (C04's parcel argument)." · (3) "Raise Re_τ from 180 to 5200: what happens to the share of direct
  dissipation in the whole-channel budget?" — "It falls: more of the pressure work goes through the turbulence." · (4) "Turn
  the damping off: where does the peak go?" — "Closer to the wall: s = ½ needs κ²y⁺² = 2, i.e. y⁺ = √2/κ ≈ 3.4."
- **Selftest parity rows:** `{name: 'production y+=12', js: budgetAt(12, 1000, 0.41, 26).production, py:
  'ch12.channel_energy_budget_at(12.0, 1000.0, 0.41, 26.0)["production"]', rtol: 1e-6}` · `{name: 'mean dissipation y+=100', js:
  budgetAt(100, 1000, 0.41, 26).mean_dissipation, py: 'ch12.channel_energy_budget_at(100.0, 1000.0, 0.41, 26.0)["mean_dissipation"]',
  rtol: 1e-6}` · `{name: 'U+ at 30', js: budgetAt(30, 1000, 0.41, 26).Uplus, py: 'ch12.channel_energy_budget_at(30.0, 1000.0, 0.41,
  26.0)["Uplus"]', rtol: 1e-4}` · invariant `{name: 'peak production', js: peakProduction(5200, 0.41, 26), expect: 0.25, rtol: 0.01}`.
- **Fit plan:** 360×640: `prof` (45 %) over `bars` (55 %); `whole` hidden (its two percentages in the Explain tab and a
  readout); bar labels shortened ("prod.", "diss."). 844×390: `prof` | `bars`. Desktop: `whole` under both.

### E6 · turbulent_jet_similarity
- **Title:** "How does a jet spread without a turbulence model?" · **Summary:** "One conserved flux and one shape: toggle raw ↔
  rescaled profiles, break the exponents, switch to wakes and plumes." · **CORE:** C09 (also N101, N109, N116, N122–N125, N130) ·
  **Reference:** `angular_frequency_explorer_1.html` (modes and linked views on one state) with the ch09 `free_jet_similarity`
  stage (raw ↔ rescaled, invariant bars).
- **meta:** `viz:order 6` · `viz:sections 12.8` · `viz:equations 12.56 12.62 12.63 12.65 12.66 12.68` · `viz:fluidpy
  ch12.plane_jet_mean_velocity ch12.jet_momentum_flux_per_span ch12.plane_jet_volume_flux ch12.free_shear_exponents
  ch12.wrong_exponent_fluxes ch12.profile_integrals ch12.general_similarity_check` · `viz:derivations D13 D14 D15 D16`.
- **Physics (JS ↔ Python):** `F(xi, xh)` = exp(−ln2 ξ²/ξ½²) ↔ `ch12.gaussian_profile`; `integrals(xh)` → [I₁ = ξ½√(π/ln2), I₂ =
  ξ½√(π/2ln2)] ↔ `ch12.profile_integrals(xh)["I1"]`, `["I2"]`; `jetU(x, y, JsRho, xh)` = I₂^{−1/2}(J_s/ρ)^{1/2}x^{−1/2}F(y/x) ↔
  `ch12.plane_jet_mean_velocity(x, y, Js, rho, C5="from_invariant", xi_half=xh)`; `volFlux(x, JsRho, xh)` ↔
  `ch12.plane_jet_volume_flux`; `exponents(flow)` → {width, velocity, scalar, reynolds} (exact fractions as [num, den]) ↔
  `ch12.free_shear_exponents(flow)`; `wrongFluxes(x, n, m)` = [x^{2n+m}, x^{n+m}] ↔ `ch12.wrong_exponent_fluxes`; generic flow
  profiles U_CL ∝ x^{a}, δ ∝ x^{b} with a Gaussian in the flow's similarity variable ↔ `ch12.free_shear_profile`; the laminar
  ghost (δ ∝ x^{2/3}, U ∝ x^{−1/3}) ↔ `core.jets.free_jet` (shape only; no parity).
- **Views** (rows [1, 1]): 1. `field` "The flow" (row 0, full width): x 0…10 (units of d × 20), y ±3; colour = mean speed
  (purple scale); half-width lines ±ξ½x (white); entrainment arrows at the edges with length ∝ ½ dV̇/dx; the laminar jet's
  half-width as a muted dashed ghost; five station lines, the selected one bold. For wakes the colour is the deficit; for a
  plume the picture is turned upright. Pointer: drag the bold station. 2. `prof` "Profiles at five stations" (row 1, flex
  1.2): raw U(y) (five curves, spreading and sinking) or rescaled U/U_CL vs y/δ (one curve when the exponents are right;
  fanned out when the trial exponent is wrong). 3. `inv` "What is conserved" (row 1, flex 0.8): bars along x for the flow's
  invariant (flat), the volume flux (rising) and the scalar flux (flat, amber), each normalised by its value at the first
  station; with a wrong exponent the invariant bar slopes and turns rose.
- **Controls:** `flow` chips plane jet / round jet / plane wake / plane plume (mode) · `resc` toggle raw / rescaled · `n` "Trial
  decay exponent $n$" (−1.5…0, step 0.05; default = the correct one for the flow) · `xs` "Station $x$" (1…10) · `scalar`
  toggle (optional).
- **Presets:** "correct exponents" {n: correct} · "wrong decay" {flow: plane jet, n: −0.4} · "laminar jet (Ch. 9)" {ghost on,
  resc: false} · "round jet" · "plane wake".
- **Status:** "✅ plane jet: δ ∝ x, U_CL ∝ x^{−1/2}: momentum flux flat, volume flux ∝ x^{1/2}" · "⚠️ n = −0.4: momentum flux
  grows as x^{0.2} — not allowed by J_s = ρ∫U²dy = const (12.62)" · "round jet: U_CL ∝ x^{−1}; local Reynolds number constant".
- **Readouts:** "U_CL at x" · "half-width" · "J_s(x)/J_s(1)" · "V̇(x)/V̇(1)" · "entrainment v_e" · "Re exponent".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (four flows) · **presets** · **status** ·
  **terms** (the three invariant bars; click one to see its integral written out).
- **Explain:**
  0. *What the views show* — "Top: the mean flow (<b class=c-accent>purple</b> = fast); white lines mark where the speed is
     half the centreline value; arrows are the still fluid being drawn in. Bottom left: velocity profiles at five distances.
     Bottom right: quantities integrated across the flow at each distance."
  1. *The invariant* — "J_s = ρ∫U²dy = ρU_CL²δ∫F²dξ (12.65): constant ⇒ U_CL²δ = const. With U_CL ∝ x^n, δ ∝ x^m: **2n + m = 0**."
  2. *The growth law* — "shape-keeping in the similarity equation {δU′_CL/U_CL}F² − {δU′_CL/U_CL + δ′}F′∫₀^ξF dξ = {Ψ/U_CL²}G′ (12.63)
     needs δ′ = const: **m = 1**. (Wakes: dδ/dx ~ ΔU/U_∞ instead ⇒ m − 1 = n.)"
  3. *Solve the two equations* — "plane jet: m = 1, n = −½ ⇒ U(x, y) = C₅(J_s/ρ)^{1/2}x^{−1/2}F(y/x) (12.66)." Table of the four
     flows (invariant · growth law · m · n) with the current row lit.
  4. *At your station* — "x = 4: U_CL = 2.577/√4 = **1.289**; half-width 0.40; V̇ ∝ x^{1/2} (12.68) = **2.0** × its value at x = 1;
     v_e = ½ dV̇/dx = …".
  5. *With your trial exponent* — "n = −0.4, m = 1: momentum flux ∝ x^{2n+m} = x^{0.2}; volume flux ∝ x^{n+m} = x^{0.6}."
  6. *Local Reynolds number* — "U_CLδ/ν ∝ x^{n+m}: plane jet +½ (grows — it stays turbulent), round jet 0, plane wake 0, round
     wake −⅓ (falls)."
  7. *Reading the current setting* — correct: "No turbulence model was used: conservation fixed one exponent and
     shape-keeping the other. The constants (the spreading rate ξ½, C₅) are the only things experiments must supply." ·
     wrong n: "The profiles still look plausible at one station — but the bar shows the jet creating momentum from nothing.
     The collapse and the invariant must hold together." · wake: "A wake's deficit is carried at the free-stream speed, so
     it spreads more slowly (√x) than a jet (x)." · plume: "Buoyancy keeps adding momentum: the centreline speed does not
     decay at all, while the scalar is diluted as 1/x."
- **Derivation tab:** **D13** (10 steps) `view: 'inv'`; step 9 `watch` "the momentum bar is flat at every station". **D14** (14
  steps, all shown; long lines split with `aligned`) `view: 'prof'`, `set` {resc: true}; step 4 `live` "∂U/∂x = U′_CL F −
  U_CL F′ ξ δ′/δ at ξ = 0.1, x = 2: …"; step 11 `watch` "the two ξFF′ terms are equal and opposite"; step 14 `live` "the three
  brackets with δ = x, U_CL ∝ x^{−1/2}: −½, +½, C₃". **D15** (10 steps) `view: 'field'`; step 6 `set` {n: −0.4} then step 7 `set`
  {n: −0.5}, `watch` "the rose bar turns flat exactly at n = −½"; step 10 `live` "V̇(4)/V̇(1) = 2". **D16** (12 steps) `view:
  'inv'`; steps 3, 6, 9, 11 `set` {flow: round jet / plane wake / plane plume / …}; `live` the two equations for the selected
  flow from `free_shear_exponents(flow, return_equations=True)`.
- **Code:**
  ```python
  ex = ch12.free_shear_exponents("{{flow}}")       # width {{m}}, velocity {{n}}
  U = ch12.plane_jet_mean_velocity({{x}}, y, 1.0, 1.0, C5="from_invariant",
                                   xi_half=0.10)   # (12.66): U_CL = {{ucl}}
  J = ch12.jet_momentum_flux_per_span(y, U, 1.0)   # (12.62): {{J}} (same at every x)
  V = ch12.plane_jet_volume_flux({{x}}, 1.0, 1.0, C5="from_invariant",
                                 xi_half=0.10)     # (12.68): {{V}} ~ x^(1/2)
  mom, vol = ch12.wrong_exponent_fluxes({{x}}, {{ntrial}}, 1.0)   # {{mom}}, {{vol}}
  ```
- **Walkthrough (7 steps):** 1. "No model needed?" — "A jet widens and slows. Can we say how fast without knowing the
  turbulence?" · 2. "Something is conserved" — "No force acts on the jet: its momentum flux J_s = ρ∫U²dy is the same at every
  station." `highlight: ['view:inv']`, `eq: 'inv'` · 3. "One shape" — "Toggle 'rescaled': five profiles become one,
  U = U_CL(x)F(y/δ)." `set` {resc: true}, `controls: ['resc']` · 4. "Two equations" — "Conservation: 2n + m = 0. Shape-keeping:
  m = 1. So n = −½." `derive: {id: 'D15', step: 7}` · 5. "Break it" — "Set n = −0.4: the momentum bar climbs. Physics says
  no." `set` {n: -0.4}, `controls: ['n']` · 6. "Entrainment" — "The volume flux grows as x^{1/2}: the jet swallows still
  fluid from the sides." `eq: 'vol'` · 7. "Your turn" — "Predict the exponents of a round jet (invariant ∝ U²δ²), then switch."
  `controls: ['flow']`.
- **Equations:** `sim` ref 'Eq. (12.56)' U(x, y) = U_CL(x)F(y/δ(x)) · `inv` ref 'Eq. (12.62)' J_s ≡ ρ_s∫[U²]_{x=0}dy ≅ ρ∫U²dy = const.,
  live · `simeq` ref 'Eq. (12.63)' {δU′_CL/U_CL}F² − {δU′_CL/U_CL + δ′}F′∫₀^ξ F dξ = {Ψ/U_CL²}G′ (aligned) · `exp` ref 'Eq. (12.65)' J_s =
  ρC₄²x^{2γ+1}∫F²dξ ⇒ 2γ + 1 = 0 · `far` ref 'Eq. (12.66)' U = C₅(J_s/ρ)^{1/2}x^{−1/2}F(y/x), live · `vol` ref 'Eq. (12.68)' V̇(x) =
  C₅(J_s/ρ)^{1/2}x^{+1/2}∫F dξ, live. (The scalar law Ȳ = C₆(Ṁ_s/√(ρJ_s))x^{−1/2}H(y/x) (12.71) appears when `scalar` is on.)
- **Check yourself:** (1) "From x = 1 to x = 4 the centreline speed of a plane jet falls by what factor, and the volume flux
  rises by what factor?" — "Both 2: x^{−1/2} and x^{+1/2}." · (2) "Which bar tells you a trial exponent is wrong?" — "The
  invariant: it must be flat." `set {n: -0.7}` · (3) "Why does a round jet's Reynolds number not change downstream?" — "U_CL ∝
  1/x and δ ∝ x: the product is constant." `set {flow: round jet}` · (4) "A plane wake widens as x^{1/2}: which assumption
  differs from the jet?" — "The deficit is advected at U_∞, so dδ/dx ~ ΔU/U_∞ is not constant."
- **Selftest parity rows:** `{name: 'U_CL x=4', js: jetU(4, 0, 1, 0.1), py: 'ch12.plane_jet_mean_velocity(4.0, 0.0, 1.0, 1.0,
  C5="from_invariant", xi_half=0.1)', rtol: 1e-10}` (1.28868) · `{name: 'volume flux x=1', js: volFlux(1, 1, 0.1), py:
  'ch12.plane_jet_volume_flux(1.0, 1.0, 1.0, C5="from_invariant", xi_half=0.1)', rtol: 1e-10}` (0.548705) · `{name: 'wrong momentum',
  js: wrongFluxes(4, -0.4, 1)[0], py: 'ch12.wrong_exponent_fluxes(4.0, -0.4, 1.0)[0]', rtol: 1e-12}` (1.3195) · `{name: 'round plume
  velocity exponent', js: exponents('round_plume').velocity, py: 'ch12.free_shear_exponents("round_plume")["velocity"]*1.0', rtol:
  1e-12}` (−0.33333) · invariant `{name: 'J flat', js: momentumFlux(7, 1, 0.1)/momentumFlux(1, 1, 0.1), expect: 1, rtol: 1e-6}`.
- **Fit plan:** 360×640: `field` (40 %) over `prof` (60 %); `inv` hidden — its two numbers ("J_s ratio", "V̇ ratio") are
  readouts and the status carries the verdict; chips wrap to two lines. 844×390: `prof` | `inv`, `field` as a strip. Desktop:
  rows [1, 1].

### E7 · law_of_the_wall
- **Title:** "Where does the logarithm come from?" · **Summary:** "Wall units collapse every profile near the wall, outer units
  far from it; where both hold, y dU/dy can depend on neither — slide Re_τ and watch the log layer open." · **CORE:** C10, C11
  (also N136, N137, N139, N142, N144, N147, N149, N153, N154, N156, N159) · **Reference:**
  `amplitude_phase_second_order_II_3.html` (crosshair readouts across linked windows; numbered live derivation).
- **meta:** `viz:order 7` · `viz:sections 12.9` · `viz:equations 12.80 12.81 12.82 12.84 12.87 12.88 12.93` · `viz:fluidpy
  ch12.wall_units ch12.viscous_sublayer ch12.log_law ch12.spalding_uplus ch12.composite_profile_plus ch12.stress_partition
  ch12.layer_name ch12.log_law_indicator ch12.rough_wall_log_law ch12.log_law_crossing ch12.friction_law_from_overlap` ·
  `viz:derivations D17 D18 D19 D20`.
- **Physics (JS ↔ Python):** `spaldingU(yp, kappa, B)` (Newton, C.5-2) ↔ `ch12.spalding_uplus(yplus, kappa=, B=)`; `composite(yp,
  Ret, kappa, B, Pi)` = Spalding + (2Π/κ)W(y⁺/Re_τ), W = 3ξ² − 2ξ³ ↔ `ch12.composite_profile_plus(yplus, Re_tau, kappa=, B=, Pi=)`;
  `logLaw(yp, kappa, B)` ↔ `ch12.log_law`; `crossing(kappa, B)` ↔ `ch12.log_law_crossing`; `partition(yp, Ret, kappa, B)` ↔
  `ch12.stress_partition(…)["viscous"]`, `["reynolds"]`; `layerName(yp, yOverDelta)` ↔ `ch12.layer_name` (exact strings);
  `indicator` = y⁺dU⁺/dy⁺ by differentiating the composite ↔ `ch12.log_law_indicator`; `rough(yOverY0, kappa)` ↔
  `ch12.rough_wall_log_law(y, 1.0, y0, kappa=)`; `frictionLaw(Ret, kappa, A, B)` ↔ `ch12.friction_law_from_overlap`; `wallUnits
  (tau0, rho, nu)` → [u_*, l_ν] ↔ `ch12.friction_velocity`, `ch12.viscous_length`. DNS dots: three short columns from
  `reference/ch12/` (Lee & Moser 2015), embedded as a 30-point table per case with the citation in the view title.
- **Views** (rows [1.3, 1]): 1. `prof` "U⁺ against y⁺" (row 0, full width; semi-log, y⁺ 0.5…10⁵): the composite profile bold
  purple; U⁺ = y⁺ rose dashed; the log line teal dashed; coloured bands sublayer / buffer / logarithmic / wake with their
  names; DNS dots (muted) when Re_τ matches a case; a crosshair at the cursor; in **outer** mode the axes become
  (U_∞ − U)/u_* against y/δ and curves at four Re_τ are overlaid (they collapse at the right, fan at the left). Pointer: drag
  the crosshair. 2. `stress` "Who carries the stress" (row 1, flex 1): total (straight line), viscous (rose), Reynolds
  (orange) against y⁺ (log) with the crossing ◆. 3. `ind` "y⁺ dU⁺/dy⁺" (row 1, flex 1, `hidePortrait`): the indicator with the
  level 1/κ dashed; the plateau shaded where within 5 %.
- **Controls:** `lret` "$\log_{10}\mathrm{Re}_\tau$" (2…6, step 0.05, default 3) · `kappa` "$\kappa$" (0.36…0.43, default 0.41) · `B` "$B$"
  (3.5…6, default 5.0) · `scaling` chips inner / outer (mode) · `wall` chips smooth / rough + `y0p` "$y_0^+$" (rough only;
  optional) · `Pi` "Wake strength $\Pi$" (0…1, default 0.2 — illustrative; optional) · `yp` cursor (pointer).
- **Presets:** "Re_τ = 180" · "Re_τ = 1000" · "Re_τ = 5200" · "atmospheric surface layer (10⁶)" · "rough wall" {wall: rough}.
- **Status** (exact `layer_name` strings): "y⁺ = 2.0: viscous sublayer — U⁺ = y⁺ within 1 %" · "y⁺ = 12: buffer layer — neither
  law holds (linear 12.0, log 11.1, actual 9.2)" · "y⁺ = 300: logarithmic layer" · "y/δ = 0.5: wake region — depends on δ".
- **Readouts:** "U⁺ (profile)" · "U⁺ = y⁺" · "log law" · "viscous share" · "log layer [decades]" · "C_f".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** · **modes** (inner / outer
  scaling) · **inspector** (crosshair at y⁺ = 100: "U⁺ = (1/κ) ln y⁺ + B = ln(100)/0.41 + 5.0 = 11.23 + 5.0 = **16.23**; Spalding
  16.08; viscous share of the stress 2.4 %").
- **Explain:**
  0. *What the views show* — "Top: mean velocity against distance from the wall, both in wall units, on a logarithmic
     distance axis. <b class=c-rose>Rose dashed</b>: the sublayer law. <b class=c-teal>Teal dashed</b>: the logarithmic law.
     Bottom left: how the (linear) total stress is shared between viscosity and turbulence. Bottom right (hidden on phones):
     y⁺dU⁺/dy⁺, which is flat where the profile is a logarithm."
  1. *Wall units* — "u_*² = τ₀/ρ (12.81): for τ₀ = 0.3 Pa in air u_* = 0.5 m/s, l_ν = ν/u_* = 30 µm; U⁺ = U/u_* = f(y⁺) (12.80);
     δ⁺ = Re_τ = **1000**."
  2. *At your height* — "sublayer law U⁺ = y⁺ (12.82) = …; log law U⁺ = (1/κ) ln y⁺ + B (12.88) = …; profile = …" (three boxed
     numbers).
  3. *Where the two lines cross* — "y⁺ = (1/κ) ln y⁺ + B ⇒ y⁺ = **10.8** (κ = 0.41, B = 5.0): the middle of the buffer layer."
  4. *How long is the log layer* — "from y⁺ ≈ 30 to y/δ ≈ 0.15: log₁₀(0.15 Re_τ/30) = **0.70** decades at Re_τ = 1000; 2.7 at
     10⁵. Below Re_τ ≈ 200 there is none."
  5. *The stress* — "τ̄ = τ₀(1 − 2y/h): viscous share dU⁺/dy⁺ = …, Reynolds share = …; equal near y⁺ ≈ 11."
  6. *The friction law* — "adding (U_∞ − U)/u_* = −(1/κ) ln ξ + A (12.84 with (12.89)) to the log law: U_∞⁺ = (1/κ) ln δ⁺ + A + B =
     …; C_f = 2/(U_∞⁺)² = …".
  7. *The indicator* (or hint) — "y⁺dU⁺/dy⁺ = 1/κ = **2.44** on the plateau (12.87)."
  8. *Reading the current setting* — low Re_τ: "The wall layer and the outer layer overlap nowhere: viscosity is felt across
     most of the flow and there is no logarithm to speak of." · high Re_τ: "A wide range where the flow has forgotten the
     viscous length and does not yet feel δ. There only y is left, so dU/dy = u_*/(κy)." · outer mode: "Divide by δ instead:
     now the wakes agree and the wall regions disagree. The log layer is the only place that looks the same in both
     pictures — that is what 'overlap' means." · rough: "Roughness taller than the sublayer removes viscosity from the
     problem: U⁺ = (1/κ) ln(y/y₀) (12.93); the drag coefficient no longer depends on Reynolds number. This is the wind
     profile over land and sea."
- **Derivation tab:** **D17** (9 steps) `view: 'stress'`; step 7 `live` "τ̄/τ₀ = 1 − 2y/h = … at your height". **D18** (7 steps)
  `view: 'prof'`; step 4 `live` "u_* = √(τ₀/ρ) = 0.5 m/s, y⁺ = y u_*/ν"; step 7 `set` {yp: 3}, `watch` "the profile sits on the
  rose line up to y⁺ ≈ 5". **D19** (12 steps) `view: 'prof'` (phones) / `'ind'` (wide); step 3 `set` {scaling: outer}; step 7
  `set` {scaling: inner, lret: 4}, `live` "y⁺ df/dy⁺ = … = −ξ dF/dξ"; step 8 `watch` "the plateau at 1/κ: raise Re_τ and it
  widens"; step 12 `live` "U_∞⁺ = (1/κ) ln δ⁺ + A + B = …". **D20** (5 steps) `view: 'prof'`, `set` {wall: rough}; step 3 `live`
  "U⁺ = (1/κ) ln(y/y₀)". Interpret: section 8.
- **Code:**
  ```python
  u_star = WT.friction_velocity(0.3, 1.2)             # (12.81): 0.5 m/s
  yp = {{yp}}                                          # y+ = y u*/nu
  U_sub = WT.viscous_sublayer(yp)                      # (12.82): {{usub}}
  U_log = WT.log_law(yp, kappa={{kappa}}, B={{B}})     # (12.88): {{ulog}}
  U_in = WT.spalding_uplus(yp, kappa={{kappa}}, B={{B}})   # {{uin}}
  name = WT.layer_name(yp, yp/{{Ret}})                 # "{{layer}}"
  Uinf = WT.friction_law_from_overlap({{Ret}}, kappa={{kappa}}, A=1.0, B={{B}})  # {{uinf}}
  ```
- **Walkthrough (8 steps):** 1. "One curve for every wall?" — "Different flows, different Reynolds numbers — the same curve
  near the wall. Why?" · 2. "Wall units" — "Only τ₀ and ν matter there: u_* = √(τ₀/ρ), l_ν = ν/u_*. Then U⁺ = f(y⁺)." `eq:
  'wall'`, `derive: {id: 'D18', step: 4}` · 3. "The sublayer" — "Below y⁺ ≈ 5 the stress is all viscous: U⁺ = y⁺." `set` {yp: 3},
  `highlight: ['view:stress']` · 4. "Raise Re_τ" — "The inner curve does not move; a straight piece grows between y⁺ ≈ 30 and
  0.15 δ⁺." `controls: ['lret']` · 5. "The other collapse" — "Switch to outer scaling: now the wakes agree and the wall
  regions do not." `set` {scaling: 'outer'} · 6. "Where both hold" — "There y dU/dy can depend on neither length: it is the
  constant u_*/κ." `derive: {id: 'D19', step: 8}` · 7. "A logarithm" — "Integrate dU/dy = u_*/(κy): U⁺ = (1/κ) ln y⁺ + B. Move
  the crosshair to y⁺ = 100." `inspect: true`, `eq: 'log'` · 8. "Your turn" — "Predict how many decades of log layer the
  atmosphere (Re_τ = 10⁶) has; then load it."
- **Equations:** `wall` ref 'Eq. (12.80)' U⁺ ≡ U/u_* = f(yu_*/ν) = f(y⁺) · `ustar` ref 'Eq. (12.81)' u_*² ≡ τ₀/ρ, live · `sub` ref
  'Eq. (12.82)' U⁺ = y⁺, live · `def` ref 'Eq. (12.84)' (U_∞ − U)/u_* = F(y/δ) · `match` ref 'Eq. (12.87)' −ξ dF/dξ = y⁺ df/dy⁺ (= 1/κ), live ·
  `log` ref 'Eq. (12.88)' U⁺ = (1/κ) ln(y⁺) + B, live · `rough` ref 'Eq. (12.93)' U⁺ = (1/κ) ln(y/y₀). (The linear stress τ̄ = τ₀(1 −
  2y/h) from 0 = −∂P/∂x + ∂τ̄/∂y (12.76) is in D17.)
- **Check yourself:** (1) "At y⁺ = 12, which law is right?" — "Neither: buffer layer; the profile lies below both lines." `set
  {yp: 12}` · (2) "Changing κ from 0.41 to 0.38 does what to the line?" — "Steepens it: slope per decade 2.303/κ rises from
  5.6 to 6.1." · (3) "At Re_τ = 180, how wide is the plateau of y⁺dU⁺/dy⁺?" — "There is none: 0.15 × 180 = 27 < 30." `set {lret:
  2.255}` · (4) "Why does a rough-wall law contain no viscosity?" — "The elements shed eddies directly; the sublayer — the
  only place viscosity mattered — is gone."
- **Selftest parity rows:** `{name: 'log law 100', js: logLaw(100, 0.41, 5.0), py: 'ch12.log_law(100.0, kappa=0.41, B=5.0)', rtol:
  1e-12}` (16.2321) · `{name: 'Spalding 12', js: spaldingU(12, 0.41, 5.0), py: 'ch12.spalding_uplus(12.0, kappa=0.41, B=5.0)', rtol:
  1e-8}` (9.15694) · `{name: 'crossing', js: crossing(0.41, 5.0), py: 'ch12.log_law_crossing(kappa=0.41, B=5.0)', rtol: 1e-8}` (10.8049) ·
  `{name: 'viscous share', js: partition(30, 1000, 0.41, 5.0)[0], py: 'ch12.stress_partition(30.0, 1000.0, kappa=0.41,
  B=5.0)["viscous"]', rtol: 1e-6}` · `{name: 'layer name', js: layerName(12, 0.012), py: 'ch12.layer_name(12.0, 0.012)', exact: true}`.
- **Fit plan:** 360×640: status (1 line) · `prof` (60 %) over `stress` (40 %); `ind` hidden — the readout "y⁺dU⁺/dy⁺" and the
  plateau value 1/κ sit in the `prof` title; band names shortened (sub · buf · log · wake). 844×390: `prof` | `stress`.
  Desktop: rows [1.3, 1].

### E8 · mixing_length_closure
- **Title:** "What does a closure constant do?" · **Summary:** "Replace the Reynolds stress by l_T²(dU/dy)² with l_T = κy: κ
  turns the log line, the wall damping slides it, and switching the model off returns the laminar parabola." · **CORE:** C12
  (also N132, N162, N167, N169–N172, N176, N181) · **Reference:** `stride_padding_playground.html` (classic presets, formula
  with numbers, badges) + `forced_damped_vibrations.html` (Explain).
- **meta:** `viz:order 8` · `viz:sections 12.10` · `viz:equations 12.94 12.98 12.100 12.101` · `viz:fluidpy
  ch12.mixing_length_wall_profile ch12.mixing_length_intercept ch12.channel_mixing_length ch12.mixing_length_eddy_viscosity
  ch12.mixing_length_stress` · `viz:derivations D21`.
- **Physics (JS ↔ Python):** `lplus(yp, kappa, Ap)` = κy⁺(1 − e^{−y⁺/A⁺}) (A⁺ = 0 ⇒ κy⁺); `slope(yp, kappa, Ap, tau)` = 2τ⁺/(1 +
  √(1 + 4l⁺²τ⁺)) ↔ `ch12.mixing_length_wall_profile(yplus, kappa, damping, A_plus)["slope"]` (τ⁺ = 1 in the constant-stress
  layer); `Uplus(yp, …)` by Simpson on a log grid (C.5-3) ↔ `["Uplus"]`; `intercept(kappa, Ap)` = lim [U⁺ − (1/κ) ln y⁺] — closed
  form [ln(4κ) − 1]/κ when A⁺ = 0, numerical otherwise ↔ `ch12.mixing_length_intercept(kappa, A_plus)`; `nuT(yp, …)` = l⁺²·slope ↔
  `["nuT_over_nu"]`; channel: τ⁺ = 1 − y⁺/Re_τ with l⁺ capped at 0.09 Re_τ ↔ `ch12.channel_mixing_length(Re_tau, kappa, A_plus)`
  (table from `reference/ch12/explainer_tables.json` for the four Re_τ; live integration for other values); laminar ghost U⁺ =
  y⁺(1 − y⁺/2Re_τ) ↔ `core.laminar.channel_flow` (same pressure gradient).
- **Views** (rows [1, 1]): 1. `chan` "The channel" (row 0, flex 1): U/U_centre-laminar against y/δ from the wall to the
  centreline: the model (purple), the laminar parabola at the same pressure gradient (muted, 20× faster — drawn scaled with
  the factor printed), DNS dots when available. 2. `semi` "U⁺(y⁺)" (row 0, flex 1.2; semi-log): the model bold; a tangent log
  line with its slope 1/κ and intercept B marked by a bracket on the y⁺ = 1 axis; the reference line κ = 0.41, B = 5.0 (muted
  dashed); the sublayer line. Pointer: click → inspector. 3. `nut` "Eddy viscosity and mixing length" (row 1, `hidePortrait`):
  ν_T/ν (orange) and l_T⁺ (teal) against y⁺ on log axes, κy⁺ dashed.
- **Controls:** `kappa` "$\kappa$" (0.30…0.50, default 0.41) · `Ap` "Damping $A^+$" (0…40, default 26; 0 = none) · `lret`
  "$\log_{10}\mathrm{Re}_\tau$" (2.25…4, default 3) · `model` toggle on / off (laminar) · `cap` "Core cap of $l_T/\delta$" (0.05…0.2,
  default 0.09; optional).
- **Presets:** "no damping — intercept too low" {Ap: 0} · "van Driest 26" {Ap: 26, kappa: 0.41} · "laminar" {model: off} · "Re_τ =
  5200" {lret: 3.716}.
- **Status:** "κ = 0.41, A⁺ = 26: B = 5.28 — matches smooth-wall data (≈ 5)" · "⚠️ no damping: B = −1.23 — right slope, line 6
  units too low" · "model off: laminar parabola, centreline U⁺ = Re_τ/2 = 500".
- **Readouts:** "Slope 1/κ" · "Implied B" · "ν_T/ν at cursor" · "U_bulk⁺" · "C_f" · "C_f laminar".
- **Depth features:** Explain + Code + Derivation · **linked views** · **presets** · **status** · **inspector** (click y⁺ = 10, A⁺
  = 0: "κ²y⁺²s² + s = 1: 16.81 s² + s − 1 = 0 ⇒ s = 2/(1 + √(1 + 67.24)) = **0.216**; viscous 21.6 %, Reynolds 78.4 % of τ₀; ν_T/ν =
  16.81 × 0.216 = 3.63").
- **Explain:**
  0. *What the views show* — "Left: the mean profile the model predicts (<b class=c-accent>purple</b>) against the laminar
     flow driven by the same pressure gradient (grey). Right: the same profile in wall units; the bracket shows the slope and
     the intercept of its straight part. Bottom (hidden on phones): the eddy viscosity the model uses."
  1. *The model* — "⟨u_iu_j⟩ = ⅔ē δ_ij − ν_T(∂U_i/∂x_j + ∂U_j/∂x_i) (12.94) with ν_T ~ l_T u_T (12.98), u_T = l_T dU/dy: −⟨uv⟩ =
     l_T²(dU/dy)², l_T = κy."
  2. *The equation it gives* — "0 = −(1/ρ)dP/dx + ∂/∂y(ν dU/dy + κ²y²(dU/dy)²) (12.100); near the wall the bracket is the
     constant τ₀/ρ; in wall units s + l⁺²s² = 1."
  3. *At your height* — "l⁺ = κy⁺(1 − e^{−y⁺/A⁺}) = …; s = 2/(1 + √(1 + 4l⁺²)) = **…**; ν_T/ν = l⁺²s = …".
  4. *Far from the wall* — "l⁺ ≫ 1: s → 1/(κy⁺), i.e. dU/dy ≅ √(τ₀/ρ)/(κy), U/u_* ≅ (1/κ) ln y + const. (12.101): slope 1/κ =
     **2.44**."
  5. *The intercept the model implies* — "B = lim [U⁺ − (1/κ) ln y⁺] = **…**; without damping exactly [ln(4κ) − 1]/κ = −1.23."
  6. *The whole channel* — "U_bulk⁺ = …, C_f = 2/U_bulk⁺² = …; laminar at the same pressure gradient: U_bulk⁺ = Re_τ/3 = …".
  7. *Reading the current setting* — no damping: "The eddies are allowed to mix right down to the wall, so the sublayer is
     far too thin and the whole log line sits about 6 units low. The slope is right: κ controls the slope only." · A⁺ = 26:
     "Switching the eddies off within about 26 wall units restores the sublayer and puts the line where the data are. The
     model has two constants and each does one job." · other κ: "The line pivots: 2.303/κ per decade." · model off: "No
     turbulent stress: the parabola of Ch. 8, twenty times faster for the same push — the price of turbulent friction."
- **Derivation tab:** **D21** (11 steps) `view: 'semi'`; step 3 `watch` "the orange ν_T curve: it grows with distance from the
  wall because bigger eddies fit there"; step 7 `set` {Ap: 0, yp: 10}, `live` "16.81 s² + s − 1 = 0"; step 8 `live` "s = 2/(1 +
  √(1 + 4κ²y⁺²)) = 0.216"; step 10 `live` "s → 1/(κy⁺) = 0.0244 at y⁺ = 100 (model 0.0241)"; step 11 `set` {Ap: 0} then {Ap: 26},
  `watch` "the bracket: B jumps from −1.23 to 5.28". Interpret: section 7.
- **Code:**
  ```python
  p = ch12.mixing_length_wall_profile({{yp}}, {{kappa}}, damping={{damp}},
                                      A_plus={{Ap}})
  s = p["slope"]                  # dU+/dy+ = 2/(1 + sqrt(1 + 4 l+^2)) = {{s}}
  nuT = p["nuT_over_nu"]          # l+^2 s = {{nut}}
  B = ch12.mixing_length_intercept({{kappa}}, A_plus={{ApOrNone}})   # {{B}}
  ch = ch12.channel_mixing_length({{Ret}}, {{kappa}}, A_plus={{Ap}})
  Cf = ch["Cf"]                   # {{Cf}}  (laminar: {{CfLam}})
  ```
- **Walkthrough (6 steps):** 1. "Guess the stress" — "The averaged equations lack one thing: ⟨uv⟩. What is the simplest
  guess?" · 2. "A length and a velocity" — "ν_T ~ l_T u_T, with u_T = l_T dU/dy and l_T = κy: −⟨uv⟩ = κ²y²(dU/dy)²." `eq:
  'ml'`, `highlight: ['view:nut']` · 3. "A quadratic at every height" — "Constant total stress: κ²y⁺²s² + s = 1. Click y⁺ =
  10." `set` {Ap: 0}, `inspect: true`, `derive: {id: 'D21', step: 8}` · 4. "The logarithm falls out" — "Far from the wall s →
  1/(κy⁺): a straight line of slope 1/κ. But look where it sits: B = −1.2." `readouts: ['B']` · 5. "Damp the eddies" — "Slide
  A⁺ to 26: the line rises to B = 5.3. Now slide κ: it only turns." `controls: ['Ap', 'kappa']` · 6. "Your turn" — "Switch
  the model off. Predict the ratio of laminar to turbulent centreline speed, then read it." `controls: ['model']`.
- **Equations:** `tvh` ref 'Eq. (12.94)' ⟨u_iu_j⟩ = ⅔ē δ_ij − ν_T(∂U_i/∂x_j + ∂U_j/∂x_i) · `scal` ref 'Eq. (12.98)' ν_T, κ_T, or κ_mT ~
  l_T u_T · `ml` ref 'p. 594 (no number)' −⟨uv⟩ = l_T²(dU/dy)², live · `wall` ref 'Eq. (12.100)' 0 = −(1/ρ)dP/dx + ∂/∂y(ν dU/dy +
  κ²y²(dU/dy)²) · `slope` ref 'D21 (ours)' dU⁺/dy⁺ = 2/(1 + √(1 + 4κ²y⁺²)), live · `log` ref 'Eq. (12.101)' dU/dy ≅ √(τ₀/ρ)/(κy), or U/u_*
  ≅ (1/κ) ln y + const.
- **Check yourself:** (1) "Which constant changes the slope of the log line, which its height?" — "κ the slope; A⁺ the height."
  · (2) "Without damping, what B does κ = 0.41 give, and why is it so low?" — "−1.23 = [ln(4κ) − 1]/κ: turbulent mixing at the
  wall leaves almost no sublayer." `set {Ap: 0}` · (3) "At y⁺ = 10 with no damping, what share of the stress is viscous?" —
  "21.6 %." · (4) "Is ν_T a property of the fluid?" — "No: it is κ²y²∣dU/dy∣ — zero at the wall, large in the flow, different
  in every flow."
- **Selftest parity rows:** `{name: 'slope y+=10', js: slope(10, 0.41, 0, 1), py: 'ch12.mixing_length_wall_profile(10.0, 0.41)
  ["slope"]', rtol: 1e-12}` (0.215965) · `{name: 'U+ y+=100 damped', js: Uplus(100, 0.41, 26), py: 'ch12.mixing_length_wall_profile
  (100.0, 0.41, damping="van_driest", A_plus=26.0)["Uplus"]', rtol: 1e-4}` (16.528) · `{name: 'B undamped', js: intercept(0.41, 0), py:
  'ch12.mixing_length_intercept(0.41)', rtol: 1e-9}` (−1.23245) · `{name: 'B van Driest', js: intercept(0.41, 26), py:
  'ch12.mixing_length_intercept(0.41, A_plus=26.0)', rtol: 2e-3}` (5.277) · invariant `{name: 'closed form', js: intercept(0.41, 0),
  expect: (Math.log(1.64) - 1)/0.41, rtol: 1e-9}`.
- **Fit plan:** 360×640: status · `semi` (60 %) over `chan` (40 %); `nut` hidden — "ν_T/ν at cursor" is a readout; the bracket
  labels "1/κ" and "B" stay on `semi`. 844×390: `chan` | `semi`. Desktop: `nut` under both.

### E9 · stratified_surface_layer
- **Title:** "When does stratification kill turbulence?" · **Summary:** "Drag the surface heat flux through zero: the
  Monin–Obukhov length flips sign, the wind profile bends the other way, and the height where buoyancy beats shear moves." ·
  **CORE:** C14, C15 (also R05, R06, N54, N159, N182–N189) · **Reference:** `forced_damped_vibrations.html` (Explain with a
  regime-dependent reading) with the ch01 `parcel_stability` badge (two conventions).
- **meta:** `viz:order 9` · `viz:sections 12.11` · `viz:equations 12.107 12.108 12.109 12.110 12.111` · `viz:fluidpy
  ch12.surface_layer_state ch12.monin_obukhov_from_fluxes ch12.flux_richardson_surface_layer ch12.gradient_richardson_thermal
  ch12.turbulence_regime ch12.surface_layer_wind ch12.gradient_richardson_surface_layer` · `viz:derivations D24 D25`.
- **Physics (JS ↔ Python):** `wT(H, rho, cp)` = H/(ρc_p); `LM(ustar, wT, T, kappa)` = −u_*³T/(κ g ⟨wT′⟩) (±∞ at zero flux) ↔
  `ch12.monin_obukhov_from_fluxes(tau, H, rho, cp, T, kappa=)`; `Rf(z, L)` = z/L ↔ `ch12.flux_richardson_surface_layer`; `phiM(zeta)`
  = 1 + 5ζ (ζ ≥ 0), (1 − 16ζ)^{−1/4} (ζ < 0) ↔ `ch12.dimensionless_shear`; `wind(z, ustar, z0, L, kappa)` (log-linear when stable,
  Businger–Dyer when unstable, C.5-4) ↔ `ch12.surface_layer_wind(z, u_star, z0, L_M, kappa=, unstable="businger_dyer")`; `Ri(z, L,
  PrT)` = Pr_T·Rf ↔ `ch12.gradient_richardson_surface_layer`; `regime(Rf)` ↔ `ch12.turbulence_regime` (exact strings); temperature
  profile: dθ/dz = −(⟨wT′⟩/(κu_*z)) φ_h with φ_h = Pr_T φ_m (model, stated) and T(z) = θ(z) + Γ_a z; verdict strings ↔
  `ch12.surface_layer_state(…)["verdict_kundu"]`, `["verdict_met"]`; everything in one call ↔ `ch12.surface_layer_state`.
- **Views** (rows [1.2, 1]): 1. `wind` "Wind profile" (row 0, flex 1.2; height on a log axis 0.1…200 m, wind horizontal): U(z)
  bold purple; the neutral logarithm as a muted ghost; the book's log-linear form dashed where it differs (unstable side,
  labelled "outside its range"); a blue horizontal line at z = ∣L_M∣; shading: below it "forced convection / shear-driven",
  above it "free convection" (unstable, pale orange) or "turbulence starves above z = L_M/4" (stable, pale blue, starting
  at the Rf = ¼ height ▲); the cursor height. Pointer: drag the cursor. 2. `temp` "T(z) and θ(z)" (row 0, flex 0.8): T (blue
  solid), θ (blue dashed), the adiabat through the surface value (muted); the legend line gives the local gradient in
  both conventions. 3. `bud` "Budget at the cursor" (row 1, `hidePortrait`): bars shear production (orange, +), buoyancy
  (blue, ±), dissipation as the remainder (rose, −), with Rf and Ri printed.
- **Controls:** `ustar` "Friction velocity $u_*$" (0.05…1 m/s, default 0.3) · `H` "Surface heat flux $H$" (−100…+400 W/m², step 5,
  default −30) · `z0` chips sea 2 × 10⁻⁴ m / grass 0.03 m / forest 1 m · `z` "Height $z$" (cursor, 0.5…200 m, default 10) · `conv`
  toggle Kundu ($\Gamma=dT/dz$) / meteorology ($\Gamma=-dT/dz$) · `PrT` "$\nu_T/\kappa_T$" (0.5…2, default 1; optional).
- **Presets:** "sunny afternoon" {ustar: 0.4, H: 250} · "neutral overcast" {H: 0} · "clear calm night" {ustar: 0.1, H: −20} ·
  "strong wind night" {ustar: 0.5, H: −40} · "sea surface" {z0: sea, ustar: 0.25, H: 20}.
- **Status** (three parts, both conventions; exact strings pinned): "🌙 stable · L_M = +83 m · Rf(10 m) = 0.12 shear-driven ·
  stable ⇔ dT/dz > Γa: … > −9.8 K/km · (met: Γ < Γa: … < 9.8 K/km)" — the existing `lapse_rate_stability` wording; its meteorological "Γa" is +9.8 K/km / "☀️ unstable · L_M = −17 m · Rf(10 m) = −0.60 convective ·
  …" / "neutral · L_M = ∞ · logarithmic profile".
- **Readouts:** "⟨wT′⟩" [K m/s] · "L_M" [m] · "Rf(z)" · "Ri(z)" · "U(z)" · "U neutral" · "Rf = ¼ at z".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **presets** · **status** (verdict → criterion in the
  chosen convention → the other convention's relation) · **terms** (production, buoyancy, dissipation) · **notes** ("What a
  model's surface scheme does with these numbers: from the lowest-level wind and the air–surface temperature difference it
  iterates for u_* and H using exactly this profile, then hands τ₀ = ρu_*² and H to the model").
- **Explain:**
  0. *What the views show* — "Left: the mean wind against height (log axis). Grey: what it would be with no heat flux. Blue
     line: the height ∣L_M∣. Right: thermometer temperature T and potential temperature θ. Bottom (hidden on phones): the
     turbulent energy budget at the cursor height."
  1. *From watts to the kinematic flux* — "⟨wT′⟩ = H/(ρc_p) = −30/(1.2 × 1005) = **−0.0249 K m/s** (downward)."
  2. *The Monin–Obukhov length* — "L_M ≡ −u_*³/(κ α g ⟨wT′⟩) (12.110) with α = 1/T = 1/300: −0.027/(0.4 × 0.00333 × 9.81 × (−0.0249))
     = **+83 m**." boxed; "sign: positive when the heat flux is downward (stable)."
  3. *The flux Richardson number* — "Rf = −gα⟨wT′⟩/(−⟨uw⟩ dU/dz) (12.107); in the log layer −⟨uw⟩ = u_*², dU/dz = u_*/κz, so Rf = z/L_M
     (12.111) = 10/83 = **0.12**; it reaches ¼ at z = L_M/4 = **21 m**."
  4. *The gradient Richardson number* — "Ri ≡ N²/(dU/dz)² (12.108) = (ν_T/κ_T) Rf (12.109) = 1.0 × 0.12 = …; N² = gα dθ/dz."
  5. *The wind* — "U = (u_*/κ)[ln(z/z₀) + 5z/L_M] = 0.75 × (5.81 + 0.60) = **4.81 m/s**; neutral 4.36 m/s."
  6. *The same gradient in two conventions* — "at your height dT/dz = … K/km. Kundu: stable ⇔ dT/dz > Γ_a = −9.8 K/km.
     Meteorology: Γ = −dT/dz = …; stable ⇔ Γ < Γ_d = +9.8 K/km. Same layer, same verdict: negate the number, flip the
     inequality. N² = gα(dT/dz − Γ_a) = gα(Γ_d − Γ)."
  7. *The budget* (or hint) — "production u_*³/(κz)·φ_m = …; buoyancy gα⟨wT′⟩ = …; ratio = −Rf."
  8. *Reading the current setting* — stable: "Heat flows down into the cold ground; lifting dense air costs the turbulence
     energy. Below z ≈ L_M/4 the wind shear still pays the bill; above it the turbulence cannot sustain itself and the air
     decouples from the surface — the still nights when smoke lies flat." · unstable: "The ground heats the air; buoyancy
     adds energy at every height while shear production fades as 1/z. Above ∣L_M∣ thermals do the mixing (free convection)
     and the wind profile is nearly uniform." · neutral: "No heat flux: L_M is infinite and the profile is the logarithm of
     C11 all the way up." · very small u_* at night: "L_M ∝ u_*³: halve the wind and L_M falls eightfold. The turbulent layer
     collapses to a few metres."
- **Derivation tab:** **D24** (8 steps) `view: 'bud'` (wide) / `'wind'` (phones); step 4 `live` "−⟨uw⟩ dU/dz = … ; gα⟨wT′⟩ = …";
  step 5 `live` "Rf = … "; step 8 `live` "Ri = Pr_T Rf = …". **D25** (7 steps) `view: 'wind'`; step 3 `live` "Rf = κ z gα(−⟨wT′⟩)/u_*³
  = z/L_M = …"; step 6 `set` {H: −30}, `live` "U = (u_*/κ)[ln(z/z₀) + 5z/L_M] = 4.81 m/s"; step 7 `watch` "drag H to +150: the
  bold curve crosses to the other side of the grey logarithm". Interpret: section 8.
- **Code:**
  ```python
  st = ch12.surface_layer_state({{ustar}}, {{H}}, 300.0, {{z0}}, {{z}},
                                1.2, 1005.0, 0.4, convention="{{conv}}")
  L = st["L_M"]               # (12.110): -u*^3 / (kappa alpha g <wT'>) = {{L}} m
  Rf = st["Rf"]               # (12.111): z / L_M = {{Rf}}
  Ri = st["Ri"]               # (12.109): Pr_T * Rf = {{Ri}}
  U = st["U"]                 # {{U}} m/s   (neutral {{Un}})
  print(st["regime"], "|", st["verdict_kundu"], "|", st["verdict_met"])
  ```
- **Walkthrough (7 steps):** 1. "A still night" — "At dusk the wind near the ground dies. What switched the turbulence off?"
  `set` {H: -30, ustar: 0.3} · 2. "Cost and income" — "Buoyancy removes gα⟨wT′⟩; shear supplies −⟨uw⟩dU/dz. Their ratio is
  Rf." `terms: true`, `eq: 'rf'` · 3. "It grows with height" — "Shear production falls as 1/z, the heat flux does not: Rf =
  z/L_M." `derive: {id: 'D25', step: 3}`, `controls: ['z']` · 4. "A length" — "L_M = −u_*³/(κ α g ⟨wT′⟩) = +83 m: the blue line."
  `highlight: ['readout:L']` · 5. "Through zero" — "Drag H from −30 to +150 W/m²: L_M flips sign and the profile crosses the
  neutral line." `controls: ['H']` · 6. "Two conventions, one verdict" — "Toggle the convention: the number changes sign, the
  inequality flips, the verdict stays." `controls: ['conv']` · 7. "Your turn" — "Predict L_M when the night wind halves; then
  set u_* = 0.15." `notes: true`.
- **Equations:** `rf` ref 'Eq. (12.107)' Rf = −gα⟨wT′⟩/(−⟨uw⟩(dU/dz)), live · `ri` ref 'Eq. (12.108)' Ri ≡ N²/(dU/dz)² = αg(dT̄/dz)/(dU/dz)²,
  note "T̄ is potential temperature; with the thermometer gradient N² = gα(dT/dz − Γ_a) = gα(Γ_d − Γ_met)" · `link` ref 'Eq. (12.109)'
  Ri = (ν_T/κ_T) Rf, live · `lm` ref 'Eq. (12.110)' L_M ≡ −u_*³/(κ α g ⟨wT′⟩), live · `zl` ref 'Eq. (12.111)' Rf = z/L_M, live · `prof` ref
  'p. 598 (no number)' U = (u_*/κ)[ln(z/z₀) + 5z/L_M], live.
- **Check yourself:** (1) "Halve u_* at fixed H. What happens to L_M?" — "÷ 8: L_M ∝ u_*³." · (2) "Is an atmosphere cooling at
  6.5 K/km with height stable? Say it in both conventions." — "Yes: dT/dz = −6.5 > Γ_a = −9.8 K/km (Kundu); Γ = 6.5 < Γ_d = 9.8
  K/km (meteorology)." · (3) "At what height does Rf reach ¼ for L_M = 83 m?" — "21 m." · (4) "With an upward heat flux, is Rf
  positive or negative, and what does that mean?" — "Negative: buoyancy *adds* turbulent energy (convective)." `set {H: 150}`.
- **Selftest parity rows:** `{name: 'L_M night', js: LM(0.3, -30/(1.2*1005), 300, 0.4), py: 'ch12.monin_obukhov_from_fluxes(0.108,
  -30.0, 1.2, 1005.0, 300.0, kappa=0.4)', rtol: 1e-10}` (≈ 83.0 m; the JS uses the same g = 9.80665) · `{name: 'U 10 m stable',
  js: wind(10, 0.3, 0.03, 82.98165, 0.4), py: 'ch12.surface_layer_wind(10.0, 0.3, 0.03, 82.98165, kappa=0.4)', rtol: 1e-10}` (4.8088) ·
  `{name: 'U 10 m unstable', js: wind(10, 0.3, 0.03, -16.6, 0.4), py: 'ch12.surface_layer_wind(10.0, 0.3, 0.03, -16.6, kappa=0.4,
  unstable="businger_dyer")', rtol: 1e-9}` · `{name: 'Rf', js: Rf(10, 82.98165), py: 'ch12.flux_richardson_surface_layer(10.0, 82.98165)',
  rtol: 1e-12}` · `{name: 'regime', js: regime(0.12), py: 'ch12.turbulence_regime(0.12)', exact: true}` · `{name: 'verdict kundu', js:
  verdictKundu(0.010), py: 'ch12.gradient_richardson_thermal(0.010, 0.1, 1/300, Gamma_a=ch12.adiabatic_lapse_rate())["verdict_kundu"]', exact: true}` (`Gamma_a` is required — review M1; never type −0.0098).
- **Fit plan:** 360×640: status on two short lines (verdict · L_M · Rf / criterion in the chosen convention; the other
  convention's relation moves to the Explain tab's section 6 and the `temp` title) · `wind` (60 %) over `temp` (40 %); `bud`
  hidden — Rf and Ri are readouts. 844×390: `wind` | `temp`. Desktop: `bud` under both.

### E10 · taylor_dispersion
- **Title:** "Why does a puff spread like t, then like √t?" · **Summary:** "Particles leave a source on one clock: while they
  remember their velocity the cloud is a wedge, afterwards a random walk — and the eddy diffusivity grows before it
  settles." · **CORE:** C16 (also N06, N196, N199, N208, N211–N216, slips #12, #13) · **Reference:**
  `angular_frequency_explorer_1.html` (one clock, modes, end-of-run summary, regime notes).
- **meta:** `viz:order 10` · `viz:sections 12.12` · `viz:equations 12.117 12.119 12.121 12.123 12.127 12.129` · `viz:fluidpy
  ch12.taylor_dispersion_exponential ch12.taylor_dispersion_gaussian ch12.eddy_diffusivity_exponential ch12.langevin_particles
  ch12.dispersion_regime ch12.dispersion_local_slope ch12.smoke_plume_width` · `viz:derivations D26 D27 D28`.
- **Physics (JS ↔ Python):** `X2exp(t, u2, Lam)` = 2u²Λ²[t/Λ − 1 + e^{−t/Λ}] (with `expm1` for small t/Λ) ↔
  `ch12.taylor_dispersion_exponential`; `X2gauss(t, u2, tc)` = 2u²[(√π/2)t_c t erf(t/t_c) − (t_c²/2)(1 − e^{−t²/t_c²})] ↔
  `ch12.taylor_dispersion_gaussian`; `DT(t, u2, Lam)` = u²Λ(1 − e^{−t/Λ}) ↔ `ch12.eddy_diffusivity_exponential`; `slope(t, Lam)` =
  d ln⟨X²⟩/d ln t ↔ `ch12.dispersion_local_slope`; `regime(t, Lam)` ↔ `ch12.dispersion_regime` (exact strings); `plumeWidth(x, U,
  w, Lam)` = √X2exp(x/U, w², Λ) ↔ `ch12.smoke_plume_width`; particles: exact OU velocity update + exact position update
  (`Viz.rng`) ↔ `ch12.langevin_particles` (same law; no parity on samples).
- **Views** (rows [1.1, 1]): 1. `cloud` "The cloud" (row 0, full width): **puff mode** — 300 particles leaving a point, time
  running (trails faint), the ±X_rms envelope bold teal, the ballistic wedge ±u_rms t and the diffusive parabola
  ±u_rms√(2Λ_t t) dashed, a muted Gaussian ghost of constant D = u²Λ_t; **plume mode** — the time-averaged plume behind a
  chimney in a wind U (colour = concentration, amber), the same two envelopes drawn against x = Ut, labels "∝ x" near and
  "∝ √x" far. 2. `msd` "⟨X²⟩ against time" (row 1, flex 1.2; log–log): Taylor's formula bold, particle dots, the two limits
  dashed (slopes 2 and 1 labelled), the diffusive line with its offset (thin), the constant-D ghost (muted), a vertical
  line at t = Λ_t, the current time ●. Pointer: click → set t. 3. `corr` "r(τ) and D_T(t)" (row 1, flex 0.8, `hidePortrait`):
  r(τ) with the area up to t shaded; D_T(t)/u²Λ_t rising to 1.
- **Controls:** `urms` "$u_{rms}$" (0.1…3 m/s, default 1) · `Lam` "Memory $\Lambda_t$" (1…100 s log, default 10) · `sys` chips puff
  in time / plume in a wind (mode) + `U` "Wind $U$" (1…15 m/s, default 5; plume only) · `shape` chips exponential / Gaussian
  (optional) · transport `t`.
- **Transport:** `t` 0…20 Λ_t (log-paced so both regimes get screen time), `end: 'hold'`; end card: "measured slopes: 1.97 for
  t < 0.3 Λ_t, 1.03 for t > 10 Λ_t · X_rms(20 Λ_t) = 61.6 m · D_T reached 100 % of u²Λ_t".
- **Presets:** "short memory" {Lam: 1} · "long memory" {Lam: 100} · "t = Λ_t" {t: Lam} · "chimney in a 5 m/s wind" {sys: plume, U:
  5, urms: 0.5, Lam: 20}.
- **Status** (exact `dispersion_regime` strings + numbers): "t = 0.3 Λ_t: ballistic — X_rms ≈ u_rms t; D_T still growing (26 %
  of final)" · "t = Λ_t: transition — neither limit holds (slope 1.72)" · "t = 10 Λ_t: diffusive — X_rms ≈ u_rms√(2Λ_t t); D_T =
  u²Λ_t".
- **Readouts:** "X_rms" [m] · "ballistic u t" · "diffusive √(2u²Λt)" · "local slope" · "D_T" [m²/s] · "D_T / final".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** (end card with measured slopes) ·
  **modes** (puff / plume) · **presets** · **status** · **terms** (D_T reached | still to come, a two-part bar summing to u²Λ_t).
- **Explain:**
  0. *What the views show* — "Top: marked particles released at one point; the <b class=c-teal>teal</b> envelope is their
     rms distance. Dashed: the two limiting laws. Grey: what a constant diffusivity would give. Bottom left: the
     mean-square distance on log axes — slope 2 means 'like t', slope 1 means 'like √t'. Bottom right (hidden on phones): the
     velocity correlation and the eddy diffusivity."
  1. *Taylor's formula with your numbers* — "⟨X_α²⟩ = 2⟨u_α²⟩ t ∫₀ᵗ(1 − τ/t) r_α(τ) dτ (12.119); for r = e^{−τ/Λ_t}: 2u²Λ_t²[t/Λ_t − 1 +
     e^{−t/Λ_t}] = 2 × 1 × 100 × [1 − 1 + 0.3679] = **73.6 m²** at t = 10 s; X_rms = **8.58 m**."
  2. *Against the short-time law* — "(X_α)_rms = (u_α)_rms t (12.121) = 10.0 m — 17 % too large at t = Λ_t; within 2 % for t <
     0.1 Λ_t."
  3. *Against the long-time law* — "(X_α)_rms = (u_α)_rms√(2Λ_t t) (12.123) = 14.1 m — 65 % too large at t = Λ_t; within 3 % for
     t > 20 Λ_t. The exact curve approaches 2u²Λ_t t − 2u²Λ_t²: a straight line that does not pass through the origin."
  4. *The eddy diffusivity* — "D_T ≡ ½ d⟨X_α²⟩/dt = ⟨u_α²⟩∫₀ᵗ r_α dτ (12.127) = u²Λ_t(1 − e^{−t/Λ_t}) = 10 × 0.632 = **6.3 m²/s**;
     final value D_T ≅ ⟨u_α²⟩Λ_t (12.129, for t ≫ Λ_t — the page prints ≪) = 10 m²/s."
  5. *As a random walk* — "(R_n)_rms = L√n (12.125) with steps of duration 2Λ_t = 20 s and length u_rms × 2Λ_t = 20 m: after t =
     1000 s, n = 50, 20√50 = 141 m = u_rms√(2Λ_t t) ✓."
  6. *Plume mode* (or hint) — "t = x/U: at x = 100 m in a 5 m/s wind, t = 20 s = Λ_t: Z_rms = 8.6 m. Near the chimney the
     width grows ∝ x, far away ∝ √x (the book's figure caption has the two swapped — slip #13)."
  7. *Reading the current setting* — ballistic: "Each particle still has the velocity it left with, so the cloud opens as a
     wedge. Calling this 'diffusion' would be wrong: a constant-D cloud (grey) starts far too wide." · transition: "Particles
     are forgetting. No simple law: use the full formula." · diffusive: "Each particle has changed its mind many times; the
     cloud grows like a random walk and one number, D_T = u²Λ_t, describes it — an eddy diffusivity is (rms velocity)² ×
     memory, the same product of a velocity and a length as ν_T ~ l_T u_T (12.98)." · long memory: "Slow, large eddies: the wedge lasts longer and the final diffusivity is
     larger. In the ocean Λ_t is days, which is why drifter clusters spread ballistically for days."
- **Derivation tab:** **D26** (11 steps) `view: 'corr'` (wide) / `'msd'` (phones); step 6 `live` "∫₀ᵗ r dτ = Λ_t(1 − e^{−t/Λ_t}) = …"
  with the shaded area; step 11 `live` "⟨X²⟩ = … m²". **D27** (8 steps) `view: 'msd'`; step 2 `set` {t: 0.1 Lam}, `live` "u²t² = …
  vs exact …"; step 5 `set` {t: 20 Lam}, `live` "2u²Λ_t t = … vs exact …"; step 6 `watch` "the thin line with the offset lies on
  the bold curve; the dashed one stays above it". **D28** (7 steps) `view: 'corr'` / `'msd'`; step 4 `live` "D_T = u²∫₀ᵗ r dτ = …";
  step 7 `watch` "the two-part bar: D_T reaches its final value only after several Λ_t". Interpret: section 7.
- **Code:**
  ```python
  X2 = ch12.taylor_dispersion_exponential({{t}}, {{u2}}, {{Lam}})   # (12.119): {{X2}} m2
  ball = {{u2}} * {{t}}**2                    # (12.120): {{ball}}
  diff = 2 * {{u2}} * {{Lam}} * {{t}}         # (12.122): {{diff}}
  DT = ch12.eddy_diffusivity_exponential({{t}}, {{u2}}, {{Lam}})    # (12.127): {{DT}}
  print(ch12.dispersion_regime({{t}}, {{Lam}}),
        ch12.dispersion_local_slope({{t}}, {{Lam}}))               # {{reg}}, {{slope}}
  X, u = ch12.langevin_particles(2000, t_grid, {{urms}}, {{Lam}}, seed=0)
  ```
- **Walkthrough (7 steps):** 1. "A wedge, then a parabola" — "Smoke opens like a wedge near a chimney and ever more slowly
  far away. Press ▶." `play: true` · 2. "Remembering" — "At first each particle keeps its velocity: X = ut, so X_rms = u_rms
  t. Slope 2." `set` {t: 1}, `eq: 'short'` · 3. "Forgetting" — "After a few memory times the steps are independent: a
  random walk, X_rms = u_rms√(2Λ_t t). Slope 1." `set` {t: 200}, `eq: 'long'` · 4. "One formula for both" — "⟨X²⟩ = 2⟨u²⟩ t ∫₀ᵗ
  (1 − τ/t) r dτ: only the velocity variance and its correlation enter." `derive: {id: 'D26', step: 11}` · 5. "Not a
  constant" — "D_T is the running area under r: it grows as u²t, then settles at u²Λ_t." `terms: true`, `eq: 'dt'` · 6. "Move
  the bend" — "Drag Λ_t: the change of slope follows it." `controls: ['Lam']` · 7. "Your turn" — "Switch to the plume.
  Predict the distance where the wedge ends in a 5 m/s wind with Λ_t = 20 s." `set` {sys: 'plume'}.
- **Equations:** `rate` ref 'Eq. (12.117)' d⟨X_α²⟩/dt = 2⟨u_α²⟩∫₀ᵗ r_α(τ)dτ · `tay` ref 'Eq. (12.119)' ⟨X_α²⟩(t) = 2⟨u_α²⟩ t ∫₀ᵗ(1 − τ/t)r_α(τ)dτ,
  live · `short` ref 'Eq. (12.121)' (X_α)_rms = (u_α)_rms t for t ≪ Λ_t, live · `long` ref 'Eq. (12.123)' (X_α)_rms = (u_α)_rms√(2Λ_t t) for
  t ≫ Λ_t, live · `dt` ref 'Eq. (12.127)' D_T ≡ ½ d⟨X_α²⟩/dt = ⟨u_α²⟩∫₀ᵗ r_α(τ)dτ, live · `dtl` ref 'Eq. (12.129), condition corrected' D_T ≅
  ⟨u_α²⟩Λ_t for t ≫ Λ_t, note "the page prints t ≪ Λ_t (slip #12); the short-time form is D_T ≅ ⟨u_α²⟩t (12.128)".
- **Check yourself:** (1) "At t = Λ_t, which limit is closer?" — "The ballistic one (17 % high against 65 %): the local slope
  is still 1.72." `set {t: 10, Lam: 10}` · (2) "Multiply Λ_t by 10 at fixed u_rms: what happens to the late-time X_rms and to
  D_T?" — "X_rms × √10; D_T × 10." · (3) "Why does a constant-D model fail near the source?" — "It spreads as √t from the
  start — infinitely fast at t = 0 — while real particles need a memory time before they random-walk." · (4) "Is D_T ≅ ⟨u²⟩Λ_t
  valid for t ≪ Λ_t, as printed?" — "No: there D_T ≅ ⟨u²⟩t; the constant value needs t ≫ Λ_t."
- **Selftest parity rows:** `{name: 'X2 t=10', js: X2exp(10, 1, 10), py: 'ch12.taylor_dispersion_exponential(10.0, 1.0, 10.0)', rtol:
  1e-12}` (73.5759) · `{name: 'X2 small t', js: X2exp(0.01, 1, 10), py: 'ch12.taylor_dispersion_exponential(0.01, 1.0, 10.0)', rtol:
  1e-9}` · `{name: 'X2 gaussian', js: X2gauss(10, 1, 10), py: 'ch12.taylor_dispersion_gaussian(10.0, 1.0, 10.0)', rtol: 1e-9}` (86.1528) ·
  `{name: 'D_T', js: DT(10, 1, 10), py: 'ch12.eddy_diffusivity_exponential(10.0, 1.0, 10.0)', rtol: 1e-12}` (6.32121) · `{name: 'slope',
  js: slope(10, 10), py: 'ch12.dispersion_local_slope(10.0, 10.0)', rtol: 1e-10}` (1.71828) · `{name: 'plume width', js: plumeWidth
  (1000, 5, 0.5, 20), py: 'ch12.smoke_plume_width(1000.0, 5.0, 0.5, 20.0)', rtol: 1e-10}` (42.4265) · `{name: 'regime', js: regime(3,
  10), py: 'ch12.dispersion_regime(3.0, 10.0)', exact: true}`.
- **Fit plan:** 360×640: status (1 line) · `cloud` (45 %) over `msd` (55 %); `corr` hidden — "D_T / final" is a readout and the
  terms bar sits in the Explore tab; transport without step buttons. 844×390: `cloud` | `msd`. Desktop: rows [1.1, 1].

### B1 · k_epsilon_decay
(backup — built only if one of E1–E10 fails review)
- **Title:** "Where do the k–ε constants come from?" · **Summary:** "Let model turbulence decay: the slope of ē(t) fixes C_ε2,
  and the log layer ties the rest to the von Kármán constant." · **CORE:** C13 (also N175–N178, N180) · **Reference:**
  `forced_damped_vibrations.html`.
- **meta:** `viz:order 11` · `viz:sections 12.10` · `viz:equations 12.103 12.104 12.105` · `viz:fluidpy ch12.k_epsilon_decay
  ch12.k_epsilon_loglayer_kappa ch12.k_epsilon_eddy_viscosity ch12.k_epsilon_length_scale` · `viz:derivations D22 D23`.
- **Physics (JS ↔ Python):** `decay(e0, eps0, t, C2)` → [ē, ε̄, n, t₀] with n = 1/(C_ε2 − 1), t₀ = nē₀/ε̄₀, ē = ē₀(1 + t/t₀)^{−n}, ε̄ =
  ε̄₀(1 + t/t₀)^{−(n+1)} ↔ `ch12.k_epsilon_decay`; `kappaLog(Cmu, C1, C2, sigEps)` = [√C_μ(C_ε2 − C_ε1)σ_ε]^{1/2} ↔
  `ch12.k_epsilon_loglayer_kappa`; `nuT(e, eps, Cmu)` ↔ `ch12.k_epsilon_eddy_viscosity`; `lT(e, eps)` ↔ `ch12.k_epsilon_length_scale`;
  an RK4 integration of the pair as a visual cross-check (dots).
- **Views:** 1. `decay` "ē(t) and ε̄(t)" (log–log): closed form (teal, rose), RK4 dots, the slopes −n and −(n + 1) as
  labelled triangles, a band of measured decay exponents 1.1–1.3 (muted, "grid-turbulence experiments, public range"). 2.
  `scales` "ν_T and l_T" (`hidePortrait`): both against time. 3. `dial` "Implied κ": a gauge from 0.3 to 0.5 with the needle
  at the model's κ and a mark at 0.41.
- **Controls:** `C2` "$C_{\varepsilon2}$" (1.5…2.5, default 1.92) · `C1` "$C_{\varepsilon1}$" (1.0…1.8, default 1.44) · `Cmu` "$C_\mu$" (0.05…0.15,
  default 0.09) · `sigEps` "$\sigma_\varepsilon$" (0.8…1.6, default 1.3) · `e0`, `eps0` (optional) · transport `t` (log-paced).
- **Presets:** "standard constants" · "C_ε2 = 2.0 (n = 1)" · "measured decay n = 1.3" {C2: 1.769}. **Status:** "n = 1.09 — inside
  the measured range 1.1–1.3? just below · implied κ = 0.433 (measured ≈ 0.38–0.41)".
- **Depth features:** Explain + Code + Derivation · linked views · transport · presets · status · terms (in the budget of ē:
  production 0, dissipation −ε̄, transport 0 for decaying turbulence; in the log-layer mode: production = dissipation).
- **Explain:** 0. what the views show · 1. "decay: dē/dt = −ε̄, dε̄/dt = −C_ε2 ε̄²/ē — the homogeneous form of ∂ē/∂t + U_j∂ē/∂x_j =
  ∂/∂x_j((ν_T/σ_e)∂ē/∂x_j) − ε̄ − ⟨u_iu_j⟩∂U_i/∂x_j (12.103) and of ∂ε̄/∂t + U_j∂ε̄/∂x_j = ∂/∂x_j((ν_T/σ_ε)∂ε̄/∂x_j) − C_ε1(⟨u_iu_j⟩∂U_i/∂x_j)ε̄/ē − C_ε2ε̄²/ē (12.105)" · 2. "n = 1/(C_ε2 − 1) = 1/0.92 = **1.087**; t₀ = nē₀/ε̄₀" ·
  3. "ν_T = C_μē²/ε̄ (12.104) = …; l_T = ē^{3/2}/ε̄ = … (grows: small eddies die first)" · 4. "log layer: κ² = √C_μ(C_ε2 − C_ε1)σ_ε = 0.3 ×
  0.48 × 1.3 = 0.187 ⇒ κ = **0.433**; ē = u_*²/√C_μ = 3.33 u_*²" · 5. at the current time (live) · 6. *Reading the current
  setting*: standard — "The five constants are a compromise: the decay exponent is at the low edge of the measurements and
  κ a little high. Changing one to fix a flow spoils another — which is why the standard set has survived." · C_ε2 = 2 —
  "Energy ∝ 1/t exactly: the 'final period' law, too slow for grid turbulence." · n = 1.3 — "Matches the upper measurements
  but, with the other constants fixed, the implied κ falls to 0.36."
- **Derivation tab:** **D22** (6 steps) `view: 'scales'`; **D23** (12 steps) `view: 'decay'`; step 4 `live` "ε̄ = ε̄₀(ē/ē₀)^{C_ε2}"; step
  7 `live` "n = …"; step 12 `live` "κ = …" with the dial highlighted.
- **Code:**
  ```python
  e, eps, n, t0 = ch12.k_epsilon_decay({{e0}}, {{eps0}}, {{t}}, C_eps2={{C2}})  # n = {{n}}
  nuT = ch12.k_epsilon_eddy_viscosity(e, eps, C_mu={{Cmu}})                    # (12.104): {{nuT}}
  kappa = ch12.k_epsilon_loglayer_kappa({{Cmu}}, {{C1}}, {{C2}}, {{sigEps}})   # {{kappa}}
  ```
- **Walkthrough (5 steps):** 1. "Five constants — fudge?" — "k–ε carries five numbers. Are they free knobs, or does each answer to a measurement? Press ▶ and let model turbulence decay." · 2. "Let it decay" `play: true` — "No shear, no production: ē falls
  as a power law." · 3. "The slope is C_ε2" — "n = 1/(C_ε2 − 1). Drag C_ε2 until n matches 1.3." `controls: ['C2']`, `derive:
  {id: 'D23', step: 7}` · 4. "The log layer" — "Production = dissipation gives κ² = √C_μ(C_ε2 − C_ε1)σ_ε: watch the dial."
  `highlight: ['view:dial']` · 5. "Your turn" — "Can you get n = 1.3 and κ = 0.41 together by changing σ_ε?"
- **Equations:** `k` ref 'Eq. (12.103)' ∂ē/∂t + U_j∂ē/∂x_j = ∂/∂x_j((ν_T/σ_e)∂ē/∂x_j) − ε̄ − ⟨u_iu_j⟩∂U_i/∂x_j · `nut` ref 'Eq. (12.104)' ν_T =
  C_μē²/ε̄, live · `eps` ref 'Eq. (12.105)' ∂ε̄/∂t + U_j∂ε̄/∂x_j = ∂/∂x_j((ν_T/σ_ε)∂ε̄/∂x_j) − C_ε1(⟨u_iu_j⟩∂U_i/∂x_j)ε̄/ē − C_ε2 ε̄²/ē (aligned).
- **Check yourself:** (1) "Which constant alone sets the decay exponent?" — "C_ε2." · (2) "With the standard set, what κ does
  the model imply?" — "0.433." · (3) "Does the length scale ē^{3/2}/ε̄ grow or shrink during decay?" — "It grows as
  (1 + t/t₀)^{1 − n/2}: the small eddies die first."
- **Selftest parity rows:** `{name: 'e(10)', js: decay(1, 1, 10, 1.92)[0], py: 'ch12.k_epsilon_decay(1.0, 1.0, 10.0, C_eps2=1.92)[0]',
  rtol: 1e-10}` (0.080112) · `{name: 'kappa', js: kappaLog(0.09, 1.44, 1.92, 1.3), py: 'ch12.k_epsilon_loglayer_kappa(0.09, 1.44, 1.92,
  1.3)', rtol: 1e-12}` (0.432666).
- **Fit plan:** 360×640: `decay` (65 %) over `dial` (35 %); `scales` hidden.

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

| Section | What costs time | Full | FAST |
|---|---|---|---|
| front matter, §12.1–12.2 | imports; 128² synthetic field + white noise | 4 s | 3 s |
| §12.3 C01 | 64-member ensemble × 2000 samples; sympy window integral; A1 (40 frames, dpi 80); F1 (30 steps × 4 traces × 400 pts) | 16 s | 9 s |
| §12.4 C02, C03 | one 2 × 10⁵-sample record; FFT correlations; A2 (60 frames); three-routes figure; F2 (25 steps) | 18 s | 10 s |
| §12.5 C04 | `rans_sympy`, `reynolds_stress_budget_sympy` (cached to `outputs/ch12/cache/`, P252: ≈ 15 s first run, < 1 s after); 10⁵ parcels; A3 (80 frames) | 22 s (first) / 9 s | 6 s |
| §12.6 C05 | 3-D synthetic field 64³ + f, g by FFT (cached); D07, D08 sympy checks (< 5 s each) | 14 s | 8 s (48³) |
| §12.7 C06–C08 | `tke_budget_sympy`, `mean_energy_budget_sympy` (cached); channel budget at one Re_τ; A4 (11 frames); F3 (21 steps × 3 traces × 300 pts); model spectrum quadratures | 26 s | 14 s |
| §12.8 C09 | `plane_jet_similarity_sympy` (cached, ≈ 8 s first run); F4 (24 steps); table | 14 s | 8 s |
| §12.9 C10, C11 | Spalding inversions (vectorised brentq, 400 pts); F5 (20 steps × 3 traces); `overlap_matching_sympy` | 12 s | 8 s |
| §12.10 C12, C13 | mixing-length profiles for F6 (21 × 3, cached); channel at Re_τ = 1000 (Picard, < 0.5 s); `solve_ivp` decay | 12 s | 7 s |
| §12.11 C14, C15 | closed forms; F7 (31 steps × 3 traces); `temperature_variance_sympy` (cached) | 9 s | 6 s |
| §12.12 C16 | 10⁴ Langevin particles × 10³ steps (FAST 2000 × 400); 4000 random walkers; A5 (90 frames, FAST 40); F8 (25 steps) | 20 s | 10 s |
| explainers | 10 `show_viz` cells (iframes; no computation) | 3 s | 3 s |
| **Total** | | **≈ 170 s (first run ≈ 195 s)** | **≈ 95 s** |

FAST plan: `FAST = setup_notebook()`; animations ≤ 90 frames at dpi 80 (FAST 24–40), `player="frames"` for A1 and A4 and
`"video"` for A2, A3, A5; slider figures ≤ 31 steps × ≤ 4 traces × ≤ 400 points; the nine sympy engines and the 3-D field
are cached with a parameter hash; `k_epsilon_channel` is not run. Page budget: five animations < 2 MB each, eight plotly
figures ≈ 150 kB each.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is explained.
"primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text is the primer term, used verbatim in
`nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the earlier primer; a CORE/RECAP
id alone = taught there; "Cxx (Nnn)" = the NOTE placed in that block; "Cxx (D0n)" = the derivation where it is built; "Cxx (gloss
…)" = one sentence where it is used. Earlier chapters' material that is neither a primer nor a ch12 RECAP is named by chapter
and section. New primers P280–P306 (27) in first-use order: probability vocabulary, standard error, ergodicity,
Ornstein–Uhlenbeck signal, sliding window, sinc (C01; random-phase fields in §12.1, which opens C01), covariance matrix,
never-negative quadratic, correlation by FFT (C02), Fourier pair, DFT as a Riemann sum, Welch, density change of variable
(C03), sympy averaging operator, symmetric-tensor counting, kinematic vs dynamic fluxes (C04), derivatives of functions of r
(C05), log-grid quadrature (C08), semi-log axes, inner/outer/overlap matching (C10), composite profile (C11), Picard on the
eddy viscosity (C12), dividing two ODEs (C13), stability parameter and flux–profile integration (C15), triangle double
integral, induction (C16).

| Concept | First used in | Explained by |
|---|---|---|
| the conventions block: overloaded symbols κ, k / K, e / E / ε, λ / Λ, η, f / g / F / G, R_ij, τ, δ, θ, α, S, N, L, σ, Π, U₀, h / d | C01 | front matter (A.0 row 5, before C01), repeated per block |
| tilde = total, capital or over-bar = mean, lower case or prime = fluctuation | C01 | front matter (A.0 row 5) and C04 (N42) |
| five marks of turbulence; random waves are not turbulence | C01 | C01 (N01) |
| turbulent field is divergence-free, white noise is not | C01 | C01 (N02) |
| random-phase synthetic fields (np.fft.ifftn) | C01 | primer (in C01) |
| largest eddy ~ layer thickness; turbulent / irrotational interface | C01 | C01 (N03) |
| scope: incompressible, no Coriolis; 2-D geostrophic turbulence deferred | C01 | C01 (N04) |
| historical timeline (Reynolds, Taylor, Prandtl, von Kármán, Richardson, Kolmogorov, Obukhov) | C01 | C01 (N05) |
| realization, ensemble, ensemble average, expected value | C01 | C01 (N07) |
| random variable, probability density and histogram | C01 | primer (in C01) |
| Ornstein–Uhlenbeck signal (a Langevin equation) | C01 | primer (in C01) |
| m-th moment (12.1); mean (12.10) | C01 | C01 (N19) |
| stationary; time average (12.2) | C01 | C01 (N08) |
| sliding-window average with np.cumsum | C01 | primer (in C01) |
| homogeneous; volume average (12.3) | C01 | C01 (N09) |
| ergodicity | C01 | primer (in C01) |
| which average in the field; t_c ≪ Δt ≪ drift time | C01 | C01 (N11) |
| commutation rules (12.4)–(12.9) | C01 | C01 (D01; N12–N17) |
| mean of a product = product of means + covariance | C01 | C01 (D01; N18) |
| standard error of a mean | C01 | primer (in C01) |
| scatter of an N-member mean ∝ N^{−1/2} | C01 | C01 (N20) |
| central moments (12.11); variance, skewness, kurtosis, rms (book vs usual normalisation) | C01 | C01 (N21, N22) |
| the sinc function and np.sinc | C01 | primer (in C01) |
| window factors sinh(x)/x and sin(x)/x (Example 12.1) | C01 | C01 (N23) |
| sinh | C01 | knowledge/primers.md: hyperbolic functions (ch07 P168) — reminder |
| seeded random generators (np.random.default_rng) | C01 | knowledge/primers.md: seeded generator (ch01 P10) — reminder |
| Gaussian distribution | C01 | knowledge/primers.md: Gaussian components (ch01 P11) — reminder |
| log–log slopes and np.polyfit | C01 | knowledge/primers.md: power laws on log axes (ch01 P13) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| animate / show_animation, slider_figure, show_viz | C01 | knowledge/primers.md: animate (ch01 P16), slider_figure (P17), show_viz (P18) — reminder |
| partial derivative ∂/∂t, ∂/∂x_j | C01 | knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| integral as a limit of sums; fundamental theorem of calculus | C01 | knowledge/primers.md: fundamental theorem of calculus (ch02 P84) — reminder |
| sympy (symbols, integrate, simplify) | C01 | knowledge/primers.md: sympy (ch01 P40; series P117) — reminder |
| Lorenz sensitivity to initial conditions (why statistics) | C01 | Ch. 11 §11.14 (reminder in N07) |
| correlation tensor R_ij (12.12); autocorrelation (12.13) | C02 | C02 (N24, N26) |
| covariance, the covariance matrix and its ellipse | C02 | primer (in C02) |
| uncorrelated / correlated / anticorrelated; cross- vs autocorrelation | C02 | C02 (N25) |
| correlation coefficients r_12 (12.14), r_11 (12.15) | C02 | C02 (N27, N28) |
| a quadratic that is never negative has discriminant ≤ 0 | C02 | primer (in C02) |
| Schwartz inequality (12.16), −1 ≤ r ≤ 1 | C02 | C02 (D02; N29) |
| quadratic formula | C02 | knowledge/primers.md: complex square roots and the quadratic formula (ch06 P159) — reminder |
| lag form R_11(τ) (12.17), evenness, R_ij(−τ) = R_ji(τ) | C02 | C02 (D03) |
| substitution / change of variable in an integral or average | C02 | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| even and odd functions | C02 | knowledge/primers.md: even and odd functions (ch11 P261) — reminder |
| correlation by FFT with zero padding | C02 | primer (in C02) |
| integral time scale Λ_t (12.18); equal-area rectangle | C02 | C02 (N31, N34) |
| correlation time t_c; independent samples N ≈ Δt/t_c | C02 | C02 (N32) |
| delayed copy: cross-correlation peaks at the delay | C02 | C02 (N30) |
| Taylor series to second order | C02 | knowledge/primers.md: Taylor expansion (ch01 P26; ch03 P98) — reminder |
| product rule | C02 | knowledge/primers.md: product rule (ch01 P38) — reminder |
| Taylor microscale λ_t (12.19); osculating parabola | C02 | C02 (D04; N33) |
| np.trapezoid; Simpson | C02 | knowledge/primers.md: trapezoid rule (ch01 P37; ch09 P203) — reminder |
| energy spectrum S_e(ω) (12.20) | C03 | C03 |
| Fourier-transform pair: where the 2π sits | C03 | primer (in C03) |
| Euler's formula e^{iθ} = cos θ + i sin θ | C03 | knowledge/primers.md: complex exponential (ch01 P45) — reminder |
| improper integrals | C03 | knowledge/primers.md: improper integral as a limit (ch05 P144) — reminder |
| inverse transform (12.21); variance integral (12.22); S_e(0) = variance × Λ_t/π | C03 | C03 (D05; N35–N37) |
| the discrete Fourier transform as a Riemann sum; Parseval | C03 | primer (in C03) |
| np.fft.rfft, Fourier modes | C03 | knowledge/primers.md: Fourier modes and the FFT (ch05 P142); np.fft.rfft (ch07 P181) — reminder |
| segment averaging (Welch) and leakage | C03 | primer (in C03) |
| np.corrcoef (numpy's own correlation-coefficient matrix; its off-diagonal entry is r) | C02 | C02 (gloss where `TS.correlation_coefficient` is compared with it, N27–N28 row) |
| np.isclose (the one-number twin of np.allclose) | C02 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| periodogram of a finite record (Wiener–Khinchin) | C03 | C03 (N38; gloss) |
| spatial correlation tensor R_ij(r) (12.23) | C03 | C03 (N39) |
| change of variable in a density; units of a spectral density | C03 | primer (in C03) |
| Taylor's frozen-turbulence hypothesis; k_1 = ω/U_0, S_11 = U_0 S_e | C03 | C03 (N40) |
| scipy.integrate.quad | C03 | knowledge/primers.md: quad (ch03 P87) — reminder |
| why the mean matters; range of scales too wide to resolve | C04 | C04 (N41) |
| Reynolds decomposition (12.24)–(12.26) | C04 | C04 (N42–N44) |
| Boussinesq equations in flux form | C04 | R01 (Ch. 4 §4.9) |
| summation convention, δ_ij, free vs dummy index | C04 | Ch. 2 §2.1, §2.7 (reminder in A.0 row 6) |
| a sympy averaging operator | C04 | primer (in C04) |
| mean continuity (12.27); fluctuation continuity (12.28); (12.29) | C04 | C04 (D06; N45–N48) |
| Reynolds-averaged momentum equation (12.30); mean stress tensor | C04 | C04 (D06) |
| viscous term as a stress divergence, 2μS_ij | C04 | Ch. 4 §4.5–4.6 (reminder in D06 step 11) |
| Reynolds stress tensor: symmetric, six components, normal and shear | C04 | C04 (N49) |
| counting the independent components of a symmetric tensor | C04 | primer (in C04) |
| principal axes of a symmetric tensor | C04 | Ch. 2 §2.11 (reminder in P287) |
| sign convention of a stress on a face | C04 | R02 (Ch. 2 §2.4) |
| displaced-parcel argument: ⟨uv⟩ < 0 for dU/dy > 0 | C04 | C04 (N50) |
| Reynolds stress as a momentum flux | C04 | C04 (N51) |
| mean temperature equation (12.31), heat flux (12.32), turbulent heat flux | C04 | C04 (N52–N54) |
| Fourier's and Fick's laws | C04 | Ch. 1 §1.5 (reminder in N53, N56) |
| kinematic vs dynamic fluxes | C04 | primer (in C04) |
| passive scalar, mixture density; scalar equations (12.33), (12.34) | C04 | C04 (N55–N57) |
| RANS equations; constant-density form | C04 | C04 (N58) |
| Reynolds-stress transport equation (12.35) | C04 | C04 (N59) |
| closure problem; RANS models, DNS, LES | C04 | C04 (N60) |
| no Bernoulli integral with Reynolds stresses | C04 | C04 (N61) |
| caching slow symbolic results | C04 | knowledge/primers.md: caching runs (ch10 P252) — reminder |
| homogeneous isotropic turbulence; grid turbulence; local isotropy | C05 | C05 (N62; gloss) |
| isotropy statements (12.36), (12.37) | C05 | C05 (N64, N65) |
| longitudinal and transverse correlations f, g (12.38); scales (12.39) | C05 | C05 (N66, N67, N71) |
| isotropic tensors (built from δ_ij and a vector) | C05 | knowledge/primers.md: isotropic tensors (ch04 P119, P120) — reminder |
| derivatives of functions of r = ∣r∣ | C05 | primer (in C05) |
| isotropic tensor (12.40), incompressible form (12.41), g = f + (r/2)f′ | C05 | C05 (D07; N68–N70) |
| integration by parts | C05 | knowledge/primers.md: integration by parts (ch09 P218a) — reminder |
| turbulent kinetic energy ē = ½⟨u_i²⟩; ⟨u²⟩ is one component | C05 | C05 (N72) |
| mean dissipation rate (12.42) | C05 | C05 (N73); Ch. 4 §4.8 for (4.58) |
| gradient moments from the correlation tensor; 2 : 4 : −1 | C05 | C05 (D08; N74) |
| isotropic dissipation (12.43) | C05 | C05 (D08) |
| Taylor-scale Reynolds number (12.44) | C05 | C05 (N75) |
| one-dimensional wavenumber spectrum S_11(k_1) (12.45) | C05 | C05 (N76) |
| measurement recipe: record → frozen field → R_11 → S_11 | C05 | C05 (N77) |
| np.gradient | C05 | knowledge/primers.md: np.gradient (ch01 P22) — reminder |
| mean-flow kinetic-energy budget (12.46) | C06 | C06 (D09; N78) |
| symmetric contraction τ_ij ∂U_i/∂x_j = τ_ij S_ij | C06 | Ch. 2 §2.10 (reminder in D09 step 7) |
| Gauss's theorem: divergence terms only transport | C06 | Ch. 2 §2.12 (reminder in N79) |
| turbulent kinetic-energy budget (12.47) | C06 | C06 (D10) |
| shear production, dissipation ε̄, buoyant production | C06 | C06 (N81) |
| direct dissipation of the mean flow ~ 1/Re of production | C06 | C06 (N80) |
| isotropic turbulence has no shear production | C06 | C06 (N83) |
| mixing an unstable layer releases potential energy | C06 | C06 (N82) |
| disturbance-energy equation (11.88) | C06 | Ch. 11 §11.10 (reminder in N81) |
| outer scales ΔU, L | C07 | C07 (N84) |
| order-of-magnitude scaling ("∼") | C07 | knowledge/primers.md: order-of-magnitude scaling (ch04 P130) — reminder |
| ε̄ ~ (ΔU)³/L (12.48), (12.49) | C07 | C07 (D11; N85, N86) |
| the cascade; u′(l′) ~ (ε̄ l′)^{1/3} | C07 | C07 (N87) |
| vortex stretching | C07 | Ch. 5 §5.4 (reminder in N01, N87) |
| Π theorem and exponent matching | C07 | Ch. 1 §1.11 (reminder in D11; `DIM.solve_exponents`) |
| np.linalg.solve | C07 | knowledge/primers.md: np.linalg.solve (ch01 P57) — reminder |
| Kolmogorov scales (12.50); ηu_K/ν = 1 | C07 | C07 (D11; N88) |
| scale separation (12.51); DNS cost Re^{9/4} | C07 | C07 (D11; N89) |
| real cases: millimetre η in atmosphere and ocean | C07 | C07 (N90) |
| Taylor microscale between η and L (12.52); R_λ ~ Re_L^{1/2}; ordering of scales | C07 | C07 (N91, N92) |
| pandas DataFrame for small tables | C07 | front matter (gloss in A.0 row 5: a DataFrame is a labelled table) and the gloss list of knowledge/primers.md |
| universal spectrum (12.53) | C08 | C08 (D12; N93) |
| Kolmogorov's −5/3 law (12.54), slip #1 | C08 | C08 (D12) |
| one- vs two-sided normalisation (12.55) | C08 | C08 (D12; N94) |
| three-dimensional spectrum; C_1 = (18/55) C | C08 | C08 (N95) |
| quadrature on a logarithmic grid | C08 | primer (in C08) |
| model spectrum (Pao; Pope) | C08 | C08 (N97) |
| collapse of spectra in Kolmogorov scaling | C08 | C08 (N96) |
| live widgets (kernel only) | C08 | knowledge/primers.md: live widgets (ch01 P47) — reminder |
| free vs wall-bounded shear flow; jet, wake, shear layer, plume | C09 | C09 (N100) |
| entrainment; self-preservation | C09 | C09 (N101) |
| similarity forms (12.56), (12.57) | C09 | C09 (N102, N103) |
| mean equations (12.58)–(12.60) | C09 | C09 (N104–N106) |
| two-length (thin-layer) scaling | C09 | knowledge/primers.md: anisotropic scaling (ch08 P188) — reminder |
| thin-layer equations (12.61); invariant J_s (12.62) | C09 | C09 (D13; N107–N109) |
| chain rule with a moving similarity variable | C09 | knowledge/primers.md: similarity-variable chain rule (ch09 P206) — reminder |
| Leibniz rule | C09 | knowledge/primers.md: Leibniz rule (ch08 P189) — reminder |
| similarity equation (12.63) | C09 | C09 (D14; N110, N111) |
| separation argument: a function of x equal to a function of ξ is a constant | C09 | knowledge/primers.md: separation of variables (ch07 P167) — reminder |
| exponent matching | C09 | knowledge/primers.md: exponent matching (ch08 P197) — reminder |
| (12.64), virtual origin, (12.65), far field (12.66), (12.67), volume flux (12.68) | C09 | C09 (D15; N112–N116) |
| Gaussian profile and Gaussian integrals | C09 | C09 (N122); knowledge/primers.md: Gaussian integral (ch08 P194) — reminder |
| scalar field of the jet (12.69)–(12.73) | C09 | C09 (N117–N121) |
| exponents of Table 12.1 from invariant + growth law | C09 | C09 (D16; N123, N124) |
| fractions.Fraction | C09 | knowledge/primers.md: fractions.Fraction (ch01 P60) — reminder |
| general similarity (12.74), slip #5 | C09 | C09 (N125) |
| stoichiometric mass fraction (Example 12.2, our numbers) | C09 | C09 (N126; gloss on mole vs mass fraction) |
| jet energy budget (12.75) and its three balances | C09 | C09 (N127–N129) |
| constant eddy viscosity gives sech² | C09 | C09 (N130) |
| laminar jet as a ghost (δ ∝ x^{2/3}, U ∝ x^{−1/3}) | C09 | Ch. 9 §9.10 (reminder in F4) |
| two length scales at a wall; no Re independence | C10 | C10 (N131) |
| turbulent vs laminar mean profiles | C10 | C10 (N132) |
| fully developed flow; force balance on a slug | C10 | R03 (Ch. 8 §8.2) |
| channel equations (12.76), (12.77); linear total stress | C10 | C10 (D17; N133–N136) |
| pipe balance (12.91) | C10 | R04 (Ch. 8 §8.2) |
| stress partition; constant-stress layer; boundary-layer equation (12.78) | C10 | C10 (N137, N138) |
| inner, outer and overlap: two descriptions that must agree where both hold | C10 | primer (in C10) |
| inner layer, outer layer, overlap region | C10 | C10 (N139) |
| law of the wall (12.80); u_* (12.81); l_ν, y⁺, U⁺, δ⁺ | C10 | C10 (D18; N140, N141) |
| viscous sublayer U⁺ = y⁺ (12.82) | C10 | C10 (D18; N142) |
| semi-log axes: a logarithm is a straight line | C10 | primer (in C10) |
| scaled (dimensionless) variables | C10 | knowledge/primers.md: scaled variables and the chain rule (ch04 P133) — reminder |
| defect law (12.83), (12.84) | C11 | C11 (D19; N143, N144) |
| chain rule | C11 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| overlap matching (12.85)–(12.87) | C11 | C11 (D19; N145–N147) |
| logarithmic law (12.88), (12.89); κ, B | C11 | C11 (D19; N148) |
| composite profile: inner law + outer correction | C11 | primer (in C11) |
| Spalding's formula; Coles' wake | C11 | C11 (N153, N154) |
| scipy.optimize.brentq | C11 | knowledge/primers.md: brentq (ch03 P108) — reminder |
| layer names: sublayer, buffer, logarithmic, wake | C11 | C11 (N149) |
| log-law constants by flow; (12.92) | C11 | C11 (N156, N157) |
| rough-wall law (12.93); roughness length; neutral drag coefficient | C11 | C11 (D20; N158, N159) |
| zero-pressure-gradient fits; skin-friction laws; Example 12.3 (our numbers) | C11 | C11 (N150–N152, N161) |
| displacement and momentum thickness | C11 | Ch. 9 §9.2 (reminder in N150) |
| pipe friction law from the log law; Darcy friction factor | C11 | C11 (N160) |
| open questions; universality | C11 | C11 (N155) |
| eddy viscosity, eddy diffusivity; turbulent-viscosity hypothesis (12.94) | C12 | C12 (N162) |
| gradient diffusion (12.95), (12.96); effective-viscosity RANS (12.97), slip #6 | C12 | C12 (N163–N165) |
| why the molecular analogy is imperfect (Knudsen number contrast) | C12 | C12 (N166); Ch. 1 §1.4 |
| ν_T ~ l_T u_T (12.98) | C12 | C12 (N167) |
| mixing length; (12.99)–(12.101); exact slope; intercept | C12 | C12 (D21; N168–N172) |
| van Driest damping | C12 | C12 (gloss in the code explanation; N181) |
| scipy.integrate.cumulative_trapezoid | C12 | knowledge/primers.md: cumulative_trapezoid (ch08 P190) — reminder |
| a nonlinear diffusion problem by Picard iteration on the eddy viscosity | C12 | primer (in C12) |
| mixing-length channel against DNS | C12 | C12 (N181) |
| stretched grid | C12 | Ch. 10 §10.2 (`fd.stretched_grid`; reminder in N181) |
| reading reference data files | C12 | knowledge/primers.md: reading reference data (ch10 P251) — reminder |
| convective velocity and eddy diffusivity (12.102) | C12 | C12 (N173, N174) |
| one-equation ingredients | C13 | C13 (N175) |
| modelled energy equation (12.103); ν_T = C_μ ē²/ε̄ (12.104) | C13 | C13 (D22; N176) |
| modelled dissipation equation (12.105); the five constants | C13 | C13 (N177, N178) |
| wall functions; Reynolds-stress closures | C13 | C13 (N179; gloss) |
| dividing two ODEs to eliminate time | C13 | primer (in C13) |
| decay exponent n = 1/(C_ε2 − 1); log-layer relation for κ | C13 | C13 (D23; N180) |
| solve_ivp; RK4 loop | C13 | knowledge/primers.md: solve_ivp (ch01 P31; ch03 P94), RK4 (ch03 P95) — reminder |
| potential temperature; two lapse-rate conventions; N² | C14 | R05 (Ch. 1 §1.10) |
| inequalities under a sign change | C14 | knowledge/primers.md: multiplying an inequality by −1 (ch01 P48) — reminder |
| reduced budget (12.106) | C14 | C14 (D24; N182) |
| flux Richardson number (12.107); Rf_cr ≈ ¼ | C14 | C14 (D24; N183) |
| gradient Richardson number (12.108) | C14 | R06 (Ch. 11 §11.7) |
| necessary vs sufficient (theorem vs observation) | C14 | knowledge/primers.md: necessary vs sufficient (ch11 P255) — reminder |
| Ri = (ν_T/κ_T) Rf (12.109); turbulent Prandtl number; Reynolds analogy | C14 | C14 (D24; N184, N185) |
| internal waves carry momentum but not heat | C14 | Ch. 7 §7.8 (reminder in N185) |
| Monin–Obukhov length (12.110) | C15 | C15 |
| Rf = z/L_M (12.111); forced vs free convection | C15 | C15 (D25; N186, N187) |
| the stability parameter ζ = z/L and integrating a flux–profile relation | C15 | primer (in C15) |
| log-linear wind profile; Businger–Dyer form | C15 | C15 (D25; N188, N189; gloss) |
| temperature-variance budget (12.112), slip #15 | C15 | C15 (N190) |
| temperature spectrum; Obukhov–Corrsin (12.113); Batchelor scale; (12.114) | C15 | C15 (N191–N195) |
| Lagrangian description; particle paths | C16 | C16 (N196); Ch. 3 §3.2 |
| Lagrangian autocorrelation r_α | C16 | C16 (N199) |
| a double integral over a triangle | C16 | primer (in C16) |
| iterated integrals | C16 | knowledge/primers.md: volume and surface integrals as sums (ch02 P83); quad and dblquad (ch03 P87) — reminder |
| (12.115)–(12.118); Taylor's formula (12.119) | C16 | C16 (D26; N197, N198, N200–N202) |
| ballistic and diffusive limits (12.120)–(12.123); exponential and Gaussian closed forms | C16 | C16 (D27; N203–N208) |
| erf | C16 | knowledge/primers.md: error function (ch08 P195) — reminder |
| proof by induction | C16 | primer (in C16) |
| random walk (12.124), (12.125) | C16 | C16 (N209–N211) |
| smoke plume width, t = x/U; slip #13 | C16 | C16 (N212) |
| Gaussian spreading of a line vortex, σ² = 2νt | C16 | R07 (Ch. 8; Ch. 3; Ch. 5) |
| diffusivity from the growth of a variance (12.126) | C16 | R08 (Ch. 1 §1.5) |
| differentiation under the integral sign | C16 | knowledge/primers.md: differentiation under the integral sign (ch03 P109) — reminder |
| eddy diffusivity D_T (12.127)–(12.129), slip #12 | C16 | C16 (D28; N213–N215) |
| turbulent vs molecular diffusion; Richardson's 4/3 law | C16 | C16 (N216, N06) |
| concluding remarks; pointers to Ch. 13 | C16 | C16 (N217) |
| exercises; literature | C16 | S01, S02 (pointers in §12.13) |

---

## Part F — derivation storyboards

One block per `D` row of the curation (§4b), in order. Builders copy each block word for word into `nb.derivation(id, title,
ref=…, goal=…, start=(tex, plain), plan=[…], uses=[…], steps=[dict(did=…, tex=…, why=…, plain=…), …], result=(tex, plain),
interpret=…, check=…, check_src=…)` and, where the heading names an explainer, into that explainer's `derivations: [...]` (same
`did` titles, same step count). Every step has four parts: **did** (the move) · **tex** (the new line) · **why** (why it is
allowed and why we make it) · **plain** (what the line says). An over-bar is the ensemble average; $\tilde u=U+u$. Wherever a
step uses a book equation the equation is written out with its number. Colours of terms as convention 9. ★★★ blocks carry a
**sympy check intent** with a `check_src` sketch that re-runs the derivation's own construction (ch01 lesson).

**Repeat mentions inside a block: completed (inlining pass of 2026-10-07).** Every block writes each equation out where it
is first used, and every later line of the block that names an equation again — a *why*, a Start, Tools, Result, Check or
Traps line — now either carries the compact equation beside the number or, where the equation is long, points in words to
the block's own Start or Result ("the Start's ē equation", "D14's Result"). A Result line that used to read "(number) as in
step n" now holds the equation itself, which is what `nb.derivation(result=(tex, plain))` needs. Builders copy the lines as
they stand; nothing is left for them to inline. Two conventions remain, both deliberate:

- **Own-line labels.** A *why* that ends "This is (N.M)." — or "First of (N.M).", "These are (N.M) and (N.M)." — labels the
  step's own `tex` line: that line *is* the numbered equation and sits directly above the sentence in the notebook cell and
  on the Derivation page. There are 73 of these. The equation is not repeated in the prose, so that the *why* stays inside
  its length limit.
- **Explainer copies.** In an explainer's plain-text *why* the builder writes an inlined equation in Unicode (no raw TeX).
  If that takes a phone *why* past 35 words, the first sentence stays and the number-plus-equation moves to the step's
  `tex` (or `live`) field — never the number alone.

### D01 · The rules of averaging (12.4)–(12.9), $\overline{\bar u}=\bar u$, and $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$ — ★, 8 steps, in C01 (notebook · `reynolds_averaging_window`)
- **Goal:** find out which operations may be swapped with an average — and the one that may not. Every later derivation of the
  chapter uses only these rules. · **Start:** $\overline{u^m(\mathbf x,t)}\equiv\frac1N\sum_{n=1}^N\big(u(\mathbf x,t{:}n)\big)^m$, the finite-N form of $\langle u^m(\mathbf x,t)\rangle=\lim_{N\to\infty}\frac1N\sum_{n=1}^N\big(u(\mathbf x,t{:}n)\big)^m$ (12.1) — "add the N
  realizations and divide by N". · **Plan:** • write the average as a plain finite sum • pass each linear operation through the
  sum • try a product and see what is left • take N → ∞ last. · **Tools:** finite sums (school algebra); partial derivative
  (P25, reminder); integral as a limit of sums (P84, reminder); expanding a product. · **Assumptions:** ensemble average; the
  same N realizations for every variable; the fields are differentiable.
- **Steps:**
  1. **did** Write the average as a sum · **tex** $\overline{u}=\dfrac1N\big(u_1+u_2+\dots+u_N\big)$, $u_n\equiv u(\mathbf x,t{:}n)$ · **why** This is the
     definition with m = 1; writing it out shows an average is nothing more than adding and dividing by a fixed number. ·
     **plain** The mean is the sum of the runs divided by how many there are.
  2. **did** Average a sum of two variables · **tex** $\overline{u+v}=\dfrac1N\sum_n(u_n+v_n)=\dfrac1N\sum_nu_n+\dfrac1N\sum_nv_n=\overline u+\overline v$ · **why** A finite
     sum can be added in any order (commutative, associative). This is (12.4), $\overline{u^m+v^m}=\overline{u^m}+\overline{v^m}$, for m = 1; any m works
     the same way. · **plain** The mean of a sum is the sum of the means.
  3. **did** Average a constant times a variable · **tex** $\overline{Au}=\dfrac1N\sum_nAu_n=A\,\dfrac1N\sum_nu_n=A\,\overline u$ · **why** A is the same in every
     realization, so it factors out of the sum (distributive law). This is (12.5), $\overline{Au^m}=A\overline{u^m}$; with u = 1 it gives
     $\overline A=A$. · **plain** A fixed number passes straight through an average.
  4. **did** Average a derivative · **tex** $\overline{\dfrac{\partial u}{\partial t}}=\dfrac1N\sum_n\dfrac{\partial u_n}{\partial t}=\dfrac{\partial}{\partial t}\Big(\dfrac1N\sum_nu_n\Big)=\dfrac{\partial\overline u}{\partial t}$ · **why** The derivative of a
     finite sum is the sum of the derivatives (linearity of ∂), and 1/N is a constant. This is (12.6); with $x_j$ for t it is
     (12.8), $\overline{\partial u^m/\partial x_j}=\partial\overline{u^m}/\partial x_j$. · **plain** Averaging and differentiating can be done in either order.
  5. **did** Average an integral · **tex** $\overline{\displaystyle\int_a^bu\,dt}=\dfrac1N\sum_n\int_a^bu_n\,dt=\int_a^b\Big(\dfrac1N\sum_nu_n\Big)dt=\int_a^b\overline u\,dt$ · **why** An integral is a
     limit of sums, so it is linear too and swaps with the finite sum over n. This is (12.7); the spatial version is (12.9),
     $\overline{\int u^m\,d\mathbf x}=\int\overline{u^m}\,d\mathbf x$. · **plain** Averaging and integrating can be done in either order.
  6. **did** Average something already averaged · **tex** $\overline{\bar u}=\bar u,\qquad\overline{\bar u\,v}=\bar u\,\overline v$ · **why** $\bar u$ is one number at each (x, t), the
     same in every realization — a "constant" in the sense of step 3. We need this to handle products next. · **plain** A mean
     is not random any more, so averaging it again changes nothing.
  7. **did** Average a product of two total fields · **tex** $\overline{\tilde u\tilde v}=\overline{(\bar u+u)(\bar v+v)}=\bar u\bar v+\bar u\,\overline v+\overline u\,\bar v+\overline{uv}=\bar u\bar v+\overline{uv}$ · **why** Expand
     the product, then use steps 2, 3 and 6; the two middle terms vanish because a fluctuation has zero mean, $\overline u=\overline v=0$ (C04's $\overline{u_i}=0$ (12.26)). A product is **not** linear, so nothing lets us split $\overline{uv}$. · **plain** The mean of a product is
     the product of the means plus the covariance of the fluctuations.
  8. **did** Let the number of realizations grow · **tex** $\langle\tilde u\tilde v\rangle=\lim_{N\to\infty}\overline{\tilde u\tilde v}=\langle\tilde u\rangle\langle\tilde v\rangle+\langle uv\rangle,\qquad\langle uv\rangle\neq0\ \text{in general}$ · **why** Every line
     above holds for each finite N, so it holds in the limit; taking the limit last avoids swapping a limit with a
     derivative. · **plain** However many runs we average, the leftover covariance stays.
- **Result:** $\overline{u+v}=\bar u+\bar v$, $\overline{Au}=A\bar u$, $\overline{\partial u/\partial t}=\partial\bar u/\partial t$, $\overline{\partial u/\partial x_j}=\partial\bar u/\partial x_j$, $\overline{\int u}=\int\bar u$, but
  $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$ — "an average commutes with everything linear and with nothing else".
- **Check:** units — both sides of every rule carry the units of u (or u/t, u·t, uv). Number: two runs ũ = (1, 3), ṽ = (2, 6):
  $\bar u\bar v$ = 8, $\overline{uv}$ = 2, $\overline{\tilde u\tilde v}$ = 10 ✓. Code: `TS.check_averaging_rules` returns round-off for every entry except
  "product". Special case: if u and v are uncorrelated, $\overline{uv}=0$ and the product rule "works".
- **What it means:** when we average the equations of motion every linear term keeps its form; the nonlinear term
  $\tilde u_j\tilde u_i$ leaves $\overline{u_iu_j}$ behind — the Reynolds stress. **Fails when:** the "average" is a finite time window: then $\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t$ (12.6) holds only if the window is short against the drift of the mean (C01, N14).
- **Traps:** thinking the rules are special to turbulence (they are linearity); forgetting that $\overline{u^2}\neq\bar u^2$ is the same
  failure; using different sets of runs for u and v; commuting a *time* average with ∂/∂t exactly.

### D02 · The Schwartz inequality $\lvert\overline{uv}\rvert\le\sqrt{\overline{u^2}}\sqrt{\overline{v^2}}$ (12.16) and $-1\le r\le1$ — ★★, 6 steps, in C02 (notebook)
- **Goal:** show that a correlation coefficient can never leave the range −1 … 1 — so "r = 0.9" means the same in every flow. ·
  **Start:** $\overline{(u+\lambda v)^2}\ge0$ for every real number λ — "the average of a square cannot be negative". · **Plan:** • expand
  the square into a quadratic in λ • a quadratic that is never negative has discriminant ≤ 0 • rearrange. · **Tools:** the rules
  of D01; discriminant of a quadratic (P288, primer above; quadratic formula P159). · **Assumptions:** finite variances; u and
  v are fluctuations at any two points and times.
- **Steps:**
  1. **did** Average a square that contains a free number λ · **tex** $\overline{(u+\lambda v)^2}\ge0\quad\text{for every real }\lambda$ · **why** Each realization
     contributes a square, which is ≥ 0, so their average is ≥ 0. λ is a device: one inequality for every λ is a strong
     statement. · **plain** No mixture of u and v can have a negative mean square.
  2. **did** Expand and average term by term · **tex** $\overline{v^2}\,\lambda^2+2\,\overline{uv}\,\lambda+\overline{u^2}\ge0$ · **why** $(u+\lambda v)^2=u^2+2\lambda uv+\lambda^2v^2$; λ is a constant, so D01's rules $\overline{u^m+v^m}=\overline{u^m}+\overline{v^m}$ (12.4) and $\overline{Au^m}=A\overline{u^m}$ (12.5) apply. · **plain** The inequality is a parabola in λ that never dips below zero.
  3. **did** Apply the discriminant condition · **tex** $(2\,\overline{uv})^2-4\,\overline{v^2}\;\overline{u^2}\le0$ · **why** A parabola $a\lambda^2+b\lambda+c$ with a > 0 that is
     never negative has at most one real root, which requires $b^2-4ac\le0$ (P288). · **plain** The cross term cannot be too large
     compared with the two variances.
  4. **did** Rearrange and take the square root · **tex** $\lvert\overline{uv}\rvert\le\sqrt{\overline{u^2}}\,\sqrt{\overline{v^2}}$ · **why** Divide by 4, move one term across, take the positive
     root of both (non-negative) sides; $\sqrt{x^2}=\lvert x\rvert$ — hence the absolute value, which the printed $\overline{uv}\le\sqrt{\overline{u^2}}\sqrt{\overline{v^2}}$ (12.16) omits. · **plain** A
     covariance is never bigger in size than the product of the two rms values.
  5. **did** Divide by the right-hand side · **tex** $-1\le r\equiv\dfrac{\overline{uv}}{\sqrt{\overline{u^2}}\sqrt{\overline{v^2}}}\le1$ · **why** Both rms values are positive, so dividing keeps the
     inequality; this is the coefficient of (12.14), $r_{12}=\overline{u_1u_2}/(\sqrt{\overline{u_1^2}}\sqrt{\overline{u_2^2}})$. · **plain** Every correlation coefficient
     lies between −1 and +1.
  6. **did** Find when equality holds · **tex** $r=\pm1\iff v=c\,u\ \text{in every realization}$ · **why** Equality means the parabola touches zero at
     some $\lambda_0$: $\overline{(u+\lambda_0v)^2}=0$, which forces $u+\lambda_0v=0$ in every realization. · **plain** Perfect correlation means one signal is
     a fixed multiple of the other.
- **Result:** $\lvert\overline{uv}\rvert\le\sqrt{\overline{u^2}}\sqrt{\overline{v^2}}$ (12.16), so $-1\le r\le1$ — "a correlation coefficient is a pure number between −1 and 1".
- **Check:** units — both sides u·v. Number: u = (1, −1, 2, −2), v = (2, 0, 1, −3): $\overline{uv}$ = 2.5, bound √2.5 × √3.5 = 2.958, r =
  0.845 ✓. Limits: v = u gives r = 1; v = −u gives −1. Code: `TS.correlation_coefficient` vs `np.corrcoef`.
- **What it means:** the Reynolds shear stress can never exceed $\rho\,u_{rms}v_{rms}$; in shear flows it is about 0.4 of that. The
  same inequality is the Cauchy–Schwarz inequality for vectors, $\lvert\mathbf a\cdot\mathbf b\rvert\le\lvert\mathbf a\rvert\lvert\mathbf b\rvert$, with "average of a product" as the
  dot product.
- **Traps:** writing the discriminant condition as ≥ 0 (that is the condition for two real roots — a parabola that *does*
  cross zero); dropping the absolute value, as the page does (then r ≥ −1 would not follow).

### D03 · Stationary lag forms (12.17), $R_{11}(\tau)=R_{11}(-\tau)$ and $R_{ij}(-\tau)=R_{ji}(\tau)$ — ★, 5 steps, in C02 (notebook · `reynolds_averaging_window` · `correlation_and_spectrum`)
- **Goal:** show that for a stationary signal the correlation depends only on the time *difference*, and that the
  autocorrelation is the same for a lag forward or backward. · **Start:** $R_{11}(t_1,t_2)\equiv\overline{u_1(t_1)u_1(t_2)}$ at one point — $R_{11}(\mathbf x_1,t_1,\mathbf x_2,t_2)\equiv\overline{u_1(\mathbf x_1,t_1)u_1(\mathbf x_2,t_2)}$ (12.13) with $\mathbf x_1=\mathbf x_2$. · **Plan:** • use stationarity to drop the absolute time • flip the sign of the lag • shift the time origin.
  · **Tools:** substitution / change of variable (P106); stationarity (N08). · **Assumptions:** temporally stationary; same
  point in space.
- **Steps:**
  1. **did** Name the two times by a start and a lag · **tex** $R_{11}(t,\,t+\tau)=\overline{u_1(t)\,u_1(t+\tau)}$ · **why** Any pair of times can be written as
     a time t and a difference τ = t₂ − t₁; nothing has been assumed yet. · **plain** We describe the pair by "when" and "how
     far apart".
  2. **did** Use stationarity · **tex** $R_{11}(\tau)=\overline{u_1(t)\,u_1(t+\tau)}\quad\text{(the same for every }t)$ · **why** Stationary means the statistics do not
     change when the time origin is moved, so the average cannot depend on t. This is (12.17). · **plain** Only the lag
     matters.
  3. **did** Replace τ by −τ · **tex** $R_{11}(-\tau)=\overline{u_1(t)\,u_1(t-\tau)}$ · **why** Step 2 holds for any real lag, negative ones included; we want
     to compare the two. · **plain** Looking backward by τ instead of forward.
  4. **did** Shift the time origin by τ · **tex** $\overline{u_1(t)\,u_1(t-\tau)}=\overline{u_1(t'+\tau)\,u_1(t')}=R_{11}(\tau),\quad t'=t-\tau$ · **why** Substituting t = t′ + τ is allowed
     because, by step 2, the average does not care what the start time is called; the two factors commute. · **plain** The
     autocorrelation is an even function of the lag.
  5. **did** Repeat for two different components · **tex** $R_{ij}(-\tau)=\overline{u_i(t)\,u_j(t-\tau)}=\overline{u_i(t'+\tau)\,u_j(t')}=R_{ji}(\tau)$ · **why** The same shift; but now
     swapping the factors swaps the *roles* of i and j, so we get $R_{ji}$, not $R_{ij}$. · **plain** A cross-correlation is not
     even: mirroring the lag swaps which signal leads.
- **Result:** $R_{11}(\tau)=\overline{u_1(t)u_1(t+\tau)}=R_{11}(-\tau)$ (12.17); $R_{ij}(-\tau)=R_{ji}(\tau)$ — "an autocorrelation is even; a cross-correlation
  mirrors into its transpose".
- **Check:** units u². τ = 0 gives the variance $\overline{u_1^2}$. Example: $R=\sigma^2e^{-\lvert\tau\rvert/\tau_c}$ is even ✓. Code: `TS.cross_correlation` of a
  signal and its copy delayed by 0.8 s peaks at +0.8 s one way round and −0.8 s the other.
- **What it means:** we only ever need r(τ) for τ ≥ 0 — which is why Λ_t integrates from 0 and why the spectrum (D05) is real.
  **Fails when:** the signal is not stationary (a decaying flow): then R depends on t as well.
- **Traps:** assuming a cross-correlation is even; forgetting that the shift is legitimate *only* because of stationarity.

### D04 · The Taylor microscale $\lambda_t^2\equiv-2/[d^2r_{11}/d\tau^2]_{\tau=0}$ (12.19) as the foot of the osculating parabola; $\lambda_t^2=2\overline{u^2}/\overline{(du/dt)^2}$ — ★★, 7 steps, in C02 (notebook · `correlation_and_spectrum`)
- **Goal:** give the definition $\lambda_t^2\equiv-2\big/\big[d^2r_{11}/d\tau^2\big]_{\tau=0}$ (12.19) a picture (a parabola fitted to the top of the correlation curve) and a meaning (how
  fast the signal changes). · **Start:** $r_{11}(\tau)$, even, with $r_{11}(0)=1$ (D03). · **Plan:** • Taylor-expand r about τ = 0 •
  evenness kills the linear term • define λ_t where the parabola reaches zero • relate the curvature to the mean-square time
  derivative. · **Tools:** Taylor series to second order (P26); even functions (P261); product rule (P38); D01's $\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t$ (12.6);
  stationarity. · **Assumptions:** r twice differentiable at 0 (a smooth signal — **not** an Ornstein–Uhlenbeck one).
- **Steps:**
  1. **did** Expand r about zero lag · **tex** $r_{11}(\tau)=1+r_{11}'(0)\,\tau+\tfrac12r_{11}''(0)\,\tau^2+\dots$ · **why** Taylor's theorem for a twice-differentiable
     function; $r_{11}(0)=1$ by its definition. We want the shape of the curve's top. · **plain** Near zero lag the curve is a
     constant plus a slope plus a bend.
  2. **did** Use evenness · **tex** $r_{11}'(0)=0$ · **why** An even, differentiable function has zero slope at the origin: $r(\tau)=r(-\tau)$
     differentiated gives $r'(\tau)=-r'(-\tau)$, so $r'(0)=-r'(0)$. · **plain** The top of the curve is flat.
  3. **did** Keep the first two terms · **tex** $r_{11}(\tau)\approx1+\tfrac12r_{11}''(0)\,\tau^2,\qquad r_{11}''(0)<0$ · **why** For small τ higher terms are negligible (≈). The
     curvature is negative because r has its maximum, 1, at τ = 0 (D02). · **plain** The top of the curve is a downward
     parabola.
  4. **did** Call λ_t the lag where this parabola reaches zero · **tex** $0=1+\tfrac12r_{11}''(0)\lambda_t^2\ \Rightarrow\ \lambda_t^2=-\dfrac{2}{r_{11}''(0)},\qquad r_{11}\approx1-\dfrac{\tau^2}{\lambda_t^2}$ · **why** This
     is the definition (12.19); the minus sign makes $\lambda_t^2$ positive because $r''(0)<0$. · **plain** λ_t is where the fitted
     parabola would hit the axis.
  5. **did** Differentiate the correlation twice in the lag · **tex** $R_{11}''(\tau)=\overline{u(t)\,\dfrac{d^2u}{dt^2}(t+\tau)}$ · **why** τ appears only in the second factor;
     the derivative passes inside the average by D01's $\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t$ (12.6). · **plain** The curvature of R is a correlation between the
     signal and its own acceleration.
  6. **did** Set τ = 0 and use the product rule · **tex** $\overline{u\,\ddot u}=\dfrac{d}{dt}\overline{u\,\dot u}-\overline{\dot u^{\,2}}=-\overline{\dot u^{\,2}}$ · **why** $\frac{d}{dt}(u\dot u)=\dot u^2+u\ddot u$; and
     $\overline{u\dot u}=\tfrac12\frac{d}{dt}\overline{u^2}=0$ because a stationary signal has constant variance. · **plain** The curvature at the top is minus
     the mean-square rate of change.
  7. **did** Combine with step 4 · **tex** $r_{11}''(0)=-\dfrac{\overline{(du/dt)^2}}{\overline{u^2}}\ \Rightarrow\ \lambda_t^2=\dfrac{2\,\overline{u^2}}{\overline{(du/dt)^2}}$ · **why** Divide step 6 by $R_{11}(0)=\overline{u^2}$ to get r″(0)
     and substitute. · **plain** The microscale is the signal's size divided by how fast it changes.
- **Result:** $\lambda_t^2=-2/r_{11}''(0)=2\overline{u^2}/\overline{(du/dt)^2}$ — "λ_t measures the fastest wiggles, Λ_t the memory".
- **Check:** units — s². Gaussian pair $r=e^{-\tau^2/\tau_c^2}$: $r''(0)=-2/\tau_c^2$, λ_t = τ_c ✓ (`correlation_spectrum_pair("gaussian", …)
  ["lambda_t"]`). Pure cosine $r=\cos\omega\tau$: λ_t = √2/ω.
- **What it means:** with r for the spatial correlation f, the same parabola gives λ_f, and the mean-square *velocity
  gradient* gives the dissipation (C05). **Fails when:** r has a corner at 0: $r=e^{-\lvert\tau\rvert/\tau_c}\approx1-\lvert\tau\rvert/\tau_c$ has no parabola — an
  Ornstein–Uhlenbeck signal has infinite $\overline{\dot u^2}$.
- **Traps:** forgetting the minus sign; keeping a linear term; fitting the parabola over lags that are not small.

### D05 · The spectrum–correlation pair (12.20)–(12.21), $S_e$ real and even, $\overline{u_1^2}=\int S_e\,d\omega$ (12.22) and $S_e(0)=\overline{u_1^2}\Lambda_t/\pi$ — ★★, 9 steps, in C03 (notebook · `correlation_and_spectrum`)
- **Goal:** show that the spectrum is real, that its area is the variance, and that its value at zero frequency is the
  memory time in disguise. · **Start:** $S_e(\omega)\equiv\dfrac1{2\pi}\displaystyle\int_{-\infty}^{+\infty}R_{11}(\tau)\,e^{-i\omega\tau}\,d\tau$ (12.20). · **Plan:** • split the exponential with Euler's
  formula • the odd part integrates to zero • state the inverse and set τ = 0 • set ω = 0. · **Tools:** Euler's formula (P45);
  even and odd functions (P261); improper integrals (P144); the Fourier-transform pair (P290, primer above). ·
  **Assumptions:** stationary; $R_{11}$ decays fast enough to be integrable.
- **Steps:**
  1. **did** Start from the definition · **tex** $S_e(\omega)=\dfrac1{2\pi}\displaystyle\int_{-\infty}^{+\infty}R_{11}(\tau)\,e^{-i\omega\tau}\,d\tau$ · **why** This is (12.20): the forward Fourier transform in
     the book's convention — the 1/2π here, angular frequency ω in rad/s. · **plain** The spectrum weighs the correlation
     against a wave of frequency ω.
  2. **did** Split the complex exponential · **tex** $e^{-i\omega\tau}=\cos\omega\tau-i\sin\omega\tau$ · **why** Euler's formula (P45); it separates a part that is even
     in τ from a part that is odd. · **plain** A complex wave is a cosine and a sine.
  3. **did** Split the integral the same way · **tex** $S_e=\dfrac1{2\pi}\displaystyle\int_{-\infty}^{\infty}R_{11}\cos\omega\tau\,d\tau-\dfrac{i}{2\pi}\int_{-\infty}^{\infty}R_{11}\sin\omega\tau\,d\tau$ · **why** Integration is linear. ·
     **plain** A real part and an imaginary part.
  4. **did** Drop the sine integral · **tex** $\displaystyle\int_{-\infty}^{\infty}R_{11}(\tau)\sin\omega\tau\,d\tau=0$ · **why** $R_{11}$ is even (D03) and sine is odd, so the integrand is odd;
     an odd function integrated over a symmetric range gives zero (P261). · **plain** The spectrum has no imaginary part.
  5. **did** Fold the cosine integral onto τ ≥ 0 · **tex** $S_e(\omega)=\dfrac1\pi\displaystyle\int_0^\infty R_{11}(\tau)\cos\omega\tau\,d\tau=S_e(-\omega)$ · **why** The integrand is even, so the two halves
     are equal: twice the integral from 0. Cosine is even in ω, so S is too. · **plain** The spectrum is real and symmetric:
     only τ ≥ 0 and ω ≥ 0 are needed.
  6. **did** Write the inverse transform · **tex** $R_{11}(\tau)=\displaystyle\int_{-\infty}^{+\infty}S_e(\omega)\,e^{+i\omega\tau}\,d\omega$ · **why** The Fourier inversion theorem (P290): transforming
     back needs the opposite sign in the exponent and no further 2π. This is (12.21); the notebook verifies it on
     $R=\sigma^2e^{-\lvert\tau\rvert/\tau_c}$. · **plain** The correlation can be rebuilt from the spectrum — the same information.
  7. **did** Set the lag to zero · **tex** $\overline{u_1^2}=R_{11}(0)=\displaystyle\int_{-\infty}^{+\infty}S_e(\omega)\,d\omega$ · **why** $e^{0}=1$ and $R_{11}(0)$ is the variance. This is (12.22). · **plain**
     The area under the spectrum is the variance: S says how the variance is shared among frequencies.
  8. **did** Set the frequency to zero in step 5 · **tex** $S_e(0)=\dfrac1\pi\displaystyle\int_0^\infty R_{11}(\tau)\,d\tau$ · **why** cos 0 = 1. We want the low-frequency end of the
     spectrum. · **plain** The spectrum at zero frequency is the area under the correlation.
  9. **did** Use the definition of the integral scale · **tex** $S_e(0)=\dfrac{\overline{u_1^2}}{\pi}\displaystyle\int_0^\infty r_{11}(\tau)\,d\tau=\dfrac{\overline{u_1^2}}{\pi}\,\Lambda_t$ · **why** $R_{11}=\overline{u_1^2}\,r_{11}$ and
     $\Lambda_t\equiv\int_0^\infty r_{11}\,d\tau$ (12.18). · **plain** A long memory means a tall spectrum at low frequency.
- **Result:** $S_e(\omega)=\frac1\pi\int_0^\infty R_{11}\cos\omega\tau\,d\tau$ is real and even; $\overline{u_1^2}=\int_{-\infty}^\infty S_e\,d\omega$ (12.22); $S_e(0)=\overline{u_1^2}\Lambda_t/\pi$ — "area =
  variance, height at zero = variance × memory / π".
- **Check:** units — $R_{11}$ [m²/s²] × dτ [s] ⇒ $S_e$ [m²/s]; × dω [1/s] ⇒ m²/s² ✓. Exponential pair: $S_e=\sigma^2\tau_c/[\pi(1+\omega^2\tau_c^2)]$,
  $\int S_e\,d\omega=\sigma^2$ ✓, $S_e(0)=\sigma^2\tau_c/\pi$ ✓ (σ = 1, τ_c = 0.5 s: 0.159). Code: Parseval to round-off.
- **What it means:** correlation and spectrum are one description in two forms; stretch one and the other squeezes. The
  two-sided convention halves every ordinate relative to a one-sided spectrum. **Fails when:** R does not decay (a pure
  wave): the spectrum is then a spike.
- **Traps:** putting the 2π in the wrong transform; using cyclic frequency with angular formulas; forgetting the factor 2
  when folding to one side; reading $S_e(0)$ as "energy at zero frequency" rather than a density.

### D06 · Mean continuity (12.27), $\partial u_i/\partial x_i=0$ (12.28) and the Reynolds-averaged momentum equation (12.30) — ★★, 13 steps, in C04 (notebook · `reynolds_stress_parcels`)
- **Goal:** find the equations obeyed by the mean flow, and see exactly where the fluctuations enter them. · **Start:** the
  Boussinesq set (R01): $\dfrac{\partial\tilde u_i}{\partial x_i}=0$ (4.10) and $\dfrac{\partial\tilde u_i}{\partial t}+\tilde u_j\dfrac{\partial\tilde u_i}{\partial x_j}=-\dfrac1{\rho_0}\dfrac{\partial\tilde p}{\partial x_i}-g\big[1-\alpha(\tilde T-T_0)\big]\delta_{i3}+\nu\dfrac{\partial^2\tilde u_i}{\partial x_j^2}$ (4.86). · **Plan:** • put
  the momentum equation in flux form • substitute mean + fluctuation • average with D01's rules • tidy the result into a
  stress. · **Tools:** D01's rules; Reynolds decomposition $\tilde u_i=U_i+u_i$ (12.24) with $\overline{\tilde u_i}=U_i$ (12.25) and $\overline{u_i}=0$ (12.26); product rule (P38); summation convention (ch02).
  · **Assumptions:** Boussinesq (constant ρ₀, ν, α); ensemble average; $\bar T$ is potential temperature in the buoyancy term.
- **Steps:**
  1. **did** Write the advection term as a divergence · **tex** $\tilde u_j\dfrac{\partial\tilde u_i}{\partial x_j}=\dfrac{\partial}{\partial x_j}(\tilde u_j\tilde u_i)-\tilde u_i\dfrac{\partial\tilde u_j}{\partial x_j}=\dfrac{\partial}{\partial x_j}(\tilde u_j\tilde u_i)$ · **why** Product rule, then
     continuity $\partial\tilde u_j/\partial x_j=0$ (4.10). A derivative of a product is easier to average than a product with a derivative. ·
     **plain** Advection is the divergence of a momentum flux.
  2. **did** Decompose the continuity equation · **tex** $\dfrac{\partial}{\partial x_i}(U_i+u_i)=0$ · **why** Substitute $\tilde u_i=U_i+u_i$ (12.24). · **plain** Mean plus fluctuation
     is divergence-free.
  3. **did** Average it · **tex** $\dfrac{\partial}{\partial x_i}\big(U_i+\overline{u_i}\big)=\dfrac{\partial U_i}{\partial x_i}=0$ · **why** The average passes through the derivative (D01's $\overline{\partial u^m/\partial x_j}=\partial\overline{u^m}/\partial x_j$ (12.8)) and
     $\overline{u_i}=0$ (12.26). This is (12.27). · **plain** The mean flow is divergence-free.
  4. **did** Subtract step 3 from step 2 · **tex** $\dfrac{\partial u_i}{\partial x_i}=0$ · **why** Two true equations may be subtracted. This is (12.28); we need it
     in D10. · **plain** The fluctuation is divergence-free too.
  5. **did** Substitute the decomposition into the flux-form momentum equation · **tex** $\begin{aligned}&\dfrac{\partial(U_i+u_i)}{\partial t}+\dfrac{\partial}{\partial x_j}\big((U_j+u_j)(U_i+u_i)\big)=-\dfrac1{\rho_0}\dfrac{\partial(P+p)}{\partial x_i}\\&\quad-g\big[1-\alpha(\bar T+T'-T_0)\big]\delta_{i3}+\nu\dfrac{\partial^2(U_i+u_i)}{\partial x_j^2}\end{aligned}$ ·
     **why** $\tilde u_i=U_i+u_i$, $\tilde p=P+p$, $\tilde T=\bar T+T'$ (12.24) in every field. This is (12.29); nothing has been approximated. · **plain** The exact equation, written for
     mean + fluctuation.
  6. **did** Expand the product · **tex** $(U_j+u_j)(U_i+u_i)=U_iU_j+U_iu_j+u_iU_j+u_iu_j$ · **why** Plain algebra; four terms: mean × mean, two mixed, fluctuation
     × fluctuation. · **plain** The momentum flux has four pieces.
  7. **did** Average the product · **tex** $\overline{(U_j+u_j)(U_i+u_i)}=U_iU_j+U_i\overline{u_j}+\overline{u_i}U_j+\overline{u_iu_j}=U_iU_j+\overline{u_iu_j}$ · **why** D01 step 7: a mean factors
     out, a single fluctuation averages to zero, a product of two fluctuations does not. · **plain** The mixed terms vanish;
     the covariance survives.
  8. **did** Average the linear terms · **tex** $\overline{\dfrac{\partial(U_i+u_i)}{\partial t}}=\dfrac{\partial U_i}{\partial t},\ \ \overline{\dfrac{\partial(P+p)}{\partial x_i}}=\dfrac{\partial P}{\partial x_i},\ \ \overline{T'}=0,\ \ \overline{\dfrac{\partial^2(U_i+u_i)}{\partial x_j^2}}=\dfrac{\partial^2U_i}{\partial x_j^2}$ · **why** Each is linear in the
     fields, so D01's rules $\overline{u+v}=\bar u+\bar v$ (12.4), $\overline{Au}=A\bar u$ (12.5), $\overline{\partial u/\partial t}=\partial\bar u/\partial t$ (12.6), $\overline{\partial u/\partial x_j}=\partial\bar u/\partial x_j$ (12.8) and $\overline{u_i}=0$ (12.26) leave the mean part only. · **plain** Every linear term keeps its form with
     the mean in place of the total.
  9. **did** Collect the averaged equation · **tex** $\dfrac{\partial U_i}{\partial t}+\dfrac{\partial}{\partial x_j}(U_iU_j)+\dfrac{\partial\overline{u_iu_j}}{\partial x_j}=-\dfrac1{\rho_0}\dfrac{\partial P}{\partial x_i}-g\big[1-\alpha(\bar T-T_0)\big]\delta_{i3}+\nu\dfrac{\partial^2U_i}{\partial x_j^2}$ · **why** Steps 7 and 8
     inside step 5. One term has no counterpart in the original equation. · **plain** The mean obeys the same equation plus
     the divergence of $\overline{u_iu_j}$.
  10. **did** Return the mean advection to its usual form · **tex** $\dfrac{\partial}{\partial x_j}(U_iU_j)=U_j\dfrac{\partial U_i}{\partial x_j}+U_i\dfrac{\partial U_j}{\partial x_j}=U_j\dfrac{\partial U_i}{\partial x_j}$ · **why** Product rule and the
      mean continuity equation $\partial U_j/\partial x_j=0$ (12.27) — step 1 run backwards for the mean. · **plain** The mean flow advects its
      own momentum.
  11. **did** Write the viscous term as a divergence of a stress · **tex** $\nu\dfrac{\partial^2U_i}{\partial x_j^2}=\dfrac1{\rho_0}\dfrac{\partial}{\partial x_j}\big(2\mu\bar S_{ij}\big),\quad\bar S_{ij}=\tfrac12\Big(\dfrac{\partial U_i}{\partial x_j}+\dfrac{\partial U_j}{\partial x_i}\Big)$ · **why**
      $\partial(2\bar S_{ij})/\partial x_j=\partial^2U_i/\partial x_j^2+\partial(\partial U_j/\partial x_j)/\partial x_i$ and the last term is zero by $\partial U_j/\partial x_j=0$ (12.27); μ = ρ₀ν (the mean-flow form of Ch. 4's $\mu\,\partial^2u_j/\partial x_i^2=2\mu\,\partial S_{ij}/\partial x_i$ (4.40)). ·
      **plain** Viscous friction of the mean flow is a stress, as in a laminar flow.
  12. **did** Write the pressure gradient the same way and move the new term to the right · **tex** $-\dfrac1{\rho_0}\dfrac{\partial P}{\partial x_i}-\dfrac{\partial\overline{u_iu_j}}{\partial x_j}=\dfrac1{\rho_0}\dfrac{\partial}{\partial x_j}\big(-P\delta_{ij}-\rho_0\overline{u_iu_j}\big)$
      · **why** $\partial(P\delta_{ij})/\partial x_j=\partial P/\partial x_i$; moving a term across the equals sign changes its sign; ρ₀ is constant. ·
      **plain** The covariance sits beside pressure and viscous stress — it acts as a stress.
  13. **did** Assemble · **tex** $\dfrac{\partial U_i}{\partial t}+U_j\dfrac{\partial U_i}{\partial x_j}=-g\big[1-\alpha(\bar T-T_0)\big]\delta_{i3}+\dfrac1{\rho_0}\dfrac{\partial}{\partial x_j}\big(-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}\big)$ · **why** Steps 9–12 together.
      This is (12.30) with $\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}$. · **plain** The mean flow feels pressure, viscous stress and a third
      stress made by the fluctuations.
- **Result:** $\dfrac{\partial U_i}{\partial x_i}=0$ (12.27); $\dfrac{\partial U_i}{\partial t}+U_j\dfrac{\partial U_i}{\partial x_j}=-g[1-\alpha(\bar T-T_0)]\delta_{i3}+\dfrac1{\rho_0}\dfrac{\partial\bar\tau_{ij}}{\partial x_j}$ (12.30) — "the mean flow obeys the
  Navier–Stokes equations with one extra stress, $-\rho_0\overline{u_iu_j}$".
- **Check:** units — $\rho_0\overline{u_iu_j}$ is kg/m³ × m²/s² = Pa ✓. Limit: no fluctuations ⇒ $\overline{u_iu_j}=0$, the stress is $\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}$ (12.30) and the Result is the laminar equation ✓. Number:
  $\overline{uv}=-0.1$ m²/s² in air ⇒ 0.12 Pa (C04 worked example). Code: `ch12.rans_sympy()` reproduces steps 5–9 symbolically.
- **What it means:** four equations — $\partial U_i/\partial x_i=0$ (12.27) and the three components of $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-g[1-\alpha(\bar T-T_0)]\delta_{i3}+\frac1{\rho_0}\frac{\partial\bar\tau_{ij}}{\partial x_j}$ (12.30) — for ten unknowns ($U_i$, P and six $\overline{u_iu_j}$): the closure problem. The
  rest of the chapter measures, scales or models this one term. **Fails when:** the density varies strongly
  (compressible flows need mass-weighted averages).
- **Traps:** averaging $\tilde u_j\,\partial\tilde u_i/\partial x_j$ without first writing it as a divergence (possible, but messier); thinking
  $\overline{U_iu_j}\neq0$; losing the minus sign when the new term crosses to the right; reading $\bar T$ as the thermometer
  temperature.

### D07 · The isotropic two-point tensor (12.40), its incompressible form (12.41), $g=f+\tfrac r2f'$, $\Lambda_g=\Lambda_f/2$, $\lambda_g=\lambda_f/\sqrt2$ — ★★★, 14 steps, in C05 (notebook)
- **Goal:** show that in isotropic, incompressible turbulence the whole two-point correlation tensor — nine functions of a
  vector — is fixed by one function f(r). · **Start:** $R_{ij}(\mathbf r)\equiv\overline{u_i(\mathbf x)u_j(\mathbf x+\mathbf r)}$ (12.23) with the definitions
  $f(r)\equiv\overline{u_\parallel(\mathbf x+\mathbf r)u_\parallel(\mathbf x)}/\overline{u^2}$, $g(r)\equiv\overline{u_\perp(\mathbf x+\mathbf r)u_\perp(\mathbf x)}/\overline{u^2}$ (12.38). · **Plan:** • write the most general tensor that isotropy
  allows • read its two functions off by pointing **r** along one axis • impose incompressibility • integrate and
  differentiate the result. · **Tools:** isotropic tensors (P119); derivatives of functions of r (P297, primer above);
  $\partial u_i/\partial x_i=0$ (12.28) from D06; integration by parts (P218a); Taylor series (P26); sympy (P40). · **Assumptions:** homogeneous, isotropic
  (reflection-invariant), incompressible, three-dimensional; f decays faster than 1/r.
- **Steps:**
  1. **did** Note what R can depend on · **tex** $R_{ij}=R_{ij}(\mathbf r)$ · **why** Homogeneity: the statistics are the same at every **x**, so only the separation **r** remains: $R_{ij}(\mathbf r)\equiv\overline{u_i(\mathbf x)u_j(\mathbf x+\mathbf r)}$ (12.23). · **plain** Only "how far apart and in which direction" matters.
  2. **did** Write the most general isotropic form · **tex** $R_{ij}=F_R(r)\,r_ir_j+G_R(r)\,\delta_{ij}$ · **why** With no preferred direction, a second-order tensor
     can be built only from $\delta_{ij}$ and the vector **r** itself, with scalar coefficients depending on r = ∣**r**∣; an
     $\epsilon_{ijk}r_k$ term would change sign under reflection. This is (12.40). · **plain** Two unknown functions instead of nine.
  3. **did** Point **r** along the 1-axis and take the 11-component · **tex** $\mathbf r=(r,0,0):\quad R_{11}=F_Rr^2+G_R$ · **why** We are free to choose the axes
     (isotropy); then $r_1r_1=r^2$ and $\delta_{11}=1$. · **plain** The component along the separation.
  4. **did** Recognise the longitudinal correlation · **tex** $R_{11}=\overline{u_1(\mathbf x)u_1(\mathbf x+r\mathbf e_1)}=\overline{u^2}\,f(r)$ · **why** Both velocity components lie along **r**: that is the definition $f(r)\equiv\overline{u_\parallel(\mathbf x+\mathbf r)u_\parallel(\mathbf x)}/\overline{u_\parallel^2}$ (12.38). · **plain** Along the line, the tensor is f.
  5. **did** Take the 22-component · **tex** $R_{22}=G_R=\overline{u^2}\,g(r)$ · **why** $r_2=0$, $\delta_{22}=1$; both components are perpendicular to **r**: the
     definition of g. · **plain** Across the line, the tensor is g.
  6. **did** Solve for the two functions · **tex** $G_R=\overline{u^2}\,g,\qquad F_R=\overline{u^2}\,\dfrac{f-g}{r^2}$ · **why** Subtract step 5 from steps 3–4 and divide by r². These are
     the F and G the book quotes after $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40). · **plain** The tensor is known once f and g are.
  7. **did** Impose incompressibility at the second point · **tex** $\dfrac{\partial R_{ij}}{\partial r_j}=\overline{u_i(\mathbf x)\,\dfrac{\partial u_j}{\partial x_j}(\mathbf x+\mathbf r)}=0$ · **why** Differentiating with respect to **r**
     acts only on the second factor (D01's $\overline{\partial u^m/\partial x_j}=\partial\overline{u^m}/\partial x_j$ (12.8)), and $\partial u_j/\partial x_j=0$ (12.28). · **plain** The tensor is divergence-free in **r**.
  8. **did** Differentiate the first part · **tex** $\dfrac{\partial}{\partial r_j}\big(F_Rr_ir_j\big)=F_R'\dfrac{r_j}{r}r_ir_j+F_R\big(\delta_{ij}r_j+r_i\delta_{jj}\big)=\big(rF_R'+4F_R\big)r_i$ · **why** Product rule with
     $\partial r/\partial r_j=r_j/r$, $\partial r_i/\partial r_j=\delta_{ij}$, $r_jr_j=r^2$, $\delta_{jj}=3$ (P297). · **plain** The divergence of the $r_ir_j$ part points along **r**.
  9. **did** Differentiate the second part · **tex** $\dfrac{\partial}{\partial r_j}\big(G_R\delta_{ij}\big)=G_R'\dfrac{r_i}{r}$ · **why** Chain rule; $\delta_{ij}$ picks out j = i. · **plain** So does the
     divergence of the $\delta_{ij}$ part.
  10. **did** Add them and remove the common factor · **tex** $rF_R'+4F_R+\dfrac{G_R'}{r}=0$ · **why** Step 7 must hold for every **r**, so the coefficient of
      $r_i$ vanishes. · **plain** One ordinary differential equation links F and G.
  11. **did** Substitute step 6 · **tex** $\dfrac{f'}{r}+\dfrac{2(f-g)}{r^2}=0\ \Rightarrow\ g=f+\dfrac r2\dfrac{df}{dr}$ · **why** $rF_R'=\overline{u^2}\big[(f'-g')/r-2(f-g)/r^2\big]$; the g′ terms cancel and
      $-2+4=2$. · **plain** The transverse correlation is fixed by the longitudinal one.
  12. **did** Put g back into the tensor · **tex** $R_{ij}=\overline{u^2}\Big\{f\,\delta_{ij}+\dfrac r2\dfrac{df}{dr}\Big(\delta_{ij}-\dfrac{r_ir_j}{r^2}\Big)\Big\}$ · **why** $F_R=\overline{u^2}(f-g)/r^2=-\overline{u^2}f'/(2r)$ and $G_R=\overline{u^2}(f+rf'/2)$ in
      step 2. This is (12.41). · **plain** One function f(r) gives all nine components.
  13. **did** Integrate g over all r · **tex** $\Lambda_g=\displaystyle\int_0^\infty\Big(f+\dfrac r2f'\Big)dr=\Lambda_f+\Big[\dfrac{rf}{2}\Big]_0^\infty-\dfrac12\int_0^\infty f\,dr=\dfrac{\Lambda_f}{2}$ · **why** Integration by parts on the second
      term (P218a); the boundary term vanishes because f decays faster than 1/r. Definitions (12.39). · **plain** The
      transverse integral scale is half the longitudinal one.
  14. **did** Differentiate g twice at r = 0 · **tex** $g''=2f''+\dfrac r2f'''\ \Rightarrow\ g''(0)=2f''(0)\ \Rightarrow\ \lambda_g=\dfrac{\lambda_f}{\sqrt2}$ · **why** $g'=\tfrac32f'+\tfrac r2f''$, then once more;
      $\lambda^2\equiv-2/[\text{second derivative at }0]$ (12.39). · **plain** The transverse correlation is twice as sharply curved at its
      top.
- **Result:** $R_{ij}=\overline{u^2}\{f\delta_{ij}+\frac r2f'(\delta_{ij}-r_ir_j/r^2)\}$ (12.41); $g=f+\frac r2f'$; $\Lambda_g=\Lambda_f/2$; $\lambda_g=\lambda_f/\sqrt2$ — "isotropy and
  incompressibility leave one free function".
- **Check:** units — $R_{ij}$ in m²/s². r → 0: f = g = 1, $R_{ij}=\overline{u^2}\delta_{ij}$ ✓ (equal normal stresses, no shear stress). Trace:
  $R_{ii}=\overline{u^2}(3f+rf')$, at r = 0 it is $3\overline{u^2}=2\bar e$ ✓. Gaussian $f=e^{-r^2/L^2}$: $g=(1-r^2/L^2)e^{-r^2/L^2}$ — negative beyond r = L,
  $\Lambda_g=\tfrac{\sqrt\pi}4L=\Lambda_f/2$ ✓.
- **sympy check intent:** rebuild the tensor of step 12 for a symbolic f, take its divergence, and get zero; redo steps 13–14
  for a Gaussian. `check_src` sketch:
  ```python
  import sympy as sp                                             # symbolic algebra
  r1, r2, r3, u2 = sp.symbols("r1 r2 r3 u2", positive=True)     # separation components, one-component variance
  rv = sp.Matrix([r1, r2, r3]); r = sp.sqrt(r1**2 + r2**2 + r3**2)   # the vector r and its length
  f = sp.Function("f")                                           # the longitudinal correlation, left general
  g = f(r) + r/2*sp.diff(f(r), r)                                # step 11: g = f + (r/2) f'
  R = sp.Matrix(3, 3, lambda i, j: u2*((f(r) - g)/r**2*rv[i]*rv[j] + g*sp.KroneckerDelta(i, j)))  # steps 2, 6
  div = [sp.simplify(sum(sp.diff(R[i, j], rv[j]) for j in range(3))) for i in range(3)]          # step 7
  assert div == [0, 0, 0]                                        # divergence-free for ANY f
  s, L = sp.symbols("s L", positive=True); fG = sp.exp(-s**2/L**2)   # a Gaussian f to test steps 13 and 14
  gG = fG + s/2*sp.diff(fG, s)                                   # its transverse partner
  assert sp.simplify(sp.integrate(gG, (s, 0, sp.oo)) - sp.integrate(fG, (s, 0, sp.oo))/2) == 0   # Lambda_g = Lambda_f/2
  assert sp.simplify(sp.diff(gG, s, 2).subs(s, 0) - 2*sp.diff(fG, s, 2).subs(s, 0)) == 0         # g''(0) = 2 f''(0)
  ```
- **What it means:** measure f with one probe towed along a line and you know every two-point correlation; g must go
  negative somewhere (fluid crossing a line has to come back). **Fails when:** the turbulence is anisotropic (shear, walls,
  stratification) or two-dimensional (there $g=d(rf)/dr$ and the factors change).
- **Traps:** using $\delta_{jj}=1$ (it is 3); forgetting $\partial r/\partial r_j=r_j/r$; quoting the scale definitions $\Lambda_f\equiv\int_0^\infty f(r)\,dr$, $\lambda_f^2\equiv-2/[d^2f/dr^2]_{r=0}$ (12.39) for the tensor as Exercise 12.18 does
  (slip #14: the tensor is $R_{ij}=F(r)r_ir_j+G(r)\delta_{ij}$ (12.40)); confusing these F, G with the jet profiles of $U=U_{CL}(x)F(y/\delta(x))$ (12.56) and $-\overline{uv}=\Psi(x)G(y/\delta(x))$ (12.57).

### D08 · The isotropic dissipation (12.42) → (12.43): moments 2 : 4 : −1 and $\bar\varepsilon=30\nu\overline{u^2}/\lambda_f^2=15\nu\overline{u^2}/\lambda_g^2=15\nu\overline{(\partial u_1/\partial x_1)^2}$ — ★★★, 13 steps, in C05 (notebook)
- **Goal:** reduce the dissipation rate — an average over nine velocity gradients — to one measurable number. · **Start:**
  $\bar\varepsilon=\dfrac\nu2\,\overline{\Big(\dfrac{\partial u_i}{\partial x_j}+\dfrac{\partial u_j}{\partial x_i}\Big)^2}$ (12.42). · **Plan:** • expand the square and count equal terms using isotropy • express gradient
  moments through the correlation tensor • expand f near r = 0 • read off the three moments and add. · **Tools:** summation
  convention (ch02); the isotropy statements $\overline{(\partial u_1/\partial x_1)^n}=\overline{(\partial u_2/\partial x_2)^n}=\overline{(\partial u_3/\partial x_3)^n}$ (12.36) and $\overline{(\partial u_1/\partial x_2)^n}=\overline{(\partial u_1/\partial x_3)^n}=\dots=\overline{(\partial u_3/\partial x_2)^n}$ (12.37); D07's $R_{ij}=\overline{u^2}\{f\delta_{ij}+\frac r2f'(\delta_{ij}-r_ir_j/r^2)\}$ (12.41); Taylor expansion (P26); P297; sympy (P40). ·
  **Assumptions:** homogeneous, isotropic, incompressible, three-dimensional.
- **Steps:**
  1. **did** Expand the square · **tex** $\Big(\dfrac{\partial u_i}{\partial x_j}+\dfrac{\partial u_j}{\partial x_i}\Big)^2=2\dfrac{\partial u_i}{\partial x_j}\dfrac{\partial u_i}{\partial x_j}+2\dfrac{\partial u_i}{\partial x_j}\dfrac{\partial u_j}{\partial x_i}$ · **why** $(a+b)^2=a^2+2ab+b^2$; summed over i and j, the two squares
     are the same sum (rename i ↔ j). · **plain** Two kinds of double sum.
  2. **did** Average · **tex** $\bar\varepsilon=\nu\Big[\,\overline{\dfrac{\partial u_i}{\partial x_j}\dfrac{\partial u_i}{\partial x_j}}+\overline{\dfrac{\partial u_i}{\partial x_j}\dfrac{\partial u_j}{\partial x_i}}\,\Big]$ · **why** D01's linearity; the factor ½ of $\bar\varepsilon=\frac\nu2\overline{(\partial u_i/\partial x_j+\partial u_j/\partial x_i)^2}$ (12.42) cancels the 2. · **plain** The
     dissipation is ν times two sums of nine gradient moments each.
  3. **did** Count the first sum · **tex** $\overline{\dfrac{\partial u_i}{\partial x_j}\dfrac{\partial u_i}{\partial x_j}}=3a+6b,\quad a\equiv\overline{\Big(\dfrac{\partial u_1}{\partial x_1}\Big)^2},\ b\equiv\overline{\Big(\dfrac{\partial u_1}{\partial x_2}\Big)^2}$ · **why** Three terms have i = j; by (12.36),
     $\overline{(\partial u_1/\partial x_1)^n}=\overline{(\partial u_2/\partial x_2)^n}=\overline{(\partial u_3/\partial x_3)^n}$, they are equal. Six have i ≠ j and are equal by (12.37), $\overline{(\partial u_1/\partial x_2)^n}=\overline{(\partial u_1/\partial x_3)^n}=\overline{(\partial u_2/\partial x_1)^n}=\overline{(\partial u_2/\partial x_3)^n}=\overline{(\partial u_3/\partial x_1)^n}=\overline{(\partial u_3/\partial x_2)^n}$. · **plain** Three
     "along" moments and six "across" moments.
  4. **did** Count the second sum · **tex** $\overline{\dfrac{\partial u_i}{\partial x_j}\dfrac{\partial u_j}{\partial x_i}}=3a+6c,\quad c\equiv\overline{\dfrac{\partial u_1}{\partial x_2}\dfrac{\partial u_2}{\partial x_1}}$ · **why** For i = j the product is again a square (a); the
     six i ≠ j terms are all equal to c by isotropy. · **plain** The cross sum brings a third kind of moment.
  5. **did** Add · **tex** $\bar\varepsilon=\nu(6a+6b+6c)=6\nu\Big\{\overline{\Big(\dfrac{\partial u_1}{\partial x_1}\Big)^2}+\overline{\Big(\dfrac{\partial u_1}{\partial x_2}\Big)^2}+\overline{\dfrac{\partial u_1}{\partial x_2}\dfrac{\partial u_2}{\partial x_1}}\Big\}$ · **why** Steps 2–4. This is the first equality of
     (12.43). · **plain** Three different moments — not three equal ones.
  6. **did** Differentiate the correlation once in **r** · **tex** $\dfrac{\partial R_{ij}}{\partial r_l}=\overline{u_i(\mathbf x)\,\dfrac{\partial u_j}{\partial x_l}(\mathbf x+\mathbf r)}$ · **why** **r** appears only in the second factor;
     the derivative passes through the average (D01's $\overline{\partial u^m/\partial x_j}=\partial\overline{u^m}/\partial x_j$ (12.8)). · **plain** One derivative lands on the second velocity.
  7. **did** Shift the origin, differentiate again, set r = 0 · **tex** $\dfrac{\partial^2R_{ij}}{\partial r_k\partial r_l}\Big|_{r=0}=-\overline{\dfrac{\partial u_i}{\partial x_k}\dfrac{\partial u_j}{\partial x_l}}$ · **why** By homogeneity
     $\overline{u_i(\mathbf x)\partial_lu_j(\mathbf x+\mathbf r)}=\overline{u_i(\mathbf x'-\mathbf r)\partial_lu_j(\mathbf x')}$; now **r** sits in the first factor with a minus sign. · **plain**
     Gradient moments are (minus) the curvature of the correlation at zero separation.
  8. **did** Expand f near r = 0 · **tex** $f\approx1-\dfrac{r^2}{\lambda_f^2},\qquad\dfrac r2\dfrac{df}{dr}\approx-\dfrac{r^2}{\lambda_f^2}$ · **why** Taylor series of an even function with
     $\lambda_f^2\equiv-2/f''(0)$ (12.39), exactly as in D04; ≈ because only second derivatives at 0 are needed. · **plain** Near zero
     separation f is a parabola.
  9. **did** Insert into the isotropic tensor · **tex** $R_{ij}\approx\overline{u^2}\Big[\delta_{ij}-\dfrac{2r^2}{\lambda_f^2}\delta_{ij}+\dfrac{r_ir_j}{\lambda_f^2}\Big]$ · **why** (12.41),
     $R_{ij}=\overline{u^2}\{f\delta_{ij}+\frac r2f'(\delta_{ij}-r_ir_j/r^2)\}$, with step 8: $-r^2/\lambda_f^2-r^2/\lambda_f^2=-2r^2/\lambda_f^2$. · **plain** Near the origin the tensor is
     quadratic in **r**.
  10. **did** Differentiate twice · **tex** $-\dfrac{\partial^2R_{ij}}{\partial r_k\partial r_l}=\dfrac{\overline{u^2}}{\lambda_f^2}\big[4\delta_{ij}\delta_{kl}-\delta_{ik}\delta_{jl}-\delta_{il}\delta_{jk}\big]$ · **why** $\partial^2(r^2)/\partial r_k\partial r_l=2\delta_{kl}$ and
      $\partial^2(r_ir_j)/\partial r_k\partial r_l=\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}$ (P297). · **plain** Every gradient moment, in one formula.
  11. **did** Read off the three moments · **tex** $a=\dfrac{2\overline{u^2}}{\lambda_f^2},\qquad b=\dfrac{4\overline{u^2}}{\lambda_f^2},\qquad c=-\dfrac{\overline{u^2}}{\lambda_f^2}$ · **why** a: i = j = k = l = 1 gives 4 − 1 − 1; b: i =
      j = 1, k = l = 2 gives 4; c: i = 1, k = 2, j = 2, l = 1 gives −1. · **plain** The moments are in the ratio 2 : 4 : −1.
  12. **did** Add them · **tex** $\bar\varepsilon=6\nu\,(2+4-1)\dfrac{\overline{u^2}}{\lambda_f^2}=30\nu\dfrac{\overline{u^2}}{\lambda_f^2}=-15\nu\,\overline{u^2}\Big[\dfrac{d^2f}{dr^2}\Big]_{r=0}$ · **why** Step 5 with step 11; $1/\lambda_f^2=-f''(0)/2$. The
      middle equalities of (12.43). · **plain** The dissipation is fixed by the curvature of f at the origin.
  13. **did** Rewrite with the transverse microscale and with one gradient · **tex** $\bar\varepsilon=15\nu\dfrac{\overline{u^2}}{\lambda_g^2}=15\nu\,\overline{\Big(\dfrac{\partial u_1}{\partial x_1}\Big)^2}$ · **why** $\lambda_g^2=\lambda_f^2/2$ (D07
      step 14); and $a=2\overline{u^2}/\lambda_f^2$ from step 11. The first form is the last member of (12.43); the one-gradient form $15\nu\overline{(\partial u_1/\partial x_1)^2}$ is ours (a corollary of step 11, not part of the book's numbered equation). · **plain** One velocity gradient, measured along a line, gives the whole
      dissipation.
- **Result:** $\bar\varepsilon=6\nu\{a+b+c\}=-15\nu\overline{u^2}f''(0)=30\nu\overline{u^2}/\lambda_f^2=15\nu\overline{u^2}/\lambda_g^2$ (12.43) — "symmetry turns nine gradients into
  one".
- **Check:** units — ν [m²/s] × [m²/s²]/[m²] = m²/s³ ✓. Consistency: in homogeneous incompressible turbulence the cross sum
  $\overline{\partial_ju_i\,\partial_iu_j}=3a+6c$ must vanish (it is a double divergence): 3 × 2 + 6 × (−1) = 0 ✓ — so also $\bar\varepsilon=\nu\overline{(\partial_ju_i)^2}$ = ν(3a +
  6b) = 30ν$\overline{u^2}/\lambda_f^2$ ✓. Number: ν = 1.5 × 10⁻⁵, $\overline{u^2}$ = 1, λ_f = 0.01 m ⇒ 4.5 m²/s³.
- **sympy check intent:** build the quadratic tensor of step 9, differentiate it as in step 10, and read the three moments
  and the total. `check_src` sketch:
  ```python
  import sympy as sp                                           # symbolic algebra
  r1, r2, r3, u2, lam, nu = sp.symbols("r1 r2 r3 u2 lambda_f nu", positive=True)   # separation components, one-component variance, λ_f, viscosity
  rv = [r1, r2, r3]; rr = r1**2 + r2**2 + r3**2                # separation and its square
  d = sp.KroneckerDelta                                        # delta_ij
  R = lambda i, j: u2*(d(i, j) - 2*rr/lam**2*d(i, j) + rv[i]*rv[j]/lam**2)   # step 9
  M = lambda i, k, j, l: -sp.diff(R(i, j), rv[k], rv[l])       # step 7: <d_k u_i d_l u_j> = -d2R_ij/dr_k dr_l
  a, b, c = M(0, 0, 0, 0), M(0, 1, 0, 1), M(0, 1, 1, 0)        # the three moments of step 3 and 4
  assert [sp.simplify(x*lam**2/u2) for x in (a, b, c)] == [2, 4, -1]          # step 11
  eps = nu*sum(M(i, j, i, j) + M(i, j, j, i) for i in range(3) for j in range(3))   # step 2, all 18 terms
  assert sp.simplify(eps - 30*nu*u2/lam**2) == 0               # step 12
  assert sp.simplify(sum(M(i, j, j, i) for i in range(3) for j in range(3))) == 0   # the cross sum vanishes
  ```
- **What it means:** oceanographers measure $\overline{(\partial u_1/\partial x_1)^2}$ with a shear probe and multiply by 15ν (or 7.5ν for a transverse
  gradient) to get ε̄. **Fails when:** the small scales are not isotropic (low Reynolds number, strong stratification), or
  in two dimensions (the factor is not 15).
- **Traps:** assuming the three moments in the bracket of $\bar\varepsilon=6\nu\big\{\overline{(\partial u_1/\partial x_1)^2}+\overline{(\partial u_1/\partial x_2)^2}+\overline{(\partial u_1/\partial x_2)(\partial u_2/\partial x_1)}\big\}$ (12.43) are equal; losing the minus sign in step 7; mixing λ_f and λ_g
  (a factor 2 in ε̄); using the total variance $2\bar e$ where $\overline{u^2}$ is one component.

### D09 · The kinetic-energy budget of the mean flow (12.46) — ★★, 9 steps, in C06 (notebook · `turbulent_energy_budget`)
- **Goal:** find where the mean flow's kinetic energy $\bar E=\tfrac12U_i^2$ goes — in particular, the term that hands energy to the
  turbulence. · **Start:** $\dfrac{\partial U_i}{\partial t}+U_j\dfrac{\partial U_i}{\partial x_j}=-g\big[1-\alpha(\bar T-T_0)\big]\delta_{i3}+\dfrac1{\rho_0}\dfrac{\partial\bar\tau_{ij}}{\partial x_j}$, $\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}$ (12.30). · **Plan:** •
  multiply by $U_i$ • turn each product into a derivative of $\bar E$ or a divergence • what cannot be made a divergence is a
  source or sink. · **Tools:** product rule (P38); mean continuity $\partial U_i/\partial x_i=0$ (12.27); contraction of a symmetric tensor with a velocity
  gradient (ch02 §2.10). · **Assumptions:** Boussinesq; constant ν.
- **Steps:**
  1. **did** Multiply the mean momentum equation by $U_i$ (sum over i) · **tex** $U_i\dfrac{\partial U_i}{\partial t}+U_iU_j\dfrac{\partial U_i}{\partial x_j}=-gU_3\big[1-\alpha(\bar T-T_0)\big]+\dfrac{U_i}{\rho_0}\dfrac{\partial\bar\tau_{ij}}{\partial x_j}$ · **why** Force ×
     velocity = power: multiplying a momentum equation by the velocity gives an energy equation (the move of Ch. 4's
     mechanical-energy equation). · **plain** The rate of working of each force on the mean flow.
  2. **did** Recognise the time derivative · **tex** $U_i\dfrac{\partial U_i}{\partial t}=\dfrac{\partial}{\partial t}\Big(\tfrac12U_i^2\Big)=\dfrac{\partial\bar E}{\partial t}$ · **why** Chain rule backwards: $\partial(U_i^2)/\partial t=2U_i\,\partial U_i/\partial t$. ·
     **plain** The first term is the rate of change of mean kinetic energy.
  3. **did** Do the same for advection · **tex** $U_iU_j\dfrac{\partial U_i}{\partial x_j}=U_j\dfrac{\partial\bar E}{\partial x_j}$ · **why** Same chain rule with $x_j$ for t. · **plain** Mean energy carried
     along by the mean flow.
  4. **did** Rewrite the gravity term with the mean density · **tex** $-gU_3\big[1-\alpha(\bar T-T_0)\big]=-\dfrac{g}{\rho_0}\,\bar\rho\,U_3$ · **why** Boussinesq equation of state,
     $\bar\rho=\rho_0[1-\alpha(\bar T-T_0)]$. · **plain** Rising mean motion loses kinetic energy to potential energy.
  5. **did** Split the stress term with the product rule · **tex** $\dfrac{U_i}{\rho_0}\dfrac{\partial\bar\tau_{ij}}{\partial x_j}=\dfrac1{\rho_0}\dfrac{\partial}{\partial x_j}\big(U_i\bar\tau_{ij}\big)-\dfrac{\bar\tau_{ij}}{\rho_0}\dfrac{\partial U_i}{\partial x_j}$ · **why** $\partial(U_i\bar\tau_{ij})/\partial x_j=U_i\,\partial\bar\tau_{ij}/\partial x_j+\bar\tau_{ij}\,\partial U_i/\partial x_j$.
     A divergence only moves energy; the remainder creates or destroys it. · **plain** Stress work = transport of energy −
     deformation work.
  6. **did** Write out the transport part · **tex** $\dfrac{U_i\bar\tau_{ij}}{\rho_0}=-\dfrac{U_jP}{\rho_0}+2\nu U_i\bar S_{ij}-\overline{u_iu_j}\,U_i$ · **why** Substitute $\bar\tau_{ij}$; $U_i\delta_{ij}=U_j$; μ/ρ₀ = ν. · **plain**
     Energy is carried by pressure, by viscous stress and by Reynolds stress.
  7. **did** Contract the pressure and viscous parts with the gradient · **tex** $\dfrac{-P\delta_{ij}+2\mu\bar S_{ij}}{\rho_0}\dfrac{\partial U_i}{\partial x_j}=-\dfrac{P}{\rho_0}\dfrac{\partial U_i}{\partial x_i}+2\nu\bar S_{ij}\dfrac{\partial U_i}{\partial x_j}=2\nu\bar S_{ij}\bar S_{ij}$ · **why** $\partial U_i/\partial x_i=0$
     (12.27); a symmetric tensor contracted with a gradient sees only its symmetric part (ch02 §2.10). · **plain** Pressure
     does no deformation work in an incompressible flow; viscosity always dissipates.
  8. **did** Contract the Reynolds part · **tex** $-\dfrac{(-\rho_0\overline{u_iu_j})}{\rho_0}\dfrac{\partial U_i}{\partial x_j}=+\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}$ · **why** The minus of step 5 meets the minus in
     $\bar\tau_{ij}$. For U(y): $\overline{uv}\,dU/dy$, negative when $dU/dy>0$ (C04). · **plain** The work of the mean flow against the Reynolds
     stress — a loss.
  9. **did** Assemble · **tex** $\begin{aligned}\dfrac{\partial\bar E}{\partial t}+U_j\dfrac{\partial\bar E}{\partial x_j}=&\ \dfrac{\partial}{\partial x_j}\Big(-\dfrac{U_jP}{\rho_0}+2\nu U_i\bar S_{ij}-\overline{u_iu_j}U_i\Big)\\&-2\nu\bar S_{ij}\bar S_{ij}+\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}-\dfrac{g}{\rho_0}\bar\rho U_3\end{aligned}$ · **why** Steps 2–8 in step 1. This is (12.46). ·
     **plain** Mean energy changes by transport, direct viscous loss, loss to turbulence and loss to potential energy.
- **Result:** $\frac{\partial\bar E}{\partial t}+U_j\frac{\partial\bar E}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac{U_jP}{\rho_0}+2\nu U_i\bar S_{ij}-\overline{u_iu_j}U_i\big)-2\nu\bar S_{ij}\bar S_{ij}+\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}-\frac g{\rho_0}\bar\rho U_3$ (12.46), as in step 9 — "the mean flow loses energy to the turbulence at the rate $-\overline{u_iu_j}\,\partial U_i/\partial x_j$".
- **Check:** units — every term m²/s³. No fluctuations ⇒ Ch. 4's mechanical-energy equation ✓. Size: direct dissipation /
  production ~ ν/(UL) = 1/Re where the mean gradient is ~ ΔU/L (free shear flows, the outer part of a wall flow) — not in the viscous wall layer, where the gradient is $u_*^2/\nu$ (N80). Code: `ch12.mean_energy_budget_sympy()` residual 0; integral identity in a channel: pressure
  work = direct dissipation + production (`channel_energy_budget`).
- **What it means:** at high Reynolds number the mean flow, away from walls, hardly feels viscosity directly; it is drained by the Reynolds
  stress, and D10 shows the same term arriving in the turbulence. **Fails when:** density varies beyond Boussinesq.
- **Traps:** calling $+\overline{u_iu_j}\,\partial U_i/\partial x_j$ a gain because of its plus sign (it is negative); forgetting the divergence terms
  vanish only when integrated over a closed or homogeneous region.

### D10 · The turbulent kinetic-energy budget (12.47) — ★★★, 15 steps, in C06 (notebook · `turbulent_energy_budget`)
- **Goal:** a budget for the turbulent kinetic energy $\bar e=\tfrac12\overline{u_i^2}$: what feeds it, what moves it, what destroys it. ·
  **Start:** the total (Boussinesq) momentum equation $\frac{\partial\tilde u_i}{\partial t}+\tilde u_j\frac{\partial\tilde u_i}{\partial x_j}=-\frac1{\rho_0}\frac{\partial\tilde p}{\partial x_i}-g[1-\alpha(\tilde T-T_0)]\delta_{i3}+\nu\frac{\partial^2\tilde u_i}{\partial x_j^2}$ (4.86) and the mean one, $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-g[1-\alpha(\bar T-T_0)]\delta_{i3}+\frac1{\rho_0}\frac{\partial\bar\tau_{ij}}{\partial x_j}$ with $\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}$ (12.30). · **Plan:** • subtract mean from total to get the
  equation of the fluctuation • multiply by $u_i$ and average • turn each term into a derivative of $\bar e$, a divergence, or a
  source • handle the viscous term last. · **Tools:** D01's rules; product rule on products of three factors (P38);
  $\partial u_i/\partial x_i=0$ (12.28); symmetric contraction (ch02 §2.10); the sympy averaging operator (P294). · **Assumptions:**
  Boussinesq; constant ν; $T'$ is the potential-temperature fluctuation.
- **Steps:**
  1. **did** Write the total equation with every field decomposed · **tex** $\begin{aligned}&\dfrac{\partial(U_i+u_i)}{\partial t}+(U_j+u_j)\dfrac{\partial(U_i+u_i)}{\partial x_j}=-\dfrac1{\rho_0}\dfrac{\partial(P+p)}{\partial x_i}\\&\quad-g\big[1-\alpha(\bar T+T'-T_0)\big]\delta_{i3}+\nu\dfrac{\partial^2(U_i+u_i)}{\partial x_j^2}\end{aligned}$ · **why** The Start's Boussinesq equation (4.86), $D\mathbf u/Dt=-\rho_0^{-1}\nabla p'+(\rho'/\rho_0)\mathbf g+\nu\nabla^2\mathbf u$ in Ch. 4, with $\tilde u_i=U_i+u_i$ (12.24); the advective form is more convenient here than the flux form of D06. · **plain** The exact equation for
     mean + fluctuation.
  2. **did** Write the mean equation in the same form · **tex** $\dfrac{\partial U_i}{\partial t}+U_j\dfrac{\partial U_i}{\partial x_j}=-\dfrac1{\rho_0}\dfrac{\partial P}{\partial x_i}-g\big[1-\alpha(\bar T-T_0)\big]\delta_{i3}+\nu\dfrac{\partial^2U_i}{\partial x_j^2}-\dfrac{\partial\overline{u_iu_j}}{\partial x_j}$ · **why** This is D06 step 9 with step 10 — the mean momentum equation before its stresses were grouped into $\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}$ (12.30). · **plain** The equation of the mean alone.
  3. **did** Subtract · **tex** $\dfrac{\partial u_i}{\partial t}+U_j\dfrac{\partial u_i}{\partial x_j}+u_j\dfrac{\partial U_i}{\partial x_j}+u_j\dfrac{\partial u_i}{\partial x_j}-\dfrac{\partial\overline{u_iu_j}}{\partial x_j}=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x_i}+g\alpha T'\delta_{i3}+\nu\dfrac{\partial^2u_i}{\partial x_j^2}$ · **why** Total minus mean leaves the
     fluctuation's own equation; $(U_j+u_j)\partial(U_i+u_i)/\partial x_j$ gives four terms, one of which is in the mean equation. ·
     **plain** How one fluctuation evolves: carried by the mean, distorted by the mean shear, acting on itself.
  4. **did** Multiply by $u_i$ and average · **tex** $\overline{u_i\dfrac{\partial u_i}{\partial t}}+U_j\overline{u_i\dfrac{\partial u_i}{\partial x_j}}+\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}+\overline{u_iu_j\dfrac{\partial u_i}{\partial x_j}}-\overline{u_i}\,\dfrac{\partial\overline{u_iu_j}}{\partial x_j}=-\dfrac1{\rho_0}\overline{u_i\dfrac{\partial p}{\partial x_i}}+g\alpha\overline{u_3T'}+\nu\overline{u_i\dfrac{\partial^2u_i}{\partial x_j^2}}$ · **why** Force × velocity again (D09 step 1), now for the fluctuation;
     the average is taken because single realizations are not reproducible. Each term is treated in turn below. ·
     **plain** The mean rate of working on the fluctuations.
  5. **did** Time derivative · **tex** $\overline{u_i\dfrac{\partial u_i}{\partial t}}=\dfrac{\partial}{\partial t}\Big(\tfrac12\overline{u_i^2}\Big)=\dfrac{\partial\bar e}{\partial t}$ · **why** Chain rule backwards, then D01's $\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t$ (12.6). · **plain** The rate of
     change of turbulent kinetic energy.
  6. **did** Advection by the mean flow · **tex** $U_j\,\overline{u_i\dfrac{\partial u_i}{\partial x_j}}=U_j\dfrac{\partial\bar e}{\partial x_j}$ · **why** $U_j$ is a mean and factors out; chain rule and $\overline{\partial u^m/\partial x_j}=\partial\overline{u^m}/\partial x_j$ (12.8). ·
     **plain** Turbulent energy carried along by the mean flow.
  7. **did** The mean-stress term · **tex** $\overline{u_i}\,\dfrac{\partial\overline{u_iu_j}}{\partial x_j}=0$ · **why** $\partial\overline{u_iu_j}/\partial x_j$ is already an average — a "constant" for the averaging
     — and $\overline{u_i}=0$ (12.26). · **plain** This term kept the fluctuation equation mean-free; it does no work on average.
  8. **did** The shear term · **tex** $\overline{u_iu_j}\,\dfrac{\partial U_i}{\partial x_j}\quad(\text{on the left})\ \longrightarrow\ -\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}\quad(\text{on the right})$ · **why** The mean gradient factors out
     of the average. Moved to the right it is exactly minus the term of D09 step 8. · **plain** <span style="color:#f97316">Shear production</span>: what the
     mean flow lost, the turbulence gains.
  9. **did** The triple term · **tex** $\overline{u_iu_j\dfrac{\partial u_i}{\partial x_j}}=\overline{u_j\dfrac{\partial}{\partial x_j}\Big(\tfrac12u_i^2\Big)}=\dfrac{\partial}{\partial x_j}\Big(\tfrac12\overline{u_i^2u_j}\Big)$ · **why** Chain rule; then $u_j\,\partial\phi/\partial x_j=\partial(u_j\phi)/\partial x_j$ because
     $\partial u_j/\partial x_j=0$ (12.28). · **plain** The fluctuations carry their own energy about: a divergence, so transport only.
  10. **did** The pressure term · **tex** $-\dfrac1{\rho_0}\overline{u_i\dfrac{\partial p}{\partial x_i}}=-\dfrac1{\rho_0}\dfrac{\partial}{\partial x_i}\overline{p\,u_i}$ · **why** Product rule: $\partial(pu_i)/\partial x_i=u_i\,\partial p/\partial x_i+p\,\partial u_i/\partial x_i$ and the last term is zero by $\partial u_i/\partial x_i=0$ (12.28). · **plain** Pressure fluctuations only move energy from place to place.
  11. **did** The buoyancy term · **tex** $\overline{u_i\,g\alpha T'\delta_{i3}}=g\alpha\,\overline{u_3T'}$ · **why** $\delta_{i3}$ selects the vertical component; g and α are constants. ·
      **plain** <span style="color:#3b82f6">Buoyant production</span> when warm fluid rises ($\overline{u_3T'}>0$), destruction when the heat flux is downward.
  12. **did** Rewrite the viscous term with the strain rate · **tex** $\dfrac{\partial^2u_i}{\partial x_j^2}=\dfrac{\partial}{\partial x_j}\big(2S'_{ij}\big),\quad S'_{ij}=\tfrac12\Big(\dfrac{\partial u_i}{\partial x_j}+\dfrac{\partial u_j}{\partial x_i}\Big)$ · **why** $\partial(2S'_{ij})/\partial x_j=\partial^2u_i/\partial x_j^2+\partial(\partial u_j/\partial x_j)/\partial x_i$
      and the last term is zero by $\partial u_j/\partial x_j=0$ (12.28). We want the book's form. · **plain** Viscous force = divergence of the viscous
      stress of the fluctuation.
  13. **did** Split it with the product rule · **tex** $\nu\,u_i\dfrac{\partial}{\partial x_j}\big(2S'_{ij}\big)=\dfrac{\partial}{\partial x_j}\big(2\nu\,u_iS'_{ij}\big)-2\nu\,S'_{ij}\dfrac{\partial u_i}{\partial x_j}$ · **why** The same split as D09 step 5: a divergence
      (transport) and a remainder. · **plain** Viscous transport minus viscous deformation work.
  14. **did** Contract the remainder · **tex** $S'_{ij}\dfrac{\partial u_i}{\partial x_j}=S'_{ij}S'_{ij}\ \Rightarrow\ \bar\varepsilon\equiv2\nu\,\overline{S'_{ij}S'_{ij}}\ \ge0$ · **why** A symmetric tensor contracted with a gradient
      sees only the symmetric part (ch02 §2.10); a sum of squares cannot be negative. · **plain** <span style="color:#f43f5e">Dissipation</span>: always a
      loss.
  15. **did** Assemble · **tex** $\begin{aligned}\dfrac{\partial\bar e}{\partial t}+U_j\dfrac{\partial\bar e}{\partial x_j}=&\ \dfrac{\partial}{\partial x_j}\Big(-\dfrac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\tfrac12\overline{u_i^2u_j}\Big)\\&-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}\end{aligned}$ · **why** Steps 5–14 in step 4, transport
      terms gathered in one divergence. This is (12.47). · **plain** Turbulent energy: transport − dissipation + shear
      production ± buoyancy.
- **Result:** $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47), as in step 15 — "the turbulence is fed by $-\overline{u_iu_j}\,\partial U_i/\partial x_j$ (and by buoyancy when heated from below) and
  drained by $\bar\varepsilon$".
- **Check:** units — m²/s³. Add the two budgets: $+\overline{u_iu_j}\,\partial U_i/\partial x_j$ in the mean-flow budget (D09's Result) and $-\overline{u_iu_j}\,\partial U_i/\partial x_j$ in this Result cancel — total kinetic energy is only transported,
  dissipated or exchanged with potential energy ✓. Isotropic turbulence: production = $\overline{u_1^2}\,\partial U_i/\partial x_i=0$ ✓. Homogeneous
  steady shear flow: production = ε̄. Log layer: $u_*^2\cdot u_*/(\kappa z)$ = 6.75 × 10⁻³ m²/s³ for u_* = 0.3 m/s, z = 10 m.
- **sympy check intent:** run the same moves on a manufactured two-dimensional solenoidal fluctuation and check that each
  rewriting step is an identity. `check_src` sketch:
  ```python
  import sympy as sp                                             # symbolic algebra
  x, y, nu = sp.symbols("x y nu", real=True)                     # two coordinates are enough to test the identities
  psi = sp.Function("psi")(x, y); p = sp.Function("p")(x, y)     # any stream function and pressure
  u = [sp.diff(psi, y), -sp.diff(psi, x)]; X = [x, y]            # a divergence-free fluctuation: (12.28) holds exactly
  S = lambda i, j: (sp.diff(u[i], X[j]) + sp.diff(u[j], X[i]))/2 # fluctuating strain rate S'_ij
  rng = range(2)   # two components are enough to test the identities
  # step 9: u_i u_j d_j u_i  ==  d_j( u_j u_i^2 / 2 )
  lhs9 = sum(u[i]*u[j]*sp.diff(u[i], X[j]) for i in rng for j in rng)   # step 9, left side: u_i u_j ∂u_i/∂x_j
  rhs9 = sum(sp.diff(u[j]*sum(ui**2 for ui in u)/2, X[j]) for j in rng)   # step 9, right side: ∂(½ u_i² u_j)/∂x_j
  assert sp.simplify(lhs9 - rhs9) == 0   # the triple term is a pure divergence (uses ∂u_j/∂x_j = 0)
  # step 10: u_i d_i p == d_i (p u_i)
  assert sp.simplify(sum(u[i]*sp.diff(p, X[i]) for i in rng) - sum(sp.diff(p*u[i], X[i]) for i in rng)) == 0   # step 10: u_i ∂p/∂x_i = ∂(p u_i)/∂x_i
  # steps 12-14: nu u_i lap(u_i) == d_j(2 nu u_i S_ij) - 2 nu S_ij S_ij
  lhs = nu*sum(u[i]*(sp.diff(u[i], x, 2) + sp.diff(u[i], y, 2)) for i in rng)   # steps 12–14, left side: ν u_i ∇²u_i
  rhs = sum(sp.diff(2*nu*u[i]*S(i, j), X[j]) for i in rng for j in rng) - 2*nu*sum(S(i, j)**2 for i in rng for j in rng)   # right side: ∂(2ν u_i S′_ij)/∂x_j − 2ν S′_ij S′_ij
  assert sp.simplify(lhs - rhs) == 0                             # every pointwise identity holds before averaging
  ```
  (The averaging itself — steps 5–8, 11 — is checked by `ch12.tke_budget_sympy()` with the averaging operator of P294, and
  by `reynolds_stress_budget_sympy()`: half the trace of the Reynolds-stress equation (N59 in C04) equals this Result.)
- **What it means:** the three sources and sinks are what every turbulence closure models (C13) and what the Richardson
  numbers compare (C14). **Fails when:** density fluctuations are not small (compressible turbulence has extra pressure–
  dilatation terms).
- **Traps:** forgetting the $+\partial\overline{u_iu_j}/\partial x_j$ in the fluctuation equation (it vanishes only after step 7); the sign of
  production; writing the dissipation as $\nu\overline{(\partial_ju_i)^2}$ — it differs from $2\nu\overline{S'_{ij}S'_{ij}}$ by a divergence that vanishes only in
  homogeneous turbulence; slip #2 (the page's label under the left side names Ē — it is the turbulent ē); slip #10 (use
  $\tfrac12\overline{u_i^2u_j}$ for the triple correlation everywhere).

### D11 · $\bar\varepsilon\sim(\Delta U)^3/L$ (12.48)–(12.49), the Kolmogorov scales (12.50), $\eta u_K/\nu=1$ and $\eta/L\sim\mathrm{Re}_L^{-3/4}$ (12.51) — ★, 9 steps, in C07 (notebook · `energy_cascade_spectrum`)
- **Goal:** estimate how much energy the turbulence dissipates and how small the eddies that do it are. · **Start:** the turbulent energy budget $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47), taken steady, in a shear flow of velocity difference ΔU across a width L. · **Plan:** • production
  balances dissipation • estimate production from the large eddies • the smallest eddies know only ν and ε̄: build a length
  and a velocity from them • compare with L. · **Tools:** order-of-magnitude scaling (P130); exponent matching / Π theorem
  (ch01 C69, `DIM.solve_exponents`). · **Assumptions:** stationary turbulence; high Reynolds number; Kolmogorov's first
  hypothesis (small scales depend on ν and ε̄ only).
- **Steps:**
  1. **did** Balance the budget on average · **tex** $\dot W\equiv-\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}=\bar\varepsilon$ · **why** In stationary turbulence $\partial\bar e/\partial t=0$; integrated over the
     flow the transport (divergence) terms of the Start's budget vanish, leaving production = dissipation. This is $\dot W=\bar\varepsilon$ of (12.49). · **plain**
     What the mean flow supplies is dissipated at the same rate.
  2. **did** Estimate each factor of the production · **tex** $\overline{u_iu_j}\sim(\Delta U)^2,\qquad\dfrac{\partial U_i}{\partial x_j}\sim\dfrac{\Delta U}{L}$ · **why** The large eddies are made by the velocity
     difference, so their velocities are a fraction of ΔU; the mean velocity changes by ΔU over L ("∼" = same order). ·
     **plain** Big eddies: speed like ΔU, size like L.
  3. **did** Multiply · **tex** $\bar\varepsilon=\dot W\sim(\Delta U)^2\Big[\dfrac{\Delta U}{L}\Big]=\dfrac{(\Delta U)^3}{L}$ · **why** Steps 1 and 2. This is (12.48)–(12.49); no viscosity appears. · **plain** The
     dissipation rate is set by the large scales alone.
  4. **did** Seek a length made of ν and ε̄ only · **tex** $\eta=\nu^{a}\,\bar\varepsilon^{\,b}$ · **why** Kolmogorov's hypothesis: the smallest eddies are far removed
     from the geometry; the only quantities they can know are the viscosity and the energy handed down to them. · **plain**
     The smallest scale must be a product of powers of ν and ε̄.
  5. **did** Match units · **tex** $\mathrm m=\big(\mathrm m^2\mathrm s^{-1}\big)^a\big(\mathrm m^2\mathrm s^{-3}\big)^b\ \Rightarrow\ 2a+2b=1,\quad-a-3b=0$ · **why** Both sides must have the same power of metres
     and of seconds (ch01's Π method). · **plain** Two equations for two exponents.
  6. **did** Solve · **tex** $a=\tfrac34,\ b=-\tfrac14\ \Rightarrow\ \eta=\big(\nu^3/\bar\varepsilon\big)^{1/4}$ · **why** From the second equation a = −3b; then −6b + 2b = 1. This is the first
     of (12.50). · **plain** The Kolmogorov length.
  7. **did** Repeat for a velocity and a time · **tex** $u_K=(\nu\bar\varepsilon)^{1/4},\qquad\tau_\eta=(\nu/\bar\varepsilon)^{1/2}$ · **why** Same matching with m s⁻¹ and s on the left; the second of $\eta=(\nu^3/\bar\varepsilon)^{1/4}$, $u_K=(\nu\bar\varepsilon)^{1/4}$ (12.50) (the time scale is ours). · **plain** The speed and lifetime of the smallest eddies.
  8. **did** Form their Reynolds number · **tex** $\dfrac{\eta\,u_K}{\nu}=\dfrac{(\nu^3/\bar\varepsilon)^{1/4}(\nu\bar\varepsilon)^{1/4}}{\nu}=1$ · **why** Multiply the two scales; the powers of ε̄ cancel and ν^{3/4+1/4} =
     ν. · **plain** At the Kolmogorov scale inertia and viscosity are equally strong: the cascade stops there.
  9. **did** Eliminate ε̄ with step 3 · **tex** $\dfrac\eta L\sim\dfrac1L\Big(\dfrac{\nu^3L}{(\Delta U)^3}\Big)^{1/4}=\Big(\dfrac{\nu}{\Delta U\,L}\Big)^{3/4}=\mathrm{Re}_L^{-3/4}$ · **why** Substitute $\bar\varepsilon\sim(\Delta U)^3/L$; $L^{1/4}/L=L^{-3/4}$. This is
     (12.51) with $\mathrm{Re}_L=\Delta UL/\nu$. · **plain** The higher the Reynolds number, the further apart the largest and smallest
     eddies.
- **Result:** $\bar\varepsilon\sim(\Delta U)^3/L$; $\eta=(\nu^3/\bar\varepsilon)^{1/4}$, $u_K=(\nu\bar\varepsilon)^{1/4}$ (12.50); $\eta/L\sim\mathrm{Re}_L^{-3/4}$ (12.51) — "viscosity sets where, not how
  much".
- **Check:** units — ε̄: (m/s)³/m = m²/s³ ✓; η: (m⁶s⁻³/m²s⁻³)^{1/4} = m ✓. Also (ours) $u_K/\Delta U=\mathrm{Re}_L^{-1/4}$,
  $\tau_\eta\Delta U/L=\mathrm{Re}_L^{-1/2}$. Number: ΔU = 1 m/s, L = 1 m, water: η = 32 µm, Re_L = 10⁶, η/L = 3.2 × 10⁻⁵ ✓.
- **What it means:** doubling ν does not change the dissipation rate; it makes the smallest eddies 2^{3/4} = 1.7 times
  larger. A simulation that resolves η needs $(L/\eta)^3\sim\mathrm{Re}_L^{9/4}$ points. **Fails when:** Re is low (no separation), or the
  turbulence is not stationary (then ε̄ lags the supply).
- **Traps:** reading "∼" as "="; forgetting ε̄ is per unit mass; believing ε̄ depends on ν because ν appears in its
  definition $2\nu\overline{S'_{ij}S'_{ij}}$ (the gradients adjust).

### D12 · The universal spectrum (12.53), Kolmogorov's law $S_{11}=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ (12.54, corrected) and the normalisation (12.55) — ★★, 9 steps, in C08 (notebook · `energy_cascade_spectrum`)
- **Goal:** find the shape of the spectrum between the large eddies and the smallest ones without solving any equation. ·
  **Start:** $S_{11}(k_1)=\dfrac1{2\pi}\displaystyle\int_{-\infty}^{+\infty}R_{11}(r_1)\,e^{-ik_1r_1}\,dr_1$ (12.45). · **Plan:** • find the units of $S_{11}$ • at small scales let it
  depend on ε̄, ν, k₁: two dimensionless groups • in the inertial range drop ν: one group, a power law • fix the constant's
  convention. · **Tools:** units of a spectral density (P293); Π theorem (ch01 C69); log–log slopes (P13). · **Assumptions:**
  high Reynolds number; local isotropy; the inertial range knows neither L nor ν (Kolmogorov's second hypothesis).
- **Steps:**
  1. **did** Find the units of the spectrum · **tex** $[S_{11}]=[R_{11}]\,[r_1]=\mathrm m^2\mathrm s^{-2}\times\mathrm m=\mathrm m^3\mathrm s^{-2}$ · **why** In $S_{11}(k_1)=\frac1{2\pi}\int_{-\infty}^{+\infty}R_{11}(r_1)\,e^{-ik_1r_1}\,dr_1$ (12.45) a correlation (velocity²) is
     integrated over a distance; the exponential and 1/2π have no units. · **plain** A spectrum is "variance per unit
     wavenumber".
  2. **did** List what small-scale eddies can depend on · **tex** $S_{11}=\mathrm{fn}(\bar\varepsilon,\ \nu,\ k_1)$ · **why** Far down the cascade the mean shear is forgotten
     (D11 step 4); the wavenumber says which eddy size we ask about. Four quantities, two dimensions (m, s) ⇒ two groups. ·
     **plain** One dimensionless spectrum as a function of one dimensionless wavenumber.
  3. **did** Form the two groups · **tex** $\dfrac{S_{11}(k_1)}{\nu^{5/4}\bar\varepsilon^{1/4}}=\Phi\Big(\dfrac{k_1\nu^{3/4}}{\bar\varepsilon^{1/4}}\Big)$ · **why** $\nu^{5/4}\bar\varepsilon^{1/4}$ has units m^{10/4+2/4}s^{−5/4−3/4} = m³s⁻²;
     $\nu^{3/4}\bar\varepsilon^{-1/4}$ is a length. This is the first form of (12.53). · **plain** In these units every flow's small-scale
     spectrum is the same curve Φ.
  4. **did** Recognise the Kolmogorov scales · **tex** $\dfrac{S_{11}(k_1)}{u_K^2\,\eta}=\Phi(k_1\eta)\quad\text{for }k_1\gg2\pi/L$ · **why** $u_K^2\eta=(\nu\bar\varepsilon)^{1/2}(\nu^3/\bar\varepsilon)^{1/4}=\nu^{5/4}\bar\varepsilon^{1/4}$ and
     $\eta=\nu^{3/4}\bar\varepsilon^{-1/4}$ by (12.50). The second form of (12.53). · **plain** Velocity² × length of the smallest eddies is the
     natural unit.
  5. **did** In the inertial range drop the viscosity · **tex** $S_{11}=C_1\,\bar\varepsilon^{\,a}\,k_1^{\,b}:\quad\mathrm m^3\mathrm s^{-2}=\big(\mathrm m^2\mathrm s^{-3}\big)^a\big(\mathrm m^{-1}\big)^b$ · **why** For $2\pi/L\ll k_1\ll2\pi/\eta$ eddies are too
     small to feel L and too large to feel ν. Three quantities, two dimensions ⇒ one group ⇒ a pure power law. · **plain**
     Only one combination of ε̄ and k₁ has the units of a spectrum.
  6. **did** Solve for the exponents · **tex** $-3a=-2,\ \ 2a-b=3\ \Rightarrow\ a=\tfrac23,\ b=-\tfrac53:\qquad S_{11}(k_1)=C_1\,\bar\varepsilon^{2/3}\,k_1^{-5/3}$ · **why** Seconds give a; metres then give b = 2a − 3 =
     −5/3. This is (12.54) with the exponent the page misprints as +5/3 (slip #1). · **plain** The −5/3 law.
  7. **did** Check it against the universal form · **tex** $\Phi(k_1\eta)=C_1\,(k_1\eta)^{-5/3}\ \Rightarrow\ S_{11}=C_1u_K^2\eta\,(k_1\eta)^{-5/3}=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ · **why** $u_K^2\eta^{-2/3}=(\nu\bar\varepsilon)^{1/2}(\nu^3/\bar\varepsilon)^{-1/6}=\bar\varepsilon^{2/3}$:
     ν cancels only for the exponent −5/3. · **plain** The inertial range is the part of Φ where viscosity drops out.
  8. **did** State the range of validity · **tex** $\dfrac{2\pi}{L}\ll k_1\ll\dfrac{2\pi}{\eta},\qquad\dfrac{L}{\eta}\sim\mathrm{Re}_L^{3/4}$ · **why** Both conditions of step 5; their ratio is D11's $\eta/L\sim\mathrm{Re}_L^{-3/4}$ (12.51). · **plain** The −5/3 range is ¾ log₁₀Re_L decades wide at most.
  9. **did** Fix what the constant refers to · **tex** $\displaystyle\int_{-\infty}^{+\infty}S_{11}(k_1)\,dk_1=\overline{u_1^2}=\int_0^{+\infty}2S_{11}(k_1)\,dk_1$ · **why** $S_{11}$ is even (D05), so one-sided
     plots show $2S_{11}$. This is (12.55): a quoted $C_1$ is for the two-sided form; the one-sided constant is twice as
     large. · **plain** Two conventions, a factor 2.
- **Result:** $S_{11}/(u_K^2\eta)=\Phi(k_1\eta)$ (12.53); $S_{11}(k_1)=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ for $2\pi/L\ll k_1\ll2\pi/\eta$ (12.54) — "in the inertial range
  the spectrum falls 46-fold per decade".
- **Check:** units — (m²/s³)^{2/3}(1/m)^{−5/3} = m^{4/3+5/3}s⁻² = m³/s² ✓; the printed $k_1^{+5/3}$ gives m^{−1/3}s⁻² ✗ and an infinite
  variance. Isotropy links the constants: one-sided $C_1=\tfrac{18}{55}C$ ≈ 0.49 for C ≈ 1.5 (Sreenivasan 1995); two-sided 0.245.
  Code: `fit_inertial_range` on the model spectrum gives slope −1.667.
- **What it means:** the same argument with a temperature dimension gives the Obukhov–Corrsin spectrum
  $S_T\propto\bar\varepsilon_T\bar\varepsilon^{-1/3}K^{-5/3}$ (12.113), and with a length gives Richardson's $K\sim\bar\varepsilon^{1/3}l^{4/3}$. **Fails when:** Re_L is too low
  for a range to exist, or other parameters enter (stratification: N; rotation: f — Ch. 13).
- **Traps:** the sign of the exponent (slip #1); one- vs two-sided constants; using k in cycles per metre with formulas in
  rad/m; expecting the law at the largest scales.

### D13 · The thin-layer jet equations (12.61) and the momentum-flux invariant $J_s=\rho\int U^2dy$ (12.62) — ★★, 10 steps, in C09 (notebook · `turbulent_jet_similarity`)
- **Goal:** simplify the mean equations for a slender jet and find the quantity that is the same at every distance from the
  nozzle. · **Start:** $\dfrac{\partial U}{\partial x}+\dfrac{\partial V}{\partial y}=0$ (12.58), $U\dfrac{\partial U}{\partial x}+V\dfrac{\partial U}{\partial y}=-\dfrac1\rho\dfrac{\partial P}{\partial x}+\nu\Big(\dfrac{\partial^2U}{\partial x^2}+\dfrac{\partial^2U}{\partial y^2}\Big)-\dfrac{\partial\overline{u^2}}{\partial x}-\dfrac{\partial\overline{uv}}{\partial y}$ (12.59) and
  $U\dfrac{\partial V}{\partial x}+V\dfrac{\partial V}{\partial y}=-\dfrac1\rho\dfrac{\partial P}{\partial y}+\nu\Big(\dfrac{\partial^2V}{\partial x^2}+\dfrac{\partial^2V}{\partial y^2}\Big)-\dfrac{\partial\overline{uv}}{\partial x}-\dfrac{\partial\overline{v^2}}{\partial y}$ (12.60). · **Plan:** • compare term sizes for a thin layer • the y-equation gives the pressure •
  drop the small terms • add continuity to make a conservation form • integrate across the jet. · **Tools:** two-length
  scaling (P188, ch09); product rule (P38); integrals over an infinite range (P144). · **Assumptions:** steady in the mean;
  constant density; high Reynolds number; slender (δ ≪ x); still surroundings with uniform pressure.
- **Steps:**
  1. **did** Set the scales · **tex** $U\sim U_{CL},\quad x\sim L_x,\quad y\sim\delta\ll L_x,\quad V\sim U_{CL}\,\delta/L_x$ · **why** The jet is long and thin; continuity $\partial U/\partial x+\partial V/\partial y=0$ (12.58) makes
     $\partial V/\partial y\sim\partial U/\partial x$, which fixes the size of V (Ch. 9's boundary-layer scaling). · **plain** Cross-stream changes are fast,
     stream-wise ones slow, and V is small.
  2. **did** Keep the largest terms of the y-equation · **tex** $0\cong-\dfrac1\rho\dfrac{\partial}{\partial y}\big(P+\rho\overline{v^2}\big)$ · **why** In the cross-stream equation of the Start the inertia terms are ~ $U_{CL}^2\delta/L_x^2$, the
     stress term $\partial\overline{v^2}/\partial y\sim u'^2/\delta$ is larger by $(L_x/\delta)^2(u'/U_{CL})^2$; only the pressure can balance it. This is the second of (12.61), $0\cong-\frac1\rho\frac{\partial}{\partial y}\big(P+\rho\overline{v^2}\big)$.
     · **plain** Across the jet, pressure and the normal Reynolds stress add to a constant.
  3. **did** Integrate across the jet · **tex** $P+\rho\,\overline{v^2}=P_\infty$ · **why** Integrate step 2 in y out to the still fluid, where $\overline{v^2}=0$ and the
     pressure is the uniform $P_\infty$. · **plain** Inside the jet the mean pressure is slightly below ambient.
  4. **did** Differentiate in x · **tex** $-\dfrac1\rho\dfrac{\partial P}{\partial x}=\dfrac{\partial\overline{v^2}}{\partial x}$ · **why** $P_\infty$ is constant. This term is ~ $u'^2/L_x$, as small as $\partial\overline{u^2}/\partial x$ and
     smaller than $\partial\overline{uv}/\partial y\sim u'^2/\delta$ by δ/L_x. · **plain** The stream-wise pressure gradient is negligible.
  5. **did** Drop the small terms of the x-equation · **tex** $U\dfrac{\partial U}{\partial x}+V\dfrac{\partial U}{\partial y}\cong-\dfrac{\partial\overline{uv}}{\partial y}$ · **why** Viscous terms are ~ 1/Re of the inertia terms;
     pressure and $\overline{u^2}$ gradients were shown small in step 4. This is the first of (12.61), $Urac{\partial U}{\partial x}+Vrac{\partial U}{\partial y}\cong-rac{\partial\overline{uv}}{\partial y}$. · **plain** Mean advection is balanced by
     the cross-stream gradient of the Reynolds shear stress.
  6. **did** Add U times the continuity equation · **tex** $U\dfrac{\partial U}{\partial x}+V\dfrac{\partial U}{\partial y}+U\Big(\dfrac{\partial U}{\partial x}+\dfrac{\partial V}{\partial y}\Big)=\dfrac{\partial(U^2)}{\partial x}+\dfrac{\partial(VU)}{\partial y}$ · **why** The bracket is zero by $\partial U/\partial x+\partial V/\partial y=0$ (12.58), so
     adding it changes nothing; the product rule then makes exact derivatives. · **plain** Advection as the divergence of a
     momentum flux.
  7. **did** Write the conservation form · **tex** $\dfrac{\partial}{\partial x}\big(U^2\big)+\dfrac{\partial}{\partial y}\big(VU+\overline{uv}\big)\cong0$ · **why** Steps 5 and 6, everything on one side. · **plain** Whatever
     x-momentum leaves along the jet must have come in from the sides — and nothing does.
  8. **did** Integrate across the jet · **tex** $\dfrac{d}{dx}\displaystyle\int_{-\infty}^{+\infty}U^2\,dy+\Big[VU+\overline{uv}\Big]_{y=-\infty}^{y=+\infty}\cong0$ · **why** Integrate each term over y; the x-derivative
     comes outside because the limits are fixed. · **plain** Change of momentum flux along the jet = flux through the edges.
  9. **did** Evaluate the edge term · **tex** $\Big[VU+\overline{uv}\Big]_{-\infty}^{+\infty}=0$ · **why** Outside the jet the fluid has no stream-wise velocity (U → 0) and
     no turbulence ($\overline{uv}\to0$) — although V does **not** vanish there (entrainment). · **plain** The inflowing fluid brings
     volume but no x-momentum.
  10. **did** Integrate in x and multiply by the density · **tex** $J_s\equiv\rho_s\displaystyle\int_{-\infty}^{+\infty}\big[U^2\big]_{x=0}\,dy\cong\rho\int_{-\infty}^{+\infty}U^2\,dy=\mathit{const.}$ · **why** A zero derivative means
      a constant, fixed at the nozzle (x = 0). This is (12.62); downstream the density is the ambient ρ because the jet
      fluid is diluted. · **plain** The momentum flux per unit span is the same at every station.
- **Result:** $U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}\cong-\frac{\partial\overline{uv}}{\partial y}$, $0\cong-\frac1\rho\frac{\partial}{\partial y}\big(P+\rho\overline{v^2}\big)$ (12.61) and $J_s=\rho\int_{-\infty}^{+\infty}U^2dy=\mathit{const.}$ (12.62) — "a jet conserves its momentum flux".
- **Check:** units — ρU²·y: kg m⁻³ × m² s⁻² × m = N/m ✓. A top-hat slot of width d: $J_s=\rho_sU_0^2d$. Code: `jet_momentum_flux_per_span`
  flat to 1e-6 along the similarity solution; `thin_shear_layer_terms` prints the dropped terms, of relative order (δ/x)² — small for ξ½ = 0.1: at x = 1 m, y = 0.1 m in air the viscous term is about 0.001 of the stress gradient and V/U ≈ 0.02 (printed by the cell above the derivation).
- **What it means:** with no wall and no pressure gradient nothing can change the jet's momentum; everything else (width,
  speed, volume flux) must arrange itself around this invariant (D15). **Fails when:** there is a co-flow or pressure
  gradient, buoyancy (plumes: momentum grows), or a nearby wall (Coanda effect).
- **Traps:** concluding V → 0 at the edges; dropping $\overline{v^2}$ and then keeping ∂P/∂x; forgetting that the volume flux is *not*
  conserved.

### D14 · The similarity equation of the plane jet (12.63) — ★★★, 14 steps, in C09 (notebook · `turbulent_jet_similarity`)
- **Goal:** put the self-preserving forms into the jet equation and see what conditions on the scales $U_{CL}(x)$, δ(x), Ψ(x)
  make a single profile shape possible. The book calls this "somewhat tedious" and skips it. · **Start:**
  $U\dfrac{\partial U}{\partial x}+V\dfrac{\partial U}{\partial y}\cong-\dfrac{\partial\overline{uv}}{\partial y}$ (12.61) with $U=U_{CL}(x)F(\xi)$ (12.56), $-\overline{uv}=\Psi(x)G(\xi)$ (12.57), $\xi=y/\delta(x)$. · **Plan:** • eliminate V by
  continuity • differentiate the similarity forms (ξ moves with x!) • do the y-integral by parts • add the terms: two cancel
  • divide by $U_{CL}^2/\delta$. · **Tools:** chain rule with a moving similarity variable (P206); integration by parts (P218a);
  change of variable dy = δ dξ (P106); sympy (P40). · **Assumptions:** self-preservation with one length δ(x); symmetric jet
  (V = 0 on the axis).
- **Steps:**
  1. **did** Solve continuity for V · **tex** $V(x,y)=-\displaystyle\int_0^y\dfrac{\partial U}{\partial x}\,dy$ · **why** Integrate $\partial V/\partial y=-\partial U/\partial x$ (12.58) from the axis, where V = 0 by
     symmetry. · **plain** The cross-flow is whatever is needed to make up for the slowing stream.
  2. **did** Eliminate V · **tex** $U\dfrac{\partial U}{\partial x}-\Big[\displaystyle\int_0^y\dfrac{\partial U}{\partial x}\,dy\Big]\dfrac{\partial U}{\partial y}\cong-\dfrac{\partial\overline{uv}}{\partial y}$ · **why** Substitute step 1 in $U\,\partial U/\partial x+V\,\partial U/\partial y\cong-\partial\overline{uv}/\partial y$ (12.61); now only U and $\overline{uv}$ remain. · **plain** One
     equation for the mean velocity and the shear stress.
  3. **did** Differentiate the similarity variable · **tex** $\dfrac{\partial\xi}{\partial x}=-\dfrac{y\,\delta'}{\delta^2}=-\dfrac{\xi\,\delta'}{\delta},\qquad\dfrac{\partial\xi}{\partial y}=\dfrac1\delta$ · **why** ξ = y/δ(x): at fixed y, ξ decreases as the
     jet widens (P206). A prime on δ, $U_{CL}$ means d/dx. · **plain** A fixed point in space slides inward in the scaled
     coordinate as the jet grows.
  4. **did** Differentiate U in x · **tex** $\dfrac{\partial U}{\partial x}=U_{CL}'\,F-U_{CL}\,F'\,\dfrac{\xi\,\delta'}{\delta}$ · **why** Product rule on $U_{CL}(x)F(\xi)$, then the chain rule with step 3 (F′ =
     dF/dξ). · **plain** The speed at a point changes because the centreline slows and because the profile widens.
  5. **did** Differentiate U in y · **tex** $\dfrac{\partial U}{\partial y}=\dfrac{U_{CL}}{\delta}\,F'$ · **why** Chain rule with $\partial\xi/\partial y=1/\delta$. · **plain** The cross-stream slope.
  6. **did** Form the first term · **tex** $U\dfrac{\partial U}{\partial x}=U_{CL}U_{CL}'\,F^2-U_{CL}^2\,\dfrac{\delta'}{\delta}\,\xi FF'$ · **why** Multiply step 4 by $U=U_{CL}F$. · **plain** Stream-wise advection in
     similarity form.
  7. **did** Change the integration variable · **tex** $\displaystyle\int_0^y\dfrac{\partial U}{\partial x}\,dy=\delta\,U_{CL}'\int_0^\xi F\,d\xi-U_{CL}\,\delta'\int_0^\xi\xi F'\,d\xi$ · **why** dy = δ dξ at fixed x; the factors that depend only
     on x come outside. · **plain** The integral becomes two integrals over the profile shape.
  8. **did** Integrate the second one by parts · **tex** $\displaystyle\int_0^\xi\xi F'\,d\xi=\xi F-\int_0^\xi F\,d\xi$ · **why** $\int u\,dv=uv-\int v\,du$ with u = ξ, dv = F′dξ (P218a); the
     boundary term at 0 vanishes. · **plain** Everything is now in terms of F and its running integral.
  9. **did** Collect V · **tex** $V=-\displaystyle\int_0^y\dfrac{\partial U}{\partial x}dy=U_{CL}\,\delta'\,\xi F-\big(\delta\,U_{CL}'+U_{CL}\,\delta'\big)I,\qquad I(\xi)\equiv\int_0^\xi F\,d\xi$ · **why** Steps 7 and 8 with the overall minus
     sign of step 1. · **plain** The cross-flow in similarity form: an outward part and an inward (entraining) part.
  10. **did** Form the second term · **tex** $V\dfrac{\partial U}{\partial y}=U_{CL}^2\,\dfrac{\delta'}{\delta}\,\xi FF'-\Big(U_{CL}U_{CL}'+U_{CL}^2\dfrac{\delta'}{\delta}\Big)F'I$ · **why** Multiply step 9 by step 5. · **plain** Cross-stream
      advection.
  11. **did** Add the two advection terms · **tex** $U\dfrac{\partial U}{\partial x}+V\dfrac{\partial U}{\partial y}=U_{CL}U_{CL}'\,F^2-\Big(U_{CL}U_{CL}'+U_{CL}^2\dfrac{\delta'}{\delta}\Big)F'I$ · **why** The two $\xi FF'$ terms of steps 6 and 10
      are equal and opposite and cancel. · **plain** Only two terms survive on the left.
  12. **did** Differentiate the stress · **tex** $-\dfrac{\partial\overline{uv}}{\partial y}=\dfrac{\partial}{\partial y}\big(\Psi\,G\big)=\dfrac{\Psi}{\delta}\,G'$ · **why** $-\overline{uv}=\Psi(x)G(\xi)$ (12.57) and the chain rule of step 5. · **plain** The
      stress gradient in similarity form.
  13. **did** Multiply both sides by $\delta/U_{CL}^2$ · **tex** $\dfrac{\delta\,U_{CL}'}{U_{CL}}F^2-\Big(\dfrac{\delta\,U_{CL}'}{U_{CL}}+\delta'\Big)F'I=\dfrac{\Psi}{U_{CL}^2}G'$ · **why** This makes every coefficient
      dimensionless and puts all the x-dependence into three groups. · **plain** Functions of x multiply functions of ξ.
  14. **did** Write it in the book's form · **tex** $\Big\{\dfrac{\delta U_{CL}'}{U_{CL}}\Big\}F^2-\Big\{\dfrac{\delta U_{CL}'}{U_{CL}}+\delta'\Big\}F'\displaystyle\int_0^\xi F\,d\xi=\Big\{\dfrac{\Psi}{U_{CL}^2}\Big\}G'$ · **why** Replace I by its definition. This is
      (12.63). · **plain** Three brackets that depend on x only; three shape functions that depend on ξ only.
- **Result:** $\big\{\frac{\delta U'_{CL}}{U_{CL}}\big\}F^2-\big\{\frac{\delta U'_{CL}}{U_{CL}}+\delta'\big\}F'\int_0^\xi F\,d\xi=\big\{\frac\Psi{U_{CL}^2}\big\}G'$ (12.63), as in step 14 — "self-preservation is possible only if the three brackets behave alike in x".
- **Check:** units — every bracket is dimensionless (δ × (1/x) ; Ψ/U² with Ψ in m²/s²). With δ = x, $U_{CL}\propto x^{-1/2}$: brackets
  −½, +½, $C_3$, and the left side is the exact derivative $-\tfrac12(FI)'$ (used in D15 step 8). Symmetry: left side even in ξ,
  so G′ is even and G odd ✓.
- **sympy check intent:** take concrete power-law scales and a Gaussian F, compute U, V from continuity and $U\,\partial U/\partial x+V\,\partial U/\partial y$ (the left side of the thin-layer equation in the Start) directly, and compare with $(U_{CL}^2/\delta)\big[\{\delta U'_{CL}/U_{CL}\}F^2-\{\delta U'_{CL}/U_{CL}+\delta'\}F'\int_0^\xi F\,d\xi\big]$, the left side of the Result. `check_src` sketch:
  ```python
  import sympy as sp                                              # symbolic algebra
  x, y, s, a, c4 = sp.symbols("x y s a c4", positive=True)        # s is a dummy for the y-integral
  m, n = sp.Rational(1), sp.Rational(-1, 2)                       # delta ~ x^m, U_CL ~ x^n (any m, n work)
  delta, UCL = x**m, c4*x**n                                      # the two scales
  F = lambda xi: sp.exp(-a*xi**2)                                 # a Gaussian profile shape
  U = UCL*F(y/delta)                                              # (12.56)
  V = -sp.integrate(sp.diff(UCL*F(s/delta), x), (s, 0, y))        # step 1: continuity
  lhs = sp.simplify(U*sp.diff(U, x) + V*sp.diff(U, y))            # left side of (12.61), computed directly
  xi = sp.symbols("xi", positive=True)   # the similarity variable ξ = y/δ
  I = sp.integrate(F(s), (s, 0, xi))                              # running integral of F
  A1 = delta*sp.diff(UCL, x)/UCL; A2 = A1 + sp.diff(delta, x)     # the first two brackets of (12.63)
  rhs = (UCL**2/delta)*(A1*F(xi)**2 - A2*sp.diff(F(xi), xi)*I)    # steps 11-13, times U_CL^2/delta
  assert sp.simplify(lhs.subs(y, xi*delta) - rhs) == 0            # identical for every x and xi
  assert sp.simplify(A1 - n) == 0 and sp.simplify(A2 - (m + n)) == 0   # brackets are the constants n and m + n
  ```
- **What it means:** the equation separates into x-dependent brackets and ξ-dependent shapes — the door to D15. **Fails
  when:** two lengths matter (near the nozzle; wall jets), so no single ξ exists.
- **Traps:** forgetting that ξ depends on x (dropping the $\xi\delta'/\delta$ term); sign of V; integrating $\xi F'$ wrongly; missing the
  cancellation of the two ξFF′ terms.

### D15 · $\delta\propto x$, $U_{CL}\propto x^{-1/2}$ from (12.64)–(12.65), the far field (12.66)–(12.67) and $\dot V\propto x^{1/2}$ (12.68) — ★★, 10 steps, in C09 (notebook · `turbulent_jet_similarity`)
- **Goal:** get the growth and decay laws of the plane jet from the similarity equation and the invariant. · **Start:**
  $\Big\{\dfrac{\delta U_{CL}'}{U_{CL}}\Big\}F^2-\Big\{\dfrac{\delta U_{CL}'}{U_{CL}}+\delta'\Big\}F'\displaystyle\int_0^\xi F\,d\xi=\Big\{\dfrac{\Psi}{U_{CL}^2}\Big\}G'$ (12.63) and $J_s=\rho\displaystyle\int_{-\infty}^{+\infty}U^2dy=\mathit{const.}$ (12.62). · **Plan:** • each bracket must be
  a constant • integrate the two resulting ODEs • use the invariant to fix the exponent • integrate the profile for the
  volume flux. · **Tools:** separation argument (P167); separable ODE; exponent matching (P197); Gaussian integral (P194). ·
  **Assumptions:** the simplest similarity (brackets constant); virtual origin neglected.
- **Steps:**
  1. **did** Require each bracket to be constant · **tex** $\dfrac{\delta U_{CL}'}{U_{CL}}=C_1,\qquad\dfrac{\delta U_{CL}'}{U_{CL}}+\delta'=C_2,\qquad\dfrac{\Psi}{U_{CL}^2}=C_3$ · **why** F, G depend on ξ only; if the brackets varied
     differently with x the equation could not hold at every x for one pair of shapes (P167). This is (12.64). · **plain** Two
     ODEs for the scales and one algebraic relation.
  2. **did** Subtract the first from the second · **tex** $\delta'=C_2-C_1\ \Rightarrow\ \delta(x)=(C_2-C_1)(x-x_o)$ · **why** A constant derivative integrates to a straight
     line; $x_o$ is the integration constant ("virtual origin"). · **plain** The jet widens linearly.
  3. **did** Choose the convention for δ · **tex** $C_2-C_1=1,\ x_o\approx0\ \Rightarrow\ \delta=x,\quad\xi=y/x$ · **why** δ is only a scale, so its constant factor can be
     absorbed into the shape F; $x_o$ is of the order of the slot width, negligible far downstream. · **plain** We measure y in
     units of the distance from the nozzle.
  4. **did** Integrate the first condition · **tex** $\dfrac{x\,U_{CL}'}{U_{CL}}=C_1\ \Rightarrow\ \ln U_{CL}=C_1\ln x+\text{const}\ \Rightarrow\ U_{CL}=C_4\,x^{\gamma},\ \gamma\equiv C_1$ · **why** With δ = x the equation is
     separable: $dU_{CL}/U_{CL}=C_1\,dx/x$. · **plain** The centreline speed is a power of x — which power is not yet known.
  5. **did** Put the similarity form in the invariant · **tex** $J_s=\rho\,U_{CL}^2\,\delta\displaystyle\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi=\rho\,C_4^2\,x^{2\gamma+1}\int_{-\infty}^{+\infty}F^2(\xi)\,d\xi$ · **why** $U=U_{CL}F$, dy = δ dξ, δ = x. This
     is (12.65). · **plain** The momentum flux as a power of x times a pure number.
  6. **did** Demand that it does not depend on x · **tex** $2\gamma+1=0\ \Rightarrow\ \gamma=-\tfrac12$ · **why** $J_s=\rho\int U^2dy=\mathit{const.}$ (12.62) is the same at every x, so the exponent of x
     must be zero (P197). · **plain** Conservation of momentum fixes the decay exponent.
  7. **did** Write the far-field velocity · **tex** $U(x,y)=C_5\,(J_s/\rho)^{1/2}\,x^{-1/2}\,F(y/x),\qquad C_5=\Big(\displaystyle\int_{-\infty}^{+\infty}F^2\,d\xi\Big)^{-1/2}$ · **why** From step 5 with γ = −½:
     $C_4=(J_s/\rho)^{1/2}(\int F^2)^{-1/2}$. This is (12.66); $C_5=C_4(\rho/J_s)^{1/2}$ (the text's "C₄" after $-\overline{uv}=C_3C_5^2(J_s/\rho)x^{-1}G(y/x)$ (12.67) is this C₅ — slip
     #4). · **plain** The whole mean velocity field from one invariant and one shape.
  8. **did** Integrate the similarity equation once · **tex** $-\tfrac12\big(F\,I\big)'=C_3\,G'\ \Rightarrow\ C_3\,G(\xi)=-\tfrac12F(\xi)\displaystyle\int_0^\xi F\,d\xi$ · **why** With $C_1=-\tfrac12$, $C_2=+\tfrac12$ the left side of the similarity equation (D14's Result) is $-\tfrac12F^2-\tfrac12F'I=-\tfrac12(FI)'$; G(0) = 0. Then $-\overline{uv}=C_3U_{CL}^2G=C_3C_5^2(J_s/\rho)x^{-1}G(y/x)$ (12.67). ·
     **plain** The shear-stress shape follows from the velocity shape; it is negative where dU/dy is negative.
  9. **did** Integrate the velocity across the jet · **tex** $\dot V(x)=\displaystyle\int_{-\infty}^{+\infty}U\,dy=C_5\,(J_s/\rho)^{1/2}\,x^{+1/2}\int_{-\infty}^{+\infty}F(\xi)\,d\xi$ · **why** dy = x dξ turns $x^{-1/2}$ into
     $x^{+1/2}$. This is (12.68). · **plain** The volume flux grows downstream: the jet entrains.
  10. **did** Differentiate for the entrainment velocity · **tex** $v_e\equiv-V(+\infty)=\tfrac12\dfrac{d\dot V}{dx}=\dfrac{\dot V}{4x},\qquad\dfrac{U_{CL}\,\delta}{\nu}\propto x^{1/2}$ · **why** Mass balance: the extra
      flux enters equally through both edges (D14 step 9 at ξ → ∞ gives the same). The local Reynolds number is the
      product of the two scales. · **plain** Still fluid drifts in at a speed ∝ $x^{-1/2}$, and the jet becomes *more* turbulent
      downstream.
- **Result:** δ = x, $U=C_5(J_s/\rho)^{1/2}x^{-1/2}F(y/x)$ (12.66), $\dot V\propto x^{1/2}$ (12.68) — "linear spreading, inverse-square-root decay,
  entrainment".
- **Check:** units — $(J_s/\rho)^{1/2}x^{-1/2}$: (m³/s²)^{1/2} m^{−1/2} = m/s ✓. Gaussian F with illustrative ξ½ = 0.10: $\int F^2$ = 0.1505, $C_5$ =
  2.577; $J_s$ recovered to 1e-6 at x = 1, 2, 4; $\dot V(4)/\dot V(1)=2$. Laminar jet (Ch. 9): δ ∝ x^{2/3}, U ∝ x^{−1/3} — different
  because there ν (a constant) replaces an eddy viscosity that grows as $U_{CL}\delta\propto x^{1/2}$.
- **What it means:** no turbulence model was used — only conservation and shape-keeping; experiments supply two numbers
  (the spreading rate and $C_3$). **Fails when:** near the nozzle (x ≲ 20 d), or if the constants are not universal (N125's
  more general similarity).
- **Traps:** thinking C₁, C₂ are both free (their difference is a convention, γ = C₁ is fixed by $J_s$); the sign of G;
  confusing the constants' names (slip #4); assuming the volume flux is conserved.

### D16 · The exponents of Table 12.1 from each flow's invariant and growth law — ★★, 12 steps, in C09 (notebook · `turbulent_jet_similarity`)
- **Goal:** derive, for every free shear flow of Table 12.1, how its width and velocity scale vary with distance — the book
  prints only the table. · **Start:** a velocity scale $U_s\propto x^{n}$ and a width $\delta\propto x^{m}$ (and, for scalars and buoyancy, a
  third scale ∝ $x^{b}$). · **Plan:** • one equation from what is conserved • one from how fast the layer can widen • solve two
  linear equations — flow by flow. · **Tools:** D13–D15 as the template; exponent matching (P197); exact fractions
  (`fractions.Fraction`, P60). · **Assumptions:** far field; high Re; self-preservation; small deficit for wakes; Boussinesq
  for plumes.
- **Steps:**
  1. **did** State the growth law for flows that carry themselves · **tex** $\dfrac{d\delta}{dx}\sim\dfrac{u'}{U_{\text{adv}}}=\dfrac{u'}{U_s}=\text{const}\ \Rightarrow\ m=1$ · **why** The layer widens at the
     turbulent velocity u′ ∝ $U_s$ while being carried downstream at $U_s$ itself; the ratio is a pure number (the content
     of $\delta'=C_2-C_1$ in D15). Holds for jets, plumes, shear layers. · **plain** Self-propelled flows spread linearly.
  2. **did** Recall the plane jet · **tex** $U_s^2\,\delta=\text{const}\ \Rightarrow\ 2n+m=0\ \Rightarrow\ n=-\tfrac12$ · **why** Momentum flux per span, $J_s=\rho\int U^2dy$ (12.62), with step 1.
     The template. · **plain** Plane jet: (m, n) = (1, −½).
  3. **did** Round jet: integrate the momentum flux over a disc · **tex** $\displaystyle\int_0^\infty\rho U^2\,2\pi r\,dr\propto U_s^2\,\delta^2=\text{const}\ \Rightarrow\ 2n+2m=0\ \Rightarrow\ n=-1$ · **why** The area element
     is 2πr dr, so the invariant carries δ² instead of δ; m = 1 from step 1. · **plain** A round jet slows faster, as 1/x.
  4. **did** Round jet: scalar and Reynolds number · **tex** $Y_s\,U_s\,\delta^2=\text{const}\ \Rightarrow\ Y_s\propto x^{-1};\qquad\dfrac{U_s\delta}{\nu}\propto x^{\,n+m}=x^{0}$ · **why** The flux of nozzle fluid is
     conserved (the analogue of $\dot M_s\cong\rho\int\bar YU\,dy$ (12.70)); the local Reynolds number is the product of the scales. · **plain** A round jet
     keeps the Reynolds number it started with.
  5. **did** Wakes: change the growth law · **tex** $\dfrac{d\delta}{dx}\sim\dfrac{u'}{U_\infty}\sim\dfrac{\Delta U}{U_\infty}\ \Rightarrow\ m-1=n$ · **why** A wake's deficit ΔU ≪ $U_\infty$ is carried at the
     free-stream speed (advection linearised to $U_\infty\,\partial/\partial x$) while it spreads at u′ ∝ ΔU: the ratio is no longer constant.
     · **plain** Wakes spread ever more slowly as their deficit fades.
  6. **did** Plane wake: conserve the momentum deficit · **tex** $\rho\,U_\infty\displaystyle\int\Delta U\,dy=\text{drag per span}\ \Rightarrow\ n+m=0\ \Rightarrow\ m=\tfrac12,\ n=-\tfrac12$ · **why** The deficit flux equals the
     body's drag (Ch. 4); combine with step 5: m − 1 = −m. · **plain** Plane wake: width ∝ √x, deficit ∝ 1/√x.
  7. **did** Round wake · **tex** $\Delta U\,\delta^2=\text{const}\ \Rightarrow\ n+2m=0,\ m-1=n\ \Rightarrow\ m=\tfrac13,\ n=-\tfrac23$ · **why** Disc area again; then m − 1 = −2m. Local Reynolds number
     ∝ $x^{n+m}=x^{-1/3}$. · **plain** A round wake's Reynolds number falls: far enough downstream it can relaminarise.
  8. **did** Plumes: conserve the buoyancy flux · **tex** plane plume: $g'\,U_s\,\delta=\text{const}\ \Rightarrow\ b+n+m=0,\qquad g'\equiv g\,\Delta\rho/\rho\propto x^{b}$ · **why** A plume is fed with heat (or
     light fluid), not momentum: the flux of buoyancy through every cross-section equals the source's. · **plain** What a
     plume conserves is its buoyancy flux.
  9. **did** Plumes: let buoyancy increase the momentum flux · **tex** $\dfrac{d}{dx}\big(U_s^2\,\delta\big)\sim g'\,\delta\ \Rightarrow\ 2n+m-1=b+m\ \Rightarrow\ b=2n-1$ · **why** The stream-wise (vertical)
     momentum equation gains the force g′ per unit mass, integrated across the plume; match powers of x. · **plain** The
     momentum flux of a plume grows with height.
  10. **did** Plane plume: solve · **tex** $m=1,\ b=2n-1,\ b+n+1=0\ \Rightarrow\ n=0,\ b=-1$ · **why** Three linear equations (steps 1, 8, 9). The scalar (or
      temperature excess) scales like g′. · **plain** A plane plume rises at constant speed while its excess temperature
      falls as 1/x.
  11. **did** Round plume: solve · **tex** $b+n+2m=0,\ b=2n-1,\ m=1\ \Rightarrow\ n=-\tfrac13,\ b=-\tfrac53$ · **why** Disc area in both the buoyancy flux (δ²) and the
      momentum balance ($\frac{d}{dx}(U_s^2\delta^2)\sim g'\delta^2$ gives the same b = 2n − 1). · **plain** A round plume: speed ∝ x^{−1/3},
      excess ∝ x^{−5/3}.
  12. **did** Shear layer, and collect · **tex** shear layer: $U_s=U_1-U_2=\text{const}\ \Rightarrow\ n=0,\ m=1$ · **why** The velocity difference is imposed by the
      two streams; only the width grows. All rows: `free_shear_exponents(flow)`. · **plain** A mixing layer simply widens
      linearly between two fixed speeds.
- **Result** (width m, velocity n, scalar, local Reynolds number n + m): plane jet (1, −½, −½, +½) · round jet (1, −1, −1, 0) ·
  plane wake (½, −½, −½, 0) · round wake (⅓, −⅔, −⅔, −⅓) · plane plume (1, 0, −1, +1) · round plume (1, −⅓, −5/3, +⅔) · shear layer
  (1, 0, –, +1) — "an invariant and a growth law give every row".
- **Check:** plane jet agrees with $U=C_5(J_s/\rho)^{1/2}x^{-1/2}F(y/x)$ (12.66) and $\bar Y\propto x^{-1/2}$ (12.71) ✓; similarity variables of
  the wakes: $y/\sqrt{\theta_px}$ and $r/(\theta_r^2x)^{1/3}$ ✓ (m = ½, ⅓), where $\theta_p$ = (drag per unit span)/($\rho U_\infty^2$) is a length and $\theta_r^2$ = drag/($\rho U_\infty^2$) an area — the momentum thicknesses of the two wakes. Units of each invariant: N/m (plane) or N (round). Code:
  `free_shear_exponents` returns these as exact `Fraction`s; `wrong_exponent_fluxes` is flat only when 2n + m = 0.
- **What it means:** the table is not seven facts to memorise but one recipe; the *constants* (spreading rates, amplitude
  factors) are measurements and are not reproduced here. **Fails when:** the invariant is not conserved (co-flowing jets
  with pressure gradients, stratified surroundings — a plume in a stable atmosphere stops rising).
- **Traps:** using the jet's growth law for a wake; forgetting δ² for round flows; confusing which quantity is conserved in
  a plume (buoyancy flux, not momentum flux).

### D17 · The linear total stress $\bar\tau=\tau_0(1-2y/h)$ from (12.76)–(12.77), and $dP/dx=-2\tau_0/h$ (12.90), $-4\tau_0/d$ (12.91) — ★★, 9 steps, in C10 (notebook · `law_of_the_wall`)
- **Goal:** show that in a fully developed channel the total shear stress is a straight line across the channel, whatever
  the turbulence does. · **Start:** the constant-density form of the mean momentum equation (12.30), $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=\frac1\rho\frac{\partial\bar\tau_{ij}}{\partial x_j}$ with $\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}-\rho\overline{u_iu_j}$, for flow between walls at y = 0 and
  y = h. · **Plan:** • use "fully developed" to delete terms • the wall-normal equation fixes how pressure varies with y • a
  function of y equal to a function of x is a constant • integrate and use symmetry. · **Tools:** fully developed flow
  (ch08); separation argument (P167); integration of a constant; control-volume balance (ch08 D07). · **Assumptions:** steady
  in the mean; statistics independent of x and z; smooth identical walls; constant density.
- **Steps:**
  1. **did** Apply "fully developed" · **tex** $\dfrac{\partial}{\partial x}\big(\text{any mean velocity or Reynolds stress}\big)=0,\qquad V=0$ · **why** Far from the inlet the statistics no longer change
     downstream; then continuity gives $\partial V/\partial y=0$, and V = 0 at the wall, so V = 0 everywhere. Only P may depend on x. ·
     **plain** A mean flow U(y) only.
  2. **did** Reduce the stream-wise equation · **tex** $0=-\dfrac{\partial P}{\partial x}+\dfrac{\partial\bar\tau}{\partial y},\qquad\bar\tau=\mu\dfrac{\partial U}{\partial y}-\rho_0\overline{uv}$ · **why** The left side of the Start, $\partial U_i/\partial t+U_j\,\partial U_i/\partial x_j$, vanishes (steady, V = 0,
     ∂U/∂x = 0); of the stress divergence only the y-derivative of the xy-component remains. This is the first of (12.76), $0=-\frac{\partial P}{\partial x}+\frac{\partial\bar\tau}{\partial y}$; the definition of $\bar\tau$ is the unnumbered line that follows it in the book. · **plain**
     Pressure gradient balances the gradient of the total shear stress.
  3. **did** Reduce the wall-normal equation · **tex** $0=-\dfrac{\partial}{\partial y}\big(P+\rho\,\overline{v^2}\big)$ · **why** Same reduction for the y-component: only pressure and the normal
     Reynolds stress survive. This is the second of (12.76), $0=-\frac{\partial}{\partial y}\big(P+\rho\overline{v^2}\big)$. · **plain** Across the channel, P + ρ⟨v²⟩ is uniform.
  4. **did** Integrate from the wall · **tex** $P(x,y)-P(x,0)=-\rho\,\overline{v^2}+\rho\big[\overline{v^2}\big]_{y=0}=-\rho\,\overline{v^2}$ · **why** Fundamental theorem of calculus; the
     fluctuation vanishes at the wall (no-slip, impermeable). · **plain** The pressure inside is a little lower than at the
     wall.
  5. **did** Differentiate in x · **tex** $\dfrac{\partial}{\partial x}P(x,y)-\dfrac{d}{dx}P(x,0)=-\rho\dfrac{\partial}{\partial x}\overline{v^2}=0$ · **why** Step 1: no Reynolds stress depends on x. This is (12.77). ·
     **plain** The stream-wise pressure gradient is the same at every height: a function of x only.
  6. **did** Separate · **tex** $\dfrac{d\bar\tau}{dy}(y)=\dfrac{dP}{dx}(x)=\text{constant}$ · **why** In step 2 the left side depends on y only (step 1), the right on x
     only (step 5); they can be equal for all x and y only if both are one constant (P167). · **plain** The stress gradient
     is uniform.
  7. **did** Integrate in y · **tex** $\bar\tau(y)=\tau_0+\dfrac{dP}{dx}\,y$ · **why** A constant derivative gives a straight line; at the wall the Reynolds stress is
     zero and $\bar\tau(0)=\tau_0$, the wall shear stress. · **plain** The total stress is linear in y.
  8. **did** Use symmetry at the other wall · **tex** $\bar\tau(h)=-\tau_0\ \Rightarrow\ \dfrac{dP}{dx}=-\dfrac{2\tau_0}{h}$ · **why** The flow is symmetric about the centreline, so the
     stress on the upper wall is equal and opposite. This is (12.90), h = full height. · **plain** The pressure drop pays
     for the friction on both walls.
  9. **did** Substitute back; do the pipe by a force balance · **tex** $\bar\tau(y)=\tau_0\Big(1-\dfrac{2y}{h}\Big);\qquad\dfrac{\pi d^2}{4}\,dP=-\tau_0\,\pi d\,dx\ \Rightarrow\ \dfrac{dP}{dx}=-\dfrac{4\tau_0}{d}$ · **why** Steps 7 and 8; for a pipe
     the pressure force on a slug's ends balances friction on its side (ch08 D07). The second is (12.91). · **plain** Zero
     stress on the centreline, τ₀ at the wall, a straight line between.
- **Result:** $\bar\tau(y)=\mu\,dU/dy-\rho\overline{uv}=\tau_0(1-2y/h)$; $dP/dx=-2\tau_0/h$ (12.90); pipe $dP/dx=-4\tau_0/d$ (12.91) — "the total stress is linear;
  only its split between viscous and Reynolds parts depends on the turbulence".
- **Check:** units — Pa/m both sides. Laminar limit: $\overline{uv}=0$, μ dU/dy linear ⇒ the parabola of Ch. 8 ✓. Number: τ₀ = 0.3 Pa, h
  = 0.1 m ⇒ dP/dx = −6 Pa/m. Code: `channel_momentum_residual` ≈ 0; pipe parity with `LAM.pipe_wall_stress`.
- **What it means:** very near the wall (y ≪ h) the stress is nearly the constant τ₀ — the "constant-stress layer" on which
  the law of the wall (D18) and the mixing-length model (D21) stand. **Fails when:** the flow is still developing, or the
  walls differ (the zero-stress point moves off the centreline).
- **Traps:** using the half-height in $dP/dx=-2\tau_0/h$ (12.90) (h is the full height); thinking ⟨v²⟩ changes dP/dx; citing the proof as
  Exercise 12.31 (slip #7: it is 12.32).

### D18 · The law of the wall (12.80) with $u_*$ (12.81), and the viscous sublayer $U^+=y^+$ (12.82) — ★, 7 steps, in C10 (notebook · `law_of_the_wall`)
- **Goal:** find the variables in which the mean velocity near any smooth wall is one universal curve. · **Start:**
  $U=U(\rho,\tau_0,\nu,y)$ (12.79). · **Plan:** • count dimensionless groups • build a velocity and a length from the wall
  quantities • very near the wall integrate the stress balance. · **Tools:** Π theorem (ch01 C69, `DIM.pi_groups`);
  integration of a constant. · **Assumptions:** smooth wall; close enough to it that δ and $U_\infty$ do not matter; constant total
  stress (y ≪ h, D17).
- **Steps:**
  1. **did** List what the near-wall flow depends on · **tex** $U=U(\rho,\ \tau_0,\ \nu,\ y)$ · **why** The wall acts on the fluid only through its stress τ₀;
     the fluid's properties are ρ and ν; y is where we look. The outer thickness δ is deliberately left out. This is (12.79).
     · **plain** Five quantities.
  2. **did** Count the groups · **tex** $5\ \text{quantities}-3\ \text{dimensions (kg, m, s)}=2\ \text{groups}$ · **why** Buckingham's Π theorem (ch01). · **plain** One
     dimensionless velocity is a function of one dimensionless distance.
  3. **did** Build a velocity from the wall stress · **tex** $u_*\equiv\sqrt{\tau_0/\rho}$ · **why** τ₀/ρ has units (N/m²)/(kg/m³) = m²/s²; it is the only
     velocity that can be made from τ₀ and ρ. This is (12.81), $u_*^2\equiv\tau_0/\rho$. · **plain** The friction velocity: a measure of
     the wall stress in m/s.
  4. **did** Build a length and form the groups · **tex** $l_\nu=\dfrac{\nu}{u_*};\qquad U^+\equiv\dfrac{U}{u_*}=f\Big(\dfrac{y\,u_*}{\nu}\Big)=f(y^+)$ · **why** ν/u_* is m²s⁻¹/(m s⁻¹) = m; dividing U and y by
     these scales gives the two groups. This is (12.80). · **plain** In wall units every smooth-wall flow has the same
     profile near the wall.
  5. **did** Very near the wall keep only the viscous stress · **tex** $\mu\dfrac{dU}{dy}=\tau_0$ · **why** At the wall the fluctuations vanish, so $\overline{uv}\to0$; and the
     total stress is still ≈ τ₀ because y ≪ h (D17). · **plain** In the sublayer the stress is carried by viscosity alone.
  6. **did** Integrate from the wall · **tex** $U=\dfrac{\tau_0\,y}{\mu}$ · **why** τ₀/μ is constant; U = 0 at y = 0 (no slip). · **plain** A linear profile.
  7. **did** Write it in wall units · **tex** $\dfrac{U}{u_*}=\dfrac{\rho u_*^2\,y}{\mu\,u_*}=\dfrac{y\,u_*}{\nu}\ \Rightarrow\ U^+=y^+$ · **why** τ₀ = ρu_*² and μ/ρ = ν. This is (12.82), good to about y⁺ ≈ 5.
     · **plain** In the viscous sublayer the function f is simply the identity.
- **Result:** $U^+=f(y^+)$ (12.80), $u_*^2=\tau_0/\rho$ (12.81); $U^+=y^+$ for $y^+\lesssim5$ (12.82) — "the wall stress alone sets the scales
  near a wall".
- **Check:** units — $u_*$ m/s, $l_\nu$ m. Number: τ₀ = 0.3 Pa in air ⇒ u_* = 0.5 m/s, l_ν = 30 µm; at y = 0.09 mm, y⁺ = 3, U = 1.5 m/s.
  Spalding's profile differs from U⁺ = y⁺ by 0.02 % at y⁺ = 1 and 2.7 % at y⁺ = 5.
- **What it means:** profiles from pipes, channels and boundary layers at any Reynolds number collapse near the wall when
  scaled this way; C11 extends f beyond the sublayer. **Fails when:** the wall is rough (elements taller than ~ 5 l_ν),
  or the pressure gradient is strong enough that the stress is not constant across the inner layer.
- **Traps:** including δ in the inner list; reading U⁺ = y⁺ on a semi-log plot as "curved, so not linear"; confusing u_*
  with a velocity that exists somewhere in the flow.

### D19 · The defect law (12.84), overlap matching (12.85)–(12.87) and the logarithmic law (12.88)–(12.89); friction law $U_\infty^+=\tfrac1\kappa\ln\delta^++A+B$ — ★★, 12 steps, in C11 (notebook · `law_of_the_wall`)
- **Goal:** find the velocity profile in the region where both the near-wall and the outer descriptions hold. · **Start:**
  the law of the wall $U/u_*=f(y^+)$ (12.80) and the outer list $U=U(\rho,\tau_0,\delta,y)$ (12.83). · **Plan:** • build the outer law by
  dimensional analysis • differentiate both laws • equate the gradients where both hold • separate variables • integrate. ·
  **Tools:** Π theorem (ch01); chain rule (P49); matching in the overlap (P300, primer in C10); separation argument (P167);
  $\int dy/y=\ln y$. · **Assumptions:** an overlap region exists: $y^+\gg1$ and $y/\delta\ll1$ at once, i.e. $\delta^+\gg1$.
- **Steps:**
  1. **did** List what the outer flow depends on · **tex** $U=U(\rho,\ \tau_0,\ \delta,\ y)$ · **why** Far from the wall the Reynolds stress carries the momentum
     and viscosity does not matter; the wall is still felt through τ₀, and the layer's size δ enters. This is (12.83). ·
     **plain** Outer flow: no ν, but δ.
  2. **did** Form the outer law for the velocity defect · **tex** $\dfrac{U_\infty-U}{u_*}=F_d\Big(\dfrac{y}{\delta}\Big)=F_d(\xi)$ · **why** Π theorem: one group y/δ. The *defect*
     is used because without ν the outer flow knows U only up to an additive constant (it cannot see the no-slip wall).
     This is (12.84) (the book's F). · **plain** How far the velocity falls short of the edge value, in units of u_*.
  3. **did** Differentiate the inner law · **tex** $\dfrac{dU}{dy}=u_*\dfrac{df}{dy^+}\dfrac{dy^+}{dy}=\dfrac{u_*^2}{\nu}\dfrac{df}{dy^+}$ · **why** Chain rule with $y^+=yu_*/\nu$. This is (12.85). · **plain** The mean
     shear from the inner description.
  4. **did** Differentiate the outer law · **tex** $-\dfrac{dU}{dy}=\dfrac{u_*}{\delta}\dfrac{dF_d}{d\xi}$ · **why** Chain rule with ξ = y/δ; $U_\infty$ is constant. This is (12.86). · **plain**
     The same shear from the outer description.
  5. **did** Equate them in the overlap · **tex** $\dfrac{u_*^2}{\nu}\dfrac{df}{dy^+}=-\dfrac{u_*}{\delta}\dfrac{dF_d}{d\xi}$ · **why** Where both laws are valid they describe the same flow, so the
     shear must agree (matching, P300). Gradients are matched, not velocities, because the two laws differ by a constant.
     · **plain** Two formulas for one slope.
  6. **did** Multiply by $y/u_*$ · **tex** $y^+\dfrac{df}{dy^+}=-\xi\dfrac{dF_d}{d\xi}$ · **why** $y\,u_*/\nu=y^+$ on the left and $y/\delta=\xi$ on the right: each side now
     contains only its own variable. This is (12.87). · **plain** A function of y⁺ alone equals a function of ξ alone.
  7. **did** Conclude that both sides are constant · **tex** $y^+\dfrac{df}{dy^+}=-\xi\dfrac{dF_d}{d\xi}=\dfrac1\kappa$ · **why** y⁺ and ξ can be varied independently (change δ at fixed
     y⁺), so neither side can depend on its variable (P167). κ is the von Kármán constant (≈ 0.4; not the thermal
     diffusivity). · **plain** In the overlap, y dU/dy = u_*/κ — a number, the same at every height.
  8. **did** Write the inner equation as an ODE · **tex** $\dfrac{df}{dy^+}=\dfrac{1}{\kappa\,y^+}$ · **why** Divide step 7 by y⁺. · **plain** The slope falls off as 1/y.
  9. **did** Integrate it · **tex** $U^+=f(y^+)=\dfrac1\kappa\ln(y^+)+B$ · **why** $\int dy^+/y^+=\ln y^+$; B is the integration constant, fixed by experiment
     (it remembers the sublayer and buffer layer). This is (12.88). · **plain** The logarithmic law.
  10. **did** Integrate the outer equation the same way · **tex** $F_d(\xi)=-\dfrac1\kappa\ln(\xi)+A$ · **why** $dF_d/d\xi=-1/(\kappa\xi)$. This is (12.89); A remembers
      the outer flow (the wake). · **plain** The same logarithm, seen from outside.
  11. **did** Add the two forms at one height · **tex** $\dfrac{U_\infty}{u_*}=\dfrac1\kappa\ln y^+-\dfrac1\kappa\ln\xi+A+B=\dfrac1\kappa\ln\delta^++A+B$ · **why** $U/u_*+(U_\infty-U)/u_*=U_\infty/u_*$ and
      $\ln y^+-\ln(y/\delta)=\ln(\delta u_*/\nu)$: y drops out. (Ours.) · **plain** A friction law: the edge velocity in wall units grows
      as the logarithm of the Reynolds number.
  12. **did** Turn it into a skin-friction coefficient · **tex** $C_f\equiv\dfrac{\tau_0}{\tfrac12\rho U_\infty^2}=\dfrac{2}{(U_\infty^+)^2}=\dfrac{2}{\big[\tfrac1\kappa\ln\delta^++A+B\big]^2}$ · **why** τ₀ = ρu_*². · **plain** Turbulent skin
      friction falls only logarithmically with Reynolds number.
- **Result:** $U^+=\tfrac1\kappa\ln y^++B$ (12.88); $F_d(\xi)=-\tfrac1\kappa\ln\xi+A$ (12.89); $U_\infty^+=\tfrac1\kappa\ln\delta^++A+B$ — "the logarithm is what remains
  when the shear may depend on neither length".
- **Check:** units — all dimensionless. Illustrative κ = 0.41, B = 5.0: U⁺(100) = 16.23, U⁺(1000) = 21.85 (+5.62 per decade).
  With A = 1.0, δ⁺ = 1000: U_∞⁺ = 22.85, C_f = 3.8 × 10⁻³. `overlap_matching_sympy()` reproduces steps 3–10. The same argument in
  dimensional form: dU/dy can depend only on u_* and y ⇒ dU/dy = u_*/(κy).
- **What it means:** wall functions in CFD, the roughness-length wind profile and bulk drag laws are all this logarithm.
  **Fails when:** δ⁺ ≲ 200 (no overlap), near separation (u_* → 0), or with strong pressure gradients.
- **Traps:** matching velocities instead of gradients; the sign in the defect law; ln vs log₁₀ (slope per decade is
  2.303/κ); treating κ and B as exact universal numbers (they vary a little with flow type — required keywords in code).

### D20 · The rough-wall law $U^+=\tfrac1\kappa\ln(y/y_0)$ (12.93), the equivalent $B=-\tfrac1\kappa\ln y_0^+$ and the neutral drag coefficient — ★, 5 steps, in C11 (notebook · `law_of_the_wall`)
- **Goal:** write the logarithmic law for a rough surface, where there is no viscous sublayer to anchor it. · **Start:**
  $dU/dy=u_*/(\kappa y)$ (D19 step 7 in dimensional form). · **Plan:** • integrate without using ν • name the height where the
  line extrapolates to zero • compare with the smooth law • form a drag coefficient. · **Tools:** logarithm rules;
  definition of a drag coefficient. · **Assumptions:** fully rough (elements stick out of the sublayer, so viscosity is
  irrelevant); y well above the elements.
- **Steps:**
  1. **did** Integrate the overlap shear in dimensional form · **tex** $\dfrac{U}{u_*}=\dfrac1\kappa\ln y+C$ · **why** Same integral as D19 step 9, but we may not
     use $l_\nu$ to make y dimensionless: on a rough wall viscosity has dropped out. · **plain** Still a logarithm; its constant
     must come from the roughness.
  2. **did** Define the roughness length · **tex** $U=0\ \text{at}\ y=y_0\ \Rightarrow\ C=-\dfrac1\kappa\ln y_0$ · **why** The constant is fixed by the height at which the straight
     line (extrapolated downward) crosses U = 0 — a property of the surface. · **plain** $y_0$ packages everything about the
     roughness into one length.
  3. **did** Combine · **tex** $U^+=\dfrac{U}{u_*}=\dfrac1\kappa\ln\Big(\dfrac{y}{y_0}\Big)$ · **why** $\ln y-\ln y_0=\ln(y/y_0)$. This is (12.93). · **plain** The rough-wall logarithmic law.
  4. **did** Compare with the smooth form · **tex** $\dfrac1\kappa\ln y^++B=\dfrac1\kappa\ln\dfrac{y}{y_0}\ \Rightarrow\ B=-\dfrac1\kappa\ln y_0^+,\quad y_0^+\equiv\dfrac{y_0u_*}{\nu}$ · **why** Write $\ln(y/y_0)=\ln y^+-\ln y_0^+$. (Ours.) A smooth
     wall with B = 5.0 behaves like $y_0^+=e^{-\kappa B}=0.13$. · **plain** Roughness lowers the line: a larger $y_0$ is a smaller B.
  5. **did** Form the drag coefficient at a reference height · **tex** $C_D\equiv\dfrac{\tau_0}{\rho\,U(z)^2}=\dfrac{u_*^2}{U(z)^2}=\Big[\dfrac{\kappa}{\ln(z/z_0)}\Big]^2$ · **why** τ₀ = ρu_*² and step 3 with y → z,
     $y_0\to z_0$ (the meteorological names). · **plain** The drag coefficient depends only on the measuring height and the
     roughness — not on the Reynolds number.
- **Result:** $U^+=\tfrac1\kappa\ln(y/y_0)$ (12.93); $C_D=[\kappa/\ln(z/z_0)]^2$ — "over a rough surface the wind is logarithmic in height
  with zero at $y_0$".
- **Check:** units — dimensionless. Number: z₀ = 0.03 m, u_* = 0.4 m/s, κ = 0.41: U(10 m) = 5.67 m/s; C_D(10 m) = 4.98 × 10⁻³.
  Sea (z₀ = 2 × 10⁻⁴ m): C_D = 1.4 × 10⁻³.
- **What it means:** this is the neutral wind profile of the atmospheric surface layer; C15 adds the stability correction.
  **Fails when:** y is within a few element heights (use a displacement height), or the roughness is transitional.
- **Traps:** thinking $y_0$ is the height of the roughness elements (it is typically a tenth to a thirtieth of it);
  evaluating the law below $y_0$; forgetting that C_D depends on the reference height.

### D21 · The mixing length: $-\overline{uv}=l_T^2(dU/dy)^2$, the wall model (12.100), $dU^+/dy^+=2/(1+\sqrt{1+4\kappa^2y^{+2}})$ and the logarithmic limit (12.101) — ★★, 11 steps, in C12 (notebook · `mixing_length_closure`)
- **Goal:** build the simplest closure for the Reynolds shear stress and show that, near a wall, it returns the logarithmic
  law. · **Start:** $\overline{u_iu_j}=\dfrac23\bar e\,\delta_{ij}-\nu_T\Big(\dfrac{\partial U_i}{\partial x_j}+\dfrac{\partial U_j}{\partial x_i}\Big)$ (12.94) and $0=-\dfrac1\rho\dfrac{dP}{dx}+\dfrac{\partial}{\partial y}\Big(\nu\dfrac{\partial U}{\partial y}-\overline{uv}\Big)$ (12.99). · **Plan:** • specialise the
  hypothesis to a simple shear • estimate ν_T from a length and a velocity • put l_T = κy into the momentum equation •
  integrate once • solve the quadratic • take the two limits. · **Tools:** the parcel estimate of C04 (N50); quadratic
  formula (P159); limits; `cumulative_trapezoid` (P190). · **Assumptions:** unidirectional mean flow U(y); zero pressure
  gradient near the wall (constant-stress layer, D17); the eddy size is proportional to the distance from the wall.
- **Steps:**
  1. **did** Specialise the hypothesis to U(y) · **tex** $-\overline{uv}=\nu_T\dfrac{dU}{dy}$ · **why** In $\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ (12.94) take i = 1, j = 2: $\delta_{12}=0$ and only $\partial U_1/\partial x_2=dU/dy$
     is non-zero. · **plain** The Reynolds shear stress is modelled like a viscous stress with a different viscosity.
  2. **did** Estimate the eddy viscosity · **tex** $\nu_T\sim l_T\,u_T$ · **why** A diffusivity is a length times a velocity (units m²/s); for turbulence
     the carriers are eddies of size $l_T$ and speed $u_T$. This is (12.98). · **plain** Big, fast eddies mix more.
  3. **did** Estimate the eddy velocity from the shear · **tex** $u_T\sim l_T\dfrac{dU}{dy}\ \Rightarrow\ -\overline{uv}=l_T^2\Big(\dfrac{dU}{dy}\Big)^2$ · **why** A parcel displaced by $l_T$ arrives with a
     velocity excess $l_T\,dU/dy$ (the parcel argument of C04); the proportionality constant is absorbed into $l_T$. · **plain**
     The stress grows as the square of the shear.
  4. **did** Choose the length near a wall and insert · **tex** $l_T=\kappa y:\qquad0=-\dfrac1\rho\dfrac{dP}{dx}+\dfrac{\partial}{\partial y}\Big(\nu\dfrac{dU}{dy}+\kappa^2y^2\Big(\dfrac{dU}{dy}\Big)^2\Big)$ · **why** An eddy cannot be larger than its
     distance from the wall, so $l_T\propto y$; put step 3 into $0=-\frac1\rho\frac{dP}{dx}+\frac{\partial}{\partial y}\big(\nu\frac{\partial U}{\partial y}-\overline{uv}\big)$ (12.99). This is (12.100). · **plain** One equation for U(y)
     with one constant, κ.
  5. **did** Integrate once at zero pressure gradient · **tex** $\nu\dfrac{dU}{dy}+\kappa^2y^2\Big(\dfrac{dU}{dy}\Big)^2=\text{const}$ · **why** With dP/dx = 0 the bracket has zero derivative,
     so it is constant. · **plain** The total stress — viscous plus turbulent — is the same at every height.
  6. **did** Evaluate the constant at the wall · **tex** $\text{const}=\nu\dfrac{dU}{dy}\Big|_{y=0}=\dfrac{\tau_0}{\rho}=u_*^2$ · **why** At y = 0 the mixing-length term vanishes and the
     viscous stress is the wall stress; $u_*^2\equiv\tau_0/\rho$ (12.81). · **plain** The constant is the wall stress.
  7. **did** Write it in wall units · **tex** $s+\kappa^2y^{+2}s^2=1,\qquad s\equiv\dfrac{dU^+}{dy^+}$ · **why** Divide by $u_*^2$: $\nu(dU/dy)/u_*^2=dU^+/dy^+$ and
     $\kappa^2y^2(dU/dy)^2/u_*^2=\kappa^2y^{+2}s^2$. · **plain** Viscous share + turbulent share = 1.
  8. **did** Solve the quadratic for the slope · **tex** $s=\dfrac{-1+\sqrt{1+4\kappa^2y^{+2}}}{2\kappa^2y^{+2}}=\dfrac{2}{1+\sqrt{1+4\kappa^2y^{+2}}}$ · **why** Quadratic formula with the positive root (the
     velocity increases away from the wall); multiply top and bottom by $1+\sqrt{\ }$ to avoid 0/0 at the wall. (Ours.) · **plain**
     The exact slope of the model profile at every height.
  9. **did** Take the limit near the wall · **tex** $y^+\to0:\ s\to1\ \Rightarrow\ U^+=y^+$ · **why** The square root → 1. · **plain** The model contains the viscous
     sublayer $U^+=y^+$ (12.82).
  10. **did** Take the limit far from the wall · **tex** $y^+\gg1:\ s\to\dfrac{1}{\kappa y^+}\ \Rightarrow\ \dfrac{dU}{dy}\cong\sqrt{\dfrac{\tau_0}{\rho}}\,\dfrac{1}{\kappa y},\quad\dfrac{U}{u_*}\cong\dfrac1\kappa\ln y+\text{const.}$ · **why** $\sqrt{1+4\kappa^2y^{+2}}\approx2\kappa y^+\gg1$;
      equivalently, drop the viscous term in step 5. This is (12.101). · **plain** Outside the sublayer the model gives the
      logarithmic law with slope 1/κ.
  11. **did** Integrate the exact slope · **tex** $U^+=\dfrac1\kappa\Big[\sinh^{-1}(a)-\dfrac{\sqrt{1+a^2}-1}{a}\Big],\ a=2\kappa y^+;\qquad B=\lim\Big(U^+-\dfrac1\kappa\ln y^+\Big)=\dfrac{\ln(4\kappa)-1}{\kappa}$ · **why**
      $\int da/(1+\sqrt{1+a^2})=\sinh^{-1}a-(\sqrt{1+a^2}-1)/a$; for large a, $\sinh^{-1}a\to\ln(2a)$. (Ours.) · **plain** The model also predicts
      the intercept — and gets it badly wrong: −1.23 for κ = 0.41 instead of about 5.
- **Result:** $-\overline{uv}=\kappa^2y^2(dU/dy)^2$; $dU^+/dy^+=2/(1+\sqrt{1+4\kappa^2y^{+2}})$; $U/u_*\cong\tfrac1\kappa\ln y+\text{const.}$ (12.101) — "one guessed length returns
  the log law's slope".
- **Check:** units — $l_T^2(dU/dy)^2$: m² s⁻² ✓. y⁺ = 10, κ = 0.41: s = 0.216 (log slope 0.244); y⁺ = 100: 0.0241 (0.0244).
  Intercept: [ln(1.64) − 1]/0.41 = −1.23 (numerical quadrature −1.232); with van Driest damping $l_T=\kappa y[1-e^{-y^+/26}]$: B =
  5.28.
- **What it means:** κ controls the slope only. The plain model lets eddies mix right down to the wall, so the sublayer is
  too thin and the line sits about 6 units low; damping the mixing length within ~ 26 wall units fixes the intercept. Two
  constants, two jobs. **Fails when:** no single length is obvious (separated flows, jets in cross-flow) — then C13.
- **Traps:** losing the sign for negative shear (use $l_T^2\lvert dU/dy\rvert\,dU/dy$); taking the negative root; believing that
  recovering the log law validates the model (the slope was built in by l_T ∝ y).

### D22 · The modelled energy equation (12.103) and the k–ε eddy viscosity $\nu_T=C_\mu\bar e^2/\bar\varepsilon$ (12.104) — ★, 6 steps, in C13 (notebook)
- **Goal:** turn the exact (but unclosed) energy budget into an equation a computer can solve, and build an eddy viscosity
  from its two variables. · **Start:** $\dfrac{\partial\bar e}{\partial t}+U_j\dfrac{\partial\bar e}{\partial x_j}=\dfrac{\partial}{\partial x_j}\Big(-\dfrac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\tfrac12\overline{u_i^2u_j}\Big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47) and
  $\nu_T\sim l_Tu_T$ (12.98). · **Plan:** • drop buoyancy (constant density) • replace the unknown transport by a gradient flux
  • name the dissipation • build scales from ē and ε̄. · **Tools:** gradient-diffusion hypothesis (N162); substitution; units
  check. · **Assumptions:** constant density; high Reynolds number; transport of ē behaves like diffusion.
- **Steps:**
  1. **did** Drop the buoyancy term · **tex** $g\alpha\,\overline{u_3T'}\to0$ · **why** From §12.8 to §12.10 the density is constant (the switch announced at the top of §12.8). · **plain** No
     buoyant production in this model.
  2. **did** Model the three transport terms as one gradient flux · **tex** $-\dfrac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\tfrac12\overline{u_i^2u_j}=\dfrac{\nu_T}{\sigma_e}\dfrac{\partial\bar e}{\partial x_j}$ · **why** They only move energy from
     where there is much to where there is little; the simplest such law is diffusion down the gradient, with a
     diffusivity $\nu_T/\sigma_e$ ($\sigma_e$ = a constant of order one). A model, not a theorem. · **plain** Turbulent energy
     diffuses like a scalar.
  3. **did** Name the dissipation and keep the production exact · **tex** $\dfrac{\partial\bar e}{\partial t}+U_j\dfrac{\partial\bar e}{\partial x_j}=\dfrac{\partial}{\partial x_j}\Big(\dfrac{\nu_T}{\sigma_e}\dfrac{\partial\bar e}{\partial x_j}\Big)-\bar\varepsilon-\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}$ · **why** $\bar\varepsilon\equiv2\nu\overline{S'_{ij}S'_{ij}}$ becomes an
     unknown of its own; the production needs no model once $\overline{u_iu_j}$ is given by the closure $\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ (12.94). This is (12.103). · **plain**
     A transport equation for ē with one diffusion term, one sink, one source.
  4. **did** Build scales from ē and ε̄ · **tex** $u_T=\sqrt{\bar e},\qquad t_T=\dfrac{\bar e}{\bar\varepsilon},\qquad l_T=u_T\,t_T=\dfrac{\bar e^{3/2}}{\bar\varepsilon}$ · **why** ē [m²/s²] gives a velocity; ē/ε̄ [s] is the
     time in which dissipation would use up the energy; their product is a length. · **plain** The eddies' speed, lifetime
     and size from the two model variables.
  5. **did** Form the eddy viscosity · **tex** $\nu_T=C_\mu\,l_T\,u_T=C_\mu\Big[\dfrac{\bar e^{3/2}}{\bar\varepsilon}\Big]\sqrt{\bar e}=C_\mu\dfrac{\bar e^{\,2}}{\bar\varepsilon}$ · **why** $\nu_T\sim l_Tu_T$ (12.98) with step 4; $C_\mu$ is the proportionality
     constant. This is (12.104). · **plain** More energy, or slower dissipation, means more mixing.
  6. **did** Write the production with the model stress · **tex** $-\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}=2\nu_T\,\bar S_{ij}\bar S_{ij}\ \ge0$ · **why** Insert $\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ (12.94): the $\tfrac23\bar e\delta_{ij}$ part
     contracts to $\partial U_i/\partial x_i=0$ (12.27); the rest is a symmetric contraction. · **plain** In this model production can never
     be negative.
- **Result:** $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\frac{\nu_T}{\sigma_e}\frac{\partial\bar e}{\partial x_j}\big)-\bar\varepsilon-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}$ (12.103) with $\nu_T=C_\mu\bar e^2/\bar\varepsilon$ (12.104) — "two fields, ē and ε̄, give a viscosity everywhere".
- **Check:** units — ē²/ε̄: (m²/s²)²/(m²/s³) = m²/s ✓. Number: ē = 1 m²/s², ε̄ = 1 m²/s³, $C_\mu$ = 0.09 ⇒ ν_T = 0.09 m²/s, l_T = 1 m.
- **What it means:** the model still needs an equation for ε̄ — $\frac{\partial\bar\varepsilon}{\partial t}+U_j\frac{\partial\bar\varepsilon}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\frac{\nu_T}{\sigma_\varepsilon}\frac{\partial\bar\varepsilon}{\partial x_j}\big)-C_{\varepsilon1}\big(\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}\big)\frac{\bar\varepsilon}{\bar e}-C_{\varepsilon2}\frac{\bar\varepsilon^2}{\bar e}$ (12.105), built by analogy rather than derived. **Fails
  when:** transport is not gradient-like (counter-gradient fluxes in convection), or the stress is not aligned with the
  mean strain (strong rotation, curvature).
- **Traps:** thinking the Result is exact (two of its terms are models); writing the book's printed $u_j$ in the viscous
  transport (slip #8: it is $2\nu\overline{u_iS'_{ij}}$, the viscous transport of the exact budget in the Start).

### D23 · What fixes the k–ε constants: decay $\bar e\propto(t+t_0)^{-n}$, $n=1/(C_{\varepsilon2}-1)$, and the log layer $\kappa^2=\sqrt{C_\mu}(C_{\varepsilon2}-C_{\varepsilon1})\sigma_\varepsilon$ — ★★, 12 steps, in C13 (notebook)
- **Goal:** show that two of the model's constants are tied to things one can measure — how fast grid turbulence decays,
  and the von Kármán constant. The book states the constants; this derivation is ours. · **Start:** $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(\frac{\nu_T}{\sigma_e}\frac{\partial\bar e}{\partial x_j}\big)-\bar\varepsilon-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}$ (12.103) and $\dfrac{\partial\bar\varepsilon}{\partial t}+U_j\dfrac{\partial\bar\varepsilon}{\partial x_j}=\dfrac{\partial}{\partial x_j}\Big(\dfrac{\nu_T}{\sigma_\varepsilon}\dfrac{\partial\bar\varepsilon}{\partial x_j}\Big)-C_{\varepsilon1}\Big(\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}\Big)\dfrac{\bar\varepsilon}{\bar e}-C_{\varepsilon2}\dfrac{\bar\varepsilon^{\,2}}{\bar e}$ (12.105). · **Plan:** • homogeneous decay: two ODEs, divide them, integrate
  twice • log layer: production = dissipation, uniform ē, and what is left of the ε̄ equation. · **Tools:** dividing two
  ODEs (P303, primer above); separable ODE; the log-layer values of C11. · **Assumptions:** (a) homogeneous turbulence with
  no mean shear; (b) a steady constant-stress layer in local equilibrium.
- **Steps:**
  1. **did** Reduce the ē equation for decaying homogeneous turbulence · **tex** $\dfrac{d\bar e}{dt}=-\bar\varepsilon$ · **why** Homogeneous: every spatial derivative of a
     statistic vanishes (no transport, no advection); no mean shear: no production. · **plain** The energy only decays.
  2. **did** Reduce the ε̄ equation the same way · **tex** $\dfrac{d\bar\varepsilon}{dt}=-C_{\varepsilon2}\dfrac{\bar\varepsilon^{\,2}}{\bar e}$ · **why** Only the destruction term $-C_{\varepsilon2}\bar\varepsilon^2/\bar e$ of the Start's ε̄ equation survives. · **plain** The
     dissipation rate decays too, at a rate set by one constant.
  3. **did** Divide the second by the first · **tex** $\dfrac{d\bar\varepsilon}{d\bar e}=C_{\varepsilon2}\dfrac{\bar\varepsilon}{\bar e}$ · **why** $d\bar\varepsilon/d\bar e=(d\bar\varepsilon/dt)/(d\bar e/dt)$ eliminates time (P303). · **plain** How ε̄
     changes as ē falls.
  4. **did** Integrate · **tex** $\ln\bar\varepsilon=C_{\varepsilon2}\ln\bar e+\text{const}\ \Rightarrow\ \bar\varepsilon=\bar\varepsilon_0\Big(\dfrac{\bar e}{\bar e_0}\Big)^{C_{\varepsilon2}}$ · **why** Separable: $d\bar\varepsilon/\bar\varepsilon=C_{\varepsilon2}\,d\bar e/\bar e$; the constant is fixed by the
     initial values. · **plain** Dissipation is a power of the energy.
  5. **did** Put it back into step 1 · **tex** $\dfrac{d\bar e}{dt}=-\bar\varepsilon_0\Big(\dfrac{\bar e}{\bar e_0}\Big)^{C_{\varepsilon2}}$ · **why** Now one equation for ē alone. · **plain** The decay rate
     depends on how much energy is left.
  6. **did** Separate and integrate · **tex** $\Big(\dfrac{\bar e}{\bar e_0}\Big)^{-(C_{\varepsilon2}-1)}=1+(C_{\varepsilon2}-1)\dfrac{\bar\varepsilon_0}{\bar e_0}\,t$ · **why** $\int\bar e^{-C_{\varepsilon2}}d\bar e=\bar e^{1-C_{\varepsilon2}}/(1-C_{\varepsilon2})$; apply ē = ē₀ at t = 0. ·
     **plain** A power of the energy grows linearly in time.
  7. **did** Solve for ē · **tex** $\bar e=\bar e_0\Big(1+\dfrac{t}{t_0}\Big)^{-n},\quad n=\dfrac{1}{C_{\varepsilon2}-1},\quad t_0=\dfrac{n\,\bar e_0}{\bar\varepsilon_0};\qquad\bar\varepsilon=\bar\varepsilon_0\Big(1+\dfrac{t}{t_0}\Big)^{-(n+1)}$ · **why** Raise step 6 to the power −n; ε̄ from
     step 1 by differentiation. · **plain** Power-law decay: the exponent is measured behind grids, so it fixes $C_{\varepsilon2}$.
  8. **did** In a log layer set production equal to dissipation · **tex** $\bar\varepsilon=-\overline{uv}\dfrac{dU}{dy}=u_*^2\cdot\dfrac{u_*}{\kappa y}=\dfrac{u_*^3}{\kappa y}$ · **why** Steady, with ē uniform (step 10)
     there is no transport in the Start's ē equation; the log-layer values are $-\overline{uv}=u_*^2$, $dU/dy=u_*/(\kappa y)$ from $U^+=\tfrac1\kappa\ln y^++B$ (12.88). · **plain**
     Local equilibrium: what is produced at a height is dissipated there.
  9. **did** Find the eddy viscosity there · **tex** $\nu_T=\dfrac{-\overline{uv}}{dU/dy}=\kappa\,u_*\,y$ · **why** Definition of ν_T for a simple shear (D21 step 1). · **plain**
     The eddy viscosity grows linearly with height.
  10. **did** Use the model's formula for ν_T · **tex** $\kappa u_*y=C_\mu\dfrac{\bar e^{\,2}}{\bar\varepsilon}=C_\mu\,\bar e^{\,2}\,\dfrac{\kappa y}{u_*^3}\ \Rightarrow\ \bar e=\dfrac{u_*^2}{\sqrt{C_\mu}}$ · **why** $\nu_T=C_\mu\bar e^2/\bar\varepsilon$ (12.104) with step 8; y cancels, so ē
      is the same at every height — which justifies dropping its transport. · **plain** In the log layer the turbulent
      energy is a fixed multiple of u_*².
  11. **did** Reduce the ε̄ equation in the log layer · **tex** $0=\dfrac{d}{dy}\Big(\dfrac{\nu_T}{\sigma_\varepsilon}\dfrac{d\bar\varepsilon}{dy}\Big)+\big(C_{\varepsilon1}-C_{\varepsilon2}\big)\dfrac{\bar\varepsilon^{\,2}}{\bar e}$ · **why** Steady, V = 0; the production
      $-\overline{u_iu_j}\,\partial U_i/\partial x_j$ equals ε̄ (step 8), so the source is $+C_{\varepsilon1}\bar\varepsilon^2/\bar e$. Unlike ē, ε̄ varies with y, so its
      diffusion remains. · **plain** Diffusion of ε̄ must make up the difference between its source and its sink.
  12. **did** Insert steps 8–10 · **tex** $\dfrac{u_*^4}{\sigma_\varepsilon y^2}=\big(C_{\varepsilon2}-C_{\varepsilon1}\big)\dfrac{\sqrt{C_\mu}\,u_*^4}{\kappa^2y^2}\ \Rightarrow\ \kappa^2=\sqrt{C_\mu}\,\big(C_{\varepsilon2}-C_{\varepsilon1}\big)\,\sigma_\varepsilon$ · **why** $d\bar\varepsilon/dy=-u_*^3/(\kappa y^2)$,
      $\nu_T\,d\bar\varepsilon/dy=-u_*^4/y$, its derivative $+u_*^4/y^2$; and $\bar\varepsilon^2/\bar e=u_*^4\sqrt{C_\mu}/(\kappa^2y^2)$. · **plain** The model reproduces
      the log law only if its constants obey this relation.
- **Result:** $n=1/(C_{\varepsilon2}-1)$; $\kappa^2=\sqrt{C_\mu}(C_{\varepsilon2}-C_{\varepsilon1})\sigma_\varepsilon$; $\bar e=u_*^2/\sqrt{C_\mu}$ — "the constants are tied to the decay exponent and to
  κ".
- **Check:** units — t₀ = ē₀/ε̄₀ in s ✓; κ dimensionless ✓. Standard set ($C_\mu$ = 0.09, $C_{\varepsilon1}$ = 1.44, $C_{\varepsilon2}$ = 1.92, $\sigma_\varepsilon$ = 1.3;
  Launder & Sharma 1974): n = 1.087, κ = 0.433, ē/u_*² = 3.33. ē₀ = ε̄₀ = 1: ē(10 s) = 0.0801 (closed form = `solve_ivp` = RK4 to
  1e-6). Limit $C_{\varepsilon2}\to2$: n = 1.
- **What it means:** measured decay exponents (about 1.1–1.3) and κ ≈ 0.38–0.41 bracket the standard values: the set is a
  compromise, not a fit to one flow. $C_\mu$ itself comes from the measured ratio $-\overline{uv}/\bar e\approx0.3$ in equilibrium shear layers
  ($\sqrt{C_\mu}=0.3$). **Fails when:** the flow is far from these two calibration states.
- **Traps:** the sign of the $C_{\varepsilon1}$ term (the book writes $-C_{\varepsilon1}(\overline{u_iu_j}\,\partial U_i/\partial x_j)\bar\varepsilon/\bar e$, which is positive because the
  bracket is negative); dropping the diffusion of ε̄ in the log layer (then no relation for κ follows); counting six
  constants (slip #9: there are five).

### D24 · The reduced budget (12.106), the flux Richardson number (12.107) and $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109) — ★, 8 steps, in C14 (notebook · `stratified_surface_layer`)
- **Goal:** one number that says whether a stratified shear flow can keep its turbulence alive, and its link to the number
  one can measure from mean profiles. · **Start:** the turbulent energy budget $\frac{\partial\bar e}{\partial t}+U_j\frac{\partial\bar e}{\partial x_j}=\frac{\partial}{\partial x_j}\big(-\frac1{\rho_0}\overline{pu_j}+2\nu\overline{u_iS'_{ij}}-\frac12\overline{u_i^2u_j}\big)-2\nu\overline{S'_{ij}S'_{ij}}-\overline{u_iu_j}\frac{\partial U_i}{\partial x_j}+g\alpha\overline{u_3T'}$ (12.47) and the closures
  $\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ (12.94), $\overline{u_iT'}=-\kappa_T\,\partial\bar T/\partial x_i$ (12.95). · **Plan:** • specialise the budget to a mean flow U(z) with
  horizontally uniform statistics • take the ratio of the buoyancy term to the shear production • substitute the two
  closures. · **Tools:** horizontal homogeneity; substitution; N² from the potential-temperature gradient (ch01, R05). ·
  **Assumptions:** mean flow U(z) along x; statistics independent of x and y (except $\bar e$ following the flow); Boussinesq;
  $\bar T$, $T'$ are **potential** temperature.
- **Steps:**
  1. **did** Specialise the left side and the transport · **tex** $\dfrac{\partial\bar e}{\partial t}+U\dfrac{\partial\bar e}{\partial x}=-\dfrac{\partial}{\partial z}\Big(\dfrac1{\rho_0}\overline{pw}+\tfrac12\overline{u_i^2w}\Big)+\dots$ · **why** The mean flow is (U(z), 0, 0), so
     $U_j\partial/\partial x_j=U\partial/\partial x$; transport fluxes can vary only with z; viscous transport is negligible at high Re. (Triple correlation written $\tfrac12\overline{u_i^2u_j}$, as in the Start: slip #10.) · **plain** Energy is carried along x by the wind and redistributed
     vertically.
  2. **did** Specialise the shear production · **tex** $-\overline{u_iu_j}\dfrac{\partial U_i}{\partial x_j}=-\overline{uw}\dfrac{dU}{dz}$ · **why** The only non-zero mean gradient is $\partial U_1/\partial x_3=dU/dz$. ·
     **plain** Production by the vertical shear of the wind.
  3. **did** Keep buoyancy and dissipation · **tex** $\dfrac{\partial\bar e}{\partial t}+U\dfrac{\partial\bar e}{\partial x}=-\dfrac{\partial}{\partial z}\Big(\dfrac1{\rho_0}\overline{pw}+\tfrac12\overline{u_i^2w}\Big)-\overline{uw}\dfrac{dU}{dz}+g\alpha\,\overline{wT'}-\bar\varepsilon$ · **why** $u_3=w$; $\bar\varepsilon\equiv2\nu\overline{S'_{ij}S'_{ij}}$. This is (12.106).
     · **plain** Transport, shear production, buoyancy, dissipation.
  4. **did** Read the signs of the two sources · **tex** $-\overline{uw}\dfrac{dU}{dz}>0;\qquad g\alpha\,\overline{wT'}\ \begin{cases}>0&\text{heat flux upward (unstable)}\\<0&\text{heat flux downward (stable)}\end{cases}$ · **why** $\overline{uw}<0$ when dU/dz > 0 (C04). Stable:
     dθ/dz > 0, rising parcels are cold ($T'<0$ with w > 0). · **plain** Shear always feeds turbulence; buoyancy feeds or
     drains it.
  5. **did** Take the ratio · **tex** $\mathrm{Rf}\equiv\dfrac{-g\alpha\,\overline{wT'}}{-\overline{uw}\,(dU/dz)}=\dfrac{\text{buoyant destruction}}{\text{shear production}}$ · **why** Dividing two terms of one budget gives a pure number
     that does not depend on units or on the size of the flow. This is (12.107). · **plain** The fraction of the shear
     production that buoyancy takes away (negative when buoyancy adds).
  6. **did** Model the two fluxes · **tex** $-\overline{uw}=\nu_T\dfrac{dU}{dz},\qquad\overline{wT'}=-\kappa_T\dfrac{d\bar T}{dz}$ · **why** The Start's closures: $\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ (12.94) with i = 1, j = 3 and $\overline{u_iT'}=-\kappa_T\,\partial\bar T/\partial x_i$ (12.95) with i = 3: eddy
     viscosity $\nu_T$ and eddy diffusivity $\kappa_T$ (not the molecular κ, not von Kármán's). · **plain** Fluxes run down the mean
     gradients.
  7. **did** Substitute into Rf · **tex** $\mathrm{Rf}=\dfrac{g\alpha\,\kappa_T\,(d\bar T/dz)}{\nu_T\,(dU/dz)^2}=\dfrac{\kappa_T}{\nu_T}\,\dfrac{g\alpha\,(d\bar T/dz)}{(dU/dz)^2}$ · **why** The two minus signs in the numerator cancel. · **plain** The
     flux ratio in terms of mean gradients.
  8. **did** Recognise the gradient Richardson number · **tex** $\mathrm{Ri}\equiv\dfrac{N^2}{(dU/dz)^2}=\dfrac{\alpha g\,(d\bar T/dz)}{(dU/dz)^2}\ \Rightarrow\ \mathrm{Ri}=\dfrac{\nu_T}{\kappa_T}\,\mathrm{Rf}$ · **why** $\mathrm{Ri}\equiv N^2/(dU/dz)^2$ (12.108) with $N^2=g\alpha\,d\bar T/dz$ for
     potential temperature. This is (12.109). · **plain** The measurable Ri is Rf times the turbulent Prandtl number.
- **Result:** $\mathrm{Rf}=\dfrac{-g\alpha\overline{wT'}}{-\overline{uw}(dU/dz)}$ (12.107); $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$ (12.109) — "Rf compares what buoyancy removes with what
  shear supplies".
- **Check:** units — numerator (m/s²)(1/K)(K m/s) = m²/s³; denominator (m²/s²)(1/s) = m²/s³ ✓. Night example: $\overline{wT'}$ = −0.02 K
  m/s, $\overline{uw}$ = −0.09 m²/s², dU/dz = 0.1 s⁻¹, α = 1/300: Rf = 6.54 × 10⁻⁴/9.0 × 10⁻³ = 0.073. With the thermometer gradient
  dT/dz = +10 K/km: Kundu $dT/dz-\Gamma_a$ = 10 − (−9.8) = 19.8 K/km; meteorology $\Gamma_d-\Gamma_{met}$ = 9.8 − (−10) = 19.8 K/km; Ri = 0.065;
  $\mathrm{Pr}_T$ = 0.89 = Ri/Rf ✓. Neutral: $\overline{wT'}=0$ ⇒ Rf = 0.
- **What it means:** in steady state without transport, $\bar\varepsilon=P(1-\mathrm{Rf})$: observations show turbulence cannot sustain
  itself when Rf exceeds about ¼ (dissipation needs the rest). **Fails when:** transport is important (the top of a
  convective layer), or the flow is not horizontally uniform.
- **Traps:** the sign of Rf (upward heat flux ⇒ negative); using the thermometer gradient in $\mathrm{Ri}=\alpha g(d\bar T/dz)/(dU/dz)^2$ (12.108) without subtracting
  $\Gamma_a$ (an isothermal layer would look neutral instead of strongly stable); confusing the observed $\mathrm{Rf}_{cr}\approx\tfrac14$ with
  Ch. 11's theorem Ri > ¼; three different κ's on one page.

### D25 · $\mathrm{Rf}=z/L_M$ (12.111) and the log-linear wind profile $U=\tfrac{u_*}\kappa\big[\ln\tfrac z{z_0}+5\tfrac z{L_M}\big]$ — ★, 7 steps, in C15 (notebook · `stratified_surface_layer`)
- **Goal:** show that in the surface layer the flux Richardson number is simply height divided by one length, and get the
  wind profile that bends with stability. · **Start:** $\mathrm{Rf}=\dfrac{-g\alpha\overline{wT'}}{-\overline{uw}(dU/dz)}$ (12.107) and $L_M\equiv-\dfrac{u_*^3}{\kappa\,\alpha\,g\,\overline{wT'}}$ (12.110). · **Plan:** • put
  the log-layer values of stress and shear into Rf • recognise L_M • correct the shear for stability to first order •
  integrate from the roughness height. · **Tools:** the log law of C11; substitution; the stability parameter ζ = z/L and the
  flux–profile relation (P304, primer above); $\int(1/z)\,dz=\ln z$. · **Assumptions:** surface (constant-flux) layer: $-\overline{uw}=u_*^2$
  and $\overline{wT'}$ independent of height; z well above $z_0$; mild stratification for the linear correction.
- **Steps:**
  1. **did** Take the stress and shear of the neutral log layer · **tex** $-\overline{uw}=u_*^2,\qquad\dfrac{dU}{dz}=\dfrac{u_*}{\kappa z}$ · **why** Close to the ground the stress is the
     surface stress ($u_*^2\equiv\tau_0/\rho$ (12.81)) and differentiating $U^+=\tfrac1\kappa\ln y^++B$ (12.88) gives the shear. Here κ is
     von Kármán's constant. · **plain** Near the ground: constant stress, shear ∝ 1/z.
  2. **did** Put them into Rf · **tex** $\mathrm{Rf}=\dfrac{-g\alpha\,\overline{wT'}}{u_*^2\cdot u_*/(\kappa z)}=-\dfrac{\kappa\,\alpha\,g\,\overline{wT'}}{u_*^3}\,z$ · **why** Substitution in $\mathrm{Rf}=\dfrac{-g\alpha\overline{wT'}}{-\overline{uw}(dU/dz)}$ (12.107). · **plain** Rf grows in proportion to
     height: shear production fades as 1/z, the heat flux does not.
  3. **did** Recognise the Monin–Obukhov length · **tex** $\mathrm{Rf}=\dfrac{z}{L_M},\qquad L_M\equiv-\dfrac{u_*^3}{\kappa\,\alpha\,g\,\overline{wT'}}$ · **why** The factor multiplying z is exactly $1/L_M$ by $L_M\equiv-u_*^3/(\kappa\alpha g\overline{wT'})$ (12.110); the minus sign makes $L_M>0$ when the heat flux is downward (stable). This is (12.111). · **plain** $L_M$ is the
     height at which buoyancy would equal shear production.
  4. **did** Correct the shear for stability · **tex** $\phi_m\equiv\dfrac{\kappa z}{u_*}\dfrac{dU}{dz}=1+\beta\,\dfrac{z}{L_M}$ · **why** Dimensional analysis: the dimensionless shear can
     depend only on ζ = z/L_M; it is 1 when neutral; for small ζ keep the first term of its Taylor series. β is empirical
     (≈ 5). (Ours.) · **plain** Stable air needs more shear to carry the same stress.
  5. **did** Solve for the shear · **tex** $\dfrac{dU}{dz}=\dfrac{u_*}{\kappa}\Big(\dfrac1z+\dfrac{\beta}{L_M}\Big)$ · **why** Multiply step 4 by $u_*/(\kappa z)$. · **plain** The neutral shear plus a
     height-independent extra.
  6. **did** Integrate from the roughness height · **tex** $U(z)=\dfrac{u_*}{\kappa}\Big[\ln\dfrac{z}{z_0}+\beta\,\dfrac{z-z_0}{L_M}\Big]$ · **why** U = 0 at $z=z_0$, as in the rough-wall law
     $U^+=\tfrac1\kappa\ln(y/y_0)$ (12.93); $\int dz/z=\ln z$, $\int dz=z$. · **plain** A logarithm plus a straight line.
  7. **did** Neglect $z_0$ against z and set β = 5 · **tex** $U=\dfrac{u_*}{\kappa}\Big[\ln\dfrac{z}{z_0}+5\,\dfrac{z}{L_M}\Big]$ · **why** $z_0$ is centimetres, z metres. This is the book's
     log-linear profile (unnumbered, p. 598). · **plain** Stable: more wind aloft than the logarithm; neutral ($L_M\to\infty$): the
     logarithm.
- **Result:** $\mathrm{Rf}=z/L_M$ (12.111); $U=\tfrac{u_*}\kappa[\ln(z/z_0)+5z/L_M]$ — "below ∣L_M∣ the wind makes the turbulence; above it
  buoyancy makes or kills it".
- **Check:** units — $u_*^3/(g\alpha\overline{wT'})$: (m³/s³)/((m/s²)(1/K)(K m/s)) = m ✓. Night: u_* = 0.3 m/s, H = −30 W/m² ⇒ $\overline{wT'}$ = −0.0249 K
  m/s, T = 300 K: L_M = +83 m; Rf(10 m) = 0.12; U(10 m) = 4.81 m/s (neutral 4.36). Day, H = +150 W/m²: L_M = −17 m. Neutral: H = 0
  ⇒ L_M = ∞ ⇒ the rough-wall log law ✓.
- **What it means:** every bulk surface-flux scheme iterates on this profile (and its temperature twin) to get u_* and H
  from model-level wind and temperature. **Fails when:** z ≥ ∣L_M∣/5 on the unstable side — there the log-linear shear 1 + 5z/L_M reaches zero, and beyond it the linear form even gives negative
  winds; use a commonly used form (often called Businger–Dyer), $\phi_m=(1-16\zeta)^{-1/4}$, there — and for very stable layers (ζ ≳ 1), where
  turbulence is intermittent.
- **Traps:** the sign of L_M; κ is von Kármán's constant here although the temperature-variance budget on the next page (N190 in C15) uses κ for the thermal diffusivity (our $\kappa_{th}$); α = 1/T needs T in kelvin; H in W/m² must be divided by ρc_p.

### D26 · Taylor's formula: (12.115)–(12.117), the double integral (12.118) and $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\,t\int_0^t(1-\tau/t)r_\alpha\,d\tau$ (12.119) — ★★, 11 steps, in C16 (notebook · `taylor_dispersion`)
- **Goal:** express how far marked particles have spread, on average, purely in terms of the statistics of their velocity.
  · **Start:** a particle released at the origin at t = 0: $\dfrac{dX_\alpha}{dt}=u_\alpha(t)$, where $u_\alpha$ is the velocity *of that particle*
  (Lagrangian) along one axis α (no sum). · **Plan:** • differentiate the mean-square displacement • write X as the time
  integral of u • recognise the velocity autocorrelation • integrate in time • simplify the double integral by parts. ·
  **Tools:** D01's rules $\overline{\partial u/\partial t}=\partial\bar u/\partial t$ (12.6) and $\overline{\int u\,dt}=\int\bar u\,dt$ (12.7); chain rule (P49); substitution (P106); even functions (P261); a double integral
  over a triangle (P305, primer above); integration by parts (P218a). · **Assumptions:** stationary homogeneous turbulence;
  zero mean velocity; average over many particles (or releases).
- **Steps:**
  1. **did** Write the displacement as an integral of the velocity · **tex** $X_\alpha(t)=\displaystyle\int_0^tu_\alpha(t')\,dt'$ · **why** Integrate $dX_\alpha/dt=u_\alpha$ from the release
     (X = 0 at t = 0); t′ is a dummy time. · **plain** Where a particle is = the sum of all its past velocities.
  2. **did** Differentiate the square · **tex** $\dfrac{d}{dt}\big(X_\alpha^2\big)=2X_\alpha\dfrac{dX_\alpha}{dt}$ · **why** Chain rule (P49). · **plain** The square grows at twice displacement ×
     velocity.
  3. **did** Average over particles · **tex** $\dfrac{d}{dt}\big(\overline{X_\alpha^2}\big)=2\,\overline{X_\alpha\dfrac{dX_\alpha}{dt}}$ · **why** The average commutes with d/dt, $\overline{\partial u^m/\partial t}=\partial\overline{u^m}/\partial t$ (12.6). This is
     (12.115). · **plain** The growth rate of the cloud's mean-square size.
  4. **did** Replace dX/dt by the velocity · **tex** $\dfrac{d}{dt}\big(\overline{X_\alpha^2}\big)=2\,\overline{X_\alpha u_\alpha}$ · **why** Definition of the particle's velocity. · **plain** The cloud
     grows when displaced particles are still moving outward.
  5. **did** Insert step 1 and move the average inside · **tex** $\dfrac{d}{dt}\big(\overline{X_\alpha^2}\big)=2\,\overline{\Big[\displaystyle\int_0^tu_\alpha(t')dt'\Big]u_\alpha(t)}=2\displaystyle\int_0^t\overline{u_\alpha(t')\,u_\alpha(t)}\,dt'$ · **why** $u_\alpha(t)$ does not depend
     on t′, so it goes inside the integral; the average commutes with time integration, $\overline{\int u\,dt}=\int\bar u\,dt$ (12.7). This is
     (12.116). · **plain** The growth rate is the accumulated correlation between present and past velocity.
  6. **did** Use stationarity · **tex** $\overline{u_\alpha(t')\,u_\alpha(t)}=\overline{u_\alpha^2}\;r_\alpha(t-t')$ · **why** In stationary turbulence $\overline{u_\alpha^2}$ is constant and the correlation depends
     only on the time difference: $r_\alpha(\tau)\equiv\overline{u_\alpha(t)u_\alpha(t+\tau)}/\overline{u_\alpha^2}$. The book writes $r_\alpha(t'-t)$; the two are equal because r is even
     (D03). · **plain** How well a particle remembers the velocity it had a time t − t′ ago.
  7. **did** Change variable to the lag · **tex** $\dfrac{d}{dt}\big(\overline{X_\alpha^2}\big)=2\,\overline{u_\alpha^2}\displaystyle\int_0^tr_\alpha(\tau)\,d\tau,\qquad\tau=t-t'$ · **why** dτ = −dt′; t′ = 0 → τ = t and t′ = t → τ = 0; the minus
     sign swaps the limits back. This is (12.117). · **plain** The cloud grows at twice the variance times the area under
     the correlation up to now.
  8. **did** Integrate in time from the release · **tex** $\overline{X_\alpha^2}(t)=2\,\overline{u_\alpha^2}\displaystyle\int_0^tdt'\int_0^{t'}r_\alpha(\tau)\,d\tau$ · **why** $\overline{X_\alpha^2}=0$ at t = 0; rename the upper limit of step 7
     as t′ and integrate it from 0 to t. This is (12.118). · **plain** A double integral over a triangle in the (t′, τ)
     plane.
  9. **did** Integrate by parts in t′ · **tex** $\displaystyle\int_0^tdt'\int_0^{t'}r_\alpha\,d\tau=\Big[t'\int_0^{t'}r_\alpha(\tau)\,d\tau\Big]_{t'=0}^{t}-\int_0^tt'\,r_\alpha(t')\,dt'$ · **why** $\int v'w=[vw]-\int vw'$ with v = t′ and
     $w=\int_0^{t'}r\,d\tau$, whose derivative is $r(t')$ (fundamental theorem, P84). · **plain** Trade the inner integral for a
     weight.
  10. **did** Evaluate the bracket and rename · **tex** $=t\displaystyle\int_0^tr_\alpha(\tau)\,d\tau-\int_0^t\tau\,r_\alpha(\tau)\,d\tau=t\int_0^t\Big(1-\dfrac{\tau}{t}\Big)r_\alpha(\tau)\,d\tau$ · **why** The bracket vanishes at t′ = 0; t′ in the last
      integral is a dummy, renamed τ; factor out t. Same result as counting strips of the triangle (P305). · **plain** Each
      lag τ is weighted by the fraction of the time (1 − τ/t) for which it was available.
  11. **did** Assemble · **tex** $\overline{X_\alpha^2}(t)=2\,\overline{u_\alpha^2}\;t\displaystyle\int_0^t\Big(1-\dfrac{\tau}{t}\Big)r_\alpha(\tau)\,d\tau$ · **why** Step 10 in step 8. This is (12.119). · **plain** The spread of the
      cloud needs only the velocity variance and its autocorrelation.
- **Result:** $\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\,t\int_0^t\big(1-\frac\tau t\big)r_\alpha(\tau)\,d\tau$ (12.119) — "dispersion is the velocity autocorrelation, integrated twice".
- **Check:** units — (m²/s²)(s)(s) = m² ✓. The double integral $\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\int_0^tdt'\int_0^{t'}r_\alpha(\tau)\,d\tau$ (12.118) and the Result agree numerically to 1e-8 for $r=e^{-\tau/\Lambda_t}$.
  Differentiating the Result returns $\frac{d}{dt}\overline{X_\alpha^2}=2\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$ (12.117) ✓. Particles: 10⁴ Langevin particles reproduce it within 5 standard errors.
- **What it means:** no eddy diffusivity was assumed — the result is exact kinematics; all the physics is in r(τ) and its
  integral scale Λ_t. **Fails when:** the turbulence is not homogeneous or stationary along the particle's path (a plume
  rising through a boundary layer), or there is mean shear (shear dispersion is faster).
- **Traps:** using the Eulerian (fixed-point) correlation — the formula needs the correlation *following a particle*;
  summing over α (it is a label); the limits when changing t′ → τ.

### D27 · The ballistic limit (12.120)–(12.121), the diffusive limit (12.122)–(12.123) with its offset, and the closed form for $r=e^{-\tau/\Lambda_t}$ — ★, 8 steps, in C16 (notebook · `taylor_dispersion`)
- **Goal:** read the two simple laws — spreading like t, then like √t — out of Taylor's formula, and get one formula that
  covers everything in between. · **Start:** $\overline{X_\alpha^2}(t)=2\overline{u_\alpha^2}\,t\displaystyle\int_0^t\Big(1-\dfrac{\tau}{t}\Big)r_\alpha(\tau)\,d\tau$ (12.119). · **Plan:** • short times: r ≈ 1 • long
  times: the integrals reach their limits • exponential r: do the integrals exactly. · **Tools:** limits of integrals;
  improper integrals (P144); integration by parts (P218a); sympy (P40). · **Assumptions:** as D26; r integrable with
  integral scale $\Lambda_t\equiv\int_0^\infty r_\alpha\,d\tau$ (Lagrangian).
- **Steps:**
  1. **did** For short times set r to one · **tex** $t\ll\Lambda_t:\quad r_\alpha(\tau)\approx1\ \ \text{for}\ 0\le\tau\le t$ · **why** r(0) = 1 and r changes only over times of order Λ_t:
     the particle has not yet forgotten its velocity. · **plain** Perfect memory.
  2. **did** Do the integral · **tex** $\overline{X_\alpha^2}\simeq2\,\overline{u_\alpha^2}\,t\displaystyle\int_0^t\Big(1-\dfrac{\tau}{t}\Big)d\tau=\overline{u_\alpha^2}\,t^2\ \Rightarrow\ (X_\alpha)_{rms}=(u_\alpha)_{rms}\,t$ · **why** $\int_0^t(1-\tau/t)\,d\tau=t-t/2=t/2$. These are
     (12.120) and (12.121). · **plain** Ballistic: each particle flies straight, so the cloud grows in proportion to
     time.
  3. **did** For long times split the integral · **tex** $\overline{X_\alpha^2}=2\,\overline{u_\alpha^2}\Big[t\displaystyle\int_0^tr_\alpha\,d\tau-\int_0^t\tau\,r_\alpha\,d\tau\Big]$ · **why** Multiply out (1 − τ/t); this is step 10 of D26 read
     backwards. · **plain** A part that grows with t and a correction.
  4. **did** Let both integrals reach their limits · **tex** $t\gg\Lambda_t:\quad\displaystyle\int_0^tr_\alpha\,d\tau\to\Lambda_t,\qquad\int_0^t\tau\,r_\alpha\,d\tau\to\int_0^\infty\tau\,r_\alpha\,d\tau=\text{const}$ · **why** r has decayed to zero
     well before t, so extending the upper limits to ∞ changes nothing (P144); the second integral converges if r decays
     faster than 1/τ². · **plain** The memory is used up.
  5. **did** Write the long-time law · **tex** $\overline{X_\alpha^2}\simeq2\,\overline{u_\alpha^2}\,\Lambda_t\,t-2\,\overline{u_\alpha^2}\displaystyle\int_0^\infty\tau\,r_\alpha\,d\tau\ \approx\ 2\,\overline{u_\alpha^2}\,\Lambda_t\,t\ \Rightarrow\ (X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}$ · **why** The second
     term is a constant, negligible against the first as t → ∞. These are (12.122) and (12.123). · **plain** Diffusive:
     the cloud grows like the square root of time, as in a random walk.
  6. **did** For an exponential memory do the first integral · **tex** $r_\alpha=e^{-\tau/\Lambda_t}:\quad\displaystyle\int_0^tr_\alpha\,d\tau=\Lambda_t\big(1-e^{-t/\Lambda_t}\big)$ · **why** Elementary integral; this
     shape is what a Langevin particle has (P283). · **plain** The remembered fraction of the memory time.
  7. **did** Do the second integral · **tex** $\displaystyle\int_0^t\tau\,e^{-\tau/\Lambda_t}d\tau=\Lambda_t^2\Big[1-e^{-t/\Lambda_t}\Big(1+\dfrac{t}{\Lambda_t}\Big)\Big]$ · **why** Integration by parts (P218a) with u = τ. · **plain** The
     correction, in closed form.
  8. **did** Combine · **tex** $\overline{X_\alpha^2}=2\,\overline{u_\alpha^2}\,\Lambda_t^2\Big[\dfrac{t}{\Lambda_t}-1+e^{-t/\Lambda_t}\Big]$ · **why** Step 3 with steps 6 and 7: the two $t\Lambda_te^{-t/\Lambda_t}$ terms cancel. (Ours.) ·
     **plain** One formula from the wedge to the parabola.
- **Result:** $(X_\alpha)_{rms}=(u_\alpha)_{rms}t$ for $t\ll\Lambda_t$ (12.121); $(X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}$ for $t\gg\Lambda_t$ (12.123); exponential memory:
  $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\Lambda_t^2[t/\Lambda_t-1+e^{-t/\Lambda_t}]$ — "t first, √t later, with the change at t ≈ Λ_t".
- **Check:** small t: $e^{-x}\approx1-x+x^2/2$ ⇒ $\overline{u^2}t^2$ ✓; large t: $2\overline{u^2}\Lambda_t(t-\Lambda_t)$ ✓ (offset $-2\overline{u^2}\Lambda_t^2$, since
  $\int_0^\infty\tau e^{-\tau/\Lambda_t}d\tau=\Lambda_t^2$). Numbers (u_rms = 1 m/s, Λ_t = 10 s): 0.967, 73.6, 1800 m² at t = 1, 10, 100 s (ballistic 1, 100,
  10⁴; diffusive 20, 200, 2000). Gaussian memory $r=e^{-\tau^2/t_c^2}$: $\Lambda_t=\tfrac{\sqrt\pi}2t_c$ and
  $\overline{X^2}=2\overline{u^2}[\tfrac{\sqrt\pi}2t_ct\,\mathrm{erf}(t/t_c)-\tfrac{t_c^2}2(1-e^{-t^2/t_c^2})]$ (86.15 m² at t = t_c = 10 s, checked by quadrature).
- **What it means:** near a source turbulent spreading is not diffusion at all; far from it, it is — with diffusivity
  $\overline{u^2}\Lambda_t$ (D28). On a plot of $\overline{X^2}$ against t the long-time line does not pass through the origin. **Fails when:** r has
  a long tail (anomalous diffusion), or for *relative* dispersion of particle pairs (Richardson's law instead).
- **Traps:** using the diffusive law at t ~ Λ_t (at t = Λ_t it gives an rms displacement 65 % too large — a factor 2.7 in the mean square); forgetting that Λ_t here is Lagrangian; the
  book's reference "(11.119)" (slip #11) means $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\,t\int_0^t(1-\tau/t)r_\alpha\,d\tau$ (12.119).

### D28 · The eddy diffusivity: $\nu=\tfrac12\,d\sigma^2/dt$ (12.126), $D_T=\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$ (12.127) and its limits (12.128)–(12.129, condition corrected) — ★, 7 steps, in C16 (notebook · `taylor_dispersion`)
- **Goal:** define an "effective diffusivity" for turbulent spreading and see why it is not a constant. · **Start:** the
  molecular yardstick (R07): a diffusing patch is a Gaussian whose variance per coordinate grows as $\sigma^2=2\nu t$. · **Plan:** •
  read a diffusivity off the growth of a variance • apply the same reading to the cloud of particles • take the two limits.
  · **Tools:** Gaussian spreading (ch01 σ² = 2Dt, R08); D26's $\frac{d}{dt}\overline{X_\alpha^2}=2\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$ (12.117); limits (D27). · **Assumptions:** as D26.
- **Steps:**
  1. **did** Recall how a diffusing patch grows · **tex** $\sigma^2=2\,\nu\,t$ · **why** The Gaussian solution of the diffusion equation,
     $\propto\exp(-r^2/4\nu t)$, has variance 2νt along each coordinate (Ch. 1; R07). · **plain** Under diffusion the variance grows
     linearly in time.
  2. **did** Solve for the diffusivity · **tex** $\nu=\dfrac12\dfrac{d\sigma^2}{dt}$ · **why** Differentiate step 1. This is (12.126): a diffusivity is half the growth
     rate of a variance. · **plain** A recipe to *measure* a diffusivity from a spreading cloud.
  3. **did** Apply the recipe to turbulent dispersion · **tex** $D_T\equiv\dfrac12\dfrac{d}{dt}\big(\overline{X_\alpha^2}\big)$ · **why** $\overline{X_\alpha^2}$ is the variance of the particle cloud along α; this
     defines an effective (eddy) diffusivity. First part of (12.127). · **plain** How fast the cloud's variance is growing
     right now.
  4. **did** Insert the rate from Taylor's analysis · **tex** $D_T=\overline{u_\alpha^2}\displaystyle\int_0^tr_\alpha(\tau)\,d\tau$ · **why** $\frac{d}{dt}\overline{X_\alpha^2}=2\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$ (12.117), halved. Second part of
     (12.127). · **plain** The eddy diffusivity is the velocity variance times the memory used so far.
  5. **did** Short times · **tex** $t\ll\Lambda_t:\quad r_\alpha\approx1\ \Rightarrow\ D_T\cong\overline{u_\alpha^2}\,t$ · **why** As D27 step 1. This is (12.128). · **plain** At first the "diffusivity"
     grows in proportion to time — not a diffusion process.
  6. **did** Long times · **tex** $t\gg\Lambda_t:\quad\displaystyle\int_0^tr_\alpha\,d\tau\to\Lambda_t\ \Rightarrow\ D_T\cong\overline{u_\alpha^2}\,\Lambda_t$ · **why** As D27 step 4. This is (12.129) with the condition
     corrected: the page prints $t\ll\Lambda_t$ (slip #12). · **plain** Only after several memory times does the diffusivity settle
     to a constant.
  7. **did** Read it as a length times a velocity · **tex** $D_T=u_{rms}\times\big(u_{rms}\Lambda_t\big)=u_T\,l_T$ · **why** $\overline{u_\alpha^2}=u_{rms}^2$; $u_{rms}\Lambda_t$ is the distance a particle
     travels while it remembers its velocity. This is $\nu_T\sim l_Tu_T$ (12.98), now derived. · **plain** The eddy diffusivity of
     C12, with its length identified.
- **Result:** $D_T\equiv\tfrac12\frac{d}{dt}\overline{X_\alpha^2}=\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$ (12.127); $D_T\cong\overline{u_\alpha^2}t$ for $t\ll\Lambda_t$ (12.128); $D_T\cong\overline{u_\alpha^2}\Lambda_t$ for $t\gg\Lambda_t$
  (12.129) — "an eddy diffusivity grows before it saturates".
- **Check:** units — (m²/s²)(s) = m²/s ✓. Exponential memory: $D_T=\overline{u^2}\Lambda_t(1-e^{-t/\Lambda_t})$: 0.95, 6.32, 10.0 m²/s at t = 1, 10, 100 s
  for u_rms = 1 m/s, Λ_t = 10 s. Against air's molecular ν = 1.5 × 10⁻⁵ m²/s: 10⁶ times larger.
- **What it means:** a constant-D (Fickian) model of a plume is right far downwind and wrong near the stack, where the
  plume is a wedge; and for a *patch* that keeps meeting larger eddies as it grows, the effective diffusivity keeps
  growing with its size, $K\sim\bar\varepsilon^{1/3}l^{4/3}$ (Richardson). **Fails when:** the scale of the cloud is comparable with the
  largest eddies of an inhomogeneous flow.
- **Traps:** σ² = 2νt per coordinate here, while Ch. 3's vortex core radius used 4νt; treating $D_T$ as a property of the
  fluid; the printed condition $t\ll\Lambda_t$ on $D_T\cong\overline{u_\alpha^2}\Lambda_t$ (12.129) (⚠️ slip #12: it holds for $t\gg\Lambda_t$).

---

## Closing check — design gate (audit of 2026-10-07; bullets only, so that the ledger parser of Part E is not disturbed)

- **CORE blocks:** 16 of 16 (`#### C01 —` … `#### C16 —`, each with one `nb.core`), every one with plain words, the idea,
  primers, derivation(s), a worked example, code, a from-scratch `nb.check_agree`, at least one visual with see / read /
  change, and a "What would change if…" link (C16's was missing and is now row 22 of A.12).
- **Sections:** 13 of 13 `nb.section` rows (12.1–12.13). **Curation ids:** 243 of 243 placed (16 C · 217 N · 8 R · 2 S);
  8 `nb.recap`, 2 `nb.pointer`.
- **Derivations:** 28 of 28 blocks in Part F and 28 `nb.derivation` rows in Part A; step counts equal the curation's
  (8, 6, 5, 7, 9, 13, 14, 13, 9, 15, 9, 9, 10, 14, 10, 12, 9, 7, 12, 5, 11, 6, 12, 8, 7, 11, 8, 7 = 266 steps); every step has
  did / tex / why / plain; every block has Goal, Start, Plan, Tools, Assumptions, Result, Check, What it means, Traps; the
  four ★★★ blocks (D07, D08, D10, D14) carry a sympy check intent with `check_src`.
- **Explainers:** E1–E10 and B1, each with title, summary, CORE ids, reference explainer, meta, physics, views, controls,
  presets, status, depth features, Explain (numbered sections ending in *Reading the current setting*), Derivation tab,
  Code, walkthrough (6, 7, 6, 7, 6, 7, 8, 6, 7, 7 and 5 steps), equations, ≥ 3 check questions, ≥ 2 parity rows, fit plan.
  The ★★★ derivations of blocks that have an explainer are in it (D10 in E5, D14 in E6); C05 (D07, D08) has none.
- **Primers:** 27 new (P280–P306), each with an `nb.primer` row in Part A and a "primer (in Cxx)" row in Part E.
- **Ledger:** 231 rows, none without "Explained by"; every C/R id in that column exists in the curation.
- **Part C:** every `ch12.` / `TS.` / `WT.` name used in Parts A, B, D, E, F is in Part C; every Part C name but one
  (C.6 item 1) is called or named by a storyboard row; nothing was added, renamed or re-signed.
- **Equation mentions (inlining pass of 2026-10-07):** no equation is cited by its number alone in the header, Part A,
  Part B or Part F — a scratch scanner (kept outside the repo) finds 0, down from 214 (header 48 · A 57 · B 2 · F 107). Each
  mention has the equation beside it, in TeX (Unicode in the Part B drafts), or was reworded to point at the Start, the
  Result or the note that writes the equation out. Not counted as mentions, by design: heading lines (`#### Cnn —`,
  `### Dnn ·`; parsed by tools, the `nb.core` title or the block's Start carries the equation); comments inside code
  listings and `check_src`; the book's wrong reference quoted in the row of slip #11; and the 73 own-line labels of
  Part F ("This is (N.M)." directly under the step's own `tex` line). Part C, Part D and the ledger of Part E were not
  scanned and not changed.
- **Open at hand-over:** the three builder notes of C.6 item 3; the expected pair of C11 row 10 (`WT.fit_log_law`),
  pending the verifier's ruling (`reports/ch12_verification.md`, flagged item 1).

