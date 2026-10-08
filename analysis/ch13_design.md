# Chapter 13 — Geophysical Fluid Dynamics: lesson design
(from `analysis/ch13_curation.md` (A 17 · B 154 · C 47 = 218 rows; CORE 17 · NOTE 168 · RECAP 31 · SKIP 2; 29 derivations D01–D29
(★ 2 · ★★ 18 · ★★★ 9: D06, D11, D12, D16, D17, D19, D22, D25, D26); E1–E10 + backup B1; primers from **P307**) and
`analysis/ch13.md` (§2 inventory with the LaTeX the analyst read from the page images, §2b derivations d01…d42, §4
implementation rows I01–I29, §5 the three new core modules, §9 slips #1–#13 and traps T1–T17). 2026-10-07, lesson-designer.

**Equation provenance (honest statement).** Every equation placed here is taken from the analyst's page-image transcription in
`analysis/ch13.md` §2. For this design I re-read 21 rendered pages against that transcription (PDF page = printed page + 27):
p627 (the six stresses $\tau_{xz}=\tau_{zx}=\rho\nu_v\frac{\partial u}{\partial z}+\rho\nu_H\frac{\partial w}{\partial x}$ … of (13.5) and the three friction components of (13.6), the first being $F_x=\frac{\partial\tau_{xx}}{\partial x}+\frac{\partial\tau_{xy}}{\partial y}+\frac{\partial\tau_{xz}}{\partial z}=\nu_H\big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\big)+\nu_v\frac{\partial^2u}{\partial z^2}$ — slip #2 seen: no 1/ρ);
p629 (the determinant, $2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)$ (13.7), $f=2\Omega\sin\theta$ (13.8), $T_i=2\pi/f$ — slip #5 seen — and the three members of (13.9), the first being $\frac{Du}{Dt}-fv=-\frac1{\rho_0}\frac{\partial p}{\partial x}+\nu_H\big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\big)+\nu_v\frac{\partial^2u}{\partial z^2}$);
p632 ($0=-\frac{\partial p}{\partial z}-g\rho$ (13.14), the thermal-wind pair $\frac{\partial v}{\partial z}=-\frac g{\rho_0f}\frac{\partial\rho}{\partial x}$, $\frac{\partial u}{\partial z}=\frac g{\rho_0f}\frac{\partial\rho}{\partial y}$ (13.15), the tank pair — slip #3 seen: the second is printed $-2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$ under (13.17) — and $E=\frac{\rho\nu U/L^2}{\rho fU}=\frac\nu{fL^2}$ (13.18));
p633 ($2\Omega\big(\frac{\partial v}{\partial y}+\frac{\partial u}{\partial x}\big)=0$, $\frac{\partial w}{\partial z}=0$ (13.19), $\frac{\partial v}{\partial z}=\frac{\partial u}{\partial z}=0$ (13.20), $\partial\mathbf u/\partial z=0$ (13.21));
p635 ($-fv=\nu_v\frac{d^2u}{dz^2}$, $fu=\nu_v\frac{d^2v}{dz^2}$ (13.22)–(13.23), the three conditions — slip #4 seen: $z\to\infty$ — $\frac{d^2V}{dz^2}=\frac{if}{\nu_v}V$ (13.27), $V=Ae^{(1+i)z/\delta}+Be^{-(1+i)z/\delta}$ with $\delta=\sqrt{2\nu_v/f}$ (13.28)–(13.29), $A=\frac{\tau\delta(1-i)}{2\rho\nu_v}$ and the cos / sin forms of u and v);
p636 ($\int_{-\infty}^0u\,dz=0$, $\int_{-\infty}^0v\,dz=-\frac\tau{\rho f}$ (13.30) and the remark on $-\rho fv=d\tau/dz$);
p640 ($V=Ae^{-(1+i)z/\delta}+Be^{(1+i)z/\delta}+U$ (13.40), $u=U[1-e^{-z/\delta}\cos(z/\delta)]$, $v=Ue^{-z/\delta}\sin(z/\delta)$ (13.41));
p643 ($(H+\eta)\frac{\partial u}{\partial x}+(H+\eta)\frac{\partial v}{\partial y}+w(\eta)-w(0)=0$ (13.43), $\frac{\partial\eta}{\partial t}+\frac{\partial}{\partial x}[u(H+\eta)]+\frac{\partial}{\partial y}[v(H+\eta)]=0$ (13.44), the three members of (13.45), $c^2=gH_e$ (13.46));
p645 ($\frac{d\psi_n/dz}{N^2\int_{-H}^z\psi_n\,dz}=\frac{\rho_0}g\frac{w_n}{\partial\rho_n/\partial t}\equiv-\frac1{c_n^2}$ (13.55), $\frac d{dz}\big(\frac1{N^2}\frac{d\psi_n}{dz}\big)+\frac1{c_n^2}\psi_n=0$ (13.56), the modal set (13.57)–(13.61) and $c_n^2\equiv gH_e$ (13.62));
p650 ((13.72)–(13.74), the vorticity equation $\frac\partial{\partial t}\big(\frac{\partial u}{\partial y}-\frac{\partial v}{\partial x}\big)-f_0\big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\big)-\beta v=0$, $\frac{\partial^3v}{\partial t^3}-gH\frac\partial{\partial t}\nabla_H^2v+f_0^2\frac{\partial v}{\partial t}-gH\beta\frac{\partial v}{\partial x}=0$ (13.75) and $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76));
p652 ($\omega^2-f^2=gH(k^2+l^2)$ (13.81), $\omega^2=f^2+gHK^2$ (13.82));
p657 ($\eta=\eta_0e^{-fy/c}\cos k(x-ct)$, $u=\eta_0\sqrt{g/H}\,e^{-fy/c}\cos k(x-ct)$ (13.87), $\Lambda\equiv c/f$);
p660 ($\frac{D(\zeta+f)}{Dt}=\frac{\zeta+f_0}h\frac{Dh}{Dt}$ (13.93), $\frac{Df}{Dt}=v\beta$, $\frac D{Dt}\big(\frac{\zeta+f}h\big)=0$ (13.94), $\frac f{h_0}=\frac{\zeta+f}{h_1}$);
p663 ($\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w+f^2\frac{\partial^2w}{\partial z^2}=0$ (13.96) — slip #8 seen — (13.97)–(13.100), $f<\omega<N$);
p668–p669 ($\omega^2-f^2=\frac{k^2}{m^2}(N^2-\omega^2)$ (13.112); **the form $\omega^2=f^2\sin^2\theta+N^2\cos^2\theta$ is printed on the next page without a number** — this design therefore writes it as "(unnumbered, p. 669)" and never under (13.112); the three regime forms);
p673 ((13.116), $\frac\partial{\partial t}\big(\frac{\partial^2\eta}{\partial x^2}+\frac{\partial^2\eta}{\partial y^2}-\frac{f_0^2}{c^2}\eta\big)+\beta\frac{\partial\eta}{\partial x}=0$ (13.117), $\omega=-\frac{\beta k}{k^2+l^2+f_0^2/c^2}$ (13.118));
p677 ($\frac\partial{\partial t}(\nabla^2\psi)+U\frac\partial{\partial x}(\nabla^2\psi)+\big(\beta-\frac{d^2U}{dy^2}\big)\frac{\partial\psi}{\partial x}=0$ (13.123), the normal-mode equation, $\frac d{dy}(\bar\zeta+f)=\beta-\frac{d^2U}{dy^2}$ (13.124));
p681 ((13.133)–(13.135), $\big(\frac\partial{\partial t}+U\frac\partial{\partial x}\big)\big[\nabla_H^2p'+\frac{f^2}{N^2}\frac{\partial^2p'}{\partial z^2}\big]=0$ (13.136), (13.137)–(13.140));
p682 (the two lid conditions, the 2 × 2 system, $c=\frac{U_0}2\pm\frac{U_0}{\alpha H}\sqrt{\big(\frac{\alpha H}2-\tanh\frac{\alpha H}2\big)\big(\frac{\alpha H}2-\coth\frac{\alpha H}2\big)}$ (13.141), $\frac{\alpha_cH}2=\coth\frac{\alpha_cH}2$);
p687 (the two members of (13.145), $\frac{S_1}{S_2}=\frac{K_2-K_0}{K_0-K_1}\frac{K_2+K_0}{K_1+K_0}$ and $\frac{K_1^2S_1}{K_2^2S_2}=\frac{K_1^2}{K_2^2}\frac{K_2^2-K_0^2}{K_0^2-K_1^2}$, and the two spectral ranges).
No difference from the transcription was found except the placement of that one number, described under p668–p669. Three pages of earlier chapters were read for the cross-references: ch04 p117 (the rotating-frame equation $\rho\big(\frac{D'\mathbf u'}{Dt}\big)_{O'1'2'3'}=-\nabla'p+\rho\big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\big]+\mu\nabla'^2\mathbf u'$ (4.45)), ch04 p135 ($\rho\frac{D\mathbf u}{Dt}=-\nabla p'+\rho'\mathbf g+\mu\nabla^2\mathbf u$ (4.84)) and ch01 p21 ($\frac{d\rho_\theta}{dz}=\frac{d\rho}{dz}-\frac{d\rho_a}{dz}\cong\frac{d\rho}{dz}+\frac{\rho g}{c^2}$ (1.35)); the other earlier-chapter equations quoted here are taken from the written-out forms in `analysis/ch01_design.md` … `analysis/ch12_design.md`. **Two cross-references of the curation were corrected against those pages:** the hydrostatic equation is $dp/dz=-\rho g$ (1.8), not (1.14) as recap R07 was annotated; and the number (4.89) belongs to the Boussinesq heat equation $DT/Dt=\kappa\nabla^2T$ (4.89), not to a density equation (recap R06). The other 49 pages
were **not** re-read by the designer; builders render the page (`tools/render_pages.py ch13 --eq N.M`) before setting any
equation not in that list. Numbers marked *expect* were computed for this design with a scratch script (numpy / scipy / sympy
closed forms; no `fluidpy` chapter-13 function existed yet) from the inputs stated beside them; the builder compares an
executed cell against them and **reports** any difference instead of editing the expectation. The nine `check_src` sketches of
Part F were executed as written (sympy 1.x, each under 3 s) and passed.

**Part C is written first and is the contract both the implementer and the builders keep** (first saved 2026-10-07 19:39 EDT,
before any other part).

**Binding conventions for every builder (analysis §9, curation decisions 4–11, §8, §9).**
1. **Imports and aliases.** `from fluidpy import ch13_geophysical_fluid_dynamics as ch13`; `from fluidpy.core import gfd as GFD,
   vertical_modes as VM, shallow_water as SW`. **`ch13` re-exports every public name of `GFD`, `VM` and `SW`** (the ch11 / ch12
   pattern; one alias: `SW.potential_vorticity` → `ch13.sw_potential_vorticity`), so explainer parity rows write `ch13.<name>`
   only. Reused: `from fluidpy.core import rotating as ROT, stratification as STRAT, similarity as SIM, waves as WAV, stability
   as ST, laminar as LAM, turbstats as TS, dimensional as DIM`; `from fluidpy import ch05_vorticity_dynamics as ch05,
   ch11_instability as ch11, ch12_turbulence as ch12`.
2. **Scalar-callable and parity-friendly** (the `shot.py` evaluator knows only `np`, `math`, `ch13`): every public closed form
   accepts Python floats and returns a float, a complex, a tuple or a `dict` of floats (arrays only when an array goes in).
   Parity `py:` rows use dict keys, integer indices and `np.…` only (never ndarray methods or builtins); details in C.6.
3. **Units.** SI. Gravity is g = `G0` = 9.80665 m/s² unless a cell says "g ≈ 10" for hand arithmetic. Angles in radians inside functions
   (`lat_rad`); degrees only in prose and through `np.deg2rad`. Transports per unit width [m²/s]. Frequencies in rad/s; periods
   are quoted in hours or days with the conversion shown.
4. **Symbols (overloading resolved — the conventions block before C01 shows this table; each block repeats the row it needs).**

| Book symbol | Meanings in this chapter (and earlier) | Notebook / explainer symbol | Code name |
|---|---|---|---|
| f | Coriolis parameter $f=2\Omega\sin\theta$ (13.8), signed; earlier chapters: a generic function, a frequency in Hz, the similarity function of ch09, the correlation function of ch12 | $f$ = Coriolis parameter only; $f_0$ its value at the central latitude; a frequency is always $\omega$ [rad/s]; a generic function is $F(\cdot)$ or $G(\cdot)$ | `f`, `f0` |
| Ω, ω | The capital Ω is the earth's rotation rate (vector $\boldsymbol\Omega$, magnitude Ω); ω wave frequency; $\omega_x$, $\omega_y$ horizontal vorticity components in $-f\frac{dv}{dz}=\nu_v\frac{d^2\omega_y}{dz^2}$ (13.31) | Both as in the book (Ω, ω); horizontal vorticity always with its subscript | `Omega`, `omega`, `omega_x` |
| β | planetary vorticity gradient $\beta=2\Omega\cos\theta_0/R$ in $f=f_0+\beta y$ (13.10); on p. 626 the haline contraction coefficient in $\delta\rho/\rho_0=\beta\,\delta S$ | Only df/dy is called β; the haline coefficient is written $\beta_S$ in our own lines | `beta`, `beta_S` |
| N | buoyancy frequency $N^2\equiv-\frac g{\rho_0}\frac{d\rho}{dz}$ (potential density); earlier: a count (ch12's N members) | $N$ = buoyancy frequency; counts are $n$, `nx`, `n_modes` | `N`, `N2` |
| H, h, η | H mean layer depth (§13.8), ocean depth (§13.9), the scale of variation of N in the WKB condition $Hm\gg1$ (§13.14), the distance between the lids (§13.17), scale height c²/g (§13.3); h total depth $H+\eta$ over an uneven bottom (§13.13); η surface displacement (ch12: the Kolmogorov length; p. 688 uses it that way once) | $H$ depth or lid separation; $H_N$ for the WKB scale; $H_s$ for the scale height; $H_e$ equivalent depth; $h$ total depth; η displacement ($\eta_K$ if the Kolmogorov length is meant) | `H`, `H_N`, `H_s`, `He`, `h`, `eta` |
| c, $c_n$ | sound speed in $\frac{d\rho_\theta}{dz}=\frac{d\rho}{dz}+\frac{g\rho}{c^2}$ (13.1); long-wave speed $\sqrt{gH}$ from §13.8 on; modal speed $c_n$; complex phase speed of a normal mode in §13.16–13.17; $c_x$ zonal phase speed; $\mathbf c_g$ group velocity | $c_s$ sound speed in our own lines; $c$ long-wave speed; $c_n$ mode n; $c=c_r+ic_i$ only in C15's Rayleigh–Kuo note and C16, said each time | `c_sound`, `c`, `c_n`, complex `c` |
| k, l, m, K | Wavenumbers: k eastward, l northward, m vertical; $K=\sqrt{k^2+l^2}$ horizontal magnitude (in §13.14 **K** is the three-dimensional wavevector); K perturbation kinetic energy in $\frac{dK}{dt}=-g\int w'\rho'\,dx\,dy\,dz$ (§13.17); $K_0,K_1,K_2$ wavenumbers of Fjørtoft's triad; K(z) our eddy-viscosity profile | $k,l,m$; $K$ wavenumber magnitude; $K_E$ for the Eady kinetic energy in our own lines; $K_v(z)$ for the eddy-viscosity profile; l (p. 688) the eddy length is written $L_\beta$ | `k`, `l`, `m`, `K`, `KE`, `K` (array at faces in `ekman_solve`) |
| ζ | relative vorticity $\zeta\equiv\frac{\partial v}{\partial x}-\frac{\partial u}{\partial y}$; ch07: interface displacement; ch12: the stability parameter z/L_M | Relative vorticity only (ζ); $\bar\zeta$ of the basic flow, $\zeta_g$ geostrophic | `zeta` |
| θ | latitude (§13.4); angle of the wavevector with the horizontal in $\tan\theta=m/k$ (§13.14); potential temperature (in $\rho_\theta$, Exercise 13.7, ch01) | Latitude is θ; $\theta_K$ for the wavevector angle in our own lines (the book's θ quoted inside its equations with a ⚠️ line); $\theta_p$ potential temperature | `lat_rad`, `theta_K`, `theta_p` |
| δ | Ekman thickness $\delta=\sqrt{2\nu_v/f}$ (13.29) (an e-folding scale; the oceanographic Ekman depth is πδ); earlier: boundary-layer thickness, Kronecker delta, a small increment | Ekman e-folding thickness δ; $D_E=\pi\delta$ Ekman depth | `delta`, `ekman_depth(..., convention=)` |
| Λ, λ | Rossby radius Λ: $c/f$ external, $\sqrt{g'H}/f$ or $NH/(n\pi f)$ internal, and $NH/f$ **without π** in §13.17; ch12: integral scales $\Lambda_t,\Lambda_f$. Wavelength λ (ch12: Taylor microscale; elsewhere longitude or an eigenvalue) | $\Lambda$ = c/abs(f) with the wave speed named; $\Lambda_E=NH/f$ for the Eady radius; λ wavelength only (longitude is never used; the root of a characteristic equation is written $r$) | `rossby_radius*`; `lam` |
| ψ | vertical mode $\psi_n(z)$ (§13.9); stream function with $u'=-\partial\psi/\partial y$, $v'=\partial\psi/\partial x$ (§13.5, §13.16) — the opposite sign to ch04 and ch11 | $\psi_n$ modes; ψ stream function with the sign stated beside it | `psi` (modes), `psi_sf` |
| E, Ro, R | E Ekman number $\nu/(fL^2)$ (13.18) (earlier: energy); Ro Rossby number $U/(fL)$ (13.13) (ch04's `rossby_number` is U/(2ΩL)); R earth's radius (earlier: gas constant, a correlation) | $E$ Ekman number; energies are written out (KE, PE); Ro; $R$ = 6371 km | `ekman_number`, `rossby_number`, `EARTH_RADIUS_MEAN` |
| α | thermal expansion coefficient (§13.3, our temperature form of the thermal wind); scaled wavenumber $\alpha^2\equiv\frac{N^2}{f^2}(k^2+l^2)$ (13.139); enstrophy flux in $S\propto\alpha^{2/3}K^{-3}$ (§13.18) | $\alpha_T$ expansion in our own lines where two meet; α Eady wavenumber in C16; $\alpha_Z$ enstrophy flux in our own lines | `alpha` (expansion, keyword-only), `alphaH`, `alpha_ens` |
| U, V, τ, u | U velocity scale, geostrophic interior velocity, mean current, $U(z)=U_0z/H$; V complex velocity $u+iv$ (§13.6–13.7); τ wind stress (and $\tau_{ij}$); u on p. 688 an rms speed | $U$ each time named; $V$ complex velocity; $V_g$ its geostrophic part; τ stress; $u_{rms}$ | `U`, `U_g`, `V_g`, `tau_x`, `tau_y`, `u_rms` |
| primes, i, S | primes = perturbations (dropped by the book from the thin-shell equations of §13.4 on, trap T3; ch12: lower case = fluctuation); i = √−1 and the label in $T_i$; S(K) the energy spectrum with $\overline{u^2}=\int_0^\infty S\,dK$ (no ½, one-sided; ch12 used two-sided spectra with ∫S = variance) | as in the book, with a ⚠️ line where the convention changes | — |

5. **Normalisations.** Vertical modes: ψ_n(0) = 1 (free surface or lid at z = 0). Eady: $\Lambda_E=NH/f$, αH dimensionless, growth
   rate in units of fU₀/(NH). Spectra of §13.18: one-sided in K, $\overline{u^2}=\int_0^\infty S\,dK$ (trap T16). Kelvin and
   Poincaré amplitudes are surface-height amplitudes [m]. Rossby functions return signed ω for signed k.
6. **Sign and frame conventions.** x east, y north, z up. **z = 0 is**: the sea surface with the ocean in z < 0 (§13.6 surface
   Ekman layer, §13.9 modes, §13.14); the solid surface with the fluid in z > 0 (§13.7); the flat bottom (§13.8); the lower lid
   (§13.17) — every block opens by saying which (trap T4). **Every "to the right", "clockwise", "counter-clockwise round a
   low", "coast on the right" in the book assumes the northern hemisphere, f > 0.** Every such sentence in the notebook reads
   "in the northern hemisphere (f > 0); mirror for f < 0", and every explainer where direction matters (E1, E2, E3, E4, E7, E8)
   has a hemisphere control (N / S chips = sign of f). Code takes abs(f) and sign(f). A wind is named by where it comes from, a
   current by where it goes (`ch13.wind_from_to`). **Earth's rotation**: the sidereal `OMEGA_EARTH` = 7.292115 × 10⁻⁵ rad/s
   everywhere; the book's "2π per day" (a solar day) is 0.27 % smaller (*expect* 0.274 %; trap box in C01, R11).
7. **Lapse rate (the user's standing rule).** Compute with Kundu's $\Gamma\equiv dT/dz$ ($\Gamma_a=-g/C_p\approx-9.8$ K/km; stable when
   $dT/dz>\Gamma_a$) and **always** show the meteorological $\Gamma_{met}\equiv-dT/dz$ ($\Gamma_d\approx+9.8$ K/km; stable when
   $\Gamma_{met}<\Gamma_d$) beside it: a two-row table (`ch13.lapse_rate_table(dT_dz, Gamma_a=…)`), a worked conversion (negate the
   number, flip the inequality), both inequalities in every legend, badge and slider trace (F2). `Gamma_a` is always obtained
   from `STRAT.adiabatic_lapse_rate()` and passed by keyword (`ch12.gradient_richardson_thermal(…, Gamma_a=Gamma_a)`), never typed.
   The ch01 library verdict string writes "Γ < Γa" for the meteorological form: it is **never shown raw**; the notebook builds
   its own strings and explains once, in the §13.2 recap, that the library's "Γ < Γa" means $\Gamma_{met}<\Gamma_d$. N² in this
   chapter is always from **potential** density (or potential temperature), never from in-situ temperature alone.
8. **Book slips taught in corrected form** — house form `> ⚠️ **slip #k — the book prints** … **; the correct form is** …` where
   used ("the book prints", never "the page prints"); keys `ch13.book_slips()`.

| Slip | The book prints | Correct | Taught in | Planted wrong variant a test must fail |
|---|---|---|---|---|
| #1 | the momentum member of (13.2) as $\frac{D\mathbf u}{Dt}+2\boldsymbol\Omega\times\mathbf u=+\frac1{\rho_0}\nabla p-\frac{g\rho}{\rho_0}\mathbf e_z+\mathbf F$ | $\frac{D\mathbf u}{Dt}+2\boldsymbol\Omega\times\mathbf u=-\frac1{\rho_0}\nabla p-\frac{g\rho}{\rho_0}\mathbf e_z+\mathbf F$ | C01 (R04), D02 Start | `boussinesq_rotating_sympy(printed=True)`: the rest state leaves a residual |
| #2 | $F_i=\partial\tau_{ij}/\partial x_j$ as the friction force per unit mass, above the three components of (13.6) such as $F_x=\nu_H\big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\big)+\nu_v\frac{\partial^2u}{\partial z^2}$ | $F_i=\frac1\rho\,\partial\tau_{ij}/\partial x_j$ | C01 (R10), D01 | a `pint` dimension check of the printed form |
| #3 | (13.17) as $-2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$ | $2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$ | C03 (N21), D05 | `taylor_proudman_sympy(printed=True)` reaches $\partial u/\partial x-\partial v/\partial y=0$, not $\partial w/\partial z=0$ (13.19) |
| #4 | (13.26) as $u,v\to0$ for $z\to\infty$ | $u,v\to0$ for $z\to-\infty$ (the ocean is below the surface z = 0) | C04 (N28), D06 | the root bounded upward violates the surface-stress condition |
| #5 | after $T_i=2\pi/f$: the subscript "does refer to a component" | it does not; i is a label for "inertial" | C01 (N15) | — |
| #6 | the first root of $\tan\frac{NH}{c_n}=\frac{c_nN}g$ (13.69) is said to be at NH/c_n = 1 | $NH/c_0\ll1$ (ours, computed: 0.0558 for H = 4200 m, N = 2.7 × 10⁻³ s⁻¹) | C08 (R22) | `modes_uniform_N`: the first root asserted below 0.1 |
| #7 | the low-frequency paragraph closes with ω ≫ f as the range where the first term of $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76) is negligible | the range ω ≪ f | C09 (N83), D13, E6 toggle | `dispersion_term_sizes`: the ω³ term is the largest at ω = 3f and the smallest at ω = 0.1f |
| #8 | in §13.14, N is taken "depth independent" | depth dependent, N(z), as $m^2(z)\equiv\frac{(k^2+l^2)[N^2(z)-\omega^2]}{\omega^2-f^2}$ (13.99) requires | C14 (N105), D19 | — |
| #9 | cross-references "Section 4.18", "Section 4.11", "Section 8.7", "Section 5.7", and the number 7.128 for the definition of $N^2$ | §4.9 (Boussinesq), §4.9 (buoyancy), §8.4 (impulsively started plate), §5.6 (vortex stretching), and $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ (7.127) | conventions block; each recap cites our section with the equation written out | — |
| #10 | a laminar Ekman thickness for air about a quarter below $\sqrt{2\nu/f}$ with its own stated inputs | the value of the formula $\delta=\sqrt{2\nu_v/f}$ (13.29) | C06 (N48) | the private test asserts the formula and records the gap |
| #11 | a typical internal Rossby radius about three times $NH/(\pi f)$ with the chapter's own typical N, H, f | an inconsistency of inputs, not of the formula; both shown with our inputs | C12 (R24) | the private test records both |
| #12 | an exercise answer for the long Rossby-wave speed that follows only from a round β | with β at the stated latitude the speed is about a fifth smaller | C15 (N135), as a remark on rounding | the private test reproduces the printed value only with the round β |
| #13 | the 1940s theory of baroclinic instability attributed to V. Bjerknes and others | **not asserted in the notebook until a first-hand source is read**; C16 says only "the problem solved here is known as the Eady problem" | C16 opening | — |

9. **Traps (true as printed, easy to misread) — each a `> ⚠️ Common confusion:` callout where it first bites.**

| Trap | What is easy to misread | Callout placed in |
|---|---|---|
| T1 | the book's Ω is 2π per solar day; the sidereal value is 0.27 % larger | C01 (R11, P308) |
| T2 | "close to neutral" compares with the moist adiabat; the standard troposphere is stable to dry displacements; N² is from potential density | C01, §13.2 recap (R02, N05, R03) |
| T3 | from the thin-shell equations on, p and ρ are perturbations with the primes dropped | C01 (R08), repeated in C03 |
| T4 | where z = 0 is changes from section to section | conventions block; first line of C04, C06, C07, C08, C16 |
| T5 | every formula with √f, "to the right", "clockwise", $e^{-fy/c}$ assumes f > 0 | conventions block (N02), C04 (D06), C11 (D15), every hemisphere control |
| T6 | the thickness δ is an e-folding scale; the "Ekman depth" is πδ | C04 (N30) |
| T7 | the book's "orbit" figure is a velocity hodograph; axis ratio quoted as ω/f in §13.11 and f/ω in §13.14 | C10 (N89), C14 (N120) |
| T8 | with ω > 0 a Rossby wave has k < 0; the "maximum" phase speed is a maximum of magnitude; group-velocity arrows point into the circles | C15 (D23, N132) |
| T9 | "all roots real, two superinertial" holds under β-plane scaling; the fast roots have opposite signs and only satisfy abs(ω) > abs(f) | C09 (N82, D13) |
| T10 | the parameter f is replaced by f₀ except where it is differentiated; potential-vorticity conservation itself needs no such step | C09 (N81), C13 (D17, R26) |
| T11 | three Rossby radii: external, internal (with π), and the Eady radius NH/f (no π) | C12 (R24), C16 (N161) |
| T12 | the modal amplitudes have different units: u_n, v_n [m/s], p_n [m²/s²], w_n [1/s], ρ_n [kg/m²] | C08 (N58) |
| T13 | overloaded letters (the table of convention 4) | conventions block; repeated at C14 (θ), C16 (α), C17 (α, η, l) |
| T14 | here u = −∂ψ/∂y, v = ∂ψ/∂x; ch04 and ch11 use the opposite sign | C02 (N18), C15 (R30) |
| T15 | approximations that enter silently (continuity and uniform coefficients in the friction force; uniform density in the no-shear step; no vertical advection of momentum in the Eady set; √m frozen in the WKB velocity; the lee-wave k is a magnitude) | D01, D05, D25, C14 (N116, N125) |
| T16 | §13.18's spectrum is one-sided with no factor ½; ch12's was two-sided | C17 (N164) |
| T17 | two Rossby numbers: U/(fL) here, U/(2ΩL) in ch04 — they differ by sin θ | C02 (R13) |

10. **Colour code (curation §5; the same in figures, derivations, bars and explainers):** pressure-gradient force = orange ·
    Coriolis force = teal · friction = rose · acceleration / tendency = purple (accent) · density and temperature = blue ·
    vorticity and potential vorticity = amber · references and ghosts = muted.
11. **Honest labels.** (a) Closed forms are "analytic"; the C-grid shallow-water model, the exact-spectral QG evolution used on
    a periodic box, the forward–backward 1-D march and the pseudo-spectral barotropic model are **our numerical choices** — each
    gets a `> 🔧 **Our choice** — …` box where first used (P334, P335). (b) **"Not in the book — ours"** opens every appearance of:
    the thermal wind in temperature form (N20), Ekman pumping and the one line on Sverdrup balance (N36), the Ekman layer with
    $K_v(z)$ (N37), geostrophic adjustment (N19, D16), the Eady growth rate in physical units (N162), the two-dimensional
    turbulence run (N168); also the step-flow streamline solution (D18), the Poincaré and Rossby group-velocity components, the
    finite-depth Ekman solution and the third branch on the Poincaré–Kelvin figure. (c) **"Ours, computed" with no citation**:
    the attribution of slip #13, the Eady coefficient 0.30982, the −3 range, the Rhines scale, the closed form of the adjusted
    step and its energy ratio 1/3, and the C-grid scheme's attribution. No paper or textbook is cited for these anywhere.
12. **Public-repo rule.** No book number is reproduced: worked examples use `ch13.illustrative_inputs()` (35° N and 35° S; 60° N
    for the Ekman layers; 12° N for Rossby waves; ocean 4200 m, N = 2.7 × 10⁻³ s⁻¹; atmosphere 9 km, N = 1.1 × 10⁻² s⁻¹, U₀ = 27 m/s;
    τ = 0.07 N/m², ν_v = 0.03 m²/s; U_g = 12 m/s, ν_v = 7 m²/s; U = 14 m/s, L = 1400 km; two-layer 120 m, Δρ = 3.1 kg/m³; mean
    current 17 m/s; wavelengths 730 km and 3100 km; inertial current 0.23 m/s; Fjørtoft pair K₁ = K₀/3, K₂ = 3K₀/2; rms speeds
    13 m/s and 0.08 m/s). The Eady constants are quoted to five digits (2.3994, 1.6061, 2.6187, 3.9120, 0.30982) as ours,
    computed, never in a two-digit rounding beside Λ or k. Observed figures (the Oregon current profile, the mid-troposphere
    height map) are not reproduced; our model output replaces them.
13. **Sentences and table cells never start with a symbol** (a capitalisation pass once turned κ into Κ): write "The thickness δ
    …", not "δ …".

---

## Part C — functions the builders will call (the implementer's contract)

Names and signatures are those of `analysis/ch13.md` §4 (rows I01–I29) and curation §8 (items 1–11) and §9. "Eq." = the book
equation implemented (written out in the docstring). Items marked **(+)** are not in analysis §4 (C.5 lists them). Return shapes
are part of the contract (convention 2): scalars in → Python floats (or complex) out; a `dict` has exactly the keys listed; a
tuple has exactly the members listed, in that order. Everything after a bare `*` is keyword-only. `f` may have either sign in
every function unless the row says "f > 0 only"; `s = sign(f)`; a function with `f` in a denominator raises `ValueError` at
`f == 0` (never returns inf or NaN). Latitudes in radians (`lat_rad`). Grids are `[j, i]` = (y, x). The default of g is
`G0 = 9.80665` m/s² (`core.thermo.G0`).

### C.0 Reused (existing; called, not changed)
`ROT.coriolis_parameter(lat_rad, Omega=OMEGA_EARTH)` (ch04; $f=2\Omega\sin\theta$ (13.8)), `ROT.OMEGA_EARTH = 7.292115e-5` rad/s
(sidereal), `ROT.coriolis_acceleration(Omega, u)` (the full $2\boldsymbol\Omega\times\mathbf u$) · `SIM.rossby_number(U, Omega, l,
factor=2.0)` (ch04's U/(2Ωl)) · `STRAT.lapse_rate_stability(dT_dz, Gamma_a=None, convention="kundu" or "meteorology")` (returns
`.verdict`, `.text`), `STRAT.adiabatic_lapse_rate(T=None, cp=CP_AIR, alpha=None, g=G0)`, `STRAT.lapse_rate_convention`,
`STRAT.brunt_vaisala_sq(rho0, drho_dz, drho_a_dz, g)`, `STRAT.brunt_vaisala_sq_from_theta(theta, dtheta_dz, g)`,
`STRAT.brunt_vaisala_sq_from_lapse(T, dT_dz, cp, g)`, `STRAT.ocean_potential_density_gradient(drho_dz, rho, c, g)`,
`STRAT.isentropic_density_gradient(rho, c, g)`, `STRAT.stability_timescale(N2)` (ch01) · `core.statics` USSA-1976 temperature
profile (ch01) · `ch01.seawater_density_linear` · `ch05.column_relative_vorticity(h, h0, zeta0=0.0, f=1e-4)` ·
`WAV.internal_wave_omega(k, m, N, l=0.0)`, `WAV.group_velocity_vector(omega_fn, K)`, `WAV.doppler_frequency(omega, U, K)`,
`WAV.two_layer_long_wave_speed(H, rho1, rho2, g)`, `WAV.reduced_gravity_book(rho1, rho2, g, ref="lower")` (ch07) ·
`LAM.stokes_layer(nu, omega)` (ch08; the same $e^{-(1+i)y/\delta}$ structure) · `ST.rayleigh_eigs_contour(k, U, Up, Upp, domain, N,
…)`, `ST.cheb_grid`, `ST.generalized_eigs` (ch11), `ch11.rayleigh_criterion(y, U=None, Upp=None)`, `ch11.gradient_richardson` ·
`ch12.gradient_richardson_thermal(dTdz, dUdz, alpha, g=G0, *, Gamma_a, convention="kundu")`, `TS.shell_spectrum(components, L,
n_bins=None)`, `TS.synthetic_solenoidal_field(n, L, spectrum, seed=0, dim=2)` (ch12) · `DIM.pi_groups` (ch01) ·
`core.operators.gradient`, `curl` · `core.potential` (cylinder flow for the Taylor-column sketch) · `core.mac.MacGrid` (layout
only) · machinery: `setup_notebook`, `show_viz`, `animate`, `show_animation`, `slider_figure`, `animate_figure`, `live`,
`savefig`, `COLORS`.

### C.1 `fluidpy/core/gfd.py` (NEW, `GFD`; re-exported by `ch13`) — rotation, balances, Ekman layers, waves, instability

**C.1a Constants and geometry**

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `OMEGA_EARTH` (re-export of `ROT.OMEGA_EARTH`), `OMEGA_SOLAR_DAY = 2*np.pi/86400.0`, `EARTH_RADIUS_MEAN = 6.371e6`, `SIDEREAL_DAY = 2*np.pi/OMEGA_EARTH` **(+)** | floats [rad/s], [rad/s], [m], [s] | R11, trap T1; `OMEGA_SOLAR_DAY` is for the private tests and the trap box only |
| `coriolis_parameter(lat_rad, Omega=OMEGA_EARTH)` (re-export of `ROT.coriolis_parameter`) | value f [1/s], signed | $f=2\Omega\sin\theta$ (13.8) |
| `hemisphere(lat_rad)` | dict: "sign" (+1.0, −1.0 or 0.0), "name" ("northern", "southern", "equator"), "turns" ("right", "left", "none"), "cyclonic" ("counter-clockwise", "clockwise", "none") | N02, trap T5 |
| `earth_rotation_local(lat_rad, Omega=OMEGA_EARTH)` | array (3,) = (0, Ω cos θ, Ω sin θ) [rad/s] | N13, $\boldsymbol\Omega=(0,\Omega\cos\theta,\Omega\sin\theta)$ |
| `inertial_period(f)` | $2\pi/\lvert f\rvert$ [s]; `np.inf` at f = 0 (the one function that returns inf: the period really is infinite) | N15, $T_i=2\pi/f$ |
| `f_plane(lat0_rad, Omega=OMEGA_EARTH)` | value f₀ [1/s] | N16 |
| `beta_parameter(lat_rad, Omega=OMEGA_EARTH, R=EARTH_RADIUS_MEAN)` | value β = 2Ω cos θ₀ / R [1/(m s)] | $f=f_0+\beta y$ (13.10) |
| `beta_plane(y, lat0_rad, Omega=OMEGA_EARTH, R=EARTH_RADIUS_MEAN)` | value f₀ + βy [1/s] | (13.10) |
| `beta_plane_error(y, lat0_rad, Omega=OMEGA_EARTH, R=EARTH_RADIUS_MEAN)` | relative error (f_β − f_exact)/f_exact with f_exact = 2Ω sin(θ₀ + y/R) [–] | N17 (ours) |
| `vertical_velocity_scale(U, H, L)` | W = UH/L [m/s] | N12, $W/U\sim H/L$ |
| `coriolis_acceleration_local(u, v, w, lat_rad, Omega=OMEGA_EARTH, thin=True)` | tuple (a_x, a_y, a_z) [m/s²] of $2\boldsymbol\Omega\times\mathbf u$: `thin=True` → (−fv, fu, −2Ωu cos θ); `thin=False` → (2Ω(w cos θ − v sin θ), 2Ωu sin θ, −2Ωu cos θ) | $2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)$ (13.7) |
| `thin_layer_terms(U, L, H, lat_rad, *, drho_over_rho0=1e-3, nu_H=0.0, nu_v=0.0, g=G0, Omega=OMEGA_EARTH)` | dict: "W", "aspect", "f", "Ro", "x" → dict("acceleration" U²/L, "coriolis" abs(f)U, "coriolis_w" 2Ω cos θ·W, "friction_H" ν_H U/L², "friction_v" ν_v U/H², "pressure" (= the largest of the others)), "z" → dict("acceleration" UW/L, "coriolis" 2ΩU cos θ, "buoyancy" g·drho_over_rho0, "friction_H", "friction_v") [m/s²] | N12, N14, C01 (scale sizes of every term of (13.9)) |
| `eddy_friction(lap_h, d2z, nu_H, nu_v)` | value ν_H·lap_h + ν_v·d2z [m/s²] | (13.6), one component |
| `scale_height(c, g=G0)` | value c²/g [m] | R05 |
| `buoyancy_frequency_sq(z, rho, rho0, g=G0)` | array N² = −(g/ρ₀) dρ/dz [1/s²] (second-order stencils, one-sided at the ends) | R03, $N^2\equiv-\frac{g}{\rho_0}\frac{d\rho}{dz}$ |

**C.1b Balances**

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `geostrophic_velocity(dpdx, dpdy, f, rho0)` | tuple (u, v) = (−dpdy/(ρ₀f), dpdx/(ρ₀f)) [m/s] | $-fv=-\frac1{\rho_0}\frac{\partial p}{\partial x}$, $fu=-\frac1{\rho_0}\frac{\partial p}{\partial y}$ (13.11)–(13.12) |
| `geostrophic_from_field(p, dx, dy, f, rho0)` | tuple (u, v) arrays `[j, i]` (centred differences, one-sided at edges) | (13.11)–(13.12) |
| `geostrophic_from_height(eta, dx, dy, f, g=G0)` | tuple (u, v) arrays | (13.116) $u\simeq-\frac g{f_0}\frac{\partial\eta}{\partial y}$, $v\simeq\frac g{f_0}\frac{\partial\eta}{\partial x}$ |
| `geostrophic_streamfunction(p, f, rho0)` | value ψ = p/(fρ₀) [m²/s] (u = −∂ψ/∂y, v = ∂ψ/∂x) | N18, trap T14 |
| `rossby_number(U, f, L)` | U/(abs(f)L) [–] | $\mathrm{Ro}=\frac{U^2/L}{fU}=\frac U{fL}$ (13.13) |
| `ekman_number(nu, f, L)` | value ν/(abs(f)L²) [–] | $E=\frac{\rho\nu U/L^2}{\rho fU}=\frac\nu{fL^2}$ (13.18) |
| `thermal_wind_shear(drho_dx, drho_dy, f, rho0, g=G0)` | tuple (du_dz, dv_dz) = (g·drho_dy/(ρ₀f), −g·drho_dx/(ρ₀f)) [1/s] | (13.15) |
| `thermal_wind_from_temperature(dTdx, dTdy, f, *, alpha, g=G0)` | tuple (du_dz, dv_dz) = (−gα·dTdy/f, gα·dTdx/f) [1/s] | N20 (ours): (13.15) with ρ′ = −ρ₀αT′ |
| `thermal_wind_integrate(z, drho_dy, f, rho0, g=G0, u_ref=0.0)` | array u(z) [m/s] by cumulative trapezoid from z[0], where u = u_ref (`drho_dy` scalar or array on z) | (13.15b), (13.128) |
| `taylor_proudman_residual(u_fn, x, h=1e-4)` | dict: "du_dz", "dv_dz", "dw_dz" at the point x = (x, y, z) by central differences of the callable `u_fn(x, y, z) -> (u, v, w)` | $\partial\mathbf u/\partial z=0$ (13.21) |
| `parcel_adjust(t, G, f, r=0.0, V0=0j)` **(+)** | complex V = u + iv [m/s] of a parcel under the uniform complex pressure-gradient acceleration G = (1/ρ₀)(∂p/∂x + i ∂p/∂y) [m/s²] with linear drag r [1/s]: V = V∞ + (V0 − V∞)e^{−(r+if)t}, V∞ = −G/(r + if) | curation §8 note 1: closed form of dV/dt + (r + if)V = −G (ours; E1's stage) |

**C.1c Ekman layers** (surface layer: z ≤ 0, z = 0 the sea surface; bottom layer: z ≥ 0 from the solid surface; δ = √(2ν_v/abs(f)))

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `ekman_depth(nu_v, f, convention="efold")` | value δ [m]; `convention="pi"` → πδ | $\delta=\sqrt{2\nu_v/f}$ (13.29); trap T6 |
| `ekman_surface(z, tau_x, tau_y, rho, nu_v, f, U_g=0.0, V_g=0.0, as_complex=False)` | tuple (u, v) [m/s] (or complex V if `as_complex`): V = V_g + (τ_x + iτ_y)(1 − is)/(ρ√(2ν_v abs(f)))·e^{(1+is)z/δ} | C04; (13.27)–(13.29) with both hemispheres |
| `ekman_transport(tau_x, tau_y, rho, f)` | tuple (M_x, M_y) = (τ_y/(ρf), −τ_x/(ρf)) [m²/s] | $\int_{-\infty}^0u\,dz=0$, $\int_{-\infty}^0v\,dz=-\frac\tau{\rho f}$ (13.30) |
| `ekman_transport_partial(z, tau_x, tau_y, rho, nu_v, f)` **(+)** | tuple (M_x(z), M_y(z)): transport between depth z and the surface, closed form M(z) = M(−∞)·(1 − e^{(1+is)z/δ}) in complex form | D07 (ours); E3's running transport |
| `ekman_residual(z, u, v, nu_v, f, U_g=0.0, V_g=0.0)` | tuple (r_x, r_y) arrays: −f(v − V_g) − ν_v u″ and f(u − U_g) − ν_v v″ on interior nodes [m/s²] | (13.22)–(13.23), (13.33)–(13.34) |
| `ekman_vorticity_balance(z, tau_x, rho, nu_v, f)` | dict: "tilt_x" (= −f dv/dz), "diff_x" (= ν_v d²ω_y/dz²), "tilt_y" (= −f du/dz), "diff_y" (= ν_v d²ω_x/dz²), "omega_x", "omega_y" (closed form) | (13.31) |
| `ekman_solve(z, K, f, *, tau=None, rho=None, V_g=0.0, bottom="decay", top="stress")` | complex array V at the nodes z (monotonic, non-uniform allowed); K at the faces (array of len(z) − 1, or a callable K(z_face)); `top="stress"` needs `tau` (complex τ_x + iτ_y) and `rho`; `top="geostrophic"` sets V = V_g at the top; `bottom="decay"` sets V = V_g, `"noslip"` sets V = 0 | N37 (ours): $\frac{d}{dz}\big(K\frac{dV}{dz}\big)=if(V-V_g)$ |
| `ekman_pumping(tau_x, tau_y, dx, dy, rho, f)` | array w_E `[j, i]` = (1/ρ)[∂(τ_y/f)/∂x − ∂(τ_x/f)/∂y] [m/s] (`f` scalar or array) | N36 (ours), D09 |
| `ekman_pumping_bottom(zeta_g, nu_v, f)` | value w = s·(δ/2)·ζ_g [m/s] (positive = upward out of the layer) | N36 (ours), D09 |
| `ekman_finite_depth(z, tau_x, tau_y, U_g, V_g, H, nu_v, rho, f)` | complex array V(z), −H ≤ z ≤ 0: surface stress at z = 0, no slip at z = −H, exact (two exponentials) | N34 (ours) |
| `ekman_bottom(z, U_g, V_g, nu_v, f, as_complex=False)` | tuple (u, v) (or complex): V = V_g(1 − e^{−(1+is)z/δ}), V_g = U_g + iV_g | (13.41) |
| `ekman_bottom_transport(U_g, V_g, nu_v, f)` | tuple (M_x, M_y) = the integral of (V − V_g) over 0 ≤ z < ∞ = −V_g δ/(1 + is) [m²/s]; for f > 0, V_g = (U, 0): (−½Uδ, +½Uδ) | N47, $\int_0^\infty v\,dz=\tfrac12U\delta$ |
| `ekman_force_balance(z, U_g, nu_v, f, rho=1.0)` | dict of 2-tuples (x, y) [N/m³ for the given rho; m/s² with rho = 1]: "coriolis" = ρf(v, −u), "pressure" = ρf(0, U_g), "friction" = ρν_v(u″, v″); plus "u", "v", "speed", "angle_to_isobars" [rad, ≥ 0], "sum" (2-tuple, round-off) | N49; (13.33)–(13.34) as three vectors |
| `eddy_viscosity_from_depth(delta, f)` | value ν_v = abs(f)δ²/2 [m²/s] | N48 |

**C.1d Shallow-water waves and the Rossby radius**

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `long_wave_speed(H, g=G0)` | √(gH) [m/s] | (13.86) |
| `equivalent_depth(c, g=G0)` | value c²/g [m] | $c^2=gH_e$ (13.46), (13.62) |
| `baroclinic_mode_speed(N, H, n=1)` | NH/(nπ) [m/s] | $c_n=\frac{NH}{n\pi}$ (13.71) |
| `shallow_water_omega(k, l, c, f0, beta)` | tuple (ω₋, ω_R, ω₊) in ascending order [rad/s]: the three real roots of the cubic; **(nan, nan, nan) where `shallow_water_discriminant` is negative** (the one documented exception to "never NaN": callers test the discriminant first) | $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76) |
| `shallow_water_branches(k, l, c, f0, beta)` **(+)** | dict: "poincare_plus", "poincare_minus", "rossby" (the three roots above), "kelvin" (= c·k) | curation §8 note 5 |
| `shallow_water_discriminant(k, l, c, f0, beta)` | 4(c²K² + f₀²)³ − 27(c²βk)² [1/s⁶] | N82 |
| `shallow_water_regime(omega, f0)` | dict: "ratio" abs(ω/f₀), "regime" ("high" for ratio > 3, "near-inertial" for 1 ≤ ratio ≤ 3, "low" for ratio < 1), "neglect" (the term of the cubic that may be dropped: "f0^2 omega", "c^2 beta k", "omega^3") | N83 |
| `dispersion_term_sizes(k, l, c, f0, beta, omega)` | dict: "omega3" = ω³, "gravity" = −c²K²ω, "rotation" = −f₀²ω, "beta" = −c²βk, "sum" | (13.76) term by term |
| `poincare_omega(K, f, c)` | +√(f² + c²K²) [rad/s] | $\omega^2=f^2+gHK^2$ (13.82) |
| `poincare_group_velocity(k, l, f, c)` | tuple (c_gx, c_gy) = c²(k, l)/ω [m/s] | D14 (ours) |
| `poincare_amplitudes(k, l, omega, f, g, eta_hat)` | tuple (û, v̂) complex [m/s] | (13.80) |
| `poincare_fields(x, y, t, k, l, eta_hat, H, f, g=G0)` | tuple (η, u, v) real fields; ω = +`poincare_omega` | (13.80), (13.83) |
| `poincare_orbit(t, k, eta_hat, H, f, g=G0)` | dict: "u", "v" (velocity at x = 0), "x", "y" (particle displacement about its mean position), "omega", "axis_ratio" (= ω/abs(f)), "sense" ("clockwise" for f > 0, "counter-clockwise" for f < 0) | (13.83), N89 |
| `inertial_oscillation(t, u0, v0, f)` | tuple (u, v, x, y): u + iv = (u0 + iv0)e^{−ift}; x + iy = i(u + iv − u0 − iv0)/f | N90 |
| `inertial_radius(q, f)` | value q/abs(f) [m] | N90 |
| `kelvin_wave(x, y, t, eta0, k, H, f, g=G0, direction=None)` | tuple (η, u): η = η₀e^{−y/Λ}cos k(x − d·ct), u = d·η₀√(g/H)e^{−y/Λ}cos k(x − d·ct), fluid in y ≥ 0, d = `direction` (+1 toward +x, −1 toward −x; default sign(f), the trapped one); raises `ValueError` if d·sign(f) < 0 (the solution would grow offshore) | (13.87) |
| `kelvin_omega(k, c)` **(+)** | value c·k [rad/s] | (13.86) |
| `kelvin_decay_side(f, direction)` **(+)** | dict: "trapped" (bool: direction·sign(f) > 0 for fluid in y ≥ 0), "coast_on" ("right" or "left" of the direction of travel for a trapped wave in this hemisphere) | curation §8 note 6 |
| `kelvin_residuals(x, y, t, eta0, k, H, f, g=G0, direction=None, h=1.0, ht=1.0)` | dict: "continuity", "x_momentum", "y_geostrophy" — residuals of the three equations by central differences, each divided by the size of its largest term | (13.84) |
| `rossby_radius(c, f)` | value Λ = c/abs(f) [m] | C12, $\Lambda\equiv c/f$ |
| `rossby_radius_internal(N, H, f, n=1, with_pi=True)` | NH/(nπ abs(f)); `with_pi=False` → NH/abs(f) (the Eady Λ; n ignored) [m] | R24, N161; trap T11 |
| `rossby_radius_two_layer(H1, rho1, rho2, f, g=G0)` | √(g′H₁)/abs(f), g′ = g(ρ₂ − ρ₁)/ρ₂ [m] | R24 |
| `geostrophic_adjustment_1d(x, eta0, H, f, g=G0)` | tuple (η, v): end state of the initial step η₀ sgn(x): η = η₀ sgn(x)(1 − e^{−abs(x)/Λ}), v = s·(gη₀/c)e^{−abs(x)/Λ}; at f = 0 raises (no steady end state) | N19 (ours), D16 |
| `adjustment_energy(eta0, H, f, g=G0, L=None, rho=1000.0)` **(+)** | dict per unit length of the step [J/m]: "pe_released" (= (3/2)ρgη₀²Λ for `L=None`), "ke_jet" (= (1/2)ρgη₀²Λ), "radiated" (difference), "ratio" (ke_jet/pe_released); a finite `L` integrates over abs(x) ≤ L instead (closed form) | curation §8 note 7 (ours, computed) |

**C.1e Internal waves with rotation**

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `inertia_gravity_m2(k, l, omega, N, f)` | value m² = (k² + l²)(N² − ω²)/(ω² − f²) [1/m²] (negative outside the band; raises at ω² = f²) | (13.99) |
| `inertia_gravity_band(omega, N, f)` | dict: "propagating" (bool: abs(f) < ω < N), "where" ("below f", "in band", "above N") | N110 |
| `inertia_gravity_omega(k, m, N, f, l=0.0, approx=None)` | value ω ≥ 0 [rad/s]: ω² = (N²K_h² + f²m²)/(K_h² + m²); `approx` in (None, "nonrotating", "hydrostatic", "midrange") gives N K_h/√(K_h² + m²), √(f² + N²K_h²/m²), N K_h/abs(m) | (13.112), N123 |
| `inertia_gravity_regime(omega, N, f)` | dict: "regime" ("near-inertial", "mid", "near-buoyancy"), "err_nonrotating", "err_hydrostatic", "err_midrange" (relative error in m² of each approximation at this ω) | N123 |
| `inertia_gravity_group_velocity(k, m, N, f)` | tuple (c_gx, c_gz) [m/s] | N124, D21 |
| `wkb_vertical_structure(z, m, A0=1.0, sign=+1)` | complex array ŵ = A₀ m^{−1/2} exp(sign·i∫m dz) (phase by cumulative trapezoid from z[0]) | (13.104) |
| `inertia_gravity_fields(x, z, t, k, omega, N_fn, f, A0=1.0, sign=+1)` | tuple (u, v, w) real arrays on the grid (z, x) | (13.108) |
| `inertia_gravity_hodograph(t, omega, f, sign=+1)` | tuple (u, v) = (∓cos ωt, ±(f/ω) sin ωt), upper signs for `sign=+1` | (13.110) |
| `lee_wave_m(U, N, k)` | value m = √(N²/U² − k²) [1/m]; raises `ValueError("evanescent: k > N/U")` otherwise | (13.113) with ω = kU |

**C.1f Rossby waves, potential vorticity, instability, turbulence**

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `rossby_omega(k, l, beta, f0=0.0, c=np.inf, U=0.0)` | value ω = Uk − βk/(k² + l² + f₀²/c²) [rad/s], signed (k < 0 gives ω > 0 when U = 0) | (13.118), (13.120) |
| `rossby_phase_speed(k, l, beta, f0=0.0, c=np.inf, U=0.0)` | value c_x = ω/k [m/s] | (13.119), (13.120) |
| `rossby_group_velocity(k, l, beta, f0=0.0, c=np.inf, U=0.0)` | tuple (c_gx, c_gy) = (U + β(k² − l² − F)/(k² + l² + F)², 2βkl/(k² + l² + F)²), F = f₀²/c² | N132, D23 (ours) |
| `rossby_max_frequency(beta, f0, c)` | dict: "omega_max" = βc/(2abs(f₀)), "k" = −abs(f₀)/c (l = 0) | N132 |
| `rossby_omega_circle(omega, beta, f0, c)` | dict: "center_k" = −β/(2ω), "center_l" = 0.0, "radius" (√ of (β/2ω)² − f₀²/c²), "exists" (bool: the radicand is positive) | N132 |
| `rossby_long_wave_speed(beta, f0, c)` | −βc²/f₀² [m/s] | N135 |
| `stationary_rossby_wavelength(U, beta)` | 2π√(U/β) [m] (U > 0) | (13.120) |
| `rossby_packet_spectrum(k0, sigma, n)` **(+)** | tuple (k, amplitudes): n wavenumbers evenly spaced over k0 ± 3σ with Gaussian amplitudes that sum to 1 | curation §8 note 8 |
| `potential_vorticity(zeta, f, h)` | (ζ + f)/h [1/(m s)] | (13.94) |
| `step_vorticity(f, h0, h1)` | value ζ = f(h₁ − h₀)/h₀ [1/s] | N102 |
| `absolute_vorticity_gradient(y, U, beta)` | array β − U″(y) (U samples on y, second-order stencils) | (13.124) |
| `rayleigh_kuo_criterion(y, U, beta, tol=1e-12)` | dict: "changes_sign" (bool), "y_zero" (list of zero crossings), "min", "max" of β − U″ | (13.124) |
| `eady_alpha(k, l, N, f)` | value α = N√(k² + l²)/abs(f) [1/m] | (13.139) |
| `eady_factors(alphaH)` **(+)** | dict: "tanh_factor" = x − tanh x, "coth_factor" = x − coth x, "product", x = αH/2 (series for x < 1e-3) | radicand of (13.141) |
| `eady_phase_speed(alphaH, U0)` | complex c [m/s]: U₀/2 + (U₀/αH)√(product); the root with Im c ≥ 0 (for a neutral pair, the faster one) | (13.141) |
| `eady_growth_rate(k, l, N, f, H, U0)` | value σ = k·Im c [1/s] (0 beyond the cut-off) | N161 |
| `eady_critical()` | value α_cH [–], root of x = coth x doubled (`brentq`); ours, computed: 2.3994 | N161, D27 |
| `eady_fastest()` | dict: "alphaH" (1.6061), "sigma_nd" (= σNH/(abs(f)U₀) = 0.30982), "ci_over_U0", "cr_over_U0" (0.5), for l = 0 (bounded Brent maximisation) | N161, N162 (ours, computed) |
| `eady_max_growth_rate(f, N, dUdz)` | value σ_max = `eady_fastest()["sigma_nd"]`·abs(f)·dUdz/N [1/s] | N162 (ours) |
| `eady_time_scale(f, N, dUdz)` | 1/σ_max [s] | N162 (ours) |
| `eady_mode(z, alphaH, U0, H)` | dict: "p_hat" (complex array, normalised to max abs = 1), "amplitude", "phase" [rad] (phase of p̂ against height; it increases with height for the growing mode: a westward tilt), "A", "B" | (13.140) |
| `eady_vertical_velocity(z, alphaH, U0, H, k, N, f, rho0)` | complex array ŵ(z) for the mode above with unit pressure amplitude | (13.135) |
| `eady_fluxes(z, alphaH, U0, H, N, f, rho0, g=G0, k=None)` | dict: "w_rho" (mean of w′ρ′ over a wavelength, array on z), "v_rho", "v_T_sign" (+1 poleward for the growing mode), "phase_tilt" [rad] between top and bottom | N163 |
| `fjortoft_transfer(K0, K1, K2, S0=1.0)` | dict: "S1", "S2", "energy_ratio" (S₁/S₂), "enstrophy_ratio" (K₁²S₁/(K₂²S₂)) | (13.145) |
| `enstrophy_spectrum(K, S)` | K²S | N164 |
| `two_d_cascade_spectrum(K, K0, eps, alpha_ens, C_E=1.0, C_Z=1.0)` | array S(K): C_E ε^{2/3}K^{−5/3} for K < K₀, C_Z α^{2/3}K^{−3} for K > K₀ (shape only; constants are inputs) | N166 |
| `rhines_length(u_rms, beta)` | √(u/β) [m] | N167 |

### C.2 `fluidpy/core/vertical_modes.py` (NEW, `VM`; re-exported by `ch13`)

`Modes` is a `typing.NamedTuple` with fields `(c, psi, He, z, weights)`: `c` array (n_modes,) [m/s] in decreasing order (index 0
= barotropic with a free surface; with a rigid lid index 0 is the first baroclinic mode); `psi` array (n_modes, len(z)), each
normalised to ψ(z = 0) = 1; `He = c²/g`; `weights` the quadrature weights on z. Integer indexing works for parity rows:
`ch13.vertical_modes(...)[0][1]` is c₁ with a free surface. **With `lid="rigid"` there is no barotropic entry: index 0 is
the first baroclinic mode, so mode n is at index n − 1** (measured on the implementer's code: rigid c = 3.6096, 1.8048,
1.2032 m/s; free c = 203.05, 3.6085, 1.8047 m/s for uniform N = 2.7 × 10⁻³ s⁻¹, H = 4200 m, 401 nodes).

| Function (signature) | Returns [units] | Book eq. / source |
|---|---|---|
| `vertical_modes(z, N2, g=G0, n_modes=4, lid="free", method="fd")` | `Modes`; z nodes from −H to 0 (increasing), N2 at the nodes (> 0); `lid` in ("free", "rigid"); `method` in ("fd", "cheb") | (13.56) with (13.64), (13.65) |
| `vertical_modes_shooting(z, N2_fn, g=G0, n=1, lid="free")` **(+)** | tuple (c_n, psi_n): integrate from the bottom, root-find the surface condition (the method E5's JS uses) | curation §8 note 4 |
| `modes_uniform_N(N, H, g=G0, n_modes=4, lid="free", nz=201)` | `Modes` with exact speeds (`brentq` on tan(NH/c) = cN/g, one root per branch) and closed-form ψ_n on a uniform z grid of nz points | (13.69), (13.70), (13.71) |
| `rigid_lid_error(N, H, g=G0, n=1)` | (c_n(free) − c_n(rigid))/c_n(rigid) [–] | N75 |
| `w_structure(modes, n)` | array: the integral of ψ_n from −H to z | (13.53) |
| `rho_structure(modes, n)` | array dψ_n/dz | (13.54) |
| `project(modes, profile)` | array of coefficients a_n = ∫profile·ψ_n dz / ∫ψ_n² dz (weight 1 for both lids); with rigid-lid modes the depth mean of the profile is not represented (a constant projects to zero) and must be added back by the caller | N56, N60 |
| `reconstruct(modes, coeffs)` | array Σ a_nψ_n(z) | (13.52) |
| `orthogonality_matrix(modes, kind="psi", N2=None, g=G0, lid="free")` | array (n, n). `kind="psi"` (default; the orchestrator's "weight-1" form): ∫ψ_mψ_n dz/√(∫ψ_m²∫ψ_n²), weight 1 and **no surface term for either lid**. `kind="energy"` (needs `N2`): the second relation, ∫ψ_m′ψ_n′/N² dz + ψ_m(0)ψ_n(0)/g for a free surface (the surface term belongs here), normalised likewise | N60 |
| `modal_amplitudes(c_n, p_n_t, p_n=None, rho0=1.0, g=G0)` | dict: "w_n" = p_n_t/c_n²; "rho_n" = −(ρ₀/g)p_n when `p_n` is given, otherwise `None` | (13.60), (13.61) |
| `wkb_mode_speed(z, N, n=1)` | ∫N dz/(nπ) [m/s] | ours (the WKB estimate) |

### C.3 `fluidpy/core/shallow_water.py` (NEW, `SW`; re-exported by `ch13`) — numerical models (our choices)

| Function (signature) | Returns [units] | Source |
|---|---|---|
| `ShallowWater(nx, ny, Lx, Ly, H, f0, beta=0.0, g=G0, bc="periodic", bottom=None)` | grid object (Arakawa C layout as `core.mac.MacGrid`): attributes `x_c`, `y_c`, `dx`, `dy`, `H`, `f_u`, `f_v`, `f_q`; `bc` in ("periodic", "channel", "closed") | our choice of scheme (P335) |
| `gaussian_bump(model, amplitude, x0, y0, radius)` | array η `[ny, nx]` | helper |
| `geostrophic_state(model, eta)` | dict state {"eta", "u", "v"} in geostrophic balance with η | helper |
| `continuity_tendency(model, state, linear=True)` | array ∂η/∂t (flux form: sums to zero over a closed or periodic domain) | (13.44), (13.45a) |
| `momentum_tendencies(model, state, linear=True)` | tuple (∂u/∂t, ∂v/∂t) | (13.45b, c); (13.88)–(13.89) for `linear=False` |
| `relative_vorticity(model, state)` | array ζ at the corners | R25 |
| `potential_vorticity(model, state)` | array q = (ζ + f)/h at the corners (re-exported by `ch13` as `sw_potential_vorticity`, because `GFD.potential_vorticity` owns the bare name) | (13.94) |
| `energy(model, state)` | dict "kinetic", "potential", "total" [J per unit density] | diagnostic |
| `dt_limit(model, courant=0.5)` | largest stable time step [s] (gravity-wave CFL and inertial limit) | P335 |
| `step(model, state, dt, linear=True)` | new state | SSP-RK3 |
| `run(model, state, t_end, dt=None, linear=True, save_every=1)` | dict: "t" (n,), "eta" (n, ny, nx), "u", "v", "energy" (n,) | — |
| `advect_particles(model, history, xp0, yp0)` | arrays xp, yp (n, n_particles) by RK4 with bilinear interpolation | C13 figure (cached run) |
| `linear_1d_step(eta, u, v, dx, dt, H, f, g=G0)` **(+)** | tuple (eta, u, v) after one forward–backward step: η at centres (n,), u at faces (n + 1,) with u = 0 at both walls, v at centres | curation §8 note 7 |
| `linear_1d_run(eta0, dx, dt, n_steps, H, f, g=G0, save_every=1)` **(+)** | dict: "t", "eta" (n_saved, n), "u", "v", "pv" (linear potential vorticity v_x − fη/H at the interior faces, first and last saved step) | A2, E8 parity |
| `qg_linear_evolve(eta0, x, y, t, beta, f0, c)` | array η(x, y, t) by exact spectral evolution of (13.117) (`numpy.fft`; periodic) | (13.117) |
| `qg_linear_evolve_1d(amplitudes, k, l, t, beta, f0, c, U=0.0)` **(+)** | complex array a_j·exp(−iω_j t) with ω_j = `rossby_omega(k_j, l, …)` (the field is the real part of Σ a_j(t) e^{ik_jx}) | curation §8 note 8 |
| `barotropic_vorticity_rhs(zeta_hat, KX, KY, beta=0.0, nu=0.0, dealias=True)` | array: Fourier transform of −J(ψ, ζ) − β ∂ψ/∂x + ν∇²ζ | (13.122) with our viscosity |
| `barotropic_run(n, L, zeta0, *, t_end, dt, beta=0.0, nu=0.0, dealias=True, save_every=1)` | dict: "t", "zeta" (n_saved, n, n), "energy", "enstrophy", "K_E", "K_Z" (centroids) | N168 (ours); pseudo-spectral, RK4, 2/3 rule |
| `barotropic_invariants(zeta, L)` | dict "energy" [m²/s²], "enstrophy" [1/s²] (domain means) | (13.143), (13.144) |
| `spectral_centroids(zeta, L)` | dict "K_E", "K_Z" [rad/m]: energy- and enstrophy-weighted mean wavenumbers | ours |

### C.4 `fluidpy/ch13_geophysical_fluid_dynamics.py` (chapter module; re-exports every public name of C.1, C.2, C.3)

`from fluidpy.core.gfd import *`, `from fluidpy.core.vertical_modes import *`, `from fluidpy.core.shallow_water import *` (with
`__all__` in each core module), except that `SW.potential_vorticity` is re-exported as `sw_potential_vorticity`. Parity rows
write `ch13.<name>` only.

| Function (signature) | Returns | Source |
|---|---|---|
| `illustrative_inputs()` **(+)** | dict of our inputs (never the book's): "lat" 35° in rad, "lat_ekman" 60°, "lat_rossby" 12°, "ocean_H" 4200.0, "ocean_N" 2.7e-3, "atm_H" 9000.0, "atm_N" 1.1e-2, "atm_U0" 27.0, "tau" 0.07, "nu_v_ocean" 0.03, "U_g" 12.0, "nu_v_atm" 7.0, "U_syn" 14.0, "L_syn" 1.4e6, "H1" 120.0, "drho" 3.1, "rho_ocean" 1027.0, "U_mean" 17.0, "wavelengths" (7.3e5, 3.1e6), "q_inertial" 0.23, "u_rms" (13.0, 0.08) | curation decision 11 |
| `conventions_table()` **(+)** | `pandas.DataFrame`: the symbol table of convention 4 (letter, meanings, where) | conventions block |
| `book_slips()` **(+)** | dict keyed "1" … "13"; each value a dict with "id", "where", "printed", "correct", "taught_in", "corrected", "how_to_tell", "coded_in", "test" (LaTeX strings; no book numbers) | slip table |
| `traps()` **(+)** | list of 17 dicts, one per trap T1–T17 (kept apart from the slips) | trap table |
| `lapse_rate_table(dT_dz, *, Gamma_a)` **(+)** | `pandas.DataFrame` of two rows ("Kundu", "meteorology") with columns "symbol", "value_K_per_km", "adiabatic_K_per_km", "criterion", "verdict"; wraps `STRAT.lapse_rate_stability`; the criterion strings are ours ("dT/dz = −6.5 > −9.8 K/km", "Γ_met = 6.5 < Γ_d = 9.8 K/km") | lapse-rate rule |
| `wind_from_to(u, v)` | dict: "wind_name" (where it comes from, e.g. "westerly"), "current_name" (where it goes, e.g. "eastward"), "from_deg", "to_deg" (compass) | N02 |
| `ocean_density_gradient_budget(drho_dz, rho, c, g=G0)` | dict: "in_situ", "adiabatic" (= −gρ/c²), "potential" (their difference), "compression_share" | N03 with (13.1) |
| `atmosphere_layers(z)` | dict: "T" [K], "dT_dz" [K/m] (Kundu sign), "N2", "layer" (names) from USSA-1976 | N04 |
| `ocean_profile_idealized(z, *, T_surface=18.0, T_deep=2.0, h_mixed=50.0, h_thermocline=400.0, alpha=2.0e-4, rho0=1027.0)` | dict: "T" [°C], "rho_theta" [kg/m³], "N2" [1/s²] (analytic), "N" | N06 (ours) |
| `thermocline_N2(z, *, N_deep=5.0e-4, N_peak=8.0e-3, z_t=-300.0, width=150.0)` **(+)** | array N²(z) = (N_deep + (N_peak − N_deep)·exp(−((z − z_t)/width)²))² | E5, F6 (ours) |
| `anisotropic_eddy_stress(G, nu_H, nu_v, rho)` | array (3, 3) τ from G[i, j] = ∂u_i/∂x_j | (13.5) |
| `term_table_thin_layer(U, L, H, lat_rad, **kw)` **(+)** | `pandas.DataFrame` (rows = terms, columns = x and z components, size and ratio to the Coriolis term), built on `GFD.thin_layer_terms` | C01 figure |
| `pressure_centre(x, y, *, dp, R_c, lat_rad, rho0=1.2, kind="low")` | dict: "p" (Gaussian anomaly, −dp at the centre for a low) `[j, i]`, "u", "v" (geostrophic), "f" | R14, E1 |
| `jet_section(y, z, *, dT, width, lat_rad, alpha, T0=280.0, lapse=-6.5e-3, u_surface=0.0)` **(+)** | dict: "T" (z, y) [K] = T0 + lapse·z − (dT/2)·tanh(y/width)·sign(lat), "U" (z, y) [m/s] in exact thermal-wind balance, "dTdy" (y,) | curation §8 note 2 (ours) |
| `viscous_layer_thicknesses(nu, *, t=None, x=None, U=None, f=None)` | dict: "diffusive" √(νt), "boundary_layer" √(νx/U), "ekman" √(2ν/abs(f)) (None for the ones whose inputs are missing) | R16 |
| `ekman_surface_printed(z, tau, rho, nu_v, f)` | tuple (u, v): the book's cos / sin form, f > 0 only (raises for f ≤ 0) | C04 as printed |
| `coastal_upwelling(tau_alongshore, coast_side, lat_rad, rho=1027.0)` | dict: "transport_offshore" [m²/s] (positive = away from the coast), "upwelling" (bool), "kelvin_direction" ("poleward") ; `coast_side` in ("east", "west") = which side of the water the land is on; τ positive northward | N96 |
| `flow_over_step(x, U, beta, f0, h0, h1)` | dict: "Y" (streamline displacement [m]), "zeta" [1/s], "wavelength" (2π√(U/β) for U > 0, None for U < 0), "decay_length" (√(abs(U)/β)) | N102, N103, D18 (ours) |
| `vertical_structure_solve(z, N_fn, k, omega, f, w0=0.0, dw0=1.0)` | array ŵ(z) from `solve_ivp` (DOP853, rtol 1e-10) | (13.100) |
| `wkb_error(z, N_fn, k, omega, f)` | dict: "max_rel_error", "Hm" (smallest value of m·N/abs(dN/dz) on z); raises if the wave does not propagate everywhere on z. Measured: the error falls as 1/Hm (about 0.15/Hm), not as its square | N114 |
| `w_equation_rotating_residual(w_fn, x, y, z, t, N, f, h=1.0, ht=1.0)` | float residual of (13.96) by central differences | (13.96) |
| `lee_wave_field(x, z, U, N, k, h0)` | dict: "psi" (streamfunction of the total flow), "w", "m", "tilt" ("upstream") | N126 |
| `basin_crossing_time(L, lat_rad, c)` | time [s] = L/abs(`rossby_long_wave_speed`) | N135 |
| `rayleigh_kuo_eigs(k, U, Up, Upp, beta, domain=(-1.0, 1.0), N=120, delta=0.2, **kw)` | dict: "c" (complex leading eigenvalue; nan when no growing mode is found), "growth_rate" (k·Im c; 0.0 then), "stable" (bool) — passes β to `ST.rayleigh_eigs_contour(..., beta=beta)` as its own keyword (β is **not** folded into `Upp`); further keywords (`bc="decay"`, `y_max`, `parity`) go through. Non-dimensional: k in 1/L, β in U₀/L² | R30 |
| `eady_basic_state(y, z, *, N, f, H, U0, rho0, g=G0)` | dict: "rho" (z, y), "U" (z,), "slope" (isopycnal slope = fU₀/(N²H)), "drho_dy" | N143, (13.128) |
| `eady_matrix(c, alphaH, U0, H)` | complex array (2, 2): the coefficients of (A, B) in the two boundary conditions | N160 |
| `eady_numeric_eigs(k, l, N, f, H, U0, n=64)` | dict: "c" (complex, largest Im), "growth_rate" (Chebyshev generalised eigenproblem of (13.136) with the two lid conditions) | independent route |
| `load_reference_run(name)` **(+)** | dict of arrays from `reference/ch13/<name>.npz` ("kelvin_basin", "turbulence_f", "turbulence_beta", "pv_particles"), or `None` when the file is absent | caches (Part D) |
| sympy engines `boussinesq_rotating_sympy(printed=False)`, `perturbation_form_sympy()`, `eddy_friction_force_sympy()`, `thin_layer_sympy()`, `taylor_proudman_sympy(printed=False)`, `hydrostatic_linear_set_sympy()`, `v_equation_sympy()`, `rotating_internal_wave_set_sympy()`, `w_equation_rotating_sympy()`, `pv_conservation_sympy()`, `qg_vorticity_sympy()`, `eady_qg_sympy()` | each a dict with at least "residual" (sympy expression, 0 when the result holds) and "ok" (bool); `printed=True` builds the slip and returns "ok" False | slips #1, #3; D01, D02, D05, D12, D17, D19, D22, D25 |
| figure helpers `fig_ekman_surface(...)`, `fig_ekman_bottom(...)`, `fig_mode_roots(...)`, `fig_vertical_modes(...)`, `fig_poincare_kelvin_dispersion(...)`, `fig_kelvin_sections(...)`, `fig_inertia_gravity_orbit(...)`, `fig_rossby_dispersion(...)` | each returns a matplotlib `Figure` built from the functions above (all keyword arguments have defaults from `illustrative_inputs()`) | our analogues of the book's figures |

### C.4b `scripts/ch13_drawings.py` and `scripts/ch13_make_caches.py`
Sketch helpers imported by the notebook with `from ch13_drawings import *` (pure matplotlib, no physics; each returns a `Figure`):
`draw_tangent_plane()`, `draw_parcel_balance(hemisphere=+1)`, `draw_thermal_wind_wedge()`, `draw_taylor_column()`,
`draw_ekman_layer_sketch(kind="surface")`, `draw_shallow_layer()`, `draw_mode_stack()`, `draw_kelvin_sections()`,
`draw_pv_column()`, `draw_step_flow_sketch()`, `draw_eady_wedge()`, `draw_cascade_arrows()`. `scripts/ch13_make_caches.py
--which all` regenerates `reference/ch13/kelvin_basin.npz`, `turbulence_f.npz`, `turbulence_beta.npz`, `pv_particles.npz`
(float16 frames; < 2 MB together; our own model output only) and prints their sizes and a checksum.

### C.5 Not in `analysis/ch13.md` §4 (flag list for the implementer)
From curation §8: `parcel_adjust`, `ekman_transport_partial`, `shallow_water_branches`, `kelvin_omega`, `kelvin_decay_side`,
`adjustment_energy`, `rossby_packet_spectrum`, `eady_factors` (GFD); `vertical_modes_shooting`, `wkb_mode_speed` (VM; the
second is in §4 row I16); `linear_1d_step`, `linear_1d_run`, `qg_linear_evolve_1d` (SW); `jet_section`, `thermocline_N2`,
`term_table_thin_layer`, `lapse_rate_table`, `conventions_table`, `illustrative_inputs`, `book_slips`, `traps` (ch13). Added by this
design: `SIDEREAL_DAY`, `load_reference_run`, the alias `sw_potential_vorticity`, `scripts/ch13_make_caches.py`. Clarified
against §4 (not renamed): `thermal_wind_from_temperature` takes `alpha` keyword-only (curation §9); `barotropic_run` takes
`t_end`, `dt` keyword-only (§4 wrote them after defaulted arguments, which Python does not allow); `adjustment_energy` gains
`rho=1000.0` and a default `L=None`; `modes_uniform_N` gains `nz=201`; `kelvin_residuals` gains the difference steps `h`, `ht`;
`thin_layer_terms` has the explicit signature above; `pressure_centre`, `eady_basic_state`, `jet_section` take their physical
parameters keyword-only; `eady_fluxes` gains `k=None` (default: the fastest wave); `lee_wave_m` raises when evanescent.

### C.6 Parity conventions (the explainers reproduce these exactly)
1. `py:` expressions use only `ch13.<name>(…)`, `np.…`, `math.…`, numbers, strings, dict keys and integer indices. Latitudes
   are written `np.deg2rad(35.0)`; f is obtained in the row itself (`ch13.coriolis_parameter(np.deg2rad(35.0))`), never typed.
2. Every explainer with a hemisphere control has one parity row in each hemisphere.
3. Constants the page could compute are computed in the page and pinned by a row (`ch13.eady_critical()`,
   `ch13.eady_fastest()["sigma_nd"]`, `ch13.adjustment_energy(...)["ratio"]`); none is typed.
4. Complex results are compared through `np.real(...)` and `np.imag(...)` in two rows.
5. The JS mirrors use `OMEGA = 7.292115e-5` and `R = 6.371e6` (the two typed constants allowed: they are definitions) and are
   pinned by a row against `ch13.coriolis_parameter` and `ch13.beta_parameter`.

### C.7 Test-only and notebook-only functions
**Test-only (called by no notebook cell and no explainer):** `OMEGA_SOLAR_DAY` (private book-value tests; shown once as a
number in the trap box, computed inline as 2π/86400), `ekman_surface_printed` (the notebook shows the printed form as LaTeX
and one comparison cell — so it is called once; listed here because no explainer mirrors it), `SW.step`,
`SW.continuity_tendency`, `SW.momentum_tendencies`, `SW.energy`, `SW.dt_limit`, `SW.advect_particles`, `SW.geostrophic_state`, `SW.relative_vorticity`,
`SW.gaussian_bump`, `SW.run` and the `ShallowWater` class (used by `scripts/ch13_make_caches.py` and the slow tests; the
notebook loads their output through `load_reference_run`, and one optional `slow`-tagged cell under `if not FAST` calls
`SW.run` for 20 steps to show the API), `vertical_modes(method="cheb")` (the cross-check route), `eady_numeric_eigs` (called in
one check cell of C16; otherwise a test route), `w_equation_rotating_residual`, `rotating_internal_wave_set_sympy`,
`hydrostatic_linear_set_sympy`, `perturbation_form_sympy`, `thin_layer_sympy` (each called once in a one-line "the engine
agrees" cell under its derivation or note; no explainer). **Every other function of C.1–C.4 is called by a storyboard row
of Part A or mirrored by an explainer of Part B.**

## Part A — notebook storyboard (`notebooks/build_ch13.py` → `notebooks/ch13_geophysical_fluid_dynamics.ipynb`)

**One line per book section** (cell numbers are estimates; ≈ 700 cells in all):
- §13.1 → N01 N02 (open the C01 story; P307; the conventions block) — cells ≈ 10–20
- §13.2 → R01 R02 R03, N03 N04 N05 N06 N07 (stratification in two conventions; F2) — cells ≈ 21–45
- §13.3 → R04 R05 R06 R07 R08 R09 R10, N08 N09 N10 N11 (the rotating Boussinesq set and the eddy friction) — cells ≈ 46–66
- §13.4 → C01 (R11 R12 · N12–N17 · D01 D02 · P308 · F1 · plotly tangent plane · term bars) — cells ≈ 67–115
- §13.5 → C02 (R13 R14 · N18 N22 · D03 · P309 · **E1**), C03 (R15 · N20 N21 N23–N26 · D04 D05 · F3 · **E2**) — cells ≈ 116–190
- §13.6 → C04 (R16 · N27–N31 N33 N35 · D06 · P310 P311 · F4), C05 (N32 N36 N37 N50 · D07 D09 · P312 P313 · live · **E3**) — cells ≈ 191–265
- §13.7 → C06 (N34 N38–N49 · D08 · F5 · **E4**) — cells ≈ 266–305
- §13.8 → C07 (R17 R18 · N51–N54 · D10 · P335) — cells ≈ 306–335
- §13.9 → C08 (R19–R22 · N55–N77 · D11 · P314–P318 · F6 · **E5**) — cells ≈ 336–395
- §13.10 → C09 (N78–N83 · D12 D13 · P319 P320 · F7) — cells ≈ 396–430
- §13.11 → C10 (N84–N90 · D14 · P321 P322 · A1 · **E6**) — cells ≈ 431–465
- §13.12 → C11 (R23 · N91–N96 · D15 · P323 · A3 · **E7**), C12 (R24 · N19 · D16 · P324 · A2 · **E8**) — cells ≈ 466–520
- §13.13 → C13 (R25 R26 R27 · N97–N104 · D17 D18 · P325 · figs) — cells ≈ 521–560
- §13.14 → C14 (R28 R29 · N105–N126 · D19 D20 D21 · P326 P327 · plotly helix · figs) — cells ≈ 561–610
- §13.15 → C15 (N127–N138 · D22 D23 · P328 P329 · A4 · F8 · **E9**) — cells ≈ 611–650
- §13.16 → **C15 (continued)**: R30, N139–N142 · D24 — cells ≈ 651–665
- §13.17 → C16 (N143–N163 · D25 D26 D27 · P330 P331 P332 · A5 · F9 · live · **E10**) — cells ≈ 666–705
- §13.18 → C17 (R31 · N164–N168 · D28 D29 · P333 P334 · A6); S01, S02; summary — cells ≈ 706–735

Every CORE block follows: problem in plain words → idea → primers → maths (notes and derivations, Part F) → tiny example →
code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would change if…". Code
drafts give intent + exact calls; the builder comments every line (novice grade, units, the equation written out next to its
number). *expect* = numbers the executed cell must print (computed by the designer's scratch script from the stated inputs —
hypotheses to be reported against, not edited). *see / read / change* = the three figure notes. Note ids open every note in
bold (**N09 [C]**). Wherever a draft names an equation, the equation is written out with its number, and the builder keeps
it so. `inp = ch13.illustrative_inputs()` supplies every worked number; `lat = np.deg2rad(35.0)`, `f35 =
GFD.coriolis_parameter(lat)`, `f60`, `f12` likewise. **Every sentence that says "to the right", "clockwise" or "coast on the
right" carries "in the northern hemisphere (f > 0); mirror for f < 0".**

**Re-entry rule (binding; `tools/nbkit.py`).** `nb.recap(...)` and `nb.section(...)` close the current CORE block (they set
`nb.current_core = None`), and `nb.derivation`, the code count and the visual count all need an open block. This storyboard
places recaps *where they are needed*, often inside a block. The builder therefore writes `nb.current_core = "Cnn"` on the
line after every `nb.recap` that sits inside block Cnn (the pattern already used in `build_ch02.py`, `build_ch09.py` and
`build_ch12.py`), and after `nb.section("13.16", …)`, whose rows continue block C15. Rows marked **(re-enter Cnn)** below say
where. Recaps of §13.2 and §13.3 come before `nb.core("C01", …)` and need no re-entry.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…17…], prerequisites=[…])`. **Big idea (draft):** "Look at a weather map: the wind does not
   blow from high pressure to low, it blows *around* the highs and lows. Look at a map of ocean currents: the water piles up in
   the middle of each basin and circles it. Two things make the atmosphere and ocean behave so unlike water in a sink — the
   earth turns, and the fluid is layered, light over heavy. This chapter adds those two ingredients to everything the book has
   built. Rotation first: on a thin shell only the local vertical part of the earth's spin matters, and slow flows settle into
   a stand-off between the Coriolis force and the pressure gradient (geostrophy), with a vertical shear wherever temperature
   changes horizontally (thermal wind) and thin friction layers at the sea surface and the ground (Ekman layers) that move
   water at right angles to the wind. Then the waves that a rotating layer supports — fast ones that cannot go slower than the
   earth's own turning rate, a coast-hugging one, and one slow planetary wave that exists only because the Coriolis parameter
   changes with latitude. One conserved quantity, potential vorticity, organises all the slow motion. Finally, why weather
   exists at all (baroclinic instability) and why its eddies grow into jets instead of breaking down (two-dimensional
   turbulence)." **Roadmap (one line per CORE):** C01 the equations on a thin rotating shell · C02 geostrophic balance · C03
   thermal wind · C04 the surface Ekman spiral · C05 Ekman transport and pumping · C06 the bottom Ekman layer · C07 the
   shallow-water equations · C08 vertical normal modes · C09 one cubic for all the waves · C10 Poincaré waves · C11 the Kelvin
   wave · C12 the Rossby radius and geostrophic adjustment · C13 potential vorticity · C14 inertia–gravity waves · C15 Rossby
   waves (and barotropic instability) · C16 the Eady problem · C17 two conserved quantities and the inverse cascade.
   **Prerequisites:** the rotating-frame momentum equation $\rho\big(\frac{D'\mathbf u'}{Dt}\big)_{O'1'2'3'}=-\nabla'p+\rho\big[\mathbf g-\frac{d\mathbf U}{dt}-2\boldsymbol\Omega\times\mathbf u'-\frac{d\boldsymbol\Omega}{dt}\times\mathbf x'-\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf x')\big]+\mu\nabla'^2\mathbf u'$ (4.45) (for the steadily turning earth only the Coriolis and centrifugal terms of the bracket survive) and the Boussinesq
   approximation (Ch. 4), static stability, potential temperature and the two lapse-rate conventions (Ch. 1), vorticity,
   stretching and the column argument (Ch. 5), surface and internal gravity waves, group velocity (Ch. 7), the Stokes layer
   (Ch. 8), boundary-layer scaling (Ch. 9), Rayleigh's criterion and normal modes with a complex phase speed (Ch. 11), eddy
   viscosity and the energy cascade (Ch. 12).
2. `nb.explainer_index([...])` — 10 rows: ("geostrophic_balance", "Why doesn't air flow straight from high to low pressure?",
   "the Coriolis force turns it until the only motion left runs along the isobars") · ("thermal_wind", "Why is there a jet
   stream, and why is it westerly in both hemispheres?", "a horizontal temperature gradient sets the change of wind with
   height") · ("ekman_spiral", "The wind blows east — why does the water go south?", "the eddy viscosity sets how deep the
   spiral reaches, not how much water it carries") · ("ekman_force_balance", "Why does air spiral into a low?", "friction
   weakens the Coriolis force near the ground and the pressure force wins") · ("vertical_modes", "How can one layer stand for a
   stratified ocean?", "each vertical mode is a shallow-water system of its own, much shallower, depth") ·
   ("shallow_water_dispersion", "Are Poincaré, Kelvin and Rossby waves separate theories?", "three roots of one cubic, decades
   apart in frequency") · ("kelvin_wave", "Why does this wave run only one way along a coast?", "only one direction gives a
   slope that decays offshore") · ("geostrophic_adjustment", "Why doesn't a pile of water flatten out on a rotating planet?",
   "rotation holds up whatever lies beyond a Rossby radius") · ("rossby_waves", "The crests go west — so how does the energy go
   east?", "phase and group velocity part company at the Rossby radius") · ("eady_instability", "Where do storms come from?",
   "two boundary waves that hold each other in place and grow, feeding on the tilt of the density surfaces").
3. `nb.setup()`.
4. `nb.code` — **chapter imports**: `import numpy as np`, `import sympy as sp`, `import matplotlib.pyplot as plt`, `import pandas
   as pd`, `import plotly.graph_objects as go`, `from scipy.integrate import solve_ivp, cumulative_trapezoid`, `from
   scipy.optimize import brentq, minimize_scalar`, `from scipy.linalg import eigh_tridiagonal`, `from fluidpy import
   ch13_geophysical_fluid_dynamics as ch13`, `from fluidpy.core import gfd as GFD, vertical_modes as VM, shallow_water as SW`,
   `from fluidpy.core import rotating as ROT, stratification as STRAT, similarity as SIM, waves as WAV, stability as ST, laminar
   as LAM, turbstats as TS, dimensional as DIM`, `from fluidpy import ch05_vorticity_dynamics as ch05, ch11_instability as ch11,
   ch12_turbulence as ch12`, `from fluidpy.core.interact import slider_figure, animate_figure, live`, `from fluidpy.core.anim
   import animate`, `from fluidpy.core.style import COLORS, savefig`, `from fluidpy.core.thermo import G0`, `import sys;
   sys.path.insert(0, "scripts"); from ch13_drawings import *`; then `inp = ch13.illustrative_inputs()`, `lat =
   inp["lat"]`, `f35 = GFD.coriolis_parameter(lat)`, `f60 = GFD.coriolis_parameter(inp["lat_ekman"])`, `f12 =
   GFD.coriolis_parameter(inp["lat_rossby"])`, `beta35 = GFD.beta_parameter(lat)`, `beta12 =
   GFD.beta_parameter(inp["lat_rossby"])`. *explain:* one line per import ("`GFD`, `VM`
   and `SW` are the three new core modules; `ch13` re-exports all of them, so `ch13.rossby_omega` and `GFD.rossby_omega` are the
   same function"). *expect:* f35 = 8.365 × 10⁻⁵ s⁻¹, f60 = 1.263 × 10⁻⁴ s⁻¹, f12 = 3.032 × 10⁻⁵ s⁻¹, beta35 = 1.875 × 10⁻¹¹ m⁻¹ s⁻¹.
5. `nb.md` — **⚠️ Conventions in this chapter** (the conventions block of curation decision 8, **before C01**): (a) axes: x
   east, y north, z up, (u, v, w); (b) the table "where is z = 0?" (sea surface in §13.6, §13.9, §13.14 · solid surface in §13.7
   · flat bottom in §13.8 · lower lid in §13.17) — trap **T4**; (c) hemisphere: "every 'right', 'clockwise' and 'coast on the
   right' in the book is the northern-hemisphere case f > 0; mirror for f < 0" — trap **T5**; (d) primes dropped from the
   thin-shell equations on — trap **T3**; (e) `ch13.conventions_table()` displayed (the symbol table of header convention 4:
   f, Ω / ω, β, N, H / h / η, c / c_n, k / l / m / K, ζ, θ, δ, Λ / λ, ψ, E / Ro / R, α, U / V / τ) — trap **T13**; (f) three
   Rossby radii (trap **T11**), two Rossby numbers (trap **T17**), the stream-function sign (trap **T14**), the spectrum
   normalisation (trap **T16**), each in one line with the block that uses it; (g) > ⚠️ **slip #9 — the book prints**
   cross-references to "Section 4.18", "Section 4.11", "Section 8.7", "Section 5.7" and to equation number 7.128**; the
   correct form is** §4.9, §4.9, §8.4, §5.6 and the definition $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ (7.127) — we cite our section
   numbers with the equation written out; (h)
   `pd.DataFrame(ch13.book_slips()).T` (13 rows; columns where, printed, correct, taught_in, corrected, how_to_tell, coded_in,
   test) and `pd.DataFrame(ch13.traps())` (17 rows) with the sentence "slips are false as printed; traps are true but easy to
   misread — they are kept apart". (i) "**Which inputs?** All worked numbers use our own illustrative inputs
   (`ch13.illustrative_inputs()`), never the book's."

### A.1 §13.1 Introduction — N01 N02 (open the C01 story)
1. `nb.section("13.1", "Introduction", intro="**What is this section about?** What makes the atmosphere and the ocean a subject
   of their own: the earth's rotation and the layering of density. And the words — easterly, eastward, cyclonic — that the
   rest of the chapter uses.")`
2. `nb.note` — **N01 [C]** "**The picture.** Geophysical fluid dynamics is the dynamics of the atmosphere and the ocean. Two
   features set it apart: rotation and vertical density stratification. Their most visible signature is that large-scale
   flow runs *along* lines of constant pressure, not across them." + `nb.figure`: two panels from `ch13.pressure_centre(x, y,
   dp=400.0, R_c=6e5, lat_rad=lat, kind="low")` (isobars + geostrophic quiver) and the same with `kind="high"`. *see:* arrows
   circling the centres; *read:* "no arrow points at the low — C02 explains why"; *change:* "…at 35° S both circulations
   reverse (C02)". Pointer: rotation → C01–C06; stratification → §13.2 recap, C08, C14, C16.
3. `nb.primer("geographic vocabulary: zonal / meridional, easterly wind vs eastward current, poleward / equatorward, cyclonic
   / anticyclonic", "Zonal = along a latitude circle (east–west); meridional = north–south. A wind is named by where it comes
   FROM (a westerly blows toward the east); a current by where it goes TO (an eastward current). Poleward and equatorward
   avoid saying north or south. Cyclonic = turning in the same sense as the earth below: counter-clockwise in the northern
   hemisphere, clockwise in the southern.", code="print(ch13.wind_from_to(10.0, 0.0))     # u = +10 m/s: a 'westerly' wind, an 'eastward'
   current\nprint(GFD.hemisphere(np.deg2rad(-35.0)))  # sign -1, turns 'left', cyclonic 'clockwise'")` (**P307**)
4. `nb.note` — **N02 [B]** the conventions box in words (z up; winds by origin, currents by destination; heights and depths;
   "all the book's figures are for the northern hemisphere; the sign of every rotational effect follows sign(f)") +
   `GFD.hemisphere(lat)` and `GFD.hemisphere(-lat)` printed side by side as a two-row table (hemisphere · sign of f · deflection
   · sense round a low).

### A.2 §13.2 Vertical Variation of Density in the Atmosphere and Ocean — R01 R02 R03, N03–N07
1. `nb.section("13.2", "Vertical Variation of Density in the Atmosphere and Ocean", intro="**What is this section about?**
   A reminder of static stability from Chapter 1, in the form this chapter needs: which density gradient decides stability,
   how the lapse rate says the same thing in two sign conventions, and the buoyancy frequency N(z) that later sets the speed
   of internal modes (C08), the frequency range of internal waves (C14) and the growth of storms (C16).")`
2. `nb.recap("R01", "Static stability is decided by the potential density", "A parcel moved up expands and its density falls
   even if nothing else changes. What decides whether it sinks back is the density it would have at a common pressure — the
   potential density $\\rho_\\theta$: $\\dfrac{d\\rho_\\theta}{dz}=\\dfrac{d\\rho}{dz}+\\dfrac{g\\rho}{c^2}$ (13.1), where $c$ is the speed of sound ($c_s$ in our own
   lines). Stable when $d\\rho_\\theta/dz<0$.", where="Ch. 1 §1.10, where the same statement is $\\frac{d\\rho_\\theta}{dz}=\\frac{d\\rho}{dz}-\\frac{d\\rho_a}{dz}\\cong\\frac{d\\rho}{dz}+\\frac{\\rho g}{c^2}$ (1.35)")` + `nb.code`:
   `STRAT.ocean_potential_density_gradient(drho_dz, rho, c_sound)`, `STRAT.isentropic_density_gradient(rho, c_sound)`.
3. `nb.note` — **N03 [B]** "Most of the increase of in-situ density with depth is compression." `b =
   ch13.ocean_density_gradient_budget(drho_dz=-5.24e-3, rho=1027.0, c=1500.0)`. *expect* (inputs ρ = 1027 kg m⁻³, c_s = 1500
   m s⁻¹, in-situ gradient −5.240 × 10⁻³ kg m⁻⁴): adiabatic −4.476 × 10⁻³, potential −7.63 × 10⁻⁴ kg m⁻⁴, compression share
   85.4 %. "So N² must be built from the potential density — the in-situ gradient would overstate it almost sevenfold."
   (Builder: write the unit as kg m⁻⁴, never with a slash.)
4. `nb.recap("R02", "The adiabatic gradient, in two sign conventions", "The adiabatic density gradient is $-g\\rho/c^2$. For a
   gas the same idea is the adiabatic temperature gradient. **Kundu:** $\\Gamma\\equiv dT/dz$, adiabatic value $\\Gamma_a=-g/C_p\\approx-9.8$
   K/km, stable when $dT/dz>\\Gamma_a$. **Meteorology:** $\\Gamma_{met}\\equiv-dT/dz$, dry adiabatic value $\\Gamma_d=+g/C_p\\approx+9.8$ K/km, stable when
   $\\Gamma_{met}<\\Gamma_d$. Same physics: negate the number, flip the inequality. This section of the book writes no symbol Γ; it
   quotes rates of decrease in words, which are the meteorological magnitudes.", where="Ch. 1 §1.10 (C52–C55) and Ch. 12
   (C14)")` + `nb.code`: `Gamma_a = STRAT.adiabatic_lapse_rate()`; `ch13.lapse_rate_table(-6.5e-3, Gamma_a=Gamma_a)` displayed.
   *expect:* Γ_a = −9.77 K/km; row "Kundu": dT/dz = −6.5 > −9.77 K/km → stable; row "meteorology": Γ_met = 6.5 < Γ_d = 9.77 K/km
   → stable. + `nb.md` **⚠️ Common confusion (trap T2, the library string):** "`STRAT.lapse_rate_stability(...,
   convention='meteorology').text` prints 'Γ < Γa'. Read it as $\Gamma_{met}<\Gamma_d$: in that convention the library's Γa is the
   positive number +9.77 K/km. We do not show that raw string again."
5. `nb.worked_example("is a −8 K/km layer stable?", "1. Kundu: dT/dz = −8 K/km; is −8 > −9.8? Yes → stable. 2. Meteorology:
   negate: Γ_met = +8 K/km; flip the inequality: is 8 < 9.8? Yes → stable. 3. Margin: 1.8 K/km in both. 4. A −11 K/km layer:
   −11 > −9.8 is false; 11 < 9.8 is false → unstable in both.")`
6. `nb.note` — **N04 [B]** our Fig. 13.1 (`ch13.atmosphere_layers(z)`, z = 0…50 km): left T(z) with troposphere, tropopause,
   stratosphere, stratopause named; right N²(z). Legend carries both forms: "troposphere: dT/dz = −6.5 > −9.8 K/km ⇔ Γ_met =
   6.5 < 9.8 K/km". *expect:* N² = 1.11 × 10⁻⁴ s⁻² at the ground (N = 1.05 × 10⁻² s⁻¹, period 9.9 min); about four times larger
   in the isothermal lower stratosphere. *see / read* ("N² jumps at the tropopause: the stratosphere is a lid") */ change*
   ("…the troposphere followed the dry adiabat: N² = 0 there and the Eady growth rate of C16 would be unbounded").
7. `nb.note` — **N05 [B]** **⚠️ Common confusion (trap T2):** "'The troposphere is close to neutral' compares the observed
   lapse rate with the **moist** adiabat (named only: saturated air cools more slowly because condensation releases heat).
   To dry displacements the standard troposphere is statically stable in both conventions: $\Gamma_a<dT/dz<0$ (Kundu), i.e.
   $0<\Gamma_{met}<\Gamma_d$ (meteorology); N² > 0." + `nb.code`: `ch12.gradient_richardson_thermal(-6.5e-3, 3.0e-3, 1/288.15,
   Gamma_a=Gamma_a)["Ri"]` *expect:* about +12 for a shear of 3 m/s per km (N²/shear² = 1.11 × 10⁻⁴/9 × 10⁻⁶). Badge text built
   by us: "stable · dT/dz = −6.5 > −9.8 K/km · Γ_met = 6.5 < 9.8 K/km".
8. `nb.plotly` — **F2** (`slider_figure`, slider = environmental dT/dz from −12 to +2 K/km, 29 steps): trace 1 the buoyancy
   acceleration of a parcel lifted 100 m, trace 2 N² (negative values dashed), trace names carry both conventions and the
   verdict from `ch13.lapse_rate_table` ("dT/dz = −6.5 > −9.8 K/km ⇔ Γ_met = 6.5 < 9.8 K/km → stable"). *see / read* ("the
   verdict flips where the two lines cross zero — at the adiabatic value, whichever sign you write it with") */ change*.
9. `nb.note` — **N06 [B]** our Fig. 13.2 (`prof = ch13.ocean_profile_idealized(z)`, z from −4200 m to 0): T, potential density
   and N against depth; mixed layer, thermocline, abyss named. "Salinity is secondary and left out of our profile." The
   same profile feeds C08.
10. `nb.recap("R03", "The buoyancy frequency", "$N^2\\equiv-\\dfrac{g}{\\rho_0}\\dfrac{d\\rho}{dz}$ with ρ the **potential** density and ρ₀ a constant
    reference value: the squared frequency at which a displaced parcel bobs. Where it enters this chapter: C08 (mode speeds
    $c_n\\approx NH/n\\pi$), C14 (internal waves need ω < N), C16 (the Eady radius NH/f and growth rate).", where="Ch. 1, $N^2=-\\frac g\\rho\\big(\\frac{d\\rho}{dz}-\\frac{d\\rho_a}{dz}\\big)$
    (1.29); Ch. 7, $N^2\\equiv-\\frac g{\\rho_0}\\frac{d\\bar\\rho}{dz}$ (7.127)")` + `nb.code`: `GFD.buoyancy_frequency_sq(z, prof["rho_theta"], 1027.0)` against
    `prof["N2"]` (`assert np.allclose(..., rtol=2e-2)` away from the ends).
11. `nb.note` — **N07 [C]** "The largest N of our thermocline profile and its period: `STRAT.stability_timescale(N2.max())`."
    *expect* (N_peak-like value 8 × 10⁻³ s⁻¹): about 13 min. "This is the upper frequency limit of C14 and the profile maximum
    of C08."

### A.3 §13.3 Equations of Motion — R04–R10, N08–N11
1. `nb.section("13.3", "Equations of Motion", intro="**What is this section about?** The starting equations — Chapter 4's
   Boussinesq set written in a rotating frame — and the one new modelling choice: friction by eddies, with a large
   horizontal and a small vertical eddy viscosity.")`
2. `nb.recap("R04", "The rotating Boussinesq equations", "Continuity $\\nabla\\cdot\\mathbf u=0$; momentum
   $\\dfrac{D\\mathbf u}{Dt}+2\\boldsymbol\\Omega\\times\\mathbf u=-\\dfrac1{\\rho_0}\\nabla p-\\dfrac{g\\rho}{\\rho_0}\\mathbf e_z+\\mathbf F$; density $\\dfrac{D\\rho}{Dt}=0$ — together (13.2).
   Term by term: acceleration seen from the turning earth, the Coriolis acceleration, the pressure-gradient force, weight,
   friction per unit mass.", where="Ch. 4: the rotating-frame equation $\\rho\\big(\\frac{D'\\mathbf u'}{Dt}\\big)_{O'1'2'3'}=-\\nabla'p+\\rho\\big[\\mathbf g-\\frac{d\\mathbf U}{dt}-2\\boldsymbol\\Omega\\times\\mathbf u'-\\frac{d\\boldsymbol\\Omega}{dt}\\times\\mathbf x'-\\boldsymbol\\Omega\\times(\\boldsymbol\\Omega\\times\\mathbf x')\\big]+\\mu\\nabla'^2\\mathbf u'$ (4.45) and the
   Boussinesq momentum equation $\\frac{D\\mathbf u}{Dt}=-\\frac1{\\rho_0}\\nabla p'+\\frac{\\rho'}{\\rho_0}\\mathbf g+\\nu\\nabla^2\\mathbf u$ (4.86)")` + > ⚠️ **slip #1 — the book prints** the pressure
   term of the momentum equation as $+\frac1{\rho_0}\nabla p$**; the correct form is** $-\frac1{\rho_0}\nabla p$ (a fluid at rest must satisfy
   $0=-\nabla p-g\rho\,\mathbf e_z$). + `nb.code`: `ch13.boussinesq_rotating_sympy()["ok"]` → True;
   `ch13.boussinesq_rotating_sympy(printed=True)["residual"]` → non-zero ($2\nabla\bar p/\rho_0$).
3. `nb.recap("R05", "When Boussinesq holds", "The layer must be thin against the scale height $c^2/g$ (written $c_s^2/g$ in
   our lines), over which density changes by its own size through compression.", where="Ch. 4 §4.9")` + `GFD.scale_height(1500.0)`,
   `GFD.scale_height(340.0)` *expect:* 229 km (ocean) and 11.8 km (air): "true for any ocean depth; only roughly true for the
   troposphere — which is why meteorology often uses pressure coordinates (not in the book)".
4. `nb.recap("R06", "Density is carried with the fluid", "$D\\rho/Dt=0$ follows from $DT/Dt=0$ (or $DS/Dt=0$) and a linear
   equation of state, $\\delta\\rho/\\rho_0=-\\alpha\\,\\delta T$ and $\\delta\\rho/\\rho_0=\\beta_S\\,\\delta S$. ⚠️ The book writes the haline coefficient as β; from
   §13.4 on β is df/dy, so we write $\\beta_S$ here. The expansion coefficient α comes back in C03.", where="Ch. 4, the Boussinesq heat equation $DT/Dt=\\kappa\\nabla^2T$ (4.89) with its diffusion term
   dropped; Ch. 1 linear equation of state")`
5. `nb.recap("R07", "The state of rest", "$\\dfrac{d\\bar p}{dz}=-\\bar\\rho g$ (13.3): the hydrostatic balance of the motionless reference
   state.", where="Ch. 1, $dp/dz=-\\rho g$ (1.8)")`
6. `nb.recap("R08", "Rest state plus perturbation", "$\\rho(\\mathbf x,t)=\\bar\\rho(z)+\\rho'(\\mathbf x,t)$ and $p(\\mathbf x,t)=\\bar p(z)+p'(\\mathbf x,t)$ (13.4).
   ⚠️ **Trap T3:** from the thin-shell equations of C01 on, the book drops the primes: p and ρ there are *perturbations*.",
   where="Ch. 4, where $p'=p-p_s$, $\\rho'=\\rho-\\rho_s$ give $\\rho\\frac{D\\mathbf u}{Dt}=-\\nabla p'+\\rho'\\mathbf g+\\mu\\nabla^2\\mathbf u$ (4.84); Ch. 7, $\\rho=\\bar\\rho(z)+\\rho'$ (7.124)")`
7. `nb.recap("R09", "The pressure and gravity terms keep their form", "Subtract the rest state: $-\\dfrac1{\\rho_0}\\nabla p-\\dfrac{g\\rho}{\\rho_0}\\mathbf e_z=
   -\\dfrac1{\\rho_0}\\nabla p'-\\dfrac{g\\rho'}{\\rho_0}\\mathbf e_z$, because $\\nabla\\bar p=(d\\bar p/dz)\\mathbf e_z=-\\bar\\rho g\\,\\mathbf e_z$ cancels the weight of the rest state.",
   where="Ch. 4, derivation D27")` + `ch13.perturbation_form_sympy()["ok"]` → True.
8. `nb.recap("R10", "Friction per unit mass, and the eddy-viscosity idea", "The friction force per unit mass is the divergence
   of the stress divided by density, $F_i=\\dfrac1\\rho\\dfrac{\\partial\\tau_{ij}}{\\partial x_j}$. For large-scale flow the stress is carried by turbulent eddies
   and is modelled as an eddy viscosity times a velocity gradient.", where="Ch. 4 §4.4 and Ch. 12, the eddy-viscosity hypothesis $\\overline{u_iu_j}=\\frac23\\bar e\\delta_{ij}-\\nu_T\\big(\\frac{\\partial U_i}{\\partial x_j}+\\frac{\\partial U_j}{\\partial x_i}\\big)$ (12.94)")` +
   > ⚠️ **slip #2 — the book prints** $F_i=\partial\tau_{ij}/\partial x_j$ for a force per unit mass**; the correct form is**
   $F_i=\frac1\rho\,\partial\tau_{ij}/\partial x_j$ (N/m³ divided by kg/m³ gives m/s²).
9. `nb.note` — **N08 [B]** the six stresses, all of (13.5): $\tau_{xz}=\tau_{zx}=\rho\nu_v\frac{\partial u}{\partial z}+\rho\nu_H\frac{\partial w}{\partial x}$, $\tau_{yz}=\tau_{zy}=\rho\nu_v\frac{\partial v}{\partial z}+\rho\nu_H\frac{\partial w}{\partial y}$,
   $\tau_{xy}=\tau_{yx}=\rho\nu_H\big(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\big)$, $\tau_{xx}=2\rho\nu_H\frac{\partial u}{\partial x}$, $\tau_{yy}=2\rho\nu_H\frac{\partial v}{\partial y}$, $\tau_{zz}=2\rho\nu_v\frac{\partial w}{\partial z}$. "Two coefficients because the
   eddies are wide and flat: horizontal exchange is strong, vertical exchange weak, $\nu_H\gg\nu_v$." `tau =
   ch13.anisotropic_eddy_stress(G, nu_H, nu_v, rho)`; check: with `nu_H == nu_v` the result equals $2\rho\nu S_{ij}$ of Ch. 4
   (`assert np.allclose`).
10. `nb.note` — **N09 [C]** "This stress is not frame-indifferent: a rigid rotation in a vertical plane, $u=\omega_rz$,
    $w=-\omega_rx$, has no deformation yet gives $\tau_{xz}=\rho\omega_r(\nu_v-\nu_H)\neq0$ (one line of code). The book accepts the defect for
    simplicity; so do we. Pointer: Ch. 4 §4.5."
11. `nb.note` — **N10 [B]** the friction force, all of (13.6) per unit mass: $F_x=\nu_H\big(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\big)+\nu_v\frac{\partial^2u}{\partial z^2}$,
    $F_y=\nu_H\big(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\big)+\nu_v\frac{\partial^2v}{\partial z^2}$, $F_z=\nu_H\big(\frac{\partial^2w}{\partial x^2}+\frac{\partial^2w}{\partial y^2}\big)+\nu_v\frac{\partial^2w}{\partial z^2}$. "The book says the
    stresses 'become' this. The six differentiations and the use of continuity are written out as derivation D01 inside the
    C01 block below." `GFD.eddy_friction(lap_h, d2z, nu_H, nu_v)`.
12. `nb.note` — **N11 [C]** a small table of **our own illustrative** eddy coefficients against molecular values (labelled
    illustrative): lower atmosphere ν_v ≈ 7 m² s⁻¹, ν_H ≈ 3 × 10⁴ m² s⁻¹; upper ocean ν_v ≈ 0.03 m² s⁻¹, ν_H ≈ 300 m² s⁻¹;
    molecular: air 1.5 × 10⁻⁵, water 1.0 × 10⁻⁶ m² s⁻¹. "Used where ν_v sets the Ekman thickness (C04, C06)."

---

### A.4 §13.4 Approximate Equations for a Thin Layer on a Rotating Sphere — C01
#### C01 — The thin-shell equations (13.9), with the Coriolis parameter (13.8)
1. `nb.section("13.4", "Approximate Equations for a Thin Layer on a Rotating Sphere", intro="**What is this section about?**
   Standing at one latitude, we replace the sphere by its tangent plane and ask which parts of the earth's rotation a thin
   layer of fluid can feel. The answer is one number per latitude, f, and a set of equations that every later section
   simplifies further.")`
2. `nb.core("C01", "The equations of motion on a thin rotating shell: $\\dfrac{Du}{Dt}-fv=-\\dfrac1{\\rho_0}\\dfrac{\\partial p}{\\partial x}+F_x$,
   $\\dfrac{Dv}{Dt}+fu=-\\dfrac1{\\rho_0}\\dfrac{\\partial p}{\\partial y}+F_y$, $\\dfrac{Dw}{Dt}=-\\dfrac1{\\rho_0}\\dfrac{\\partial p}{\\partial z}-\\dfrac{g\\rho}{\\rho_0}+F_z$ (13.9), with $f=2\\Omega\\sin\\theta$ (13.8)",
   question="Which part of the earth's rotation does a thin layer of air or water actually feel?")`
3. `nb.md` — **Plain words:** "The atmosphere is about 10 km deep and its weather systems are 1000 km and more across: a sheet
   thinner, in proportion, than the paper of a wall map. In such a sheet the fluid can hardly move up or down. Of the earth's
   spin, the part that matters is the part about the *local vertical* — a turntable under your feet that turns once a day at
   the pole and not at all on the equator. Climate models call dropping the rest the 'traditional approximation'; every
   dynamical core makes it."
4. `nb.md` — **The idea** (ASCII): the rotation vector seen from a point at latitude θ splits into a vertical part Ω sin θ
   ("spins the ground under you") and a northward horizontal part Ω cos θ ("tips the ground; only couples to vertical motion").
   Two-column table: pole (θ = 90°: all vertical, f = 2Ω) · equator (θ = 0: all horizontal, f = 0).
5. **(re-enter C01 after this recap)** `nb.recap("R11", "How fast the earth turns", "One turn against the stars takes 86 164 s (a sidereal day), so
   $\\Omega=2\\pi/86\\,164\\ \\mathrm{s}=7.292\\times10^{-5}$ rad/s (`ROT.OMEGA_EARTH`). ⚠️ **Trap T1:** the book takes $\\Omega=2\\pi$ rad/day with a
   24-hour solar day, 0.27 % smaller. Every number of ours uses the sidereal value.", where="Ch. 4 §4.7")` +
   `nb.primer("sidereal day vs solar day", "In one year the earth turns 366.25 times against the stars but the sun crosses
   the sky only 365.25 times, because one turn is used up going round the sun. So one true rotation takes 24 h × 365.25/366.25
   = 23 h 56 min 4 s. Dynamics cares about the true rotation.", code="T_sid = 86400.0 * 365.25 / 366.25       # one
   rotation against the stars [s]\nprint(T_sid, 2*np.pi/T_sid)             # 86164.1 s, 7.2921e-05 rad/s\nprint(GFD.OMEGA_EARTH /
   (2*np.pi/86400.0) - 1)   # 0.00274: the solar-day value is 0.27 % smaller")` (**P308**) *expect:* 86164.1; 7.2921e-05; 0.00274.
6. `nb.note` — **N13 [B]** $\boldsymbol\Omega=(0,\ \Omega\cos\theta,\ \Omega\sin\theta)$ on the tangent plane + `nb.plotly`: a 3-D sphere with the
   rotation axis, a tangent plane at a dropdown of latitudes (0°, 35°, 60°, 90°) and the two components of Ω as arrows
   (`GFD.earth_rotation_local(lat)`); `height=520`. *see / read* ("the teal vertical arrow is f/2") */ change*.
7. `nb.note` — **N12 [B]** $W/U\sim H/L$ from continuity: `GFD.vertical_velocity_scale(14.0, 9000.0, 1.4e6)` *expect:* W = 0.09
   m/s, aspect ratio 6.4 × 10⁻³. "Step 3 of D02."
8. `nb.derivation("D01", …)` — Part F D01 (6 steps), ref "13.6": the friction force from the anisotropic stress (N10; the book
   never writes it out). + `nb.code`: `ch13.eddy_friction_force_sympy()["ok"]` → True.
9. `nb.derivation("D02", …)` — Part F D02 (9 steps), ref "13.9": **N14 [B]** is its steps 4–6: before the approximation the
   x-component of the Coriolis acceleration is $2\Omega(w\cos\theta-v\sin\theta)$; dropping $w\cos\theta$ against $v\sin\theta$ and the vertical
   component against gravity gives $2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)$ (13.7) — the "traditional approximation" (the
   book does not name it). + `ch13.thin_layer_sympy()["ok"]` → True.
10. **(re-enter C01)** `nb.recap("R12", "The Coriolis parameter", "$f=2\\Omega\\sin\\theta$ (13.8): twice the local vertical component of the earth's
    rotation, also called the planetary vorticity. Positive in the northern hemisphere, negative in the southern, zero on the
    equator, odd in latitude.", where="Ch. 4 §4.7, Ch. 5 §5.6")` + `nb.code`: `GFD.coriolis_parameter(lat)`,
    `GFD.coriolis_parameter(-lat)`, at 90°. *expect:* +8.365 × 10⁻⁵, −8.365 × 10⁻⁵, 1.458 × 10⁻⁴ s⁻¹.
11. `nb.note` — **N15 [B]** inertial period $T_i=2\pi/f$: `GFD.inertial_period(f35)/3600` *expect:* 20.86 h (13.82 h at 60°, 57.56 h
    at 12°, 11.97 h at the pole). > ⚠️ **slip #5 — the book prints** that the subscript of $T_i$ "does refer to a
    component"**; the correct form is** "does not": i stands for "inertial". Used in C10.
12. `nb.worked_example("f, β and the size of what we dropped at 30° N", "Take Ω ≈ 7.3 × 10⁻⁵ s⁻¹, R ≈ 6400 km, sin 30° = 0.5,
    cos 30° ≈ 0.87. 1. $f=2\\Omega\\sin\\theta$ = 2 × 7.3 × 10⁻⁵ × 0.5 = 7.3 × 10⁻⁵ s⁻¹. 2. Inertial period 2π/f = 86 000 s ≈ 24 h (at 30°
    the pendulum day is one day). 3. $\\beta=2\\Omega\\cos\\theta/R$ = 1.46 × 10⁻⁴ × 0.87/6.4 × 10⁶ ≈ 2.0 × 10⁻¹¹ m⁻¹ s⁻¹. 4. A storm with U = 20
    m/s, L = 2000 km, H = 10 km: W ≈ UH/L = 0.1 m/s. 5. Dropped term over kept term: (W cos θ)/(U sin θ) = (0.1 × 0.87)/(20 ×
    0.5) ≈ 0.009 — under one per cent.")`
13. `nb.code` — `terms = GFD.thin_layer_terms(14.0, 1.4e6, 9000.0, lat, nu_H=3e4, nu_v=7.0)`; `ch13.term_table_thin_layer(14.0,
    1.4e6, 9000.0, lat, nu_H=3e4, nu_v=7.0)`; `GFD.coriolis_acceleration_local(14.0, 0.0, 0.09, lat, thin=False)` against
    `thin=True`. *expect:* Ro = 0.120; x-equation: acceleration 1.4 × 10⁻⁴, Coriolis fU = 1.17 × 10⁻³, dropped 2Ω cos θ·W = 1.08 ×
    10⁻⁵ m/s² (0.9 % of fU); z-equation: Coriolis 2ΩU cos θ = 1.67 × 10⁻³ m/s² against buoyancy g ρ′/ρ₀ ≈ 9.8 × 10⁻³ (for ρ′/ρ₀ =
    10⁻³) and against g itself: 1.7 × 10⁻⁴ of g. *explain:* 4 numbered points.
14. `nb.check_agree` — **from scratch:** `f_mine = 2*GFD.OMEGA_EARTH*np.sin(lat)`; `beta_mine =
    2*GFD.OMEGA_EARTH*np.cos(lat)/GFD.EARTH_RADIUS_MEAN`; the full Coriolis acceleration `2*np.cross(GFD.earth_rotation_local(lat),
    [14.0, 0.0, 0.09])` against `GFD.coriolis_acceleration_local(..., thin=False)`; `assert np.allclose` for each; then the
    relative size of the difference to `thin=True` printed.
15. `nb.figure` — the term bars (log axis): one group per thin-shell equation, bars coloured by convention 10 (acceleration
    purple, Coriolis teal, pressure orange, friction rose, the dropped pieces hatched grey). *see:* in the horizontal, Coriolis
    and pressure tower over everything; in the vertical, pressure and weight do. *read:* "the two tall pairs are the two
    balances of the next blocks: geostrophic (C02) and hydrostatic". *change:* "…L shrinks to 10 km (a thunderstorm): Ro = 17,
    the acceleration bar overtakes Coriolis and rotation stops mattering".
16. `nb.note` — **N16 [C]** "**f-plane:** f held at its value at the central latitude, $f=f_0=2\Omega\sin\theta_0$ (`GFD.f_plane(lat)`).
    Used in C04–C12, C14, C16." **N17 [B]** "**β-plane:** keep the first term of a Taylor expansion in the northward distance
    y: $f=f_0+\beta y$ with $\beta\equiv(df/dy)_{\theta_0}=(df/d\theta\cdot d\theta/dy)_{\theta_0}=2\Omega\cos\theta_0/R$ (13.10), since $dy=R\,d\theta$." Reminder of P26
    (Taylor expansion, ch01–ch02) in one sentence. `GFD.beta_parameter(lat)`, `GFD.beta_parameter(inp["lat_rossby"])` *expect:*
    1.875 × 10⁻¹¹ and 2.239 × 10⁻¹¹ m⁻¹ s⁻¹. "The result is stated, not derived; its error is measured next."
17. `nb.plotly` — **F1** (`slider_figure`, slider = central latitude θ₀ from 5° to 75° in 15 steps): f exact (`2Ω sin(θ₀ +
    y/R)`), the β-plane line `GFD.beta_plane(y, lat0)` and 100 × `GFD.beta_plane_error(y, lat0)` over y = ±2000 km; trace names
    carry f₀ and β. *expect* at 35°: error +0.29 % at y = +500 km, +1.1 % at +1000 km, +4.0 % at +2000 km, +8.1 % at −2000 km.
    *see / read* ("the error grows as y², and faster toward the equator, where f₀ itself is small") */ change*.
18. `nb.md` — **Reading the three equations** (the three thin-shell equations once more, each term labelled and coloured), then
    **What would change if…** "…the flow is slow and wide, so that the acceleration and friction bars are negligible? Only
    the two tall bars are left in each horizontal equation — geostrophic balance (C02)."


---

### A.5 §13.5 Geostrophic Flow — C02, C03
#### C02 — Geostrophic balance (13.11)–(13.12)
1. `nb.section("13.5", "Geostrophic Flow", intro="**What is this section about?** The balance that rules every slow,
   large-scale flow: Coriolis force against pressure gradient. Then its two consequences — a wind that changes with height
   wherever temperature changes horizontally, and, when it does not, fluid that moves in rigid columns.")`
2. `nb.core("C02", "Geostrophic balance: $-fv=-\\dfrac1{\\rho_0}\\dfrac{\\partial p}{\\partial x}$ and $fu=-\\dfrac1{\\rho_0}\\dfrac{\\partial p}{\\partial y}$ (13.11)–(13.12)",
   question="Why doesn't the air simply flow from high to low pressure?")`
3. `nb.md` — **Plain words:** "Open any weather map. The wind arrows run *along* the isobars, with low pressure on their
   left in the northern hemisphere. A ball on a hill rolls downhill; air on a pressure 'hill' goes round it. The reason is
   that on a turning earth anything that moves is pushed sideways, and the push grows with speed. The air speeds up toward
   low pressure, is turned, and keeps being turned until the sideways push exactly cancels the pressure force. From then on
   it has no reason to change. Oceanographers use the same balance backwards: an altimeter measures the slope of the sea
   surface, and the slope gives the current."
4. `nb.md` — **The idea** (ASCII, three frames, northern hemisphere): `L ← ● (pressure force only: starts toward L)` →
   `● turning right, Coriolis grows with speed` → `● moving along the isobar: Coriolis (right of motion) = − pressure force`.
   Caption: "in the northern hemisphere (f > 0); mirror for f < 0".
5. `nb.primer("reading a pressure map: isobars, the pressure-gradient force points from high to low, tight spacing = strong
   force", "An isobar joins points of equal pressure, like a height contour on a hiking map. The pressure-gradient force
   per unit mass is $-\\nabla p/\\rho_0$: it points straight across the isobars from high to low, and it is large where the lines
   are close.", code="dp, dn, rho0 = 400.0, 3.0e5, 1.2      # 4 hPa between isobars 300 km apart; air density\nprint(dp / dn / rho0)
   # force per unit mass: 1.1e-3 m/s^2 (compare g = 9.8)")` (**P309**)
6. **(re-enter C02)** `nb.recap("R13", "The Rossby number", "$\\mathrm{Ro}=\\dfrac{U^2/L}{fU}=\\dfrac{U}{fL}$ (13.13): the acceleration of the flow divided by the Coriolis
   acceleration. Small Ro = rotation rules. ⚠️ **Trap T17:** Chapter 4's `SIM.rossby_number(U, Omega, l)` is $U/(2\\Omega l)$; the
   two differ by sin θ.", where="Ch. 4 §4.7")` + `nb.code`: `GFD.rossby_number(14.0, f35, 1.4e6)`, `SIM.rossby_number(14.0,
   GFD.OMEGA_EARTH, 1.4e6)` *expect:* 0.120 and 0.0686 (ratio sin 35° = 0.574).
7. `nb.note` — **N22 [B]** $E=\dfrac{\rho\nu U/L^2}{\rho fU}=\dfrac{\nu}{fL^2}$ (13.18): the Ekman number, friction over Coriolis. "Geostrophy needs
   both Ro ≪ 1 and E ≪ 1." `GFD.ekman_number(7.0, f35, 1000.0)` (vertical friction over a 1 km depth) and
   `GFD.ekman_number(3e4, f35, 1.4e6)` (horizontal) *expect:* 0.084 and 1.8 × 10⁻⁴. "The first is not small: within the lowest
   kilometre friction matters — the Ekman layer of C06. In C04 we meet E again as $(\delta/L)^2/2$."
8. `nb.derivation("D03", …)` — Part F D03 (6 steps), ref "13.11": includes **N18 [B]** (steps 5–6): $\mathbf u\cdot\nabla p=0$, and for
   constant f the pressure is a stream function, $\psi=p/(f\rho_0)$ with $u=-\partial\psi/\partial y$, $v=\partial\psi/\partial x$. `> ⚠️ Common confusion (trap T14):`
   "Chapters 4 and 11 used $u=+\partial\psi/\partial y$. Here the sign is the other way, so that ψ is high where p is high (for f > 0)." `GFD.geostrophic_streamfunction(p, f35, 1.2)` returns it.
9. `nb.worked_example("the wind from two isobars", "Isobars 4 hPa (400 Pa) apart lie 400 km apart; ρ₀ ≈ 1.25 kg/m³; f ≈ 10⁻⁴ s⁻¹
   (about 43° N). 1. Pressure-gradient force per unit mass: 400/(4 × 10⁵ × 1.25) = 8 × 10⁻⁴ m/s². 2. Balance $f\\,U=\\frac1{\\rho_0}\\lvert\\nabla p\\rvert$:
   U = 8 × 10⁻⁴/10⁻⁴ = 8 m/s. 3. Direction: along the isobars, low pressure on the left (northern hemisphere, f > 0; on the
   right for f < 0). 4. Halve the spacing → 16 m/s. 5. The same isobars at 20° N (f ≈ 5 × 10⁻⁵): 16 m/s — and at the equator
   the formula gives infinity: geostrophy fails there.")`
10. `nb.code` — `u, v = GFD.geostrophic_velocity(0.0, 400.0/3.0e5, f35, 1.2)` (pressure rising northward by 4 hPa in 300 km);
    the same with `-f35`; `GFD.geostrophic_from_height` for a sea-surface slope of 0.1 m per 100 km. *expect:* u = −13.28 m/s
    (an easterly: high pressure to the north, low on the left of the motion) for f > 0 and +13.28 m/s for f < 0; ocean
    current 0.117 m/s. *explain:* 3 numbered points.
11. `nb.check_agree` — **from scratch:** `pc = ch13.pressure_centre(x, y, dp=400.0, R_c=6e5, lat_rad=lat, kind="low")`; centred
    differences written with slices (`dpdx[:, 1:-1] = (p[:, 2:] - p[:, :-2])/(2*dx)`); `u_mine = -dpdy/(rho0*f35)`, `v_mine =
    dpdx/(rho0*f35)`; `assert np.allclose(u_mine[1:-1, 1:-1], GFD.geostrophic_from_field(p, dx, dy, f35, rho0)[0][1:-1, 1:-1])`;
    then `np.abs(u*dpdx + v*dpdy).max()` printed (round-off: the wind is along the isobars).
12. **(re-enter C02)** `nb.recap("R14", "Highs and lows", "Round a low the geostrophic wind turns counter-clockwise in the northern hemisphere
    (f > 0) and clockwise in the southern (f < 0); round a high the other way. The pressure force $-\\nabla p/\\rho_0$ points
    inward for a low; the Coriolis force $-f\\,\\mathbf e_z\\times\\mathbf u$ must point outward.", where="Ch. 4 §4.7")` + `nb.figure`: 2 × 2 panels (low /
    high × 35° N / 35° S) from `ch13.pressure_centre`, isobars orange, wind quiver teal, the sense of rotation written in each
    title. *see / read* ("low on the left of the arrows in the top row, on the right in the bottom row") */ change* ("…the
    pressure difference doubles: every arrow doubles, no arrow turns").
13. `nb.md` — **Where geostrophy fails, and how it is set up** (the verbal part of **N19**; its quantitative part is taught in
    C12): "It fails near the equator (f → 0), where friction matters (E not small), and when the flow changes in less than
    about a day (the acceleration is not small). How it is reached: release a parcel from rest in a uniform pressure
    gradient and it does not go to the low — it circles." + `nb.figure`: parcel paths from `GFD.parcel_adjust(t, G, f35,
    r=0.0)` and `r = 0.3*f35` (integrated to positions by `cumulative_trapezoid`): cycloid-like inertial loops drifting
    along the isobars at the geostrophic speed; with drag, a spiral that settles to a steady drift at an angle to the
    isobars. *expect:* mean drift of the r = 0 path = the geostrophic velocity `G/f` to 1 % over 10 inertial periods. *see /
    read / change* ("…f < 0: the loops turn the other way and the drift reverses").
14. `nb.explainer("geostrophic_balance", heading="Why doesn't air flow straight from high to low pressure?", why="Releasing
    a parcel and watching the Coriolis arrow grow and turn with the velocity until the bars cancel shows the balance being
    reached; a static map shows only the end state.", tries=["Press play with the default low and watch the teal Coriolis
    arrow swing round until it opposes the orange pressure arrow.", "Switch to S: the circulation reverses at once.", "Drag
    the latitude toward 5°: Ro climbs past 1 and the parcel cuts across the isobars.", "Turn on drag: the wind settles at
    an angle to the isobars, toward low pressure — a preview of C06."])`
15. `nb.md` — **What would change if…** "…the density, and so the pressure pattern, differs from one height to the next? Then
    the geostrophic wind differs too: the thermal wind (C03)."

#### C03 — Thermal wind (13.15), and the Taylor–Proudman theorem (13.21)
1. `nb.core("C03", "Thermal wind: $\\dfrac{\\partial v}{\\partial z}=-\\dfrac{g}{\\rho_0f}\\dfrac{\\partial\\rho}{\\partial x}$ and $\\dfrac{\\partial u}{\\partial z}=\\dfrac{g}{\\rho_0f}\\dfrac{\\partial\\rho}{\\partial y}$ (13.15)",
   question="Why is there a jet stream high up, and why is it westerly in both hemispheres?")`
2. `nb.md` — **Plain words:** "Airliners flying east across the Atlantic ride a river of air moving at 150 km/h or more; the
   wind at the ground below is far gentler. The jet sits above the place where the temperature changes fastest from south
   to north. That is no coincidence. Cold air is dense, so pressure falls faster with height in the cold column than in the
   warm one; the pressure difference between the columns therefore grows with height, and with it the geostrophic wind.
   The temperature gradient does not set the wind — it sets how the wind *changes with height*."
3. `nb.md` — **The idea** (ASCII section, pole on the left): a warm tall column and a cold short column between the same
   two pressure surfaces; the upper surface slopes down toward the pole; arrow ⊙ (out of the page = eastward) growing with
   height. Two-row table: northern hemisphere (cold to the north, f > 0 → westerly shear) · southern hemisphere (cold to
   the south, f < 0 → westerly shear again).
4. **(re-enter C03)** `nb.recap("R15", "The perturbation is hydrostatic", "$0=-\\dfrac{\\partial p}{\\partial z}-g\\rho$ (13.14). ⚠️ **Trap T3 again:** p and ρ are the
   perturbations from the state of rest (primes dropped).", where="Ch. 1 hydrostatics $dp/dz=-\\rho g$ (1.8); Ch. 7 $p'=\\rho g\\eta$ (7.52)")`
5. `nb.md` — one-line reminders: P121 (mixed partial derivatives commute, ch04) and P177 (eliminating a variable by
   cross-differentiation, ch07) — "the only two tools D04 needs".
6. `nb.derivation("D04", …)` — Part F D04 (7 steps), ref "13.15": the last two steps are **N20 [B]**, **not in the book —
   ours** (the book states it in words): with $\rho'/\rho_0=-\alpha T'$, $\dfrac{\partial u}{\partial z}=-\dfrac{g\alpha}{f}\dfrac{\partial T}{\partial y}$ and $\dfrac{\partial v}{\partial z}=\dfrac{g\alpha}{f}\dfrac{\partial T}{\partial x}$, where α = 1/T₀ for a perfect
   gas and T is potential temperature in a deep layer.
7. `nb.worked_example("a jet from a temperature gradient", "Take g ≈ 10 m/s², T₀ = 250 K so gα = g/T₀ = 0.04 m s⁻² K⁻¹, f = 10⁻⁴
   s⁻¹, and let temperature fall northward by 1 K per 100 km: ∂T/∂y = −10⁻⁵ K/m. 1. Shear: $\\frac{\\partial u}{\\partial z}=-\\frac{g\\alpha}{f}\\frac{\\partial T}{\\partial y}$ = −(0.04/10⁻⁴) ×
   (−10⁻⁵) = +4 × 10⁻³ s⁻¹ = 4 m/s per km. 2. With calm air at the ground, the wind at 10 km is 40 m/s from the west. 3.
   Southern hemisphere: ∂T/∂y = +10⁻⁵ (cold toward the south pole) and f = −10⁻⁴: the two sign changes cancel, again +4 m/s
   per km. 4. No gradient → no shear.")`
8. `nb.code` — `dudz, dvdz = GFD.thermal_wind_from_temperature(0.0, -7.0e-6, f35, alpha=1/280.0)`; `U9 = dudz*9000.0`; southern
   twin `GFD.thermal_wind_from_temperature(0.0, +7.0e-6, -f35, alpha=1/280.0)`; density form
   `GFD.thermal_wind_shear(0.0, 1.0e-5, f35, 1027.0)` (an ocean front: density rising northward by 1 kg m⁻³ per 100 km).
   *expect:* dudz = 2.931 × 10⁻³ s⁻¹; U at 9 km = 26.4 m/s; southern twin identical; ocean shear 1.14 × 10⁻³ s⁻¹ (1.14 m/s over
   1000 m). *explain:* 4 numbered points, one on why `alpha` must be passed by name.
9. `nb.check_agree` — **from scratch:** shear in one line `-G0*(1/280.0)*(-7.0e-6)/f35`; then U(z) by a hand-written
   cumulative trapezoid loop over a density gradient that decays with height; `assert np.isclose(...)` against
   `GFD.thermal_wind_from_temperature` and `assert np.allclose(U_loop, GFD.thermal_wind_integrate(z, drho_dy, f35, 1.2))`.
10. `nb.figure` — `sec = ch13.jet_section(y, z, dT=28.0, width=2.0e6, lat_rad=lat, alpha=1/280.0)`: isotherms (blue contours)
    and the zonal wind (filled contours), jet core marked. *expect:* core wind at z = 9 km ≈ 26 m/s above y = 0 (where
    ∂T/∂y = −7 × 10⁻⁶ K/m). *see:* the wind maximum sits above the tightest isotherms; *read:* "follow one isotherm: where it
    slopes most steeply the wind contours are closest in the vertical"; *change:* "…the surface were already blowing at 5
    m/s: every contour shifts by 5 m/s — the thermal wind fixes the shear, not the wind".
11. `nb.plotly` — **F3** (`slider_figure`, slider = temperature contrast dT from 0 to 40 K, 21 steps; dropdown hemisphere):
    U(z) at the front's centre from `GFD.thermal_wind_integrate` with the isotherm slope in the trace name. *see / read /
    change*.
12. `nb.note` — **N21 [B]** the rotating-tank pair: $-2\Omega v=-\dfrac1\rho\dfrac{\partial p}{\partial x}$ (13.16) and $2\Omega u=-\dfrac1\rho\dfrac{\partial p}{\partial y}$ (13.17), i.e. geostrophy
    with f = 2Ω in a homogeneous fluid. > ⚠️ **slip #3 — the book prints** $-2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$**; the correct form is**
    $2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$ (it is the second geostrophic equation with f = 2Ω; with the printed sign the next step of
    the book cannot be reached). `ch13.taylor_proudman_sympy()["ok"]` → True; `ch13.taylor_proudman_sympy(printed=True)["ok"]` →
    False.
13. `nb.derivation("D05", …)` — Part F D05 (7 steps), ref "13.21": **N23 [B]** $\dfrac{\partial w}{\partial z}=0$ (13.19) is its step 4; **N24 [B]**
    $\dfrac{\partial v}{\partial z}=\dfrac{\partial u}{\partial z}=0$ (13.20) its step 6; **N25 [B]** the theorem $\partial\mathbf u/\partial z=0$ (13.21) with its four hypotheses (steady,
    Ro ≪ 1, E ≪ 1, uniform density): "thermal wind with nothing to drive a shear".
14. `nb.note` — **N26 [B]** Taylor's experiment as our sketch: `draw_taylor_column()` (2-D flow round a cylinder from
    `core.potential` repeated at every height; the column above the short obstacle shaded). "A short obstacle on the floor of
    a rapidly rotating tank acts as if it reached the surface: dye divides ahead of the column and goes round it." +
    `GFD.taylor_proudman_residual(lambda x, y, z: (-y, x, 0.0), (0.3, 0.2, 0.5))` → all three derivatives 0. (The tank's
    dimensions are not reproduced.)
15. `nb.explainer("thermal_wind", heading="Why is there a jet stream, and why is it westerly in both hemispheres?",
    why="Dragging the temperature contrast tilts the isotherms and grows the shear at once; flipping the hemisphere flips
    both f and the gradient and the jet stays westerly — three coupled changes a single section cannot show.",
    tries=["Drag the contrast to zero: every wind profile collapses to a vertical line (the Taylor–Proudman limit).",
    "Switch N ↔ S and watch which two signs change.", "Click a point in the section to see its local gradient and shear
    worked out.", "Switch to the ocean mode: the same relation in density, with a front instead of a jet."])`
16. `nb.md` — **What would change if…** "…friction is not negligible? Near the sea surface and near the ground it is not;
    the balance gains a third force and the flow turns — the Ekman layers (C04–C06). (The sheared flow of this block is also
    the starting state of the Eady problem, C16.)"

---

### A.6 §13.6 Ekman Layer at a Free Surface — C04, C05
#### C04 — The surface Ekman spiral (13.27)–(13.29)
1. `nb.section("13.6", "Ekman Layer at a Free Surface", intro="**What is this section about?** What a steady wind does to
   the top of the ocean. The answer is a thin layer in which the current turns and weakens with depth — and whose total
   transport is at right angles to the wind. ⚠️ In this section z = 0 is the sea surface and the ocean is z < 0.")`
2. `nb.core("C04", "The surface Ekman spiral: $\\dfrac{d^2V}{dz^2}=\\dfrac{if}{\\nu_v}V$, $V\\equiv u+iv$ (13.27), with $V=A\\,e^{(1+i)z/\\delta}+B\\,e^{-(1+i)z/\\delta}$,
   $\\delta=\\sqrt{2\\nu_v/f}$ (13.28)–(13.29)", question="A steady wind blows over the sea. Which way does the water go, and how deep
   does the motion reach?")`
3. `nb.md` — **Plain words:** "Nansen noticed that Arctic ice drifts not downwind but well to the right of the wind. His
   student Ekman explained it in 1905. The wind drags the top layer of water; that layer is turned by the Coriolis force
   and drags the layer below, which is turned a little more, and so on downward. The result is a staircase of currents
   that rotates and fades with depth: a spiral."
4. `nb.md` — **The idea:** a stack of thin slabs; each slab feels the drag of the slab above (pushing it), the drag of the
   slab below (holding it back) and the Coriolis force (turning it). Steady state: the net drag on each slab balances its
   Coriolis force. Table "what can balance friction": → next row.
5. **(re-enter C04)** `nb.recap("R16", "Three ways to balance friction", "Friction can be balanced by unsteadiness (the layer thickens as
   $\\delta\\sim\\sqrt{\\nu t}$), by advection (it thickens downstream, $\\delta\\sim\\sqrt{\\nu x/U}$), or — new here — by the Coriolis force, and then the
   thickness $\\delta=\\sqrt{2\\nu_v/f}$ does **not** grow at all.", where="Ch. 8 §8.4 (impulsively started plate), Ch. 9 §9.3 (Blasius layer)")` +
   `ch13.viscous_layer_thicknesses(0.03, t=86400.0, x=1.0e5, U=0.1, f=f60)` *expect:* 50.9 m after one day, 173 m after 100 km
   at 0.1 m/s, 21.8 m for the Ekman layer — a three-row table.
6. `nb.note` — **N27 [B]** the balance: $-fv=\nu_v\dfrac{d^2u}{dz^2}$ and $fu=\nu_v\dfrac{d^2v}{dz^2}$ (13.22)–(13.23) — steady, horizontally uniform, no
   pressure gradient: only Coriolis (teal) and vertical friction (rose) remain. **N28 [B]** the conditions: $\rho\nu_v\dfrac{du}{dz}=\tau$ at
   z = 0 (13.24), $\dfrac{dv}{dz}=0$ at z = 0 (13.25), and $u,v\to0$ as $z\to-\infty$ (13.26, corrected). > ⚠️ **slip #4 — the book
   prints** $z\to\infty$**; the correct form is** $z\to-\infty$: the ocean lies below the surface z = 0.
7. `nb.primer("hodograph: the tip of the velocity vector traced as depth (or time) changes", "Plot v against u for every
   depth and join the points. The curve is the hodograph: each point is the tip of the velocity arrow at one depth. A
   straight hodograph means the current keeps its direction; a spiral means it turns.", code="s = np.linspace(0, 3, 4)
   # four depths, in units of delta\nV = np.exp(-s) * np.exp(-1j*s)               # speed decays, direction turns clockwise\nprint(np.round(V.real,
   2), np.round(V.imag, 2))  # the points of a spiral")` (**P310**)
8. `nb.primer("complex velocity V = u + iv: multiplying by i turns a vector 90° to the left; e^{(1+i)s} decays and turns at
   the same rate", "Pack the two velocity components into one complex number V = u + iv. Then iV is the same arrow
   turned 90° counter-clockwise — exactly what the Coriolis terms do — so two coupled real equations become one complex
   one. A factor $e^{(1+i)s}=e^{s}e^{is}$ changes the length by $e^s$ and the direction by s radians together.",
   code="V = 3 + 4j                                # u = 3, v = 4\nprint(1j*V)                                 # (-4+3j): turned 90 deg to the left\nprint(abs(np.exp((1+1j)*(-1.0))),
   np.angle(np.exp((1+1j)*(-1.0))))  # 0.368, -1.0 rad\nprint(np.sqrt(1j))                          # (0.707+0.707j) = (1+i)/sqrt(2)")` (**P311**) — reminder
   in one sentence of P45/P159 (complex square roots) and P44 (try $e^{rz}$ in a linear ODE).
9. `nb.derivation("D06", …)` — Part F D06 (14 steps), ref "13.27": **N29 [B]** is its step 2 ($\dfrac{d^2V}{dz^2}=\dfrac{if}{\nu_v}V$ (13.27)); **N30
   [B]** its steps 3–6 ($V=Ae^{(1+i)z/\delta}+Be^{-(1+i)z/\delta}$ with $\delta=\sqrt{2\nu_v/f}$ (13.28)–(13.29)). `> ⚠️ Common confusion (trap T6):`
   "δ is the depth over which the current falls by e (and turns by one radian). Oceanographers' 'Ekman depth' is πδ, where
   the current first points opposite to the surface current. `GFD.ekman_depth(nu_v, f, convention='efold' or 'pi')`. The
   laminar cousin is Chapter 8's Stokes layer, $u=Ue^{-y\sqrt{\omega/2\nu}}\cos(\omega t-y\sqrt{\omega/2\nu})$ (8.38), with thickness $\sqrt{2\nu/\omega}$: replace the
   oscillation frequency by f."
10. `nb.worked_example("thickness and surface speed", "Take ν_v = 0.05 m²/s, f = 10⁻⁴ s⁻¹, τ = 0.1 N/m², ρ = 1000 kg/m³. 1.
    $\\delta=\\sqrt{2\\nu_v/f}$ = √(0.1/10⁻⁴) = √1000 ≈ 31.6 m. 2. Ekman depth πδ ≈ 99 m. 3. Surface speed $V_0=\\dfrac{\\tau/\\rho}{\\sqrt{f\\nu_v}}$ = 10⁻⁴/√(5 × 10⁻⁶) =
    10⁻⁴/2.24 × 10⁻³ ≈ 0.045 m/s. 4. Direction: 45° to the right of the wind (northern hemisphere, f > 0; to the left for
    f < 0), so u = v-magnitude = 0.045/√2 ≈ 0.032 m/s. 5. At z = −δ: speed × e⁻¹ = 0.016 m/s, turned a further 57°.")`
11. `nb.code` — `delta = GFD.ekman_depth(0.03, f60)`; `GFD.ekman_depth(0.03, f60, convention="pi")`; `u, v = GFD.ekman_surface(z,
    0.07, 0.0, 1027.0, 0.03, f60)` on `z = np.linspace(-5*delta, 0, 201)`; the southern twin with `-f60`;
    `ch13.ekman_surface_printed(z, 0.07, 1027.0, 0.03, f60)` compared. *expect* (60° N, τ = 0.07 N/m² eastward, ν_v = 0.03 m²/s,
    ρ = 1027 kg/m³): δ = 21.80 m, πδ = 68.47 m; surface current (u, v) = (+0.02476, −0.02476) m/s, speed 0.03502 m/s, 45° to
    the right of the wind; at 60° S (+0.02476, +0.02476): 45° to the left; speed at z = −πδ is 4.3 % of the surface speed and
    points opposite to it; printed form equal to round-off for f > 0. *explain:* 5 numbered points.
12. `nb.check_agree` — **from scratch:** `A = 0.07*delta*(1 - 1j)/(2*1027.0*0.03)`; `V = A*np.exp((1 + 1j)*z/delta)`; `assert
    np.allclose(V.real, u) and np.allclose(V.imag, v)`; then `GFD.ekman_residual(z, u, v, 0.03, f60)` max abs printed (≈ 10⁻⁹
    relative: finite-difference error only).
13. `nb.figure` — **N31 [B]** our Fig. 13.6 (`ch13.fig_ekman_surface()`): (a) hodograph with depth marks at −z/δ = 0, π/8,
    π/4, 3π/8, π/2, π, the wind arrow, the 45° angle; (b) u and v against −z/δ. + `nb.plotly`: the spiral in 3-D (`go.Cone` or
    line + arrows at 16 depths, dropdown N / S). *see:* the arrow shortens and swings clockwise going down (northern
    hemisphere); *read:* "the surface arrow is 45° to the right of the wind; one δ down it has turned one more radian";
    *change:* "…ν_v were four times larger: the spiral reaches twice as deep and the surface current halves (F4)".
14. `nb.plotly` — **F4** (`slider_figure`, slider = ν_v from 0.005 to 0.3 m²/s, 25 log steps): hodograph (trace 1) and u(z),
    v(z) (traces 2–3); trace names carry δ, the surface speed and "transport τ/(ρf) = 0.540 m²/s" — which does not change.
    *see / read / change.*
15. `nb.note` — **N33 [B]** "A horizontal pressure gradient merely adds a depth-independent geostrophic velocity (the problem
    is linear): `GFD.ekman_surface(z, 0.07, 0.0, 1027.0, 0.03, f60, U_g=0.05)` shifts the whole hodograph by 0.05 m/s."
    (small figure: two hodographs.)
16. `nb.note` — **N35 [B]** why the layer does not grow: with $\omega_x=-\dfrac{dv}{dz}$, $\omega_y=\dfrac{du}{dz}$, the balance is $-f\dfrac{dv}{dz}=\nu_v\dfrac{d^2\omega_y}{dz^2}$
    and $-f\dfrac{du}{dz}=\nu_v\dfrac{d^2\omega_x}{dz^2}$ (13.31): tilting of the planetary vorticity (the $(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u$ term of Chapter 5's
    $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\frac1{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega$ (5.30)) balances the downward diffusion of horizontal vorticity. The result is
    stated; `b = GFD.ekman_vorticity_balance(z, 0.07, 1027.0, 0.03, f60)`, figure of `b["tilt_x"]` and `b["diff_x"]`
    (identical curves; difference at round-off).
17. `nb.md` — **What would change if…** "…we do not care about the shape of the spiral, only about how much water moves?
    Then we should add the arrows up — and the eddy viscosity drops out (C05)."

#### C05 — Ekman transport (13.30), pumping and upwelling
1. `nb.core("C05", "Ekman transport: $\\displaystyle\\int_{-\\infty}^{0}u\\,dz=0$ and $\\displaystyle\\int_{-\\infty}^{0}v\\,dz=-\\dfrac{\\tau}{\\rho f}$ (13.30)",
   question="Summed over the whole layer, where does the wind-driven water go — and what happens at a coast?")`
2. `nb.md` — **Plain words:** "Off Peru and California the wind blows toward the equator along the coast, and the water at
   the shore is cold and full of fish. The wind is not pulling the cold water up; it is moving the surface water
   *offshore*, at right angles to itself, and deeper water has to rise to replace it. The same sideways transport, driven
   by the trade winds and the westerlies, piles water up in the middle of each ocean basin and sets the great gyres
   turning."
3. `nb.md` — **The idea:** "Add up all the arrows of the spiral. The along-wind parts cancel; what is left points 90° to
   the right of the wind (northern hemisphere, f > 0; to the left for f < 0). Reason in one line: over the whole layer the
   only horizontal forces are the wind stress on top and the Coriolis force on the total transport — they must cancel, and
   the Coriolis force is at right angles to the transport."
4. `nb.primer("depth-integrated quantities: transport per unit width [m²/s] and its divergence as a vertical velocity",
   "Integrating a velocity over depth gives a transport per unit width, in m²/s: the volume crossing a one-metre-wide
   gate each second. If more leaves a column sideways than enters, the difference must come through its top or bottom:
   the divergence of the transport [m/s] is a vertical velocity.", code="z = np.linspace(-100, 0, 1001)            # depth
   [m]\nu = 0.1*np.exp(z/20)                       # a current that decays with depth [m/s]\nprint(np.trapezoid(u, z))                  #
   2.0 m^2/s (exactly 0.1*20*(1 - e^-5))")` (**P313**) + `nb.primer("first integral: integrating an ODE once across a layer to
   get a transport without solving it", "If an equation says 'something = d(flux)/dz', integrating across the layer
   gives 'total something = flux at the top − flux at the bottom'. You need the flux only at the two ends, not the
   solution in between.", code="tau_top, tau_bottom, rho, f = 0.07, 0.0, 1027.0, 1.0e-4\nprint(-(tau_top - tau_bottom)/(rho*f))
   # the y-transport, with no profile in sight")` (**P312**)
5. `nb.derivation("D07", …)` — Part F D07 (8 steps), ref "13.30": its second route (steps 6–8) is **N32 [B]**: integrate
   $-\rho fv=\dfrac{d\tau}{dz}$ once; only "the stress vanishes at depth" is used, so the result does not depend on the eddy-viscosity
   model.
6. `nb.worked_example("how much water, and which way", "τ = 0.1 N/m² toward the east, ρ = 1000 kg/m³, f = 10⁻⁴ s⁻¹. 1. Transport
    $\\int v\\,dz=-\\tau/(\\rho f)$ = −0.1/(1000 × 10⁻⁴) = −1 m²/s: one cubic metre per second through every metre of an east–west
    line, toward the **south** (to the right of the wind). 2. Along a 1000 km line: 10⁶ m³/s = 1 Sverdrup. 3. Southern
    hemisphere (f = −10⁻⁴): +1 m²/s, toward the north — to the left of the wind. 4. Double the eddy viscosity: same answer.")`
7. `nb.code` — `Mx, My = GFD.ekman_transport(0.07, 0.0, 1027.0, f60)`; `GFD.ekman_transport_partial(-delta, 0.07, 0.0, 1027.0,
   0.03, f60)`; the same at `-np.pi*delta`. *expect:* (0, −0.5397) m²/s; the layer above z = −δ carries 86 % of the transport's
   magnitude (abs(1 − e^{−(1+i)}) = 0.859) but not yet in the final direction; above z = −πδ, 104 % (a small overshoot that
   the deeper, opposed currents remove). *explain:* 3 numbered points.
8. `nb.check_agree` — **from scratch:** for ν_v = 0.03 and 0.12 m²/s: `u, v = GFD.ekman_surface(z, …)` on z from −40δ to 0,
   `np.trapezoid(v, z)`, `np.trapezoid(u, z)`; `assert np.isclose(My_num, My, rtol=1e-4)` for both, `abs(Mx_num) < 1e-6`.
   *expect:* δ doubles (21.80 → 43.59 m), the surface speed halves (0.03502 → 0.01751 m/s), the transport stays −0.5397 m²/s.
9. `nb.note` — **N37 [B]** **Not in the book — ours** (the book remarks that the constant-viscosity spiral is rarely
   observed): the layer with a depth-dependent eddy viscosity, $\dfrac{d}{dz}\Big(K_v\dfrac{dV}{dz}\Big)=if\,(V-V_g)$, solved numerically. `> 🔧 **Our
   choice** — a second-order finite-volume scheme on a stretched grid with $K_v$ at the cell faces, one complex tridiagonal
   solve (the pattern of Chapter 12's eddy-viscosity solver).` `V = GFD.ekman_solve(z, K, f60, tau=0.07+0j, rho=1027.0)` for
   three profiles (constant 0.03; 0.03(1 + abs(z)/20); 0.003 + 0.03 e^{z/30} m²/s). + `nb.figure`: three hodographs (very
   different shapes, surface angles away from 45°) and a bar chart of the three transports. *expect* (hypothesis from D07's
   second route; the designer did not run this solver): the three transports agree with −0.5397 m²/s to better than 1 %; the
   constant case matches `GFD.ekman_surface` to 10⁻³ relative. *see / read* ("the shape belongs to the turbulence model,
   the transport to Newton's law") */ change*.
10. `nb.derivation("D09", …)` — Part F D09 (7 steps), ref "": **N36 [B]** **Not in the book — ours** (the book gives only the
    consequence, rising air in a low): Ekman pumping $w_E=\dfrac1\rho\,\mathbf e_z\cdot\nabla\times\Big(\dfrac{\boldsymbol\tau}{f}\Big)$ under a surface layer, and
    $w=\dfrac\delta2\zeta_g$ above a bottom layer. One line, labelled beyond the book: "**Sverdrup balance** — in the interior below the
    surface layer, stretching by this pumping is balanced by moving columns north or south, $\beta V=f\,w_E$ per unit width
    (depth-integrated, with no vertical velocity at depth): the wind's curl sets the north–south transport of a whole gyre."
11. `nb.code` — `wE = GFD.ekman_pumping(tau_x, tau_y, dx, dy, 1027.0, f35)` for a zonal wind `tau_x = -0.07*np.cos(np.pi*Y/Ly)`
    (easterlies in the south, westerlies in the north of a 2000 km band, the pattern over a subtropical gyre); coastal case
    `ch13.coastal_upwelling(-0.07, "east", lat)`. *expect:* pumping downward everywhere in the band (anticyclonic stress curl
    for f > 0), peak magnitude π × 0.07/(2 × 10⁶ × 1027 × f35) = 1.28 × 10⁻⁶ m/s ≈ 40 m per year; coastal: offshore transport
    0.815 m²/s (35° N) → upwelling of 2.0 × 10⁻⁵ m/s ≈ 1.8 m per day if it is fed over a 40 km wide strip. + `nb.figure`: map
    of τ arrows and w_E colour. *see / read / change* ("…in the southern hemisphere the same wind band pumps upward").
12. `nb.note` — **N50 [C]** "The wind-driven circulation and Stommel's explanation of strong western boundary currents are
    named by the book, not presented. Pointer: the Sverdrup line above; any dynamical-oceanography text for the gyre."
13. `nb.live` — free τ, ν_v, latitude (both signs): hodograph + transport arrow (`live(fn, tau=(0.01, 0.3, 0.01), nu_v=(0.005,
    0.3, 0.005), lat_deg=(-80, 80, 5))`); paired with F4 for the published page.
14. `nb.explainer("ekman_spiral", heading="The wind blows east — why does the water go south?", why="Orbiting the spiral in
    3-D, sliding a depth cursor down the hodograph and watching the running transport swing round to 90° connects three
    pictures of one solution; changing ν_v changes δ and the surface speed and leaves the transport untouched.",
    tries=["Press play: the depth cursor descends and the running transport arrow swings to 90° from the wind.", "Use the
    preset 'four times ν_v': twice as deep, half as fast at the surface, same transport.", "Switch to S: surface current and
    transport flip to the left of the wind.", "Choose 'wind along a coast' and read the status: upwelling or downwelling?"])`
15. `nb.md` — **What would change if…** "…the boundary is solid and the fluid above it is already moving? The same
    equation with a different boundary condition: the bottom Ekman layer (C06)."

---

### A.7 §13.7 Ekman Layer on a Rigid Surface — C06
#### C06 — The bottom Ekman layer (13.41)
1. `nb.section("13.7", "Ekman Layer on a Rigid Surface", intro="**What is this section about?** The friction layer under a
   geostrophic flow — the lowest kilometre of the atmosphere, or the bottom of the ocean. ⚠️ Here z = 0 is the solid
   surface and the fluid is z > 0.")`
2. `nb.core("C06", "The bottom Ekman layer: $u=U\\big[1-e^{-z/\\delta}\\cos(z/\\delta)\\big]$, $v=Ue^{-z/\\delta}\\sin(z/\\delta)$ (13.41)",
   question="If the wind follows the isobars, why does air spiral into a low?")`
3. `nb.md` — **Plain words:** "Satellite pictures show clouds spiralling *into* a depression, not circling it. Aloft the
   wind does follow the isobars. But near the ground friction slows the air; slower air feels a weaker Coriolis force;
   the pressure force, which does not care about speed, is no longer cancelled, and the air drifts toward low pressure.
   Air converging on a low has nowhere to go but up — clouds and rain. Stir a cup of tea and the leaves collect in the
   middle for the same reason (Ch. 9)."
4. `nb.md` — **The idea** (force triangle, ASCII): aloft: pressure (orange, toward low) = − Coriolis (teal). Near the
   ground: Coriolis shorter and turned, friction (rose) closes the triangle; the wind has a component toward low pressure.
5. `nb.note` — **N38 [C]** the interior is geostrophic, $fU=-\dfrac1\rho\dfrac{dp}{dy}$ (13.32): it enters as the constant U. **N39 [B]** so the
   balance in the layer has three forces: $-fv=\nu_v\dfrac{d^2u}{dz^2}$ and $fu=\nu_v\dfrac{d^2v}{dz^2}+fU$ (13.33)–(13.34) — the pressure gradient is hidden
   in the term fU. **N40 [C]** far-field condition $u=U,\ v=0$ as $z\to\infty$ (13.35); **N41 [C]** no slip $u=0,\ v=0$ at z = 0 (13.36).
6. `nb.derivation("D08", …)` — Part F D08 (9 steps), ref "13.41": **N42 [B]** $\dfrac{d^2V}{dz^2}=\dfrac{if}{\nu_v}(V-U)$ (13.37) is its step 2; **N43 [C]**
   $V=U$ as $z\to\infty$ (13.38) and **N44 [C]** $V=0$ at z = 0 (13.39) are its conditions; **N45 [B]** $V=Ae^{-(1+i)z/\delta}+Be^{(1+i)z/\delta}+U$ (13.40)
   its step 3; **N47 [B]** the transport $\displaystyle\int_0^\infty v\,dz=U\Big[\dfrac{\nu_v}{2f}\Big]^{1/2}=\dfrac12U\delta$ its steps 8–9.
7. `nb.worked_example("the layer under a 10 m/s wind", "U = 10 m/s, ν_v = 5 m²/s, f = 10⁻⁴ s⁻¹. 1. $\\delta=\\sqrt{2\\nu_v/f}$ = √(10/10⁻⁴) = √10⁵
   ≈ 316 m. 2. At z = δ: e⁻¹ = 0.368, cos 1 = 0.540, sin 1 = 0.841 → u = 10(1 − 0.368 × 0.540) = 8.0 m/s, v = 10 × 0.368 ×
   0.841 = 3.1 m/s: the wind points 21° to the left of the isobars' direction (toward low pressure; northern hemisphere).
   3. Very near the ground u ≈ v: 45°. 4. Transport toward low pressure: ½Uδ = 0.5 × 10 × 316 ≈ 1580 m²/s.")`
8. `nb.code` — `delta = GFD.ekman_depth(7.0, f60)`; `u, v = GFD.ekman_bottom(z, 12.0, 0.0, 7.0, f60)`;
   `GFD.ekman_bottom_transport(12.0, 0.0, 7.0, f60)`; `fb = GFD.ekman_force_balance(0.5*delta, 12.0, 7.0, f60)`. *expect* (60° N,
   U = 12 m/s, ν_v = 7 m²/s): δ = 332.9 m (πδ = 1046 m); (u, v) = (5.61, 3.49) m/s at z = 0.5δ (31.9° from the isobars), (9.61,
   3.71) at δ (21.1°), (12.00, 2.49) at πδ/2 (11.7°), (12.52, 0.00) at πδ; transport (M_x, M_y) = (−1998, +1998) m²/s;
   `fb["sum"]` at round-off. *explain:* 4 numbered points.
9. `nb.check_agree` — **from scratch:** at z = 0.5δ: `cor = f60*np.array([v0, -u0])`, `pgf = f60*np.array([0.0, 12.0])`, `fric =
   7.0*np.array([d2u, d2v])` with second differences of the profile on a fine grid; `assert np.allclose(cor + pgf + fric, 0,
   atol=1e-7)`; `assert np.allclose(cor, fb["coriolis"])`.
10. `nb.figure` — **N46 [B]** our Fig. 13.9 (`ch13.fig_ekman_bottom()`): (a) hodograph from the origin to U with marks at z/δ =
    π/4, π/2, π; (b) u and v against z/δ. **Stated precisely (the orchestrator's ruling on the overshoot):** "The component u overshoots U.
    Its largest value is at z = 3πδ/4, where $u/U=1+e^{-3\pi/4}/\sqrt2=1.067$ (6.7 % above U); at z = πδ, where v first returns to
    zero, $u/U=1+e^{-\pi}=1.043$." *expect:* argmax of u at z/δ = 2.356; 12.80 m/s; 12.52 m/s at πδ. *see / read* ("the hodograph
    leaves the origin at 45° to the left of U and winds into the point (U, 0)") */ change*.
11. `nb.note` — **N48 [B]** numbers with our inputs: laminar air (ν = 1.5 × 10⁻⁵ m²/s) would give δ = 0.49 m at 60° N, against
    an observed layer of the order of a kilometre: `GFD.eddy_viscosity_from_depth(1000.0, f60)` *expect:* 63 m²/s — "eddies,
    not molecules, carry the stress". > ⚠️ **slip #10 — the book prints** a laminar thickness about a quarter below what
    its own formula $\delta=\sqrt{2\nu/f}$ gives for its stated inputs**; the correct form is** the value of the formula (our numbers above
    are computed, not quoted).
12. `nb.note` — **N49 [B]** the force triangle at three heights (figure: `GFD.ekman_force_balance` at z = 0.1δ, δ, 3δ; arrows
    pressure orange, Coriolis teal, friction rose, closing to zero; the wind vector drawn in grey): $-f\,\mathbf e_z\times\mathbf u-\dfrac1\rho\nabla p+\nu_v\dfrac{d^2\mathbf u}{dz^2}=0$.
    "Inflow and rising air in a low, outflow and sinking in a high — in both hemispheres." + bottom pumping (the ours-item of
    C05 applied here): `GFD.ekman_pumping_bottom(1.0e-5, 7.0, f60)` *expect:* 1.66 mm/s for ζ_g = 10⁻⁵ s⁻¹; with
    `-f60` and ζ_g = −10⁻⁵ (a southern low) the same +1.66 mm/s. "Spin-down: this pumping squashes the columns above and
    removes the vorticity of a 9 km deep vortex in about 2H/(fδ) ≈ 5 days (*expect* 4.95 days)."
13. `nb.note` — **N34 [B]** "Observed current profiles in shallow water look like a surface layer on top of a bottom layer.
    The observation is not reproduced; **our** exact finite-depth solution shows the same structure." `V =
    GFD.ekman_finite_depth(z, 0.07, 0.0, 0.05, 0.0, 100.0, 0.03, 1027.0, f60)` (ours, labelled): figure of u, v against depth
    and the hodograph, for H = 100 m ≈ 4.6δ. *see / read / change* ("…H < 2δ: the two layers merge and the current runs
    nearly downwind").
14. `nb.plotly` — **F5** (`slider_figure`, slider = U from 2 to 30 m/s; dropdown ν_v = 2, 7, 20 m²/s): hodograph and the
    cross-isobar angle against z/δ; trace name carries δ and ½Uδ. *read:* "the angle curve does not depend on U at all".
15. `nb.explainer("ekman_force_balance", heading="If the wind follows the isobars, why does air spiral into a low?",
    why="Sliding the height cursor from the ground to the top of the layer turns the force triangle continuously from
    'pressure against friction' to 'pressure against Coriolis'; the triangle always closes.", tries=["Slide the height to
    0: the wind is 45° from the isobars.", "Jump to the preset z = πδ: v is zero and u is 4 % above U.", "Switch the
    centre from low to high and read the pumping sign.", "Switch to S and check that air still flows *into* the low."])`
16. `nb.md` — **What would change if…** "…we forget friction again and ask how such a layer of fluid *moves in time* — its
    waves? For that we need equations for a thin layer with a free surface (C07)."


---

### A.8 §13.8 Shallow-Water Equations — C07
#### C07 — The linear shallow-water equations (13.45)
1. `nb.section("13.8", "Shallow-Water Equations", intro="**What is this section about?** The simplest model that still has
   gravity waves and rotation: one thin layer of uniform density with a free surface. Three equations for the surface
   height and the two horizontal velocities. Sections 13.10–13.15 are all about its solutions. ⚠️ Here z = 0 is the flat
   bottom and the surface is at z = H + η.")`
2. `nb.core("C07", "The linear shallow-water equations: $\\dfrac{\\partial\\eta}{\\partial t}+H\\Big(\\dfrac{\\partial u}{\\partial x}+\\dfrac{\\partial v}{\\partial y}\\Big)=0$, $\\dfrac{\\partial u}{\\partial t}-fv=-g\\dfrac{\\partial\\eta}{\\partial x}$,
   $\\dfrac{\\partial v}{\\partial t}+fu=-g\\dfrac{\\partial\\eta}{\\partial y}$ (13.45)", question="What is the least we must keep to describe a tide or a tsunami crossing an
   ocean on a turning earth?")`
3. `nb.md` — **Plain words:** "A tsunami in mid-ocean is a few hundred kilometres long in water four kilometres deep. To
   such a wave the ocean is a puddle: the water moves almost horizontally, the same at every depth, and the pressure at
   any point is just the weight of water above it. Then the whole three-dimensional problem collapses to a map: how high is
   the surface, and which way is the column moving. Tides, storm surges, the adjustment of the ocean to a change of wind —
   and, with one substitution, each internal mode of a stratified ocean — all obey these three equations."
4. `nb.md` — **The idea** + `nb.note` **N54 [C]** the geometry sketch `draw_shallow_layer()` (H, η, z from the bottom; trap T4
   line: "which z = 0 does each section use?" as a four-row table). ASCII: `convergence of columns → surface rises →
   slope pushes columns apart again`.
5. **(re-enter C07)** `nb.recap("R17", "Hydrostatic pressure under a displaced surface", "With zero pressure at the surface z = H + η, the
   pressure at height z is the weight of the water above: $p=\\rho g(H+\\eta-z)$.", where="Ch. 7, long waves are hydrostatic:
   $p'=\\rho g\\eta$ (7.52)")`
6. **(re-enter C07)** `nb.recap("R18", "The pressure gradient is the surface slope", "$\\dfrac{\\partial p}{\\partial x}=\\rho g\\dfrac{\\partial\\eta}{\\partial x}$ and $\\dfrac{\\partial p}{\\partial y}=\\rho g\\dfrac{\\partial\\eta}{\\partial y}$ (13.42): neither depends
   on z, so every level of a column is pushed alike and a flow that starts depth-independent stays so.", where="Ch. 7 §7.2")`
7. `nb.md` — one-line reminders: P313 (transport per unit width and its divergence, primed in C05), P38 (product rule), and
   the kinematic surface condition of Ch. 7, $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}\cong\frac{\partial\eta}{\partial t}$ (7.17), "here in its exact form w(η) = Dη/Dt".
8. `nb.derivation("D10", …)` — Part F D10 (9 steps), ref "13.45": **N51 [B]** is its step 3,
   $(H+\eta)\dfrac{\partial u}{\partial x}+(H+\eta)\dfrac{\partial v}{\partial y}+w(\eta)-w(0)=0$ (13.43); **N52 [B]** its steps 4–6, the nonlinear continuity equation
   $\dfrac{\partial\eta}{\partial t}+\dfrac{\partial}{\partial x}\big[u(H+\eta)\big]+\dfrac{\partial}{\partial y}\big[v(H+\eta)\big]=0$ (13.44): "the divergence of the transport lowers the surface".
9. `nb.worked_example("how fast, and how strong a current", "H = 4000 m, g ≈ 10 m/s², no rotation. 1. Try η = F(x − ct): the two
   one-dimensional equations give c² = gH, so c = √(40 000) = 200 m/s (720 km/h — a jet aircraft). 2. An ocean 6000 km wide
   is crossed in 30 000 s ≈ 8 h. 3. A crest 1 m high carries a current u = cη/H = 200 × 1/4000 = 0.05 m/s: the water barely
   moves; the *shape* travels. 4. With H = 1.3 m instead (the equivalent depth of an internal mode, C08): c = √13 ≈ 3.6
   m/s.")`
10. `nb.primer("a C-grid shallow-water step and its CFL limit with c = √(gH) (our choice of scheme)", "To step the equations
    in time we store η at the centres of grid cells and the velocities on the cell faces (a 'C-grid'), so every difference
    is taken across exactly one cell. Forward–backward stepping: update η with the old velocities, then the velocities
    with the new η. It is stable only if a wave cannot cross more than one cell per step: $c\\,\\Delta t/\\Delta x\\le1$ with $c=\\sqrt{gH}$ (the CFL
    limit of Ch. 10). This scheme is OUR choice — the book has no numerical model.", code="dx, H = 1.0e4, 1.33
    # 10 km cells; a layer 1.33 m deep\nc = np.sqrt(G0*H)                            # 3.61 m/s\nprint(c, 0.5*dx/c)                             # time step at
    Courant number 0.5: 1385 s")` (**P335**) + `> 🔧 **Our choice** — forward–backward time stepping on a staggered grid in one
    dimension here; an energy-conserving C-grid scheme with third-order Runge–Kutta steps for the two-dimensional model
    `SW.ShallowWater`, used only to make the cached runs of C11 and C13.`
11. `nb.code` — `out = SW.linear_1d_run(eta0, dx=1.0e4, dt=1385.0, n_steps=144, H=1.3286, f=0.0)` with `eta0 =
    np.exp(-(x/1.0e5)**2)` on 400 cells (−2000…2000 km); `GFD.long_wave_speed(1.3286)`; `GFD.long_wave_speed(4200.0)`.
    *expect* (designer's own scratch march, same scheme): after 2 × 10⁵ s the bump has split into two pulses of height 0.500
    at x = ±715 km (speed 3.58 m/s on the 10 km grid against c = 3.610 m/s); √(gH) = 202.95 m/s for 4200 m. *explain:* 4
    numbered points.
12. `nb.check_agree` — **from scratch** (the 25-line march): `for n in range(n_steps): eta -= dt*H*np.diff(u)/dx; u[1:-1] +=
    dt*(f*v_face - G0*np.diff(eta)/dx); v -= dt*f*u_centre`, every line commented; `assert np.allclose(eta_mine,
    out["eta"][-1])` (same scheme → same numbers to round-off); then with `f=f35` both again.
13. `nb.figure` — x–t diagrams (Hovmöller) of η: left without rotation (two straight rays of slope ±c, nothing left behind),
    right with f = f35 (rays that disperse into a wave train, and a part of the bump that **stays**). Dashed lines x = ±ct.
    *see / read* ("slope of a ray = speed; whatever sits on the line x = 0 at late time never left") */ change* ("…H is
    quadrupled: the rays are twice as steep"). Forward pointer: "why something stays is C12".
14. `nb.note` — **N53 [B]** equivalent depth: "A stratified mode with long-wave speed c behaves like a homogeneous layer of
    depth $H_e$ defined by $c^2=gH_e$ (13.46)." `GFD.equivalent_depth(3.6096)` *expect:* 1.329 m. "Made precise in C08."
15. `nb.code` (tag `slow`, under `if not FAST:`) — the two-dimensional model's API in six lines: `m = SW.ShallowWater(64, 64,
    2.0e6, 2.0e6, 1.3286, f35, bc="closed")`; `st = {"eta": SW.gaussian_bump(m, 0.05, 1.0e6, 1.0e6, 1.5e5), "u": …, "v": …}`;
    `hist = SW.run(m, st, t_end=20*SW.dt_limit(m))`; volume `hist["eta"].sum(axis=(1, 2))` constant to round-off
    (`SW.continuity_tendency` sums to zero). (Everything else two-dimensional comes from the caches of Part D.)
16. `nb.md` — **What would change if…** "…the fluid is not one layer but continuously stratified? Surprisingly little: it
    splits into a stack of such layers, each with its own depth (C08)."

---

### A.9 §13.9 Normal Modes in a Continuously Stratified Layer — C08
#### C08 — Vertical normal modes (13.56) and the equivalent depth (13.62)
1. `nb.section("13.9", "Normal Modes in a Continuously Stratified Layer", intro="**What is this section about?** How a
   continuously stratified ocean or atmosphere can be replaced by a handful of shallow-water systems: the vertical
   structure separates out as an eigenvalue problem set by N(z). ⚠️ Here z = 0 is the surface and the flat bottom is at
   z = −H.")`
2. `nb.core("C08", "Vertical normal modes: $\\dfrac{d}{dz}\\Big(\\dfrac{1}{N^2}\\dfrac{d\\psi_n}{dz}\\Big)+\\dfrac{1}{c_n^2}\\psi_n=0$ (13.56), each a shallow-water system with
   $c_n^2\\equiv gH_e$ (13.62)", question="How can a one-layer model say anything about an ocean whose density changes all the way
   down?")`
3. `nb.md` — **Plain words:** "When El Niño begins, a bulge in the thermocline crosses the Pacific in about two months —
   far slower than a tsunami (hours) yet it is the same kind of wave. The ocean has many ways to wobble: all together, top
   to bottom (fast), or with the upper ocean moving against the deep ocean (slow), or in three alternating layers (slower
   still). Each of these *modes* has a fixed vertical shape and moves horizontally exactly like a single shallow layer —
   of a different, much smaller, depth. For the first internal mode that depth is about a metre."
4. `nb.md` — **The idea** (ASCII, three columns of arrows): mode 0 `→ → → →` (barotropic) · mode 1 `→ → | ← ←` · mode 2
   `→ | ← | →`; "mode n has n levels where the horizontal velocity changes sign". Table: mode · speed · equivalent depth.
5. **(re-enter C08 after each recap of this row)** `nb.recap("R19", "Continuity", "$\\dfrac{\\partial u}{\\partial x}+\\dfrac{\\partial v}{\\partial y}+\\dfrac{\\partial w}{\\partial z}=0$ (13.47).", where="Ch. 4")` · `nb.note` **N55 [B]** the
   horizontal momentum equations $\dfrac{\partial u}{\partial t}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$ and $\dfrac{\partial v}{\partial t}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$ (13.48)–(13.49) ·
   `nb.recap("R20", "Hydrostatic balance and the density equation", "$0=-\\dfrac{\\partial p}{\\partial z}-g\\rho$ (13.50) and $\\dfrac{\\partial\\rho}{\\partial t}-\\dfrac{\\rho_0N^2}{g}w=0$ (13.51):
   density at a point changes only because the background gradient is carried up or down, $w\\,d\\bar\\rho/dz=-\\rho_0N^2w/g$.",
   where="Ch. 7, $\\frac{d\\bar\\rho}{dz}=-\\frac{\\rho_0N^2}{g}\\Rightarrow\\frac{\\partial\\rho'}{\\partial t}-\\frac{N^2\\rho_0}gw=0$ (7.131), there without rotation")` +
   `ch13.hydrostatic_linear_set_sympy()["ok"]` → True. "These five are the Start of D11."
6. `nb.primer("Sturm–Liouville problem: a ladder of eigenvalues, eigenfunctions with n zero crossings, orthogonality", "A
   string fixed at both ends can only vibrate in certain shapes — one hump, two, three… — each with its own frequency.
   Any equation of the form (p ψ′)′ + λ w ψ = 0 with conditions at two ends behaves the same way: only special values λ₀
   < λ₁ < … allow a solution, the n-th solution crosses zero n times, and two different solutions are orthogonal (their
   weighted product integrates to zero). Here λ = 1/c² and p = 1/N².", code="z = np.linspace(-1, 0, 2001)
   # a layer of unit depth\np1, p2 = np.cos(np.pi*z), np.cos(2*np.pi*z)        # the first two shapes for uniform N, rigid lid\nprint(round(np.trapezoid(p1*p2, z),
   12), np.trapezoid(p1*p1, z))   # 0.0 and 0.5: orthogonal, not normalised")` (**P314**) — one-line reminders of P167
   (separation of variables), P257 (eigenvalue problem for a differential operator), P218a (integration by parts).
7. `nb.derivation("D11", …)` — Part F D11 (13 steps), ref "13.56": **N56 [B]** step 1,
   $[u,v,p/\rho_0]=\sum_{n=0}^\infty[u_n,v_n,p_n]\,\psi_n(z)$ (13.52); **N57 [B]** step 2, $w=\sum_{n=0}^\infty w_n\int_{-H}^z\psi_n\,dz$ (13.53); **N58 [B]** step 3,
   $\rho=\sum_{n=0}^\infty\rho_n\dfrac{d\psi_n}{dz}$ (13.54); **N59 [B]** step 5, the separation constant $\dfrac{d\psi_n/dz}{N^2\int_{-H}^z\psi_n\,dz}=\dfrac{\rho_0}{g}\dfrac{w_n}{\partial\rho_n/\partial t}\equiv-\dfrac1{c_n^2}$ (13.55); **N61 [B]**
   step 8, $\dfrac{\partial u_n}{\partial x}+\dfrac{\partial v_n}{\partial y}+\dfrac1{c_n^2}\dfrac{\partial p_n}{\partial t}=0$ (13.57); **N62 [B]** step 9, $\dfrac{\partial u_n}{\partial t}-fv_n=-\dfrac{\partial p_n}{\partial x}$ and
   $\dfrac{\partial v_n}{\partial t}+fu_n=-\dfrac{\partial p_n}{\partial y}$ (13.58)–(13.59); **N63 [B]** steps 3 and 7, $p_n=-\dfrac g{\rho_0}\rho_n$ and $w_n=\dfrac1{c_n^2}\dfrac{\partial p_n}{\partial t}$ (13.60)–(13.61);
   **N64 [B]** step 10, the identification $p_n\leftrightarrow g\eta$, $c_n^2\leftrightarrow gH$ and the definition $c_n^2\equiv gH_e$ (13.62); **N60 [B]** steps 11–13,
   orthogonality with weight 1, $\int_{-H}^0\psi_m\psi_n\,dz=0$ for m ≠ n, **for the rigid lid and for the free surface alike** (the
   book only states it; the surface term $\psi_m(0)\psi_n(0)/g$ belongs to the second, "energy" relation, not to this one). `> ⚠️ Common confusion (trap
   T12):` "The amplitudes do not share units: ψ_n is a pure number, so u_n and v_n are velocities [m/s] and p_n is p/ρ₀
   [m²/s²]; but w_n multiplies an integral of ψ_n over z [m], so w_n is in s⁻¹, and ρ_n multiplies dψ_n/dz [1/m], so ρ_n is in
   kg/m²." `VM.modal_amplitudes(c_n, p_n_t, p_n=p_n)` (without `p_n` the entry "rho_n" is `None`).
8. `nb.note` — **N65 [B]** the link the boundary conditions need: $w=\dfrac{g(\partial\rho/\partial t)}{\rho_0N^2}=-\dfrac1{\rho_0N^2}\dfrac{\partial^2p}{\partial z\,\partial t}=-\dfrac1{N^2}\sum_{n=0}^\infty\dfrac{\partial p_n}{\partial t}\dfrac{d\psi_n}{dz}$ (13.63). **N66
   [B]** flat bottom, w = 0: $\dfrac{d\psi_n}{dz}=0$ at $z=-H$ (13.64). **(re-enter C08)** `nb.recap("R21", "The linearised free surface", "$w=\\partial\\eta/\\partial t$ and
   $p=\\rho_0g\\eta$ at z = 0 (the book labels this pair (13.65′)), which combine to $\\partial p/\\partial t=\\rho_0gw$ at z = 0.", where="Ch. 7: the
   kinematic condition $\\big(\\frac{\\partial\\phi}{\\partial z}\\big)_{z=0}\\cong\\frac{\\partial\\eta}{\\partial t}$ (7.18) and the dynamic condition $\\big(\\frac{\\partial\\phi}{\\partial t}\\big)_{z=0}\\cong-g\\eta$ (7.21)")`
   **N67 [B]** hence the surface condition on the modes, $\dfrac{d\psi_n}{dz}+\dfrac{N^2}{g}\psi_n=0$ at $z=0$ (13.65).
9. `nb.primer("Robin (mixed) boundary condition, between Dirichlet and Neumann", "A Dirichlet condition fixes the value
   (ψ = 0), a Neumann condition fixes the slope (ψ′ = 0), a Robin condition ties them together (ψ′ + aψ = 0). With a = 0 it
   is Neumann. Here a = N²/g is tiny — about 10⁻⁶ per metre in the ocean — so for internal modes the free surface is
   almost a rigid lid.", code="N, g = 2.7e-3, 9.80665\nprint(N**2/g, N**2/g*4200)        # 7.4e-07 per metre; 0.0031 over the
   whole depth")` (**P315**)
10. `nb.note` — the uniform-N worked case: **N68 [C]** $\dfrac{d^2\psi_n}{dz^2}+\dfrac{N^2}{c_n^2}\psi_n=0$ (13.66); **N69 [C]** its solution
    $\psi_n=A_n\cos\dfrac{Nz}{c_n}+B_n\sin\dfrac{Nz}{c_n}$ (13.67); **N70 [C]** the surface condition gives $B_n=-\dfrac{c_nN}{g}A_n$ (13.68); **N71 [B]** the bottom
    condition then gives the eigenvalue condition $\tan\dfrac{NH}{c_n}=\dfrac{c_nN}{g}$ (13.69). The result is stated, not derived (the book
    writes it out).
11. `nb.primer("graphical roots of a transcendental equation such as tan x = εx", "Some equations have the unknown both
    inside and outside a function and cannot be solved by algebra. Plot both sides against x: the roots are where the
    curves cross. Reading the picture also tells you roughly where they are — here, just above 0, π, 2π, … — which is all
    a root-finder such as `brentq` needs.", code="eps = 0.0031                                   # N^2 H / g for our ocean\nF =
    lambda X: np.tan(X) - eps/X                # zero where tan X = eps/X\nprint(brentq(F, 1e-6, 1.0), brentq(F, np.pi + 1e-9, np.pi + 1.0) -
    np.pi)   # 0.0558, 0.00099")` (**P316**; reminder of P108 `brentq`)
12. `nb.figure` — **N72 [B]** our Fig. 13.12 (`ch13.fig_mode_roots()`): tan X and (N²H/g)/X against X = NH/c_n from 0 to
    2.5π, roots marked. *see / read* ("the hyperbola is so low that it meets each tan branch almost at its foot: X ≈ nπ")
    */ change* ("…g were a thousand times smaller: the roots would climb the branches and the lid approximation would fail").
13. **(re-enter C08)** `nb.recap("R22", "The barotropic mode", "The first root has $NH/c_0\\ll1$; there $\\tan X\\approx X$ gives $c_0=\\sqrt{gH}$ (13.70), the
    long-wave speed of a homogeneous ocean, with a nearly uniform structure $\\psi_0\\simeq1-N^2z/g\\simeq1$.", where="Ch. 7, long
    waves")` + > ⚠️ **slip #6 — the book prints** that this root occurs "for NH/c_n = 1"**; the correct form is** $NH/c_0\ll1$
    (ours, computed: 0.0558 for our ocean).
14. `nb.note` — **N73 [B]** the baroclinic roots: since $c_nN/g\ll1$, $\tan\dfrac{NH}{c_n}=0$, so $c_n=\dfrac{NH}{n\pi}$, $n=1,2,3,\dots$ (13.71).
    `GFD.baroclinic_mode_speed(2.7e-3, 4200.0, n)`. **N74 [B]** numbers with our inputs (H = 4200 m, N = 2.7 × 10⁻³ s⁻¹): *expect*
    c₁ = 3.610 m/s, H_e = 1.329 m, against c₀ = 202.95 m/s — "a factor 56 in speed, 3200 in equivalent depth".
15. `nb.worked_example("the first baroclinic mode by hand", "Take N = π × 10⁻³ ≈ 3.14 × 10⁻³ s⁻¹ and H = 4000 m. 1. $c_1=NH/\\pi$ =
    (π × 10⁻³ × 4000)/π = 4 m/s. 2. Equivalent depth $H_e=c_1^2/g$ ≈ 16/10 = 1.6 m. 3. Mode 2 travels at half that speed, c₂ = c₁/2; mode 3 at a third, c₁/3 ≈ 1.33 m/s. 4. Barotropic: √(gH) = √40 000 = 200 m/s. 5. Size of the surface term: N²H/g = 10⁻⁵ ×
    4000/10 = 0.004 ≪ 1 — the surface is almost a lid for these modes.")`
16. `nb.primer("scipy.linalg.eigh_tridiagonal for a discretised eigenproblem", "Replace ψ by its values on a grid and the
    second derivative by differences: the differential eigenproblem becomes a matrix one, with only the diagonal and its
    two neighbours filled. `eigh_tridiagonal(d, e)` returns all eigenvalues (sorted) and eigenvectors of a symmetric
    tridiagonal matrix from its diagonal d and off-diagonal e.", code="n = 200; h = 1.0/n                            #
    interior points of a unit string\nlam, vec = eigh_tridiagonal(2*np.ones(n-1)/h**2, -np.ones(n-2)/h**2)\nprint(np.sqrt(lam[:3])/np.pi)
    # 1, 2, 3 (to 4 digits): the modes of a string")` (**P317**; reminders of P80 eigenvalues, P247 generalised symmetric
    eigenproblem)
17. `nb.code` — `z = np.linspace(-4200.0, 0.0, 401)`; `m = VM.vertical_modes(z, np.full_like(z, 2.7e-3**2), n_modes=4,
    lid="free")`; `mu = VM.modes_uniform_N(2.7e-3, 4200.0, n_modes=4)`; `VM.rigid_lid_error(2.7e-3, 4200.0, n=1)`;
    `GFD.rossby_radius(m.c[1], f35)`. *expect:* numerical c = (203.05, 3.6085, 1.8047, 1.2032) m/s; exact (`brentq`) the same
    to 5 digits; exact barotropic speed 0.052 % above √(gH); rigid-lid error of c₁ = +3.2 × 10⁻⁴ (free-surface c₁ = 3.6085
    against NH/π = 3.6096), 7.9 × 10⁻⁵ for n = 2, 3.5 × 10⁻⁵ for n = 3; first root X₀ = 0.05585; Λ₁ = 43.1 km at 35° N. *explain:*
    5 numbered points.
18. `nb.check_agree` — **from scratch:** the mode equation for uniform N with a rigid lid as a small symmetric matrix
    (`A = (2I − shifted)/h²` with the two Neumann end rows halved), `np.linalg.eigh`, `c_mine = N/np.sqrt(lam[1])`; `assert
    np.isclose(c_mine, GFD.baroclinic_mode_speed(2.7e-3, 4200.0), rtol=1e-4)` and `assert np.isclose(c_mine,
    VM.vertical_modes(z, N2, lid="rigid").c[0], rtol=1e-6)` — **index 0**: a rigid-lid `Modes` has no barotropic entry, so its
    first element is the first baroclinic mode (free-surface modes keep it at index 1).
19. `nb.note` — **N75 [B]** the rigid-lid approximation: "Replace the surface condition by w = 0 at z = 0, that is, the slope of ψ_n vanishes there.
    Then $\psi_n=\cos\dfrac{n\pi z}{H}$, n = 0, 1, 2, …: the baroclinic speeds change by 3 parts in 10⁴ or less (cell above) and the
    barotropic mode disappears (c₀ → ∞). What the lid does **not** mean: the surface pressure still varies under it — the lid
    pushes back." *(Reported difference from the curation's wording "a few parts in 10⁵": that holds from n = 3 on; for n = 1
    the computed value is 3.2 × 10⁻⁴.)*
20. `nb.primer("projecting a profile on modes: the inner product of two functions", "Two vectors are orthogonal when
    their dot product is zero; two functions when the integral of their product is zero. Because the modes are
    orthogonal, the amount of mode n in any profile F(z) is found by one integral, a_n = ∫Fψ_n dz / ∫ψ_n² dz — like
    reading off a component along an axis.", code="z = np.linspace(-1, 0, 2001); F = 1 + z                       # a linear
    profile\np1 = np.cos(np.pi*z)\nprint(np.trapezoid(F*p1, z)/np.trapezoid(p1*p1, z))            # 0.405 = 4/pi^2: its mode-1
    content")` (**P318**) + `nb.code`: `VM.orthogonality_matrix(m)` for the free-surface modes `m` (default `kind="psi"`:
    weight 1, no surface term — identity to 10⁻⁸ off the diagonal, *expect* confirmed on the implementer's code and by the
    designer's own finite-volume solve), `VM.orthogonality_matrix(m, kind="energy", N2=N2)` (the second relation, which
    carries the surface term), `a = VM.project(m, profile)`, `VM.reconstruct(m, a)`, `VM.w_structure(m, 1)`,
    `VM.rho_structure(m, 1)`; then the same projection with the rigid-lid modes `mr = VM.vertical_modes(z, N2, lid="rigid")`:
    "`VM.project(mr, np.ones_like(z))` is zero — rigid-lid modes cannot hold a depth mean; add `profile.mean()` back".
21. `nb.figure` — **N76 [B]** our Fig. 13.13 (`ch13.fig_vertical_modes()`): left, the first three modes for uniform N
    (cosines); right, the same for the thermocline profile `ch13.thermocline_N2(z)` of §13.2, with N(z) in blue behind.
    *expect* (thermocline profile, defaults N_deep = 5 × 10⁻⁴, N_peak = 8 × 10⁻³ s⁻¹, z_t = −300 m, width 150 m; designer's
    finite-volume scratch solve, 801 nodes): c₁ = 1.772, c₂ = 0.618, c₃ = 0.466 m/s; H_e = 0.320 m for mode 1; Λ₁ = 21.2 km at
    35° N; mode 1 changes sign at z = −446 m; `VM.wkb_mode_speed` gives 1.30 m/s for n = 1 — **27 % low**, because a sharp
    thermocline is not a slowly varying medium (a reported surprise; the estimate is good to a few per cent only for
    smooth profiles). *see / read* ("mode n crosses zero n times; with a thermocline the crossings crowd into it") */
    change*.
22. `nb.plotly` — **F6** (`slider_figure`, slider = thermocline depth −z_t from 100 to 800 m, 15 steps): N(z) and the first
    three modes; trace names carry c_n and H_e. *expect:* c₁ = 1.245, 1.772, 2.474 m/s for z_t = −150, −300, −600 m. *read:*
    "a deeper thermocline is a faster first mode — the El Niño thermocline signal in one slider".
23. `nb.note` — **N77 [B]** the fine print: flat bottom; no sheared mean current; hydrostatic, so ω ≪ N and the shapes do
    not depend on frequency; valid with or without f and β.
24. `nb.explainer("vertical_modes", heading="How can one layer stand for a stratified ocean?", why="Reshaping N(z) and
    watching the shapes, speeds and equivalent depths respond shows what the eigenproblem does; stepping n shows the zero
    crossings appear one by one; the lid toggle removes the barotropic mode and barely moves the others.", tries=["Step n
    from 0 to 3 and count the zero crossings.", "Deepen the thermocline and watch c₁ and the Rossby radius grow.", "Toggle
    the rigid lid: the top rung of the ladder disappears, the others do not move.", "Click a depth to see the two terms
    of the mode equation cancel there."])`
25. `nb.md` — **What would change if…** "…we now ask what waves one of these shallow-water systems supports when the earth
    turns? One cubic equation answers for all of them (C09)."

---

### A.10 §13.10 High- and Low-Frequency Regimes in Shallow-Water Equations — C09
#### C09 — The complete dispersion relation (13.76)
1. `nb.section("13.10", "High- and Low-Frequency Regimes in Shallow-Water Equations", intro="**What is this section about?**
   One equation that contains every wave of a rotating shallow layer, and how to tell, from the frequency alone, which
   terms matter.")`
2. `nb.core("C09", "The complete dispersion relation of rotating shallow water: $\\omega^3-c^2\\omega K^2-f_0^2\\omega-c^2\\beta k=0$, with
   $K^2=k^2+l^2$, $c=\\sqrt{gH}$ (13.76)", question="Poincaré, Kelvin, Rossby — how many different kinds of wave does one layer of
   water really have?")`
3. `nb.md` — **Plain words:** "A weather forecast must not let fast gravity waves swamp the slow weather; a tide model
   wants exactly those fast waves. Both run the same equations. The waves separate by frequency: two fast ones (period
   hours) that are gravity waves bent by rotation, and one slow one (period days to years) that exists only because the
   Coriolis parameter changes with latitude. Knowing which is which tells a modeller which terms may be dropped."
4. `nb.md` — **The idea:** a table of the four terms of the cubic — ω³ (pure time change) · −c²K²ω (gravity) · −f₀²ω
   (rotation) · −c²βk (variation of f) — with "who balances whom" for fast and slow waves.
5. `nb.md` — one-line reminders: P177 (operator elimination for linear PDEs, ch07), P176 (complex amplitudes, plane waves),
   P198 (dominant balance), P256 (`np.roots`).
6. `nb.derivation("D12", …)` — Part F D12 (10 steps), ref "13.75": **N78 [C]** its steps 2–3,
   $\dfrac{\partial^2u}{\partial t^2}-f\dfrac{\partial v}{\partial t}=gH\dfrac{\partial}{\partial x}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ and $\dfrac{\partial^2v}{\partial t^2}+f\dfrac{\partial u}{\partial t}=gH\dfrac{\partial}{\partial y}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ (13.72)–(13.73); **N79 [C]** step 4,
   $\dfrac{\partial^3v}{\partial t^3}+f\Big[f\dfrac{\partial v}{\partial t}+gH\dfrac{\partial}{\partial x}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)\Big]=gH\dfrac{\partial^2}{\partial y\,\partial t}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ (13.74); **N80 [B]** step 6, the linear vorticity equation on the β-plane,
   $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial u}{\partial y}-\dfrac{\partial v}{\partial x}\Big)-f_0\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)-\beta v=0$ ("the first appearance of the idea C13 turns into potential vorticity");
   **N81 [B]** the result, $\dfrac{\partial^3v}{\partial t^3}-gH\dfrac{\partial}{\partial t}\nabla_H^2v+f_0^2\dfrac{\partial v}{\partial t}-gH\beta\dfrac{\partial v}{\partial x}=0$ (13.75). `> ⚠️ Common confusion (trap T10):` "f is
   replaced by the constant f₀ everywhere except where it is differentiated (that derivative is β)." +
   `ch13.v_equation_sympy()["ok"]` → True.
7. `nb.primer("three real roots of a cubic: the discriminant and the trigonometric form", "A cubic ω³ + pω + q = 0 with
   real p, q has three real roots exactly when its discriminant −4p³ − 27q² is positive (this needs p < 0). They are then
   given without complex arithmetic by $\\omega_j=2\\sqrt{-p/3}\\,\\cos\\big[\\tfrac13\\arccos\\big(\\tfrac{3q}{2p}\\sqrt{-3/p}\\big)-\\tfrac{2\\pi j}3\\big]$, j = 0, 1, 2. When one root is thousands of times
   smaller than the others this form keeps its digits; a general polynomial solver may not.", code="p, q = -7.0, 6.0
   # w^3 - 7w + 6 = (w - 1)(w - 2)(w + 3)\nprint(-4*p**3 - 27*q**2)                    # 400 > 0: three real roots\nm = 2*np.sqrt(-p/3); th =
   np.arccos(3*q/(p*m))/3\nprint(np.sort(m*np.cos(th - 2*np.pi*np.arange(3)/3)))   # [-3.  1.  2.]")` (**P319**; reminder of P262
   cube roots)
8. `nb.derivation("D13", …)` — Part F D13 (7 steps), ref "13.76": **N82 [B]** its steps 4–5: the discriminant
   $4(c^2K^2+f_0^2)^3-27(c^2\beta k)^2$ is positive, so all three roots are real (the book says only "it can be shown"). `> ⚠️
   Common confusion (trap T9):` "'Two superinertial, one subinertial' is true under β-plane scaling; the two fast roots
   have opposite signs and only satisfy abs(ω) > abs(f₀), not ω ≫ f." **N83 [B]** steps 6–7, the three regimes with the two
   scale ratios $\dfrac{c^2\beta k}{c^2\omega K^2}\sim\dfrac\beta{\omega K}$ and $\dfrac{\omega^3}{c^2\beta k}\ll1$. > ⚠️ **slip #7 — the book prints** ω ≫ f as the range in which the
   first term of $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76) is negligible**; the correct form is** ω ≪ f.
9. `nb.primer("reading a dispersion diagram with several branches: signed ω for signed k, a logarithmic frequency axis",
   "A dispersion diagram plots frequency against wavenumber; each curve ('branch') is one kind of wave. The slope of the
   line from the origin to a point is the phase speed, the slope of the curve itself the group velocity. We let k carry
   the direction (k > 0 eastward) and keep ω ≥ 0 for the plot; when branches differ by factors of a thousand, the frequency
   axis must be logarithmic or the slow branch hides on the axis.", code="k = np.array([1e-7, 1e-6, 1e-5])\nprint(np.sqrt(1e-8 +
   4e4*k**2))               # a fast branch: flattens to f = 1e-4 as k -> 0\nprint(2e-11*k/(k**2 + 2.5e-13))               # a slow
   branch: rises, peaks, falls")` (**P320**)
10. `nb.worked_example("fast and slow root by hand", "f₀ = 10⁻⁴ s⁻¹, c = 200 m/s, β = 2 × 10⁻¹¹ m⁻¹ s⁻¹, k = 10⁻⁶ m⁻¹ (wavelength 6300
    km), l = 0. 1. Fast: drop the β term → ω² = f₀² + c²k² = 10⁻⁸ + 4 × 10⁻⁸ → ω = ±2.24 × 10⁻⁴ s⁻¹ (period 7.8 h). 2. Slow: drop
    ω³ → ω = −c²βk/(c²k² + f₀²) = −(4 × 10⁴ × 2 × 10⁻¹¹ × 10⁻⁶)/(5 × 10⁻⁸) = −1.6 × 10⁻⁵ s⁻¹ (period 4.5 days; minus = westward for
    k > 0). 3. Was dropping ω³ fair? (1.6 × 10⁻⁵)³ = 4 × 10⁻¹⁵ against c²βk = 8 × 10⁻¹³: half a per cent. 4. Sum of the three roots
    must be 0 (no ω² term): +2.24 − 2.24 + the small one — so the fast pair is not exactly symmetric; the slow root is the
    difference.")`
11. `nb.code` — `k = 2*np.pi/3.1e6`; `w = GFD.shallow_water_omega(k, 0.0, 202.95, f35, beta35)`;
    `GFD.dispersion_term_sizes(k, 0.0, 202.95, f35, beta35, w[1])` and `(…, w[2])`; `GFD.shallow_water_discriminant(...)`;
    `GFD.shallow_water_regime(w[1], f35)`; `GFD.shallow_water_branches(...)`; then the low-latitude case `k20 = 2*np.pi/2.0e7`,
    `GFD.shallow_water_discriminant(k20, 0.0, 202.95, f12, beta12)` (negative) and `GFD.shallow_water_omega(k20, 0.0, 202.95, f12,
    beta12)` → `(nan, nan, nan)`: "test the discriminant first; NaN is the function's way of saying 'no three real roots'".
    *expect* (35° N, c = 202.95 m/s, wavelength 3100
    km): roots (−4.1525 × 10⁻⁴, −8.888 × 10⁻⁶, +4.2414 × 10⁻⁴) s⁻¹; periods 4.12 h (fast) and 8.18 days (slow); slow/fast = 0.021;
    discriminant > 0; for the first baroclinic mode (c = 3.610 m/s): (−8.3936 × 10⁻⁵, −7.023 × 10⁻⁸, +8.4006 × 10⁻⁵) s⁻¹ — 20.8 h
    and 2.8 years; the slow root differs from the Rossby formula of C15 by 4.5 × 10⁻⁴ (external) and 7 × 10⁻⁷ (baroclinic),
    relative. *explain:* 5 numbered points.
12. `nb.check_agree` — **from scratch:** `np.sort(np.roots([1.0, 0.0, -(c**2*k**2 + f35**2), -c**2*beta35*k]).real)`; `assert
    np.allclose(roots_mine, w, rtol=1e-9)`; and the sum of the roots divided by the largest is below 10⁻¹².
13. `nb.figure` — four bars per frequency for the three cases ω = 3f, ω = f (taken on the fast branch) and ω = 0.1f (slow
    branch): `GFD.dispersion_term_sizes`, bars ω³ purple, −c²K²ω orange, −f₀²ω teal, −c²βk amber, log axis. *see:* at 3f the
    ω³ bar is the tallest and the β bar invisible; at 0.1f the ω³ bar is the smallest. *read:* "each regime of the book is
    'drop the shortest bar'". *change:* "…β = 0: the amber bar vanishes and the slow root becomes ω = 0 — a steady
    geostrophic flow, the state C12's adjustment ends in". + `nb.md` **An honest caveat (ours, computed):** "all three roots
    are real for every wavelength exactly when $\beta c<f_0^2$, i.e. when the Rossby radius c/f₀ is smaller than R tan θ₀. For
    the external mode of a 4200 m ocean that fails equatorward of 26.3°; there the discriminant is negative for a band of
    very long waves and the cubic has no three real roots — an artefact of freezing f at f₀, not a real instability.
    `GFD.shallow_water_discriminant` is the flag; `GFD.shallow_water_omega` returns NaN there and the figure leaves the band
    blank with a label." *expect:* βc/f₀² = 0.544 at 35° (external), 0.0097 (baroclinic), 4.94 at 12° (external); at 12° N
    with the external speed the discriminant is negative for all wavelengths above 12 500 km (verified: 20 000 km gives
    NaN, and `np.roots` a complex pair there).
14. `nb.plotly` — **F7** (`slider_figure`, slider = latitude from 5° to 75°, 15 steps; dropdown external / first baroclinic
    c): abs(ω) of the three branches of `GFD.shallow_water_branches` and the Kelvin line `GFD.kelvin_omega(k, c)` against k
    on log–log axes, ω = f dashed. *see / read* ("the gap between the fast and slow branches is the reason filtered
    'quasi-geostrophic' models work") */ change*. (The Kelvin line is explained in C11.)
15. `nb.md` — **What would change if…** "…we look at the fast roots alone, on an f-plane? They are the gravity waves of Ch. 7
    with a floor under their frequency (C10)."

---

### A.11 §13.11 Gravity Waves with Rotation — C10
#### C10 — Poincaré waves (13.82)
1. `nb.section("13.11", "Gravity Waves with Rotation", intro="**What is this section about?** What rotation does to a long
   gravity wave: it cannot oscillate more slowly than f, it becomes dispersive, and the water under it moves in ellipses.")`
2. `nb.core("C10", "Poincaré waves: $\\omega^2=f^2+gHK^2$, $K=\\sqrt{k^2+l^2}$ (13.82)", question="What does the earth's rotation do to a
   long gravity wave?")`
3. `nb.md` — **Plain words:** "After a storm passes, current meters in the open ocean show the water going round in
   circles a few kilometres across, once every 17 to 20 hours at mid-latitudes, for days. Nothing is pushing it; it is
   coasting, and the Coriolis force keeps bending its path. That circle is the longest possible gravity wave. Shorter
   waves — the tides among them — are mixtures: partly a gravity wave sloshing back and forth, partly this turning."
4. `nb.md` — **The idea:** two restoring forces, two frequencies: gravity alone gives ω = cK; rotation alone gives ω = f;
   together $\omega^2=f^2+c^2K^2$ — "add the squares, like the sides of a right triangle".
5. `nb.note` — **N84 [B]** plane-wave amplitudes on the f-plane: $-i\omega\hat u-f\hat v=-ikg\hat\eta$, $-i\omega\hat v+f\hat u=-ilg\hat\eta$,
   $-i\omega\hat\eta+iH(k\hat u+l\hat v)=0$ (13.77)–(13.79). "These are the Start of D14."
6. `nb.primer("polarisation relations: solving a 2 × 2 complex system for the velocity amplitudes", "For a plane wave
   every field is a complex amplitude times the same $e^{i(kx+ly-\\omega t)}$. The momentum equations become two linear equations
   for the two velocity amplitudes in terms of the height amplitude. Solving them (Cramer's rule) tells how the current
   is oriented and timed relative to the surface — the wave's 'polarisation'. A factor i means a quarter-period shift.",
   code="w, f, k, g = 2.0, 1.0, 1.0, 1.0         # easy numbers, l = 0\nM = np.array([[-1j*w, -f], [f, -1j*w]])       # unknowns
   (u_hat, v_hat)\nprint(np.linalg.solve(M, [-1j*k*g, 0.0]))   # [0.667, -0.333j]: v lags u by a quarter period")` (**P321**)
7. `nb.derivation("D14", …)` — Part F D14 (8 steps), ref "13.82": **N85 [B]** its step 3, $\hat u=\dfrac{g\hat\eta}{\omega^2-f^2}(\omega k+ifl)$ and
   $\hat v=\dfrac{g\hat\eta}{\omega^2-f^2}(-ifk+\omega l)$ (13.80); **N86 [C]** step 5, $\omega^2-f^2=gH(k^2+l^2)$ (13.81), the same relation before K is introduced;
   steps 7–8 (ours): the group velocity $\mathbf c_g=c^2\mathbf K/\omega$ and $c_pc_g=c^2$.
8. `nb.worked_example("a long wave at mid-latitude", "f = 10⁻⁴ s⁻¹, c = 200 m/s, wavelength 6300 km (K = 10⁻⁶ m⁻¹). 1. Frequency: ω² = 10⁻⁸ +
   (200 × 10⁻⁶)² = 10⁻⁸ + 4 × 10⁻⁸ → ω = 2.24 × 10⁻⁴ s⁻¹ = 2.24 f. 2. Phase speed ω/K = 224 m/s — faster than c. 3. Group velocity
   c²K/ω = 4 × 10⁴ × 10⁻⁶/2.24 × 10⁻⁴ = 179 m/s — slower than c; product 224 × 179 = 4.0 × 10⁴ = c². 4. Current ellipse: axis ratio
   ω/f = 2.24, turning clockwise (northern hemisphere, f > 0; counter-clockwise for f < 0). 5. Ten times longer wave: ω =
   1.02 f — almost a pure inertial circle.")`
9. `nb.code` — `K = 2*np.pi/3.1e6`; `w = GFD.poincare_omega(K, f35, 202.95)`; `GFD.poincare_group_velocity(K, 0.0, f35, 202.95)`;
   `uh, vh = GFD.poincare_amplitudes(K, 0.0, w, f35, G0, 0.5)`; `orb = GFD.poincare_orbit(t, K, 0.5, 4200.0, f35)`. *expect:*
   ω = 4.198 × 10⁻⁴ s⁻¹ = 5.018 f; c_p = 207.10 m/s, c_g = 198.88 m/s, product = c²; velocity amplitudes for η̂ = 0.5 m: u 0.0247
   m/s, v 0.0049 m/s (ratio ω/f); `orb["sense"]` = "clockwise"; for `-f35` "counter-clockwise". *explain:* 4 points.
10. `nb.check_agree` — **from scratch:** `w_mine = np.sqrt(f35**2 + G0*4200.0*K**2)`; `u_mine = w_mine*0.5/(K*4200.0)*np.cos(-w_mine*t)`,
    `v_mine = f35*0.5/(K*4200.0)*np.sin(-w_mine*t)`; `assert np.allclose` against `orb["u"]`, `orb["v"]`.
11. `nb.figure` — **N87 [B]** our Fig. 13.14 (`ch13.fig_poincare_kelvin_dispersion()`): ω against K, the hyperbola from ω = f
    and its asymptote ω = cK; the Kelvin line is drawn dashed "(C11)"; inset (**ours**, labelled): the same on a
    logarithmic axis with the slow Rossby branch three decades below. *see / read / change* ("…f = 0: the hyperbola
    collapses onto its asymptote — Chapter 7's long waves").
12. `nb.note` — **N88 [B]** the real fields for a wave along x (l = 0), with $\eta=\hat\eta\cos(kx-\omega t)$: $u=\dfrac{\omega\hat\eta}{kH}\cos(kx-\omega t)$ and
    $v=\dfrac{f\hat\eta}{kH}\sin(kx-\omega t)$ (13.83) — stated (real parts of D14's amplitudes, using $\omega^2-f^2=gHk^2$). `GFD.poincare_fields(x, y,
    t, K, 0.0, 0.5, 4200.0, f35)`.
13. `nb.primer("sense of rotation of (a cos ωt, b sin ωt) from the sign of the swept area", "The point (a cos ωt, b sin
    ωt) runs round an ellipse. Which way? Its 'angular momentum' x·dy/dt − y·dx/dt = abω is positive for counter-clockwise
    motion. So (a cos ωt, b sin ωt) with a, b, ω > 0 is counter-clockwise — and (a cos(−ωt), b sin(−ωt)), which is what a
    wave gives at a fixed place, is clockwise.", code="t = np.linspace(0, 1, 5); a, b, w = 2.0, 1.0, -2*np.pi   # the wave's
    phase is -wt\nx, y = a*np.cos(w*t), b*np.sin(w*t)\nprint(np.sign((x[:-1]*np.diff(y) - y[:-1]*np.diff(x)).sum()))   # -1.0:
    clockwise")` (**P322**; reminder of P104, a linear map of a circle is an ellipse)
14. `nb.note` — **N89 [B]** our Fig. 13.15: the velocity ellipse with ωt = 0, π/2, π marked, axes $2\omega\hat\eta/(kH)$ and $2f\hat\eta/(kH)$.
    `> ⚠️ Common confusion (trap T7):` "The book's 'orbit' figure is a velocity hodograph. The particle's path has the same
    shape and sense, with both axes divided by ω. And §13.11 quotes the axis ratio as ω/f (long over short), §13.14 as f/ω
    (short over long) — the same ellipse."
15. `nb.animation` — **A1** (video, 80 frames, FAST 40): left, the current vector and the particle path of one parcel while
    the wavelength grows from 300 km to 30 000 km (`GFD.poincare_orbit`): the ellipse fattens into a circle; right, the
    same in the southern hemisphere, turning the other way; the frequency ratio ω/f in the title falls from 21 toward 1.
    *see / read / change*.
16. `nb.note` — **N90 [B]** inertial motion, the K → 0 limit: $\partial u/\partial t-fv=0$, $\partial v/\partial t+fu=0$, solved by $u=q\cos ft$,
    $v=-q\sin ft$: a circle of radius $r=q/f$ in one inertial period. `GFD.inertial_oscillation(t, 0.23, 0.0, f35)`,
    `GFD.inertial_radius(0.23, f35)` *expect:* radius 2.75 km, period 20.86 h at 35° N. "Seen from space this is Chapter 4's
    particle moving in a straight line while the earth turns under it."
17. `nb.explainer("shallow_water_dispersion", heading="Are Poincaré, Kelvin and Rossby waves separate theories?",
    why="Moving a wavenumber cursor along the diagram while the four term bars of the cubic rebalance shows which terms
    each wave is made of; sliding f to zero or β to zero removes a whole branch.", tries=["Select the slow root and read
    the status: which term is negligible?", "Use the preset 'f → 0': the hyperbola collapses onto ω = cK.", "Switch β off:
    the slow root drops to ω = 0.", "Toggle 'printed vs corrected' to see slip #7 in the regime text."])`
18. `nb.md` — **What would change if…** "…a coast stands in the way, so that the water cannot move across it? A wave with
    no cross-shore velocity at all becomes possible — and it can have any frequency, even below f (C11)."

---

### A.12 §13.12 Kelvin Wave — C11, C12
#### C11 — The Kelvin wave (13.87)
1. `nb.section("13.12", "Kelvin Wave", intro="**What is this section about?** The wave that leans on a coast: it travels
   at the ordinary long-wave speed, only in one direction, and dies away offshore over a distance that turns out to be the
   most important length scale of rotating fluids.")`
2. `nb.core("C11", "The Kelvin wave: $\\eta=\\eta_0e^{-fy/c}\\cos k(x-ct)$, $u=\\eta_0\\sqrt{g/H}\\,e^{-fy/c}\\cos k(x-ct)$ (13.87)",
   question="How can a wave exist below the frequency f, and why does it run only one way along a coast?")`
3. `nb.md` — **Plain words:** "The tide in the North Sea does not slosh in and out; it travels round the basin,
   counter-clockwise, as a wave whose range is largest at the shore. In the Pacific, a relaxation of the trade winds sends
   a signal along the equator and then up and down the coast of the Americas. Both are Kelvin waves: gravity waves that
   use a boundary to dodge the Coriolis force."
4. `nb.md` — **The idea** + **N91 [B]**: "At a wall the water cannot move across the shore, so v = 0 there — try v = 0
   everywhere. Then the Coriolis force on the along-shore current has nothing to turn; it must be *held* by a pressure
   force: the surface tilts across the shore, $fu=-g\dfrac{\partial\eta}{\partial y}$. Under a crest the current flows with the wave, so the surface must
   be high at the wall and fall offshore — and that fixes the direction of travel: coast on the right (northern
   hemisphere, f > 0; on the left for f < 0)." ASCII cross-section.
5. `nb.note` — **N94 [B]** the equations with v ≡ 0: $\dfrac{\partial\eta}{\partial t}+H\dfrac{\partial u}{\partial x}=0$, $\dfrac{\partial u}{\partial t}=-g\dfrac{\partial\eta}{\partial x}$, $fu=-g\dfrac{\partial\eta}{\partial y}$ (13.84). **N93 [C]** "Compare
   the cross-shore equation of the two waves, $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$: for a Poincaré wave the Coriolis term is partly balanced
   by the acceleration ∂v/∂t; for a Kelvin wave v = 0 and it is balanced entirely by the slope — geostrophically."
6. `nb.primer("trapped solutions: keeping the exponential that decays away from a boundary", "A first- or second-order
   equation in the offshore coordinate often has two exponential solutions, one growing and one decaying away from the
   wall. In a half-space only the decaying one is physical (finite energy). Which one decays depends on signs in the
   problem — here on the sign of f and on the direction of travel.", code="y = np.array([0.0, 1.0, 2.0])                 #
   distance from the wall in units of c/|f|\nfor s in (+1, -1):                             # s = sign(f) * direction of travel\n
   print(s, np.exp(-s*y))                     # +1 decays (kept), -1 grows (rejected)")` (**P323**; reminder of P210 first-order
   linear ODE)
7. `nb.derivation("D15", …)` — Part F D15 (8 steps), ref "13.87": **N95 [B]** its steps 2–4: with
   $[u,\eta]=[\hat u(y),\hat\eta(y)]e^{i(kx-\omega t)}$, $-i\omega\hat\eta+iHk\hat u=0$, $-i\omega\hat u=-igk\hat\eta$, $f\hat u=-g\dfrac{d\hat\eta}{dy}$ (13.85), hence $\hat\eta\,[\omega^2-gHk^2]=0$.
8. **(re-enter C11)** `nb.recap("R23", "The Kelvin wave is not dispersive", "$c=\\sqrt{gH}$ (13.86), i.e. $\\omega=\\pm k\\sqrt{gH}$: the non-rotating
   long-wave speed, for every frequency — including below f, where no Poincaré wave exists.", where="Ch. 7, long-wave speed")`
9. `nb.worked_example("a Kelvin wave on a deep ocean", "H = 4000 m, g ≈ 10 m/s², f = 10⁻⁴ s⁻¹, crest height η₀ = 1 m at the
   coast. 1. Speed: c = √(gH) = 200 m/s. 2. Trapping width Λ = c/f = 200/10⁻⁴ = 2 × 10⁶ m = 2000 km. 3. Current under the crest at
   the coast: u = η₀√(g/H) = 1 × √(10/4000) = 0.05 m/s. 4. At y = Λ offshore both are e⁻¹ = 37 % of that. 5. Check the
   cross-shore balance at the coast: fu = 10⁻⁴ × 0.05 = 5 × 10⁻⁶; −g ∂η/∂y = g η₀/Λ = 10/2 × 10⁶ = 5 × 10⁻⁶ ✓.")`
10. `nb.code` — `eta, u = GFD.kelvin_wave(x, y, t, 0.5, 2*np.pi/3.1e6, 4200.0, f35)`; `GFD.kelvin_decay_side(f35, +1)`;
    `GFD.kelvin_decay_side(-f35, +1)`; `GFD.rossby_radius(202.95, f35)`; `GFD.kelvin_omega(2*np.pi/3.1e6, 202.95)`; a call
    with `direction=-1` inside `try/except ValueError`. *expect:* Λ = 2426 km; u = 0.0242 m/s at the coast under the crest;
    period 4.24 h; "trapped": True with the coast on the right for f > 0, travel toward +x; for f < 0 the +x wave is not
    trapped (coast must be on the left); the wrong direction raises. *explain:* 4 points.
11. `nb.check_agree` — **from scratch:** the three equations checked on the formula by finite differences (`deta_dt + H*du_dx`,
    `du_dt + G0*deta_dx`, `f35*u + G0*deta_dy`, each divided by its largest term); `assert` each below 10⁻⁵; compare with
    `GFD.kelvin_residuals(x0, y0, t0, 0.5, k, 4200.0, f35)`.
12. `nb.figure` — **N92 [B]** our Figs. 13.16–13.17 (`ch13.fig_kelvin_sections()`): cross-shore sections of η through a crest
    and a trough, at a single coast and in a channel 2Λ wide (two waves, one on each wall, travelling in opposite
    directions); + `nb.plotly`: `go.Surface` of η(x, y) hugging the coast, dropdown N / S. *see / read* ("the slope across
    the shore reverses between crest and trough, with the current") */ change*.
13. `nb.animation` — **A3** (video, the 24 cached frames): a Kelvin wave travelling round a closed square basin with the wall on its
    right, from the cached C-grid run `ch13.load_reference_run("kelvin_basin")` (64², equivalent depth chosen so that Λ is a
    fifth of the basin; `> 🔧 **Our choice** — our numerical model, cached`). If the cache is absent
    (`load_reference_run` returns None) the straight-coast closed form `GFD.kelvin_wave` is animated instead and the
    caption says so. *see:* the bulge goes round counter-clockwise (f > 0) and turns the corners; *read:* "at every wall
    the coast is on the right of the direction of travel"; *change:* "…f < 0: clockwise".
14. `nb.note` — **N96 [B]** application joining C05 and C11 (`ch13.coastal_upwelling(-0.07, "east", lat)`): "An equatorward
    wind along an eastern ocean boundary drives Ekman transport offshore; deeper water rises; the lifted thermocline is a
    disturbance that leaves **poleward** along the coast as an internal Kelvin wave (coast on its right in the northern
    hemisphere)." Speed and width with our two-layer inputs → C12.
15. `nb.explainer("kelvin_wave", heading="Why does this wave run only one way along a coast?", why="Reversing the
    direction of travel (or the hemisphere) makes the offshore profile grow instead of decay and the status badge reject
    it; widening the channel separates the two coastal waves; switching to the internal mode shrinks the trapping width
    from thousands of kilometres to tens.", tries=["Reverse the direction: the status says 'would grow offshore — not a
    solution'.", "Switch to S and see which way the trapped wave now runs.", "Choose the channel 2Λ wide and watch the two
    walls carry opposite waves.", "Switch to the internal mode and read the new Λ."])`
16. `nb.md` — **What would change if…** "…there is no coast? The width c/f still means something: it is how far a gravity
    wave gets before rotation takes over (C12)."

#### C12 — The Rossby radius of deformation, and geostrophic adjustment
1. `nb.core("C12", "The Rossby radius of deformation $\\Lambda\\equiv c/f$", question="If I pile water up and let go, why doesn't
   it flatten out as it would without rotation — and how wide is what stays?")`
2. `nb.md` — **Plain words:** "Pour a bucket of water into a bath and the bump spreads until the surface is flat. Do the
   same on a planet-sized turntable and it does not: the water starts to spread, the Coriolis force turns the outflow
   into a current running *along* the edge of the bump, and that current's own Coriolis force holds the rest of the bump
   up. The ocean is full of such fronts and eddies standing in slopes that gravity cannot flatten. Their width is the
   Rossby radius — tens of kilometres in the ocean, about a thousand in the atmosphere — and it is why an ocean model
   needs a grid of a few kilometres to 'resolve eddies'."
3. `nb.md` — **The idea:** "In a time 1/f a gravity wave travels a distance c/f. Disturbances narrower than that spread
   before rotation notices them; wider ones are caught by rotation before gravity can flatten them." Table: scale ≪ Λ:
   behaves as non-rotating · scale ≫ Λ: stays, geostrophically balanced.
4. **(re-enter C12)** `nb.recap("R24", "Internal Kelvin waves and the internal radius", "On the interface between a thin upper layer and a
   deep lower one the long-wave speed is $c=\\sqrt{g'H}$ with the reduced gravity $g'=g(\\rho_2-\\rho_1)/\\rho_2$; in a continuously
   stratified layer the speeds are $c=NH/n\\pi$ (C08). The corresponding radius $\\Lambda=NH/(\\pi f)$ for n = 1 is far smaller than
   the external one; the interface moves far more than the surface, in the opposite sense.", where="Ch. 7, $g'=g(\\rho_2-\\rho_1)/\\rho_2$
   (7.117)")` + `> ⚠️ Common confusion (trap T11):` "Three Rossby radii share one name: external √(gH)/f; internal
   √(g′H)/f or NH/(nπf); and in C16 the Eady radius NH/f **without** the π." + > ⚠️ **slip #11 — the book prints** a typical
   internal radius about three times what $NH/(\pi f)$ gives with its own typical N, H and f**; the correct form is** the
   value of the formula — the discrepancy is in the inputs, not the formula (both recorded privately).
5. `nb.worked_example("three radii by hand", "f = 10⁻⁴ s⁻¹. 1. External, H = 4000 m: c = 200 m/s → Λ = 2000 km. 2. Two-layer,
   upper layer 100 m, Δρ/ρ = 0.003: g′ = 0.03 m/s², c = √(0.03 × 100) = √3 ≈ 1.7 m/s → Λ ≈ 17 km. 3. Uniform N = π × 10⁻³ s⁻¹
   over 4000 m: c₁ = NH/π = 4 m/s → Λ = 40 km. 4. Ratio external / internal: 50 to 100 — one ocean, two utterly different
   scales.")`
6. `nb.code` — `GFD.rossby_radius(GFD.long_wave_speed(4200.0), f35)`; `GFD.rossby_radius_internal(2.7e-3, 4200.0, f35)`;
   `GFD.rossby_radius_internal(2.7e-3, 4200.0, f35, with_pi=False)`; `GFD.rossby_radius_two_layer(120.0, 1027.0 - 3.1, 1027.0,
   f35)`; the same at 12° N. *expect:* 2426 km; 43.15 km; 135.6 km (no π — the Eady radius of this ocean); 22.5 km (g′ =
   0.0296 m/s², c = 1.885 m/s ≈ 163 km per day); at 12° N the n = 1 radius is 119.0 km. *explain:* 4 points.
7. `nb.primer("an adjustment problem: what a steady end state can remember (a conserved quantity pins it)", "Release an
   unbalanced state and wait. Waves carry away whatever can travel. What is left must be steady — but there are
   infinitely many steady states. The one nature picks is singled out by a quantity that cannot change at any point
   while the waves pass. Find that quantity, demand it has its initial value, and the end state follows without solving
   for the transient at all.", code="# a bank analogy: transfers move money between branches (waves), the total cannot
   change\nstart = np.array([10.0, 0.0, 0.0]); end = np.full(3, start.sum()/3)\nprint(end, end.sum() == start.sum())")`
   (**P324**)
8. `nb.derivation("D16", …)` — Part F D16 (11 steps), ref "": **N19 [B]** **Not in the book — ours** (the book gives only
   a verbal account of how geostrophy is set up, restated in C02): adjustment of a step $\eta_0\,\mathrm{sgn}(x)$ to the steady
   state $\eta=\eta_0\,\mathrm{sgn}(x)\,(1-e^{-\lvert x\rvert/\Lambda})$ with a jet $v=(g\eta_0/c)\,e^{-\lvert x\rvert/\Lambda}$ along the step. Labelled "analytic (ours)", no
   citation.
9. `nb.code` — `eta_end, v_end = GFD.geostrophic_adjustment_1d(x, 0.05, 1.3286, f35)`; `en = GFD.adjustment_energy(0.05,
   1.3286, f35)`. *expect* (first baroclinic mode as a layer of equivalent depth 1.3286 m, step ±0.05 m, 35° N): Λ = 43.15
   km; jet maximum gη₀/c = 0.136 m/s at the step, along +y for f > 0 and along −y for f < 0; `en["ratio"]` = 1/3 exactly
   ("of the potential energy released, one third stays as kinetic energy of the jet; two thirds leave with the waves").
10. `nb.check_agree` — **from scratch:** Λ three ways by hand (external, two-layer, uniform N) against the three functions;
    the adjusted step from its two exponentials `eta0*np.sign(x)*(1 - np.exp(-abs(x)/Lam))`; `assert np.allclose`; and the
    energy integrals by `np.trapezoid` over ±20Λ: `assert np.isclose(KE/PE, en["ratio"], rtol=1e-3)`.
11. `nb.animation` — **A2** (frames; stops at t = 0, 1, 3, 5 inertial periods; 2 rows): `SW.linear_1d_run` from the step on
    800 cells of 5 km (±2000 km), top row f = 0 (two fronts run off at ±c, the surface between them flat), bottom row
    f = f35 (Poincaré waves leave; a front of width Λ and a jet remain); the closed-form end state as a dashed ghost; Λ
    marked. *expect* (designer's scratch march, same scheme): the instantaneous η within 3Λ of the step differs from the
    end state by 12 %, 6 % and 4 % of η₀ after 1, 3 and 5 inertial periods (near-inertial oscillations die away only
    slowly), while the **mean over the fifth inertial period** matches the closed form to 0.2 % in η and 0.4 % in v. The
    animation's last panel therefore shows the period mean and says so. *see / read / change* ("…the step is replaced by
    a bump much narrower than Λ: almost all of it radiates away").
12. `nb.figure` — the linear potential vorticity $\zeta-f\eta/H$ against x at the first and last step (`out["pv"]`): two curves
    on top of each other; and energy bars from `GFD.adjustment_energy`. *read:* "this is the quantity that pinned the end
    state — C13 makes it exact and nonlinear".
13. `nb.explainer("geostrophic_adjustment", heading="Why doesn't a pile of water flatten out on a rotating planet?",
    why="Pressing play with f = 0 (everything radiates, nothing stays) and then with rotation (waves leave, a front of
    width Λ stays) is the whole idea; the linked bars show how much potential energy was released and how much stayed.",
    tries=["Play the preset 'no rotation', then the mid-latitude preset, and compare the end cards.", "Halve f: the front
    is twice as wide.", "Switch from the step to the narrow bump: what fraction survives?", "Watch the potential-vorticity
    curve: it does not move."])`
14. `nb.md` — **What would change if…** "…the motion is not small and f varies with latitude? The quantity that pinned the
    end state here survives all of that: potential vorticity (C13)."


---

### A.13 §13.13 Potential Vorticity Conservation in Shallow-Water Theory — C13
#### C13 — Potential vorticity: (13.93) and its compact form (13.94)
1. `nb.section("13.13", "Potential Vorticity Conservation in Shallow-Water Theory", intro="**What is this section about?**
   The one conservation law behind all slow, large-scale motion: every column of fluid carries the ratio of its absolute
   vorticity to its depth unchanged. Stretch it, squash it or move it to another latitude, and its spin must respond.")`
2. `nb.core("C13", "Potential vorticity: $\\dfrac{D(\\zeta+f)}{Dt}=\\dfrac{\\zeta+f_0}{h}\\dfrac{Dh}{Dt}$ (13.93), hence $\\dfrac{D}{Dt}\\Big(\\dfrac{\\zeta+f}{h}\\Big)=0$ with $f=f_0+\\beta y$ (13.94)",
   question="What does a column of fluid keep as it moves across an ocean of changing depth and latitude?")`
3. `nb.md` — **Plain words:** "A skater who pulls her arms in spins faster. A column of ocean that is stretched taller
   gets thinner (its volume is fixed) and spins faster too — except that the column already carries a spin it never
   chose: the earth's, f, which depends on latitude. So there are two ways to change a column's own spin ζ: change its
   height, or carry it north or south. Westerlies crossing the Rockies are squashed, turn toward the equator, overshoot
   and meander downstream for thousands of kilometres. Dynamicists treat potential vorticity as the dye that marks
   large-scale air and water masses."
4. `nb.md` — **The idea** + **N104 [C]** the geometry sketch `draw_pv_column()` (a column of depth h over an uneven bottom,
   η measured from a reference level; a thin tall column beside a short fat one of the same volume). Table: stretch
   (h ↑) → ζ + f ↑ · move poleward (f ↑ in the north) → ζ ↓ · both hemispheres in separate rows.
5. `nb.note` — **N97 [B]** $\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}-fv=-g\dfrac{\partial\eta}{\partial x}$ (13.88); **N98 [B]** $\dfrac{\partial v}{\partial t}+u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.89);
   **N99 [B]** $\dfrac{\partial h}{\partial t}+\dfrac{\partial}{\partial x}(uh)+\dfrac{\partial}{\partial y}(vh)=0$ (13.90), with h the total depth over an uneven bottom and $f=f_0+\beta y$. "These three
   are the Start of D17; unlike C07 they are nonlinear."
6. **(re-enter C13)** `nb.recap("R25", "Relative and absolute vorticity", "$\\zeta\\equiv\\dfrac{\\partial v}{\\partial x}-\\dfrac{\\partial u}{\\partial y}$ is the vertical vorticity measured
   by someone turning with the earth (relative); adding the planet's own, ζ + f, gives the absolute vorticity.", where="Ch. 3
   and Ch. 5 §5.6")` + `nb.primer("materially conserved: Dq/Dt = 0 labels a parcel; it does not mean q is steady at a
   point", "D/Dt follows one parcel. Dq/Dt = 0 says: each parcel keeps its own value of q for ever, like a dye. At a
   fixed place q can still change, because parcels with different values pass by. Steady means ∂q/∂t = 0 — a different
   statement.", code="x = np.linspace(0, 10, 11); q0 = lambda x: x**2      # each parcel's label\nU, t = 2.0, 1.5\nprint(q0(x -
   U*t)[5], q0(x)[5])   # at the fixed point x = 5: 4.0 now, 25.0 before — not steady, yet every parcel kept its q")` (**P325**;
   reminders of P38 product rule, P265 quotient rule, the material derivative $\frac{D}{Dt}=\frac\partial{\partial t}+u\frac\partial{\partial x}+v\frac\partial{\partial y}$ of Ch. 3).
7. `nb.derivation("D17", …)` — Part F D17 (12 steps), ref "13.93": **N100 [B]** its steps 1–3, the cross-differentiated
   momentum equation $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial v}{\partial x}-\dfrac{\partial u}{\partial y}\Big)+\dfrac{\partial}{\partial x}\Big[u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}\Big]-\dfrac{\partial}{\partial y}\Big[u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}\Big]+f_0\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)+\beta v=0$ (13.91); **N101 [B]** step 6,
   $\dfrac{D\zeta}{Dt}+(\zeta+f_0)\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)+\beta v=0$ (13.92). `> ⚠️ Common confusion (trap T10):` "The book replaces f by f₀ in the
   coefficient and then restores it. The exact shallow-water result needs no such step — our check cell keeps f = f₀ + βy
   throughout."
8. **(re-enter C13)** `nb.recap("R26", "Potential vorticity is conserved following the motion", "$\\dfrac{D}{Dt}\\Big(\\dfrac{\\zeta+f}{h}\\Big)=0$, $f=f_0+\\beta y$
   (13.94) — the boxed result of D17. Chapter 5 obtained the f-plane version from Kelvin's theorem for a column; here it
   comes from the full nonlinear equations, with f varying.", where="Ch. 5 §5.6, derivation D20: $(\\omega_z+2\\Omega)/h$ constant
   for a column")` + `ch13.pv_conservation_sympy()["ok"]` → True.
9. **(re-enter C13)** `nb.recap("R27", "How to read it", "Stretching creates relative vorticity even from ζ = 0, because the
   column carries f. There is no tilting term: a shallow layer moves in vertical columns. An increase of h makes ζ + f
   more positive in the northern hemisphere (f > 0) and more negative in the southern (f < 0) — larger in magnitude in
   both.", where="Ch. 5 §5.6 (C11)")`
10. `nb.worked_example("a column crossing a ridge", "f = 10⁻⁴ s⁻¹; a column with no relative vorticity (ζ = 0) and depth h₀ =
    4000 m moves onto a plateau where h₁ = 3600 m, without changing latitude. 1. Before: q = (0 + f)/h₀ = 10⁻⁴/4000 = 2.5 ×
    10⁻⁸ m⁻¹ s⁻¹. 2. After: (ζ + f)/h₁ = q → ζ + f = 2.5 × 10⁻⁸ × 3600 = 0.9 × 10⁻⁴. 3. So ζ = −0.1 × 10⁻⁴ = −0.1 f: clockwise
    (anticyclonic) spin in the northern hemisphere. 4. Same formula in one step: $\\zeta=f(h_1-h_0)/h_0$ = 10⁻⁴ × (−400/4000). 5.
    Instead keep the depth and move the column 500 km north where f is larger by βΔy = 2 × 10⁻¹¹ × 5 × 10⁵ = 10⁻⁵: ζ = −10⁻⁵
    again. A 10 % squash equals a 500 km trip north.")`
11. `nb.code` — `q0 = GFD.potential_vorticity(0.0, f35, 9000.0)`; `z1 = GFD.step_vorticity(f35, 9000.0, 8550.0)`;
    `ch05.column_relative_vorticity(8550.0, 9000.0, 0.0, f35)` (the Chapter 5 function, same number); the southern twin with
    `-f35`. *expect* (a 9000 m layer squashed by 5 % at 35° N): ζ = −4.183 × 10⁻⁶ s⁻¹ = −0.05 f (clockwise); at 35° S, +4.183 ×
    10⁻⁶ s⁻¹ (counter-clockwise, which is anticyclonic there too). *explain:* 3 points.
12. `nb.check_agree` — **from scratch:** `q_before = (0.0 + f35)/9000.0`; `zeta_after = q_before*8550.0 - f35`; `assert
    np.isclose(zeta_after, GFD.step_vorticity(f35, 9000.0, 8550.0))`; `assert np.isclose(GFD.potential_vorticity(zeta_after,
    f35, 8550.0), q_before)`.
13. `nb.derivation("D18", …)` — Part F D18 (7 steps), ref "": **N102 [B]** eastward flow over a step — the vorticity just
    downstream from $\dfrac{f}{h_0}=\dfrac{\zeta+f}{h_1}$, i.e. $\zeta=\dfrac{f(h_1-h_0)}{h_0}<0$, and a standing meander of wavelength $\lambda=2\pi\sqrt{U/\beta}$; **N103 [B]**
    westward flow — no oscillation, exponential on both sides. **Ours** (the book argues in words): the streamline
    equation $Y''+(\beta/U)Y=\text{const}$ and its solutions.
14. `nb.figure` — `e = ch13.flow_over_step(x, 17.0, beta35, f35, 9000.0, 8550.0)` and `w = ch13.flow_over_step(x, -17.0, …)`:
    two panels of streamlines (several starting latitudes, each displaced by Y(x)), the step shaded; the relative
    vorticity ζ(x) below each. *expect* (ours, computed: U = ±17 m/s, 5 % step, 35° N): eastward — the stream swings between
    Y = 0 and Y = −446 km about the mean −223 km, wavelength 5983 km, steepest slope 0.23; westward — deflection begins
    *before* the step with e-folding length √(|U|/β) = 952 km, is −112 km at the step and tends to −223 km far
    downstream. *see / read* ("eastward: a wake of stationary waves behind the step; westward: the step is felt
    upstream and there is no wake") */ change* ("…U were four times larger: the meander is twice as long").
    `> ⚠️ **Where our solution and the book's sketch differ:**` "The book draws the westward stream back at its original
    latitude far downstream. With a step, potential-vorticity conservation does not allow that: a column with ζ = 0 over
    the shallower depth must sit where f is smaller by the same 5 %, i.e. 223 km equatorward. It returns to its latitude
    only if the depth returns too (a ridge instead of a step)." The callout is titled, in the notebook, **"where our
    potential-vorticity solution and the book's description differ"** with $Y_p=f_0(h_1-h_0)/(\beta h_0)$ written out; by the
    orchestrator's ruling it is a labelled remark, **not** a slip, unless the math-verifier confirms it independently.
15. `nb.figure` — nonlinear check from the cache: `pv = ch13.load_reference_run("pv_particles")` (our 48² nonlinear C-grid
    run with 24 particles at 13 saved times, key `"q"`; `> 🔧 **Our choice** — our numerical model, cached`): potential vorticity of each particle
    against time (flat lines) beside its relative vorticity and local depth (both wandering). If the cache is absent the
    cell prints one sentence and the closed-form figure above stands. *read:* "ζ and h each change by tens of per cent;
    their combination (ζ + f)/h does not".
16. `nb.md` — **What would change if…** "…the stratification is continuous and we ask for waves of any slope, fast or slow?
    First the fast ones (C14); the slow ones, which are this conservation law in motion, are Rossby waves (C15)."

---

### A.14 §13.14 Internal Waves — C14
#### C14 — Inertia–gravity waves (13.112)
1. `nb.section("13.14", "Internal Waves", intro="**What is this section about?** Chapter 7's internal gravity waves with
   the earth's rotation added. The frequency is still set by the direction of the wavevector alone, but now it is caught
   between f and N. ⚠️ In this section the book's θ is the angle of the wavevector with the horizontal, not latitude, and
   its H is the scale over which N changes, not a depth.")`
2. `nb.core("C14", "Inertia–gravity waves: $\\omega^2-f^2=\\dfrac{k^2}{m^2}(N^2-\\omega^2)$ (13.112), equivalently $\\omega^2=f^2\\sin^2\\theta+N^2\\cos^2\\theta$ (unnumbered, p. 669)",
   question="What does rotation change about the internal waves of Chapter 7?")`
3. `nb.md` — **Plain words:** "Lower a current meter into the thermocline and the record is full of oscillations with
   periods from about ten minutes to most of a day. The short end is the buoyancy period 2π/N; the long end is the
   inertial period 2π/f. In between, every period belongs to a wave whose crests have a definite slope: steep crests,
   fast wave; nearly flat crests, slow wave, almost an inertial circle. These waves carry the energy of wind and tide
   into the deep ocean, where their breaking does much of the mixing that ocean models have to parameterise."
4. `nb.md` — **The idea:** two springs again: buoyancy (stiffness N²) resists vertical motion, rotation (stiffness f²)
   resists horizontal motion. Fluid moves along the crests; if the crests are vertical the motion is vertical and only
   buoyancy acts (ω = N); if they are horizontal only rotation acts (ω = f); in between, a weighted mix. Table θ_K = 0°,
   45°, 90° → ω.
5. **(re-enter C14)** `nb.recap("R28", "Internal waves without rotation", "Anisotropic; ω ≤ N; frequency set by the angle
   of the wavevector, $\\omega=N\\cos\\theta$; phase and group velocity perpendicular, with opposite vertical components.",
   where="Ch. 7 §7.8 (C15–C16)")` · **(re-enter C14)** `nb.recap("R29", "The linear rotating Boussinesq set", "Continuity
   $\\dfrac{\\partial u}{\\partial x}+\\dfrac{\\partial v}{\\partial y}+\\dfrac{\\partial w}{\\partial z}=0$; $\\dfrac{\\partial u}{\\partial t}-fv=-\\dfrac1{\\rho_0}\\dfrac{\\partial p}{\\partial x}$; $\\dfrac{\\partial v}{\\partial t}+fu=-\\dfrac1{\\rho_0}\\dfrac{\\partial p}{\\partial y}$; $\\dfrac{\\partial w}{\\partial t}=-\\dfrac1{\\rho_0}\\dfrac{\\partial p}{\\partial z}-\\dfrac{\\rho g}{\\rho_0}$;
   $\\dfrac{\\partial\\rho}{\\partial t}-\\dfrac{\\rho_0N^2}{g}w=0$ — the five members of (13.95). New against Chapter 7: the two Coriolis terms. Not hydrostatic.",
   where="Ch. 7, the same set with f = 0, ending in $\\frac{\\partial^2}{\\partial t^2}\\nabla^2w+N^2\\nabla_H^2w=0$ (7.134)")` +
   `ch13.rotating_internal_wave_set_sympy()["ok"]` → True.
6. `nb.derivation("D19", …)` — Part F D19 (11 steps), ref "13.96": **N105 [B]** the single equation for w,
   $\dfrac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w+f^2\dfrac{\partial^2w}{\partial z^2}=0$ (13.96), written out because the book refers to its §7.8 and never shows the steps
   with rotation. > ⚠️ **slip #8 — the book prints** that N is taken "depth independent"**; the correct form is** depth
   dependent, N(z), as its own definition $m^2(z)\equiv\frac{(k^2+l^2)[N^2(z)-\omega^2]}{\omega^2-f^2}$ (13.99) requires. +
   `ch13.w_equation_rotating_sympy()["ok"]` → True; `ch13.w_equation_rotating_residual(w_fn, …)` on a plane wave → 10⁻⁹
   relative.
7. `nb.derivation("D20", …)` — Part F D20 (7 steps), ref "13.112": **N106 [C]** its step 1, the trial solution
   $[u,v,w]=[\hat u(z),\hat v(z),\hat w(z)]\,e^{i(kx+ly-\omega t)}$ (13.97); **N107 [B]** step 2, $\dfrac{d^2\hat w}{dz^2}+\dfrac{(N^2-\omega^2)(k^2+l^2)}{\omega^2-f^2}\hat w=0$ (13.98); **N108 [B]** step 3,
   $m^2(z)\equiv\dfrac{(k^2+l^2)[N^2(z)-\omega^2]}{\omega^2-f^2}$ (13.99); **N109 [B]** the oscillator form $\dfrac{d^2\hat w}{dz^2}+m^2\hat w=0$ (13.100); **N119 [C]** step 4, the
   l = 0 form $m^2=\dfrac{k^2(N^2-\omega^2)}{\omega^2-f^2}$ (13.109). `GFD.inertia_gravity_m2(k, 0.0, omega, N, f35)`.
8. `nb.note` — **N110 [B]** "Internal waves exist only for $f<\omega<N$: there m² > 0 and the structure oscillates in z; outside
   the band m² < 0 and it decays (evanescent)." `GFD.inertia_gravity_band(omega, 2.7e-3, f35)` for ω = 0.5f, 5f, 1.2N →
   "below f", "in band", "above N" (badge text built from the dict, never a raw library string).
9. `nb.worked_example("frequency from the slope of the crests", "f = 10⁻⁴ s⁻¹, N = 10⁻² s⁻¹ (N/f = 100). 1. Wavevector 60°
   from the horizontal: ω² = f² sin² 60° + N² cos² 60° = 0.75 × 10⁻⁸ + 0.25 × 10⁻⁴ ≈ 0.25 × 10⁻⁴ → ω = 5 × 10⁻³ s⁻¹ = N/2 (period 21
   min): rotation is invisible. 2. Wavevector almost vertical, m/k = 100 (tan θ = 100): cos² θ ≈ 10⁻⁴, sin² θ ≈ 1 → ω² = 10⁻⁸ +
   10⁻⁴ × 10⁻⁴ = 2 × 10⁻⁸ → ω = 1.41 f (period 12.3 h): half rotation, half buoyancy. 3. With m/k = 1000: ω = 1.005 f — an inertial
   oscillation.")`
10. `nb.code` — `for th in (5, 45, 85, 89): GFD.inertia_gravity_omega(np.cos(np.deg2rad(th)), np.sin(np.deg2rad(th)), 2.7e-3,
    f35)`; `k, m = 2*np.pi/1.0e4, 2*np.pi/200.0`; `w = GFD.inertia_gravity_omega(k, m, 2.7e-3, f35)`; `cg =
    GFD.inertia_gravity_group_velocity(k, m, 2.7e-3, f35)`; the f = 0 call against `WAV.internal_wave_omega(k, m, 2.7e-3)`.
    *expect* (N = 2.7 × 10⁻³ s⁻¹, 35° N; N/f = 32.3): ω/N = 0.996 at 5°; ω = 22.8 f at 45°; 2.98 f at 85°; 1.148 f at 89°; for a
    wave 10 km long and 200 m tall ω = 9.955 × 10⁻⁵ s⁻¹ = 1.190 f (period 17.5 h), c_g = (+0.0465, −9.31 × 10⁻⁴) m/s, phase
    speed 3.17 mm/s; the dot product of phase and group velocity is zero to round-off. *explain:* 4 points.
11. `nb.primer("group velocity as the gradient of ω in wavenumber space, read off a contour plot", "In one dimension the
    energy of a wave packet travels at dω/dk. In two or three, at the vector (∂ω/∂k, ∂ω/∂l, ∂ω/∂m): the gradient of ω in
    wavenumber space. On a contour plot of ω over the (k, m) plane it points straight uphill, at right angles to the
    contours — and its length is how crowded the contours are.", code="w = lambda k, m: np.sqrt((k**2 + 0.01*m**2)/(k**2 +
    m**2))   # N = 1, f = 0.1\nh = 1e-6; k0, m0 = 1.0, 2.0\nprint((w(k0+h, m0)-w(k0-h, m0))/(2*h), (w(k0, m0+h)-w(k0, m0-h))/(2*h))
    # (0.354, -0.177): energy goes forward and DOWN while crests go forward and up")` (**P326**)
12. `nb.derivation("D21", …)` — Part F D21 (7 steps), ref "": **N124 [B]** the group velocity
    $[c_{gx},c_{gz}]=\dfrac{(N^2-f^2)\,km}{(m^2+k^2)^{3/2}(m^2f^2+k^2N^2)^{1/2}}\,[m,\,-k]$ (the book leaves it to an exercise): phase velocity along the
    wavevector, group velocity at right angles to it with the opposite vertical sign, fluid motion along the group
    velocity.
13. `nb.check_agree` — **from scratch:** `w_mine = np.sqrt(f35**2*np.sin(th)**2 + N**2*np.cos(th)**2)` with `th =
    np.arctan2(m, k)`; group velocity by central differences of `GFD.inertia_gravity_omega` in k and m; `assert
    np.isclose(w_mine, w)`, `assert np.allclose(cg_fd, cg, rtol=1e-6)`.
14. `nb.figure` — **N123 [B]** our Fig. 13.25: ω against k for three values of m, rising from f to N on logarithmic axes;
    the three regimes marked with their approximate relations: high frequency $m^2\simeq\dfrac{k^2(N^2-\omega^2)}{\omega^2}$ (non-rotating, ω = N cos θ);
    low frequency $\omega^2\simeq f^2+\dfrac{k^2N^2}{m^2}$ (hydrostatic); mid frequency $m^2\simeq\dfrac{k^2N^2}{\omega^2}$ (both). + a three-row table of
    `GFD.inertia_gravity_regime(omega, N, f35)` errors at ω = 2f, √(fN), N/2. *see / read / change*.
15. `nb.primer("slowly varying medium (WKB): amplitude and phase ansatz, valid when the medium changes little in one
    wavelength", "If the 'stiffness' of an oscillator equation ŵ″ + m²(z)ŵ = 0 changes slowly, the solution still looks
    locally like a wave, with a local wavenumber m(z) and a slowly changing amplitude. Write ŵ = A(z)e^{iφ(z)}, and demand
    that A changes little over one wavelength. It fails where m → 0 (a turning point: the wave reflects).", code="z =
    np.linspace(0, 50, 5001); m = 1 + 0.02*z                     # wavenumber doubles over 50 units\nphase = np.cumsum(m)*(z[1]-z[0])\nw_wkb
    = np.cos(phase)/np.sqrt(m)                                 # amplitude ~ m^(-1/2)\nprint(w_wkb[0], np.abs(w_wkb[-200:]).max())
    # 1.0 at the start, 0.71 at the end")` (**P327**)
16. `nb.note` — the WKB solution, stated (the book prints its steps): **N111 [B]** the condition $H_Nm\gg1$ (the book writes
    $Hm\gg1$, with H the scale over which N varies); **N112 [C]** substituting $\hat w=A(z)e^{i\phi(z)}$ gives
    $\dfrac{d^2A}{dz^2}+A\Big[m^2-\Big(\dfrac{d\phi}{dz}\Big)^2\Big]=0$ and $2\dfrac{dA}{dz}\dfrac{d\phi}{dz}+A\dfrac{d^2\phi}{dz^2}=0$ (13.101)–(13.102); **N113 [B]** neglecting A″ gives the eikonal
    $\dfrac{d\phi}{dz}=\pm m$, $\phi=\pm\displaystyle\int^zm\,dz$ (13.103): the phase gradient is the local wavenumber; **N114 [B]** the second equation then
    integrates to $\hat w=\dfrac{A_0}{\sqrt m}\,e^{\pm i\int^zm\,dz}$ (13.104) — "the amplitude grows where the wave is long because the vertical
    energy flux must be the same at every level". + `nb.figure`: `ch13.vertical_structure_solve(z, N_fn, k, omega, f35)`
    (numerical reference) against `GFD.wkb_vertical_structure(z, m)` for the thermocline profile, and `ch13.wkb_error` for
    three wavelengths. *expect* (measured by the implementer and re-run by the designer on its test profile N = N₀(1 +
    0.5 z/D), N₀ = 2.7 × 10⁻³ s⁻¹, k = 10⁻³ m⁻¹, ω = 4 × 10⁻⁴ s⁻¹, 35° N, D = 0.5, 1, 2, 4 km): $H_Nm$ = 1.65, 3.30, 6.59, 13.2 and
    errors 0.104, 0.045, 0.020, 0.0094 — the error falls as $1/(H_Nm)$ (error × $H_Nm$ ≈ 0.17 → 0.12), **not** as its square.
    The figure plots error against $H_Nm$ on log axes with a slope −1 guide.
17. `nb.note` — the velocity field (stated): **N115 [C]** continuity for l = 0, $ik\hat u+\dfrac{d\hat w}{dz}=0$ (13.105); **N116 [C]**
    $\hat u=\mp\dfrac{A_0\sqrt m}{k}e^{\pm i\int^zm\,dz}$ (13.106) ("⚠️ trap T15: √m is treated as constant when ŵ is differentiated"); **N117 [C]**
    $\hat v=\pm\dfrac{if}{\omega}\dfrac{A_0\sqrt m}{k}e^{\pm i\int^zm\,dz}$ (13.107), from $\hat u/\hat v=i\omega/f$; **N118 [B]** the real fields
    $u=\mp\dfrac{A_0\sqrt m}{k}\cos\Big(kx\pm\displaystyle\int^zm\,dz-\omega t\Big)$, $v=\mp\dfrac{A_0f\sqrt m}{\omega k}\sin\Big(kx\pm\displaystyle\int^zm\,dz-\omega t\Big)$, $w=\dfrac{A_0}{\sqrt m}\cos\Big(kx\pm\displaystyle\int^zm\,dz-\omega t\Big)$ (13.108), upper signs for
    upward phase propagation. `GFD.inertia_gravity_fields(x, z, t, k, omega, N_fn, f35)`; continuity residual at round-off.
18. `nb.note` — **N120 [B]** the horizontal hodograph at a point, $u=\mp\cos\omega t$, $v=\pm\dfrac f\omega\sin\omega t$ (13.110): a clockwise
    ellipse in the northern hemisphere (f > 0; counter-clockwise for f < 0), long axis along the direction of
    propagation, axis ratio f/ω (trap T7: the same ellipse as C10's ω/f, quoted the other way up).
    `GFD.inertia_gravity_hodograph(t, omega, f35)`. **N121 [B]** the motion lies along the phase lines:
    $\dfrac uw=\mp\dfrac mk=\mp\tan\theta$ (13.111), with $\theta=\tan^{-1}(m/k)$ the angle of the wavevector with the horizontal (**not**
    latitude — trap T13); checked: `k*u + m*w` at round-off.
19. `nb.figure` — **N122 [B]** our Figs. 13.22–13.24 (`ch13.fig_inertia_gravity_orbit()`): (a) hodograph and the tilted
    plane of the orbit; (b) a vertical section with phase lines, phase velocity and group velocity at right angles, u
    along the crests; + `nb.plotly`: the helix traced by the velocity vectors with depth (`go.Scatter3d`), dropdown
    upward / downward phase: clockwise with depth for upward phase (northern hemisphere). *see / read / change*.
20. `nb.note` — **N125 [B]** lee waves: with f negligible, $\omega^2=\dfrac{N^2k^2}{m^2+k^2}$ (13.113); a wave that stands still over the
    ground has its intrinsic frequency Doppler-shifted to zero, $\omega_0=\omega+\mathbf K\cdot\mathbf U=0$ (Chapter 7's $\omega_0=\omega+\mathbf U\cdot\mathbf K$ (7.9)),
    so ω = kU and $U=\dfrac{N}{\sqrt{k^2+m^2}}$: only waves with k < N/U exist. `GFD.lee_wave_m(17.0, 1.1e-2, 2*np.pi/2.0e4)` *expect* (U =
    17 m/s, N = 1.1 × 10⁻² s⁻¹): shortest stationary wavelength 2πU/N = 9.71 km; a 20 km wave has m = 5.66 × 10⁻⁴ m⁻¹ (vertical
    wavelength 11.1 km); a 5 km wave raises "evanescent". (Trap T15: k here is a magnitude; the wavevector points
    upstream.) **N126 [B]** our Fig. 13.26: `ch13.lee_wave_field(x, z, 17.0, 1.1e-2, 2*np.pi/2.0e4, 300.0)` streamlines; the
    line through the crests tilts upstream with height. *see / read / change*.
21. `nb.md` — **What would change if…** "…the frequency is far below f? Then none of these waves is possible — but the
    slow root of C09's cubic is. It needs β (C15)."

---

### A.15 §13.15 Rossby Wave — C15
#### C15 — Rossby waves (13.118)
1. `nb.section("13.15", "Rossby Wave", intro="**What is this section about?** The slow wave that owes its existence to the
   change of the Coriolis parameter with latitude. Its crests always drift west; its energy can go either way.")`
2. `nb.core("C15", "Rossby waves: $\\omega=-\\dfrac{\\beta k}{k^2+l^2+f_0^2/c^2}$ (13.118)", question="The crests go west — so how can a
   storm track's energy go east, and how long does the ocean take to hear about a change in the wind?")`
3. `nb.md` — **Plain words:** "The jet stream meanders in four or five great waves round the hemisphere; when one of
   them stalls, a region gets weeks of the same weather. In the ocean, a change of wind in the east Pacific is felt in
   the west a year later at low latitudes, a decade later at mid-latitudes. Both are Rossby waves. They are not held up
   by gravity like the waves of C10 but by the conservation of potential vorticity on a planet where f grows toward the
   pole."
4. `nb.md` — **The idea** (the mechanism, with C13's law): a line of columns along a latitude circle; push one north: f
   is larger there, so its ζ must drop — it spins clockwise; push its neighbour south: it spins counter-clockwise. The
   two spins between them push the fluid *west of the northern bulge* north and the fluid east of it south: the whole
   pattern shifts **west**. ASCII row of ↻ ↺.
5. `nb.note` — **N127 [B]** "Rossby (planetary) waves exist only because of β, at ω ≪ f. The motion is
   *quasi-geostrophic*: geostrophic to lowest order, with its slow evolution set by the small departures from
   geostrophy." The observed height map of the book is replaced by a synthetic wavenumber-5 height field with its
   geostrophic wind (`GFD.geostrophic_from_height`) — figure.
6. `nb.primer("ordering in a small parameter: lowest order gives the balance, next order gives the evolution", "When a
   small number ε (here the Rossby number, or ω/f) multiplies some terms, sort every term by how many factors of ε it
   carries. The largest terms must balance among themselves — that gives a *diagnostic* relation with no time
   derivative (geostrophy). To learn how things change you must go to the next size of terms. So the small terms cannot
   all be thrown away: the right ones decide the future.", code="eps = 0.1\nprint(1.0, eps, eps**2)        # sizes of:
   Coriolis and pressure | acceleration and the divergent wind | products of wave amplitudes")` (**P328**)
7. `nb.derivation("D22", …)` — Part F D22 (10 steps), ref "13.117": **N128 [B]** its step 1,
   $(H+\eta)\Big(\dfrac{\partial\zeta}{\partial t}+u\dfrac{\partial\zeta}{\partial x}+v\dfrac{\partial\zeta}{\partial y}+\beta v\Big)-(\zeta+f_0)\Big(\dfrac{\partial\eta}{\partial t}+u\dfrac{\partial\eta}{\partial x}+v\dfrac{\partial\eta}{\partial y}\Big)=0$ (13.114); **N129 [B]** step 4, $H\dfrac{\partial\zeta}{\partial t}+H\beta v-f_0\dfrac{\partial\eta}{\partial t}=0$ (13.115); **N130
   [B]** steps 5–6, $u\simeq-\dfrac g{f_0}\dfrac{\partial\eta}{\partial y}$, $v\simeq\dfrac g{f_0}\dfrac{\partial\eta}{\partial x}$ (13.116), hence $\zeta=\dfrac g{f_0}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}\Big)$; **N131 [B]** the result,
   $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}-\dfrac{f_0^2}{c^2}\eta\Big)+\beta\dfrac{\partial\eta}{\partial x}=0$, $c=\sqrt{gH}$ (13.117). + `ch13.qg_vorticity_sympy()["ok"]` → True.
8. `nb.derivation("D23", …)` — Part F D23 (9 steps), ref "13.118": **N134 [B]** its step 3, $c_x=\dfrac\omega k=-\dfrac{\beta}{k^2+l^2+f_0^2/c^2}$ (13.119) —
   westward for every wave; **N132 [B]** steps 4–8: the circles $\Big(k+\dfrac\beta{2\omega}\Big)^2+l^2=\Big(\dfrac\beta{2\omega}\Big)^2-\dfrac{f_0^2}{c^2}$, the group velocity
   $\mathbf c_g=\mathbf e_x\dfrac{\partial\omega}{\partial k}+\mathbf e_y\dfrac{\partial\omega}{\partial l}$ with (ours) $c_{gx}=\dfrac{\beta(k^2-l^2-f_0^2/c^2)}{(k^2+l^2+f_0^2/c^2)^2}$, $c_{gy}=\dfrac{2\beta kl}{(k^2+l^2+f_0^2/c^2)^2}$, and the maximum frequency
   $\omega_{max}=\beta c/(2f_0)$ at $kc/f_0=-1$; **N136 [B]** step 9, a mean current: $c_x=U-\dfrac{\beta}{k^2+l^2+f_0^2/c^2}$ (13.120), stationary for
   $\lambda=2\pi\sqrt{U/\beta}$. `> ⚠️ Common confusion (trap T8):` "With ω taken positive the zonal wavenumber is negative. The
   'maximum' phase speed is a maximum of magnitude. The group velocity is westward for long waves and eastward for
   short ones; drawn on the constant-ω circles its arrows point inward."
9. `nb.worked_example("two Rossby waves by hand", "β = 2 × 10⁻¹¹ m⁻¹ s⁻¹. **(a) A barotropic atmospheric wave**, wavelength 6300
   km (k = −10⁻⁶ m⁻¹), l = 0, f₀²/c² negligible: 1. $\\omega=-\\beta k/k^2$ = β/abs(k) = 2 × 10⁻⁵ s⁻¹ → period 3.6 days. 2. $c_x=-\\beta/k^2$ =
   −20 m/s: westward at 20 m/s relative to the air. 3. In a westerly wind of 20 m/s it stands still — a stationary wave
   ($U=\\beta/k^2$). **(b) A long baroclinic ocean wave** with Rossby radius Λ = 40 km: 4. $c_x\\simeq-\\beta\\Lambda^2$ = −2 × 10⁻¹¹ × 1.6 × 10⁹ = −0.032
   m/s ≈ −2.8 km per day. 5. An ocean 8000 km wide is crossed in 8 × 10⁶/0.032 = 2.5 × 10⁸ s ≈ 8 years.")`
10. `nb.code` — `k = -2*np.pi/3.1e6`; `GFD.rossby_omega(k, 0.0, beta35, f35, 202.95)`; `GFD.rossby_phase_speed(...)`;
    `GFD.rossby_group_velocity(...)`; the same with `c=3.6096`; `GFD.rossby_max_frequency(beta35, f35, 3.6096)`;
    `GFD.rossby_long_wave_speed(beta35, f35, 3.6096)`; `GFD.rossby_long_wave_speed(GFD.beta_parameter(inp["lat_rossby"]), f12,
    3.6096)`; `ch13.basin_crossing_time(1.0e7, inp["lat_rossby"], 3.6096)`; `GFD.stationary_rossby_wavelength(17.0, beta35)`;
    `GFD.rossby_omega_circle(1.0e-6, beta35, f35, 202.95)`. *expect* (35° N, wavelength 3100 km): external mode ω = 8.884 × 10⁻⁶
    s⁻¹ (period 8.19 days), c_x = −4.38 m/s, c_gx = **+4.04 m/s** (eastward: this wave is shorter than 2πΛ = 15 240 km);
    first baroclinic mode ω = 7.02 × 10⁻⁸ s⁻¹ (2.8 years), c_x = −3.47 cm/s, c_gx = **−3.41 cm/s** (westward: longer than 2πΛ =
    271 km); ω_max = 4.05 × 10⁻⁷ s⁻¹ (shortest period 180 days) at 35° N; long-wave speed −3.49 cm/s at 35° N and −31.7 cm/s
    at 12° N; a 10 000 km basin is crossed in 1.00 years at 12° N and 9.07 years at 35° N; stationary wavelength 5983 km in
    a 17 m/s westerly. *explain:* 6 points. **N135 [B]** is the long-wave pair of this cell: $c_x\simeq-\dfrac{\beta c^2}{f_0^2}=-\beta\Lambda^2$,
    non-dispersive and westward. > ⚠️ **slip #12 — the book prints** an exercise answer for this speed that follows only
    from a round value of β**; the correct form is** about a fifth smaller with β evaluated at the stated latitude (a
    remark on rounding; the exercise is not reproduced).
11. `nb.check_agree` — **from scratch:** `w_mine = -beta35*k/(k**2 + l**2 + f35**2/c**2)`; group velocity by complex-step
    differences (`np.imag(omega(k + 1j*h, l))/h`, h = 10⁻²⁰); `assert np.isclose(w_mine, …)`, `assert np.allclose(cg_cs, cg,
    rtol=1e-12)`.
12. `nb.figure` — **N133 [B]** our Fig. 13.28 (`ch13.fig_rossby_dispersion()`): upper, $\omega f_0/(\beta c)$ against $kc/f_0$ for l = 0
    (maximum ½ at −1, the sign of c_gx marked on either side, the non-dispersive corner near k = 0); lower, three circles
    of constant ω in the $(kc/f_0,\ lc/f_0)$ plane with group-velocity arrows pointing inward. *see / read* ("left of the
    maximum the curve rises toward the origin: short waves, energy eastward; right of it: long waves, energy westward")
    */ change* ("…the Rossby radius were infinite (barotropic, rigid lid): no maximum — every wave sends its energy
    east").
13. `nb.plotly` — **F8** (`slider_figure`, slider = Rossby radius from 20 km to 3000 km, log, 20 steps): ω(k) and c_gx(k) for
    l = 0 with the maximum and the sign change marked; trace names carry 2πΛ. *see / read / change*.
14. `nb.note` — **N137 [B]** "The cubic of C09 again: $\omega^3-c^2\omega(k^2+l^2)-f_0^2\omega-c^2\beta k=0$ (13.121). For ω ≪ f its first term
    drops and the relation of this block follows." + `nb.figure`: relative difference between the slow root of
    `GFD.shallow_water_omega` and `GFD.rossby_omega` against ω/f for a sweep of wavelengths. *expect:* 4.5 × 10⁻⁴ at ω/f =
    0.106 (external mode, 3100 km), 7 × 10⁻⁷ at ω/f = 8 × 10⁻⁴ (baroclinic) — the difference scales as (ω/f)².
15. `nb.primer("two-dimensional FFT wavenumber grids (np.fft.fftfreq, np.fft.fft2)", "A field on a periodic grid is a sum
    of plane waves. `np.fft.fft2` returns their complex amplitudes; `2π·np.fft.fftfreq(n, d)` gives the wavenumber that
    belongs to each array index (positive first, then negative). A linear equation with constant coefficients moves each
    plane wave independently — multiply its amplitude by $e^{-i\\omega t}$ and transform back: an exact solution with no time
    stepping.", code="k = 2*np.pi*np.fft.fftfreq(8, d=1.0)       # wavenumbers of an 8-point periodic grid\nprint(np.round(k,
    2))                         # [0, 0.79, 1.57, 2.36, -3.14, -2.36, -1.57, -0.79]")` (**P329**; reminder of P142)
16. `nb.animation` — **A4** (video, 90 frames, FAST 40): a wave packet from `SW.qg_linear_evolve_1d(a, kk, 0.0, t, beta35,
    f35, 3.6096)` with `kk, a = GFD.rossby_packet_spectrum(k0, sigma, 48)` (`> 🔧 **Our choice** — exact spectral evolution
    on a periodic line`): top, a short wave (150 km): crests move west at 0.82 cm/s while the envelope moves **east** at
    0.43 cm/s; bottom, a long wave (1000 km): crests west at 3.25 cm/s, envelope **west** at 2.81 cm/s. A row of marked
    columns coloured by relative vorticity (amber). *expect* (35° N, first baroclinic mode): the four speeds above; at
    2πΛ = 271 km the envelope stands still. *see / read / change*. + one cell with the two-dimensional
    `SW.qg_linear_evolve(eta0, x, y, t, beta35, f35, 3.6096)` on a Gaussian eddy: it drifts west and sheds a wake to its
    east (static 3-panel figure).
17. `nb.note` — **N138 [C]** "This relation fails within a few degrees of the equator, where geostrophy itself breaks
    down. Equatorial waves — the equatorial Kelvin wave and the equatorial radius √(c/β) — are beyond the book; pointer:
    any text on equatorial dynamics."
18. `nb.explainer("rossby_waves", heading="The crests go west — so how does the energy go east?", why="A packet whose
    crests march west while its envelope moves east cannot be drawn still; marked columns show each one's vorticity
    change as it is displaced; a mean flow slows the crests to a standstill.", tries=["Play the short-wave preset and
    follow one crest and the envelope separately.", "Drag the wavelength through 2πΛ and watch the status flip from
    'energy east' to 'energy west'.", "Add a mean flow until the crests stand still; read the wavelength.", "Click a
    column to see its displacement, βy and the vorticity it must have."])`

### A.16 §13.16 Barotropic Instability — C15 (continued): R30, N139–N142, D24
1. `nb.section("13.16", "Barotropic Instability", intro="**What is this section about?** Chapter 11's question — when is a
   shear flow unstable? — asked again on the β-plane. The answer is Rayleigh's criterion with one change: what must
   change sign is the gradient of *absolute* vorticity.")` → **(re-enter C15: `nb.current_core = "C15"`)**.
2. `nb.note` — **N139 [B]** constant depth turns potential-vorticity conservation into conservation of absolute vorticity:
   $\Big(\dfrac{\partial}{\partial t}+\mathbf u\cdot\nabla\Big)(\zeta+f)=0$ (13.122). `SW.barotropic_vorticity_rhs` is its spectral right-hand side (used in C17). **N140
   [B]** linearised about a zonal current U(y), with $u'=-\partial\psi/\partial y$, $v'=\partial\psi/\partial x$, $\bar\zeta=-dU/dy$:
   $\dfrac{\partial}{\partial t}(\nabla^2\psi)+U\dfrac{\partial}{\partial x}(\nabla^2\psi)+\Big(\beta-\dfrac{d^2U}{dy^2}\Big)\dfrac{\partial\psi}{\partial x}=0$ (13.123) (steps 1–3 of D24).
3. **(re-enter C15)** `nb.recap("R30", "Rayleigh's equation, now with β", "For normal modes $\\psi=\\hat\\psi(y)e^{ik(x-ct)}$:
   $(U-c)\\Big[\\dfrac{d^2}{dy^2}-k^2\\Big]\\hat\\psi+\\Big[\\beta-\\dfrac{d^2U}{dy^2}\\Big]\\hat\\psi=0$ — Chapter 11's Rayleigh equation with $-U''$ replaced by $\\beta-U''$. ⚠️ Trap T14: the
   stream function has the opposite sign to Chapter 11's; the eigenvalue c does not care.", where="Ch. 11,
   $(U-c)\\big(\\frac{d^2\\phi}{dy^2}-k^2\\phi\\big)-\\frac{d^2U}{dy^2}\\phi=0$ (11.81)")`
4. `nb.derivation("D24", …)` — Part F D24 (8 steps), ref "13.124": **N141 [B]** the Rayleigh–Kuo criterion: a necessary
   (not sufficient) condition for instability is that $\dfrac{d}{dy}(\bar\zeta+f)=\beta-\dfrac{d^2U}{dy^2}$ (13.124) changes sign somewhere in the flow
   (the book says the analysis carries over and never shows it). Reminder of P255 (necessary vs sufficient).
5. `nb.code` — our easterly jet at 12° N: `U = -14.0/np.cosh(y/4.0e5)**2`; `g1 = GFD.absolute_vorticity_gradient(y, U, beta12)`;
   `GFD.rayleigh_kuo_criterion(y, U, beta12)`; the same jet 1200 km wide; `ch11.rayleigh_criterion(y, U=U)` (the β = 0
   case); `ch13.rayleigh_kuo_eigs(1.4, U_fn, Up_fn, Upp_fn, beta_nd, bc="decay", y_max=16, parity="even")` for the
   non-dimensional **westerly** jet U = sech² y at ≤ 12 (k, β) points (cached; β is passed as its own argument and reaches
   `ST.rayleigh_eigs_contour(..., beta=…)` — never folded into `Upp`). *expect:* largest U″ of
   this jet is 2U₀/L² = 1.75 × 10⁻¹⁰ m⁻¹ s⁻¹ against β = 2.24 × 10⁻¹¹: the gradient changes sign → may be unstable; with L =
   1200 km (> √(2U₀/β) = 1118 km) it does not → stable by the criterion. A westerly jet of the same shape needs L < 646 km
   ("easterly jets are the easier ones to destabilise"). Growth rates (measured by the implementer; re-run by the
   designer with the options above): for U = sech² y at k = 1.4 the growth rate k c_i is 0.1191, 0.0744, 0.0151 at β =
   0, 0.3, 0.6 (in units of U₀/L and U₀/L²) and no growing mode is found at β = 0.7 — beyond the bound max U″ = 2/3 of this
   westerly jet, as the criterion requires. (With the default wall box instead of `bc="decay"` the β = 0.6 mode is not
   found: the builder keeps the stated options.) For scale: our 400 km, 14 m/s jet at 12° N has βL²/U₀ = 0.26. + hand-written sign test `np.any(np.diff(np.sign(beta12 - Upp)) != 0)`
   with `assert` against `GFD.rayleigh_kuo_criterion(...)["changes_sign"]`.
6. `nb.figure` — **N142 [B]** our Fig. 13.29: U(y), $\bar\zeta$, f and $\bar\zeta+f$ against latitude for the 400 km jet, the extremum of
   $\bar\zeta+f$ marked; right panel: growth rate against β for three wavenumbers. *see / read* ("where the dashed ζ̄ + f curve
   turns back, its gradient changes sign: waves there can exchange places and feed on the jet") */ change*.
7. `nb.md` — **What would change if…** "…the flow has no horizontal shear at all, only a vertical one — the thermal wind
   of C03? There is no inflection point and Ri is far above ¼, yet it is unstable: baroclinic instability (C16)."

---

### A.17 §13.17 Baroclinic Instability — C16
#### C16 — The Eady problem (13.141)
1. `nb.section("13.17", "Baroclinic Instability", intro="**What is this section about?** Why mid-latitude weather exists.
   A wind that increases with height over sloping density surfaces is unstable to waves a few thousand kilometres long;
   they grow by carrying warm air poleward and upward and cold air equatorward and downward. ⚠️ Here z = 0 is the lower
   lid and z = H the upper one; α is a scaled wavenumber, not the expansion coefficient; Λ = NH/f has no π.")`
2. `nb.core("C16", "The Eady problem: $c=\\dfrac{U_0}{2}\\pm\\dfrac{U_0}{\\alpha H}\\sqrt{\\Big(\\dfrac{\\alpha H}{2}-\\tanh\\dfrac{\\alpha H}{2}\\Big)\\Big(\\dfrac{\\alpha H}{2}-\\coth\\dfrac{\\alpha H}{2}\\Big)}$ (13.141)",
   question="The jet stream has no inflection point and a Richardson number far above ¼ — so where do storms come from?")`
3. `nb.md` — **Plain words:** "The tropics are warm, the poles cold, and in between the surfaces of constant density tilt
   — a little, about one part in five hundred. That tilt is stored energy: if cold air could slide down the slope under
   warm air sliding up it, the centre of mass would fall. Rotation forbids the direct slide (C03: the tilted state is in
   balance). But a wave of the right size can do it sideways, and once it starts, it feeds itself. The waves are the
   highs and lows of the weather map; their size, about four thousand kilometres, and their growth time, a day or two,
   come out of one formula. *The problem solved here is known as the Eady problem.*"
4. `nb.md` — **The idea:** two lids, and on each lid a wave that exists because temperature varies along it (an 'edge
   wave'). The lower one, left alone, drifts east slowly; the upper one sits in a fast eastward wind but itself
   propagates west relative to it. If they can feel each other through the layer, they lock: each one's flow strengthens
   the other. They can feel each other only if the wavelength is long compared with the layer's Rossby radius.
5. `nb.note` — **N143 [B]** the basic state (`bs = ch13.eady_basic_state(y, z, N=1.1e-2, f=f35, H=9000.0, U0=27.0,
   rho0=1.2)`): density surfaces sloping up toward the pole; by the thermal wind (C03) the eastward flow increases with
   height. + `nb.primer("available potential energy and the wedge of sloping convection", "Only part of a fluid's
   potential energy can ever be released: the part that would be freed by flattening the density surfaces. A parcel
   exchange releases energy only if the heavier parcel ends up lower *and* was the denser one at the same level: its
   path must lie inside the wedge between the horizontal and the sloping density surface. Steeper than the surface: it
   costs energy (ordinary static stability). Flatter than horizontal: impossible.", code="slope_rho = 2.07e-3
   # slope of our density surfaces (computed below)\nfor s in (0.5, 1.0, 1.5):                      # path slope / surface slope\n
   print(s, 'releases energy' if 0 < s < 1 else 'does not')")` (**P332**) + figure `draw_eady_wedge()`. *expect:* slope
   fU₀/(N²H) = 2.07 × 10⁻³ (about 1 in 480); Ri = N²H²/U₀² = 13.4.
6. `nb.note` — the starting set, named: **N144 [C]** the f-plane, hydrostatic, inviscid equations for the total flow
   (13.125): $\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$, $\dfrac{\partial v}{\partial t}+u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$, $0=-\dfrac{\partial p}{\partial z}-\rho g$,
   $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$, $\dfrac{\partial\rho}{\partial t}+u\dfrac{\partial\rho}{\partial x}+v\dfrac{\partial\rho}{\partial y}+w\dfrac{\partial\rho}{\partial z}=0$ (⚠️ trap T15: no $w\,\partial u/\partial z$ in the momentum equations); **N145 [C]**
   basic state plus perturbation (13.126): $u=U(z)+u'$, $v=v'$, $w=w'$, $\rho=\bar\rho(y,z)+\rho'$, $p=\bar p(y,z)+p'$; **N146 [C]** the basic state is
   geostrophic and hydrostatic (13.127): $fU=-\dfrac1{\rho_0}\dfrac{\partial\bar p}{\partial y}$, $0=-\dfrac{\partial\bar p}{\partial z}-\bar\rho g$.
7. `nb.primer("non-trivial solutions of a homogeneous 2 × 2 system: the determinant must vanish", "Two equations
   a·A + b·B = 0 and c·A + d·B = 0 always have the solution A = B = 0. They have another one only if the two equations
   say the same thing, i.e. the determinant ad − bc is zero. When the coefficients contain an unknown (a frequency, a
   phase speed), 'determinant = 0' is the equation that fixes it.", code="c = sp.symbols('c')\nM = sp.Matrix([[c, -1], [-4,
   c]])               # a toy system with an unknown c\nprint(sp.solve(M.det(), c))                  # [-2, 2]: only these c
   allow A, B != 0")` (**P330**) + `nb.primer("hyperbolic half-angle identities and coth", "coth x = cosh x / sinh x =
   1/tanh x: large near 0, → 1 for large x. Two identities turn products at the half argument into functions of the full
   argument: 2 sinh x cosh x = sinh 2x and cosh² x + sinh² x = cosh 2x; hence tanh x + coth x = 2 coth 2x.", code="x =
   0.8\nprint(np.tanh(x) + 1/np.tanh(x), 2/np.tanh(2*x))     # equal: 2.1696\nprint(brentq(lambda x: x - 1/np.tanh(x), 0.5,
   3))    # 1.19968: where x = coth x")` (**P331**; reminders of P168 cosh, sinh, tanh and P53 determinants)
8. `nb.derivation("D25", …)` — Part F D25 (13 steps), ref "13.136": **N147 [B]** its step 3, the thermal wind of the basic
   state $\dfrac{dU}{dz}=\dfrac g{f\rho_0}\dfrac{\partial\bar\rho}{\partial y}$ (13.128), so $U=\dfrac{U_0z}H$; **N148 [C]** step 4, $\dfrac{\partial\zeta}{\partial t}+u\dfrac{\partial\zeta}{\partial x}+v\dfrac{\partial\zeta}{\partial y}-(\zeta+f)\dfrac{\partial w}{\partial z}=0$ (13.129); **N149 [B]**
   step 5, $\dfrac{\partial\zeta'}{\partial t}+U\dfrac{\partial\zeta'}{\partial x}-f\dfrac{\partial w'}{\partial z}=0$ (13.130); **N150 [C]** step 6, $u'\simeq-\dfrac1{\rho_0f}\dfrac{\partial p'}{\partial y}$, $v'\simeq\dfrac1{\rho_0f}\dfrac{\partial p'}{\partial x}$ (13.131); **N151 [B]** step
   7, $\zeta'=\dfrac1{\rho_0f}\nabla_H^2p'$ (13.132); **N152 [B]** step 8, $\dfrac{\partial\rho'}{\partial t}+U\dfrac{\partial\rho'}{\partial x}+v'\dfrac{\partial\bar\rho}{\partial y}-\dfrac{\rho_0N^2w'}{g}=0$ (13.133); **N153 [C]** step 9,
   $0=-\dfrac{\partial p'}{\partial z}-\rho'g$ (13.134); **N154 [B]** step 10, $w'=-\dfrac1{\rho_0N^2}\Big[\Big(\dfrac\partial{\partial t}+U\dfrac\partial{\partial x}\Big)\dfrac{\partial p'}{\partial z}-\dfrac{dU}{dz}\dfrac{\partial p'}{\partial x}\Big]$ (13.135)
   (`GFD.eady_vertical_velocity`); **N155 [B]** the result, $\Big(\dfrac\partial{\partial t}+U\dfrac\partial{\partial x}\Big)\Big[\nabla_H^2p'+\dfrac{f^2}{N^2}\dfrac{\partial^2p'}{\partial z^2}\Big]=0$ (13.136): "the bracket is the
   perturbation potential vorticity; it is carried by the basic flow and is zero in the interior for a normal mode — so
   everything happens at the two lids". + `ch13.eady_qg_sympy()["ok"]` → True.
9. `nb.derivation("D26", …)` — Part F D26 (12 steps), ref "13.141": **N156 [C]** its step 1, $p'=\hat p(z)e^{i(kx+ly-\omega t)}$ (13.137);
   **N157 [B]** step 2, $\dfrac{d^2\hat p}{dz^2}-\alpha^2\hat p=0$ (13.138); **N158 [B]** step 3, $\alpha^2\equiv\dfrac{N^2}{f^2}(k^2+l^2)$ (13.139) (`GFD.eady_alpha`); **N159 [C]**
   step 4, $\hat p=A\cosh\alpha\Big(z-\dfrac H2\Big)+B\sinh\alpha\Big(z-\dfrac H2\Big)$ (13.140); **N160 [B]** steps 5–8: the two lid conditions
   $\dfrac{\partial^2p'}{\partial t\,\partial z}-\dfrac{U_0}H\dfrac{\partial p'}{\partial x}=0$ at z = 0 and $\dfrac{\partial^2p'}{\partial t\,\partial z}-\dfrac{U_0}H\dfrac{\partial p'}{\partial x}+U_0\dfrac{\partial^2p'}{\partial x\,\partial z}=0$ at z = H, and the 2 × 2 system
   $A\Big[\alpha c\sinh\dfrac{\alpha H}2-\dfrac{U_0}H\cosh\dfrac{\alpha H}2\Big]+B\Big[-\alpha c\cosh\dfrac{\alpha H}2+\dfrac{U_0}H\sinh\dfrac{\alpha H}2\Big]=0$,
   $A\Big[\alpha(U_0-c)\sinh\dfrac{\alpha H}2-\dfrac{U_0}H\cosh\dfrac{\alpha H}2\Big]+B\Big[\alpha(U_0-c)\cosh\dfrac{\alpha H}2-\dfrac{U_0}H\sinh\dfrac{\alpha H}2\Big]=0$ (`ch13.eady_matrix(c, alphaH, U0, H)`).
10. `nb.worked_example("storm size and growth time by hand", "f = 10⁻⁴ s⁻¹, N = 10⁻² s⁻¹, H = 10 km, U₀ = 30 m/s. 1. Eady radius
    $\\Lambda_E=NH/f$ = 10⁻² × 10⁴/10⁻⁴ = 10⁶ m = 1000 km. 2. Fastest wave (ours, computed): αH = 1.6061, so k = 1.6061/Λ_E and the
    wavelength is 2π/k = 3.9120 × 1000 km ≈ 3912 km. 3. Growth rate σ = 0.30982 × fU₀/(NH) = 0.30982 × 10⁻⁴ × 30/100 = 9.3 × 10⁻⁶
    s⁻¹. 4. e-folding time 1/σ = 1.08 × 10⁵ s ≈ 1.2 days. 5. Shortest unstable wave: 2.6187 × 1000 km ≈ 2619 km; anything
    shorter is two neutral edge waves.")`
11. `nb.derivation("D27", …)` — Part F D27 (6 steps), ref "13.142": **N161 [B]** the marginal condition
    $\dfrac{\alpha_cH}2=\coth\Big(\dfrac{\alpha_cH}2\Big)$; instability for $\alpha H<\alpha_cH$, i.e. $\dfrac{HN}f<\dfrac{\alpha_cH}{\sqrt{k^2+l^2}}$, and for l = 0 $\dfrac{HN}f<\dfrac{\alpha_cH}k$ (13.142; the book
    prints $\alpha_cH$ as a rounded number); with $\Lambda\equiv\dfrac{HN}f$ the unstable wavelengths are $\lambda>(2\pi/\alpha_cH)\,\Lambda$. `> ⚠️ Common confusion
    (trap T11):` "this Λ has no π". **N162 [B]** **Not in the book — ours** (the book stops at the wavelength): the
    maximum growth rate $\sigma_{max}=0.30982\,\dfrac fN\dfrac{dU}{dz}$ at αH = 1.6061 (l = 0) and its e-folding time — ours, computed, no
    citation.
12. `nb.code` — `GFD.eady_critical()`; `GFD.eady_fastest()`; `GFD.eady_factors(1.6061)`; `c = GFD.eady_phase_speed(1.0, 27.0)`;
    `GFD.eady_phase_speed(3.0, 27.0)`; `GFD.rossby_radius_internal(1.1e-2, 9000.0, f35, with_pi=False)`;
    `GFD.eady_max_growth_rate(f35, 1.1e-2, 27.0/9000.0)`; `GFD.eady_time_scale(f35, 1.1e-2, 27.0/9000.0)/86400`;
    `GFD.eady_growth_rate(1.6061/Lam, 0.0, 1.1e-2, f35, 9000.0, 27.0)`; the ocean twin (N = 5 × 10⁻³ s⁻¹, H = 1000 m, U₀ = 0.1
    m/s). *expect:* α_cH = 2.3994 (to five digits: 2.39936); fastest αH = 1.6061, σ NH/(fU₀) = 0.30982, c_i/U₀ = 0.1929; the
    two factors at the fastest wave +0.1373 and −0.6990; c = (0.5 + 0.2511i) U₀ at αH = 1 → 13.5 + 6.78i m/s; at αH = 3 two
    real speeds, the faster 0.6616 U₀; atmosphere (35° N, N = 1.1 × 10⁻² s⁻¹, H = 9 km, U₀ = 27 m/s): Λ_E = 1183 km, σ_max =
    7.068 × 10⁻⁶ s⁻¹, e-folding 1.637 days, fastest wavelength 3.9120 Λ_E = 4630 km, cut-off 2.6187 Λ_E = 3099 km; ocean
    twin: Λ_E = 59.8 km, e-folding 22.3 days, fastest wavelength 234 km. *explain:* 6 points.
13. `nb.check_agree` — **from scratch:** at αH = 1.4: `M = lambda c: ch13.eady_matrix(c, 1.4, 27.0, 9000.0)`; the determinant
    is quadratic in c, so three evaluations (c = 0, 1, 2) give its coefficients and `np.roots` its zeros; `assert
    np.isclose(roots.max_imag, GFD.eady_phase_speed(1.4, 27.0))` *expect:* c = (0.5 + 0.2158i) U₀; and
    `ch13.eady_numeric_eigs(k, 0.0, 1.1e-2, f35, 9000.0, 27.0)["c"]` (the Chebyshev route) agrees to 10⁻⁶.
14. `nb.figure` — our Fig. 13.31: upper, coth x, tanh x and the line x against αH (x = αH/2), the crossing x = coth x at
    αH = 2.3994 marked; lower, the growth rate $kc_i$ in units of fU₀/(NH) against αH with its maximum 0.30982 at 1.6061.
    `> ⚠️ Common confusion:` "the growth rate is $kc_i$, not $c_i$ — $c_i$ is largest as k → 0, where nothing grows."
    *see / read / change* ("…l ≠ 0: the same curve read at α = (N/f)√(k² + l²); the growth rate carries the factor k/K and
    is smaller").
15. `nb.plotly` — **F9** (`slider_figure`, slider = U₀ from 5 to 50 m/s, 19 steps; dropdown atmosphere / ocean): growth rate
    in day⁻¹ against wavelength in km; trace name carries the e-folding time. + `nb.live`: free N, H, U₀, latitude.
16. `nb.note` — **N163 [B]** energetics (stated; the integrations are Chapter 11's energy-equation moves, whose result
    there was $\frac d{dt}\int\frac12u_i^2\,dV=-\int u_iu_j\frac{\partial U_i}{\partial x_j}\,dV-\Lambda$ (11.88) with Λ the dissipation): here the perturbation kinetic energy
    $K_E\equiv\dfrac{\rho_0}2\displaystyle\int(u'^2+v'^2)\,dx\,dy\,dz$ grows through the buoyancy flux alone, $\dfrac{dK_E}{dt}=-g\displaystyle\int w'\rho'\,dx\,dy\,dz$ — not through the
    Reynolds stress against the shear. `fl = GFD.eady_fluxes(z, 1.6061, 27.0, 9000.0, 1.1e-2, f35, 1.2)`; `md =
    GFD.eady_mode(z, 1.6061, 27.0, 9000.0)`: figure of `fl["w_rho"]` (negative at every height: light fluid rises, heavy
    sinks), `fl["v_rho"]` (heat goes poleward) and the phase of p̂ against height. *expect:* `fl["w_rho"] <= 0` everywhere;
    `md["phase"]` increases monotonically with height (the pressure pattern tilts westward with height); the total tilt
    reported by the builder.
17. `nb.animation` — **A5** (video, 80 frames, FAST 40): the fastest mode in a longitude–height section
    (`GFD.eady_mode`), pressure contours tilting westward with height and growing by e every 1.64 days (amplitude
    renormalised each e-folding, a clock in the title), temperature anomalies on the two lids, the basic-state wind
    arrows at the side. *see / read / change* ("…αH = 3: no tilt, no growth; two edge waves pass each other").
18. `nb.explainer("eady_instability", heading="Where do storms come from?", why="Sliding the wavenumber through the cut-off
    turns a growing, westward-tilting mode into two separate neutral edge waves; changing the shear, the stratification
    or the latitude rescales the growth curve and prints the e-folding time in days and the preferred wavelength in
    kilometres.", tries=["Use the preset 'fastest-growing wave' and read the e-folding time.", "Drag the wavelength
    shorter until the status changes to 'two neutral edge waves'.", "Double N: what happens to the growth rate and to the
    preferred wavelength?", "Switch to the ocean mode: mesoscale eddies instead of weather systems."])`
19. `nb.md` — **What would change if…** "…the growing eddies become strong enough to interact with each other? Then we
    are in turbulence — but of an unusual, nearly two-dimensional kind (C17)."

---

### A.18 §13.18 Geostrophic Turbulence — C17; S01, S02; summary
#### C17 — Fjørtoft's argument (13.145): energy to large scales, enstrophy to small
1. `nb.section("13.18", "Geostrophic Turbulence", intro="**What is this section about?** What turbulence does when
   rotation and stratification keep it nearly two-dimensional. Chapter 12's cascade runs backwards: energy moves to
   *larger* scales. ⚠️ In this section α is the enstrophy flux, η (once) the Kolmogorov length, l an eddy size.")`
2. `nb.core("C17", "Fjørtoft's argument: with $S_0=S_1+S_2$ and $K_0^2S_0=K_1^2S_1+K_2^2S_2$, $\\dfrac{S_1}{S_2}=\\dfrac{K_2-K_0}{K_0-K_1}\\,\\dfrac{K_2+K_0}{K_1+K_0}$ and
   $\\dfrac{K_1^2S_1}{K_2^2S_2}=\\dfrac{K_1^2}{K_2^2}\\,\\dfrac{K_2^2-K_0^2}{K_0^2-K_1^2}$ (13.145)", question="Why do the eddies of the atmosphere and ocean merge into bigger
   eddies and jets instead of breaking down into smaller ones?")`
3. `nb.md` — **Plain words:** "Stir a cup of coffee and the swirls break into smaller swirls until viscosity erases them
   (Ch. 12). Jupiter's atmosphere does the opposite: small storms merge into large ones and into bands that have lasted
   centuries. So do ocean eddies. The difference is one extra conservation law. In flat, two-dimensional flow a vortex
   cannot be stretched, so the total squared vorticity — the enstrophy — is conserved along with the energy. Two
   conserved quantities are too many for a simple downhill cascade."
4. `nb.md` — **The idea** (a see-saw): energy and enstrophy sit on the same wavenumbers but enstrophy weighs them by K².
   Move some energy to higher K: enstrophy goes up; to keep it fixed you must move *more* energy to lower K. Sketch
   `draw_cascade_arrows()`: energy ←, enstrophy →.
5. **(re-enter C17)** `nb.recap("R31", "Where the energy finally goes", "Energy must still reach the Kolmogorov scale to be
   dissipated; three-dimensional turbulence does that through the cascade with spectrum $S_{11}=C_1\\bar\\varepsilon^{2/3}k_1^{-5/3}$ (12.54, in
   its corrected form). Internal waves are a suggested route from the large eddies to it.", where="Ch. 12 §12.7 (C07–C08)")`
6. `nb.primer("enstrophy: mean-square vorticity, and its spectrum K²S(K)", "Enstrophy is to vorticity what kinetic
   energy is to velocity: its mean square. A Fourier mode of wavenumber K with velocity amplitude û has vorticity
   amplitude Kû (a derivative multiplies by K), so its share of the enstrophy is K² times its share of the energy.",
   code="K = np.array([1.0, 2.0, 4.0]); S = np.array([4.0, 2.0, 1.0])   # energy in three modes\nprint(S.sum(), (K**2*S).sum(),
   K**2*S)                    # energy 7; enstrophy 28, most of it in the SMALLEST scale")` (**P333**)
7. `nb.note` — **N164 [B]** geostrophic turbulence: nearly two-dimensional (rotation, stratification and thinness
   suppress w), no vortex stretching; with the isotropic spectra $\overline{u^2}=\displaystyle\int_0^\infty S(K)\,dK$ and $\overline{\zeta^2}=\displaystyle\int_0^\infty K^2S(K)\,dK$.
   `> ⚠️ Common confusion (trap T16):` "This S is one-sided in K with no factor ½; Chapter 12's spectra were two-sided
   with ∫S = the variance of one component." `GFD.enstrophy_spectrum(K, S)`; `TS.shell_spectrum` on a synthetic field.
   **N165 [B]** the two conservation statements $\dfrac d{dt}\displaystyle\int_0^\infty S(K)\,dK=0$ and $\dfrac d{dt}\displaystyle\int_0^\infty K^2S(K)\,dK=0$ (13.143)–(13.144) (inviscid,
   two-dimensional). `SW.barotropic_invariants(zeta, L)`.
8. `nb.derivation("D28", …)` — Part F D28 (5 steps), ref "13.145".
9. `nb.worked_example("where the energy of wavenumber 6 goes", "Energy S₀ = 1 at K₀ = 6 is moved to K₁ = 2 and K₂ = 9 (K₁ =
   K₀/3, K₂ = 3K₀/2). 1. Energy: S₁ + S₂ = 1. 2. Enstrophy: 4S₁ + 81S₂ = 36. 3. Subtract 4 × (1) from (2): 77S₂ = 32 → S₂ =
   0.416, S₁ = 0.584. 4. Energy ratio S₁/S₂ = 1.41: more energy went to the **larger** scale. 5. Enstrophy shares: 4 × 0.584 =
   2.34 against 81 × 0.416 = 33.66, ratio 0.069: fourteen times more enstrophy went to the **smaller** scale. 6. Formula
   check: $\\frac{K_2-K_0}{K_0-K_1}\\frac{K_2+K_0}{K_1+K_0}$ = (3/4)(15/8) = 1.406 ✓.")`
10. `nb.code` — `r = GFD.fjortoft_transfer(1.0, 1/3, 1.5)`; a sweep of K₂/K₀ for fixed K₁. *expect:* S₁ = 0.5844, S₂ = 0.4156,
    energy ratio 1.40625, enstrophy ratio 0.06944 (= 1/14.4). *explain:* 3 points.
11. `nb.check_agree` — **from scratch:** `A = np.array([[1.0, 1.0], [K1**2, K2**2]])`; `S1, S2 = np.linalg.solve(A, [S0,
    K0**2*S0])`; `assert np.allclose([S1, S2], [r["S1"], r["S2"]])`.
12. `nb.figure` — energy centroid before and after for many triads (K₁, K₂): a map of the energy ratio S₁/S₂ over the
    (K₁/K₀, K₂/K₀) plane with the line S₁ = S₂. *see / read* ("almost everywhere above 1") */ change*.
13. `nb.derivation("D29", …)` — Part F D29 (7 steps), ref "": **N166 [B]** the two inertial ranges, $S(K)\propto\varepsilon^{2/3}K^{-5/3}$
    (energy flux ε toward small K) and $S(K)\propto\alpha^{2/3}K^{-3}$ (enstrophy flux α toward large K) — the exponents by
    dimensional analysis (`DIM.pi_groups`; reminder of Chapter 1's Π theorem); **N167 [B]** the Rhines length
    $l\sim\sqrt{u/\beta}$. Both presented as **ours, computed** — no citation. `GFD.two_d_cascade_spectrum(K, K0, eps,
    alpha_ens)` (shape only), `GFD.rhines_length(13.0, beta35)`, `GFD.rhines_length(0.08, beta35)` *expect:* 833 km
    (atmosphere, u_rms = 13 m/s) and 65 km (ocean, 0.08 m/s) at 35° N. + figure: our Fig. 13.33 (log S against log K,
    injection at K₀, the two slopes, arrows for the two fluxes).
14. `nb.primer("the Jacobian J(ψ, ζ) and a pseudo-spectral step with 2/3 de-aliasing", "The advection of vorticity by
    the flow it induces is u·∇ζ = ψ_xζ_y − ψ_yζ_x ≡ J(ψ, ζ). A pseudo-spectral model takes derivatives in Fourier space
    (multiply by ik — exact), multiplies fields on the grid (cheap), and zeroes the top third of the wavenumbers before
    each product so that the product's new short waves do not fold back onto long ones (aliasing). This model is OUR
    choice; the book describes such calculations only in words.", code="n = 8; k = np.fft.fftfreq(n, 1/n)                  #
    integer wavenumbers -4..3\nprint(np.abs(k) <= n/3)                              # the 2/3 rule keeps |k| <= 2")` (**P334**)
15. `nb.animation` — **A6** (frames, log-spaced times; 2 rows) — **N168 [B]** **Not in the book — ours**, a demo labelled
    *qualitative*: decaying two-dimensional turbulence from the caches `ch13.load_reference_run("turbulence_f")` and
    `("turbulence_beta")` (**64²** pseudo-spectral runs, 10 saved frames each, of $\dfrac{\partial\zeta}{\partial t}+J(\psi,\zeta)+\beta\dfrac{\partial\psi}{\partial x}=\nu\nabla^2\zeta$, $\zeta=\nabla^2\psi$; `> 🔧 **Our choice** — our numerical
    model, cached`): vorticity field (vortices merging) beside the energy and enstrophy spectra and the two invariants
    against time; second row with β: the field bands into zonal jets. If a cache is absent and `not FAST`, a 64² run
    `SW.barotropic_run(64, L, zeta0, t_end=…, dt=…, nu=…)` is made live (≈ 4 s); in FAST the cell prints one sentence.
    *expect* (read by the designer from the two caches): energy falls to 0.73 of its initial value while enstrophy falls
    to 0.21 (0.22 with β); the energy centroid K_E (`SW.spectral_centroids`) moves from 8.6 to 3.2 (3.6 with β) — energy moves to larger scales; the
    zonal (k_x = 0) share of the energy at the last frame is 0.11 without β and 0.29 with β. **Every statement about this
    run is qualitative: at 64² the K⁻³ range is not resolved, and no spectral slope is read off or quoted from it** — the
    exponents come from the dimensional argument of D29 only. No jet-spacing number is claimed either. *see / read /
    change*.
16. `nb.pointer` — **S01** "Exercises 13.1–13.8 are not reproduced. The one result the text relies on — the group velocity
    of inertia–gravity waves — is derivation D21 in C14." · **S02** "Literature and supplemental reading: see the book's
    list; this notebook cites only what was read (the book itself)."
17. `nb.md` — **What would change if…** "…the fluid is a gas moving so fast that its own compressibility matters, or a
    wing rather than a planet? Those are Chapters 14 and 15. For climate dynamics, the next steps beyond this book are
    the equatorial waves named in C15 and the wind-driven gyre named in C05."
18. `nb.summary(clicked=[…17…], feeds_forward=[…], left_out=[…])` — **clicked** (one sentence per CORE the reader could now
    explain to a friend): C01 "On a thin shell only the vertical part of the earth's spin, f = 2Ω sin θ, turns the wind."
    · C02 "Slow flow runs along the isobars because that is the only direction in which Coriolis can cancel the pressure
    force." · C03 "A horizontal temperature gradient sets how the wind changes with height." · C04 "Friction against
    rotation makes a layer of fixed thickness in which the current spirals." · C05 "The wind-driven transport is at right
    angles to the wind and does not depend on the eddy viscosity; its divergence pumps water up or down." · C06 "Near the
    ground friction lets the wind cross isobars toward low pressure: air rises in lows." · C07 "Three equations for
    surface height and two velocities hold all the dynamics of a thin layer." · C08 "A stratified ocean is a stack of
    such layers; the first internal one is about a metre deep." · C09 "One cubic: two fast gravity waves and one slow
    wave that needs β." · C10 "Rotation puts a floor f under gravity-wave frequencies and turns currents into ellipses."
    · C11 "A coast replaces the missing cross-shore flow with a slope; the wave runs with the coast on its right (f >
    0)." · C12 "Beyond a Rossby radius c/f rotation holds a slope up against gravity." · C13 "Each column keeps (ζ + f)/h."
    · C14 "Internal waves live between f and N; the slope of the crests sets the frequency." · C15 "Rossby crests drift
    west because displaced columns must change their spin; energy goes east for short waves." · C16 "Storms grow by
    sliding fluid along the wedge between level and density surfaces, at a rate 0.30982 (f/N) dU/dz." · C17 "Two
    conserved quantities send energy to large scales and enstrophy to small." **feeds_forward**: Ch. 15 (the
    shallow-water ↔ gas-dynamics analogy uses √(gH) as the 'sound speed'); the user's climate-dynamics work (`GFD`, `VM`
    are importable on their own). **left_out**: equatorial waves; the wind-driven gyre and western boundary currents; the
    moist adiabat; pressure coordinates; quasi-geostrophic theory in a continuously stratified fluid beyond the Eady
    problem; the Charney problem; observed data figures.

### A.19 Placement check (every curation id has exactly one home)
- **CORE (17):** C01 A.4 · C02, C03 A.5 · C04, C05 A.6 · C06 A.7 · C07 A.8 · C08 A.9 · C09 A.10 · C10 A.11 · C11, C12 A.12 ·
  C13 A.13 · C14 A.14 · C15 A.15–A.16 · C16 A.17 · C17 A.18.
- **RECAP (31):** R01, R02, R03 A.2 · R04–R10 A.3 · R11, R12 A.4 (C01) · R13, R14 C02 · R15 C03 · R16 C04 · R17, R18 C07 ·
  R19, R20, R21, R22 C08 · R23 C11 · R24 C12 · R25, R26, R27 C13 · R28, R29 C14 · R30 A.16 (C15) · R31 C17.
- **NOTE (168):** N01, N02 A.1 · N03–N07 A.2 · N08–N11 A.3 · N12–N17 C01 · N18, N22 C02 (N19's verbal part is restated in
  C02 row 13; its home is C12 row 8) · N20, N21, N23, N24, N25, N26 C03 · N27, N28, N29, N30, N31, N33, N35 C04 · N32, N36,
  N37, N50 C05 · N34, N38–N49 C06 · N51–N54 C07 · N55–N77 C08 · N78–N83 C09 · N84–N90 C10 · N91–N96 C11 · N19 C12 ·
  N97–N104 C13 · N105–N126 C14 · N127–N138 C15 (A.15) · N139–N142 C15 (A.16) · N143–N163 C16 · N164–N168 C17.
- **SKIP (2):** S01, S02 A.18 row 16.
- **Derivations (29 `nb.derivation` rows):** D01, D02 C01 · D03 C02 · D04, D05 C03 · D06 C04 · D07, D09 C05 · D08 C06 · D10
  C07 · D11 C08 · D12, D13 C09 · D14 C10 · D15 C11 · D16 C12 · D17, D18 C13 · D19, D20, D21 C14 · D22, D23, D24 C15 · D25,
  D26, D27 C16 · D28, D29 C17.
- **Primers (29):** P307 A.1 · P308 C01 · P309 C02 · P310, P311 C04 · P312, P313 C05 · P335 C07 · P314, P315, P316, P317,
  P318 C08 · P319, P320 C09 · P321, P322 C10 · P323 C11 · P324 C12 · P325 C13 · P326, P327 C14 · P328, P329 C15 · P330,
  P331, P332 C16 · P333, P334 C17.
- **Animations:** A1 C10 · A2 C12 · A3 C11 · A4 C15 · A5 C16 · A6 C17. **Slider figures:** F1 C01 · F2 A.2 (§13.2; the
  two-convention lapse-rate figure — it precedes `nb.core("C01")`, so C01's own visuals are the 3-D tangent plane, the term
  bars and F1) · F3 C03 · F4 C04 · F5 C06 · F6 C08 · F7 C09 · F8 C15 · F9 C16. **Live cells:** C05 (with F4), C16 (with F9).
  **Explainers (10, each embedded once):** E1 C02 · E2 C03 · E3 C05 · E4 C06 · E5 C08 · E6 C10 · E7 C11 · E8 C12 · E9 C15 ·
  E10 C16. **Blocks with no explainer and their visuals:** C01 (plotly 3-D, term bars, F1) · C04 (figure, plotly 3-D, F4;
  its explainer E3 is embedded in C05) · C07 (x–t figure) · C09 (term-bar figure, F7; its explainer E6 is embedded in C10)
  · C13 (two figures) · C14 (four figures, plotly helix) · C17 (A6, two figures).
- **From-scratch `nb.check_agree` rows:** one in each of C01–C17 (17).


---

## Part B — explainer storyboards

Common to all ten (and the backup): created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch13">`; tabs
Walkthrough · Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function
that mirrors a `ch13` callable and is proved by `selftest()` parity rows (`py:` expressions use only `ch13.…`, `np.…`,
`math.…`, numbers, strings, dict keys and integer indices — C.6; **never a sampled random value and never a typed constant the
page could compute**). Explain is "Explanation & interpretation" in numbered sections built with `Viz.work.step / line / box
/ table / hint / interpret`, modelled on `forced_damped_vibrations.html` (the regime-dependent reading) and
`amplitude_phase_second_order_II_3.html` (a numbered derivation with live numbers, one section per optional view): **0** what
the views show and what each colour means · **1…n** every displayed quantity from the controls ("formula = substituted =
result — why", results boxed) · a section or hint per view hidden on phones · the values at the current time (live) ·
**Reading the current setting** (regime-dependent, with the threshold of each regime). Derivation steps are copied from Part
F (same `did` titles, same step count, same order; phones shorten *why* to its first sentence; plain-text *why* and *watch*
never contain raw TeX). Every tour, Explain, Derivation, notes, status, equation and quiz text that names a book equation
**writes it out** next to its number. Colours as header convention 10 (pressure-gradient orange, Coriolis teal, friction
rose, acceleration purple, density / temperature blue, vorticity amber, ghosts muted). Walkthrough texts ≤ 45 words, step 1
≤ 24 words, ≤ 2 extras per step, `play: false` on steps that quote numbers. **Hygiene (chapter 12 review):** every power of a
formatted number is bracketed (`(${Viz.tnum(x)})^{2}`, never scientific notation followed by a bare ^2); "the book prints",
never "the page prints"; slips in the house form `⚠️ slip #k — the book prints … ; the correct form is …`; no library status
string is shown raw; each control's extreme values are named in its help line, and where a formula stops being valid the
curve ends with a mark and a label (no NaN text anywhere). **Hemisphere:** E1, E2, E3, E4, E7, E8 carry `hemi` chips (N / S =
sign of f); every "right / left / clockwise" string is produced by one function `side(s)` that mirrors `ch13.hemisphere`.
**Fit:** at 360×640 never three stacked views — two at a time, the third is `hidePortrait` and its key number is repeated in
a visible title or readout; at 844×345 (the short landscape frame inside the chapter page) the stage is one row of two views
beside the text panel, view titles shortened, transport without step buttons. Parallel builders use private scratch
subfolders (`<scratchpad>/<slug>/`). JS constants allowed as definitions: `OMEGA = 7.292115e-5`, `R_EARTH = 6.371e6`, `G =
9.80665`.

### E1 · geostrophic_balance
- **Title:** "Why doesn't air flow straight from high to low?" · **Summary:** "Release a parcel on a pressure map and watch the
  Coriolis force turn it until it runs along the isobars." · **CORE:** C02 (also R13, N18, R14, N22 and the term names of C01) ·
  **Reference:** `forced_damped_vibrations.html` (system + graph on one time slider; the Explain panel).
- **meta:** `viz:order 1` · `viz:sections 13.5` · `viz:equations 13.8 13.11 13.12 13.13 13.18` · `viz:fluidpy
  ch13.coriolis_parameter ch13.geostrophic_velocity ch13.rossby_number ch13.ekman_number ch13.parcel_adjust ch13.pressure_centre
  ch13.hemisphere` · `viz:derivations D03`.
- **Physics (JS ↔ Python):** `fCor(latDeg, hemi)` ↔ `ch13.coriolis_parameter(lat_rad)`; `geoWind(dpdx, dpdy, f, rho0)` → [u, v] ↔
  `ch13.geostrophic_velocity(dpdx, dpdy, f, rho0)`; `rossby(U, f, L)` ↔ `ch13.rossby_number`; `ekmanNo(nu, f, L)` ↔
  `ch13.ekman_number`; `parcel(t, Gx, Gy, f, r, u0, v0)` → [u, v] (closed form $V=V_\infty+(V_0-V_\infty)e^{-(r+if)t}$, $V_\infty=-G/(r+if)$) ↔
  `ch13.parcel_adjust(t, G, f, r, V0)`; `pField(x, y, dp, Lc, kind)` (Gaussian anomaly, analytic gradient) ↔
  `ch13.pressure_centre(...)["p"]`; the parcel on the curved map is advanced with `Viz.num.rk4Step` on
  $\dot u=fv-\frac1{\rho_0}p_x-ru$, $\dot v=-fu-\frac1{\rho_0}p_y-rv$ (the same right-hand side; no parity — the uniform-gradient closed form carries it).
- **Views** (rows [1.25, 1]): 1. `map` "Pressure map" (row 0, flex 1.3, `equal`): 3000 km square; isobars every 1 hPa (orange,
  thin), the centre marked L or H, geostrophic wind arrows (muted teal quiver), one released parcel (dot) with its trail —
  bold so far, the geostrophic circle through its start as a faint ghost. Pointer: click to release the parcel there.
  2. `forces` "Forces on the parcel" (row 0, flex 1): arrows from the parcel: pressure-gradient force (orange), Coriolis force
  (teal), drag (rose, if on), their sum = acceleration (purple); the velocity as a grey arrow; scale bar in 10⁻⁴ m/s². 3. `time`
  "Speed and angle to the isobars" (row 1, `hidePortrait`): speed (purple) and cross-isobar angle (rose) against t, the
  geostrophic speed as a dashed ghost, bold up to the transport time.
- **Controls:** `centre` chips "low / high" (default low) · `dp` "Pressure difference $\Delta p$" (1…20 hPa, step 0.5, default 4;
  help: "1 hPa: a weak ripple · 20 hPa: a deep storm") · `Lc` "Size of the centre" (300…1500 km, default 600) · `lat` "Latitude"
  (3…80°, step 1, default 35; help: "3°: f is tiny, Ro large · 80°: f near its maximum") with `hemi` chips N / S · `drag`
  "Drag $r/\lvert f\rvert$" (0…1, default 0; optional; help: "0: no friction · 1: friction as strong as Coriolis") · transport `t`.
- **Transport:** `t` 0…72 h, rate 6 h/s, `end: 'hold'`; end card: "after 72 h: mean speed 5.7 m/s = the geostrophic speed · the
  parcel circled the low 0.4 times, low on its left".
- **Presets:** "mid-latitude low" {lat 35, N, low, dp 4, drag 0} · "same low, 35° S" {hemi S} · "near the equator" {lat 5} · "with
  surface drag" {drag 0.4} · "a high" {centre high}.
- **Status** (built by the page, three parts — verdict · numbers · side): Ro < 0.3: "✅ geostrophic · Ro = 0.11 · wind along
  isobars, low on the left" (f > 0) / "… low on the right" (f < 0); 0.3 ≤ Ro < 1: "≈ nearly geostrophic · Ro = 0.5 · the
  acceleration bar is visible"; Ro ≥ 1: "⚠️ Ro = 1.9 — not geostrophic: the parcel cuts across the isobars"; drag > 0 adds "·
  friction turns the wind 22° toward low pressure". Ro uses the largest geostrophic wind of the centre,
  $U_g=0.858\,\Delta p/(\rho_0\lvert f\rvert L_c)$ (the factor $\sqrt2e^{-1/2}$ is computed in the page).
- **Readouts:** "Coriolis $f$" · "Geostrophic wind" · "Rossby number" · "Speed now" · "Angle to isobars".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **terms**
  ("x-components" and "y-components": pressure + Coriolis + drag − acceleration = 0, four bars each; click a bar to
  isolate its arrow).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**Map:** orange lines join points of equal pressure; the dot is one parcel of air, the line
     behind it its path. **Forces:** <b class=c-orange>orange</b> = pressure-gradient force (toward low pressure); <b
     class=c-teal>teal</b> = Coriolis force (at right angles to the motion, to its right for f > 0, to its left for f < 0);
     <b class=c-rose>rose</b> = drag; <b class=c-accent>purple</b> = what is left over = the acceleration."
  1. *The Coriolis parameter at your latitude* — "$f=2\Omega\sin\theta$ (13.8) $=2\times7.292\times10^{-5}\times\sin35°=$ **8.37 × 10⁻⁵ s⁻¹** (negative in
     the southern hemisphere). Inertial period 2π/abs(f) = 20.9 h: the time the parcel needs for one loop."
  2. *The pressure force where the parcel is* — "$\frac1{\rho_0}\lvert\nabla p\rvert$ = (0.858 × 400 Pa / 600 km)/1.2 kg m⁻³ = **4.8 × 10⁻⁴ m/s²** at
     the steepest ring of the centre."
  3. *The wind that balances it* — "$-fv=-\frac1{\rho_0}\frac{\partial p}{\partial x}$, $fu=-\frac1{\rho_0}\frac{\partial p}{\partial y}$ (13.11)–(13.12): speed = force/abs(f) = 4.8 × 10⁻⁴/8.37 × 10⁻⁵ =
     **5.7 m/s**, along the isobars."
  4. *Is the balance a good approximation?* — "$\mathrm{Ro}=\frac{U^2/L}{fU}=\frac U{fL}$ (13.13) = 5.7/(8.37 × 10⁻⁵ × 6 × 10⁵) = **0.11**: the acceleration is
     about a tenth of the Coriolis force."
  5. *With drag* — "steady wind = $-G/(r+if)$: turned by arctan(r/abs(f)) = arctan 0.4 = **21.8°** toward low pressure and
     slowed to cos(21.8°) = 0.93 of the geostrophic speed. (Friction over Coriolis is the Ekman number,
     $E=\frac{\rho\nu U/L^2}{\rho fU}=\frac\nu{fL^2}$ (13.18).)"
  6. *The time view* (or the hint "turn the phone sideways for the time series") — "the speed oscillates about the
     geostrophic value with the inertial period; without drag the oscillation never dies — the parcel loops for ever
     about the balanced path."
  7. *At the current time* — "now t = `Viz.live('t')` h: speed `Viz.live('spd')` m/s, angle to the isobars `Viz.live('ang')`°,
     leftover acceleration `Viz.live('acc')` × 10⁻⁴ m/s²."
  8. *Reading the current setting* — Ro < 0.3: "The two big arrows nearly cancel. The parcel's loops are small wobbles on
     a path along the isobars: this is the geostrophic wind of the weather map, low pressure on its left (northern
     hemisphere; on its right in the southern)." · 0.3 ≤ Ro < 1: "The balance is only approximate; the purple leftover is a
     sizeable fraction of the others. Tight, fast systems — hurricanes — live here and need the curvature term the book
     leaves out." · Ro ≥ 1: "The Coriolis force is too weak (low latitude) or the system too small. The parcel accelerates
     toward low pressure almost as it would on a non-rotating earth: near the equator the wind does blow across the
     isobars." · drag > 0: "Friction slows the parcel, the Coriolis force weakens with it, and the pressure force wins a
     little: the wind crosses the isobars toward low pressure. That is the surface wind spiralling into a low (C06)."
- **Derivation tab:** **D03** (6 steps) `view: 'forces'` on wide screens, `'map'` on phones; goal `set` {preset mid-latitude
  low, t: 0}; step 2 `live` "U²/L = 5.7 × 5.7/600 000 = 5.4 × 10⁻⁵ against fU = 4.8 × 10⁻⁴ m/s²"; step 4 `set` {t: 72}, `watch` "the
  orange and teal arrows end up equal and opposite"; step 5 `live` "u·∇p = `Viz.live('udotgradp')` (zero when balanced)";
  step 6 `watch` "the path is an isobar: pressure is the stream function". Interpret: `s => "With Ro = " + … + " the
  acceleration is " + … + " % of the Coriolis force."`
- **Code:**
  ```python
  f = ch13.coriolis_parameter(np.deg2rad({{lat}}))        # (13.8): f = {{f}} 1/s
  u, v = ch13.geostrophic_velocity({{dpdx}}, {{dpdy}}, f, 1.2)   # (13.11)-(13.12): ({{ug}}, {{vg}}) m/s
  Ro = ch13.rossby_number({{U}}, f, {{L}})                # (13.13): Ro = {{Ro}}
  G = ({{dpdx}} + 1j*{{dpdy}})/1.2                        # pressure-gradient acceleration, complex
  V = ch13.parcel_adjust({{t}}*3600, G, f, r={{r}})        # parcel velocity now: {{V}} m/s
  ```
- **Walkthrough (6 steps):** 1. "High to low?" — "A ball rolls downhill. Press ▶ and see whether air does the same on a
  pressure map." `play: true` · 2. "Two forces" — "Orange pushes toward low pressure. Teal, the Coriolis force, is always at
  right angles to the motion and grows with speed." `set` {t: 2}, `terms: true` · 3. "The stand-off" — "After a few hours teal
  cancels orange. From then on nothing changes the wind: it runs along the isobars." `set` {t: 72}, `eq: 'geo'`, `derive: {id:
  'D03', step: 4}` · 4. "How good is it?" — "The Rossby number compares what we dropped with what we kept. Here it is 0.11."
  `readouts: ['Ro']`, `eq: 'ro'`, `code: {id: 'geo', lines: [3, 3]}` · 5. "The other hemisphere" — "Switch to S. The Coriolis
  force now points to the left of the motion and the low is circled the other way." `controls: ['hemi']` · 6. "Your turn" —
  "Predict the latitude below which this low stops being geostrophic. Then drag the latitude down." `controls: ['lat',
  'dp']`.
- **Equations:** `geo` "Geostrophic balance" ref 'Eq. (13.11)–(13.12)' $-fv=-\frac1{\rho_0}\frac{\partial p}{\partial x},\ fu=-\frac1{\rho_0}\frac{\partial p}{\partial y}$, live with the two
  sides at the parcel · `ro` "Rossby number" ref 'Eq. (13.13)' $\mathrm{Ro}=\frac{U^2/L}{fU}=\frac U{fL}$, live · `f` "Coriolis parameter" ref 'Eq.
  (13.8)' $f=2\Omega\sin\theta$, live · `psi` "Isobars are streamlines" ref 'D03 (no book number)' $\mathbf u\cdot\nabla p=0$, $\psi=p/(f\rho_0)$, note "only for
  constant f; here u = −∂ψ/∂y" · `ek` "Ekman number" ref 'Eq. (13.18)' $E=\frac{\rho\nu U/L^2}{\rho fU}=\frac\nu{fL^2}$.
- **Check yourself:** (1) "Halve the isobar spacing (double Δp). What happens to the wind speed and to its direction?" —
  "Speed doubles; direction does not change: the wind still follows the isobars." `set {dp: 8}` · (2) "Which way does air
  circle a low at 35° S?" — "Clockwise: f < 0 puts low pressure on the right of the wind." `set {hemi: 'S'}` · (3) "At what
  drag does the wind blow 45° across the isobars?" — "When r = abs(f): arctan 1 = 45°." `set {drag: 1}` · (4) "Why does the
  balance fail at 5°?" — "There f is seven times smaller than at 35°, so the wind needed to balance the same pressure force is
  seven times larger and Ro passes 1."
- **Selftest parity rows:** `{name: 'f 35N', js: fCor(35, 'N'), py: 'ch13.coriolis_parameter(np.deg2rad(35.0))', rtol: 1e-12}`
  (8.36517e-5) · `{name: 'f 35S', js: fCor(35, 'S'), py: 'ch13.coriolis_parameter(np.deg2rad(-35.0))', rtol: 1e-12}` · `{name: 'u_g',
  js: geoWind(0, 400/3e5, fCor(35, 'N'), 1.2)[0], py: 'ch13.geostrophic_velocity(0.0, 400.0/3.0e5, ch13.coriolis_parameter(np.deg2rad(35.0)),
  1.2)[0]', rtol: 1e-12}` (−13.2826) · `{name: 'parcel Re', js: parcel(7200, 5e-4, 0, fCor(35, 'N'), 2e-5, 0, 0)[0], py:
  'np.real(ch13.parcel_adjust(7200.0, 5.0e-4+0j, ch13.coriolis_parameter(np.deg2rad(35.0)), 2.0e-5, 0j))', rtol: 1e-10}` · the same
  with `np.imag` and index [1] · `{name: 'parcel S', js: parcel(7200, 5e-4, 0, fCor(35, 'S'), 0, 0, 0)[1], py:
  'np.imag(ch13.parcel_adjust(7200.0, 5.0e-4+0j, ch13.coriolis_parameter(np.deg2rad(-35.0)), 0.0, 0j))', rtol: 1e-10}` · `{name: 'Ro',
  js: rossby(14, fCor(35, 'N'), 1.4e6), py: 'ch13.rossby_number(14.0, ch13.coriolis_parameter(np.deg2rad(35.0)), 1.4e6)', rtol:
  1e-12}` (0.11954) · invariant `{name: 'forces close', js: forceSum(state)[0], expect: 0, atol: 1e-12}` at t = 72 h, drag 0.4.
- **Fit plan:** 360×640: status (1 line) · `map` (58 %) over `forces` (42 %); `time` hidden (its two numbers are the readouts
  "Speed now" and "Angle to isobars", and sit in the map's title); transport without step buttons; walkthrough card
  paged. 844×345: `map` beside `forces` in one row, text panel on the right, titles shortened to "Map" and "Forces".
  Desktop / 1000×700: rows [1.25, 1].

### E2 · thermal_wind
- **Title:** "Why is there a jet stream — and why westerly?" · **Summary:** "Drag the pole-to-equator temperature contrast and
  watch the isotherms tilt and the wind grow with height, in either hemisphere." · **CORE:** C03 (also N20, N24, N25, N147) ·
  **Reference:** `amplitude_phase_second_order_II_3.html` (linked windows; a numbered Explain with live numbers).
- **meta:** `viz:order 2` · `viz:sections 13.5 13.17` · `viz:equations 13.14 13.15 13.21` · `viz:fluidpy ch13.thermal_wind_shear
  ch13.thermal_wind_from_temperature ch13.thermal_wind_integrate ch13.jet_section ch13.coriolis_parameter` · `viz:derivations D04`.
- **Physics (JS ↔ Python):** `shearT(dTdx, dTdy, f, alpha)` → [du/dz, dv/dz] ↔ `ch13.thermal_wind_from_temperature(dTdx, dTdy, f,
  alpha=alpha)`; `shearRho(drdx, drdy, f, rho0)` ↔ `ch13.thermal_wind_shear`; `jetT(y, z, dT, w, sgn)` and `jetU(y, z, dT, w, lat,
  alpha, u0)` (analytic: $T=T_0+\gamma z-\tfrac{\Delta T}2\tanh(y/w)\,\mathrm{sgn}(\theta)$, $U=u_0+\frac{g\alpha}{f}\frac{\Delta T}{2w}\mathrm{sech}^2(y/w)\,\mathrm{sgn}(\theta)\,z$) ↔ `ch13.jet_section(...)["T"]`,
  `["U"]`; `profileU(z, dTdy, f, alpha, u0)` ↔ `ch13.thermal_wind_integrate` (in temperature form: linear in z here).
- **Views** (rows [1.3, 1]): 1. `section` "Latitude–height section" (row 0, flex 1.5): y from −3000 to 3000 km (pole on the
  left for N, on the right for S; the axis is labelled "toward the pole"), z 0…12 km (atmosphere) or 0…−1 km (ocean mode);
  isotherms (blue contours every 4 K), zonal wind as a filled field (purple = westerly, muted = easterly) with contours
  every 5 m/s, the jet core marked ✕ with its speed; a vertical cursor line at the chosen y. Pointer: drag the cursor;
  click a point for the inspector. 2. `profile` "Wind against height here" (row 0, flex 0.8): U(z) at the cursor (purple),
  the surface wind as an anchor dot, the no-gradient profile as a muted vertical ghost. 3. `bars` "The two sides of the
  relation" (row 1, `hidePortrait`): two bars that must be equal — the shear ∂u/∂z, and −(gα/f) ∂T/∂y — plus the three
  factors g α/f, ∂T/∂y as small bars with signs.
- **Controls:** `dT` "Temperature contrast $\Delta T$" (0…50 K, step 1, default 28; help: "0: no gradient, no shear · 50: a
  fierce winter contrast") · `w` "Half-width of the front" (500…3000 km, default 2000) · `lat` "Latitude" (15…70°, default 35)
  with `hemi` chips N / S · `u0` "Surface wind" (−10…10 m/s, default 0; optional) · `sys` modes "atmosphere (T) / ocean (ρ)".
- **Presets:** "winter jet" {dT 40, w 1500} · "weak summer gradient" {dT 12} · "southern hemisphere" {hemi S} · "no gradient
  (Taylor–Proudman)" {dT 0} · "an ocean front" {sys ocean, Δρ 1 kg m⁻³ over 100 km}.
- **Status:** dT > 0: "westerly shear +2.9 m/s per km · jet 26 m/s at 9 km" (the same words in both hemispheres — and the
  page adds "(f and ∂T/∂y both changed sign)" when S is chosen); dT = 0: "no horizontal gradient → no shear: the wind is the
  same at every height (Taylor–Proudman)"; ocean mode: "current shear 1.1 × 10⁻³ s⁻¹ across the front".
- **Readouts:** "$\partial T/\partial y$ here" · "Shear $\partial u/\partial z$" · "Wind at 9 km" · "Isotherm slope".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (atmosphere / ocean) · **presets** ·
  **status** · **inspector** (click a point of the section: "∂T/∂y = −(ΔT/2w) sech²(y/w) = … K/m; ∂u/∂z = −(gα/f)(∂T/∂y) = … s⁻¹;
  U = u₀ + shear × z = … m/s").
- **Explain:**
  0. *What the views show* — "**Section:** <b class=c-blue>blue</b> lines are isotherms — cold toward the pole; shading is
     the eastward wind. **Profile:** the wind above the dashed cursor. **Bars** (hidden on phones): the two sides of the
     thermal-wind relation."
  1. *The temperature gradient under the cursor* — "∂T/∂y = −ΔT/(2w) at the centre of the front = −28/(2 × 2 × 10⁶) =
     **−7.0 × 10⁻⁶ K/m** (7 K colder per 1000 km toward the pole)."
  2. *The shear it demands* — "$\frac{\partial u}{\partial z}=\frac g{\rho_0f}\frac{\partial\rho}{\partial y}$ (the second of (13.15)); with ρ′ = −ρ₀αT′ (ours) $\frac{\partial u}{\partial z}=-\frac{g\alpha}f\frac{\partial T}{\partial y}$ =
     (9.81/280)/(8.37 × 10⁻⁵) × 7.0 × 10⁻⁶ = **2.93 × 10⁻³ s⁻¹** = 2.9 m/s per km."
  3. *The wind aloft* — "U(z) = u₀ + shear × z = 0 + 2.93 × 10⁻³ × 9000 = **26 m/s** at 9 km — the jet."
  4. *Why the isotherms slope* — "slope = −(∂T/∂y)/(∂T/∂z); with the background lapse rate dT/dz = −6.5 K/km (Kundu's sign;
     meteorologists write Γ_met = +6.5 K/km) the slope is 7.0 × 10⁻⁶/6.5 × 10⁻³ ≈ **1.1 × 10⁻³**: isotherms dip toward the equator
     by about 1 km in 1000 km."
  5. *The bars* (or the hint "turn the phone sideways for the bars") — "left bar: shear; right bar: −(gα/f) ∂T/∂y. They are
     always equal; flip the hemisphere and both factors on the right change sign."
  6. *At the cursor* — "at y = `Viz.live('y')` km: ∂T/∂y = `Viz.live('dTdy')` K/m, shear `Viz.live('sh')` s⁻¹, wind at 9 km
     `Viz.live('U9')` m/s."
  7. *Reading the current setting* — dT > 0, N: "Cold air toward the pole: pressure falls faster with height on the cold
     side, so the poleward pressure force — and the westerly wind that balances it — grows with height. The jet sits
     where the isotherms are steepest." · S: "Now cold is toward the south pole (∂T/∂y > 0) and f < 0: two sign changes, the
     same westerly shear. The jet stream is westerly in both hemispheres." · dT = 0: "No horizontal gradient, nothing to
     drive a shear: $\partial\mathbf u/\partial z=0$ (13.21), the Taylor–Proudman state — fluid moves in vertical columns." · ocean: "The same
     relation in density, $\frac{\partial u}{\partial z}=\frac g{\rho_0f}\frac{\partial\rho}{\partial y}$: oceanographers measure a density section and integrate upward from a
     level assumed at rest to get the current."
- **Derivation tab:** **D04** (7 steps) `view: 'section'`; goal `set` {preset winter jet}; step 3 `watch` "the isobars
  (toggle them on) slope more steeply at each higher level"; step 5 `live` "∂u/∂z = (g/ρ₀f) ∂ρ/∂y = …"; step 7 `set` {sys:
  atmosphere}, `live` "−(gα/f) ∂T/∂y = −(0.0350/8.37 × 10⁻⁵) × (−7.0 × 10⁻⁶) = 2.93 × 10⁻³ s⁻¹".
- **Code:**
  ```python
  f = ch13.coriolis_parameter(np.deg2rad({{lat}}))          # f = {{f}} 1/s (negative in the south)
  dTdy = {{dTdy}}                                            # K/m, from the contrast and the width
  dudz, dvdz = ch13.thermal_wind_from_temperature(0.0, dTdy, f, alpha=1/280)   # ours: {{dudz}} 1/s
  U9 = {{u0}} + dudz*9000.0                                  # wind at 9 km: {{U9}} m/s
  sec = ch13.jet_section(y, z, dT={{dT}}, width={{w}}, lat_rad=np.deg2rad({{lat}}), alpha=1/280)
  ```
- **Walkthrough (6 steps):** 1. "A river of air" — "High above the calm ground blows a jet. What sets its speed? Drag the
  temperature contrast." `controls: ['dT']` · 2. "Tilted isotherms" — "Colder toward the pole means the blue lines slope.
  Where they slope most, the wind changes fastest with height." `set` {dT: 28} · 3. "The relation" — "Shear equals −(gα/f)
  times the temperature gradient: 2.9 m/s per kilometre here." `eq: 'tw'`, `derive: {id: 'D04', step: 7}` · 4. "Up to the
  jet" — "Add up the shear over 9 km: 26 m/s from the west." `readouts: ['U9']`, `code: {id: 'tw', lines: [3, 4]}` · 5.
  "Southern hemisphere" — "Switch to S: the cold is on the other side and f is negative. Two sign changes — still
  westerly." `controls: ['hemi']` · 6. "Your turn" — "Predict the wind profile when the contrast is zero. Then set it."
  `controls: ['dT', 'sys']`.
- **Equations:** `tw` "Thermal wind" ref 'Eq. (13.15)' $\frac{\partial v}{\partial z}=-\frac g{\rho_0f}\frac{\partial\rho}{\partial x},\ \frac{\partial u}{\partial z}=\frac g{\rho_0f}\frac{\partial\rho}{\partial y}$, live · `twT` "In
  temperature (ours)" ref 'not in the book — ours' $\frac{\partial u}{\partial z}=-\frac{g\alpha}f\frac{\partial T}{\partial y}$, live · `hyd` "Hydrostatic perturbation" ref 'Eq. (13.14)'
  $0=-\frac{\partial p}{\partial z}-g\rho$ · `tp` "Taylor–Proudman" ref 'Eq. (13.21)' $\partial\mathbf u/\partial z=0$, note "the no-gradient limit".
- **Check yourself:** (1) "Double the width of the front at the same contrast. What happens to the jet?" — "It halves: the
  gradient ΔT/2w halves." `set {w: 3000}` · (2) "Why is the jet westerly in the southern hemisphere too?" — "Both f and
  ∂T/∂y change sign." `set {hemi: 'S'}` · (3) "Give the surface a 5 m/s easterly. What changes aloft?" — "Every level shifts by
  −5 m/s; the shear is unchanged." `set {u0: -5}` · (4) "At the same contrast, is the jet stronger at 30° or at 60°?" —
  "At 30°: the shear goes as 1/f."
- **Selftest parity rows:** `{name: 'shear N', js: shearT(0, -7e-6, fCor(35, 'N'), 1/280)[0], py:
  'ch13.thermal_wind_from_temperature(0.0, -7.0e-6, ch13.coriolis_parameter(np.deg2rad(35.0)), alpha=1/280)[0]', rtol: 1e-12}` (2.9308e-3)
  · `{name: 'shear S', js: shearT(0, 7e-6, fCor(35, 'S'), 1/280)[0], py: 'ch13.thermal_wind_from_temperature(0.0, 7.0e-6,
  ch13.coriolis_parameter(np.deg2rad(-35.0)), alpha=1/280)[0]', rtol: 1e-12}` (the same number) · `{name: 'density form', js:
  shearRho(0, 1e-5, fCor(35, 'N'), 1027)[0], py: 'ch13.thermal_wind_shear(0.0, 1.0e-5, ch13.coriolis_parameter(np.deg2rad(35.0)),
  1027.0)[0]', rtol: 1e-12}` (1.1415e-3) · `{name: 'dv/dz', js: shearT(3e-6, 0, fCor(35, 'N'), 1/280)[1], py:
  'ch13.thermal_wind_from_temperature(3.0e-6, 0.0, ch13.coriolis_parameter(np.deg2rad(35.0)), alpha=1/280)[1]', rtol: 1e-12}` ·
  invariant `{name: 'bars equal', js: barsGap(state), expect: 0, atol: 1e-14}`.
- **Fit plan:** 360×640: status · `section` (60 %) over `profile` (40 %); `bars` hidden (the shear and its two factors are in
  the profile's title: "∂u/∂z = 2.9 m/s per km = −(gα/f) ∂T/∂y"); no transport in this explainer. 844×345: `section` beside
  `profile`. Desktop: rows [1.3, 1].

### E3 · ekman_spiral
- **Title:** "The wind blows east — why does the water go south?" · **Summary:** "Orbit the spiral, slide a depth cursor down
  it and watch the summed transport swing to 90° from the wind." · **CORE:** C04, C05 (also N30, N31, N32, N33, N35) ·
  **Reference:** `angular_frequency_explorer_1.html` (system animation + pointer + graph on one clock, presets, live working).
- **meta:** `viz:order 3` · `viz:sections 13.6` · `viz:equations 13.22 13.23 13.27 13.28 13.29 13.30` · `viz:fluidpy ch13.ekman_depth
  ch13.ekman_surface ch13.ekman_transport ch13.ekman_transport_partial ch13.ekman_residual ch13.coriolis_parameter` ·
  `viz:derivations D06 D07`.
- **Physics (JS ↔ Python):** `ekDepth(nu, f, conv)` ↔ `ch13.ekman_depth(nu_v, f, convention=)`; `ekSurface(z, tx, ty, rho, nu, f,
  Ug, Vg)` → [u, v] (complex arithmetic by hand: $V=V_g+\frac{(\tau_x+i\tau_y)(1-is)}{\rho\sqrt{2\nu_v\lvert f\rvert}}e^{(1+is)z/\delta}$) ↔ `ch13.ekman_surface`; `ekTransport(tx, ty,
  rho, f)` ↔ `ch13.ekman_transport`; `ekPartial(z, tx, ty, rho, nu, f)` ↔ `ch13.ekman_transport_partial`; `ekResidual(z, …)` (the
  two terms of each component, analytic second derivative) ↔ `ch13.ekman_residual`.
- **Views** (rows [1.3, 1]): 1. `spiral` "The spiral in 3-D" (row 0, flex 1.2; `Viz.three`, drag to orbit): the wind arrow
  above a translucent sea surface, 24 current arrows hanging at depths down to 4δ (teal → muted with depth), the tips
  joined by a line; the running transport as a thick purple arrow on the surface plane; the cursor depth as a
  translucent plate. Falls back to an oblique 2-D drawing if WebGL is missing (with a one-line note). 2. `hodo`
  "Hodograph" (row 0, flex 1, `equal`): v against u; the spiral with marks at −z/δ = π/8, π/4, π/2, π; the wind direction
  (grey), the cursor's velocity vector (teal, bold), the part of the spiral above the cursor bold; the transport
  direction as a dashed purple ray. 3. `prof` "Profiles and running transport" (row 1, `hidePortrait`): u, v against −z/δ
  (left axis), and the running transport components from the cursor to the surface (right axis, purple) with the final
  value as a ghost line.
- **Controls:** `tau` "Wind stress $\tau$" (0.01…0.3 N/m², default 0.07) · `dir` "Wind toward" (0…360°, step 15, default 90 =
  east; optional) · `nu` "Eddy viscosity $\nu_v$" (0.005…0.3 m²/s, log, default 0.03; help: "0.005: a thin, fast layer · 0.3:
  deep and slow") · `lat` "Latitude" (5…80°, default 60) with `hemi` chips N / S · `Ug` "Interior flow" (0…0.2 m/s, default 0;
  optional) · transport `zc` (the depth cursor).
- **Transport:** `zc` from 0 down to −5δ ("play = integrate downward"), rate 1 δ/s, `end: 'hold'`; end card: "summed over
  the layer: 0.540 m²/s, 90° to the right of the wind — whatever ν_v is".
- **Presets:** "four times ν_v" {nu 0.12}: twice as deep, same transport · "35° S" {hemi S, lat 35}: to the left · "wind along a
  coast" {dir 180 (toward the south), a coast drawn on the east side}: transport offshore → upwelling · "with interior
  flow" {Ug 0.05}.
- **Status:** "transport 0.540 m²/s, 90° to the right of the wind · surface current 45° to the right · δ = 21.8 m" (f > 0) /
  "… to the left …" (f < 0); coast preset adds "→ offshore: upwelling" or "→ onshore: downwelling".
- **Readouts:** "Thickness $\delta$" · "Ekman depth $\pi\delta$" · "Surface speed" · "Transport $\tau/\rho\lvert f\rvert$" · "Summed so far" · "Angle
  so far".
- **Depth features:** Explain + Code + Derivation · **3-D view** · **linked views** (3) · **transport** · **presets** ·
  **status** · **terms** (at the cursor depth, per component: Coriolis (teal) + friction (rose) = 0 — two bars for x, two
  for y; click to isolate).
- **Explain:**
  0. *What the views show* — "grey arrow = wind; <b class=c-teal>teal</b> arrows = the current at each depth; <b
     class=c-accent>purple</b> = the transport, i.e. all the currents above the cursor added up; <b class=c-rose>rose</b>
     bars = friction, the pull of the layers above and below."
  1. *How deep* — "$\delta=\sqrt{2\nu_v/f}$ (13.29) = √(2 × 0.03/1.263 × 10⁻⁴) = **21.8 m**: the current falls by e and turns one radian
     every δ. ⚠️ Oceanographers' 'Ekman depth' is πδ = **68.5 m**, where the current first opposes the surface current."
  2. *How fast at the surface* — "$V_0=\frac{\tau/\rho}{\sqrt{\lvert f\rvert\nu_v}}$ = (0.07/1027)/√(1.263 × 10⁻⁴ × 0.03) = **0.035 m/s**, pointing 45° to the right
     of the wind (to the left for f < 0)."
  3. *At the cursor* — "at z = −δ: speed × e⁻¹ = 0.0129 m/s, turned a further 57.3°."
  4. *The transport* — "$\int_{-\infty}^0u\,dz=0$, $\int_{-\infty}^0v\,dz=-\frac\tau{\rho f}$ (13.30): 0.07/(1027 × 1.263 × 10⁻⁴) = **0.540 m²/s**, at 90° to the
     right. No ν_v in it."
  5. *Summed so far* — "above the cursor the transport is (1 − e^{(1+i)z/δ}) of the total: `Viz.live('frac')` of its
     size, `Viz.live('ang')`° from the wind."
  6. *The profile view* (or the hint "turn the phone sideways for the profiles").
  7. *Reading the current setting* — default: "Each layer is dragged by the one above and turned by the Coriolis force,
     so the arrows rotate and shrink going down. Their along-wind parts cancel in the sum; only the part at right angles
     survives." · ν_v changed: "A larger eddy viscosity spreads the same momentum over a deeper layer: slower at the
     surface, deeper spiral, identical transport. The transport is fixed by the stress and by f alone." · S: "With f < 0 every
     'right' becomes 'left'." · coast: "The transport is offshore, so water must rise at the coast to replace it:
     upwelling — the cold, fertile water off Peru and California."
- **Derivation tab:** **D06** (14 steps) `view: 'hodo'`; goal `set` {zc: 0}; step 2 `watch` "multiplying by i turns an arrow
  a quarter turn to the left — exactly what the Coriolis terms do"; step 7 `set` {zc: −3δ}, `watch` "the kept root decays
  downward"; step 10 `live` "A = τδ(1 − i)/(2ρν_v) = (0.0248 − 0.0248 i) m/s"; step 13 `live` "u = 0.035 e^{z/δ} cos(−z/δ + π/4)
  = `Viz.live('u')`". **D07** (8 steps) `view: 'hodo'` (phones) / `'prof'` (wide); step 4 `live` "∫V dz = −iτ/(ρf) = −0.540 i
  m²/s"; step 7 `set` {nu: 0.12}, `watch` "the purple arrow does not move".
- **Code:**
  ```python
  f = ch13.coriolis_parameter(np.deg2rad({{lat}}))               # {{f}} 1/s
  delta = ch13.ekman_depth({{nu}}, f)                            # (13.29): {{delta}} m; pi*delta = {{pid}} m
  u, v = ch13.ekman_surface({{zc}}, {{tx}}, {{ty}}, 1027.0, {{nu}}, f)   # current at the cursor: ({{u}}, {{v}}) m/s
  Mx, My = ch13.ekman_transport({{tx}}, {{ty}}, 1027.0, f)       # (13.30): ({{Mx}}, {{My}}) m^2/s
  mx, my = ch13.ekman_transport_partial({{zc}}, {{tx}}, {{ty}}, 1027.0, {{nu}}, f)   # summed so far
  ```
- **Walkthrough (7 steps):** 1. "Downwind?" — "The wind blows east over the sea. Which way does the water go? Drag the 3-D
  view to look." · 2. "Turning with depth" — "Each layer is dragged by the one above and turned by Coriolis: the arrows
  rotate and shrink." `set` {zc: 0}, `play: true` · 3. "One complex equation" — "Write V = u + iv. Friction against
  Coriolis becomes one line." `eq: 'cplx'`, `derive: {id: 'D06', step: 2}` · 4. "How deep" — "Thickness δ = 21.8 m; the
  surface current is 45° to the right of the wind." `readouts: ['delta', 'V0']`, `eq: 'sol'` · 5. "Add them up" — "Play: the
  purple sum swings round and settles at 90° from the wind." `play: true`, `terms: true` · 6. "What viscosity does" — "Four
  times the eddy viscosity: twice as deep, half as fast — same transport." `set` {nu: 0.12}, `code: {id: 'ek', lines: [4, 4]}`
  · 7. "Your turn" — "Predict the transport's direction at 35° S for a wind toward the north. Then set it." `controls:
  ['hemi', 'dir']`.
- **Equations:** `bal` "Coriolis against friction" ref 'Eq. (13.22)–(13.23)' $-fv=\nu_v\frac{d^2u}{dz^2},\ fu=\nu_v\frac{d^2v}{dz^2}$, live at the cursor ·
  `cplx` "One complex equation" ref 'Eq. (13.27)' $\frac{d^2V}{dz^2}=\frac{if}{\nu_v}V,\ V\equiv u+iv$ · `sol` "Solution and thickness" ref 'Eq.
  (13.28)–(13.29)' $V=Ae^{(1+i)z/\delta}+Be^{-(1+i)z/\delta},\ \delta=\sqrt{2\nu_v/f}$, live, note "B = 0: nothing may grow downward; f > 0 as printed" ·
  `tr` "Transport" ref 'Eq. (13.30)' $\int_{-\infty}^0u\,dz=0,\ \int_{-\infty}^0v\,dz=-\frac\tau{\rho f}$, live.
- **Check yourself:** (1) "Quadruple ν_v. Which of δ, surface speed, transport change?" — "The thickness δ doubles, surface speed
  halves, transport stays." `set {nu: 0.12}` · (2) "How deep must you integrate to get 86 % of the transport's size?" — "To
  z = −δ: abs(1 − e^{−(1+i)}) = 0.859." `set {zc: -1}` · (3) "A wind blows toward the equator along an eastern ocean boundary
  in the northern hemisphere. Upwelling or downwelling?" — "Upwelling: the transport is to the right of the wind =
  offshore." · (4) "Where is the current exactly opposite to the surface current?" — "At z = −πδ, with 4.3 % of its speed."
- **Selftest parity rows:** `{name: 'delta', js: ekDepth(0.03, fCor(60, 'N'), 'efold'), py: 'ch13.ekman_depth(0.03,
  ch13.coriolis_parameter(np.deg2rad(60.0)))', rtol: 1e-12}` (21.7956) · `{name: 'u surface N', js: ekSurface(0, 0.07, 0, 1027, 0.03,
  fCor(60, 'N'), 0, 0)[0], py: 'ch13.ekman_surface(0.0, 0.07, 0.0, 1027.0, 0.03, ch13.coriolis_parameter(np.deg2rad(60.0)))[0]', rtol:
  1e-12}` (0.024760) · `{name: 'v surface S', js: ekSurface(0, 0.07, 0, 1027, 0.03, fCor(60, 'S'), 0, 0)[1], py: 'ch13.ekman_surface(0.0,
  0.07, 0.0, 1027.0, 0.03, ch13.coriolis_parameter(np.deg2rad(-60.0)))[1]', rtol: 1e-12}` (+0.024760) · `{name: 'v at -delta', js:
  ekSurface(-21.7956, …)[1], py: 'ch13.ekman_surface(-21.7956, 0.07, 0.0, 1027.0, 0.03, ch13.coriolis_parameter(np.deg2rad(60.0)))[1]',
  rtol: 1e-10}` · `{name: 'transport', js: ekTransport(0.07, 0, 1027, fCor(60, 'N'))[1], py: 'ch13.ekman_transport(0.07, 0.0, 1027.0,
  ch13.coriolis_parameter(np.deg2rad(60.0)))[1]', rtol: 1e-12}` (−0.53965) · `{name: 'partial', js: ekPartial(-21.7956, 0.07, 0, 1027,
  0.03, fCor(60, 'N'))[1], py: 'ch13.ekman_transport_partial(-21.7956, 0.07, 0.0, 1027.0, 0.03,
  ch13.coriolis_parameter(np.deg2rad(60.0)))[1]', rtol: 1e-10}` · invariant `{name: 'terms cancel', js: ekResidual(-10, …)[0],
  expect: 0, atol: 1e-15}`.
- **Fit plan:** 360×640: status · `spiral` (55 %) over `hodo` (45 %); `prof` hidden (the running transport and its angle
  are in the hodograph's title: "summed so far 0.46 m²/s at 61°"); if WebGL is missing the phone shows `hodo` over `prof`
  instead. 844×345: `spiral` beside `hodo`; transport slider under them without step buttons. Desktop: rows [1.3, 1].

### E4 · ekman_force_balance
- **Title:** "Why does air spiral into a low?" · **Summary:** "Slide from the ground to the top of the friction layer and watch
  three forces keep closing a triangle while the wind turns toward low pressure." · **CORE:** C06 (also N39, N47, N48, N49 and
  the bottom pumping of N36) · **Reference:** `forced_damped_vibrations.html`.
- **meta:** `viz:order 4` · `viz:sections 13.7` · `viz:equations 13.28 13.29 13.33 13.34 13.41` · `viz:fluidpy ch13.ekman_bottom
  ch13.ekman_force_balance ch13.ekman_bottom_transport ch13.ekman_pumping_bottom ch13.ekman_depth ch13.eddy_viscosity_from_depth` ·
  `viz:derivations D08`.
- **Physics (JS ↔ Python):** `ekBottom(z, Ug, Vg, nu, f)` → [u, v] ($V=V_g(1-e^{-(1+is)z/\delta})$) ↔ `ch13.ekman_bottom`; `forces(z, Ug,
  nu, f, rho)` → {coriolis, pressure, friction, angle, sum} ↔ `ch13.ekman_force_balance`; `bottomTransport(Ug, Vg, nu, f)` ↔
  `ch13.ekman_bottom_transport`; `pumpBottom(zeta, nu, f)` ↔ `ch13.ekman_pumping_bottom`; `ekDepth` as in E3; `nuFromDepth(delta,
  f)` ↔ `ch13.eddy_viscosity_from_depth`.
- **Views** (rows [1.25, 1]): 1. `tri` "Forces at this height" (row 0, flex 1.2, `equal`): faint isobars (parallel lines,
  low pressure labelled on one side — the left for N, the right for S), the wind vector at the cursor height (grey, from
  the origin), and the three forces head-to-tail closing a triangle: pressure gradient (orange, always the same),
  Coriolis (teal, perpendicular to the wind), friction (rose); the cross-isobar angle drawn as an arc with its value.
  2. `hodo` "Hodograph" (row 0, flex 1, `equal`): from the origin to (U, 0) with height marks z/δ = π/4, π/2, π, the cursor
  dot, the part below the cursor bold; the point of largest u marked ▲ "3πδ/4: 1.067 U". 3. `plan` "A low and a high from
  above" (row 1, `hidePortrait`): two circular centres with geostrophic arrows aloft (muted) and surface-layer arrows at
  the cursor height turned inward (low) and outward (high); the pumping velocity printed at each centre with ↑ or ↓.
- **Controls:** `zc` (transport) "Height $z/\delta$" · `U` "Geostrophic wind $U$" (2…30 m/s, default 12) · `nu` "Eddy viscosity
  $\nu_v$" (0.5…40 m² s⁻¹, log, default 7; help: "0.5: a layer 90 m thick · 40: nearly 800 m") · `lat` "Latitude" (10…80°,
  default 60) with `hemi` chips N / S · `centre` chips "low / high".
- **Transport:** `zc` 0…5 (in units of δ), rate 0.5 /s, `end: 'hold'`; end card: "top of the layer: wind along the isobars,
  friction gone · summed inflow ½Uδ = 1998 m²/s toward low pressure".
- **Presets:** "at the ground" {zc 0.02} · "at z = πδ: v = 0, u = 1.043 U" {zc π} · "largest u" {zc 3π/4} · "top of the layer"
  {zc 5} · "35° S" {hemi S, lat 35} · "a high" {centre high}.
- **Status:** "at z = 0.5 δ: wind 32° across the isobars, toward low pressure (to the left of the geostrophic wind)" (f > 0) /
  "… (to the right …)" (f < 0); at the top: "geostrophic: wind along the isobars"; pumping: "low: air rises at 1.7 mm/s" /
  "high: air sinks".
- **Readouts:** "Thickness $\delta$" · "Wind here" · "Angle to isobars" · "Inflow $\tfrac12U\delta$" · "Pumping $w$".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **terms**
  (paired budget: x-components and y-components, three bars each — pressure + Coriolis + friction = 0).
- **Explain:**
  0. *What the views show* — "<b class=c-orange>orange</b> = pressure-gradient force (toward low pressure, the same at
     every height); <b class=c-teal>teal</b> = Coriolis force (perpendicular to the wind, proportional to its speed); <b
     class=c-rose>rose</b> = friction; grey = the wind. The three coloured arrows always close a triangle."
  1. *Thickness of the layer* — "$\delta=\sqrt{2\nu_v/f}$ (13.29) = √(2 × 7/1.263 × 10⁻⁴) = **333 m** (and πδ = 1046 m)."
  2. *The wind at the cursor* — "$u=U[1-e^{-z/\delta}\cos(z/\delta)]$, $v=Ue^{-z/\delta}\sin(z/\delta)$ (13.41); at z = 0.5δ: u = 12(1 − 0.6065 × 0.8776) =
     **5.61 m/s**, v = 12 × 0.6065 × 0.4794 = **3.49 m/s**; angle = arctan(v/u) = **31.9°**."
  3. *The three forces (per unit mass)* — "pressure: (0, fU) = (0, 1.52 × 10⁻³) m/s²; Coriolis: (fv, −fu) = (0.44, −0.71) ×
     10⁻³; friction: what closes the triangle, (−0.44, −0.81) × 10⁻³. Sum: zero."
  4. *The inflow* — "$\int_0^\infty v\,dz=\tfrac12U\delta$ = 0.5 × 12 × 333 = **1998 m²/s** toward low pressure."
  5. *The pumping* (plan view, or the readout on phones) — "(ours) $w=\tfrac\delta2\zeta_g$: for ζ_g = 10⁻⁵ s⁻¹, w = 166 × 10⁻⁵ = **1.7
     mm/s** upward out of the layer; a high (ζ_g < 0 for f > 0) pumps downward."
  6. *At the cursor* — "at z/δ = `Viz.live('zc')`: wind (`Viz.live('u')`, `Viz.live('v')`) m/s, angle `Viz.live('ang')`°."
  7. *Reading the current setting* — near the ground (z < 0.3δ): "The air is almost stopped, so the Coriolis force is
     almost gone; pressure is held by friction alone and the little wind there is blows 45° across the isobars." · middle
     (0.3δ–2δ): "A three-way balance. The wind has a component toward low pressure: this is the inflow that feeds a
     depression." · z ≈ 3πδ/4: "The wind is 6.7 % *faster* than geostrophic — friction here is pushing forward, handing
     down momentum from above." · top (z > 3δ): "Friction has faded; Coriolis alone balances pressure; the wind follows
     the isobars." · S: "Left and right swap; the inflow is still toward low pressure." · high: "The surface wind leaves
     the centre, air sinks from above and clouds dissolve."
- **Derivation tab:** **D08** (9 steps) `view: 'tri'` (phones: `'hodo'`); goal `set` {zc: 0.5}; step 2 `watch` "the constant
  term fU is the orange arrow"; step 5 `set` {zc: 0.02}, `watch` "V = 0 at the ground"; step 6 `live` "u = `Viz.live('u')`, v
  = `Viz.live('v')`"; step 7 `set` {zc: 0.02}, `live` "angle → 45°"; step 9 `live` "½Uδ = `Viz.live('M')` m²/s".
- **Code:**
  ```python
  f = ch13.coriolis_parameter(np.deg2rad({{lat}}))            # {{f}} 1/s
  delta = ch13.ekman_depth({{nu}}, f)                         # (13.29): {{delta}} m
  u, v = ch13.ekman_bottom({{z}}, {{U}}, 0.0, {{nu}}, f)      # (13.41): ({{u}}, {{v}}) m/s at z = {{z}} m
  fb = ch13.ekman_force_balance({{z}}, {{U}}, {{nu}}, f)      # three forces; angle = {{ang}} deg
  Mx, My = ch13.ekman_bottom_transport({{U}}, 0.0, {{nu}}, f) # cross-isobar inflow = {{My}} m^2/s
  w = ch13.ekman_pumping_bottom({{zeta}}, {{nu}}, f)          # ours: {{w}} m/s
  ```
- **Walkthrough (6 steps):** 1. "Into the low" — "Aloft the wind follows the isobars. Near the ground it does not. Slide the
  height and watch." `play: true` · 2. "Three forces" — "Orange never changes. Teal shrinks as the wind slows. Rose,
  friction, closes the triangle." `set` {zc: 0.5}, `terms: true` · 3. "The solution" — "The same complex equation as the
  ocean spiral, now with a constant forcing fU." `eq: 'sol'`, `derive: {id: 'D08', step: 2}` · 4. "At the ground" — "Very
  near the surface the wind blows 45° across the isobars." `set` {zc: 0.02}, `readouts: ['ang']` · 5. "What it adds up to" —
  "Summed over the layer, ½Uδ flows toward low pressure — and has to rise." `eq: 'tr'`, `code: {id: 'eb', lines: [5, 6]}` · 6.
  "Your turn" — "Predict where u is largest and by how much it exceeds U. Then find it." `controls: ['zc', 'hemi']`.
- **Equations:** `bal` "Three-way balance" ref 'Eq. (13.33)–(13.34)' $-fv=\nu_v\frac{d^2u}{dz^2},\ fu=\nu_v\frac{d^2v}{dz^2}+fU$, live · `sol` "The layer"
  ref 'Eq. (13.41)' $u=U[1-e^{-z/\delta}\cos(z/\delta)],\ v=Ue^{-z/\delta}\sin(z/\delta)$, live · `d` "Thickness" ref 'Eq. (13.29)' $\delta=\sqrt{2\nu_v/f}$ · `tr` "Inflow"
  ref 'p. 641 (no number)' $\int_0^\infty v\,dz=\tfrac12U\delta$ · `w` "Pumping (ours)" ref 'not in the book — ours' $w=\tfrac\delta2\zeta_g$.
- **Check yourself:** (1) "Double U. What happens to the cross-isobar angle at z = δ?" — "Nothing: the angle depends only on
  z/δ." `set {U: 24}` · (2) "Where does u exceed U most, and by how much?" — "At z = 3πδ/4, by 6.7 %." `set {zc: 2.356}` · (3)
  "In the southern hemisphere, does surface air flow into a low or out of it?" — "Into it: toward low pressure in both
  hemispheres; only the turning sense changes." `set {hemi: 'S'}` · (4) "Which force balances pressure right at the ground?"
  — "Friction: the wind, and with it the Coriolis force, vanish there."
- **Selftest parity rows:** `{name: 'u at 0.5 delta', js: ekBottom(0.5*332.933, 12, 0, 7, fCor(60, 'N'))[0], py:
  'ch13.ekman_bottom(0.5*332.933, 12.0, 0.0, 7.0, ch13.coriolis_parameter(np.deg2rad(60.0)))[0]', rtol: 1e-10}` (5.6126) · `{name: 'v S',
  js: ekBottom(0.5*332.933, 12, 0, 7, fCor(60, 'S'))[1], py: 'ch13.ekman_bottom(0.5*332.933, 12.0, 0.0, 7.0,
  ch13.coriolis_parameter(np.deg2rad(-60.0)))[1]', rtol: 1e-10}` (−3.4894) · `{name: 'angle', js: forces(166.47, 12, 7, fCor(60, 'N'),
  1).angle, py: 'ch13.ekman_force_balance(166.47, 12.0, 7.0, ch13.coriolis_parameter(np.deg2rad(60.0)))["angle_to_isobars"]', rtol:
  1e-9}` · `{name: 'transport', js: bottomTransport(12, 0, 7, fCor(60, 'N'))[1], py: 'ch13.ekman_bottom_transport(12.0, 0.0, 7.0,
  ch13.coriolis_parameter(np.deg2rad(60.0)))[1]', rtol: 1e-12}` (1997.60) · `{name: 'pumping S', js: pumpBottom(-1e-5, 7, fCor(60,
  'S')), py: 'ch13.ekman_pumping_bottom(-1.0e-5, 7.0, ch13.coriolis_parameter(np.deg2rad(-60.0)))', rtol: 1e-12}` (+1.6647e-3) ·
  invariants `{name: 'triangle closes', js: forces(100, 12, 7, fCor(60, 'N'), 1).sum[0], expect: 0, atol: 1e-15}` and `{name:
  'overshoot', js: maxU(12, 7, fCor(60, 'N'))/12, py: '1 + math.exp(-0.75*math.pi)/math.sqrt(2)', rtol: 1e-6}` (1.06702).
- **Fit plan:** 360×640: status · `tri` (55 %) over `hodo` (45 %); `plan` hidden (the pumping velocity is the readout
  "Pumping w", and the status names the centre). 844×345: `tri` beside `hodo`. Desktop: rows [1.25, 1].

### E5 · vertical_modes
- **Title:** "How can one layer stand for a stratified ocean?" · **Summary:** "Reshape N(z) and watch the mode shapes, their
  speeds and their equivalent depths respond." · **CORE:** C08 (also N53, N64, N71–N73, N75, N76, R24) · **Reference:**
  `amplitude_phase_second_order_II_3.html`.
- **meta:** `viz:order 5` · `viz:sections 13.9 13.12` · `viz:equations 13.56 13.62 13.64 13.65 13.69 13.71` · `viz:fluidpy
  ch13.vertical_modes ch13.vertical_modes_shooting ch13.modes_uniform_N ch13.equivalent_depth ch13.baroclinic_mode_speed
  ch13.rossby_radius ch13.thermocline_N2 ch13.wkb_mode_speed` · `viz:derivations D11`.
- **Physics (JS ↔ Python):** `N2prof(z, p)` (uniform, or the thermocline shape $N=N_d+(N_p-N_d)e^{-((z-z_t)/w)^2}$) ↔
  `ch13.thermocline_N2(z, N_deep=, N_peak=, z_t=, width=)`; `shootMode(n, p, lid)` → {c, psi[]}: integrate $S'=-\psi/c^2$,
  $\psi'=N^2S$ (with $S=\psi'/N^2$) from the bottom by RK4 on 400 steps, bisection on c so that the surface condition holds and ψ
  has n sign changes ↔ `ch13.vertical_modes_shooting(z, N2_fn, n=n, lid=)` (rtol 1e-8) and `ch13.vertical_modes(...)[0][n]` (rtol
  1e-4; index n for a free surface, **index n − 1 for `lid="rigid"`**, whose `Modes` has no barotropic entry); `uniformRoot(N, H, n)` (bisection on tan X = (N²H/g)/X) ↔ `ch13.modes_uniform_N(N, H)[0][n]`; `He(c)` ↔
  `ch13.equivalent_depth`; `cRigid(N, H, n)` ↔ `ch13.baroclinic_mode_speed`; `Lam(c, f)` ↔ `ch13.rossby_radius`; `cWKB(p, n)` ↔
  `ch13.wkb_mode_speed`.
- **Views** (rows [1.3, 1]): 1. `shape` "N(z) and the mode" (row 0, flex 1.2): depth axis −4200…0 m; N(z) as a filled blue
  curve (top axis, s⁻¹); the selected mode ψ_n(z) bold purple, its w-structure (the running integral of ψ_n) thin teal,
  the other modes faint; zero crossings marked ○; a depth cursor. Pointer: click a depth for the inspector; drag the
  thermocline handle. 2. `ladder` "The ladder of speeds" (row 0, flex 0.8): c_n on a logarithmic axis for n = 0…4, each
  rung labelled with c_n, H_e and Λ_n = c_n/abs(f); the selected rung bold; ghost rungs = the uniform-N values NH/(nπ) and
  the WKB estimate. 3. `wave` "The selected mode as a wave" (row 1, `hidePortrait`): displacement of five density surfaces
  against x on one clock, travelling at c_n (transport).
- **Controls:** `prof` chips "uniform N / thermocline" · `zt` "Thermocline depth" (100…800 m, default 300) · `wd` "Thermocline
  thickness" (60…400 m, default 150; help: "60: a sharp step · 400: nearly smooth") · `n` "Mode number $n$" (0…4, default 1) ·
  `lid` toggle "rigid lid" (default off = free surface) · `lat` "Latitude" (5…70°, default 35; optional — only the radius
  uses it).
- **Transport:** `t` 0…2 periods of the selected mode's wave (wavelength 20Λ_n), loop.
- **Presets:** "uniform N: exact cosines" · "sharp shallow thermocline" {zt 150, wd 80} · "deep weak thermocline" {zt 600, wd
  300} · "rigid lid" {lid on, n 1}.
- **Status:** "mode 1: c = 3.61 m/s · H_e = 1.33 m · Λ = 43 km at 35°" (uniform N defaults) / thermocline defaults: "mode 1:
  c = 1.77 m/s · H_e = 0.32 m · Λ = 21 km" / n = 0: "barotropic: c = 203 m/s ≈ √(gH)" / lid on with n = 0: "the rigid lid
  removes the barotropic mode (c → ∞): pick n ≥ 1".
- **Readouts:** "Speed $c_n$" · "Equivalent depth $H_e$" · "Rossby radius $\Lambda_n$" · "Zero crossings" · "WKB estimate".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** ·
  **inspector** (click a depth: "here ψ_n = …, (1/N²)ψ_n′ = …; the two terms of the mode equation: d/dz[(1/N²)ψ_n′] = …, ψ_n/c_n² = …;
  sum = 0").
- **Explain:**
  0. *What the views show* — "<b class=c-blue>blue</b> = the buoyancy frequency N(z); <b class=c-accent>purple</b> = the
     shape ψ_n(z) shared by the horizontal velocity and the pressure of mode n; <b class=c-teal>teal</b> = the shape of
     its vertical velocity. The ladder lists the speed of each mode."
  1. *The eigenvalue for your profile* — "we search for the c that lets $\frac d{dz}\big(\frac1{N^2}\frac{d\psi_n}{dz}\big)+\frac1{c_n^2}\psi_n=0$ (13.56) meet both
     boundary conditions with n sign changes: **c₁ = 3.61 m/s** (uniform N = 2.7 × 10⁻³ s⁻¹, H = 4200 m)."
  2. *The uniform-N check* — "$c_n=\frac{NH}{n\pi}$ (13.71) = 2.7 × 10⁻³ × 4200/π = **3.610 m/s**; the free surface lowers it by 3 parts
     in 10⁴ (the exact root of $\tan\frac{NH}{c_n}=\frac{c_nN}g$ (13.69) is 3.6085)."
  3. *Equivalent depth* — "$c_n^2\equiv gH_e$ (13.62): H_e = 3.61²/9.81 = **1.33 m**. This mode moves horizontally exactly like
     a homogeneous layer 1.33 m deep."
  4. *Its Rossby radius* — "radius Λ₁ = c₁/abs(f) = 3.61/8.37 × 10⁻⁵ = **43 km** at 35°."
  5. *A quick estimate* — "(ours) WKB: c_n ≈ ∫N dz/(nπ). For uniform N it is exact; for the default thermocline it gives
     1.30 against the true 1.77 m/s — the thermocline is too sharp to count as 'slowly varying'."
  6. *The wave view* (or the hint "turn the phone sideways to see the mode travel").
  7. *Reading the current setting* — n = 0: "All levels move together; the whole depth acts: c ≈ √(gH) ≈ 203 m/s. This is
     the tide and the tsunami." · n ≥ 1, free surface: "The upper ocean moves against the deep ocean; the surface barely
     moves (by about N²H/g ≈ 0.3 % of the interface). Mode n has n levels of no horizontal motion." · thermocline: "The
     zero crossing sits in the thermocline; a deeper or stronger thermocline is a faster mode — that is the signal El
     Niño forecasts track." · lid on: "With a lid the surface cannot rise at all: the barotropic mode is gone and the
     internal speeds change by less than 0.1 %."
- **Derivation tab:** **D11** (13 steps) `view: 'shape'`; goal `set` {prof uniform, n 1}; step 2 `watch` "teal is the running
  integral of purple: it is largest where purple crosses zero"; step 5 `live` "(ψ′/N²)/∫ψ at the cursor = −1/c² =
  `Viz.live('sep')` s²/m²"; step 6 `set` {n: 2}; step 10 `live` "c₁² = gH_e → H_e = `Viz.live('He')` m"; step 13 `watch` "the
  overlap integral of modes 1 and 2 (weight 1) printed in the title: 0 — with the lid on or off".
- **Code:**
  ```python
  z = np.linspace(-4200.0, 0.0, 401)                         # depth nodes [m]
  N2 = ch13.thermocline_N2(z, z_t={{zt}}, width={{wd}})       # or np.full_like(z, 2.7e-3**2) for uniform N
  modes = ch13.vertical_modes(z, N2, n_modes=5, lid="{{lid}}")   # solves (13.56) with (13.64), (13.65)
  c = modes[0][{{idx}}]                                      # c_n = {{c}} m/s (index n; n - 1 with a rigid lid)
  He = ch13.equivalent_depth(c)                              # (13.62): {{He}} m
  Lam = ch13.rossby_radius(c, ch13.coriolis_parameter(np.deg2rad({{lat}})))   # {{Lam}} km
  ```
- **Walkthrough (6 steps):** 1. "One layer?" — "The ocean's density changes all the way down. Can a single layer describe
  it? Step the mode number." `controls: ['n']` · 2. "Shapes" — "Each mode has a fixed vertical shape. Mode n crosses zero n
  times." `set` {prof uniform, n 2} · 3. "An eigenvalue problem" — "The shape and its speed come together: only special
  speeds fit the top and bottom." `eq: 'sl'`, `derive: {id: 'D11', step: 6}` · 4. "A metre deep" — "Mode 1 travels at 3.6 m/s,
  like a layer only 1.33 m deep." `readouts: ['c', 'He']`, `eq: 'he'` · 5. "The thermocline decides" — "Switch to the
  thermocline and deepen it: the mode speeds up." `set` {prof thermocline}, `controls: ['zt']`, `code: {id: 'vm', lines: [2,
  4]}` · 6. "Your turn" — "Predict what the rigid lid does to modes 0 and 1. Then toggle it." `controls: ['lid', 'n']`.
- **Equations:** `sl` "Vertical structure" ref 'Eq. (13.56)' $\frac d{dz}\big(\frac1{N^2}\frac{d\psi_n}{dz}\big)+\frac1{c_n^2}\psi_n=0$, live at the cursor · `bc`
  "Bottom and surface" ref 'Eq. (13.64), (13.65)' $\frac{d\psi_n}{dz}=0$ at $z=-H$; $\frac{d\psi_n}{dz}+\frac{N^2}g\psi_n=0$ at $z=0$ · `he` "Equivalent depth" ref
  'Eq. (13.62)' $c_n^2\equiv gH_e$, live · `tan` "Uniform N" ref 'Eq. (13.69)' $\tan\frac{NH}{c_n}=\frac{c_nN}g$, live with the root · `cn` "Baroclinic
  speeds" ref 'Eq. (13.71)' $c_n=\frac{NH}{n\pi}$, n = 1, 2, 3, …
- **Check yourself:** (1) "How many times does mode 3 cross zero?" — "Three." `set {n: 3}` · (2) "Deepen the thermocline
  from 150 m to 600 m. What happens to c₁?" — "It doubles, 1.24 → 2.47 m/s: more of the stratified water lies above the
  node." `set {zt: 600}` · (3) "Toggle the rigid lid. Which rung of the ladder disappears, and do the remaining shapes stay orthogonal?" — "Only the
  barotropic one; and yes — the plain overlap ∫ψ_mψ_n dz is zero for both lids." ·
  (4) "Why is H_e so small?" — "Because the restoring force is buoyancy, weaker than gravity by Δρ/ρ ~ 10⁻³."
- **Selftest parity rows:** `{name: 'c1 uniform', js: shootMode(1, uni, 'free').c, py: 'ch13.modes_uniform_N(2.7e-3, 4200.0)[0][1]',
  rtol: 1e-6}` (3.60849) · `{name: 'c0 uniform', js: shootMode(0, uni, 'free').c, py: 'ch13.modes_uniform_N(2.7e-3, 4200.0)[0][0]',
  rtol: 1e-6}` (203.054) · `{name: 'c1 thermocline vs matrix', js: shootMode(1, thermo, 'free').c, py:
  'ch13.vertical_modes(np.linspace(-4200.0, 0.0, 801), ch13.thermocline_N2(np.linspace(-4200.0, 0.0, 801)))[0][1]', rtol: 1e-4}` (1.7718)
  · `{name: 'c1 thermocline vs shooting', js: shootMode(1, thermo, 'free').c, py: 'ch13.vertical_modes_shooting(np.linspace(-4200.0,
  0.0, 401), ch13.thermocline_N2, n=1)[0]', rtol: 1e-8}` · `{name: 'rigid', js: cRigid(2.7e-3, 4200, 1), py:
  'ch13.baroclinic_mode_speed(2.7e-3, 4200.0, 1)', rtol: 1e-12}` (3.60963) · `{name: 'He', js: He(3.6096), py:
  'ch13.equivalent_depth(3.6096)', rtol: 1e-12}` · `{name: 'radius', js: Lam(3.6096, fCor(35, 'N')), py: 'ch13.rossby_radius(3.6096,
  ch13.coriolis_parameter(np.deg2rad(35.0)))', rtol: 1e-12}`.
- **Fit plan:** 360×640: status · `shape` (60 %) over `ladder` (40 %); `wave` hidden (c_n, H_e and Λ_n stay in the ladder's
  title and the readouts); transport hidden on phones (the wave view is its only client). 844×345: `shape` beside
  `ladder`. Desktop: rows [1.3, 1].


### E6 · shallow_water_dispersion
- **Title:** "Poincaré, Kelvin, Rossby — separate theories?" · **Summary:** "Move a wavenumber cursor over one diagram and watch
  the four terms of one cubic rebalance for the fast and the slow waves." · **CORE:** C09, C10 (also N82, N83, N87, N90, N137
  and the Kelvin line of C11) · **Reference:** `amplitude_phase_second_order_II_3.html` (response curves with a crosshair; a
  numbered Explain).
- **meta:** `viz:order 6` · `viz:sections 13.10 13.11 13.12 13.15` · `viz:equations 13.75 13.76 13.82 13.86 13.118` · `viz:fluidpy
  ch13.shallow_water_omega ch13.shallow_water_branches ch13.dispersion_term_sizes ch13.poincare_omega ch13.poincare_orbit
  ch13.rossby_omega ch13.kelvin_omega ch13.shallow_water_regime ch13.shallow_water_discriminant ch13.beta_parameter` ·
  `viz:derivations D12 D13`.
- **Physics (JS ↔ Python):** `cubicRoots(k, l, c, f0, beta)` → [ω₋, ω_R, ω₊] by the trigonometric form (P319) ↔
  `ch13.shallow_water_omega`; `termSizes(k, l, c, f0, beta, w)` → {omega3, gravity, rotation, beta, sum} ↔
  `ch13.dispersion_term_sizes`; `poincare(K, f, c)` ↔ `ch13.poincare_omega`; `kelvin(k, c)` ↔ `ch13.kelvin_omega`; `rossby(k, l, beta,
  f0, c)` ↔ `ch13.rossby_omega`; `orbit(t, k, eta, H, f)` ↔ `ch13.poincare_orbit(...)["u"]`, `["v"]`, `["axis_ratio"]`; `regime(w, f0)`
  ↔ `ch13.shallow_water_regime(...)["regime"]`; `disc(k, l, c, f0, beta)` ↔ `ch13.shallow_water_discriminant`; `betaOf(lat)` ↔
  `ch13.beta_parameter`.
- **Views** (rows [1.3, 1]): 1. `disp` "Frequency against wavenumber" (row 0, flex 1.4): log–log, k from 10⁻⁸ to 10⁻³ m⁻¹
  (wavelength labels on the top axis), abs(ω) from 10⁻⁹ to 10⁻² s⁻¹; the two fast branches (purple, they coincide to the eye;
  the legend says "±"), the slow branch (amber), the Kelvin line (dashed teal), ω = abs(f₀) (muted horizontal), the
  inertial-period and one-year levels labelled; the cursor's three roots as dots, the selected one ringed. **No-real-roots band (binding for the builder):** `ch13.shallow_water_omega` returns (NaN, NaN, NaN) where
  `ch13.shallow_water_discriminant` < 0 (possible only at low latitude with the external speed; at 12° N every wavelength
  above 12 500 km). The JS mirror `cubicRoots` tests `disc(...)` **first** and returns `null` there — it never computes or
  prints a NaN. In that band the three branches are not drawn: each ends with a ◆ at the band's edge, the band is
  shaded muted with the label "no three real roots here — the β-plane is stretched too far", the cursor's dots and the
  term bars are replaced by the same sentence, and the readouts show "—". Pointer: drag the cursor in k. 2. `terms` "The four terms of the cubic" (row 0,
  flex 0.8): signed bars ω³ (purple), −c²K²ω (orange), −f₀²ω (teal), −c²βk (amber) for the selected root, on a
  symmetric-log axis, with their sum (≈ 0) printed. 3. `plan` "The selected wave from above" (row 1, `hidePortrait`): crests
  moving (transport), one parcel with its current ellipse and sense of rotation; for the slow root, the crests drift
  west.
- **Controls:** `lam` (cursor) "Wavelength" (10…600 000 km, log; default 3100 km) · `ln` "North–south wavenumber $l/k$" (0…3,
  default 0; optional) · `lat` "Latitude" (5…80°, default 35) with `hemi` chips N / S (optional; only the ellipse's sense
  uses it — the cubic contains f₀²) · `mode` chips "external √(gH) / first baroclinic" (c = 202.95 or 3.610 m/s) · `root`
  chips "fast + / fast − / slow" · `betaOn` toggle "β" (default on) · transport `t`.
- **Transport:** `t` 0…3 periods of the selected root, `end: 'loop'`; rate chosen per root so that one period takes 4 s.
- **Presets:** "no rotation (f → 0)" {lat 0.01} · "inertial limit (very long wave)" {lam 300 000, root fast +} · "K = 1/Λ:
  slow branch peaks" {lam 2πΛ, root slow} · "external vs baroclinic" (toggles `mode`).
- **Status** (built from `regime`): in the no-real-roots band: "⚠️ no three real roots at this wavelength (discriminant below zero): shorten the wave or move poleward"; otherwise fast root: "fast root: ω = 5.0 f — a gravity wave bent by rotation; the β term is 10⁻³
  of the others"; ω within 5 % of f: "near-inertial: the water circles, the surface hardly moves"; slow root: "slow root: ω
  = 0.11 f ≪ f — the ω³ term is negligible; this wave exists only because of β"; β off: "with β = 0 the slow root is ω = 0, a
  steady geostrophic flow".
- **Readouts:** "$\omega/f_0$" · "Period" (`Viz.fmtTime`) · "Phase speed" · "Slow / fast" · "Discriminant sign".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **terms** · **presets** · **status** · **transport**
  · a `printed / corrected` toggle for slip #7 inside the regime text of the Explain tab.
- **Explain:**
  0. *What the views show* — "<b class=c-accent>purple</b> = the two fast roots (gravity waves with rotation); <b
     class=c-amber>amber</b> = the slow root (the planetary wave); dashed <b class=c-teal>teal</b> = the Kelvin wave of a
     coast, which is not a root of the cubic (it needs a wall). Bars: the four terms of the cubic for the root you
     selected."
  1. *The parameters* — "Coriolis parameter f₀ = 2Ω sin 35° = 8.37 × 10⁻⁵ s⁻¹; β = 2Ω cos 35°/R = 1.875 × 10⁻¹¹ m⁻¹ s⁻¹; c = √(gH) = 202.9 m/s; k = 2π/3100
     km = 2.03 × 10⁻⁶ m⁻¹."
  2. *The three roots* — "$\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76) has the roots **−4.15 × 10⁻⁴, −8.89 × 10⁻⁶, +4.24 × 10⁻⁴ s⁻¹**: periods
     4.1 h (fast) and 8.2 days (slow)."
  3. *The fast roots without β* — "$\omega^2=f^2+gHK^2$ (13.82): √((8.37 × 10⁻⁵)² + (202.9 × 2.03 × 10⁻⁶)²) = **4.20 × 10⁻⁴ s⁻¹** — the
     average of the two fast magnitudes; β splits them by ± half the slow root."
  4. *The slow root without ω³* — "$\omega=-\frac{\beta k}{k^2+l^2+f_0^2/c^2}$ (13.118) = **−8.884 × 10⁻⁶ s⁻¹**; the full cubic gives −8.888 × 10⁻⁶: the
     dropped term changes it by 0.04 %."
  5. *Term sizes for the selected root* — the four bars as numbers, with the largest pair named.
  6. *Are all three real?* — "discriminant $4(c^2K^2+f_0^2)^3-27(c^2\beta k)^2$ > 0 here. It is positive for every wavelength
     when βc < f₀² (ours, computed): 0.54 at this setting." In the no-real-roots band this section reads instead: "the
     discriminant is negative: the cubic has one real root and a complex pair, which only means that f can no longer
     be frozen at f₀ over such a long wave. Nothing is plotted."
  7. *The plan view* (or the hint "turn the phone sideways to see the wave from above").
  8. *Reading the current setting* — fast, ω > 3f: "A gravity wave that hardly notices rotation: ω ≈ cK." · fast, f < ω <
     3f: "Rotation matters: the frequency cannot go below f and the current traces an ellipse with axis ratio ω/f." ·
     slow: "Three decades below the fast waves. Drop ω³ and you have the Rossby wave. ⚠️ slip #7 — the book prints 'ω ≫ f'
     for the range in which the first term is negligible; the correct form is ω ≪ f [toggle shows the printed sentence's
     claim against the bars]." · β off: "Without β the slow root sits at ω = 0: a steady geostrophic flow that neither
     propagates nor decays — what geostrophic adjustment leaves behind."
- **Derivation tab:** **D12** (10 steps) `view: 'terms'` (phones: `'disp'`); goal `set` {root slow}; step 6 `watch` "the amber
  bar is the βv term of the vorticity equation"; step 9 `live` the four operators evaluated for the cursor wave. **D13** (7
  steps) `view: 'disp'`; step 1 `live` "(−iω)³ … → ω³ − c²K²ω − f₀²ω − c²βk = `Viz.live('sum')`"; step 5 `set` {lat 12, mode
  external}, `watch` "the ◆ marks where the discriminant turns negative: the β-plane has been stretched too far"; step 6
  `set` {root fast +}; step 7 `set` {root slow}.
- **Code:**
  ```python
  f0 = ch13.coriolis_parameter(np.deg2rad({{lat}})); beta = ch13.beta_parameter(np.deg2rad({{lat}}))
  k = 2*np.pi/{{lam}}                                        # wavenumber [1/m]
  w = ch13.shallow_water_omega(k, {{l}}, {{c}}, f0, beta)    # (13.76): {{w0}}, {{w1}}, {{w2}} 1/s
  terms = ch13.dispersion_term_sizes(k, {{l}}, {{c}}, f0, beta, w[{{i}}])   # the four bars; sum = {{sum}}
  fast = ch13.poincare_omega(k, f0, {{c}})                   # (13.82): {{fast}} 1/s
  slow = ch13.rossby_omega(k, {{l}}, beta, f0, {{c}})        # (13.118): {{slow}} 1/s
  ```
- **Walkthrough (7 steps):** 1. "Three names" — "Poincaré, Kelvin, Rossby: three theories, or one? Drag the cursor along
  the diagram." · 2. "One cubic" — "Eliminate u and η and a plane wave must satisfy one cubic in ω. It has three roots."
  `eq: 'cubic'`, `derive: {id: 'D13', step: 1}` · 3. "The fast pair" — "Two roots are gravity waves that cannot go below f.
  Their β bar is invisible." `set` {root: 'fast +'}, `terms: true` · 4. "The slow one" — "The third root is a hundred to a
  thousand times slower. Its ω³ bar is the one that vanishes." `set` {root: 'slow'}, `terms: true` · 5. "Switch β off" — "The
  slow root drops to zero: without β there is no planetary wave." `set` {betaOn: false} · 6. "No rotation" — "Slide f to
  zero: the hyperbola collapses onto ω = cK." `set` {lat: 0.01}, `code: {id: 'sw', lines: [3, 5]}` · 7. "Your turn" — "Predict
  the period of the slow wave for the first baroclinic mode. Then switch the mode." `controls: ['mode', 'lam']`.
- **Equations:** `v` "One equation for v" ref 'Eq. (13.75)' $\frac{\partial^3v}{\partial t^3}-gH\frac\partial{\partial t}\nabla_H^2v+f_0^2\frac{\partial v}{\partial t}-gH\beta\frac{\partial v}{\partial x}=0$ · `cubic` "Complete dispersion
  relation" ref 'Eq. (13.76)' $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$, live with the four terms · `poin` "Poincaré" ref 'Eq. (13.82)'
  $\omega^2=f^2+gHK^2$, live · `kel` "Kelvin" ref 'Eq. (13.86)' $c=\sqrt{gH}$, i.e. $\omega=\pm k\sqrt{gH}$ · `ros` "Rossby" ref 'Eq. (13.118)'
  $\omega=-\frac{\beta k}{k^2+l^2+f_0^2/c^2}$, live.
- **Check yourself:** (1) "For which root is the ω³ term the smallest of the four?" — "The slow one (ω ≪ f)." `set {root:
  'slow'}` · (2) "What is the lowest frequency a fast wave can have at 35°?" — "f: 8.37 × 10⁻⁵ s⁻¹, a period of 20.9 h." `set
  {lam: 300000}` · (3) "Switch from external to baroclinic. By what factor does the slow wave's period change at 3100
  km?" — "From 8.2 days to 2.8 years, about 125 times." · (4) "Why is the Kelvin line dashed?" — "It is not a root of this
  cubic: it needs a coast (C11)."
- **Selftest parity rows:** `{name: 'slow root', js: cubicRoots(2*Math.PI/3.1e6, 0, 202.95, fCor(35, 'N'), betaOf(35))[1], py:
  'ch13.shallow_water_omega(2*np.pi/3.1e6, 0.0, 202.95, ch13.coriolis_parameter(np.deg2rad(35.0)), ch13.beta_parameter(np.deg2rad(35.0)))[1]',
  rtol: 1e-9}` (−8.8883e-6) · the same with index [2] (4.24135e-4) and [0] · `{name: 'baroclinic slow', js: cubicRoots(…, 3.6096,
  …)[1], py: 'ch13.shallow_water_omega(2*np.pi/3.1e6, 0.0, 3.6096, ch13.coriolis_parameter(np.deg2rad(35.0)),
  ch13.beta_parameter(np.deg2rad(35.0)))[1]', rtol: 1e-9}` (−7.0231e-8) · `{name: 'beta term', js: termSizes(…).beta, py:
  'ch13.dispersion_term_sizes(2*np.pi/3.1e6, 0.0, 202.95, ch13.coriolis_parameter(np.deg2rad(35.0)),
  ch13.beta_parameter(np.deg2rad(35.0)), 4.0e-4)["beta"]', rtol: 1e-12}` · `{name: 'poincare', js: poincare(2*Math.PI/3.1e6, fCor(35,
  'N'), 202.95), py: 'ch13.poincare_omega(2*np.pi/3.1e6, ch13.coriolis_parameter(np.deg2rad(35.0)), 202.95)', rtol: 1e-12}`
  (4.19762e-4) · `{name: 'beta', js: betaOf(35), py: 'ch13.beta_parameter(np.deg2rad(35.0))', rtol: 1e-12}` · `{name: 'axis ratio S',
  js: orbit(0, 2*Math.PI/3.1e6, 0.5, 4200, fCor(35, 'S')).axis_ratio, py: 'ch13.poincare_orbit(0.0, 2*np.pi/3.1e6, 0.5, 4200.0,
  ch13.coriolis_parameter(np.deg2rad(-35.0)))["axis_ratio"]', rtol: 1e-10}` · `{name: 'discriminant sign 12N', js: disc(2*Math.PI/2.0e7, 0, 202.95, fCor(12, 'N'), betaOf(12)) < 0 ? 1 : 0, py:
  '1.0*(ch13.shallow_water_discriminant(2*np.pi/2.0e7, 0.0, 202.95, ch13.coriolis_parameter(np.deg2rad(12.0)),
  ch13.beta_parameter(np.deg2rad(12.0))) < 0)', rtol: 0}` (1; the Python roots are NaN there, so the parity row compares the
  flag, never the roots) · invariant `{name: 'roots sum to 0', js: sum3(roots)/roots[2],
  expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: status · `disp` (58 %) over `terms` (42 %); `plan` hidden (the selected root's period and ω/f are
  in the diagram's title); transport hidden on phones. 844×345: `disp` beside `terms`. Desktop: rows [1.3, 1].

### E7 · kelvin_wave
- **Title:** "Why does this wave run only one way along a coast?" · **Summary:** "Reverse the wave or the hemisphere and watch
  the offshore profile turn from decaying to growing — and be rejected." · **CORE:** C11 (also N91, N94, R23, R24, N96 and the
  trapping scale of C12) · **Reference:** `angular_frequency_explorer_1.html`.
- **meta:** `viz:order 7` · `viz:sections 13.12` · `viz:equations 13.84 13.86 13.87` · `viz:fluidpy ch13.kelvin_wave ch13.kelvin_residuals
  ch13.kelvin_decay_side ch13.long_wave_speed ch13.rossby_radius ch13.rossby_radius_two_layer ch13.kelvin_omega` · `viz:derivations
  D15`.
- **Physics (JS ↔ Python):** `kelvin(x, y, t, eta0, k, H, f, dir)` → [η, u] ↔ `ch13.kelvin_wave(..., direction=dir)`; `decaySide(f,
  dir)` → {trapped, coast_on} ↔ `ch13.kelvin_decay_side`; `cLong(H)` ↔ `ch13.long_wave_speed`; `Lam(c, f)` ↔ `ch13.rossby_radius`;
  `LamTwoLayer(H1, r1, r2, f)` ↔ `ch13.rossby_radius_two_layer`; `omegaK(k, c)` ↔ `ch13.kelvin_omega`; `residuals(...)` (analytic
  derivatives) ↔ `ch13.kelvin_residuals`.
- **Views** (rows [1.3, 1]): 1. `plan` "Surface height from above" (row 0, flex 1.4): x along the coast (3 wavelengths), y
  offshore 0…4Λ; η as a heat map (orange high, blue low), current arrows u (teal), the coast as a thick line at y = 0
  (and a second wall at y = W in channel mode); an arrow "direction of travel" and a label "coast on the right" /
  "coast on the left". Crests move with the transport. 2. `sect` "Across the shore" (row 0, flex 1): η(y) through the
  crest (bold) and the trough (thin), the envelope $\pm\eta_0e^{-y/\Lambda}$ dashed, Λ marked on the axis; in the rejected case
  the growing profile is drawn in muted rose up to the frame edge with a ✕ and the label "grows offshore". 3. `bal`
  "Cross-shore balance at the cursor" (row 1, `hidePortrait`): two bars, Coriolis on the along-shore current $fu$ (teal)
  and the surface slope $-g\,\partial\eta/\partial y$ (orange), equal; a third thin bar for ∂v/∂t = 0.
- **Controls:** `mode` chips "external (H) / internal two-layer" · `H` "Depth $H$" (100…6000 m, default 4200; internal mode:
  upper layer 120 m, Δρ = 3.1 kg m⁻³) · `lat` "Latitude" (5…80°, default 35; help: "5°: Λ is huge · 80°: Λ is smallest") with
  `hemi` chips N / S · `lam` "Wavelength" (0.3…10 Λ, default 1.28 Λ = 3100 km) · `geo` chips "one coast / channel" (channel
  width 2Λ) · `dir` chips "travel → / ←" · transport `t`.
- **Transport:** `t` 0…2 periods, loop.
- **Presets:** "deep-ocean tide" {external, H 4200, lat 35} · "internal wave on the thermocline" {internal} · "channel 2Λ
  wide" {geo channel} · "35° S" {hemi S}.
- **Status:** `decaySide` true: "✅ trapped: coast on the right of the direction of travel · Λ = 2426 km · c = 203 m/s" (N) /
  "… coast on the left …" (S); false: "✕ this direction would grow offshore — not a solution in this hemisphere".
- **Readouts:** "Speed $c$" · "Rossby radius $\Lambda$" · "Period" · "Current at the coast" · "$fu$ vs $g\eta_y$".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** · **modes**
  (external / internal).
- **Explain:**
  0. *What the views show* — "Plan: orange = raised surface, blue = lowered; <b class=c-teal>teal</b> arrows = the
     along-shore current. Section: the surface along a line running offshore. Bars: the two forces that balance across
     the shore."
  1. *Speed* — "$c=\sqrt{gH}$ (13.86) = √(9.81 × 4200) = **203 m/s** — the same as without rotation, for every frequency."
  2. *Trapping width* — "$\Lambda\equiv c/f$ = 203/8.37 × 10⁻⁵ = **2426 km**: the height falls by e every Λ offshore."
  3. *Period* — "wavelength/c = 3100 km/203 m/s = **4.2 h**."
  4. *Current under the crest* — "$u=\eta_0\sqrt{g/H}\,e^{-fy/c}\cos k(x-ct)$ (the second of (13.87)): at the coast 0.5 × √(9.81/4200) =
     **0.024 m/s**."
  5. *The balance across the shore* — "$fu=-g\frac{\partial\eta}{\partial y}$ (the third of (13.84)): 8.37 × 10⁻⁵ × 0.024 = 2.0 × 10⁻⁶ m/s² = g η₀/Λ ✓."
  6. *The bar view* (or the hint "turn the phone sideways for the force bars").
  7. *Reading the current setting* — trapped, N: "Under a crest the current runs with the wave; the Coriolis force on it
     points to the right, toward the coast, and is held by the surface sloping down offshore. So the coast must be on
     the right." · S: "With f < 0 the Coriolis force points left, and so must the coast." · rejected: "For this direction the
     same balance needs a surface that *rises* offshore without limit. In a half-plane that is not a solution; in a
     channel it is the wave on the opposite wall." · internal: "The thermocline's Kelvin wave is a hundred times slower
     (1.9 m/s) and a hundred times narrower (Λ = 23 km): it hugs the coast — the poleward signal after coastal upwelling
     or an El Niño." · channel: "Two walls, two waves, each with its own wall on its right; where both are present the
     tide rotates about a point of no rise."
- **Derivation tab:** **D15** (8 steps) `view: 'sect'`; goal `set` {preset deep-ocean tide}; step 4 `live` "ω² = gHk² → c =
  `Viz.live('c')` m/s"; step 6 `set` {dir: '←'}, `watch` "the profile grows offshore: the status rejects it"; step 7 `set`
  {dir: '→'}, `live` "Λ = c/f = `Viz.live('Lam')` km"; step 8 `watch` "the two bars of the balance view are equal at every
  y".
- **Code:**
  ```python
  f = ch13.coriolis_parameter(np.deg2rad({{lat}}))             # {{f}} 1/s
  c = ch13.long_wave_speed({{H}})                              # (13.86): {{c}} m/s
  Lam = ch13.rossby_radius(c, f)                               # trapping width: {{Lam}} km
  side = ch13.kelvin_decay_side(f, {{dir}})                    # trapped: {{trapped}}; coast on the {{coast}}
  eta, u = ch13.kelvin_wave(x, y, t, 0.5, {{k}}, {{H}}, f, direction={{dir}})   # (13.87)
  ```
- **Walkthrough (6 steps):** 1. "Leaning on the coast" — "This wave travels along a wall and fades offshore. Press ▶." `play:
  true` · 2. "No flow across the shore" — "With v = 0 the Coriolis force on the along-shore current must be held by a
  slope." `eq: 'eqs'`, `derive: {id: 'D15', step: 1}` · 3. "The same speed" — "It moves at √(gH) = 203 m/s whatever its
  frequency — even below f." `readouts: ['c']`, `eq: 'c'` · 4. "How far offshore" — "The slope decays over Λ = c/f = 2426 km."
  `eq: 'sol'`, `code: {id: 'kw', lines: [2, 3]}` · 5. "Wrong way" — "Reverse the direction: the profile would grow offshore.
  Not allowed." `set` {dir: '←'}, `derive: {id: 'D15', step: 6}` · 6. "Your turn" — "Predict the direction of travel at 35° S
  with the same coast. Then switch." `controls: ['hemi', 'dir', 'mode']`.
- **Equations:** `eqs` "Kelvin-wave equations" ref 'Eq. (13.84)' $\frac{\partial\eta}{\partial t}+H\frac{\partial u}{\partial x}=0,\ \frac{\partial u}{\partial t}=-g\frac{\partial\eta}{\partial x},\ fu=-g\frac{\partial\eta}{\partial y}$, live residuals ·
  `c` "Speed" ref 'Eq. (13.86)' $c=\sqrt{gH}$, live · `sol` "The wave" ref 'Eq. (13.87)' $\eta=\eta_0e^{-fy/c}\cos k(x-ct),\ u=\eta_0\sqrt{g/H}\,e^{-fy/c}\cos k(x-ct)$,
  live, note "as printed: f > 0, travel toward +x" · `lam` "Rossby radius" ref 'p. 657 (no number)' $\Lambda\equiv c/f$, live.
- **Check yourself:** (1) "Halve the depth. What happens to c and to Λ?" — "Both fall by √2." `set {H: 2100}` · (2) "At 35° S,
  with the coast at y = 0 and the sea at y > 0, which way does the trapped wave travel?" — "Toward −x: the coast must be
  on its left." `set {hemi: 'S'}` · (3) "How wide is the internal wave?" — "23 km against 2426 km." `set {mode: 'internal'}` ·
  (4) "Can a Kelvin wave have a period longer than the inertial period?" — "Yes — any period; that is what
  distinguishes it from a Poincaré wave."
- **Selftest parity rows:** `{name: 'eta N', js: kelvin(1e5, 1e6, 600, 0.5, 2*Math.PI/3.1e6, 4200, fCor(35, 'N'), 1)[0], py:
  'ch13.kelvin_wave(1.0e5, 1.0e6, 600.0, 0.5, 2*np.pi/3.1e6, 4200.0, ch13.coriolis_parameter(np.deg2rad(35.0)), direction=1)[0]', rtol:
  1e-10}` · `{name: 'u S', js: kelvin(1e5, 1e6, 600, 0.5, 2*Math.PI/3.1e6, 4200, fCor(35, 'S'), -1)[1], py: 'ch13.kelvin_wave(1.0e5,
  1.0e6, 600.0, 0.5, 2*np.pi/3.1e6, 4200.0, ch13.coriolis_parameter(np.deg2rad(-35.0)), direction=-1)[1]', rtol: 1e-10}` · `{name:
  'radius', js: Lam(cLong(4200), fCor(35, 'N')), py: 'ch13.rossby_radius(ch13.long_wave_speed(4200.0),
  ch13.coriolis_parameter(np.deg2rad(35.0)))', rtol: 1e-12}` (2.42611e6) · `{name: 'two-layer radius', js: LamTwoLayer(120, 1023.9,
  1027, fCor(35, 'N')), py: 'ch13.rossby_radius_two_layer(120.0, 1023.9, 1027.0, ch13.coriolis_parameter(np.deg2rad(35.0)))', rtol:
  1e-12}` (2.2531e4) · `{name: 'trapped flag', js: decaySide(fCor(35, 'S'), 1).trapped ? 1 : 0, py:
  '1.0*ch13.kelvin_decay_side(ch13.coriolis_parameter(np.deg2rad(-35.0)), 1)["trapped"]', rtol: 0}` (0) · invariant `{name: 'geostrophy
  residual', js: residuals(…).y_geostrophy, expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: status · `plan` (58 %) over `sect` (42 %); `bal` hidden (Λ and c are in the section's title, the
  balance check is the readout "fu vs gη_y"). 844×345: `plan` beside `sect`. Desktop: rows [1.3, 1].

### E8 · geostrophic_adjustment
- **Title:** "Why doesn't a pile of water flatten out here?" · **Summary:** "Release a step in surface height with and without
  rotation: waves leave, and a front one Rossby radius wide stays." · **CORE:** C12 (also N19, the fast waves of C10, the
  conserved potential vorticity of R26) · **Reference:** `forced_damped_vibrations.html` (transient + steady state on one time
  slider).
- **meta:** `viz:order 8` · `viz:sections 13.5 13.8 13.12` · `viz:equations 13.45 13.82` · `viz:fluidpy ch13.geostrophic_adjustment_1d
  ch13.adjustment_energy ch13.linear_1d_run ch13.linear_1d_step ch13.rossby_radius ch13.poincare_omega ch13.long_wave_speed` ·
  `viz:derivations D16`.
- **Physics (JS ↔ Python):** `stepFB(eta, u, v, dx, dt, H, f)` (forward–backward on the staggered grid: η at 200 centres, u at
  201 faces with u = 0 at both ends, v at centres) ↔ `ch13.linear_1d_step` (one step) and `ch13.linear_1d_run(...)["eta"][-1]`
  (after a fixed number of steps); `endState(x, eta0, H, f)` → [η, v] ↔ `ch13.geostrophic_adjustment_1d`; `energy(eta0, H, f)` →
  {pe_released, ke_jet, radiated, ratio} ↔ `ch13.adjustment_energy`; `Lam(c, f)` ↔ `ch13.rossby_radius`; for the two bump
  shapes the end state is the Green's-function convolution $\eta_\infty(x)=\frac1{2\Lambda}\int\eta_{init}(x')e^{-\lvert x-x'\rvert/\Lambda}dx'$ (ours; labelled in the Explain tab and in the view title **"computed in the explainer; the step case is checked
  against Python"** — by the orchestrator's ruling there is no general-shape Python function; the step case is pinned by
  the row 'convolution = closed form'); `linPV(v, eta, dx, f, H)` ↔ `ch13.linear_1d_run(...)["pv"]`.
- **Views** (rows [1.2, 1]): 1. `eta` "Surface height η(x, t)" (row 0, full width): x from −25Λ to 25Λ in km; the initial
  shape (muted), the evolving surface (orange, bold), the closed-form end state (dashed ghost), Λ marked on both sides
  of the step; departing wave fronts labelled "±c t". 2. `jet` "Velocity along the step v(x, t)" (row 1, flex 1): evolving v
  (purple), end-state jet (dashed), its maximum gη₀/c labelled. 3. `budget` "Energy and potential vorticity" (row 1, flex
  0.8, `hidePortrait`): three bars — potential energy released (orange), kinetic energy of the jet (purple), energy
  carried off by waves (muted) — and, above them, the linear potential vorticity $\zeta-f\eta/H$ against x now (amber) on top
  of its initial curve (muted): they coincide.
- **Controls:** `lat` "Latitude" (0…80°, default 35; help: "0°: no rotation — nothing stays · 80°: the narrowest front")
  with `hemi` chips N / S · `mode` chips "first baroclinic (H_e = 1.33 m) / external (H = 4200 m)" · `shape` chips "step /
  narrow bump (0.3 Λ) / wide bump (5 Λ)" · transport `t` (in inertial periods; in units of Λ/c when f = 0).
- **Transport:** `t` 0…6 inertial periods, rate 0.5 period/s, `end: 'hold'`; end card (rotating): "waves gone: a front of
  width Λ = 43 km remains · jet 0.136 m/s · one third of the released energy stayed"; (f = 0): "everything radiated away:
  the surface is flat between the two fronts".
- **Presets:** "no rotation" {lat 0} · "mid-latitude, baroclinic" {lat 35, baroclinic, step} · "external mode" {external} ·
  "bump narrower than Λ" {shape narrow} · "bump wider than Λ" {shape wide}.
- **Status:** f = 0: "no rotation: gravity flattens everything"; rotating: "adjusted front: width Λ = 43 km · jet along +y
  (the high side on its right)" (N) / "… along −y (the high side on its left)" (S); while t < 1 period: "adjusting:
  Poincaré waves leaving at up to 3.6 m/s".
- **Readouts:** "Rossby radius $\Lambda$" · "Jet maximum" · "KE kept / PE released" · "$fv$ vs $g\eta_x$ at the step" · "Inertial
  period".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** (play / step / scrub, end card) ·
  **presets** · **status** · **terms** (two groups: the energy bars; and the cross-step balance $fv$ against $g\,\partial\eta/\partial x$
  converging to equal).
- **Explain:**
  0. *What the views show* — "<b class=c-orange>orange</b> = the surface; dashed = where it ends up; <b
     class=c-accent>purple</b> = the current along the step; <b class=c-amber>amber</b> = the linear potential
     vorticity, which never changes."
  1. *The Rossby radius* — "$\Lambda\equiv c/f$ with c = √(gH_e) = √(9.81 × 1.329) = 3.61 m/s: Λ = 3.61/8.37 × 10⁻⁵ = **43 km**."
  2. *What cannot change* — "(ours, from the three shallow-water equations (a) $\frac{\partial\eta}{\partial t}+H(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y})=0$, (b)
     $\frac{\partial u}{\partial t}-fv=-g\frac{\partial\eta}{\partial x}$, (c) $\frac{\partial v}{\partial t}+fu=-g\frac{\partial\eta}{\partial y}$ of (13.45)) $\frac\partial{\partial t}\big(\zeta-\frac{f\eta}H\big)=0$ at every x."
  3. *The end state* — "surface η = η₀ sgn(x)(1 − e^{−abs(x)/Λ}); at x = Λ: 0.05 × (1 − e⁻¹) = **0.032 m**."
  4. *The jet* — "jet v = (gη₀/c) e^{−abs(x)/Λ}: at the step 9.81 × 0.05/3.61 = **0.136 m/s**, falling to 0.050 at x = Λ."
  5. *The energy split* — "released potential energy (3/2)ρgη₀²Λ; kept as kinetic energy (1/2)ρgη₀²Λ: **one third**; two
     thirds left with the waves (computed for your shape — for a bump the fraction differs and is shown on the bar)."
  6. *The budget view* (or the hint "turn the phone sideways for the energy bars").
  7. *At the current time* — "now t = `Viz.live('t')` inertial periods: slope force gη_x = `Viz.live('pgf')`, Coriolis fv =
     `Viz.live('cor')` at the step."
  8. *Reading the current setting* — f = 0: "Two fronts run off at ±c and leave a flat surface: all the potential energy
     is carried away." · rotating step: "The water starts toward the low side, is turned by the Coriolis force, and ends
     as a jet *along* the step whose Coriolis force holds the remaining slope. Only the water within about a Rossby
     radius of the step ever moved across it." · narrow bump: "Narrower than Λ: gravity wins, most of the bump radiates
     away as waves, a faint balanced remnant stays." · wide bump: "Wider than Λ: rotation wins, the bump barely changes;
     only its edges adjust, each over one Λ." · external: "Here Λ = 2426 km: on the scale of a weather system the external
     mode adjusts almost as if there were no rotation."
- **Derivation tab:** **D16** (11 steps) `view: 'eta'`; goal `set` {preset mid-latitude baroclinic, t: 0}; step 3 `set` {t: 2},
  `watch` "the amber curve has not moved"; step 6 `live` "η″ − η/Λ² = −η₀ sgn(x)/Λ² with Λ = `Viz.live('Lam')` km"; step 9
  `set` {t: 6}, `watch` "the orange curve settles onto the dashed one"; step 10 `live` "v(0) = gη₀/c = `Viz.live('v0')` m/s";
  step 11 `live` "KE/PE = `Viz.live('ratio')`".
- **Code:**
  ```python
  f = ch13.coriolis_parameter(np.deg2rad({{lat}}))                  # {{f}} 1/s
  Lam = ch13.rossby_radius(ch13.long_wave_speed({{H}}), f)          # {{Lam}} km
  out = ch13.linear_1d_run(eta0, {{dx}}, {{dt}}, {{n}}, {{H}}, f)   # our forward-backward march
  eta_end, v_end = ch13.geostrophic_adjustment_1d(x, {{eta0}}, {{H}}, f)   # ours: the end state
  en = ch13.adjustment_energy({{eta0}}, {{H}}, f)                   # ratio kept = {{ratio}}
  ```
- **Walkthrough (7 steps):** 1. "Let go" — "A step in the surface, released. First without rotation: press ▶." `set` {lat:
  0}, `play: true` · 2. "Now on a turning earth" — "The same step at 35°. Waves leave again — but something stays." `set`
  {lat: 35}, `play: true` · 3. "What stays" — "A front one Rossby radius wide, Λ = c/f = 43 km, and a current along it."
  `set` {t: 6}, `eq: 'lam'` · 4. "Why it stays" — "The jet's Coriolis force balances the slope: fv = g ∂η/∂x." `terms: true` ·
  5. "What pins it" — "At every point, ζ − fη/H keeps its starting value." `derive: {id: 'D16', step: 3}`, `eq: 'pv'` · 6.
  "The bill" — "One third of the released energy is kept in the jet; two thirds left as waves." `code: {id: 'adj', lines:
  [5, 5]}` · 7. "Your turn" — "Predict what survives of a bump much narrower than Λ. Then try it." `controls: ['shape',
  'lat']`.
- **Equations:** `sw` "Linear shallow water" ref 'Eq. (13.45)' $\frac{\partial\eta}{\partial t}+H\big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\big)=0,\ \frac{\partial u}{\partial t}-fv=-g\frac{\partial\eta}{\partial x},\ \frac{\partial v}{\partial t}+fu=-g\frac{\partial\eta}{\partial y}$ · `pv`
  "What is conserved (ours)" ref 'not in the book — ours' $\frac\partial{\partial t}\big(\zeta-\frac{f\eta}H\big)=0$, live · `lam` "Rossby radius" ref 'p. 657 (no
  number)' $\Lambda\equiv c/f$, live · `end` "End state (ours)" ref 'D16 — ours' $\eta=\eta_0\,\mathrm{sgn}(x)(1-e^{-\lvert x\rvert/\Lambda}),\ v=\frac{g\eta_0}c e^{-\lvert x\rvert/\Lambda}$, live · `poin`
  "The departing waves" ref 'Eq. (13.82)' $\omega^2=f^2+gHK^2$.
- **Check yourself:** (1) "Halve f. What happens to the width of the front and to the jet's peak speed?" — "The width
  doubles; the peak gη₀/c does not change." `set {lat: 16.7}` · (2) "What fraction of the released potential energy stays
  for a step?" — "One third." · (3) "Which survives better, a bump of width 0.3 Λ or 5 Λ?" — "The wide one." `set {shape:
  'wide'}` · (4) "In the southern hemisphere, which way does the jet flow?" — "The other way: the high side is on its
  left."
- **Selftest parity rows:** `{name: 'march', js: marchEta(50)[110], py: 'ch13.linear_1d_run(0.05*np.sign(np.arange(200.0)-99.5),
  1.0e4, 1385.0, 50, 1.3286, ch13.coriolis_parameter(np.deg2rad(35.0)))["eta"][-1][110]', rtol: 1e-9}` · `{name: 'march S v', js:
  marchV(50, 'S')[105], py: 'ch13.linear_1d_run(0.05*np.sign(np.arange(200.0)-99.5), 1.0e4, 1385.0, 50, 1.3286,
  ch13.coriolis_parameter(np.deg2rad(-35.0)))["v"][-1][105]', rtol: 1e-9}` · `{name: 'end eta', js: endState(4.315e4, 0.05, 1.3286,
  fCor(35, 'N'))[0], py: 'ch13.geostrophic_adjustment_1d(4.315e4, 0.05, 1.3286, ch13.coriolis_parameter(np.deg2rad(35.0)))[0]', rtol:
  1e-12}` (0.031606) · `{name: 'end v S', js: endState(0, 0.05, 1.3286, fCor(35, 'S'))[1], py: 'ch13.geostrophic_adjustment_1d(0.0,
  0.05, 1.3286, ch13.coriolis_parameter(np.deg2rad(-35.0)))[1]', rtol: 1e-12}` (−0.13584) · `{name: 'energy ratio', js: energy(0.05,
  1.3286, fCor(35, 'N')).ratio, py: 'ch13.adjustment_energy(0.05, 1.3286, ch13.coriolis_parameter(np.deg2rad(35.0)))["ratio"]', rtol:
  1e-12}` (1/3) · `{name: 'convolution = closed form', js: convEnd(stepShape)[150], py: 'ch13.geostrophic_adjustment_1d(5.05e5, 0.05,
  1.3286, ch13.coriolis_parameter(np.deg2rad(35.0)))[0]', rtol: 1e-3}`.
- **Fit plan:** 360×640: status · `eta` (55 %) over `jet` (45 %); `budget` hidden (the kept fraction is in the η view's
  title: "KE kept / PE released = 0.33"); transport with play and scrub only. 844×345: `eta` beside `jet`. Desktop: rows
  [1.2, 1].

### E9 · rossby_waves
- **Title:** "Crests go west — so how does the energy go east?" · **Summary:** "Watch a packet whose crests march west while
  its envelope moves east, and columns whose spin changes as they are carried north and south." · **CORE:** C15 (also N132,
  N134, N135, N136 and the recap R26) · **Reference:** `angular_frequency_explorer_1.html` (system + curve with a moving dot on
  one clock).
- **meta:** `viz:order 9` · `viz:sections 13.13 13.15` · `viz:equations 13.94 13.117 13.118 13.119 13.120` · `viz:fluidpy
  ch13.rossby_omega ch13.rossby_group_velocity ch13.rossby_phase_speed ch13.rossby_max_frequency ch13.rossby_omega_circle
  ch13.rossby_long_wave_speed ch13.stationary_rossby_wavelength ch13.qg_linear_evolve_1d ch13.rossby_packet_spectrum
  ch13.beta_parameter` · `viz:derivations D22 D23`.
- **Physics (JS ↔ Python):** `omega(k, l, beta, f0, c, U)` ↔ `ch13.rossby_omega`; `cg(k, l, beta, f0, c, U)` → [c_gx, c_gy] ↔
  `ch13.rossby_group_velocity`; `cx(...)` ↔ `ch13.rossby_phase_speed`; `wmax(beta, f0, c)` ↔ `ch13.rossby_max_frequency(...)["omega_max"]`;
  `circle(w, beta, f0, c)` ↔ `ch13.rossby_omega_circle`; `cLong(beta, f0, c)` ↔ `ch13.rossby_long_wave_speed`; `lamStat(U, beta)` ↔
  `ch13.stationary_rossby_wavelength`; `packetSpectrum(k0, sigma, n)` ↔ `ch13.rossby_packet_spectrum`; `evolve(a, k, l, t, …)`
  (each of ≤ 64 modes multiplied by $e^{-i\omega_jt}$) ↔ `ch13.qg_linear_evolve_1d`.
- **Views** (rows [1.2, 1]): 1. `plan` "The height field from above" (row 0, full width): x (east) over 12 packet widths,
  y (north) one wavelength of $e^{ily}$; η as a heat map with the packet's envelope outlined (dashed); a row of nine marked
  columns that move north–south with the wave, each coloured by its relative vorticity (amber = counter-clockwise, blue
  = clockwise) with a small ↺ / ↻; two tick marks below the map: ▲ "a crest" and ◆ "the envelope's centre", with the
  distances they have travelled. 2. `curve` "Frequency against wavenumber" (row 1, flex 1): $\omega f_0/(\beta c)$ against $kc/f_0$
  from −4 to 0 for the chosen l; the cursor dot, the chord from the origin (slope = phase speed), the tangent (slope =
  group velocity, coloured by its sign), the maximum marked "energy stands still". 3. `kplane` "Wavenumber plane" (row 1,
  flex 0.8, `hidePortrait`): the circle of constant ω through the cursor, its centre at k = −β/2ω, the group-velocity
  arrow at the cursor pointing inward.
- **Controls:** `lamx` (cursor) "East–west wavelength" (0.2…20 × 2πΛ, log, default 150 km) · `l` "North–south wavenumber $l\Lambda$"
  (0…3, default 0; optional) · `lat` "Latitude" (5…70°, default 35) · `rad` chips "first baroclinic (Λ = 43 km) /
  barotropic (Λ = 2426 km)" · `U` "Mean eastward flow" (0…30 m/s barotropic, 0…0.1 m/s baroclinic; default 0) · transport `t`
  (logarithmically paced: days to years).
- **Transport:** `t` from 0 to 5 periods of the cursor wave, log-paced (the ch12 pattern), `end: 'hold'`; end card: "the
  crest moved 5.0 wavelengths west; the envelope moved 2.7 wavelengths **east**".
- **Presets:** "long baroclinic wave at 12° N" {lat 12, lamx 2000 km} · "short barotropic wave" {barotropic, lamx 3100 km} ·
  "K = 1/Λ: energy stands still" {lamx 2πΛ} · "stationary wave in a westerly" {barotropic, U 17, lamx 5983 km}.
- **Status:** KΛ > 1: "phase west (−0.82 cm/s), energy **east** (+0.43 cm/s): shorter than 2πΛ"; KΛ < 1: "phase west,
  energy **west**: longer than 2πΛ"; KΛ = 1 (± 2 %): "energy stands still"; with U such that c_x = 0: "stationary: the
  current cancels the westward drift".
- **Readouts:** "Period" (`Viz.fmtTime`) · "Phase speed $c_x$" · "Group velocity $c_{gx}$" · "$2\pi\Lambda$" · "Basin crossing (10⁴
  km)".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** (log-paced) · **presets** · **status**
  · **inspector** (click a column: "displaced north by Y = … km → f larger by βY = … s⁻¹ → to keep (ζ + f)/h it must have
  ζ = −βY + (f₀/H)η = … s⁻¹").
- **Explain:**
  0. *What the views show* — "Map: the surface height of a group of waves; the dashed outline is the group (its
     energy). Columns: <b class=c-amber>amber</b> spins counter-clockwise, blue clockwise. Curve: every Rossby wave of
     this north–south scale; the dot is yours."
  1. *The parameters* — "Here β = 1.875 × 10⁻¹¹ m⁻¹ s⁻¹, f₀ = 8.37 × 10⁻⁵ s⁻¹, c = 3.61 m/s, so Λ = 43 km and 2πΛ = 271 km; your wave: 150
     km, k = −4.19 × 10⁻⁵ m⁻¹."
  2. *Frequency and period* — "$\omega=-\frac{\beta k}{k^2+l^2+f_0^2/c^2}$ (13.118) = 1.875 × 10⁻¹¹ × 4.19 × 10⁻⁵/(1.755 × 10⁻⁹ + 5.37 × 10⁻¹⁰) = **3.43 × 10⁻⁷
     s⁻¹** → period **212 days**."
  3. *Phase speed* — "$c_x=\frac\omega k=-\frac\beta{k^2+l^2+f_0^2/c^2}$ (13.119) = **−0.82 cm/s**: westward, as for every Rossby wave."
  4. *Group velocity* — "(ours) $c_{gx}=\frac{\beta(k^2-l^2-f_0^2/c^2)}{(k^2+l^2+f_0^2/c^2)^2}$ = **+0.43 cm/s**: eastward, because k² > f₀²/c²."
  5. *The turning point* — "The group velocity c_gx = 0 at kΛ = −1 (wavelength 2πΛ = 271 km), where ω is largest: ω_max = βc/(2f₀) = 4.05 × 10⁻⁷
     s⁻¹ — no Rossby wave of this mode has a period shorter than **180 days**."
  6. *Long waves and basins* — "long-wave speed −βΛ² = −3.49 cm/s; 10⁴ km in **9.1 years** (at 12° N: 1.0 year)."
  7. *The circle view* (or the hint "turn the phone sideways for the wavenumber plane") — "all waves of one frequency lie
     on $\big(k+\frac\beta{2\omega}\big)^2+l^2=\big(\frac\beta{2\omega}\big)^2-\frac{f_0^2}{c^2}$; the group velocity points to the circle's centre."
  8. *Reading the current setting* — short: "Shorter than 2πΛ: the crests go west, the group goes east. Atmospheric
     wave trains do this — a disturbance over the Pacific shows up downstream, to the east, days later." · long: "Longer
     than 2πΛ: crests and energy both go west, almost without dispersion. This is how the ocean interior learns of a
     change of wind at its eastern side — and why the western boundary is where currents pile up." · stationary: "The
     westerly flow carries the pattern east exactly as fast as it drifts west: $c_x=U-\frac\beta{k^2+l^2+f_0^2/c^2}=0$ (13.120). Standing
     waves behind mountain ranges are of this kind." · columns: "A column moved north finds a larger f; to keep
     $\frac D{Dt}\big(\frac{\zeta+f}h\big)=0$ (13.94) it must spin clockwise. Its neighbours' spins push the pattern west."
- **Derivation tab:** **D22** (10 steps) `view: 'plan'`; goal `set` {preset long baroclinic}; step 4 `watch` "the columns'
  colour changes as they move north and south: that is the βv term"; step 7 `live` the three terms of the linearised
  equation for the cursor wave. **D23** (9 steps) `view: 'curve'`; step 1 `live` "ω = `Viz.live('w')` s⁻¹"; step 3 `watch` "the
  chord from the origin always slopes down to the left: phase speed westward"; step 7 `set` {lamx: 2πΛ}, `watch` "the
  tangent is flat: the group stands still"; step 8 `set` {lamx: 150 km}; step 9 `set` {preset stationary}.
- **Code:**
  ```python
  beta = ch13.beta_parameter(np.deg2rad({{lat}})); f0 = ch13.coriolis_parameter(np.deg2rad({{lat}}))
  k = -2*np.pi/{{lamx}}                                          # westward-propagating: k < 0 for omega > 0
  w = ch13.rossby_omega(k, {{l}}, beta, f0, {{c}}, U={{U}})       # (13.118): {{w}} 1/s, period {{T}}
  cx = ch13.rossby_phase_speed(k, {{l}}, beta, f0, {{c}}, U={{U}})   # (13.119): {{cx}} m/s
  cgx, cgy = ch13.rossby_group_velocity(k, {{l}}, beta, f0, {{c}}, U={{U}})   # ours: ({{cgx}}, {{cgy}}) m/s
  ```
- **Walkthrough (7 steps):** 1. "Which way?" — "Press ▶ and follow one crest (▲) and the group as a whole (◆)." `play: true`
  · 2. "Why west" — "A column carried north meets a larger f and must spin clockwise. Click one." `inspect: true`,
  `derive: {id: 'D22', step: 4}` · 3. "The relation" — "Every wave of this family obeys one formula; the minus sign is the
  westward drift." `eq: 'disp'`, `code: {id: 'rw', lines: [3, 3]}` · 4. "Phase and group" — "The chord gives the crests'
  speed, the tangent the group's. Here they have opposite signs." `set` {lamx: 150} · 5. "The turning point" — "At 271 km,
  one Rossby radius times 2π, the tangent is flat." `set` {lamx: 271}, `derive: {id: 'D23', step: 7}` · 6. "Long waves" —
  "Longer waves send their energy west, at up to 3.5 cm/s: nine years across an ocean." `set` {lamx: 2000}, `readouts:
  ['cross']` · 7. "Your turn" — "Predict the mean flow that makes a 6000 km barotropic wave stand still. Then find it."
  `controls: ['rad', 'U', 'lamx']`.
- **Equations:** `qg` "Quasi-geostrophic vorticity equation" ref 'Eq. (13.117)' $\frac\partial{\partial t}\big(\frac{\partial^2\eta}{\partial x^2}+\frac{\partial^2\eta}{\partial y^2}-\frac{f_0^2}{c^2}\eta\big)+\beta\frac{\partial\eta}{\partial x}=0$ · `disp`
  "Dispersion relation" ref 'Eq. (13.118)' $\omega=-\frac{\beta k}{k^2+l^2+f_0^2/c^2}$, live · `cx` "Phase speed" ref 'Eq. (13.119)' $c_x=\frac\omega k=-\frac\beta{k^2+l^2+f_0^2/c^2}$, live
  · `circ` "Circles of constant ω" ref 'p. 674 (no number)' $\big(k+\frac\beta{2\omega}\big)^2+l^2=\big(\frac\beta{2\omega}\big)^2-\frac{f_0^2}{c^2}$, live · `U` "With a mean current" ref 'Eq.
  (13.120)' $c_x=U-\frac\beta{k^2+l^2+f_0^2/c^2}$, live · `pv` "Why" ref 'Eq. (13.94)' $\frac D{Dt}\big(\frac{\zeta+f}h\big)=0$, $f=f_0+\beta y$.
- **Check yourself:** (1) "Is there any Rossby wave whose crests move east (with U = 0)?" — "No: c_x < 0 for every k and
  l." · (2) "At which wavelength does the energy change direction for the baroclinic mode at 35°?" — "2πΛ = 271 km." `set
  {lamx: 271}` · (3) "Why does the tropical ocean adjust faster than the mid-latitude ocean?" — "Because Λ² ∝ 1/f²: the long-wave
  speed βΛ² is nine times larger at 12° than at 35°." `set {lat: 12}` · (4) "What is the wavelength of a stationary wave in
  a 17 m/s westerly?" — "2π√(U/β) = 5983 km." `set {rad: 'barotropic', U: 17}`
- **Selftest parity rows:** `{name: 'omega', js: omega(-2*Math.PI/1.5e5, 0, betaOf(35), fCor(35, 'N'), 3.6096, 0), py:
  'ch13.rossby_omega(-2*np.pi/1.5e5, 0.0, ch13.beta_parameter(np.deg2rad(35.0)), ch13.coriolis_parameter(np.deg2rad(35.0)), 3.6096)',
  rtol: 1e-12}` (3.4275e-7) · `{name: 'cgx short', js: cg(-2*Math.PI/1.5e5, 0, …)[0], py: 'ch13.rossby_group_velocity(-2*np.pi/1.5e5, 0.0,
  ch13.beta_parameter(np.deg2rad(35.0)), ch13.coriolis_parameter(np.deg2rad(35.0)), 3.6096)[0]', rtol: 1e-12}` (+4.347e-3) · `{name:
  'cgx long', js: cg(-2*Math.PI/1e6, 0, …)[0], py: 'ch13.rossby_group_velocity(-2*np.pi/1.0e6, 0.0, …)[0]', rtol: 1e-12}` (−2.807e-2) ·
  `{name: 'cgy', js: cg(-2e-5, 1e-5, …)[1], py: 'ch13.rossby_group_velocity(-2.0e-5, 1.0e-5, …)[1]', rtol: 1e-12}` · `{name: 'omega max',
  js: wmax(betaOf(35), fCor(35, 'N'), 3.6096), py: 'ch13.rossby_max_frequency(ch13.beta_parameter(np.deg2rad(35.0)),
  ch13.coriolis_parameter(np.deg2rad(35.0)), 3.6096)["omega_max"]', rtol: 1e-12}` (4.0457e-7) · `{name: 'stationary', js: lamStat(17,
  betaOf(35)), py: 'ch13.stationary_rossby_wavelength(17.0, ch13.beta_parameter(np.deg2rad(35.0)))', rtol: 1e-12}` (5.98252e6) ·
  `{name: 'evolve', js: evolveRe(…)[5], py: 'np.real(ch13.qg_linear_evolve_1d(np.array([1.0+0j, 0.5+0j]), np.array([-4.0e-5, -4.4e-5]),
  0.0, 1.0e6, ch13.beta_parameter(np.deg2rad(35.0)), ch13.coriolis_parameter(np.deg2rad(35.0)), 3.6096)[1])', rtol: 1e-10}` (the
  builder writes the elided arguments out in full in every row).
- **Fit plan:** 360×640: status · `plan` (50 %) over `curve` (50 %); `kplane` hidden (the sign and size of the group
  velocity are in the curve's title: "group velocity c_gx = +0.43 cm/s → east"); transport with play and scrub only. 844×345: `plan`
  beside `curve`. Desktop: rows [1.2, 1].

### E10 · eady_instability
- **Title:** "Where do storms come from?" · **Summary:** "Slide the wavelength through the cut-off and watch a growing, tilted
  wave fall apart into two neutral edge waves; read the growth time in days." · **CORE:** C16 (also N143, N155, N160, N161,
  N162, N163) · **Reference:** `amplitude_phase_second_order_II_3.html` (response curve + system window + numbered derivation
  with live numbers).
- **meta:** `viz:order 10` · `viz:sections 13.17` · `viz:equations 13.136 13.139 13.141 13.142` · `viz:fluidpy ch13.eady_alpha
  ch13.eady_phase_speed ch13.eady_growth_rate ch13.eady_critical ch13.eady_fastest ch13.eady_max_growth_rate ch13.eady_time_scale
  ch13.eady_mode ch13.eady_fluxes ch13.eady_factors ch13.rossby_radius_internal` · `viz:derivations D25 D26`.
- **Physics (JS ↔ Python):** `alphaOf(k, l, N, f)` ↔ `ch13.eady_alpha`; `factors(aH)` → {tanh_factor, coth_factor, product} ↔
  `ch13.eady_factors` (series for αH/2 < 10⁻³); `cEady(aH, U0)` → [c_r, c_i] ↔ `ch13.eady_phase_speed`; `sigma(k, l, N, f, H, U0)` ↔
  `ch13.eady_growth_rate`; `aCrit()` (bisection on x = coth x, **computed in the page**) ↔ `ch13.eady_critical()`; `fastest()`
  (golden-section search on the growth curve, computed in the page) ↔ `ch13.eady_fastest()`; `sigmaMax(f, N, dUdz)` ↔
  `ch13.eady_max_growth_rate`; `tScale(...)` ↔ `ch13.eady_time_scale`; `mode(z, aH, U0, H)` → {amplitude, phase} ↔
  `ch13.eady_mode`; `fluxes(z, …)` ↔ `ch13.eady_fluxes(...)["w_rho"]`; `LamE(N, H, f)` ↔ `ch13.rossby_radius_internal(N, H, f,
  with_pi=False)`.
- **Views** (rows [1.25, 1]): 1. `sect` "The wave in a longitude–height section" (row 0, flex 1.4): x over two wavelengths,
  z from 0 to H; pressure-perturbation contours (orange high, blue low) from `mode`, growing with the transport and
  renormalised each e-folding (a counter "× e" in the title); the phase line through the highs drawn bold (tilting
  westward with height when unstable, vertical when neutral); at the left edge the basic wind profile U(z) as arrows and
  three sloping density surfaces (blue, slope fU₀/(N²H)). 2. `growth` "Growth rate against wavelength" (row 0, flex 1):
  σ in day⁻¹ against wavelength in km (top axis: αH); cursor dot; the maximum ▲ and the cut-off ◆ marked with their
  computed values; above it a thin strip with the two factors αH/2 − tanh(αH/2) (always positive) and αH/2 − coth(αH/2)
  (negative left of ◆) and the sign of their product. Beyond ◆ the growth curve ends on the axis with the label "two
  neutral waves". 3. `flux` "Heat flux and the wedge" (row 1, `hidePortrait`): left, the profile of −w′ρ′ (upward buoyancy
  flux, positive) and of v′T′ (poleward); right, a parcel path inside the wedge between the horizontal and a density
  surface.
- **Controls:** `lam` (cursor) "Wavelength" (0.3…12 Λ_E, default the fastest wave) · `shear` "Shear $U_0/H$" (0.5…6 m/s per
  km, default 3.0) · `N` "Buoyancy frequency $N$" (0.4…2 × 10⁻² s⁻¹, default 1.1 × 10⁻²; help: "low N: weak stratification,
  fast growth · high N: slow growth") · `lat` "Latitude" (15…75°, default 35) · `sys` modes "atmosphere (H = 9 km) / ocean (H
  = 1 km, N = 5 × 10⁻³ s⁻¹, U₀ = 0.1 m/s)" · transport `t` (in e-folding times).
- **Transport:** `t` 0…3 e-folding times (unstable) or 0…2 periods (neutral), `end: 'hold'`; end card: "after 3 e-foldings
  (4.9 days): amplitude × 20 · the crest moved at U₀/2 = 13.5 m/s".
- **Presets:** "fastest-growing wave" {lam = fastest} · "just inside the cut-off" {αH = 0.97 α_cH} · "just beyond the
  cut-off" {αH = 1.03 α_cH} · "a long wave" {lam 10 Λ_E} · "ocean mesoscale" {sys ocean}.
- **Status:** unstable: "🌀 unstable · e-folding 1.64 days · wavelength 4630 km = 3.9120 Λ_E · moves at U₀/2" (atmosphere
  defaults); beyond the cut-off: "αH = 2.47 above the cut-off 2.3994: two neutral edge waves, speeds 0.44 U₀ and 0.56 U₀"
  (the cut-off number is `aCrit()` computed in the page and formatted to five digits).
- **Readouts:** "Eady radius $\Lambda_E$" · "$\alpha H$" · "Growth rate" · "e-folding time" · "$c_r$, $c_i$" · "Tilt top − bottom".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** (atmosphere / ocean) ·
  **presets** · **status**.
- **Explain:**
  0. *What the views show* — "Section: highs (orange) and lows (blue) of pressure in a west–east slice; the bold line
     joins the highs. Left margin: the basic wind, stronger aloft, and the tilted density surfaces it is in balance with.
     Growth: how fast each wavelength grows; your wave is the dot."
  1. *The Eady radius* — "$\Lambda\equiv HN/f$ (⚠️ no π in this one) = 9000 × 1.1 × 10⁻²/8.37 × 10⁻⁵ = **1183 km**."
  2. *Your wave* — "$\alpha^2\equiv\frac{N^2}{f^2}(k^2+l^2)$ (13.139), so αH = kΛ_E = 2π × 1183/4630 = **1.606**."
  3. *The two factors* — "with x = αH/2 = 0.803: x − tanh x = **+0.137**; x − coth x = **−0.699**; product −0.096 < 0 → the square
     root is imaginary."
  4. *The phase speed* — "$c=\frac{U_0}2\pm\frac{U_0}{\alpha H}\sqrt{\big(\frac{\alpha H}2-\tanh\frac{\alpha H}2\big)\big(\frac{\alpha H}2-\coth\frac{\alpha H}2\big)}$ (13.141) = 13.5 ± 27 × 0.1929 i = **13.5 ± 5.2 i m/s**."
  5. *Growth* — "growth rate σ = k c_i = (1.606/1183 km) × 5.2 = **7.07 × 10⁻⁶ s⁻¹**; e-folding time 1/σ = **1.64 days**. (Ours, computed:
     the largest value over all wavelengths is 0.30982 fU₀/(NH).)"
  6. *The cut-off* — "the product changes sign where x = coth x: αH = **2.3994** (computed here by bisection), i.e.
     wavelengths shorter than 2.6187 Λ_E = 3099 km do not grow."
  7. *The flux view* (or the hint "turn the phone sideways for the heat flux").
  8. *At the current time* — "now t = `Viz.live('t')` e-folding times: amplitude × `Viz.live('amp')`."
  9. *Reading the current setting* — unstable: "The pattern leans westward with height. That lean lets the wave on the
     lower lid and the wave on the upper lid push each other along: poleward-moving air is warm and rising,
     equatorward-moving air cold and sinking. Heat goes poleward and upward, the density surfaces flatten, and the
     wave pays for its growth with the potential energy of their tilt." · beyond the cut-off: "Too short: each boundary
     wave decays away from its lid before it reaches the other (its vertical reach is Λ_E-scaled, 1/α). Unlinked, they
     travel at different speeds and neither grows." · long wave: "Unstable but slow: the growth rate falls like k." ·
     ocean: "The same instability with Λ_E = 60 km: ocean eddies a couple of hundred kilometres across, growing over
     weeks."
- **Derivation tab:** **D25** (13 steps) `view: 'sect'`; goal `set` {preset fastest}; step 3 `watch` "the arrows of the basic
  wind grow linearly with height: that is the thermal wind of the tilted blue surfaces"; step 10 `live` "the two terms
  of w′ at mid-height"; step 13 `watch` "the bracket — the interior potential vorticity — is zero everywhere between the
  lids". **D26** (12 steps) `view: 'growth'`; step 3 `live` "αH = `Viz.live('aH')`"; step 9 `live` the 2 × 2 determinant for
  your c: `Viz.live('det')` ≈ 0; step 11 `set` {αH = 1.03 α_cH}, `watch` "the product under the root turns positive: two real
  speeds"; step 12 `set` {preset fastest}, `live` "c = `Viz.live('c')` m/s".
- **Code:**
  ```python
  f = ch13.coriolis_parameter(np.deg2rad({{lat}}))                 # {{f}} 1/s
  LamE = ch13.rossby_radius_internal({{N}}, {{H}}, f, with_pi=False)   # NH/f = {{LamE}} km (no pi)
  aH = 2*np.pi/{{lam}}*LamE                                        # alpha*H for l = 0: {{aH}}
  c = ch13.eady_phase_speed(aH, {{U0}})                            # (13.141): {{c}} m/s
  sig = ch13.eady_growth_rate(2*np.pi/{{lam}}, 0.0, {{N}}, f, {{H}}, {{U0}})   # k*c_i = {{sig}} 1/s
  print(1/sig/86400, ch13.eady_critical(), ch13.eady_fastest()["sigma_nd"])   # {{days}} d, 2.3994, 0.30982
  ```
- **Walkthrough (7 steps):** 1. "No inflection point" — "This wind profile passes every stability test of Chapter 11. Press
  ▶ anyway." `play: true` · 2. "The stored energy" — "The blue density surfaces tilt. Flattening them would release
  energy; rotation forbids the direct route." `set` {t: 0} · 3. "Two edge waves" — "Each lid carries a wave. Linked through
  the layer, they amplify each other." `eq: 'qg'`, `derive: {id: 'D25', step: 13}` · 4. "The formula" — "Two factors under a
  root. When their product is negative the wave grows." `eq: 'c'`, `derive: {id: 'D26', step: 12}` · 5. "How fast, how big" —
  "Fastest wave: 4630 km long, e-folding in 1.64 days." `set` {preset fastest}, `readouts: ['efold']`, `code: {id: 'eady',
  lines: [4, 6]}` · 6. "Too short" — "Shorten the wave past the ◆: the tilt goes, the growth stops." `set` {αH = 1.03 α_cH} ·
  7. "Your turn" — "Predict what doubling N does to the growth time and to the preferred wavelength. Then try." `controls:
  ['N', 'shear', 'sys']`.
- **Equations:** `qg` "Perturbation equation" ref 'Eq. (13.136)' $\big(\frac\partial{\partial t}+U\frac\partial{\partial x}\big)\big[\nabla_H^2p'+\frac{f^2}{N^2}\frac{\partial^2p'}{\partial z^2}\big]=0$ · `al` "Scaled wavenumber" ref
  'Eq. (13.139)' $\alpha^2\equiv\frac{N^2}{f^2}(k^2+l^2)$, live · `c` "Phase speed" ref 'Eq. (13.141)'
  $c=\frac{U_0}2\pm\frac{U_0}{\alpha H}\sqrt{\big(\frac{\alpha H}2-\tanh\frac{\alpha H}2\big)\big(\frac{\alpha H}2-\coth\frac{\alpha H}2\big)}$, live · `cut` "Unstable band" ref 'Eq. (13.142)' $\frac{HN}f<\frac{\alpha_cH}k$, note "the book
  prints the constant rounded; we compute it from $\frac{\alpha_cH}2=\coth\frac{\alpha_cH}2$", live · `sig` "Growth rate (ours)" ref 'not in the book —
  ours, computed' $\sigma_{max}=0.30982\,\frac fN\frac{dU}{dz}$, live.
- **Check yourself:** (1) "Double the shear. What happens to the e-folding time and to the fastest wavelength?" — "The
  time halves; the wavelength does not change (it is set by Λ_E)." `set {shear: 6}` · (2) "Double N." — "Growth rate
  halves and the preferred wavelength doubles." `set {N: 0.022}` · (3) "Which waves are stable?" — "Those shorter than
  2.6187 Λ_E." · (4) "Is the growth rate largest where c_i is largest?" — "No: c_i is largest as k → 0; the growth rate
  k c_i peaks at αH = 1.6061."
- **Selftest parity rows:** `{name: 'cut-off', js: aCrit(), py: 'ch13.eady_critical()', rtol: 1e-10}` (2.39936) · `{name: 'fastest
  alphaH', js: fastest().alphaH, py: 'ch13.eady_fastest()["alphaH"]', rtol: 1e-6}` (1.60612) · `{name: 'fastest sigma', js:
  fastest().sigma_nd, py: 'ch13.eady_fastest()["sigma_nd"]', rtol: 1e-9}` (0.309817) · `{name: 'c_i', js: cEady(1.0, 27)[1], py:
  'np.imag(ch13.eady_phase_speed(1.0, 27.0))', rtol: 1e-10}` (6.7788) · `{name: 'c neutral', js: cEady(3.0, 27)[0], py:
  'np.real(ch13.eady_phase_speed(3.0, 27.0))', rtol: 1e-10}` (17.864) · `{name: 'factors', js: factors(1.6061).product, py:
  'ch13.eady_factors(1.6061)["product"]', rtol: 1e-10}` · `{name: 'sigma max', js: sigmaMax(fCor(35, 'N'), 1.1e-2, 3e-3), py:
  'ch13.eady_max_growth_rate(ch13.coriolis_parameter(np.deg2rad(35.0)), 1.1e-2, 3.0e-3)', rtol: 1e-9}` (7.0682e-6) · `{name: 'growth
  at a wavelength', js: sigma(2*Math.PI/4.0e6, 0, 1.1e-2, fCor(35, 'N'), 9000, 27), py: 'ch13.eady_growth_rate(2*np.pi/4.0e6, 0.0,
  1.1e-2, ch13.coriolis_parameter(np.deg2rad(35.0)), 9000.0, 27.0)', rtol: 1e-9}` · `{name: 'small-x series', js: factors(1e-5).product, py: 'ch13.eady_factors(1.0e-5)["product"]', rtol: 1e-8}`
  (−8.3333e-12 = −x²/3 with x = 5 × 10⁻⁶; the naive product loses six digits there, which is why both sides use the series).
- **Fit plan:** 360×640: status (one line; the second clause moves to the readouts) · `sect` (55 %) over `growth` (45 %);
  `flux` hidden (the e-folding time and the wavelength are in the growth curve's title). 844×345: `sect` beside `growth`,
  the factor strip dropped (its two numbers are in the Explain tab). Desktop: rows [1.25, 1].

### B1 · inertia_gravity_beams
- **Backup — built only if an explainer above fails review.** **Title:** "What does rotation change about internal waves?" ·
  **Summary:** "Turn the wavevector from vertical to horizontal and watch the frequency run from f to N while the currents
  trace ellipses." · **CORE:** C14 (also N108, N110, N120, N123, N124, N125) · **Reference:** chapter 7's internal-wave
  explainer and `angular_frequency_explorer_1.html`.
- **meta:** `viz:order 11` · `viz:sections 13.14` · `viz:equations 13.99 13.110 13.112` · `viz:fluidpy ch13.inertia_gravity_omega
  ch13.inertia_gravity_m2 ch13.inertia_gravity_group_velocity ch13.inertia_gravity_hodograph ch13.inertia_gravity_band` ·
  `viz:derivations D19 D20 D21`.
- **Physics (JS ↔ Python):** `igOmega(k, m, N, f)` ↔ `ch13.inertia_gravity_omega`; `igM2(k, l, w, N, f)` ↔ `ch13.inertia_gravity_m2`;
  `igCg(k, m, N, f)` ↔ `ch13.inertia_gravity_group_velocity`; `igHodo(t, w, f, sign)` ↔ `ch13.inertia_gravity_hodograph`; `igBand(w,
  N, f)` ↔ `ch13.inertia_gravity_band`.
- **Views:** 1. `sect` vertical section: phase lines moving (transport), phase-velocity arrow (purple) and group-velocity
  arrow (teal) at right angles, a parcel sliding along the crests; 2. `freq` ω against the wavevector angle from 0° to
  90°, between f and N on a log axis, with the cursor and the f = 0 curve N cos θ as a ghost; 3. `hodo` (`hidePortrait`) the
  horizontal hodograph, an ellipse of axis ratio f/ω with its sense of rotation.
- **Controls:** `th` "Angle of the wavevector" (0.5…89.9°) · `fN` "$f/N$" (0…0.2; 0 = Chapter 7) · `lam` "Wavelength" · `up` chips
  "phase up / down" · `hemi` chips N / S · transport `t`. **Presets:** near-inertial (89°) · hydrostatic mid-range (80°) · near N
  (5°) · f = 0 (Chapter 7). **Status:** "in band: ω = 2.98 f = 0.092 N" / "near-inertial: currents almost circular" / "f = 0:
  ω = N cos θ, currents along a line".
- **Explain:** 0 views and colours · 1 "$\omega^2-f^2=\frac{k^2}{m^2}(N^2-\omega^2)$ (13.112) solved for ω: with tan θ = m/k, ω² = f² sin² θ + N² cos² θ
  (unnumbered, p. 669) = … = **…**" · 2 "group velocity (D21): size and direction; phase · group = 0" · 3 "the ellipse: $u=\mp\cos\omega t$,
  $v=\pm\frac f\omega\sin\omega t$ (13.110), ratio f/ω = …" · 4 "local vertical wavenumber $m^2(z)\equiv\frac{(k^2+l^2)[N^2(z)-\omega^2]}{\omega^2-f^2}$ (13.99)" · 5 *Reading the current setting*: near N — "almost vertical motion, buoyancy alone; rotation invisible"; mid-range — "both matter; the
  hydrostatic approximation holds"; near f — "almost horizontal circles: an inertial oscillation with a slow vertical
  phase"; f = 0 — "Chapter 7's beam".
- **Derivation tab:** D19 (11 steps), D20 (7), D21 (7), `view: 'sect'`. **Code:** five lines calling the five functions.
  **Walkthrough (5 steps):** the question · the angle sets the frequency · between f and N · phase against group · your
  turn. **Check yourself:** "At what angle is ω = √(fN)?" · "Which way does the energy go when the crests move up?" · "What
  happens to the ellipse as ω → f?" **Selftest parity rows:** `{name: 'omega', js: igOmega(6.283e-4, 3.1416e-2, 2.7e-3, fCor(35,
  'N')), py: 'ch13.inertia_gravity_omega(6.283e-4, 3.1416e-2, 2.7e-3, ch13.coriolis_parameter(np.deg2rad(35.0)))', rtol: 1e-12}` ·
  `{name: 'cgz', js: igCg(…)[1], py: 'ch13.inertia_gravity_group_velocity(6.283e-4, 3.1416e-2, 2.7e-3,
  ch13.coriolis_parameter(np.deg2rad(35.0)))[1]', rtol: 1e-12}` · `{name: 'f = 0', js: igOmega(1, 2, 1, 0), py:
  'ch13.inertia_gravity_omega(1.0, 2.0, 1.0, 0.0)', rtol: 1e-12}` (0.447214). **Fit plan:** 360×640: `sect` over `freq`, `hodo`
  hidden (axis ratio in the title); 844×345 `sect` beside `freq`.


---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU; target 170–210 s)

Per-cell costs are estimates from the designer's scratch runs where one exists (marked ✓: the 1-D march of 800 cells × 6
inertial periods took 0.5 s including Python start-up; each sympy check under 3 s; closed forms negligible) and from the
analysis (§4 cost column) otherwise. The builder reports the measured numbers.

| Section | What costs time (per cell) | Full | FAST |
|---|---|---|---|
| front matter, §13.1 | imports (sympy, plotly, scipy); `illustrative_inputs`; two `pressure_centre` panels | 5 s | 5 s |
| §13.2 | USSA profile, idealised ocean profile; **F2** (29 steps × 2 traces × 200 pts) ≈ 2 s; `lapse_rate_table` | 5 s | 4 s |
| §13.3 | `boussinesq_rotating_sympy` ×2, `perturbation_form_sympy`, `eddy_friction_force_sympy` (≤ 2 s each, cached to `outputs/ch13/cache/`, P252) | 7 s (first) / 2 s | 2 s |
| §13.4 C01 | plotly 3-D tangent plane (4 dropdown cases) 2 s; `thin_layer_sympy` 2 s; term bars; **F1** (15 steps × 3 traces × 200 pts) 2 s | 8 s | 6 s |
| §13.5 C02, C03 | 2 × 2 pressure-centre figure; parcel paths (closed form); `jet_section` on 121 × 61; **F3** (21 steps × 2 hemispheres) 2 s; `taylor_proudman_sympy` ×2 (3 s); E1, E2 iframes | 10 s | 8 s |
| §13.6 C04, C05 | `fig_ekman_surface`; plotly 3-D spiral; **F4** (25 steps × 3 traces × 200 pts) 2 s; `ekman_solve` ×3 at n = 400 (< 0.2 s); `ekman_pumping` on 81 × 81; D06 sympy check (2 s ✓); live cell (not executed on the page); E3 iframe | 11 s | 9 s |
| §13.7 C06 | `fig_ekman_bottom`; force-triangle figure; `ekman_finite_depth`; **F5** (15 steps × 3 dropdown cases) 2 s; E4 iframe | 6 s | 5 s |
| §13.8 C07 | `linear_1d_run` 400 cells × 144 steps ×2 (< 0.2 s ✓); from-scratch march (same size, Python loop over steps only) 0.3 s; Hovmöller figure; the `slow` cell `SW.run` 64² × 20 steps ≈ 2 s (skipped in FAST) | 5 s | 2 s |
| §13.9 C08 | `vertical_modes` n = 401 for two profiles and two lids (< 1 s); `modes_uniform_N` (`brentq`); `fig_mode_roots`, `fig_vertical_modes`; **F6** (15 steps × eigen-solve at n = 201 ≈ 0.05 s each) 2 s; D11 sympy check (2 s ✓); E5 iframe | 9 s | 7 s |
| §13.10 C09 | closed-form cubic on 400 wavenumbers; term-bar figure; **F7** (15 latitudes × 2 modes × 4 traces × 300 pts) 3 s; D12 sympy check (1.5 s ✓); `v_equation_sympy` (≤ 5 s, cached) | 10 s (first) / 6 s | 5 s |
| §13.11 C10 | `fig_poincare_kelvin_dispersion`; **A1** (80 frames, dpi 80, closed form; FAST 40) ≈ 6 s / 3 s; E6 iframe | 9 s | 5 s |
| §13.12 C11, C12 | `fig_kelvin_sections`; plotly surface; **A3** from `reference/ch13/kelvin_basin.npz` (24 frames) ≈ 3 s / 3 s — never recomputed in the notebook; **A2** march 800 cells × 6 inertial periods ×2 (0.5 s ✓) + 4 stop-frames ×2 rows (frames player) 3 s; D16 sympy check (1.7 s ✓); E7, E8 iframes | 14 s | 9 s |
| §13.13 C13 | `flow_over_step` (closed form); PV-particle figure from `reference/ch13/pv_particles.npz` (< 1 s); D17 sympy check (1.2 s ✓); `pv_conservation_sympy` (≤ 5 s, cached) | 9 s (first) / 4 s | 4 s |
| §13.14 C14 | closed forms; `vertical_structure_solve` (DOP853, rtol 1e-10) ×3 ≈ 0.6 s; `wkb_error` ×3; plotly helix; lee-wave streamlines 201 × 101; D19 sympy check (2.2 s ✓); `w_equation_rotating_sympy` (≤ 5 s, cached) | 12 s (first) / 7 s | 6 s |
| §13.15–13.16 C15 | `fig_rossby_dispersion`; **F8** (20 steps × 2 traces × 300 pts) 2 s; **A4** (90 frames, 48 modes × 600 pts, exact; FAST 40) ≈ 7 s / 3 s; `qg_linear_evolve` 128 × 64 × 3 times (< 0.1 s); D22 sympy check (1.9 s ✓); `rayleigh_kuo_eigs` at ≤ 12 (k, β) points ≈ 10 s (cached to `outputs/ch13/cache/`; FAST 4 points); E9 iframe | 26 s (first) / 16 s | 10 s |
| §13.17 C16 | closed forms on 400 wavenumbers; `eady_numeric_eigs` at 8 wavenumbers (≈ 1 s); **F9** (19 steps × 2 systems) 2 s; **A5** (80 frames, closed form; FAST 40) ≈ 6 s / 3 s; D25, D26 sympy checks (1.9 s + 2.0 s ✓); `eady_qg_sympy` (≤ 5 s, cached); live cell; E10 iframe | 20 s (first) / 15 s | 10 s |
| §13.18 C17 | triad map 200 × 200; `two_d_cascade_spectrum`; **A6** from `reference/ch13/turbulence_f.npz` and `turbulence_beta.npz` (2 × 10 saved frames at 64², frames player) ≈ 4 s; live 64² re-run only if a cache is missing and `not FAST` (≈ 4 s per run) | 8 s | 6 s |
| explainers | 10 `show_viz` cells (iframes; no computation) | 3 s | 3 s |
| **Total** | | **≈ 180 s on a first run (cold sympy and eigenvalue caches); ≈ 145 s afterwards** | **≈ 105 s** |

**FAST plan.** `FAST = setup_notebook()`; animations ≤ 90 frames at dpi 80 (FAST 40): A1, A4, A5 `player="video"`; A2 and A6
`player="frames"` (the reader should stop at 0, 1, 3, 5 inertial periods, and at the log-spaced times); A3 `"video"`. Slider
figures ≤ 29 steps × ≤ 4 traces × ≤ 300 points (≈ 150 kB each; nine of them). The twelve sympy engines are cached with a
parameter hash in `outputs/ch13/cache/` (git-ignored); the nine ★★★ `check_src` cells always run live (17 s together ✓).

**Cached model output under `reference/ch13/` (our own output only; committed; 405 kB together as measured; regenerated by
`scripts/ch13_make_caches.py --which all`, which prints sizes and a checksum):**

| File | What | Made by | Size (target) | Cost to regenerate |
|---|---|---|---|---|
| `kelvin_basin.npz` | surface height η(t, y, x) of a Kelvin wave going round a closed 64 × 64 basin, equivalent depth chosen so that Λ = 1/5 of the basin; 24 frames | `SW.run` (C-grid, linear, closed) | 185 kB (measured) | ≈ 8 s |
| `turbulence_f.npz` | vorticity ζ(t, y, x) of decaying two-dimensional turbulence, **64²**, β = 0; 10 saved frames + energy, enstrophy, centroids K_E, K_Z at every saved step (the K⁻³ range is not resolved at this size: qualitative use only) | `SW.barotropic_run` | 78 kB (measured) | ≈ 4 s |
| `turbulence_beta.npz` | the same with β ≠ 0 (zonal share of the energy at the last frame 0.29 against 0.11 without β) | `SW.barotropic_run` | 78 kB (measured) | ≈ 4 s |
| `explainer_constants.json` | the constants the explainers pin by parity rows (Eady cut-off and fastest wave, the adjustment ratio, Ω, R, our inputs) | the implementer's cache script | 1 kB | < 1 s |
| `pv_particles.npz` | 24 particle tracks at 13 times with (ζ + f)/h along each (key `q`) and the η field, from a 48² nonlinear run | `SW.run(linear=False)`, `SW.advect_particles` | 61 kB (measured) | ≈ 10 s |

Every cell that reads a cache goes through `ch13.load_reference_run(name)` and has a stated fallback when it returns
`None` (A3: the straight-coast closed form; A6: a live 64² run unless FAST; the C13 particle figure: skipped with one
sentence). No book data, no observation and no third-party model output is stored. Page budget: six animations < 2 MB each
(A3 and A6 are the largest), nine plotly figures ≈ 150 kB each, three plotly 3-D figures ≈ 300 kB each: ≈ 11 MB < 15 MB.

---


## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is explained.
"primer (in Cxx)" = a 📎 primer placed in that block before first use (the Concept text is the primer term, used verbatim in
`nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line reminder naming the earlier primer; a CORE/RECAP
id alone = taught there; "Cxx (Nnn)" = the NOTE placed in that block; "Cxx (D0n)" = the derivation where it is built; "Cxx (gloss
…)" = one sentence where it is used. Items of §13.1–§13.3 are first used in the story that opens C01 and are entered as C01.
New primers P307–P335 (29) in first-use order: geographic vocabulary (§13.1), sidereal day (C01), pressure map (C02),
hodograph, complex velocity (C04), depth-integrated quantities, first integral (C05), C-grid step (C07), Sturm–Liouville,
Robin condition, graphical roots, `eigh_tridiagonal`, projection on modes (C08), three real roots of a cubic, reading a
dispersion diagram (C09), polarisation relations, sense of rotation (C10), trapped solutions (C11), adjustment problem (C12),
materially conserved (C13), group velocity as a gradient, WKB (C14), ordering in a small parameter, two-dimensional FFT grids
(C15), homogeneous 2 × 2 systems, hyperbolic identities and coth, available potential energy (C16), enstrophy, the Jacobian
and a pseudo-spectral step (C17).

| Concept | First used in | Explained by |
|---|---|---|
| the conventions block: axes x east, y north, z up; where z = 0 is in each section; the hemisphere rule; the overloaded letters f, Ω / ω, β, N, H / h / η, c / c_n, k / l / m / K, ζ, θ, δ, Λ / λ, ψ, E / Ro / R, α, U / V / τ | C01 | front matter (A.0 row 5, before C01), repeated per block |
| slips kept apart from traps; `ch13.book_slips()`, `ch13.conventions_table()`, `ch13.illustrative_inputs()` | C01 | front matter (A.0 row 5) |
| geophysical fluid dynamics: rotation + stratification; flow along isobars | C01 | C01 (N01) |
| geographic vocabulary: zonal / meridional, easterly wind vs eastward current, poleward / equatorward, cyclonic / anticyclonic | C01 | primer (in C01) |
| northern-hemisphere convention of the book; sign(f); `GFD.hemisphere`, `ch13.wind_from_to` | C01 | C01 (N02) |
| potential density and the stability test dρ_θ/dz < 0 | C01 | R01 |
| adiabatic density gradient −gρ/c²; compression share of the in-situ gradient | C01 | R02 and C01 (N03) |
| lapse rate in two conventions: Kundu Γ ≡ dT/dz (Γ_a ≈ −9.8 K/km) and meteorological Γ_met ≡ −dT/dz (Γ_d ≈ +9.8 K/km) | C01 | R02 |
| negate the number, flip the inequality | C01 | knowledge/primers.md: inequalities under a sign change (ch01 P48) — reminder |
| `ch13.lapse_rate_table`, `STRAT.lapse_rate_stability`, `STRAT.adiabatic_lapse_rate`; the library string "Γ < Γa" read as Γ_met < Γ_d | C01 | R02 (the trap box under it) |
| troposphere, tropopause, stratosphere; USSA-1976 profile | C01 | C01 (N04) |
| "close to neutral" = close to the moist adiabat; dry stability of the standard troposphere | C01 | C01 (N05) |
| `ch12.gradient_richardson_thermal(…, Gamma_a=…)` and the gradient Richardson number | C01 | C01 (N05), recalling the Richardson-number block of ch12 |
| mixed layer, thermocline, abyss; idealised ocean profile | C01 | C01 (N06) |
| buoyancy frequency N² from potential density | C01 | R03 |
| period 2π/N of the thermocline | C01 | C01 (N07) |
| rotating Boussinesq equations; Coriolis acceleration 2Ω × u | C01 | R04 |
| scale height c²/g and the validity of Boussinesq | C01 | R05 |
| linear equation of state, expansion coefficient α, haline coefficient β_S | C01 | R06 |
| hydrostatic rest state; rest state + perturbation; primes dropped | C01 | R07, R08 |
| perturbation form of the pressure and gravity terms | C01 | R09 |
| friction force per unit mass; eddy-viscosity hypothesis | C01 | R10 |
| anisotropic eddy stress with ν_H ≫ ν_v; frame indifference | C01 | C01 (N08, N09) |
| friction force ν_H∇_H² + ν_v∂_z² | C01 | C01 (N10, D01) |
| index notation and the summation convention | C01 | C01 (D01 tools; ch02 reminder) |
| second partial derivatives; mixed partials commute | C01 | knowledge/primers.md: Schwarz's theorem (ch04 P121) — reminder |
| partial derivative ∂/∂x | C01 | knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| illustrative eddy coefficients | C01 | C01 (N11) |
| sidereal day vs solar day | C01 | primer (in C01) |
| earth's rotation rate Ω (sidereal), `OMEGA_EARTH` | C01 | R11 |
| latitude θ and the local vertical | C01 | knowledge/primers.md: latitude, Earth's rotation rate and the local vertical (ch04 P126) — reminder |
| tangent plane; Ω = (0, Ω cos θ, Ω sin θ) | C01 | C01 (N13) |
| plotly 3-D figures with a dropdown | C01 | knowledge/primers.md: plotly 3-D arrows, lines and meshes (ch02 P64) — reminder |
| aspect ratio; W/U ~ H/L | C01 | C01 (N12) |
| scaling with two length scales | C01 | knowledge/primers.md: anisotropic scaling with two length scales (ch08 P188) — reminder |
| cross product as a determinant | C01 | knowledge/primers.md: matrices, determinants and minors (ch01 P53) — reminder |
| traditional approximation; (−fv, fu, −2Ωu cos θ) | C01 | C01 (N14, D02) |
| Coriolis parameter f = 2Ω sin θ | C01 | R12 |
| inertial period 2π/f | C01 | C01 (N15) |
| thin-shell equations | C01 | C01 |
| `GFD.thin_layer_terms`, `ch13.term_table_thin_layer`, `GFD.coriolis_acceleration_local`, `np.cross` | C01 | C01 (code row 13, from-scratch row 14) |
| logarithmic bar chart | C01 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| f-plane | C01 | C01 (N16) |
| β-plane, β = 2Ω cos θ₀/R | C01 | C01 (N17) |
| first-order Taylor expansion | C01 | knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder |
| `slider_figure` (plotly) | C01 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| `assert np.allclose` / `np.isclose` | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| sympy engines returning a residual and "ok" | C01 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| reading a pressure map: isobars, the pressure-gradient force points from high to low, tight spacing = strong force | C02 | primer (in C02) |
| Rossby number U/(fL) against ch04's U/(2ΩL) | C02 | R13 |
| Ekman number ν/(fL²) | C02 | C02 (N22) |
| order-of-magnitude scaling | C02 | knowledge/primers.md: order-of-magnitude scaling (ch04 P130) — reminder |
| geostrophic balance | C02 | C02 |
| u·∇p = 0; pressure as a stream function ψ = p/(fρ₀); sign of ψ | C02 | C02 (N18, D03) |
| highs and lows, cyclonic and anticyclonic sense in both hemispheres | C02 | R14 |
| `np.meshgrid`, the [j, i] grid layout, array slicing for centred differences | C02 | knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76) — reminder |
| contour, quiver | C02 | knowledge/primers.md: plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminder |
| inertial loops of a released parcel; `GFD.parcel_adjust`; complex numbers in numpy | C02 | C02 (row 13) and knowledge/primers.md: the complex plane in numpy (ch06 P153) — reminder |
| `show_viz` explainers | C02 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| hydrostatic balance of the perturbation | C03 | R15 |
| cross-differentiation to eliminate a variable | C03 | knowledge/primers.md: operator elimination for linear PDEs (ch07 P177) — reminder |
| thermal wind | C03 | C03 |
| thermal wind in temperature form (ours); α as a required keyword | C03 | C03 (N20, D04) |
| baroclinic vs barotropic | C03 | C03 (gloss in the plain-words row: surfaces of constant p and ρ do or do not coincide) |
| `scipy.integrate.cumulative_trapezoid` | C03 | knowledge/primers.md: `scipy.integrate.cumulative_trapezoid` (ch08 P190) — reminder |
| rotating-tank geostrophy; no stretching; no shear | C03 | C03 (N21, N23, N24) |
| Taylor–Proudman theorem and its hypotheses; Taylor column | C03 | C03 (N25, N26, D05) |
| hypotheses of a theorem (necessary conditions) | C03 | knowledge/primers.md: necessary vs sufficient conditions, "for every k", and proof by contradiction (ch11 P255) — reminder |
| three ways to balance friction; δ ~ √(νt), √(νx/U), √(2ν_v/f) | C04 | R16 |
| Ekman balance and its boundary conditions; wind stress τ | C04 | C04 (N27, N28) |
| hodograph: the tip of the velocity vector traced as depth (or time) changes | C04 | primer (in C04) |
| complex velocity V = u + iv: multiplying by i turns a vector 90° to the left; e^{(1+i)s} decays and turns at the same rate | C04 | primer (in C04) |
| trying e^{rz} in a linear ODE with constant coefficients | C04 | knowledge/primers.md: linear second-order ODE (ch01 P44) — reminder |
| square root of i; conjugates; polar form; Euler's formula | C04 | knowledge/primers.md: complex square roots and the quadratic formula (ch06 P159) — reminder |
| Ekman thickness δ and Ekman depth πδ | C04 | C04 (N29, N30, D06) |
| surface Ekman spiral | C04 | C04 |
| Stokes layer as the laminar cousin | C04 | C04 (trap T6 box; ch08 reminder) |
| added geostrophic interior flow | C04 | C04 (N33) |
| vortex-tilting balance | C04 | C04 (N35) |
| depth-integrated quantities: transport per unit width [m²/s] and its divergence as a vertical velocity | C05 | primer (in C05) |
| first integral: integrating an ODE once across a layer to get a transport without solving it | C05 | primer (in C05) |
| integral of an exponential over a half-line | C05 | knowledge/primers.md: improper integral as a limit (ch05 P144) — reminder |
| `np.trapezoid` | C05 | knowledge/primers.md: trapezoid rule (ch01 P37) — reminder |
| Ekman transport | C05 | C05 |
| transport independent of the eddy-viscosity model | C05 | C05 (N32, D07) |
| Ekman layer with K_v(z); finite-volume tridiagonal solve; `GFD.ekman_solve` | C05 | C05 (N37, the "our choice" box) and knowledge/primers.md: Crank–Nicolson with `scipy.linalg.solve_banded` (ch08 P193) — reminder |
| Ekman pumping; curl of the wind stress; Sverdrup balance (one line) | C05 | C05 (N36, D09) |
| divergence and curl of a horizontal vector | C05 | C05 (D09 tools; ch02 reminder) |
| coastal upwelling | C05 | C05 (row 11) and C11 (N96) |
| wind-driven gyres, western intensification (named) | C05 | C05 (N50) |
| live widgets (`live`) | C05 | C05 (the live note written by `nb.live`) |
| geostrophic interior above a rigid surface; three-way balance; conditions | C06 | C06 (N38–N41) |
| particular solution of a forced linear ODE | C06 | C06 (D08 step 2) and knowledge/primers.md: linear second-order ODE (ch01 P44) — reminder |
| bottom Ekman layer; cross-isobar angle; overshoot at 3πδ/4 | C06 | C06 |
| cross-isobar transport ½Uδ | C06 | C06 (N47, D08) |
| laminar thickness vs observed; eddy viscosity from depth | C06 | C06 (N48) |
| force triangle; inflow into lows; bottom pumping; spin-down | C06 | C06 (N49) |
| finite-depth Ekman solution (ours) | C06 | C06 (N34) |
| shallow layer geometry; which z = 0 | C07 | C07 (N54) |
| hydrostatic pressure under a displaced surface; surface-slope pressure gradient | C07 | R17, R18 |
| kinematic surface condition w(η) = Dη/Dt | C07 | C07 (D10 step 4; ch07 reminder) |
| product rule | C07 | knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| shallow-water continuity (flux form) | C07 | C07 (N51, N52, D10) |
| linearisation (drop products of small quantities) | C07 | C07 (D10 steps 7–9) |
| linear shallow-water equations | C07 | C07 |
| a C-grid shallow-water step and its CFL limit with c = √(gH) (our choice of scheme) | C07 | primer (in C07) |
| staggered arrays (η at centres, u at faces) | C07 | knowledge/primers.md: half-index notation and staggered array shapes (ch10 P242) — reminder |
| x–t (Hovmöller) diagram | C07 | C07 (figure row 13, the "How to read it" note) |
| equivalent depth | C07 | C07 (N53) |
| five linear hydrostatic equations of a stratified layer | C08 | R19, R20, C08 (N55) |
| Sturm–Liouville problem: a ladder of eigenvalues, eigenfunctions with n zero crossings, orthogonality | C08 | primer (in C08) |
| separation of variables | C08 | knowledge/primers.md: separation of variables for a PDE (ch07 P167) — reminder |
| eigenvalue problem for a differential operator | C08 | knowledge/primers.md: eigenvalue problem for a differential equation (ch11 P257) — reminder |
| integration by parts | C08 | knowledge/primers.md: integration by parts (ch09 P218a) — reminder |
| vertical modes ψ_n, their w- and ρ-structures, the separation constant, modal equations | C08 | C08 (N56–N64, D11) |
| units of the modal amplitudes | C08 | C08 (trap T12 box) |
| bottom and free-surface conditions on the modes | C08 | C08 (N65, N66, N67), R21 |
| Robin (mixed) boundary condition, between Dirichlet and Neumann | C08 | primer (in C08) |
| Dirichlet and Neumann conditions | C08 | knowledge/primers.md: boundary conditions (ch01 P20) — reminder |
| uniform-N eigenproblem and its roots | C08 | C08 (N68–N71) |
| graphical roots of a transcendental equation such as tan x = εx | C08 | primer (in C08) |
| `scipy.optimize.brentq` | C08 | knowledge/primers.md: scipy.optimize.brentq (ch03 P108) — reminder |
| barotropic mode | C08 | R22 |
| baroclinic modes c_n = NH/(nπ); first baroclinic mode | C08 | C08 (N73, N74) |
| scipy.linalg.eigh_tridiagonal for a discretised eigenproblem | C08 | primer (in C08) |
| eigenvalues and eigenvectors; generalised symmetric eigenproblem | C08 | knowledge/primers.md: eigenvalues and eigenvectors (ch02 P80) — reminder |
| rigid-lid approximation | C08 | C08 (N75) |
| projecting a profile on modes: the inner product of two functions | C08 | primer (in C08) |
| `VM.vertical_modes`, the `Modes` named tuple, `VM.modes_uniform_N`, `VM.project`, `VM.reconstruct`, `VM.orthogonality_matrix`, `VM.wkb_mode_speed` | C08 | C08 (code rows 17–21) |
| limits of the modal decomposition | C08 | C08 (N77) |
| plane wave e^{i(kx+ly−ωt)}; wavenumbers k, l; K | C09 | knowledge/primers.md: complex amplitudes (ch07 P176) — reminder |
| elimination to one equation for v; linear vorticity equation on the β-plane | C09 | C09 (N78–N81, D12) |
| three real roots of a cubic: the discriminant and the trigonometric form | C09 | primer (in C09) |
| `np.roots` | C09 | knowledge/primers.md: np.roots and np.lib.scimath.sqrt (ch11 P256) — reminder |
| complete dispersion relation (the cubic) | C09 | C09 |
| all roots real under β-plane scaling; superinertial and subinertial | C09 | C09 (N82, D13) |
| dominant balance; three frequency regimes | C09 | C09 (N83) and knowledge/primers.md: dominant balance (ch08 P198) — reminder |
| reading a dispersion diagram with several branches: signed ω for signed k, a logarithmic frequency axis | C09 | primer (in C09) |
| phase speed ω/k | C09 | knowledge/primers.md: phase of a wave (ch07 P165) — reminder |
| plane-wave amplitude equations | C10 | C10 (N84) |
| polarisation relations: solving a 2 × 2 complex system for the velocity amplitudes | C10 | primer (in C10) |
| `np.linalg.solve` | C10 | knowledge/primers.md: np.linalg.solve (ch01 P57) — reminder |
| Poincaré (Sverdrup, rotational gravity) waves | C10 | C10 |
| chain rule (implicit differentiation of ω²) | C10 | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| group velocity (one-dimensional dω/dk) | C10 | C10 (D14 step 7; ch07 reminder) |
| real parts of complex wave fields; velocity ellipse | C10 | C10 (N88, N89) |
| sense of rotation of (a cos ωt, b sin ωt) from the sign of the swept area | C10 | primer (in C10) |
| a linear map of a circle is an ellipse | C10 | knowledge/primers.md: linear map of a circle is an ellipse (ch03 P104) — reminder |
| `animate` and `show_animation` | C10 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| inertial oscillation; inertial circle of radius q/f | C10 | C10 (N90) |
| Kelvin-wave idea: v = 0, cross-shore geostrophy | C11 | C11 (N91, N93, N94) |
| trapped solutions: keeping the exponential that decays away from a boundary | C11 | primer (in C11) |
| first-order linear ODE | C11 | knowledge/primers.md: first-order linear ODE and the integrating factor (ch09 P210) — reminder |
| Kelvin wave; non-dispersive at √(gH) | C11 | C11 and R23 |
| `GFD.kelvin_decay_side`; `try / except ValueError` | C11 | C11 (code row 10, gloss in "What does the code above do?") |
| cached model runs; `ch13.load_reference_run` | C11 | C11 (the "our choice" box of A3) and knowledge/primers.md: caching expensive runs (np.savez and a parameter key) (ch10 P252) — reminder |
| Rossby radius of deformation Λ = c/f | C12 | C12 |
| reduced gravity; internal Kelvin waves; internal Rossby radius; three radii | C12 | R24 and knowledge/primers.md: reduced gravity and buoyancy (ch04 P131) — reminder |
| an adjustment problem: what a steady end state can remember (a conserved quantity pins it) | C12 | primer (in C12) |
| linear potential vorticity ζ − fη/H | C12 | C12 (D16 steps 2–3) |
| geostrophic adjustment; the adjusted front and jet; the energy split 1/3 | C12 | C12 (N19, D16) |
| potential and kinetic energy per unit length of a layer | C12 | C12 (D16 step 11) |
| frames player for stop-and-look animations | C12 | C12 (A2 note) and knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| nonlinear shallow-water equations over an uneven bottom; total depth h | C13 | C13 (N97, N98, N99, N104) |
| relative and absolute vorticity | C13 | R25 |
| materially conserved: Dq/Dt = 0 labels a parcel; it does not mean q is steady at a point | C13 | primer (in C13) |
| material derivative D/Dt | C13 | C13 (reminder row 6; ch03) |
| quotient rule | C13 | knowledge/primers.md: quotient rule, and differentiating with respect to K² as the variable (ch11 P265) — reminder |
| vorticity equation with divergence; absolute vorticity and stretching | C13 | C13 (N100, N101, D17) |
| potential vorticity (ζ + f)/h and its conservation | C13 | C13 and R26, R27 |
| streamline displacement Y(x); flow over a step; upstream influence | C13 | C13 (N102, N103, D18) |
| internal waves without rotation | C14 | R28 |
| linear rotating Boussinesq set (non-hydrostatic) | C14 | R29 |
| the w-equation with rotation | C14 | C14 (N105, D19) |
| vertical-structure equation; local vertical wavenumber m(z); oscillator form | C14 | C14 (N106–N109, N119, D20) |
| the band f < ω < N; evanescent | C14 | C14 (N110) |
| inertia–gravity dispersion relation; wavevector angle θ (not latitude) | C14 | C14 |
| group velocity as the gradient of ω in wavenumber space, read off a contour plot | C14 | primer (in C14) |
| group velocity of inertia–gravity waves; phase ⟂ group | C14 | C14 (N124, D21) |
| three frequency regimes (non-rotating, hydrostatic, mid-range) | C14 | C14 (N123) |
| slowly varying medium (WKB): amplitude and phase ansatz, valid when the medium changes little in one wavelength | C14 | primer (in C14) |
| WKB solution; eikonal; amplitude ∝ m^{−1/2} | C14 | C14 (N111–N114) |
| `scipy.integrate.solve_ivp` | C14 | knowledge/primers.md: scipy.integrate.solve_ivp (ch01 P31) — reminder |
| velocity field, hodograph and tilted orbit plane of an inertia–gravity wave | C14 | C14 (N115–N118, N120, N121, N122) |
| Doppler shift; lee waves; cut-off k < N/U | C14 | C14 (N125, N126; ch07 reminder) |
| quasi-geostrophic motion; planetary waves | C15 | C15 (N127) |
| ordering in a small parameter: lowest order gives the balance, next order gives the evolution | C15 | primer (in C15) |
| quasi-geostrophic vorticity equation | C15 | C15 (N128–N131, D22) |
| Rossby-wave dispersion relation | C15 | C15 |
| westward phase speed; circles of constant ω; group velocity; maximum frequency | C15 | C15 (N132, N133, N134, D23) |
| completing the square into a circle | C15 | knowledge/primers.md: completing the square into a circle (x − a)² + y² ≤ R² (ch11 P272) — reminder |
| complex-step derivative | C15 | C15 (gloss in the from-scratch row 11: perturb by ih, take the imaginary part, divide by h) |
| long-wave speed −βΛ²; basin-crossing time | C15 | C15 (N135) |
| mean-current (Doppler-shifted) Rossby waves; stationary wavelength | C15 | C15 (N136) |
| cubic reduces to the Rossby relation | C15 | C15 (N137) |
| two-dimensional FFT wavenumber grids (np.fft.fftfreq, np.fft.fft2) | C15 | primer (in C15) |
| Fourier modes | C15 | knowledge/primers.md: Fourier modes and the FFT Poisson solver (ch05 P142) — reminder |
| wave packet as a sum of modes; `SW.qg_linear_evolve_1d`, `GFD.rossby_packet_spectrum` | C15 | C15 (A4 row and its "our choice" box) |
| equatorial waves (named) | C15 | C15 (N138) |
| barotropic vorticity equation; perturbation equation about U(y) | C15 | C15 (N139, N140) |
| Rayleigh equation with β; normal modes with complex c | C15 | R30 and knowledge/primers.md: linear stability of a steady configuration (perturb, linearise, eigenvalues) (ch09 P214) — reminder |
| real and imaginary parts of a complex equation | C15 | knowledge/primers.md: real and imaginary parts of a complex identity (ch11 P260) — reminder |
| Rayleigh–Kuo criterion (necessary, not sufficient) | C15 | C15 (N141, N142, D24) |
| sloping density surfaces; Eady basic state | C16 | C16 (N143) |
| available potential energy and the wedge of sloping convection | C16 | primer (in C16) |
| f-plane hydrostatic set; decomposition; basic-state balance; thermal wind U = U₀z/H | C16 | C16 (N144–N147) |
| perturbation vorticity, density and hydrostatic equations; w′ from the pressure | C16 | C16 (N148–N154, D25) |
| Eady perturbation equation; interior potential vorticity zero | C16 | C16 (N155) |
| scaled wavenumber α; vertical structure about mid-depth | C16 | C16 (N156–N159) |
| non-trivial solutions of a homogeneous 2 × 2 system: the determinant must vanish | C16 | primer (in C16) |
| hyperbolic half-angle identities and coth | C16 | primer (in C16) |
| cosh, sinh, tanh | C16 | knowledge/primers.md: hyperbolic functions cosh, sinh, tanh (ch07 P168) — reminder |
| lid conditions; the 2 × 2 system | C16 | C16 (N160, D26) |
| Eady phase speed; complex c means growth | C16 | C16 |
| cut-off, fastest wave, Eady radius NH/f (no π) | C16 | C16 (N161, D27) |
| `scipy.optimize.minimize_scalar` | C16 | knowledge/primers.md: `scipy.optimize.minimize_scalar` (ch07 P170) — reminder |
| Eady growth rate in physical units (ours) | C16 | C16 (N162) |
| Chebyshev eigen-solve as an independent route | C16 | C16 (gloss in the from-scratch row 13) and knowledge/primers.md: Chebyshev–Gauss–Lobatto points and the differentiation matrix (ch11 P258) — reminder |
| energetics: buoyancy flux, poleward heat flux, westward tilt | C16 | C16 (N163) |
| mean of a product of two wave fields | C16 | knowledge/primers.md: mean of a product of real parts (ch07 P178) — reminder |
| three-dimensional cascade and the Kolmogorov scale | C17 | R31 |
| enstrophy: mean-square vorticity, and its spectrum K²S(K) | C17 | primer (in C17) |
| geostrophic turbulence; one-sided spectrum S(K) | C17 | C17 (N164) |
| conservation of energy and enstrophy in two dimensions | C17 | C17 (N165) |
| Fjørtoft's argument | C17 | C17 |
| two inertial ranges; −5/3 and −3 by dimensional analysis | C17 | C17 (N166, D29) |
| Π theorem; `DIM.pi_groups` | C17 | C17 (D29 tools; reminder of the Π-theorem block of ch01) |
| Rhines length | C17 | C17 (N167) |
| the Jacobian J(ψ, ζ) and a pseudo-spectral step with 2/3 de-aliasing | C17 | primer (in C17) |
| decaying two-dimensional turbulence run (ours); spectral centroids; `TS.shell_spectrum` | C17 | C17 (N168 and its "our choice" box) |
| exercises and literature | C17 | C17 (pointers S01, S02 in A.18 row 16) |


---

## Part F — derivation storyboards

One block per `D` row of the curation (§4b), in order. Builders copy each block word for word into `nb.derivation(id, title,
ref=…, goal=…, start=(tex, plain), plan=[…], uses=[…], steps=[dict(did=…, tex=…, why=…, plain=…), …], result=(tex, plain),
interpret=…, check=…, check_src=…)` and, where the heading names an explainer, into that explainer's `derivations: [...]` (same
`did` titles, same step count). Every step has four parts: **did** (the move) · **tex** (the new line) · **why** (why it is
allowed and why we make it) · **plain** (what the line says). Subscripts x, y, z, t on u, v, w, η, p, ψ, ζ denote partial
derivatives only where a line would otherwise not fit a phone, and the first such use in each block says so. Wherever a
step produces a numbered book equation, the number stands directly after that step's **tex** — beside the equation it
names. Colours of terms as header convention 10. The nine ★★★ blocks carry a **sympy check** with a `check_src` listing
that re-runs the derivation's own construction; **each listing was executed as written for this design and passed** (the
printed line is quoted under it). Hemisphere: steps are written for f > 0 as the book does; the last step or the Check of
every block where direction matters gives the f < 0 form.

### D01 · The friction force (13.6) from the anisotropic eddy stress (13.5) — ★★, 6 steps, in C01 (notebook)
- **Goal:** turn the six eddy stresses into a friction force that can stand in the momentum equation — and see why it comes
  out so simple. The book says only that the stresses "become" it. · **Start:** $F_i=\dfrac1\rho\dfrac{\partial\tau_{ij}}{\partial x_j}$ (slip #2 corrected: per unit mass)
  with the stresses $\tau_{xx}=2\rho\nu_H\dfrac{\partial u}{\partial x}$, $\tau_{xy}=\rho\nu_H\Big(\dfrac{\partial u}{\partial y}+\dfrac{\partial v}{\partial x}\Big)$, $\tau_{xz}=\rho\nu_v\dfrac{\partial u}{\partial z}+\rho\nu_H\dfrac{\partial w}{\partial x}$, three of the six members of (13.5) —
  "the force on a parcel is the net stress on its faces". · **Plan:** • write out the sum over the three faces • insert the
  stresses and differentiate • collect the terms that form the divergence of the velocity • remove them with continuity. ·
  **Tools:** index notation and the summation convention (Ch. 2, reminder); second partial derivatives and Schwarz's theorem
  (P121, reminder); continuity $\nabla\cdot\mathbf u=0$ (R04). · **Assumptions:** ρ, ν_H, ν_v uniform (step 2); incompressible (step 5).
- **Steps:**
  1. **did** Write the x-component of the force · **tex** $F_x=\dfrac1\rho\Big(\dfrac{\partial\tau_{xx}}{\partial x}+\dfrac{\partial\tau_{xy}}{\partial y}+\dfrac{\partial\tau_{xz}}{\partial z}\Big)$ · **why** The repeated index j is
     summed over x, y and z; dividing by ρ turns a force per volume into a force per mass, which is what the momentum
     equation needs. · **plain** The net push on a parcel is the difference of the stresses on opposite faces, face by face.
  2. **did** Insert the three stresses · **tex** $F_x=\dfrac{\partial}{\partial x}\Big(2\nu_H\dfrac{\partial u}{\partial x}\Big)+\dfrac{\partial}{\partial y}\Big[\nu_H\Big(\dfrac{\partial u}{\partial y}+\dfrac{\partial v}{\partial x}\Big)\Big]+\dfrac{\partial}{\partial z}\Big(\nu_v\dfrac{\partial u}{\partial z}+\nu_H\dfrac{\partial w}{\partial x}\Big)$ ·
     **why** Substitution of the Start; ρ is uniform, so it passes through the derivatives and cancels the 1/ρ in front. ·
     **plain** Each stress is an eddy viscosity times a velocity gradient.
  3. **did** Differentiate term by term · **tex** $F_x=2\nu_H\dfrac{\partial^2u}{\partial x^2}+\nu_H\dfrac{\partial^2u}{\partial y^2}+\nu_H\dfrac{\partial^2v}{\partial y\,\partial x}+\nu_v\dfrac{\partial^2u}{\partial z^2}+\nu_H\dfrac{\partial^2w}{\partial z\,\partial x}$ · **why** The two eddy
     viscosities are constants, so only the velocities are differentiated (linearity of the derivative). · **plain** Five
     second derivatives: three of u, one of v, one of w.
  4. **did** Split the first term and regroup · **tex** $F_x=\nu_H\Big(\dfrac{\partial^2u}{\partial x^2}+\dfrac{\partial^2u}{\partial y^2}\Big)+\nu_v\dfrac{\partial^2u}{\partial z^2}+\nu_H\dfrac{\partial}{\partial x}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}\Big)$ · **why**
     Mixed partial derivatives may be taken in either order (Schwarz), so the three leftover terms share an outer ∂/∂x;
     we group them because their bracket is the divergence. · **plain** Everything except a Laplacian-like part is the
     x-derivative of the velocity divergence.
  5. **did** Use continuity · **tex** $F_x=\nu_H\Big(\dfrac{\partial^2u}{\partial x^2}+\dfrac{\partial^2u}{\partial y^2}\Big)+\nu_v\dfrac{\partial^2u}{\partial z^2}$ (13.6), first member · **why** The fluid is incompressible,
     so the bracket of step 4 is zero everywhere and so is its derivative. · **plain** Horizontal mixing with the large
     coefficient, vertical mixing with the small one.
  6. **did** Repeat for the other components · **tex** $F_z=\nu_H\Big(\dfrac{\partial^2w}{\partial x^2}+\dfrac{\partial^2w}{\partial y^2}\Big)+\nu_v\dfrac{\partial^2w}{\partial z^2}+\nu_v\dfrac{\partial}{\partial z}(\nabla\cdot\mathbf u)$, last term zero ·
     **why** The same four moves with the stresses of the y- and z-faces; in the z-component the leftover terms carry
     ν_v, and again form the divergence. · **plain** All three components have the same form: ν_H on horizontal second
     derivatives, ν_v on vertical ones.
- **Result:** $F_x=\nu_H\Big(\dfrac{\partial^2u}{\partial x^2}+\dfrac{\partial^2u}{\partial y^2}\Big)+\nu_v\dfrac{\partial^2u}{\partial z^2}$, and the same with v and with w — the three members of (13.6), per unit mass:
  "anisotropic diffusion of momentum".
- **Check:** units — ν [m²/s] × velocity/length² = m/s² ✓ (the printed form without 1/ρ would be N/m³ on one side and m/s² on
  the other). Limit — ν_H = ν_v = ν gives $\nu\nabla^2u$, the Newtonian term of Chapter 4. Code —
  `ch13.eddy_friction_force_sympy()["ok"]` is True.
- **What it means:** friction acts like diffusion of each velocity component, fast sideways and slow vertically. **Fails when:**
  the eddy viscosities vary in space (extra terms with their gradients — the Ekman layer with $K_v(z)$ of C05 keeps
  them) or the flow is compressible.
- **Traps:** forgetting the 1/ρ (slip #2); dropping the cross terms without saying that continuity removes them (trap T15);
  writing ν_v on the horizontal derivatives of w.

### D02 · The Coriolis components (13.7), the Coriolis parameter (13.8) and the thin-shell equations (13.9) — ★★, 9 steps, in C01 (notebook)
- **Goal:** reduce the rotating Boussinesq momentum equation to the form every later section uses, and see exactly which
  terms are thrown away and why. · **Start:** $\dfrac{D\mathbf u}{Dt}+2\boldsymbol\Omega\times\mathbf u=-\dfrac1{\rho_0}\nabla p-\dfrac{g\rho}{\rho_0}\mathbf e_z+\mathbf F$, the momentum member of (13.2)
  with its sign corrected (slip #1), p and ρ being perturbations — "acceleration + Coriolis = pressure + buoyancy +
  friction". · **Plan:** • resolve Ω on the tangent plane • take the cross product • use thinness twice (once on each
  dropped term) • write the three components. · **Tools:** components of a vector; the cross product as a determinant (P53,
  reminder); scaling with two lengths (P188, reminder); N12's $W/U\sim H/L$. · **Assumptions:** thin layer, H ≪ L (steps 3–4, 6);
  not at the equator (step 4).
- **Steps:**
  1. **did** Resolve the rotation vector locally · **tex** $\boldsymbol\Omega=(0,\ \Omega\cos\theta,\ \Omega\sin\theta)$ · **why** At latitude θ the earth's axis lies in
     the north–vertical plane, at angle θ above the northward horizontal; there is no eastward part. We need components
     along x (east), y (north), z (up). · **plain** Part of the spin is about the local vertical, part about the northward
     horizontal.
  2. **did** Take the cross product · **tex** $2\boldsymbol\Omega\times\mathbf u=2\Omega\big[\mathbf e_x(w\cos\theta-v\sin\theta)+\mathbf e_yu\sin\theta-\mathbf e_zu\cos\theta\big]$ · **why** Expand the
     determinant with rows (e_x, e_y, e_z), (0, 2Ω cos θ, 2Ω sin θ), (u, v, w); this is the exact Coriolis acceleration. ·
     **plain** Four Coriolis terms: two in the east equation, one each in the north and vertical equations.
  3. **did** Compare w with v · **tex** $\dfrac{w\cos\theta}{v\sin\theta}\sim\dfrac{H}{L}\cot\theta\ll1$ · **why** Continuity makes the vertical velocity smaller
     than the horizontal by the aspect ratio, $W/U\sim H/L$ (N12), about 10⁻² or less; cot θ is of order one away from the
     equator. · **plain** In a thin layer the Coriolis term with w is a hundred times smaller than the one with v.
  4. **did** Drop the small term · **tex** $2\boldsymbol\Omega\times\mathbf u\cong(-2\Omega v\sin\theta,\ 2\Omega u\sin\theta,\ -2\Omega u\cos\theta)$ · **why** Step 3 justifies neglecting
     w cos θ beside v sin θ (≅ marks the approximation); it fails within a few degrees of the equator, where sin θ → 0. ·
     **plain** Horizontal motion is turned only by the vertical part of the earth's spin.
  5. **did** Name the recurring factor · **tex** $f=2\Omega\sin\theta$ (13.8), so that $2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)$ (13.7) · **why** A
     definition; it saves writing and names the one number through which rotation enters horizontal dynamics. ·
     **plain** The Coriolis parameter is twice the local vertical rotation rate; it changes sign across the equator.
  6. **did** Compare the vertical Coriolis term with buoyancy · **tex** $\dfrac{2\Omega u\cos\theta}{g\rho/\rho_0}\sim\dfrac{10^{-3}\ \mathrm{m\,s^{-2}}}{10^{-2}\ \mathrm{m\,s^{-2}}}\ll1$ · **why** For winds of
     order ten metres per second the term is about 10⁻³ m/s², against a buoyancy of 10⁻² m/s² for a density
     perturbation of one part in a thousand. So it is dropped. · **plain** In the vertical, rotation is negligible next to weight
     and pressure.
  7. **did** Write the x-component · **tex** $\dfrac{Du}{Dt}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}+\nu_H\Big(\dfrac{\partial^2u}{\partial x^2}+\dfrac{\partial^2u}{\partial y^2}\Big)+\nu_v\dfrac{\partial^2u}{\partial z^2}$ (13.9), first member · **why** Take the
     x-component of the Start with step 5 and the friction force of D01. · **plain** Eastward acceleration, minus f times
     the northward velocity, equals the pressure push plus friction.
  8. **did** Write the y-component · **tex** $\dfrac{Dv}{Dt}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}+\nu_H\Big(\dfrac{\partial^2v}{\partial x^2}+\dfrac{\partial^2v}{\partial y^2}\Big)+\nu_v\dfrac{\partial^2v}{\partial z^2}$ (13.9), second member · **why** The same,
     with the y-component fu of the Coriolis acceleration; note the sign opposite to the x-equation. · **plain** The
     Coriolis terms couple the two horizontal equations with opposite signs — that is what makes things turn.
  9. **did** Write the z-component · **tex** $\dfrac{Dw}{Dt}=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial z}-\dfrac{g\rho}{\rho_0}+\nu_H\Big(\dfrac{\partial^2w}{\partial x^2}+\dfrac{\partial^2w}{\partial y^2}\Big)+\nu_v\dfrac{\partial^2w}{\partial z^2}$ (13.9), third member · **why** The
     z-component with the Coriolis term dropped by step 6; buoyancy stays because it is one of the two largest terms. ·
     **plain** No rotation in the vertical equation at all.
- **Result:** the three members of (13.9): $\dfrac{Du}{Dt}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}+F_x$, $\dfrac{Dv}{Dt}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}+F_y$, $\dfrac{Dw}{Dt}=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial z}-\dfrac{g\rho}{\rho_0}+F_z$ —
  "rotation enters through f alone, and only in the horizontal".
- **Check:** units — f [1/s] × velocity = m/s² ✓. Limits — at the pole f = 2Ω and the layer turns like a turntable; at the
  equator f = 0 and step 3 fails. Number — 35° N: f = 8.365 × 10⁻⁵ s⁻¹; dropped/kept in step 3 = 0.9 % for U = 14 m/s, H = 9 km,
  L = 1400 km. Code — `ch13.thin_layer_sympy()["ok"]`.
- **What it means:** every later block is this set with more terms switched off. Dropping the Ω cos θ terms is the
  "traditional approximation" (a name the book does not use). **Fails when:** near the equator, or for deep convection
  where w is not small.
- **Traps:** the sign of fu in the y-equation; which term is dropped against which; from here on p and ρ are perturbations
  (trap T3).

### D03 · Geostrophic balance (13.11)–(13.12) from the thin-shell equations — ★, 6 steps, in C02 (notebook · `geostrophic_balance`)
- **Goal:** find what is left of the horizontal momentum equations for slow, large-scale flow, and what that says about
  the direction of the wind. · **Start:** $\dfrac{Du}{Dt}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}+F_x$ and $\dfrac{Dv}{Dt}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}+F_y$, the horizontal members of (13.9) —
  "acceleration + Coriolis = pressure force + friction". · **Plan:** • estimate the size of each term • form two ratios •
  drop what is small • read the direction of the wind from what is left. · **Tools:** order-of-magnitude scaling (P130,
  reminder); the dot product of two vectors. · **Assumptions:** Ro ≪ 1 and time scale ≫ 1/f (step 3); E ≪ 1 (step 3);
  constant f (step 6).
- **Steps:**
  1. **did** Estimate the size of each term · **tex** $\dfrac{Du}{Dt}\sim\dfrac{U^2}{L},\qquad fv\sim fU,\qquad F_x\sim\dfrac{\nu_vU}{H^2}$ · **why** Replace each derivative by
     a typical change over a typical distance: velocity U, horizontal scale L, vertical scale H. We need sizes, not
     values. · **plain** Three candidate terms to balance the pressure force.
  2. **did** Divide by the Coriolis term · **tex** $\mathrm{Ro}=\dfrac{U^2/L}{fU}=\dfrac{U}{fL}$ (13.13), and $E=\dfrac{\nu_v}{fH^2}$ · **why** Ratios are dimensionless, so
     "small" has a meaning; these two are the Rossby number and the Ekman number. · **plain** Ro compares acceleration
     with Coriolis, E compares friction with Coriolis.
  3. **did** Assume both are small · **tex** $\mathrm{Ro}\ll1,\quad E\ll1$ · **why** For weather systems and ocean gyres away from boundaries
     Ro is about 0.1 and E far smaller; the flow must also change slowly compared with 1/f, or ∂u/∂t would not be
     small. · **plain** Slow, wide, frictionless motion.
  4. **did** Drop the small terms · **tex** $-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$ and $fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$ (13.11)–(13.12) · **why** What is left must balance
     to the accuracy of Ro and E: only the Coriolis term can match the pressure-gradient force. · **plain** Coriolis
     force and pressure force cancel.
  5. **did** Dot the velocity with the pressure gradient · **tex** $\mathbf u\cdot\nabla p=\dfrac1{\rho_0f}\Big(-\dfrac{\partial p}{\partial y}\dfrac{\partial p}{\partial x}+\dfrac{\partial p}{\partial x}\dfrac{\partial p}{\partial y}\Big)=0$ · **why** Insert u and v
     from step 4; a zero dot product means the two vectors are perpendicular. · **plain** The wind blows along the
     isobars, never across them.
  6. **did** Define a stream function · **tex** $\psi=\dfrac{p}{f\rho_0}:\quad u=-\dfrac{\partial\psi}{\partial y},\quad v=\dfrac{\partial\psi}{\partial x}$ · **why** For constant f, step 4 has exactly this
     form; such a flow has no horizontal divergence, and lines of constant ψ — the isobars — are its streamlines. ·
     **plain** A pressure map is a map of the flow.
- **Result:** $-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$, $fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$ — "the wind follows the isobars with low pressure on its left for f > 0 (on its
  right for f < 0), at a speed proportional to their spacing".
- **Check:** units — fU [m/s²] against Δp/(ρ₀L) [Pa/(kg m⁻²)] = m/s² ✓. Number — 4 hPa over 300 km at 35° N, ρ₀ = 1.2 kg/m³:
  13.3 m/s. Limits — f → 0 gives infinite wind: the balance cannot hold at the equator. Code —
  `GFD.geostrophic_from_field` and the dot product at round-off.
- **What it means:** the balance is *diagnostic*: it gives the wind from the pressure but cannot say how either changes
  (that needs the small terms — C15, C16). **Fails when:** near the equator; in boundary layers (C06); in fast or small
  systems (Ro of order 1).
- **Traps:** the stream-function sign is opposite to Chapters 4 and 11 (trap T14); ψ exists only for constant f; "low on
  the left" flips with the sign of f; the Rossby number of Chapter 4 used 2Ω, not f (trap T17).

### D04 · Thermal wind (13.15), and its temperature form (ours) — ★★, 7 steps, in C03 (notebook · `thermal_wind`)
- **Goal:** show that a horizontal density (temperature) gradient forces the geostrophic wind to change with height. The
  book says "eliminating p". · **Start:** $-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$ and $fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$ (13.11)–(13.12), with $0=-\dfrac{\partial p}{\partial z}-g\rho$ (13.14) —
  "geostrophic in the horizontal, hydrostatic in the vertical". · **Plan:** • differentiate the geostrophic equations in z •
  differentiate the hydrostatic equation in x and in y • equate the mixed derivatives of p • (ours) replace density by
  temperature. · **Tools:** cross-differentiation (P177) and Schwarz's theorem (P121), reminders; the linear equation of
  state (R06). · **Assumptions:** f does not depend on z (steps 1, 4); Boussinesq (ρ₀ constant); steps 6–7 need
  $\rho'=-\rho_0\alpha T'$.
- **Steps:**
  1. **did** Differentiate the first balance in z · **tex** $-f\dfrac{\partial v}{\partial z}=-\dfrac1{\rho_0}\dfrac{\partial^2p}{\partial z\,\partial x}$ · **why** Both sides are differentiable and f and
     ρ₀ do not depend on z, so they pass through; we want the shear of v on the left. · **plain** The change of wind with
     height is tied to the change with height of the pressure gradient.
  2. **did** Differentiate the hydrostatic balance in x · **tex** $\dfrac{\partial^2p}{\partial x\,\partial z}=-g\dfrac{\partial\rho}{\partial x}$ · **why** Gravity g is constant; this expresses the
     same mixed derivative of p through the density. · **plain** Where the fluid is denser, pressure falls faster with
     height.
  3. **did** Equate the mixed derivatives · **tex** $\dfrac{\partial v}{\partial z}=-\dfrac{g}{\rho_0f}\dfrac{\partial\rho}{\partial x}$ (13.15), first member · **why** The order of differentiation does
     not matter (Schwarz), so steps 1 and 2 share $\partial^2p/\partial x\partial z$; eliminate it and divide by −f. · **plain** Density
     increasing eastward makes the northward wind decrease with height (f > 0).
  4. **did** Do the same with the second balance · **tex** $f\dfrac{\partial u}{\partial z}=-\dfrac1{\rho_0}\dfrac{\partial^2p}{\partial z\,\partial y},\qquad\dfrac{\partial^2p}{\partial y\,\partial z}=-g\dfrac{\partial\rho}{\partial y}$ · **why** The same two
     differentiations applied to the y-equation and to the hydrostatic balance in y. · **plain** The same link in the
     north–south direction.
  5. **did** Eliminate the pressure again · **tex** $\dfrac{\partial u}{\partial z}=\dfrac{g}{\rho_0f}\dfrac{\partial\rho}{\partial y}$ (13.15), second member · **why** Equate the mixed derivatives
     and divide by f; note the plus sign, inherited from the plus sign of fu. · **plain** Density increasing northward
     makes the eastward wind increase with height (f > 0).
  6. **did** Express density through temperature (ours) · **tex** $\dfrac{\partial\rho}{\partial y}=-\rho_0\alpha\dfrac{\partial T}{\partial y}$ · **why** For small departures the
     equation of state is linear, $\rho'=-\rho_0\alpha T'$ (R06), with α = 1/T₀ for a perfect gas; meteorologists measure
     temperature, not density. · **plain** Cold means dense.
  7. **did** Substitute (ours) · **tex** $\dfrac{\partial u}{\partial z}=-\dfrac{g\alpha}{f}\dfrac{\partial T}{\partial y},\qquad\dfrac{\partial v}{\partial z}=\dfrac{g\alpha}{f}\dfrac{\partial T}{\partial x}$ · **why** Insert step 6 into steps 5 and 3;
     the minus sign of step 6 flips both signs. · **plain** Temperature falling toward the pole makes the westerly wind
     grow with height.
- **Result:** $\dfrac{\partial v}{\partial z}=-\dfrac{g}{\rho_0f}\dfrac{\partial\rho}{\partial x}$, $\dfrac{\partial u}{\partial z}=\dfrac{g}{\rho_0f}\dfrac{\partial\rho}{\partial y}$ — the two members of (13.15); in temperature (ours)
  $\dfrac{\partial u}{\partial z}=-\dfrac{g\alpha}{f}\dfrac{\partial T}{\partial y}$ — "a horizontal temperature gradient cannot coexist with a wind that is the same at every height".
- **Check:** units — (g/ρ₀f) ∂ρ/∂y = (m s⁻²)(kg m⁻⁴)/(kg m⁻³ s⁻¹) = s⁻¹ ✓. Number — ∂T/∂y = −7 K per 1000 km, α = 1/280 K⁻¹, 35° N:
  2.93 × 10⁻³ s⁻¹, i.e. 26.4 m/s over 9 km. Hemisphere — f < 0 with ∂T/∂y > 0 gives the same westerly shear. Limit — no
  gradient, no shear (D05).
- **What it means:** the jet stream sits above the strongest temperature contrast; an oceanographer gets currents from a
  density section. Only the *shear* is fixed — one reference wind must come from elsewhere. **Fails when:** geostrophy or
  hydrostatics fail (equator, small scales).
- **Traps:** which equation is differentiated with respect to which variable; the minus of the first member and the plus
  of the second; ρ is the perturbation (trap T3); in temperature form the signs flip again and α must be supplied.

### D05 · The Taylor–Proudman theorem (13.21) — ★★, 7 steps, in C03 (notebook)
- **Goal:** show that slow, steady, frictionless motion of a *homogeneous* rotating fluid cannot vary along the rotation
  axis. The book prints one of its starting equations with the wrong sign and does not write the cross-differentiation
  out. · **Start:** $-2\Omega v=-\dfrac1\rho\dfrac{\partial p}{\partial x}$ (13.16) and $2\Omega u=-\dfrac1\rho\dfrac{\partial p}{\partial y}$ (13.17, sign corrected — slip #3), with $0=-\dfrac{\partial p}{\partial z}-g\rho$
  (13.14) and $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$ — "geostrophy with f = 2Ω, hydrostatics, continuity". · **Plan:** • cross-differentiate
  the two horizontal equations to kill p • use continuity • differentiate in z and use uniform density. · **Tools:**
  cross-differentiation (P177, reminder); continuity; the hypotheses of a theorem (P255, reminder). · **Assumptions:**
  steady, Ro ≪ 1, E ≪ 1 (the Start); uniform density (step 5).
- **Steps:**
  1. **did** Differentiate the first equation in y · **tex** $-2\Omega\dfrac{\partial v}{\partial y}=-\dfrac1\rho\dfrac{\partial^2p}{\partial y\,\partial x}$ · **why** Both Ω and ρ are constants here; we aim at
     the mixed derivative of p so that it can be eliminated. · **plain** One equation now contains ∂²p/∂x∂y.
  2. **did** Differentiate the second in x · **tex** $2\Omega\dfrac{\partial u}{\partial x}=-\dfrac1\rho\dfrac{\partial^2p}{\partial x\,\partial y}$ · **why** The same move on the other equation
     produces the same mixed derivative (Schwarz). · **plain** So does the other.
  3. **did** Subtract step 1 from step 2 · **tex** $2\Omega\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$ · **why** The right-hand sides are equal and cancel.
     With the book's printed sign of the second equation this step would give $\partial u/\partial x-\partial v/\partial y=0$ instead. · **plain** The
     horizontal flow has no divergence.
  4. **did** Use continuity · **tex** $\dfrac{\partial w}{\partial z}=0$ (13.19) · **why** Continuity says the three terms sum to zero; the first two already
     do, so the third vanishes by itself. · **plain** Columns are neither stretched nor squashed.
  5. **did** Differentiate the first equation in z · **tex** $-2\Omega\dfrac{\partial v}{\partial z}=-\dfrac1\rho\dfrac{\partial}{\partial x}\Big(\dfrac{\partial p}{\partial z}\Big)=\dfrac g\rho\dfrac{\partial\rho}{\partial x}=0$ · **why** Swap the order of
     the derivatives, insert the hydrostatic balance, and use the uniform density: its horizontal gradient is zero. ·
     **plain** No density gradient, no thermal wind.
  6. **did** Do the same for the second · **tex** $\dfrac{\partial v}{\partial z}=\dfrac{\partial u}{\partial z}=0$ (13.20) · **why** The identical move on the second equation gives
     ∂u/∂z = 0. · **plain** The horizontal velocity is the same at every height.
  7. **did** Combine · **tex** $\dfrac{\partial\mathbf u}{\partial z}=0$ (13.21) · **why** Steps 4 and 6 together cover all three velocity components. · **plain**
     The whole velocity field is independent of the coordinate along the rotation axis.
- **Result:** $\partial\mathbf u/\partial z=0$ — "the flow is two-dimensional: fluid moves in columns parallel to the rotation axis".
- **Check:** limit — it is the thermal wind (D04) with the density gradient set to zero. Code —
  `ch13.taylor_proudman_sympy()["ok"]` True and `(printed=True)["ok"]` False. Number — `GFD.taylor_proudman_residual` on rigid
  rotation returns zeros.
- **What it means:** a short obstacle on the floor of a rapidly rotating tank carries a whole column of fluid with it (a
  Taylor column). In the stratified atmosphere and ocean the density gradient prevents this — which is why D04, not D05,
  describes them. **Fails when:** any of the four hypotheses fails: unsteady, Ro or E not small, density not uniform.
- **Traps:** the printed sign (slip #3); "does not vary along the axis" is not "no motion along the axis" (w may be a
  non-zero constant in z); uniform density enters only at step 5 (trap T15).


### D06 · The surface Ekman layer: (13.22)–(13.23) → (13.27) → (13.28)–(13.29) → the constant A and u(z), v(z) — ★★★, 14 steps, in C04 (notebook · `ekman_spiral`)
- **Goal:** find the current under a steady wind when the only thing that can balance friction is the Coriolis force — and
  with it the thickness of the layer and the angle of the surface current. · **Start:** $-fv=\nu_v\dfrac{d^2u}{dz^2}$ and $fu=\nu_v\dfrac{d^2v}{dz^2}$
  (13.22)–(13.23), with $\rho\nu_v\dfrac{du}{dz}=\tau$ and $\dfrac{dv}{dz}=0$ at z = 0, and $u,v\to0$ as $z\to-\infty$ (slip #4 corrected) — "Coriolis against
  vertical friction; a stress τ along x on top; nothing moving at depth". · **Plan:** • pack u and v into one complex
  velocity • solve the resulting equation with exponentials • keep the one that dies away downward • fix its size with
  the wind stress • unpack real and imaginary parts. · **Tools:** complex velocity (P311); trying $e^{rz}$ in a linear ODE
  (P44, reminder); √i (P159, reminder); polar form and Euler's formula $e^{i\phi}=\cos\phi+i\sin\phi$ (P153, reminder). ·
  **Assumptions:** steady, horizontally uniform, no pressure gradient (the Start); constant ν_v (step 3); infinitely deep
  ocean (step 7); f > 0 until step 14.
- **Steps:**
  1. **did** Multiply the second equation by i · **tex** $ifu=\nu_v\dfrac{d^2(iv)}{dz^2}$ · **why** Multiplying an equation by a constant is always
     allowed; i is chosen so that adding the first equation will produce u + iv on both sides. · **plain** A preparation
     step: it lines the two equations up.
  2. **did** Add the first equation · **tex** $\dfrac{d^2V}{dz^2}=\dfrac{if}{\nu_v}V$, $V\equiv u+iv$ (13.27) · **why** The left sides add to $-fv+ifu=if(u+iv)$,
     the right sides to $\nu_v(u+iv)''$; one complex equation now carries both real ones. · **plain** Friction on the velocity
     arrow equals f times the arrow turned by 90°.
  3. **did** Try an exponential · **tex** $V=e^{rz}\ \Rightarrow\ r^2=\dfrac{if}{\nu_v}$ · **why** The equation is linear with constant coefficients, so
     exponentials solve it; the exponential cancels because it is never zero. · **plain** The problem becomes: find the
     square root of an imaginary number.
  4. **did** Take the square root · **tex** $r=\pm\sqrt{i}\,\sqrt{\dfrac f{\nu_v}}=\pm\dfrac{1+i}{\sqrt2}\sqrt{\dfrac f{\nu_v}}$ · **why** Because $(1+i)^2=2i$, the square root
     of i is $(1+i)/\sqrt2$; this needs f > 0. · **plain** The exponent has equal real and imaginary parts.
  5. **did** Name the length scale · **tex** $\delta=\sqrt{\dfrac{2\nu_v}{f}}\ \Rightarrow\ r=\pm\dfrac{1+i}{\delta}$ (13.29) · **why** A definition that absorbs the factor
     $\sqrt{f/2\nu_v}$; δ will turn out to be the thickness of the layer. · **plain** One length sets both the decay and the
     turning.
  6. **did** Write the general solution · **tex** $V=A\,e^{(1+i)z/\delta}+B\,e^{-(1+i)z/\delta}$ (13.28) · **why** A second-order linear equation has two
     independent solutions; any solution is a combination with complex constants A and B. · **plain** One part shrinks
     going down, the other grows.
  7. **did** Discard the part that grows downward · **tex** $B=0\ \Rightarrow\ V=A\,e^{(1+i)z/\delta}$ · **why** As $z\to-\infty$ the factor $e^{-z/\delta}$ grows
     without bound, but the current must vanish at depth; so B must be zero. · **plain** The wind's influence dies away
     with depth.
  8. **did** Combine the two surface conditions · **tex** $\rho\nu_v\dfrac{dV}{dz}=\tau$ at $z=0$ · **why** Add the first condition to i times the
     second, $\rho\nu_v(u'+iv')=\tau+i\cdot0$; the stress, like the velocity, becomes one complex number. · **plain** The wind
     pushes along x only.
  9. **did** Apply it to the solution · **tex** $\rho\nu_vA\,\dfrac{1+i}{\delta}=\tau$ · **why** Differentiate step 7, $V'=A\frac{1+i}\delta e^{(1+i)z/\delta}$, and set
     z = 0, where the exponential is 1. · **plain** One equation for the one remaining constant.
  10. **did** Solve for A · **tex** $A=\dfrac{\tau\delta}{\rho\nu_v(1+i)}=\dfrac{\tau\delta(1-i)}{2\rho\nu_v}$ · **why** Multiply numerator and denominator by the conjugate
      1 − i, since $(1+i)(1-i)=2$; this is the constant the book quotes. · **plain** The surface velocity has equal parts
      along the wind and to its right.
  11. **did** Write A in polar form · **tex** $A=\dfrac{\tau/\rho}{\sqrt{f\nu_v}}\,e^{-i\pi/4}$ · **why** Use $1-i=\sqrt2e^{-i\pi/4}$ and $\dfrac{\delta\sqrt2}{2\nu_v}=\dfrac1{\sqrt{f\nu_v}}$ (from the
      definition of δ); a length and an angle are easier to read than two components. · **plain** The surface current
      has speed (τ/ρ)/√(fν_v) and points 45° to the right of the wind.
  12. **did** Insert A into the solution · **tex** $V=\dfrac{\tau/\rho}{\sqrt{f\nu_v}}\,e^{z/\delta}\,e^{i(z/\delta-\pi/4)}$ · **why** Split $e^{(1+i)z/\delta}=e^{z/\delta}e^{iz/\delta}$ and add the
      two phase angles. · **plain** Going down, the speed falls by e and the direction turns by one radian every δ.
  13. **did** Take real and imaginary parts · **tex** $u=\dfrac{\tau/\rho}{\sqrt{f\nu_v}}e^{z/\delta}\cos\Big(-\dfrac z\delta+\dfrac\pi4\Big),\quad v=-\dfrac{\tau/\rho}{\sqrt{f\nu_v}}e^{z/\delta}\sin\Big(-\dfrac z\delta+\dfrac\pi4\Big)$ · **why**
      Euler's formula with $\cos(-x)=\cos x$ and $\sin(-x)=-\sin x$ gives the book's printed form. · **plain** At the
      surface u = −v > 0: 45° to the right of the wind; deeper, the arrow turns clockwise.
  14. **did** Repeat for f < 0 · **tex** $V=\dfrac{\tau/\rho}{\sqrt{\lvert f\rvert\nu_v}}\,e^{z/\delta}\,e^{-i(z/\delta-\pi/4)},\quad\delta=\sqrt{\dfrac{2\nu_v}{\lvert f\rvert}}$ · **why** With f negative, $r^2=-i\lvert f\rvert/\nu_v$
      and $\sqrt{-i}=(1-i)/\sqrt2$: every i in steps 4–12 changes sign, i.e. v changes sign. · **plain** In the southern
      hemisphere the spiral is the mirror image: surface current 45° to the left, turning counter-clockwise downward.
- **Result:** $V=u+iv=\dfrac{\tau\delta(1-i)}{2\rho\nu_v}e^{(1+i)z/\delta}$ with $\delta=\sqrt{2\nu_v/f}$ — "a spiral of currents that fades and turns with depth, 45°
  to the right of the wind at the surface (f > 0)".
- **Check:** units — (τ/ρ)/√(fν_v) = (m²/s²)/(m/s) = m/s ✓. Limits — ν_v → 0: δ → 0 and the surface speed → ∞ (a thinner
  layer must move faster to carry the same transport, D07); f → 0: δ → ∞, no steady layer. Number — 60° N, τ = 0.07 N/m², ν_v
  = 0.03 m²/s, ρ = 1027 kg/m³: δ = 21.80 m, surface current (0.02476, −0.02476) m/s.
- **sympy check** (`check_src`; executed as written — every assertion passes and the final line is printed):
  ```python
  import sympy as sp                                            # symbolic algebra
  z = sp.symbols("z", real=True)                                # depth coordinate, z <= 0 in the ocean
  f, nu, tau, rho = sp.symbols("f nu_v tau rho", positive=True) # northern hemisphere: f > 0
  delta = sp.sqrt(2 * nu / f)                                   # step 5: the Ekman thickness
  lam = (1 + sp.I) / delta                                      # steps 4-5: the root that decays downward
  assert sp.simplify(lam**2 - sp.I * f / nu) == 0               # step 3: r^2 = i f / nu_v
  A = tau * delta * (1 - sp.I) / (2 * rho * nu)                 # step 10: the constant from the stress condition
  V = A * sp.exp(lam * z)                                       # step 7: the bounded solution V = u + i v
  assert sp.simplify(sp.diff(V, z, 2) - sp.I * f / nu * V) == 0 # it solves the complex equation of step 2
  assert sp.simplify(rho * nu * sp.diff(V, z).subs(z, 0) - tau) == 0   # step 8: stress tau along x, none along y
  amp = (tau / rho) / sp.sqrt(f * nu) * sp.exp(z / delta)       # step 11: the amplitude of the printed form
  u_book = amp * sp.cos(-z / delta + sp.pi / 4)                 # the book's u(z)
  v_book = -amp * sp.sin(-z / delta + sp.pi / 4)                # the book's v(z)
  gap = sp.simplify(sp.expand_complex(V - (u_book + sp.I * v_book)))   # step 13: real and imaginary parts agree?
  assert gap == 0
  print("D06 checks passed: V'' = (i f/nu) V, surface stress, and the printed cos/sin form")
  ```
- **What it means:** rotation, not time or distance, sets the thickness of this boundary layer; the eddy viscosity sets
  how deep and how fast, the wind stress how strong. **Fails when:** ν_v varies with depth (real mixed layers: C05, N37),
  the water is shallower than a few δ (N34), or the wind changes within an inertial period.
- **Traps:** i × (second) + (first), not the other way round; keeping the root that grows downward (the book's printed
  $z\to\infty$, slip #4, invites it); $1/(1+i)=(1-i)/2$; the printed cos and sin hold for f > 0 only (trap T5); δ is not the
  "Ekman depth" πδ (trap T6).

### D07 · Ekman transport (13.30), two ways — ★★, 8 steps, in C05 (notebook · `ekman_spiral`)
- **Goal:** add up the currents of the whole layer — and discover that the answer does not depend on the eddy viscosity. ·
  **Start:** $V=A\,e^{(1+i)z/\delta}$ with $A=\dfrac{\tau\delta(1-i)}{2\rho\nu_v}$ (the Result of D06) — "the spiral". · **Plan:** • integrate the complex
  velocity over depth • simplify the constants • read off the two components • do it again without the solution, by
  integrating the momentum balance once. · **Tools:** the integral of an exponential over a half-line (P144, reminder);
  transport per unit width (P313); first integral of an ODE (P312). · **Assumptions:** the stress vanishes at depth (step
  7) — nothing else in the second route.
- **Steps:**
  1. **did** Define the transport as an integral · **tex** $M\equiv\displaystyle\int_{-\infty}^0V\,dz=A\displaystyle\int_{-\infty}^0e^{(1+i)z/\delta}\,dz$ · **why** Transport per unit width is
     velocity summed over depth; A is constant and leaves the integral. · **plain** Add all the arrows of the spiral.
  2. **did** Integrate the exponential · **tex** $M=\dfrac{A\delta}{1+i}$ · **why** An antiderivative of $e^{az}$ is $e^{az}/a$ with $a=(1+i)/\delta$; it is
     1/a at z = 0 and tends to zero as z → −∞ because the real part of a is positive. · **plain** The sum is finite: the
     arrows shrink fast enough.
  3. **did** Insert A · **tex** $M=\dfrac{\tau\delta^2(1-i)}{2\rho\nu_v(1+i)}$ · **why** Substitution of the constant found in D06. · **plain** Everything is
     now in terms of the wind stress and the layer's properties.
  4. **did** Simplify · **tex** $M=-\dfrac{i\,\tau}{\rho f}$ · **why** Two facts: $(1-i)/(1+i)=-i$ and $\delta^2=2\nu_v/f$; the eddy viscosity cancels
     exactly. · **plain** The transport is purely imaginary: all of it is in the y-direction.
  5. **did** Separate the components · **tex** $\displaystyle\int_{-\infty}^0u\,dz=0,\qquad\displaystyle\int_{-\infty}^0v\,dz=-\dfrac{\tau}{\rho f}$ (13.30) · **why** Real part and imaginary
     part of M = ∫u dz + i∫v dz. · **plain** No net transport along the wind; τ/(ρf) at right angles to it, to the right
     for f > 0.
  6. **did** Start again from the momentum balance in stress form · **tex** $-\rho fv=\dfrac{d\tau_x}{dz},\qquad\rho fu=\dfrac{d\tau_y}{dz}$ · **why** This is the
     Start of D06 with the stress components $\tau_x=\rho\nu_v\,du/dz$, $\tau_y=\rho\nu_v\,dv/dz$ left unexpanded — true for *any*
     relation between stress and shear. · **plain** The Coriolis force on a slab equals the difference of the stresses
     on its top and bottom.
  7. **did** Integrate the first over the layer · **tex** $-\rho f\displaystyle\int_{-\infty}^0v\,dz=\tau_x(0)-\tau_x(-\infty)=\tau$ · **why** The integral of a derivative is
     the difference of its end values (a first integral); the stress is the wind stress on top and zero at depth. ·
     **plain** The whole layer's Coriolis force balances the wind stress.
  8. **did** Integrate the second · **tex** $\rho f\displaystyle\int_{-\infty}^0u\,dz=\tau_y(0)-\tau_y(-\infty)=0$ · **why** There is no stress in the y-direction at
     either end. No eddy-viscosity model was used in steps 6–8. · **plain** The same answer without ever solving for the
     spiral.
- **Result:** $\displaystyle\int_{-\infty}^0u\,dz=0$, $\displaystyle\int_{-\infty}^0v\,dz=-\dfrac{\tau}{\rho f}$ — "the wind-driven transport is τ/(ρf), 90° to the right of the wind in the
  northern hemisphere, 90° to the left in the southern (f < 0 changes its sign), whatever the eddy viscosity".
- **Check:** units — (N/m²)/(kg/m³ · s⁻¹) = m²/s ✓. Number — τ = 0.07 N/m², ρ = 1027 kg/m³, 60° N: 0.540 m²/s; with ν_v four
  times larger: the same. Partial sums — above z = −δ the transport has 86 % of its final size. Code —
  `np.trapezoid` over `GFD.ekman_surface` for two ν_v.
- **What it means:** a wind along a coast moves water offshore or onshore (upwelling, downwelling); a wind with curl makes
  the transport converge or diverge (D09). **Fails when:** the layer feels the bottom (stress not zero at depth), or f → 0.
- **Traps:** direction flips with the sign of f; this is a volume transport per unit width [m²/s] — the mass transport is ρ
  times it; the transport relative to an interior geostrophic flow is unchanged by that flow.

### D08 · The bottom Ekman layer: (13.33)–(13.34) → (13.41), the 45° surface angle and the transport ½Uδ — ★★, 9 steps, in C06 (notebook · `ekman_force_balance`)
- **Goal:** find the wind inside the friction layer under a geostrophic flow, and how much air it carries toward low
  pressure. · **Start:** $-fv=\nu_v\dfrac{d^2u}{dz^2}$ and $fu=\nu_v\dfrac{d^2v}{dz^2}+fU$ (13.33)–(13.34), with $u=U$, $v=0$ as $z\to\infty$ and $u=v=0$ at z = 0 —
  "Coriolis, friction and a pressure gradient written as fU; geostrophic far above, no slip at the ground". · **Plan:** •
  form the complex equation • add a particular solution to the homogeneous one of D06 • apply the two conditions •
  separate components • integrate v. · **Tools:** as D06 (P311, P44); particular solution of a forced linear ODE (P44,
  reminder); $\int_0^\infty e^{-s}\sin s\,ds=\tfrac12$ (P144, reminder). · **Assumptions:** constant ν_v; steady; uniform geostrophic
  interior; flat surface; f > 0 until the Check.
- **Steps:**
  1. **did** Combine the two equations · **tex** $\dfrac{d^2V}{dz^2}=\dfrac{if}{\nu_v}(V-U)$ (13.37) · **why** As in D06: i times the second plus the first; the
     extra term ifU comes from the pressure gradient. · **plain** The same equation as for the ocean layer, now forced by
     the constant U.
  2. **did** Find a particular solution · **tex** $V_p=U$ · **why** A constant has zero second derivative and makes the right side
     zero: it is the geostrophic flow itself, which needs no friction. · **plain** Far from the ground the wind is
     simply U.
  3. **did** Add the homogeneous solution · **tex** $V=A\,e^{-(1+i)z/\delta}+B\,e^{(1+i)z/\delta}+U$ (13.40) · **why** A linear forced equation is solved by a
     particular solution plus the general solution of the unforced one (D06 step 6), with $\delta=\sqrt{2\nu_v/f}$. · **plain**
     Geostrophic wind plus a correction that can decay or grow upward.
  4. **did** Apply the far-field condition · **tex** $B=0$ · **why** As $z\to+\infty$ the term $e^{z/\delta}$ grows; V must tend to U, so its
     coefficient vanishes. Note it is the *other* exponential than in D06, because the fluid is now above z = 0. ·
     **plain** The correction is confined near the ground.
  5. **did** Apply no slip · **tex** $A=-U\ \Rightarrow\ V=U\big[1-e^{-(1+i)z/\delta}\big]$ · **why** At z = 0 the velocity is zero: A + U = 0. · **plain**
     The correction exactly cancels the geostrophic wind at the ground.
  6. **did** Take real and imaginary parts · **tex** $u=U\big[1-e^{-z/\delta}\cos(z/\delta)\big],\qquad v=Ue^{-z/\delta}\sin(z/\delta)$ (13.41) · **why** Write
     $e^{-(1+i)s}=e^{-s}(\cos s-i\sin s)$ with $s=z/\delta$ (Euler's formula). · **plain** A component v across the isobars appears
     — toward low pressure, to the left of U for f > 0.
  7. **did** Look very near the ground · **tex** $V\approx U(1+i)\dfrac z\delta\quad(z\ll\delta)$ · **why** First-order Taylor expansion $e^{-(1+i)s}\approx1-(1+i)s$
     (P26, reminder). · **plain** Close to the surface u = v: the wind blows 45° to the left of the geostrophic wind.
  8. **did** Integrate the cross-isobar component · **tex** $\displaystyle\int_0^\infty v\,dz=U\delta\displaystyle\int_0^\infty e^{-s}\sin s\,ds$ · **why** Change the variable to
     $s=z/\delta$, so $dz=\delta\,ds$; U and δ are constants. · **plain** The inflow toward low pressure, summed over the layer.
  9. **did** Evaluate the integral · **tex** $\displaystyle\int_0^\infty v\,dz=\dfrac12U\delta=U\Big[\dfrac{\nu_v}{2f}\Big]^{1/2}$ · **why** The integral of $e^{-s}\sin s$ over the
     half-line is ½ (integrate by parts twice, or take the imaginary part of $\int e^{-(1-i)s}ds=1/(1-i)$). · **plain**
     Half the geostrophic speed times the layer thickness flows toward low pressure.
- **Result:** $u=U\big[1-e^{-z/\delta}\cos(z/\delta)\big]$, $v=Ue^{-z/\delta}\sin(z/\delta)$ — the two members of (13.41); transport $\tfrac12U\delta$ toward low
  pressure — "friction lets the wind cross the isobars".
- **Check:** units — Uδ [m²/s] ✓. Limits — z → ∞: (U, 0); z = 0: (0, 0). **The overshoot (stated correctly):** du/dz = 0 where
  cos s + sin s = 0, first at s = 3π/4, so the largest u is at z = 3πδ/4 with $u/U=1+e^{-3\pi/4}/\sqrt2=1.0670$; at z = πδ, where v
  first returns to zero, $u/U=1+e^{-\pi}=1.0432$. Number — 60° N, U = 12 m/s, ν_v = 7 m²/s: δ = 332.9 m, transport 1998 m²/s.
  Hemisphere — f < 0: v changes sign (to the right of U), still toward low pressure.
- **What it means:** air converges into lows and rises (weather), diverges from highs and sinks; the same inflow, squashing
  the columns above, spins down every geostrophic vortex. **Fails when:** ν_v varies strongly with height (the real surface
  layer is logarithmic near the ground, Ch. 12) or the surface is sloping.
- **Traps:** A = −U and B = 0 — the dropped exponential is the one that grows *upward* here; "to the left of the
  geostrophic wind" is f > 0 only; the pressure gradient is hidden in the term fU; the maximum of u is at 3πδ/4, not at πδ.

### D09 · Ekman pumping: $w_E=\frac1\rho\,\mathbf e_z\cdot\nabla\times(\boldsymbol\tau/f)$ and $w=\frac\delta2\zeta_g$ — ★★, 7 steps, in C05 (notebook) — ours
- **Goal:** find the vertical velocity that a non-uniform Ekman transport forces at the edge of the layer. **Not in the
  book — ours** (the book gives only the consequence: rising air in a low). · **Start:** the transport of D07 for a stress
  in any direction, $\mathbf M=\dfrac1{\rho f}(\tau_y,\ -\tau_x)$ (the vector form of $\int_{-\infty}^0u\,dz=0$, $\int_{-\infty}^0v\,dz=-\dfrac\tau{\rho f}$ (13.30)), and continuity
  $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$ — "transport at right angles to the stress; mass is conserved". · **Plan:** • integrate continuity
  over the layer • recognise the divergence of the transport • insert the transport • repeat for the bottom layer. ·
  **Tools:** depth integration (P313); divergence and curl of a horizontal vector (Ch. 2, reminder). · **Assumptions:** the
  layer is thin, so the transport formula holds locally (step 1); f varies slowly across the wind pattern (step 5); the
  mean sea surface is level, w(0) = 0 (step 3).
- **Steps:**
  1. **did** Write the transport as a vector · **tex** $\mathbf M=\dfrac1{\rho f}\big(\tau_y,\ -\tau_x\big)$ · **why** D07 for a stress along x gave (0, −τ/ρf); a
     stress along y gives (τ/ρf, 0) by rotating the axes; the problem is linear, so the two add. · **plain** The
     transport is the stress turned 90° to the right (f > 0), divided by ρf.
  2. **did** Integrate continuity over the layer · **tex** $\displaystyle\int_{-D}^0\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)dz+w(0)-w(-D)=0$ · **why** The integral of ∂w/∂z is
     the difference of w between top and bottom of the layer, taken at a depth D of several δ. · **plain** What leaves
     sideways must be replaced from below.
  3. **did** Take the derivatives outside · **tex** $\dfrac{\partial M_x}{\partial x}+\dfrac{\partial M_y}{\partial y}=w(-D)$ · **why** The limits do not depend on x or y, so differentiation
     and integration swap; and w(0) = 0 at the level sea surface. · **plain** The divergence of the transport is the
     upward velocity at the base of the layer.
  4. **did** Name it · **tex** $w_E\equiv w(-D)=\nabla\cdot\mathbf M$ · **why** A definition: the Ekman pumping velocity (positive = upward, "suction"). ·
     **plain** A diverging surface layer sucks water up.
  5. **did** Insert the transport · **tex** $w_E=\dfrac1\rho\Big[\dfrac{\partial}{\partial x}\Big(\dfrac{\tau_y}{f}\Big)-\dfrac{\partial}{\partial y}\Big(\dfrac{\tau_x}{f}\Big)\Big]=\dfrac1\rho\,\mathbf e_z\cdot\nabla\times\Big(\dfrac{\boldsymbol\tau}{f}\Big)$ · **why** Step 1 into step 4;
     the combination is the vertical component of a curl. · **plain** Pumping is set by the curl of the wind stress.
  6. **did** Write the bottom layer's transport · **tex** $\mathbf M'=\displaystyle\int_0^\infty(\mathbf u-\mathbf u_g)\,dz=\dfrac\delta2\big(-U_g-V_g,\ U_g-V_g\big)$ · **why** From D08, $\int(V-V_g)dz=
     -V_g\delta/(1+i)=\tfrac\delta2(-1+i)V_g$ with $V_g=U_g+iV_g$; only this frictional part can diverge (the geostrophic part does
     not, D03 step 6). · **plain** Half a layer thickness of flow is diverted toward low pressure and against the wind.
  7. **did** Take its divergence · **tex** $w(\text{top})=-\nabla\cdot\mathbf M'=\dfrac\delta2\Big(\dfrac{\partial V_g}{\partial x}-\dfrac{\partial U_g}{\partial y}\Big)=\dfrac\delta2\,\zeta_g$ · **why** Continuity over the bottom layer
     with w = 0 at the ground; the terms $\partial U_g/\partial x+\partial V_g/\partial y$ vanish because the geostrophic flow is non-divergent. ·
     **plain** A cyclone pumps air up out of the friction layer at half a layer thickness times its vorticity.
- **Result:** $w_E=\dfrac1\rho\,\mathbf e_z\cdot\nabla\times\Big(\dfrac{\boldsymbol\tau}{f}\Big)$ (surface layer) and $w=\dfrac\delta2\zeta_g$ (bottom layer, f > 0; in general
  $w=\mathrm{sgn}(f)\dfrac\delta2\zeta_g$ with $\delta=\sqrt{2\nu_v/\lvert f\rvert}$) — "cyclonic forcing pumps upward, in both hemispheres".
- **Check:** units — (N/m²)/(kg/m³ · s⁻¹ · m) = m/s ✓; δζ = m/s ✓. Number — a stress changing by 0.07 N/m² over 1000 km at 35° N:
  8.1 × 10⁻⁷ m/s ≈ 26 m per year; bottom layer with δ = 333 m, ζ_g = 10⁻⁵ s⁻¹: 1.66 mm/s. Hemisphere — a southern low has ζ_g <
  0 and f < 0: w > 0 again.
- **What it means:** the winds of a subtropical gyre (easterlies equatorward, westerlies poleward) push surface water
  toward its middle and pump it down; the interior responds by moving equatorward (Sverdrup balance, one line in C05).
  **Fails when:** f → 0 (the equator: upwelling there needs its own argument), or the wind varies on the scale of δ.
- **Traps:** the sign — cyclonic stress curl gives upwelling for f > 0; keeping df/dy inside the curl is the first step
  toward Sverdrup balance and is not pursued; surface and bottom layers pump with opposite roles (one drives the
  interior, the other is driven by it).

### D10 · Shallow-water continuity (13.44) and the linear set (13.45) — ★★, 9 steps, in C07 (notebook)
- **Goal:** collapse the three-dimensional equations of a thin homogeneous layer to three equations in x, y, t. ·
  **Start:** $p=\rho g(H+\eta-z)$ (hydrostatic, zero pressure at the surface), $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$, and the inviscid horizontal
  momentum equations $\dfrac{Du}{Dt}-fv=-\dfrac1\rho\dfrac{\partial p}{\partial x}$, $\dfrac{Dv}{Dt}+fu=-\dfrac1\rho\dfrac{\partial p}{\partial y}$ (the horizontal members of (13.9) without friction) —
  "weight sets the pressure; mass is conserved; Coriolis and pressure accelerate the columns". · **Plan:** • get the
  pressure gradient from the surface slope • integrate continuity over the column • use what the surface and the bottom
  do • write it as a flux • linearise. · **Tools:** depth integration (P313); the kinematic surface condition (Ch. 7,
  reminder); product rule (P38, reminder); linearisation. · **Assumptions:** wavelength ≫ depth, so hydrostatic (the Start);
  homogeneous, inviscid; flat bottom at z = 0 (step 4); small amplitude (steps 7–9).
- **Steps:**
  1. **did** Differentiate the pressure horizontally · **tex** $\dfrac{\partial p}{\partial x}=\rho g\dfrac{\partial\eta}{\partial x},\qquad\dfrac{\partial p}{\partial y}=\rho g\dfrac{\partial\eta}{\partial y}$ (13.42) · **why** In
     $p=\rho g(H+\eta-z)$ only η depends on x and y; H is constant and z is the independent coordinate. · **plain** The
     horizontal push is the slope of the surface, the same at every depth.
  2. **did** Conclude that the columns stay upright · **tex** $\dfrac{\partial u}{\partial z}=\dfrac{\partial v}{\partial z}=0$ · **why** The force of step 1 and the Coriolis force act
     alike at all depths; a flow that starts independent of z has no reason to become dependent. · **plain** Each
     column moves as a whole.
  3. **did** Integrate continuity over the column · **tex** $(H+\eta)\dfrac{\partial u}{\partial x}+(H+\eta)\dfrac{\partial v}{\partial y}+w(\eta)-w(0)=0$ (13.43) · **why** By step 2 the
     horizontal divergence is constant over the depth H + η and leaves the integral; the integral of ∂w/∂z is w at the
     top minus w at the bottom. · **plain** Horizontal convergence must be taken up by the surface rising.
  4. **did** State what the two ends do · **tex** $w(0)=0,\qquad w(\eta)=\dfrac{D\eta}{Dt}=\dfrac{\partial\eta}{\partial t}+u\dfrac{\partial\eta}{\partial x}+v\dfrac{\partial\eta}{\partial y}$ · **why** No flow through the flat
     bottom; a parcel on the surface stays on it — the exact kinematic condition, of which Chapter 7 used the linearised
     form. · **plain** The surface moves with the fluid.
  5. **did** Substitute · **tex** $(H+\eta)\dfrac{\partial u}{\partial x}+(H+\eta)\dfrac{\partial v}{\partial y}+\dfrac{\partial\eta}{\partial t}+u\dfrac{\partial\eta}{\partial x}+v\dfrac{\partial\eta}{\partial y}=0$ · **why** Insert step 4 into step 3. · **plain** One
     equation linking the surface height and the horizontal flow.
  6. **did** Recognise two product rules · **tex** $\dfrac{\partial\eta}{\partial t}+\dfrac{\partial}{\partial x}\big[u(H+\eta)\big]+\dfrac{\partial}{\partial y}\big[v(H+\eta)\big]=0$ (13.44) · **why** Since H is constant,
     $\partial[u(H+\eta)]/\partial x=(H+\eta)\,\partial u/\partial x+u\,\partial\eta/\partial x$, and likewise in y. · **plain** The surface falls where the transport
     diverges.
  7. **did** Linearise continuity · **tex** $\dfrac{\partial\eta}{\partial t}+H\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$ · **why** For small waves η ≪ H, so H + η ≈ H, and
     products of two small quantities (uη, vη) are neglected. · **plain** For small waves the depth in the transport is
     just H.
  8. **did** Insert the pressure gradient in the x-momentum equation and linearise · **tex** $\dfrac{\partial u}{\partial t}-fv=-g\dfrac{\partial\eta}{\partial x}$ · **why** Step 1 turns
     the pressure term into −g ∂η/∂x; the advective part of Du/Dt is a product of small quantities (U ≪ c) and is
     dropped. · **plain** A surface slope accelerates the column; Coriolis turns it.
  9. **did** Do the same in y · **tex** $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$ · **why** The identical two moves on the y-equation; steps 7–9 together
     are the three members of (13.45), $\frac{\partial\eta}{\partial t}+H(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y})=0$, $\frac{\partial u}{\partial t}-fv=-g\frac{\partial\eta}{\partial x}$, $\frac{\partial v}{\partial t}+fu=-g\frac{\partial\eta}{\partial y}$. · **plain** Three linear
     equations for η, u, v.
- **Result:** $\dfrac{\partial\eta}{\partial t}+H\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$, $\dfrac{\partial u}{\partial t}-fv=-g\dfrac{\partial\eta}{\partial x}$, $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.45) — "the linear shallow-water
  equations".
- **Check:** units — H × (1/s) = m/s = ∂η/∂t ✓; g × slope = m/s² ✓. Limit — f = 0, one dimension: $\eta_{tt}=gH\,\eta_{xx}$ (subscripts =
  derivatives), waves at $c=\sqrt{gH}$, Chapter 7's long-wave speed. Conservation — the flux form of step 6 keeps the total
  volume exactly. Number — H = 4200 m: c = 203 m/s.
- **What it means:** all of §13.10–§13.15 are solutions of this set; with H replaced by an equivalent depth it also
  governs each vertical mode of a stratified fluid (D11). **Fails when:** the wavelength is not long against the depth
  (non-hydrostatic), or the amplitude is not small (then keep the flux form and the advection: C13).
- **Traps:** why u, v are independent of z (the pressure gradient is); the flux form needs the product rule; what
  "linear" drops (η/H ≪ 1 and U/c ≪ 1); z = 0 is the bottom in this section (trap T4).


### D11 · Vertical normal modes: (13.52) → (13.56) → the modal equations and (13.62); orthogonality — ★★★, 13 steps, in C08 (notebook · `vertical_modes`)
- **Goal:** show that a continuously stratified, hydrostatic layer splits into independent "modes", each with a fixed
  vertical shape and each obeying the shallow-water equations with its own speed. · **Start:** the five linear equations
  $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$ (13.47), $\dfrac{\partial u}{\partial t}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$ and $\dfrac{\partial v}{\partial t}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$ (13.48)–(13.49), $0=-\dfrac{\partial p}{\partial z}-g\rho$ (13.50),
  $\dfrac{\partial\rho}{\partial t}-\dfrac{\rho_0N^2}{g}w=0$ (13.51) — "mass, two momentum balances, hydrostatics, and density carried up and down". · **Plan:** • give
  u, v and p one common vertical shape • let continuity and hydrostatics fix the shapes of w and ρ • separate the
  density equation into a z-part and an (x, y, t)-part • read off the vertical eigenproblem and the horizontal equations
  • prove orthogonality (weight 1, for a rigid lid and for a free surface). · **Tools:** separation of variables (P167, reminder); Sturm–Liouville problems (P314);
  integration by parts (P218a, reminder). · **Assumptions:** linear; hydrostatic, so ω ≪ N; flat bottom at z = −H; no mean
  shear; N = N(z) only.
- **Steps:**
  1. **did** Give u, v and p one vertical shape · **tex** $[u,v,p/\rho_0]=\displaystyle\sum_{n=0}^\infty[u_n,v_n,p_n]\,\psi_n(z)$ (13.52) · **why** The two momentum equations contain
     no z-derivative, so if p has the shape ψ_n(z), so must u and v; the amplitudes depend on x, y, t only. · **plain**
     In each mode the horizontal flow and the pressure rise and fall together with depth.
  2. **did** Find the shape of w from continuity · **tex** $w=\displaystyle\sum_{n=0}^\infty w_n\displaystyle\int_{-H}^z\psi_n\,dz,\qquad w_n=-\Big(\dfrac{\partial u_n}{\partial x}+\dfrac{\partial v_n}{\partial y}\Big)$ (13.53) · **why** Continuity gives
     $\partial w/\partial z=-(\partial u/\partial x+\partial v/\partial y)$; integrate upward from the flat bottom, where w = 0. · **plain** The vertical velocity's
     shape is the running integral of the horizontal one's.
  3. **did** Find the shape of ρ from hydrostatics · **tex** $\rho=\displaystyle\sum_{n=0}^\infty\rho_n\dfrac{d\psi_n}{dz},\qquad\rho_n=-\dfrac{\rho_0}{g}p_n$ (13.54) · **why** Hydrostatic balance gives
     $\rho=-\frac1g\,\partial p/\partial z$; differentiate step 1 in z. · **plain** The density's shape is the slope of the pressure's
     shape.
  4. **did** Substitute into the density equation · **tex** $\displaystyle\sum_{n=0}^\infty\Big[\dfrac{\partial\rho_n}{\partial t}\dfrac{d\psi_n}{dz}-\dfrac{\rho_0N^2}{g}w_n\displaystyle\int_{-H}^z\psi_n\,dz\Big]=0$ · **why** Insert steps 2 and 3 into the
     fifth Start equation; it is the only one not yet used. · **plain** One condition is left to tie everything
     together.
  5. **did** Separate the variables · **tex** $\dfrac{d\psi_n/dz}{N^2\int_{-H}^z\psi_n\,dz}=\dfrac{\rho_0}{g}\dfrac{w_n}{\partial\rho_n/\partial t}\equiv-\dfrac1{c_n^2}$ (13.55) · **why** For independent modes each
     bracket vanishes; dividing puts everything that depends on z on the left and on (x, y, t) on the right, so both
     equal a constant. · **plain** The constant, named −1/c_n², will turn out to be a speed.
  6. **did** Differentiate the z-part · **tex** $\dfrac{d}{dz}\Big(\dfrac1{N^2}\dfrac{d\psi_n}{dz}\Big)+\dfrac1{c_n^2}\psi_n=0$ (13.56) · **why** The z-side of step 5 reads
     $\frac1{N^2}\psi_n'=-\frac1{c_n^2}\int_{-H}^z\psi_n\,dz$; one z-derivative removes the integral. · **plain** An eigenvalue problem: only special
     values of c_n allow a shape that fits top and bottom.
  7. **did** Read the (x, y, t)-part · **tex** $w_n=-\dfrac{g}{\rho_0c_n^2}\dfrac{\partial\rho_n}{\partial t}=\dfrac1{c_n^2}\dfrac{\partial p_n}{\partial t}$ (13.61) · **why** The other side of step 5, then
     $\rho_n=-\rho_0p_n/g$ from step 3. · **plain** The mode's vertical velocity follows the rate of change of its pressure.
  8. **did** Combine with step 2 · **tex** $\dfrac{\partial u_n}{\partial x}+\dfrac{\partial v_n}{\partial y}+\dfrac1{c_n^2}\dfrac{\partial p_n}{\partial t}=0$ (13.57) · **why** Step 2 says $w_n=-(\partial u_n/\partial x+\partial v_n/\partial y)$; equate
     with step 7. · **plain** A continuity equation for the mode, with p_n in the role of surface height.
  9. **did** Substitute the shapes into the momentum equations · **tex** $\dfrac{\partial u_n}{\partial t}-fv_n=-\dfrac{\partial p_n}{\partial x},\qquad\dfrac{\partial v_n}{\partial t}+fu_n=-\dfrac{\partial p_n}{\partial y}$ (13.58)–(13.59) ·
     **why** Every term carries the same ψ_n(z); for independent modes the coefficient of each ψ_n must vanish. ·
     **plain** Each mode has its own pair of momentum equations.
  10. **did** Compare with shallow water · **tex** $p_n\leftrightarrow g\eta,\qquad c_n^2\leftrightarrow gH,\qquad c_n^2\equiv gH_e$ (13.62) · **why** With these two replacements
      steps 8 and 9 are exactly the three members of the shallow-water set of D10's Result. · **plain** Each mode is a
      shallow-water layer of "equivalent depth" H_e.
  11. **did** Write the structure equation for two modes and cross-multiply · **tex**
      $\psi_n\Big(\dfrac{\psi_m'}{N^2}\Big)'-\psi_m\Big(\dfrac{\psi_n'}{N^2}\Big)'+\Big(\dfrac1{c_m^2}-\dfrac1{c_n^2}\Big)\psi_m\psi_n=0$ · **why** Multiply the equation of mode m by ψ_n and that of mode n by
      ψ_m, then subtract (a prime is d/dz); the aim is to isolate the product ψ_mψ_n. · **plain** A relation between any
      two modes.
  12. **did** Integrate over the depth by parts · **tex** $\Big[\dfrac{\psi_n\psi_m'-\psi_m\psi_n'}{N^2}\Big]_{-H}^0+\Big(\dfrac1{c_m^2}-\dfrac1{c_n^2}\Big)\displaystyle\int_{-H}^0\psi_m\psi_n\,dz=0$ · **why** The first two terms
      are the derivative of the bracket (the cross terms $\psi_n'\psi_m'/N^2$ cancel), so their integral is its end values. ·
      **plain** Everything reduces to values at the top and the bottom.
  13. **did** Use the boundary conditions · **tex** $\displaystyle\int_{-H}^0\psi_m\psi_n\,dz=0\qquad(m\neq n,\ \text{either lid})$ · **why** At the flat bottom w = 0 gives ψ′ = 0 (note N66 above). At the top either
      ψ′ = 0 (lid) or ψ′ = −(N²/g)ψ for both modes (the free-surface condition of note N67 above), so the numerator cancels; with c_m ≠ c_n the integral
      vanishes. · **plain** Different modes are orthogonal.
- **Result:** $\dfrac{d}{dz}\Big(\dfrac1{N^2}\dfrac{d\psi_n}{dz}\Big)+\dfrac1{c_n^2}\psi_n=0$ with $c_n^2\equiv gH_e$, and for each n the three shallow-water equations
  $\dfrac{\partial u_n}{\partial x}+\dfrac{\partial v_n}{\partial y}+\dfrac1{c_n^2}\dfrac{\partial p_n}{\partial t}=0$, $\dfrac{\partial u_n}{\partial t}-fv_n=-\dfrac{\partial p_n}{\partial x}$, $\dfrac{\partial v_n}{\partial t}+fu_n=-\dfrac{\partial p_n}{\partial y}$ — "a stratified layer is a stack of shallow-water
  systems".
- **Check:** units — ψ_n dimensionless, so $(1/N^2)\psi''$ is s²/m² and $\psi/c_n^2$ is s²/m² ✓. Uniform N with a rigid lid — $\psi_n=\cos(n\pi z/H)$,
  $c_n=NH/(n\pi)$: 3.61 m/s for H = 4200 m, N = 2.7 × 10⁻³ s⁻¹. Free surface — the bracket of step 12 vanishes at z = 0 too:
  $\psi_n\psi_m'-\psi_m\psi_n'=-\frac{N^2}g(\psi_n\psi_m-\psi_m\psi_n)=0$ by $\dfrac{d\psi_n}{dz}+\dfrac{N^2}{g}\psi_n=0$ at $z=0$ (13.65), so the weight is 1 for both lids
  (`VM.orthogonality_matrix(modes)`, default `kind="psi"`; verified by the designer's own finite-volume solve of the
  thermocline profile with a free surface: off-diagonal overlaps below 10⁻⁸). The surface term appears only in the second
  relation, obtained by multiplying the structure equation by ψ_m and integrating once by parts:
  $\dfrac1{c_n^2}\displaystyle\int_{-H}^0\psi_n\psi_m\,dz=\displaystyle\int_{-H}^0\dfrac{\psi_n'\psi_m'}{N^2}\,dz+\dfrac{\psi_n(0)\psi_m(0)}{g}$ (`kind="energy"`).
- **sympy check** (`check_src`; executed as written — every assertion passes and the final line is printed):
  ```python
  import sympy as sp                                            # symbolic algebra
  z, t = sp.symbols("z t", real=True)                           # height (z = -H bottom, 0 top) and time
  N, H, g, rho0 = sp.symbols("N H g rho_0", positive=True)      # uniform N for the check
  n, m = sp.symbols("n m", integer=True, positive=True)         # two mode numbers
  c = N * H / (n * sp.pi)                                       # the rigid-lid speed we expect for mode n
  psi = sp.cos(n * sp.pi * z / H)                               # its vertical structure
  ode = sp.diff(sp.diff(psi, z) / N**2, z) + psi / c**2         # step 6: the structure equation
  assert sp.simplify(ode) == 0                                  # ... is satisfied
  assert sp.diff(psi, z).subs(z, -H) == 0 and sp.diff(psi, z).subs(z, 0) == 0   # w = 0 at bottom and lid
  W = sp.integrate(psi, (z, -H, z))                             # step 2: the shape of w is the integral of psi
  ratio = sp.simplify(sp.diff(psi, z) / (N**2 * W))             # step 5: the z-side of the separation
  assert sp.simplify(ratio + 1 / c**2) == 0                     # ... equals the constant -1/c_n^2
  pn = sp.Function("p_n")(t)                                    # modal pressure amplitude p_n(t) (x, y left out)
  rho_n = -rho0 / g * pn                                        # step 3: hydrostatics, rho_n = -(rho0/g) p_n
  w_n = sp.diff(pn, t) / c**2                                   # step 7: w_n = (1/c_n^2) dp_n/dt
  dens = sp.diff(rho_n, t) * sp.diff(psi, z) - rho0 * N**2 / g * w_n * W   # step 4: the density equation
  assert sp.simplify(dens) == 0                                 # every mode satisfies it on its own
  psi_m = sp.cos(m * sp.pi * z / H)                             # a second mode
  overlap = sp.integrate(psi * psi_m, (z, -H, 0))               # step 13: the overlap integral
  assert sp.simplify(overlap.subs({n: 1, m: 2})) == 0 and sp.simplify(overlap.subs({n: 2, m: 3})) == 0
  print("D11 checks passed: structure equation, separation constant, density equation, orthogonality")
  ```
- **What it means:** every wave result for one shallow layer (C09–C15) holds for each mode with its own c_n — the first
  baroclinic mode's c₁ sets the internal Rossby radius and the speed of thermocline signals. **Fails when:** the bottom is
  not flat or there is a sheared mean current (the modes couple); the motion is not hydrostatic (ω comparable to N).
- **Traps:** w's shape is the integral of ψ_n and ρ's is its derivative, not the reverse; the separation constant is
  negative; the modal amplitudes have different units (trap T12); the surface term $\psi_m(0)\psi_n(0)/g$ belongs to the
  energy relation, not to the weight-1 orthogonality, which holds for both lids.

### D12 · One equation for v, (13.75), from the linear shallow-water set on a β-plane — ★★★, 10 steps, in C09 (notebook · `shallow_water_dispersion`)
- **Goal:** eliminate u and η from the three shallow-water equations so that a single equation for v remains, valid at
  all frequencies. · **Start:** $\dfrac{\partial\eta}{\partial t}+H\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$, $\dfrac{\partial u}{\partial t}-fv=-g\dfrac{\partial\eta}{\partial x}$, $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.45) with
  $f=f_0+\beta y$ (13.10) — "the linear shallow-water set, now with f varying northward". · **Plan:** • remove η by one
  time-derivative of each momentum equation • remove ∂²u/∂t² between the two • build a vorticity equation to remove what
  is left of u • add. · **Tools:** operator elimination for linear PDEs (P177, reminder); cross-differentiation;
  Schwarz's theorem (P121). · **Assumptions:** linear; β-plane: f is treated as the constant f₀ except where it is
  differentiated (steps 6, 8; trap T10).
- **Steps:**
  1. **did** Differentiate the x-momentum equation in time · **tex** $\dfrac{\partial^2u}{\partial t^2}-f\dfrac{\partial v}{\partial t}=-g\dfrac{\partial^2\eta}{\partial x\,\partial t}$ · **why** The parameter f does not depend on t; the
     aim is to make ∂η/∂t appear, because continuity can replace it. · **plain** How the acceleration of u changes.
  2. **did** Replace ∂η/∂t by continuity · **tex** $\dfrac{\partial^2u}{\partial t^2}-f\dfrac{\partial v}{\partial t}=gH\dfrac{\partial}{\partial x}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ (13.72) · **why** Continuity gives
     $\partial\eta/\partial t=-H(\partial u/\partial x+\partial v/\partial y)$; the two minus signs cancel. · **plain** The surface height η is gone from the x-equation.
  3. **did** Do the same to the y-momentum equation · **tex** $\dfrac{\partial^2v}{\partial t^2}+f\dfrac{\partial u}{\partial t}=gH\dfrac{\partial}{\partial y}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ (13.73) · **why** The identical
     two moves; H is constant and passes through ∂/∂y. · **plain** Now η is gone from both.
  4. **did** Differentiate step 3 in time and remove ∂²u/∂t² · **tex**
     $\dfrac{\partial^3v}{\partial t^3}+f\Big[f\dfrac{\partial v}{\partial t}+gH\dfrac{\partial}{\partial x}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)\Big]=gH\dfrac{\partial^2}{\partial y\,\partial t}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)$ (13.74) · **why** The time derivative of step 3 contains $f\,\partial^2u/\partial t^2$;
     step 2 gives that quantity as the square bracket. · **plain** Now u appears only inside the divergence.
  5. **did** Cross-differentiate the original momentum equations · **tex** $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial u}{\partial y}-\dfrac{\partial v}{\partial x}\Big)-\dfrac{\partial(fv)}{\partial y}-\dfrac{\partial(fu)}{\partial x}=0$ · **why** ∂/∂y of the
     x-equation minus ∂/∂x of the y-equation; the pressure terms $g\,\partial^2\eta/\partial x\partial y$ cancel (Schwarz). · **plain** A
     vorticity equation with no η in it.
  6. **did** Differentiate the products · **tex** $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial u}{\partial y}-\dfrac{\partial v}{\partial x}\Big)-f_0\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)-\beta v=0$ · **why** Product rule: $\partial(fv)/\partial y=f\,\partial v/\partial y+\beta v$
     because df/dy = β, while f does not depend on x; then f → f₀ in the undifferentiated factor. · **plain** Vorticity
     changes through divergence and through northward motion.
  7. **did** Take ∂/∂x and multiply by gH · **tex** $gH\dfrac{\partial^2}{\partial t\,\partial x}\Big(\dfrac{\partial u}{\partial y}-\dfrac{\partial v}{\partial x}\Big)-f_0\,gH\dfrac{\partial}{\partial x}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)-gH\beta\dfrac{\partial v}{\partial x}=0$ · **why**
     Chosen so that its middle term is minus the term $f\,gH\,\partial(\ldots)/\partial x$ of step 4. · **plain** The vorticity equation,
     dressed to match step 4.
  8. **did** Add step 7 to step 4 (with f → f₀) · **tex**
     $\dfrac{\partial^3v}{\partial t^3}+f_0^2\dfrac{\partial v}{\partial t}+gH\dfrac{\partial}{\partial t}\Big[\dfrac{\partial^2u}{\partial x\,\partial y}-\dfrac{\partial^2v}{\partial x^2}-\dfrac{\partial^2u}{\partial y\,\partial x}-\dfrac{\partial^2v}{\partial y^2}\Big]-gH\beta\dfrac{\partial v}{\partial x}=0$ · **why** The two $f_0\,gH\,\partial(\ldots)/\partial x$ terms cancel; the
     remaining divergence terms are collected under one ∂/∂t. · **plain** Only time and space derivatives of u and v
     remain.
  9. **did** Cancel the mixed derivatives of u · **tex** $\dfrac{\partial^2u}{\partial x\,\partial y}-\dfrac{\partial^2u}{\partial y\,\partial x}=0\ \Rightarrow\ [\ \cdot\ ]=-\Big(\dfrac{\partial^2v}{\partial x^2}+\dfrac{\partial^2v}{\partial y^2}\Big)=-\nabla_H^2v$ · **why**
     Schwarz's theorem again; this is the step that finally removes u. · **plain** The velocity u has dropped out completely.
  10. **did** Write the result · **tex** $\dfrac{\partial^3v}{\partial t^3}-gH\dfrac{\partial}{\partial t}\nabla_H^2v+f_0^2\dfrac{\partial v}{\partial t}-gH\beta\dfrac{\partial v}{\partial x}=0$ (13.75) · **why** Insert step 9 into step 8 and
      order the terms as the book does. · **plain** One equation for the northward velocity: gravity, rotation and β
      each contribute one term.
- **Result:** $\dfrac{\partial^3v}{\partial t^3}-gH\dfrac{\partial}{\partial t}\nabla_H^2v+f_0^2\dfrac{\partial v}{\partial t}-gH\beta\dfrac{\partial v}{\partial x}=0$, with $\nabla_H^2=\partial^2/\partial x^2+\partial^2/\partial y^2$ — "every wave of the rotating
  shallow-water layer satisfies this".
- **Check:** units — every term is (m/s)/s³ ✓ (gHβ ∂v/∂x: m²/s² · 1/(m s) · 1/s). Limits — f₀ = β = 0: $v_{tt}=gH\nabla_H^2v$ after
  one time integration, the wave equation of Chapter 7; β = 0 and no time dependence: any steady v is allowed
  (geostrophic flow).
- **sympy check** (`check_src`; executed as written — every assertion passes and the final line is printed):
  ```python
  import sympy as sp                                            # symbolic algebra
  s, a, b = sp.symbols("s a b")                                 # stand-ins for d/dt, d/dx, d/dy (they commute)
  u, v = sp.symbols("u v")                                      # the two velocity components (any plane wave)
  f0, beta, gH = sp.symbols("f_0 beta gH", positive=True)       # f0, beta and the product g*H
  div = a * u + b * v                                           # horizontal divergence du/dx + dv/dy
  e72 = s**2 * u - f0 * s * v - gH * a * div                    # step 2: (13.72) with everything on the left
  e73 = s**2 * v + f0 * s * u - gH * b * div                    # step 3: (13.73)
  e74 = sp.expand(s * e73 - f0 * e72)                           # step 4: d/dt of (13.73), then remove d2u/dt2
  book74 = s**3 * v + f0 * (f0 * s * v + gH * a * div) - gH * b * s * div   # (13.74) as printed
  assert sp.expand(e74 - book74) == 0                           # our step 4 is the book's (13.74)
  vort = s * (b * u - a * v) - f0 * div - beta * v              # step 6: the linear vorticity equation
  single = sp.expand(e74 + gH * a * vort)                       # steps 7-8: add gH d/dx of the vorticity equation
  target = s**3 * v - gH * s * (a**2 + b**2) * v + f0**2 * s * v - gH * beta * a * v   # step 10: (13.75)
  assert sp.expand(single - target) == 0                        # step 9: u has dropped out, one equation for v
  k, l, w = sp.symbols("k l omega", real=True)                  # plane wave exp(i(kx + ly - wt))
  disp = sp.expand((target / v).subs({s: -sp.I * w, a: sp.I * k, b: sp.I * l}) / sp.I)   # D13 steps 1-2
  assert sp.expand(disp - (w**3 - gH * w * (k**2 + l**2) - f0**2 * w - gH * beta * k)) == 0   # (13.76)
  print("D12 checks passed: (13.74), the elimination of u, (13.75) and the cubic (13.76)")
  ```
- **What it means:** the three time derivatives are why there are three waves per wavenumber (D13). **Fails when:** βy is
  not small against f₀ (large north–south extent, or near the equator), or the waves are not small.
- **Traps:** f is replaced by f₀ everywhere except where it is differentiated; the two $f_0\,gH\,\partial(\ldots)/\partial x$ terms cancel
  only because the vorticity equation is brought in; keep track of the order of the time derivatives.

### D13 · The cubic (13.76): three real roots, which is which, three regimes — ★★, 7 steps, in C09 (notebook · `shallow_water_dispersion`)
- **Goal:** turn the single equation for v into a relation between frequency and wavenumber, show that it has three real
  solutions, and learn which terms matter for each. · **Start:** $\dfrac{\partial^3v}{\partial t^3}-gH\dfrac{\partial}{\partial t}\nabla_H^2v+f_0^2\dfrac{\partial v}{\partial t}-gH\beta\dfrac{\partial v}{\partial x}=0$ (13.75) and the trial
  wave $v=\hat v\,e^{i(kx+ly-\omega t)}$ — "the Result of D12 and a plane wave". · **Plan:** • substitute the wave • recognise a cubic with
  no ω² term • test its discriminant • find the fast pair and the slow root by dominant balance. · **Tools:** plane-wave
  substitution (P176, reminder); the discriminant and trigonometric form of a cubic (P319); dominant balance (P198,
  reminder). · **Assumptions:** β-plane scaling for step 5 ($\beta c<f_0^2$).
- **Steps:**
  1. **did** Substitute the plane wave · **tex** $(-i\omega)^3\hat v-gH(-i\omega)\big(-k^2-l^2\big)\hat v+f_0^2(-i\omega)\hat v-gH\beta(ik)\hat v=0$ · **why** For this wave each
     ∂/∂t is a factor −iω, each ∂/∂x a factor ik, each ∂/∂y a factor il; the exponential cancels. · **plain** A
     differential equation has become algebra.
  2. **did** Divide by $i\hat v$ · **tex** $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76), with $K^2=k^2+l^2$, $c=\sqrt{gH}$ · **why** Because $(-i)^3=i$, every term carries
     one factor i; divide it out, with $\hat v\neq0$. · **plain** A cubic in ω: three waves for every wavenumber.
  3. **did** Identify the coefficients · **tex** $\omega^3+p\,\omega+q=0,\qquad p=-(c^2K^2+f_0^2),\qquad q=-c^2\beta k$ · **why** There is no ω² term, so the
     cubic is already in "depressed" form and the three roots sum to zero. · **plain** The coefficient p collects gravity and rotation,
     and q is the β term.
  4. **did** Form the discriminant · **tex** $\Delta=-4p^3-27q^2=4\big(c^2K^2+f_0^2\big)^3-27\big(c^2\beta k\big)^2$ · **why** A real cubic has three distinct real
     roots exactly when Δ > 0 (P319). The book asserts this without proof. · **plain** The sign of one number decides
     whether all three waves are real.
  5. **did** Bound it from below · **tex** $\big(c^2k^2+\tfrac12f_0^2+\tfrac12f_0^2\big)^3\ge\tfrac{27}{4}c^2k^2f_0^4\ \Rightarrow\ \Delta\ge27\,c^2k^2\big(f_0^4-c^2\beta^2\big)$ · **why** The mean of three
     positive numbers is at least their geometric mean, and $K^2\ge k^2$. So Δ > 0 whenever $\beta c<f_0^2$: the Rossby
     radius c/f₀ is smaller than f₀/β = R tan θ₀. · **plain** All three roots are real for every
     wavenumber as long as the β-plane itself makes sense.
  6. **did** Find the fast roots by dropping the β term · **tex** $\omega\big(\omega^2-c^2K^2-f_0^2\big)=0\ \Rightarrow\ \omega=\pm\sqrt{f_0^2+c^2K^2}$ · **why** For
     $\lvert\omega\rvert\ge\lvert f_0\rvert$ the ratio of the β term to the second term is of order β/(ωK) — about 0.02 for our 3100 km wave at 35° N (the code cell below prints the terms), smaller for shorter waves; neglect it. · **plain** Two
     gravity waves bent by rotation, one in each direction, never slower than f₀.
  7. **did** Find the slow root by dropping ω³ · **tex** $\omega\simeq-\dfrac{c^2\beta k}{c^2K^2+f_0^2}=-\dfrac{\beta k}{K^2+f_0^2/c^2}$ · **why** For $\lvert\omega\rvert\ll\lvert f_0\rvert$, ω³ is smaller than
     $f_0^2\omega$ by $(\omega/f_0)^2$; the remaining three terms balance. · **plain** One slow wave that exists only because of β and
     moves its crests westward.
- **Result:** $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ has three real roots when $\beta c<f_0^2$: $\omega\simeq\pm\sqrt{f_0^2+c^2K^2}$ (fast) and
  $\omega\simeq-\beta k/(K^2+f_0^2/c^2)$ (slow) — "two Poincaré waves and one Rossby wave".
- **Check:** units — every term s⁻³ ✓. Number — 35° N, c = 202.95 m/s, wavelength 3100 km: (−4.1525 × 10⁻⁴, −8.888 × 10⁻⁶, +4.2414
  × 10⁻⁴) s⁻¹; sum zero; the approximations of steps 6 and 7 are off by 1 % and 0.04 %. The bound of step 5 is sharp: Δ
  first reaches zero at $c^2k^2=f_0^2/2$ when βc = f₀² (ours, computed: for c = 203 m/s that is at latitude 26.3°; closer to
  the equator Δ < 0 for a band of very long waves — the β-plane stretched too far). Code — `GFD.shallow_water_discriminant`
  is the flag; where it is negative `GFD.shallow_water_omega` returns (nan, nan, nan), e.g. at 12° N with the external
  speed for a 20 000 km wave.
- **What it means:** a model that wants only weather may drop ω³ (filtering the fast waves); a tide model may drop β.
  **Fails when:** βc is not small against f₀² (external mode at low latitude), where the fixed-f₀ replacement is itself
  invalid.
- **Traps:** the sign of the β term (with the wrong sign the slow wave would go east); the fast roots have opposite signs
  and satisfy only abs(ω) > abs(f₀) (trap T9); slip #7 — the book prints ω ≫ f for the range where ω³ is negligible, the
  correct form is ω ≪ f; the slow root is a thousandth of the fast ones — compare roots with relative tolerances.

### D14 · Poincaré waves: (13.77)–(13.79) → (13.80) → (13.82), and the group velocity — ★★, 8 steps, in C10 (notebook)
- **Goal:** find the frequency of a long gravity wave on a rotating plane and how the current is oriented under it. ·
  **Start:** $\dfrac{\partial\eta}{\partial t}+H\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$, $\dfrac{\partial u}{\partial t}-fv=-g\dfrac{\partial\eta}{\partial x}$, $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.45) with constant f, and
  $(u,v,\eta)=(\hat u,\hat v,\hat\eta)\,e^{i(kx+ly-\omega t)}$ — "the shallow-water set on an f-plane and a plane wave". · **Plan:** • substitute • solve
  the two momentum equations for the velocity amplitudes • put them into continuity • differentiate the result for the
  group velocity. · **Tools:** complex amplitudes (P176, reminder); polarisation relations (P321); chain rule (P49,
  reminder). · **Assumptions:** f-plane; ω ≠ ±f in step 2.
- **Steps:**
  1. **did** Substitute the plane wave · **tex** $-i\omega\hat u-f\hat v=-ikg\hat\eta,\quad-i\omega\hat v+f\hat u=-ilg\hat\eta,\quad-i\omega\hat\eta+iH(k\hat u+l\hat v)=0$ (13.77)–(13.79) ·
     **why** Each derivative becomes a factor (−iω, ik or il) and the common exponential cancels. · **plain** Three
     algebraic equations for three amplitudes.
  2. **did** Solve the first two for the velocities · **tex** $\hat u=\dfrac{(-ikg\hat\eta)(-i\omega)-(-f)(-ilg\hat\eta)}{(-i\omega)^2+f^2},\qquad\hat v=\dfrac{(-i\omega)(-ilg\hat\eta)-f(-ikg\hat\eta)}{(-i\omega)^2+f^2}$ · **why** Cramer's
     rule for a 2 × 2 system; the determinant is $(-i\omega)(-i\omega)-(-f)(f)=f^2-\omega^2$, non-zero unless ω = ±f. · **plain** The
     current is fixed by the surface slope.
  3. **did** Tidy up · **tex** $\hat u=\dfrac{g\hat\eta}{\omega^2-f^2}(\omega k+ifl),\qquad\hat v=\dfrac{g\hat\eta}{\omega^2-f^2}(-ifk+\omega l)$ (13.80) · **why** Multiply out with $i^2=-1$ and change
     the sign of numerator and denominator. · **plain** Part of the current is in step with the surface (the ω terms),
     part a quarter period out of step (the terms with if).
  4. **did** Form the divergence amplitude · **tex** $k\hat u+l\hat v=\dfrac{g\hat\eta\,\omega\,(k^2+l^2)}{\omega^2-f^2}$ · **why** The terms $ifkl$ and $-ifkl$ cancel; only
     the in-phase parts of the current converge. · **plain** The out-of-phase, rotating part of the current does not
     pile water up.
  5. **did** Insert into continuity · **tex** $\omega^2-f^2=gH(k^2+l^2)$ (13.81) · **why** The third equation says $\omega\hat\eta=H(k\hat u+l\hat v)$; substitute
     step 4 and cancel $\omega\hat\eta\neq0$. · **plain** The wave can exist only at this frequency.
  6. **did** Use the magnitude of the wavenumber · **tex** $\omega^2=f^2+gHK^2,\qquad K=\sqrt{k^2+l^2}$ (13.82) · **why** A rearrangement; only K
     appears, so the wave behaves the same in every horizontal direction (isotropy). · **plain** Gravity and rotation
     add their squared frequencies.
  7. **did** Differentiate for the group velocity (ours) · **tex** $2\omega\dfrac{\partial\omega}{\partial k}=2c^2k\ \Rightarrow\ \mathbf c_g=\dfrac{c^2\mathbf K}{\omega},\qquad c^2=gH$ · **why**
     Implicit differentiation of step 6 with respect to k (and likewise l), by the chain rule; the group velocity is
     the gradient of ω in wavenumber space. · **plain** Energy travels along the wavevector, more slowly for longer
     waves.
  8. **did** Compare with the phase speed (ours) · **tex** $c_p=\dfrac\omega K,\qquad c_p\,c_g=c^2,\qquad c_g<c<c_p$ · **why** Multiply $c_p=\omega/K$ by
     $c_g=c^2K/\omega$; since ω > cK the phase speed exceeds c and the group speed falls short of it. · **plain** Crests
     outrun the non-rotating speed; energy lags behind it.
- **Result:** $\omega^2=f^2+gHK^2$ with $\hat u=\dfrac{g\hat\eta}{\omega^2-f^2}(\omega k+ifl)$, $\hat v=\dfrac{g\hat\eta}{\omega^2-f^2}(-ifk+\omega l)$ — "rotational gravity waves: dispersive, isotropic,
  never slower than f".
- **Check:** units — gHK² = (m²/s²)(1/m²) = s⁻² ✓. Limits — f = 0: ω = cK (Chapter 7); K → 0: ω → abs(f), group velocity → 0
  (inertial oscillation: energy does not travel). Number — 35° N, c = 202.95 m/s, 3100 km: ω = 5.018 f, c_p = 207.1 m/s, c_g =
  198.9 m/s.
- **What it means:** rotation matters for waves longer than about 2π times the Rossby radius; the current vector turns
  clockwise for f > 0 (counter-clockwise for f < 0) with axis ratio ω/abs(f). **Fails when:** f varies across the wave (then
  use D13's cubic) or a boundary is near (D15).
- **Traps:** dividing by ω² − f² needs ω ≠ ±f; as K → 0 the frequency tends to f, not to zero; the "orbit" figure of the
  book is a velocity hodograph (trap T7).

### D15 · The Kelvin wave: (13.84) → (13.86) → (13.87), and the trapping width — ★★, 8 steps, in C11 (notebook · `kelvin_wave`)
- **Goal:** find the wave that can run along a straight coast with no flow across the shore, its speed, and how it decays
  offshore. · **Start:** $\dfrac{\partial\eta}{\partial t}+H\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$, $\dfrac{\partial u}{\partial t}-fv=-g\dfrac{\partial\eta}{\partial x}$, $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.45), a wall along y = 0
  with the sea in y > 0 — "shallow water next to a coast". · **Plan:** • set v = 0 everywhere • look for a wave along the
  coast with an unknown offshore shape • get the speed from two equations and the shape from the third • keep the
  shape that decays. · **Tools:** plane wave with an unknown structure function; first-order linear ODE (P210, reminder);
  trapped solutions (P323). · **Assumptions:** f-plane; straight vertical wall; v ≡ 0 (step 1, justified at the end).
- **Steps:**
  1. **did** Set v = 0 everywhere · **tex** $\dfrac{\partial\eta}{\partial t}+H\dfrac{\partial u}{\partial x}=0,\qquad\dfrac{\partial u}{\partial t}=-g\dfrac{\partial\eta}{\partial x},\qquad fu=-g\dfrac{\partial\eta}{\partial y}$ (13.84) · **why** The wall forces v = 0 at
     y = 0; we try the simplest possibility, v = 0 for all y, and must check at the end that all three equations can
     still be met. · **plain** Along the shore a plain gravity wave; across it, geostrophic balance.
  2. **did** Try a wave along the coast · **tex** $[u,\eta]=[\hat u(y),\hat\eta(y)]\,e^{i(kx-\omega t)}$ · **why** The coefficients do not depend on x or t, so
     exponentials in x and t work; the dependence on y is left open. · **plain** A travelling wave whose strength may
     vary offshore.
  3. **did** Substitute · **tex** $-i\omega\hat\eta+iHk\hat u=0,\qquad-i\omega\hat u=-igk\hat\eta,\qquad f\hat u=-g\dfrac{d\hat\eta}{dy}$ (13.85) · **why** ∂/∂t → −iω, ∂/∂x → ik; the y-derivative
     stays. · **plain** Two algebraic relations and one differential equation.
  4. **did** Eliminate û between the first two · **tex** $\hat\eta\,\big[\omega^2-gHk^2\big]=0\ \Rightarrow\ \omega=\pm ck,\qquad c=\sqrt{gH}$ (13.86) · **why** The second gives
     $\hat u=(gk/\omega)\hat\eta$; insert into the first and multiply by ω. For a wave to exist the bracket must vanish. · **plain**
     The wave moves at the ordinary long-wave speed — rotation does not appear.
  5. **did** Use the cross-shore equation · **tex** $\dfrac{d\hat\eta}{dy}=-\dfrac{fk}{\omega}\hat\eta=\mp\dfrac fc\hat\eta$ · **why** Insert $\hat u=(gk/\omega)\hat\eta$ into $f\hat u=-g\,d\hat\eta/dy$;
     the upper sign goes with ω = +ck (travel toward +x). · **plain** The offshore shape obeys the simplest possible
     ODE.
  6. **did** Solve and keep the decaying solution · **tex** $\hat\eta=\eta_0\,e^{\mp fy/c};\qquad f>0\ \Rightarrow\ \text{upper sign}$ · **why** A first-order linear
     ODE has an exponential solution; in the half-plane y > 0 only the decaying one has finite energy. For f > 0 that is
     the wave travelling toward +x. · **plain** The wave must travel with the coast on its right (f > 0).
  7. **did** Take real parts · **tex** $\eta=\eta_0e^{-fy/c}\cos k(x-ct),\qquad u=\eta_0\sqrt{\dfrac gH}\,e^{-fy/c}\cos k(x-ct)$ (13.87) · **why** Real part of the complex
     wave; and $\hat u=(gk/\omega)\hat\eta=(g/c)\hat\eta=\sqrt{g/H}\,\hat\eta$. · **plain** Under a crest the current flows in the direction of
     travel; both fade offshore.
  8. **did** Read off the width and the general rule · **tex** $\Lambda\equiv\dfrac{c}{\lvert f\rvert};\qquad\text{trapped}\iff\mathrm{sgn}(f)\times(\text{direction of travel})>0$ ·
     **why** The exponent is −y/Λ; for f < 0, or a coast on the other side, the decaying root of step 6 is the other
     one, and the direction reverses. · **plain** The wave hugs the coast within one Rossby radius, keeping it on the
     right in the northern hemisphere and on the left in the southern.
- **Result:** $\eta=\eta_0e^{-fy/c}\cos k(x-ct)$, $u=\eta_0\sqrt{g/H}\,e^{-fy/c}\cos k(x-ct)$, $c=\sqrt{gH}$ — "a non-dispersive wave trapped within
  Λ = c/f of the coast".
- **Check:** units — η₀√(g/H) = m · s⁻¹ ✓. Consistency of v ≡ 0 — the y-equation becomes exact geostrophy, satisfied by
  step 5; no contradiction. Limits — f → 0: Λ → ∞, a plane wave along the wall. Number — H = 4200 m, 35° N: c = 203 m/s, Λ =
  2426 km; internal two-layer mode: c = 1.88 m/s, Λ = 22.5 km.
- **What it means:** the only shallow-water wave that exists below the frequency f; it carries tides round basins and
  thermocline signals along coasts and (in an equatorial version beyond the book) along the equator. **Fails when:** the
  coast curves on a scale shorter than Λ, or the wavelength is so long that β matters.
- **Traps:** which sign of the exponent goes with which direction; for f < 0 everything mirrors (trap T5); v ≡ 0 is an
  assumption that must be checked against the cross-shore equation, where it leaves exact geostrophy.


### D16 · Geostrophic adjustment of a step: the end state of width Λ = c/f — ★★★, 11 steps, in C12 (notebook · `geostrophic_adjustment`) — ours
- **Goal:** find what is left when a step in surface height is released on a rotating plane and the waves have gone —
  without solving for the waves. **Not in the book — ours** (the book describes in words how geostrophy is set up);
  labelled "analytic (ours)", no citation. · **Start:** $\dfrac{\partial\eta}{\partial t}+H\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$, $\dfrac{\partial u}{\partial t}-fv=-g\dfrac{\partial\eta}{\partial x}$, $\dfrac{\partial v}{\partial t}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.45) on an f-plane,
  nothing depending on y, with $\eta=\eta_0\,\mathrm{sgn}(x)$ and $u=v=0$ at t = 0 — "a step of height 2η₀ at rest". · **Plan:** • find a
  quantity that cannot change at any point • write down what a steady state must satisfy • combine the two into one
  equation for the final surface • solve it on each side and match. · **Tools:** adjustment problems (P324); second-order
  linear ODE with constant coefficients (P44, reminder); trapped solutions (P323). · **Assumptions:** linear (η₀ ≪ H);
  inviscid; constant f > 0 (the sign of f enters only in step 10); unbounded in x.
- **Steps:**
  1. **did** Drop the y-derivatives · **tex** $\dfrac{\partial\eta}{\partial t}+H\dfrac{\partial u}{\partial x}=0,\qquad\dfrac{\partial u}{\partial t}-fv=-g\dfrac{\partial\eta}{\partial x},\qquad\dfrac{\partial v}{\partial t}+fu=0$ · **why** The step is uniform along
     y and nothing can make it otherwise, so ∂/∂y = 0 for all time; v itself need not vanish. · **plain** A
     one-dimensional problem with a velocity component along the step.
  2. **did** Differentiate the third equation in x and use the first · **tex** $\dfrac{\partial}{\partial t}\dfrac{\partial v}{\partial x}-\dfrac fH\dfrac{\partial\eta}{\partial t}=0$ · **why** ∂/∂x of the third
     gives $\partial^2v/\partial x\partial t+f\,\partial u/\partial x=0$; the first equation replaces $\partial u/\partial x$ by $-\frac1H\partial\eta/\partial t$. · **plain** The vorticity of the
     jet grows exactly as the surface drops.
  3. **did** Recognise a conserved quantity · **tex** $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial v}{\partial x}-\dfrac{f\eta}{H}\Big)=0\ \Rightarrow\ \dfrac{\partial v}{\partial x}-\dfrac{f\eta}{H}=-\dfrac{f\eta_0}{H}\,\mathrm{sgn}(x)$ · **why** A quantity
     whose time derivative is zero keeps its initial value at each x; initially v = 0 and η = η₀ sgn(x). This is the
     linear form of potential vorticity (C13). · **plain** Each column remembers the surface height it started with.
  4. **did** Write what a steady end state requires · **tex** $u=0,\qquad fv=g\dfrac{\partial\eta}{\partial x}$ · **why** With ∂/∂t = 0 the first equation gives
     ∂u/∂x = 0, so u = 0 (it vanishes far away); the second is then geostrophic balance. · **plain** In the end the flow
     runs along the step, in balance with the slope across it.
  5. **did** Express the conserved quantity through η · **tex** $\dfrac gf\dfrac{\partial^2\eta}{\partial x^2}-\dfrac fH\eta=-\dfrac{f\eta_0}{H}\,\mathrm{sgn}(x)$ · **why** Differentiate step 4,
     $\partial v/\partial x=(g/f)\,\partial^2\eta/\partial x^2$, and insert into step 3. · **plain** One equation for the final surface alone.
  6. **did** Divide by g/f · **tex** $\dfrac{d^2\eta}{dx^2}-\dfrac{\eta}{\Lambda^2}=-\dfrac{\eta_0}{\Lambda^2}\,\mathrm{sgn}(x),\qquad\Lambda^2=\dfrac{gH}{f^2}=\dfrac{c^2}{f^2}$ · **why** Multiplying by f/g turns the
     coefficient f²/(gH) into 1/Λ²; the Rossby radius appears by itself. · **plain** The only length in the problem is
     Λ = c/f.
  7. **did** Solve on the right-hand side, x > 0 · **tex** $\eta=\eta_0+a\,e^{-x/\Lambda}$ · **why** A particular solution is the constant η₀;
     the homogeneous solutions are $e^{\pm x/\Lambda}$, and the growing one is excluded because η must stay finite. · **plain**
     Far from the step the surface keeps its original height.
  8. **did** Use the symmetry at x = 0 · **tex** $\eta(0)=0\ \Rightarrow\ a=-\eta_0$ · **why** The equation and the initial data are odd in x, so
     the solution is odd and vanishes at x = 0; η and its slope are then continuous there (a jump would mean infinite
     velocity). · **plain** The step is smoothed, not removed.
  9. **did** Write the solution for both sides · **tex** $\eta=\eta_0\,\mathrm{sgn}(x)\big(1-e^{-\lvert x\rvert/\Lambda}\big)$ · **why** Mirror the right-hand solution to x < 0
     using oddness. · **plain** The surface changes from −η₀ to +η₀ over a width of a few Rossby radii.
  10. **did** Get the jet from geostrophy · **tex** $v=\dfrac gf\dfrac{\partial\eta}{\partial x}=\dfrac{g\eta_0}{f\Lambda}e^{-\lvert x\rvert/\Lambda}=\dfrac{g\eta_0}{c}e^{-\lvert x\rvert/\Lambda}$ · **why** Differentiate step 9 and use step
      4 with fΛ = c for f > 0; for f < 0 the jet has the opposite sign, $v=\mathrm{sgn}(f)\,(g\eta_0/c)\,e^{-\lvert x\rvert/\Lambda}$. · **plain** A current
      along the step, fastest at the step, with the high side on its right (f > 0).
  11. **did** Compare energies · **tex** $\dfrac{\mathrm{KE}}{\mathrm{PE\ released}}=\dfrac{\tfrac12\rho H\int v^2dx}{\tfrac12\rho g\int(\eta_0^2-\eta^2)\,dx}=\dfrac{\tfrac12\rho g\eta_0^2\Lambda}{\tfrac32\rho g\eta_0^2\Lambda}=\dfrac13$ · **why** Both integrals over the
      whole line are elementary: $\int_0^\infty e^{-2s}ds=\tfrac12$ and $\int_0^\infty(2e^{-s}-e^{-2s})ds=\tfrac32$, each doubled for the two sides. ·
      **plain** A third of the released energy stays in the jet; two thirds leave with the waves.
- **Result:** $\eta=\eta_0\,\mathrm{sgn}(x)\big(1-e^{-\lvert x\rvert/\Lambda}\big)$, $v=\mathrm{sgn}(f)\,\dfrac{g\eta_0}{c}e^{-\lvert x\rvert/\Lambda}$, $u=0$, with $\Lambda=c/\lvert f\rvert$ — "rotation leaves a front one Rossby
  radius wide, held up by a jet along it".
- **Check:** units — gη₀/c = (m/s²)(m)/(m/s) = m/s ✓. Limits — f → 0: Λ → ∞, the surface flattens everywhere and no jet
  remains (the non-rotating dam break); x ≫ Λ: η → ±η₀, untouched. Number — equivalent depth 1.3286 m, η₀ = 0.05 m, 35° N:
  Λ = 43.15 km, jet 0.136 m/s. Model — our forward–backward march averaged over the fifth inertial period matches η to
  0.2 % and v to 0.4 % (the instantaneous field still oscillates at near-inertial frequency by a few per cent).
- **sympy check** (`check_src`; executed as written — every assertion passes; it prints the ratio 1/3):
  ```python
  import sympy as sp                                            # symbolic algebra
  x = sp.symbols("x", positive=True)                            # distance from the step (right-hand side, x > 0)
  eta0, g, H, f = sp.symbols("eta_0 g H f", positive=True)      # step height, gravity, depth, Coriolis (f > 0)
  c = sp.sqrt(g * H)                                            # long-wave speed
  Lam = c / f                                                   # Rossby radius of deformation
  eta = eta0 * (1 - sp.exp(-x / Lam))                           # step 9: the end-state surface for x > 0
  v = g / f * sp.diff(eta, x)                                   # step 4: geostrophic jet, f v = g d(eta)/dx
  pv_end = sp.diff(v, x) - f * eta / H                          # step 3: linear potential vorticity at the end
  pv_start = -f * eta0 / H                                      # ... and at the start (no flow, eta = eta0)
  assert sp.simplify(pv_end - pv_start) == 0                    # step 5: every column kept its value
  assert sp.simplify(sp.diff(eta, x, 2) - eta / Lam**2 + eta0 / Lam**2) == 0   # step 6: the ODE is satisfied
  assert eta.subs(x, 0) == 0                                    # step 8: eta is continuous at x = 0 (odd in x)
  assert sp.simplify(v - g * eta0 / c * sp.exp(-x / Lam)) == 0  # step 10: the jet, largest at the step
  PE = sp.integrate(sp.Rational(1, 2) * g * (eta0**2 - eta**2), (x, 0, sp.oo))   # potential energy released
  KE = sp.integrate(sp.Rational(1, 2) * H * v**2, (x, 0, sp.oo))                 # kinetic energy of the jet
  assert sp.simplify(KE / PE - sp.Rational(1, 3)) == 0          # step 11: one third stays as kinetic energy
  print("D16 checks passed: PV kept, ODE, matching, jet, KE/PE released =", sp.simplify(KE / PE))
  ```
- **What it means:** the Rossby radius is the scale beyond which rotation prevents gravity from flattening things; ocean
  fronts and eddies, and the atmosphere's highs and lows, are adjusted states of this kind. **Fails when:** the step is
  not small (nonlinear adjustment keeps the same idea with full potential vorticity); f varies (the end state then
  drifts as Rossby waves, C15); friction acts (it slowly spins the jet down, C06).
- **Traps:** the end state is not rest; the conserved quantity is fixed point by point by the initial state; η and its
  slope are continuous at x = 0; the transient waves are not part of this derivation; the jet's direction flips with the
  sign of f.

### D17 · Potential vorticity: (13.88)–(13.90) → (13.91) → (13.92) → (13.93) → (13.94) — ★★★, 12 steps, in C13 (notebook)
- **Goal:** show that each column of a shallow layer keeps its ratio of absolute vorticity to depth, starting from the
  full nonlinear equations with f varying northward. · **Start:** $\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}-fv=-g\dfrac{\partial\eta}{\partial x}$ (13.88),
  $\dfrac{\partial v}{\partial t}+u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}+fu=-g\dfrac{\partial\eta}{\partial y}$ (13.89), $\dfrac{\partial h}{\partial t}+\dfrac{\partial}{\partial x}(uh)+\dfrac{\partial}{\partial y}(vh)=0$ (13.90), with $f=f_0+\beta y$ — "nonlinear shallow water
  over an uneven bottom". · **Plan:** • cross-differentiate the momentum equations to remove the pressure • regroup the
  nonlinear terms into advection of vorticity plus vorticity times divergence • replace the divergence by the rate of
  change of depth • combine with a quotient rule. · **Tools:** cross-differentiation (P177); product rule (P38); material
  derivative (Ch. 3); quotient rule (P265) — all reminders; "materially conserved" (P325). · **Assumptions:** shallow water
  (hydrostatic, homogeneous, columnar motion); inviscid; the book's f → f₀ replacement in steps 3–9 (undone in step 10).
- **Steps:**
  1. **did** Differentiate the y-momentum equation in x · **tex** $\dfrac{\partial}{\partial t}\dfrac{\partial v}{\partial x}+\dfrac{\partial}{\partial x}\Big[u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}\Big]+\dfrac{\partial(fu)}{\partial x}=-g\dfrac{\partial^2\eta}{\partial x\,\partial y}$ · **why** Every term is
     differentiable; we aim at the mixed derivative of η so that it can be cancelled. · **plain** The first half of the
     vorticity equation.
  2. **did** Differentiate the x-momentum equation in y · **tex** $\dfrac{\partial}{\partial t}\dfrac{\partial u}{\partial y}+\dfrac{\partial}{\partial y}\Big[u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}\Big]-\dfrac{\partial(fv)}{\partial y}=-g\dfrac{\partial^2\eta}{\partial y\,\partial x}$ · **why** The same
     move on the other equation gives the same mixed derivative (Schwarz). · **plain** The second half.
  3. **did** Subtract step 2 from step 1 · **tex**
     $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial v}{\partial x}-\dfrac{\partial u}{\partial y}\Big)+\dfrac{\partial}{\partial x}\Big[u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}\Big]-\dfrac{\partial}{\partial y}\Big[u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}\Big]+f_0\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)+\beta v=0$ (13.91) · **why** The pressure terms cancel. The Coriolis terms give
     $\partial(fu)/\partial x+\partial(fv)/\partial y=f(\partial u/\partial x+\partial v/\partial y)+\beta v$, since f depends on y only; the book then writes f₀ for the undifferentiated f.
     · **plain** The term βv appears because the Coriolis parameter changes northward.
  4. **did** Expand the two brackets · **tex**
     $u_xv_x+u\,v_{xx}+v_xv_y+v\,v_{xy}-u_yu_x-u\,u_{xy}-v_yu_y-v\,u_{yy}$ (subscripts = partial derivatives) · **why** Product rule on each of the four
     products; eight terms result. The book says they rearrange "easily". · **plain** The nonlinear terms, written out.
  5. **did** Regroup them · **tex** $u\dfrac{\partial\zeta}{\partial x}+v\dfrac{\partial\zeta}{\partial y}+\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)\zeta,\qquad\zeta\equiv\dfrac{\partial v}{\partial x}-\dfrac{\partial u}{\partial y}$ · **why** Terms 2 and 6 are $u\,\partial\zeta/\partial x$, terms 4
     and 8 are $v\,\partial\zeta/\partial y$; the remaining four factor as $(u_x+v_y)(v_x-u_y)$. · **plain** Vorticity is carried by the flow
     and concentrated by convergence.
  6. **did** Collect into a vorticity equation · **tex** $\dfrac{D\zeta}{Dt}+(\zeta+f_0)\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)+\beta v=0$ (13.92) · **why** Steps 3 and 5 with
     $\frac D{Dt}\equiv\frac\partial{\partial t}+u\frac\partial{\partial x}+v\frac\partial{\partial y}$, the material derivative following the horizontal motion. · **plain** A column's relative vorticity
     changes when the flow converges and when the column moves north or south.
  7. **did** Rewrite continuity with a product rule · **tex** $\dfrac{Dh}{Dt}+h\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)=0$ · **why** In the third Start equation,
     $\partial(uh)/\partial x=u\,\partial h/\partial x+h\,\partial u/\partial x$ and likewise in y; the advective parts join ∂h/∂t. · **plain** Convergence makes the
     column taller.
  8. **did** Eliminate the divergence · **tex** $\dfrac{D\zeta}{Dt}=\dfrac{\zeta+f_0}{h}\dfrac{Dh}{Dt}-\beta v$ · **why** Step 7 gives the divergence as $-\frac1h\frac{Dh}{Dt}$; insert
     it into step 6. · **plain** Stretching a column spins it up, in proportion to the vorticity it already has —
     including the planet's.
  9. **did** Absorb the β term · **tex** $\dfrac{D(\zeta+f)}{Dt}=\dfrac{\zeta+f_0}{h}\dfrac{Dh}{Dt}$ (13.93) · **why** Since f depends only on y,
     $\frac{Df}{Dt}=\frac{\partial f}{\partial t}+u\frac{\partial f}{\partial x}+v\frac{\partial f}{\partial y}=v\beta$; move βv to the left. · **plain** The absolute vorticity changes only by
     stretching.
  10. **did** Restore f in the coefficient · **tex** $\dfrac{D(\zeta+f)}{Dt}=\dfrac{\zeta+f}{h}\dfrac{Dh}{Dt}$ · **why** The book assumes βy ≪ f₀, so f₀ and f are
      interchangeable there. In fact the replacement in step 3 was never needed: keeping f gives this line exactly
      (checked below). · **plain** The equation is exact for shallow water.
  11. **did** Apply the quotient rule to (ζ + f)/h · **tex** $\dfrac{D}{Dt}\Big(\dfrac{\zeta+f}{h}\Big)=\dfrac1{h^2}\Big[h\dfrac{D(\zeta+f)}{Dt}-(\zeta+f)\dfrac{Dh}{Dt}\Big]$ · **why** D/Dt is a
      derivative, so the quotient rule holds for it; we form the ratio because step 10 suggests it. · **plain** The rate
      of change of the ratio, in terms of things we know.
  12. **did** Insert step 10 · **tex** $\dfrac{D}{Dt}\Big(\dfrac{\zeta+f}{h}\Big)=0$, $f=f_0+\beta y$ (13.94) · **why** By step 10 the square bracket of step 11 is
      zero. · **plain** Following a column, (ζ + f)/h never changes.
- **Result:** $\dfrac{D}{Dt}\Big(\dfrac{\zeta+f}{h}\Big)=0$ — "potential vorticity is conserved following the motion".
- **Check:** units — (1/s)/m, the same on both sides of every line ✓. Limits — constant h: absolute vorticity ζ + f is
  conserved (the barotropic equation of §13.16); constant f: Chapter 5's column result. Number — a 5 % squash at 35° N:
  ζ = −0.05 f = −4.18 × 10⁻⁶ s⁻¹.
- **sympy check** (`check_src`; executed as written — every assertion passes and the final line is printed):
  ```python
  import sympy as sp                                            # symbolic algebra
  x, y, t = sp.symbols("x y t", real=True)                      # east, north, time
  g, f0, beta = sp.symbols("g f_0 beta", real=True)             # gravity, f0 and beta
  u = sp.Function("u")(x, y, t)                                 # eastward velocity
  v = sp.Function("v")(x, y, t)                                 # northward velocity
  h = sp.Function("h")(x, y, t)                                 # total depth of the layer
  b = sp.Function("b")(x, y)                                    # height of the uneven bottom
  eta = h + b                                                   # surface height above the reference level
  f = f0 + beta * y                                             # the EXACT f of the beta-plane (no f -> f0 step)
  ut = -u * sp.diff(u, x) - v * sp.diff(u, y) + f * v - g * sp.diff(eta, x)   # (13.88) solved for du/dt
  vt = -u * sp.diff(v, x) - v * sp.diff(v, y) - f * u - g * sp.diff(eta, y)   # (13.89) solved for dv/dt
  ht = -sp.diff(u * h, x) - sp.diff(v * h, y)                   # (13.90) solved for dh/dt
  zeta = sp.diff(v, x) - sp.diff(u, y)                          # relative vorticity
  zeta_t = sp.diff(vt, x) - sp.diff(ut, y)                      # steps 1-3: cross-differentiate, pressure drops out
  div = sp.diff(u, x) + sp.diff(v, y)                           # horizontal divergence
  D = lambda F, Ft: Ft + u * sp.diff(F, x) + v * sp.diff(F, y)  # material derivative with a known dF/dt
  step6 = D(zeta, zeta_t) + (zeta + f) * div + beta * v         # (13.92) with f kept exact
  assert sp.simplify(sp.expand(step6)) == 0                     # steps 4-6: the vorticity equation holds
  Dh = D(h, ht)                                                 # Dh/Dt from continuity
  assert sp.simplify(sp.expand(Dh + h * div)) == 0              # step 7: Dh/Dt = -h (du/dx + dv/dy)
  q = (zeta + f) / h                                            # potential vorticity
  q_t = (zeta_t * h - (zeta + f) * ht) / h**2                   # step 11: quotient rule (f does not depend on t)
  assert sp.simplify(sp.expand(sp.together(D(q, q_t)) * h**2)) == 0   # step 12: Dq/Dt = 0, exactly
  print("D17 checks passed: vorticity equation, continuity in D/Dt form, D/Dt[(zeta + f)/h] = 0")
  ```
- **What it means:** potential vorticity labels a column like a dye; Rossby waves, flow over topography, both instabilities
  and geostrophic turbulence are statements about it. **Fails when:** friction or forcing act (they change it slowly), or
  the layer is not shallow (stratified fluids have their own, Ertel's, version — beyond the book).
- **Traps:** where βv comes from (differentiating fu and fv); the regrouping of the four products; the book's f → f₀ and
  back — the exact result needs neither (trap T10); h is the total depth over an uneven bottom, not the mean depth.

### D18 · Flow over a step: the vorticity just downstream and the streamline equation — ★★, 7 steps, in C13 (notebook) — ours
- **Goal:** turn the book's verbal argument about a current crossing a step into an equation for the streamline, for
  eastward and westward flow. **Ours** (the book argues in words). · **Start:** $\dfrac{D}{Dt}\Big(\dfrac{\zeta+f}{h}\Big)=0$, $f=f_0+\beta y$ (13.94); a
  uniform current U along x with no relative vorticity upstream; depth h₀ before the step at x = 0 and h₁ after it. ·
  **Plan:** • apply conservation along a streamline • describe the streamline by its northward displacement Y(x) • write ζ
  in terms of Y • solve the resulting ODE for U > 0 and U < 0. · **Tools:** steady flow: streamlines are particle paths
  (Ch. 3, reminder); second-order linear ODE (P44, reminder); matching at a junction (P323). · **Assumptions:** steady;
  small deflection, $\lvert dY/dx\rvert\ll1$ (step 3); β-plane.
- **Steps:**
  1. **did** Apply conservation along a streamline · **tex** $\dfrac{\zeta+f}{h}=\dfrac{f_0}{h_0}$ · **why** In steady flow a column stays on one
     streamline, so its potential vorticity keeps the upstream value, where ζ = 0, f = f₀ and h = h₀. · **plain**
     Whatever happens downstream, this ratio is fixed.
  2. **did** Describe the streamline by its displacement · **tex** $f=f_0+\beta Y(x)$ · **why** A streamline that started at y = 0 and
     has moved north by Y sits where the Coriolis parameter is f₀ + βY. · **plain** Moving north raises f.
  3. **did** Express the vorticity through Y · **tex** $v=U\dfrac{dY}{dx},\qquad\zeta=\dfrac{\partial v}{\partial x}=U\dfrac{d^2Y}{dx^2}$ · **why** The flow is tangent to the
     streamline, so v/u = dY/dx with u ≈ U for a small deflection; and ∂u/∂y is negligible to the same order. · **plain**
     The curvature of the path is its relative vorticity (divided by U).
  4. **did** Substitute into step 1 · **tex** $\dfrac{d^2Y}{dx^2}+\dfrac{\beta}{U}Y=\dfrac{f_0}{U}\,\dfrac{h-h_0}{h_0}$ · **why** Multiply step 1 by h, insert steps 2 and 3, and
     divide by U; the right-hand side is zero before the step and a constant after it. · **plain** An oscillator
     equation if U > 0, an exponential one if U < 0.
  5. **did** Evaluate just after the step · **tex** $\zeta=U\dfrac{d^2Y}{dx^2}\Big|_{x=0^+}=\dfrac{f_0(h_1-h_0)}{h_0}$ · **why** At the step the streamline has not yet
     moved (Y = 0), so step 4 gives the curvature directly. This is the book's $\zeta=f(h_1-h_0)/h_0<0$. · **plain** A column
     squashed by the step turns clockwise (northern hemisphere).
  6. **did** Solve for eastward flow, U > 0 · **tex** $Y=Y_p\big(1-\cos\kappa x\big)\ (x>0),\qquad\kappa=\sqrt{\dfrac\beta U},\qquad Y_p=\dfrac{f_0(h_1-h_0)}{\beta h_0}$ · **why** Upstream
     Y = 0; downstream the particular solution is the constant Y_p and the homogeneous ones are cos κx and sin κx; Y
     and dY/dx are continuous at x = 0. · **plain** A standing meander of wavelength $2\pi\sqrt{U/\beta}$ about a new mean
     latitude.
  7. **did** Solve for westward flow, U < 0 · **tex** $Y=\tfrac12Y_pe^{-\mu x}\ (x>0),\qquad Y=Y_p\big(1-\tfrac12e^{\mu x}\big)\ (x<0),\qquad\mu=\sqrt{\dfrac{\beta}{\lvert U\rvert}}$ · **why**
     The fluid now arrives from x > 0. The homogeneous solutions are exponentials; keep the ones that stay bounded on
     each side and match Y and dY/dx at the step. · **plain** No oscillation: the stream starts to turn before it
     reaches the step and settles at the new latitude after it.
- **Result:** $Y''+\dfrac\beta UY=\dfrac{f_0}{U}\dfrac{h-h_0}{h_0}$; eastward flow: $Y=Y_p(1-\cos\kappa x)$ with wavelength $2\pi\sqrt{U/\beta}$; westward flow:
  exponential on both sides with e-folding length $\sqrt{\lvert U\rvert/\beta}$ — "eastward: a wake of stationary waves; westward:
  upstream influence and no wake".
- **Check:** units — β/U = 1/m² ✓; Y_p = (1/s)/(1/(m s)) = m ✓. Vorticity jump — in both cases ζ changes by f₀(h₁ − h₀)/h₀
  across the step ✓. Number (U = 17 m/s, a 5 % step, 35° N) — Y_p = −223 km, wavelength 5983 km, largest slope 0.23 (small,
  as assumed); westward: e-folding 952 km, Y = −112 km at the step. Link — $2\pi\sqrt{U/\beta}$ is the stationary Rossby
  wavelength of D23 step 9.
- **What it means:** westerlies crossing a mountain range leave a train of stationary waves downstream; easterlies do not.
  The asymmetry is the Rossby wave's: its phase can only go west, so it can stand still only in an eastward current.
  **Fails when:** the step is high (deflection not small), the flow is unsteady, or friction damps the wake.
  **Difference from the book's sketch:** for westward flow the book draws the stream back at its original latitude far
  downstream; our solution ends displaced by Y_p, as conservation of potential vorticity over a permanent change of
  depth requires (reported to the verifier).
- **Traps:** Y is a displacement, not a velocity; for U < 0 "upstream" is x > 0; the bounded solution, not the one
  starting at the step, must be chosen on each side; f in the book's jump formula is the upstream value.

### D19 · The w-equation with rotation, (13.96), from the set (13.95) — ★★★, 11 steps, in C14 (notebook)
- **Goal:** reduce five equations for u, v, w, p, ρ to one equation for the vertical velocity — the starting point for
  every internal wave in a rotating stratified fluid. The book refers to its Chapter 7 and never shows the steps with
  rotation. · **Start:** the five members of (13.95): $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$; $\dfrac{\partial u}{\partial t}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$; $\dfrac{\partial v}{\partial t}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$;
  $\dfrac{\partial w}{\partial t}=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial z}-\dfrac{\rho g}{\rho_0}$; $\dfrac{\partial\rho}{\partial t}-\dfrac{\rho_0N^2}{g}w=0$ — "linear, Boussinesq, rotating, not hydrostatic; N may depend on z". ·
  **Plan:** • take the horizontal divergence of the momentum equations • take their curl • combine to remove u and v •
  combine the vertical and density equations to remove ρ • eliminate p between the two results. · **Tools:** operator
  elimination (P177, reminder); the Chapter 7 derivation of $\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0$ (7.134), whose f = 0 pattern this
  follows. · **Assumptions:** linear; f constant; N = N(z) (it is never differentiated — see step 9).
- **Steps:**
  1. **did** Differentiate the x-momentum equation in x · **tex** $\dfrac{\partial^2u}{\partial t\,\partial x}-f\dfrac{\partial v}{\partial x}=-\dfrac1{\rho_0}\dfrac{\partial^2p}{\partial x^2}$ · **why** Both f and ρ₀ are constants; we
     prepare the horizontal divergence. · **plain** Half of a divergence equation.
  2. **did** Differentiate the y-momentum equation in y · **tex** $\dfrac{\partial^2v}{\partial t\,\partial y}+f\dfrac{\partial u}{\partial y}=-\dfrac1{\rho_0}\dfrac{\partial^2p}{\partial y^2}$ · **why** The matching move on the
     other equation. · **plain** The other half.
  3. **did** Add them · **tex** $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}\Big)-f\zeta=-\dfrac1{\rho_0}\nabla_H^2p,\qquad\zeta\equiv\dfrac{\partial v}{\partial x}-\dfrac{\partial u}{\partial y}$ · **why** The Coriolis terms combine into
     −f times the vertical vorticity; the pressure terms into the horizontal Laplacian. · **plain** Divergence changes
     because of vorticity and pressure.
  4. **did** Replace the divergence by continuity · **tex** $-\dfrac{\partial^2w}{\partial t\,\partial z}-f\zeta=-\dfrac1{\rho_0}\nabla_H^2p\quad\text{(A)}$ · **why** Continuity lets us trade the horizontal divergence for the vertical derivative of w:
     $\partial u/\partial x+\partial v/\partial y=-\partial w/\partial z$. · **plain** Now u and v appear only through ζ.
  5. **did** Form the vorticity equation · **tex** $\dfrac{\partial\zeta}{\partial t}=f\dfrac{\partial w}{\partial z}\quad\text{(B)}$ · **why** ∂/∂x of the y-equation minus ∂/∂y of the
     x-equation: the pressure cancels (Schwarz), leaving $\zeta_t+f(u_x+v_y)=0$; then continuity. This is the new step
     that rotation requires. · **plain** Stretching a column creates vorticity from the planet's spin.
  6. **did** Differentiate (A) in time and insert (B) · **tex** $-\dfrac{\partial^3w}{\partial t^2\,\partial z}-f^2\dfrac{\partial w}{\partial z}=-\dfrac1{\rho_0}\nabla_H^2\dfrac{\partial p}{\partial t}\quad\text{(C)}$ · **why** The time
     derivative of −fζ is −f² ∂w/∂z by (B); u, v and ζ are now all eliminated. · **plain** One relation between w and p.
  7. **did** Differentiate the vertical momentum equation in time · **tex** $\dfrac{\partial^2w}{\partial t^2}=-\dfrac1{\rho_0}\dfrac{\partial^2p}{\partial z\,\partial t}-\dfrac g{\rho_0}\dfrac{\partial\rho}{\partial t}$ · **why** We want
     ∂ρ/∂t to appear, because the density equation gives it in terms of w. · **plain** Preparing to remove the density.
  8. **did** Insert the density equation · **tex** $\dfrac{\partial^2w}{\partial t^2}+N^2w=-\dfrac1{\rho_0}\dfrac{\partial^2p}{\partial z\,\partial t}\quad\text{(D)}$ · **why** The fifth Start equation gives
     $\frac g{\rho_0}\frac{\partial\rho}{\partial t}=N^2w$. · **plain** A second relation between w and p: the buoyancy oscillator, forced by pressure.
  9. **did** Apply the horizontal Laplacian to (D) · **tex** $\nabla_H^2\dfrac{\partial^2w}{\partial t^2}+N^2\nabla_H^2w=-\dfrac1{\rho_0}\nabla_H^2\dfrac{\partial^2p}{\partial z\,\partial t}$ · **why** N depends on z
     only, so it passes through the horizontal derivatives untouched — which is why (D), not (C), receives $\nabla_H^2$. ·
     **plain** Dressing (D) so that its pressure term matches the one coming next.
  10. **did** Differentiate (C) in z · **tex** $-\dfrac{\partial^4w}{\partial t^2\,\partial z^2}-f^2\dfrac{\partial^2w}{\partial z^2}=-\dfrac1{\rho_0}\nabla_H^2\dfrac{\partial^2p}{\partial t\,\partial z}$ · **why** (C) contains no N, so differentiating it
      in z produces no derivative of N. · **plain** The same pressure term as in step 9.
  11. **did** Subtract step 10 from step 9 · **tex** $\dfrac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w+f^2\dfrac{\partial^2w}{\partial z^2}=0$ (13.96) · **why** The pressure terms are equal
      and cancel; $\nabla_H^2w_{tt}+w_{zztt}$ is the full Laplacian $\nabla^2=\nabla_H^2+\partial^2/\partial z^2$ of $w_{tt}$. · **plain** One equation for w: inertia,
      buoyancy acting on horizontal variations, rotation acting on vertical ones.
- **Result:** $\dfrac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w+f^2\dfrac{\partial^2w}{\partial z^2}=0$ — "the internal-wave equation with rotation".
- **Check:** units — every term is (m/s)/(s² m²) ✓. Limits — f = 0: Chapter 7's $\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0$ (7.134); N = 0: pure
  inertial waves, $\frac{\partial^2}{\partial t^2}\nabla^2w+f^2w_{zz}=0$. Plane wave — it gives the dispersion relation of D20.
- **sympy check** (`check_src`; executed as written — every assertion passes; it prints the factored 5 × 5 determinant):
  ```python
  import sympy as sp                                            # symbolic algebra
  x, y, z, t = sp.symbols("x y z t", real=True)                 # space and time
  f, rho0, g = sp.symbols("f rho_0 g", positive=True)           # Coriolis parameter, reference density, gravity
  N = sp.Function("N")(z)                                       # buoyancy frequency, allowed to vary with z
  w = sp.Function("w")(x, y, z, t)                              # vertical velocity
  p = sp.Function("p")(x, y, z, t)                              # perturbation pressure
  lapH = lambda F: sp.diff(F, x, 2) + sp.diff(F, y, 2)          # horizontal Laplacian
  C = -sp.diff(w, z, t, t) - f**2 * sp.diff(w, z) + lapH(sp.diff(p, t)) / rho0   # step 6: (C) = 0
  Dq = sp.diff(w, t, 2) + N**2 * w + sp.diff(p, z, t) / rho0    # step 8: (D) = 0
  combo = sp.expand(lapH(Dq) - sp.diff(C, z))                   # steps 9-11: lapH of (D) minus d/dz of (C)
  target = sp.diff(lapH(w) + sp.diff(w, z, 2), t, 2) + N**2 * lapH(w) + f**2 * sp.diff(w, z, 2)   # (13.96)
  assert sp.simplify(combo - sp.expand(target)) == 0            # pressure cancelled; N(z) was never differentiated
  s, a, b, d, Nc = sp.symbols("s a b d N_c")                    # constant-N cross-check: d/dt, d/dx, d/dy, d/dz
  M = sp.Matrix([[a, b, d, 0, 0],                               # continuity     (unknowns u, v, w, p/rho0, g rho/rho0)
                 [s, -f, 0, a, 0],                              # x-momentum
                 [f, s, 0, b, 0],                               # y-momentum
                 [0, 0, s, d, 1],                               # z-momentum
                 [0, 0, -Nc**2, 0, s]])                         # density equation times g/rho0
  det = sp.factor(M.det())                                      # one operator acting on any of the unknowns
  op = s**2 * (a**2 + b**2 + d**2) + Nc**2 * (a**2 + b**2) + f**2 * d**2   # the operator of (13.96)
  assert sp.simplify(det / op).free_symbols <= {s}              # the determinant is (13.96) times a power of d/dt
  print("D19 checks passed: elimination with N(z), and the 5x5 determinant =", det)
  ```
- **What it means:** rotation adds one term, with the *vertical* second derivative: it stiffens the fluid against motions
  that vary along the rotation axis, as buoyancy stiffens it against motions that vary horizontally. **Fails when:** the
  waves are not small, there is a mean shear, or f varies across the wave.
- **Traps:** the vorticity equation (step 5) is the step that brings in f²; the Laplacian is three-dimensional in the
  first term and horizontal in the second; N stays outside every z-derivative because of the order chosen in steps
  9–10; slip #8 — the book prints "depth independent" for N, the correct form is depth dependent.

### D20 · The vertical structure (13.98)–(13.99) and the dispersion relation (13.112) — ★★, 7 steps, in C14 (notebook)
- **Goal:** find which frequencies an internal wave may have in a rotating stratified fluid, and what sets them. ·
  **Start:** $\dfrac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w+f^2\dfrac{\partial^2w}{\partial z^2}=0$ (13.96) and the trial solution $[u,v,w]=[\hat u(z),\hat v(z),\hat w(z)]\,e^{i(kx+ly-\omega t)}$ (13.97) — "the
  Result of D19 and a wave that is plane in the horizontal". · **Plan:** • substitute • collect the equation for the
  vertical structure • name the local vertical wavenumber • for a locally uniform N read off the dispersion relation •
  rewrite it with the angle of the wavevector. · **Tools:** plane-wave substitution (P176, reminder); $1+\tan^2\theta=1/\cos^2\theta$. ·
  **Assumptions:** N varies slowly enough to speak of a local m (steps 4–7; made precise by the WKB note N111).
- **Steps:**
  1. **did** Substitute the trial solution · **tex** $(-i\omega)^2\Big[(ik)^2+(il)^2+\dfrac{d^2}{dz^2}\Big]\hat w+N^2\big[(ik)^2+(il)^2\big]\hat w+f^2\dfrac{d^2\hat w}{dz^2}=0$ · **why** Each ∂/∂t
     becomes −iω, ∂/∂x becomes ik, ∂/∂y becomes il; the z-dependence is kept because N may vary with z. · **plain** An
     ordinary differential equation in z for the wave's vertical structure.
  2. **did** Multiply out and collect · **tex** $\dfrac{d^2\hat w}{dz^2}+\dfrac{(N^2-\omega^2)(k^2+l^2)}{\omega^2-f^2}\hat w=0$ (13.98) · **why** Here $(-i\omega)^2=-\omega^2$ and $(ik)^2=-k^2$; the
     $\hat w''$ terms combine to $(f^2-\omega^2)\hat w''$; divide by $-(\omega^2-f^2)$. · **plain** An oscillator equation in z.
  3. **did** Name the coefficient · **tex** $m^2(z)\equiv\dfrac{(k^2+l^2)\,[N^2(z)-\omega^2]}{\omega^2-f^2}\ \Rightarrow\ \dfrac{d^2\hat w}{dz^2}+m^2\hat w=0$ (13.99)–(13.100) · **why** A
     definition: where m² > 0 the solution oscillates in z with local wavenumber m, where m² < 0 it decays. · **plain**
     Waves need m² > 0, that is, a frequency between f and N.
  4. **did** Specialise to l = 0 and a locally uniform N · **tex** $\hat w\propto e^{\pm imz},\qquad m^2=\dfrac{k^2(N^2-\omega^2)}{\omega^2-f^2}$ (13.109) · **why** With
     constant m the oscillator equation has exponential solutions; orienting x along the horizontal wavevector loses
     nothing (the problem is isotropic in the horizontal). · **plain** A plane wave with wavevector (k, m).
  5. **did** Rearrange · **tex** $\omega^2-f^2=\dfrac{k^2}{m^2}\big(N^2-\omega^2\big)$ (13.112) · **why** Multiply step 4 by $(\omega^2-f^2)/m^2$, so that the two frequency differences stand side by side. · **plain** How far the
     frequency is above f, compared with how far it is below N, is set by the slope k/m.
  6. **did** Introduce the angle of the wavevector · **tex** $\tan\theta=\dfrac mk\ \Rightarrow\ \omega^2\big(1+\tan^2\theta\big)=N^2+f^2\tan^2\theta$ · **why** Here θ is the
     angle between the wavevector and the horizontal (not latitude here); multiply step 5 by tan² θ and gather the ω²
     terms. · **plain** Only the direction of the wavevector appears, not its length.
  7. **did** Multiply by cos² θ · **tex** $\omega^2=f^2\sin^2\theta+N^2\cos^2\theta$ · **why** Use $1+\tan^2\theta=1/\cos^2\theta$ and $\tan^2\theta\cos^2\theta=\sin^2\theta$. The book prints
     this form without a number. · **plain** The squared frequency is a weighted mean of f² and N².
- **Result:** $\omega^2-f^2=\dfrac{k^2}{m^2}(N^2-\omega^2)$, equivalently $\omega^2=f^2\sin^2\theta+N^2\cos^2\theta$ — "the frequency of an inertia–gravity wave depends
  only on the tilt of its wavevector and lies between f and N".
- **Check:** units — s⁻² throughout ✓. Limits — f = 0: ω = N cos θ (Chapter 7); θ = 0 (vertical crests): ω = N; θ = 90°
  (horizontal crests): ω = f. Number — N = 2.7 × 10⁻³ s⁻¹, 35° N: θ = 45° gives ω = 22.8 f; θ = 89° gives 1.148 f. Regimes — near
  N drop f² ($m^2\simeq k^2(N^2-\omega^2)/\omega^2$); near f drop ω² against N² ($\omega^2\simeq f^2+k^2N^2/m^2$, hydrostatic); in between both
  ($m^2\simeq k^2N^2/\omega^2$).
- **What it means:** a wave-maker oscillating at frequency ω sends energy out along beams of one fixed slope; rotation
  flattens those beams as ω approaches f. **Fails when:** N changes within a vertical wavelength (then solve step 3
  numerically, or use the WKB form), or ω is outside the band (no propagation).
- **Traps:** the signs of (−iω)² and (ik)²; θ here is the angle of the wavevector with the horizontal, not latitude (trap
  T13); real m needs f < ω < N; each regime approximation has an error that the notebook tabulates.


### D21 · The group velocity of inertia–gravity waves; phase ⟂ group — ★★, 7 steps, in C14 (notebook)
- **Goal:** find where the energy of an inertia–gravity wave goes. The book leaves it to an exercise. · **Start:**
  $\omega^2-f^2=\dfrac{k^2}{m^2}(N^2-\omega^2)$ (13.112) for a locally uniform N — "the dispersion relation of D20". · **Plan:** • solve it for ω² as
  a function of (k, m) • differentiate with respect to k and to m • put the two components together • compare with the
  phase velocity and the particle motion. · **Tools:** group velocity as the gradient of ω in wavenumber space (P326);
  quotient rule (P265, reminder). · **Assumptions:** uniform N locally; l = 0; N > f.
- **Steps:**
  1. **did** Solve for ω² · **tex** $\omega^2=\dfrac{N^2k^2+f^2m^2}{k^2+m^2}$ · **why** Multiply the Start by m², gather the ω² terms and divide by
     k² + m²; an explicit function is needed before differentiating. · **plain** The frequency as a function of the two
     wavenumbers.
  2. **did** Differentiate with respect to k · **tex** $2\omega\dfrac{\partial\omega}{\partial k}=\dfrac{2N^2k\,(k^2+m^2)-2k\,(N^2k^2+f^2m^2)}{(k^2+m^2)^2}=\dfrac{2km^2(N^2-f^2)}{(k^2+m^2)^2}$ · **why** Chain rule on the left,
     quotient rule on the right; the terms $2N^2k^3$ cancel. · **plain** Differentiating ω² avoids square roots.
  3. **did** Divide by 2ω · **tex** $c_{gx}=\dfrac{\partial\omega}{\partial k}=\dfrac{km^2(N^2-f^2)}{\omega\,(k^2+m^2)^2}$ · **why** Since ω > 0 we may divide by it; the horizontal group velocity is the
     k-derivative of ω. · **plain** Energy moves horizontally in the direction of k.
  4. **did** Do the same for m · **tex** $c_{gz}=\dfrac{\partial\omega}{\partial m}=-\dfrac{k^2m\,(N^2-f^2)}{\omega\,(k^2+m^2)^2}$ · **why** The same two rules as in steps 2–3; now the numerator is
     $2f^2m(k^2+m^2)-2m(N^2k^2+f^2m^2)=-2mk^2(N^2-f^2)$. · **plain** The vertical group velocity has the opposite sign to m.
  5. **did** Insert ω and combine · **tex** $[c_{gx},c_{gz}]=\dfrac{(N^2-f^2)\,km}{(m^2+k^2)^{3/2}(m^2f^2+k^2N^2)^{1/2}}\,[m,\,-k]$ · **why** From step 1,
     $\omega=(N^2k^2+f^2m^2)^{1/2}/(k^2+m^2)^{1/2}$; both components share the prefactor. · **plain** The group velocity points along
     (m, −k).
  6. **did** Dot with the phase velocity · **tex** $\mathbf c\cdot\mathbf c_g\propto(k,\ m)\cdot(m,\ -k)=km-mk=0$ · **why** The phase velocity is along the
     wavevector (k, m). A zero dot product means perpendicular. · **plain** Energy travels along the crests, not across
     them; if the crests move up, the energy moves down.
  7. **did** Compare with the particle motion · **tex** $\dfrac uw=\mp\dfrac mk\ \Rightarrow\ (u,\ w)\parallel(m,\ -k)\parallel\mathbf c_g$ · **why** The relation
     $\dfrac uw=\mp\dfrac mk=\mp\tan\theta$ (13.111) says the motion in the x–z plane is along (m, −k) — the direction just found. ·
     **plain** In the vertical plane the fluid slides back and forth along the direction in which the energy travels.
- **Result:** $[c_{gx},c_{gz}]=\dfrac{(N^2-f^2)\,km}{(m^2+k^2)^{3/2}(m^2f^2+k^2N^2)^{1/2}}\,[m,\,-k]$, with $\mathbf c\cdot\mathbf c_g=0$ — "phase and energy move at right angles,
  with opposite vertical components".
- **Check:** units — (s⁻²)(m⁻²)/(m⁻³ · s⁻¹ m⁻¹) × m⁻¹ = m/s ✓. Limits — f = 0: Chapter 7's result; N = f: no dispersion in
  angle and zero group velocity. Number — N = 2.7 × 10⁻³ s⁻¹, 35° N, a wave 10 km long and 200 m tall: c_g = (+0.0465, −9.31 ×
  10⁻⁴) m/s. Code — central differences of `GFD.inertia_gravity_omega`.
- **What it means:** storms put near-inertial energy into the top of the ocean; it leaks downward along gently sloping
  rays while the phase is seen to move *up*. **Fails when:** N varies within a wavelength (ray tracing with the local N
  then applies).
- **Traps:** differentiate ω², then divide by 2ω; the vertical components of phase and group velocity have opposite
  signs; with f = 0 the result must reduce to Chapter 7's.

### D22 · The quasi-geostrophic vorticity equation (13.117) from potential-vorticity conservation (13.94) — ★★★, 10 steps, in C15 (notebook · `rossby_waves`)
- **Goal:** obtain one linear equation for the surface height of slow, nearly geostrophic motion on a β-plane — the
  equation whose waves are Rossby waves. · **Start:** $\dfrac{D}{Dt}\Big(\dfrac{\zeta+f}{h}\Big)=0$, $f=f_0+\beta y$ (13.94), with a flat bottom so that
  $h=H+\eta$ — "potential vorticity is conserved; depth = mean depth + surface displacement". · **Plan:** • expand the
  conservation law • keep the terms with one small factor • replace the velocities by their geostrophic values — except
  inside the term that came from the divergence • tidy up. · **Tools:** quotient rule (P265, reminder); ordering in a
  small parameter (P328); the geostrophic wind from a surface height (C02, C07). · **Assumptions:** small amplitude
  (steps 2–4); Ro ≪ 1 and ω ≪ f (step 5); βy ≪ f₀; flat bottom.
- **Steps:**
  1. **did** Expand the conservation law · **tex**
     $(H+\eta)\Big(\dfrac{\partial\zeta}{\partial t}+u\dfrac{\partial\zeta}{\partial x}+v\dfrac{\partial\zeta}{\partial y}+\beta v\Big)-(\zeta+f_0)\Big(\dfrac{\partial\eta}{\partial t}+u\dfrac{\partial\eta}{\partial x}+v\dfrac{\partial\eta}{\partial y}\Big)=0$ (13.114) · **why** Quotient rule, then multiply by h²:
     $h\,D(\zeta+f)/Dt-(\zeta+f)\,Dh/Dt=0$; with $Df/Dt=\beta v$, $Dh/Dt=D\eta/Dt$, and f → f₀ in the coefficient. · **plain** Stretching
     and northward motion both change a column's vorticity.
  2. **did** Mark the small quantities · **tex** $\zeta,\ \eta,\ u,\ v=O(\varepsilon);\qquad H,\ f_0,\ \beta=O(1)$ · **why** A wave of small amplitude: every
     wave field carries one factor of a small number ε; the background does not. · **plain** Sorting terms by size.
  3. **did** Identify the products of two small quantities · **tex** $\eta\dfrac{\partial\zeta}{\partial t},\ u\dfrac{\partial\zeta}{\partial x},\ v\dfrac{\partial\zeta}{\partial y},\ \eta\beta v,\ \zeta\dfrac{\partial\eta}{\partial t},\ u\dfrac{\partial\eta}{\partial x},\ v\dfrac{\partial\eta}{\partial y}=O(\varepsilon^2)$ ·
     **why** Each contains two wave fields, so it is smaller than the others by a factor ε and may be neglected. ·
     **plain** These are the nonlinear terms.
  4. **did** Keep the terms of first order · **tex** $H\dfrac{\partial\zeta}{\partial t}+H\beta v-f_0\dfrac{\partial\eta}{\partial t}=0$ (13.115) · **why** What remains of step 1 after removing
     the terms of step 3. · **plain** Vorticity changes by northward motion (βv) and by stretching of the column
     (f₀ ∂η/∂t).
  5. **did** Use geostrophic velocities · **tex** $u\simeq-\dfrac g{f_0}\dfrac{\partial\eta}{\partial y},\qquad v\simeq\dfrac g{f_0}\dfrac{\partial\eta}{\partial x}$ (13.116) · **why** For ω ≪ f the momentum
     equations of D10's Result reduce to geostrophic balance with pressure gρη (≃ marks the approximation). · **plain**
     To lowest order the flow follows the height contours.
  6. **did** Compute the vorticity · **tex** $\zeta=\dfrac{\partial v}{\partial x}-\dfrac{\partial u}{\partial y}=\dfrac g{f_0}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}\Big)$ · **why** Differentiate step 5; g and f₀ are
     constants. · **plain** A hill in the surface is a clockwise vortex (f > 0), a hollow a counter-clockwise one.
  7. **did** Substitute into step 4 · **tex** $\dfrac{gH}{f_0}\dfrac{\partial}{\partial t}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}\Big)+\dfrac{gH\beta}{f_0}\dfrac{\partial\eta}{\partial x}-f_0\dfrac{\partial\eta}{\partial t}=0$ · **why** Steps 5 and 6 replace ζ and v. The
     term $f_0\,\partial\eta/\partial t$ is kept: it came from the divergence, which geostrophic velocities alone would give as zero. ·
     **plain** One equation with η as the only unknown.
  8. **did** Multiply by f₀/(gH) · **tex** $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}\Big)+\beta\dfrac{\partial\eta}{\partial x}-\dfrac{f_0^2}{gH}\dfrac{\partial\eta}{\partial t}=0$ · **why** A constant factor; it makes the
     coefficient of the first term one. · **plain** The stretching term now carries f₀²/(gH).
  9. **did** Introduce the wave speed · **tex** $\dfrac{f_0^2}{gH}=\dfrac{f_0^2}{c^2}=\dfrac1{\Lambda^2},\qquad c=\sqrt{gH}$ · **why** Here $c=\sqrt{gH}$ is the long-wave speed
     (D10); the ratio c/f₀ is the Rossby radius of C12. · **plain** The Rossby radius enters as the scale at which
     stretching matters as much as relative vorticity.
  10. **did** Collect the time derivatives · **tex** $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}-\dfrac{f_0^2}{c^2}\eta\Big)+\beta\dfrac{\partial\eta}{\partial x}=0$ (13.117) · **why** Both the first and
      the last term of step 8 are time derivatives. · **plain** The bracket — relative vorticity minus stretching — can
      change only by the β-effect.
- **Result:** $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}-\dfrac{f_0^2}{c^2}\eta\Big)+\beta\dfrac{\partial\eta}{\partial x}=0$, $c=\sqrt{gH}$ — "the quasi-geostrophic form of the linearised
  potential-vorticity equation".
- **Check:** units — every term is 1/(m s) ✓. Limits — β = 0: the bracket is steady (the adjusted state of D16, where it
  equalled its initial value); c → ∞ (rigid lid): the barotropic vorticity equation. Plane wave — it gives the
  dispersion relation of D23.
- **sympy check** (`check_src`; executed as written — every assertion passes and the final line is printed):
  ```python
  import sympy as sp                                            # symbolic algebra
  x, y, t = sp.symbols("x y t", real=True)                      # east, north, time
  g, H, f0, beta = sp.symbols("g H f_0 beta", positive=True)    # gravity, mean depth, f0, beta
  eps = sp.symbols("epsilon")                                   # step 2: marks every small (wave) quantity
  eta = sp.Function("eta")(x, y, t)                             # surface displacement
  u = sp.Function("u")(x, y, t)                                 # velocities (not yet geostrophic)
  v = sp.Function("v")(x, y, t)
  zeta = sp.diff(v, x) - sp.diff(u, y)                          # relative vorticity
  D = lambda F: sp.diff(F, t) + eps * u * sp.diff(F, x) + eps * v * sp.diff(F, y)   # D/Dt, u and v are small
  full = (H + eps * eta) * (D(eps * zeta) + beta * eps * v) - (eps * zeta + f0) * D(eps * eta)   # step 1: (13.114)
  lin = sp.expand(full).coeff(eps, 1)                           # steps 3-4: keep the terms with one small factor
  assert sp.expand(lin - (H * sp.diff(zeta, t) + H * beta * v - f0 * sp.diff(eta, t))) == 0   # (13.115)
  ug = -g / f0 * sp.diff(eta, y)                                # step 5: geostrophic u, (13.116a)
  vg = g / f0 * sp.diff(eta, x)                                 # step 5: geostrophic v, (13.116b)
  zg = sp.diff(vg, x) - sp.diff(ug, y)                          # step 6: geostrophic vorticity
  assert sp.simplify(zg - g / f0 * (sp.diff(eta, x, 2) + sp.diff(eta, y, 2))) == 0
  qg = H * sp.diff(zg, t) + H * beta * vg - f0 * sp.diff(eta, t)   # step 7: substitute into (13.115)
  c2 = g * H                                                    # step 9: c^2 = gH
  book = sp.diff(sp.diff(eta, x, 2) + sp.diff(eta, y, 2) - f0**2 / c2 * eta, t) + beta * sp.diff(eta, x)   # (13.117)
  assert sp.simplify(sp.expand(qg * f0 / (g * H) - book)) == 0  # steps 8-10: multiply by f0/(gH)
  k, l, w = sp.symbols("k l omega", real=True)                  # D23 step 1: a plane wave
  wave = sp.exp(sp.I * (k * x + l * y - w * t))
  disp = sp.simplify(book.subs(eta, wave).doit() / wave)        # what the wave must satisfy
  sol = sp.solve(disp, w)[0]                                    # solve for the frequency
  assert sp.simplify(sol + beta * k / (k**2 + l**2 + f0**2 / c2)) == 0   # (13.118)
  print("D22 checks passed: (13.115), (13.117) and the dispersion relation (13.118)")
  ```
- **What it means:** "quasi-geostrophic" = geostrophic velocities everywhere *except* in the divergence, whose small
  ageostrophic part is what makes the flow evolve; the fast gravity waves have been filtered out (one time derivative
  instead of three). **Fails when:** the amplitude is large (keep the advection of the bracket by the geostrophic flow —
  the nonlinear quasi-geostrophic equation, beyond the book), near the equator, or over steep topography.
- **Traps:** which products are dropped; the geostrophic velocity may be used everywhere except in the divergence, which
  is why the stretching term survives; f → f₀ except in β; a flat bottom is assumed.

### D23 · Rossby-wave dispersion (13.118), phase speed (13.119), circles, group velocity, the mean-current form (13.120) — ★★, 9 steps, in C15 (notebook · `rossby_waves`)
- **Goal:** find the frequency of a Rossby wave and from it where its crests and its energy go. · **Start:**
  $\dfrac{\partial}{\partial t}\Big(\dfrac{\partial^2\eta}{\partial x^2}+\dfrac{\partial^2\eta}{\partial y^2}-\dfrac{f_0^2}{c^2}\eta\Big)+\beta\dfrac{\partial\eta}{\partial x}=0$ (13.117) and $\eta=\hat\eta\,e^{i(kx+ly-\omega t)}$ — "the Result of D22 and a plane wave". · **Plan:** •
  substitute • read the phase speed • rewrite as a circle in the wavenumber plane • differentiate for the group
  velocity • find where it changes sign • add a mean current. · **Tools:** plane-wave substitution (P176); completing the
  square into a circle (P272, reminder); group velocity as a gradient (P326); quotient rule (P265); Doppler shift (Ch. 7,
  reminder). · **Assumptions:** as D22; uniform mean current in step 9. Shorthand: $F\equiv f_0^2/c^2=1/\Lambda^2$.
- **Steps:**
  1. **did** Substitute the plane wave · **tex** $\omega=-\dfrac{\beta k}{k^2+l^2+f_0^2/c^2}$ (13.118) · **why** ∂/∂t → −iω, ∂/∂x → ik, ∂/∂y → il give
     $-i\omega\,(-k^2-l^2-F)+i\beta k=0$; divide by i and solve for ω. · **plain** The frequency is proportional to β and depends on the
     direction of the wavevector, not only on its length.
  2. **did** Fix the sign convention · **tex** $\omega>0\ \Rightarrow\ k<0$ · **why** The denominator is positive and β > 0, so ω and k have
     opposite signs; the book takes ω positive, the code returns a signed ω for a signed k. · **plain** With positive
     frequency the wavevector always has a westward component.
  3. **did** Form the zonal phase speed · **tex** $c_x=\dfrac\omega k=-\dfrac{\beta}{k^2+l^2+f_0^2/c^2}$ (13.119) · **why** The phase of the wave is constant
     along x = (ω/k)t. The right side is negative for every k and l. · **plain** The crests of every Rossby wave drift
     westward.
  4. **did** Rearrange into a circle · **tex** $\Big(k+\dfrac{\beta}{2\omega}\Big)^2+l^2=\Big(\dfrac{\beta}{2\omega}\Big)^2-\dfrac{f_0^2}{c^2}$ · **why** Multiply step 1 by the denominator and divide
     by ω: $k^2+(\beta/\omega)k+l^2=-F$; complete the square in k. · **plain** All waves of one frequency lie on a circle
     centred on the negative k-axis.
  5. **did** Differentiate with respect to k · **tex** $c_{gx}=\dfrac{\partial\omega}{\partial k}=\dfrac{\beta\,(k^2-l^2-F)}{(k^2+l^2+F)^2}$ · **why** Quotient rule on step 1:
     $-\beta[(k^2+l^2+F)-2k^2]/(\ldots)^2$. · **plain** The eastward speed of the energy: its sign depends on the wavelength.
  6. **did** Differentiate with respect to l · **tex** $c_{gy}=\dfrac{\partial\omega}{\partial l}=\dfrac{2\beta kl}{(k^2+l^2+F)^2}$ · **why** Only the denominator depends on l;
     the gradient $(c_{gx},c_{gy})$ of ω is at right angles to the circles of step 4 and points to their centres. · **plain**
     Energy also moves north or south unless l = 0.
  7. **did** Find where the zonal group velocity vanishes · **tex** $l=0:\quad c_{gx}=0\ \text{at}\ k=-\dfrac{f_0}{c},\qquad\omega_{max}=\dfrac{\beta c}{2f_0}$ ·
     **why** Set k² = F in step 5 (the negative root, by step 2) and insert into step 1. A zero slope of ω(k) is a
     maximum of ω. · **plain** No Rossby wave of this mode has a higher frequency; at this wavelength, 2πΛ, energy
     stands still.
  8. **did** Read the two sides of the maximum · **tex** $k^2<F:\ c_{gx}<0;\qquad k^2>F:\ c_{gx}>0;\qquad k^2\ll F:\ c_x\simeq c_{gx}\simeq-\dfrac{\beta c^2}{f_0^2}$ ·
     **why** The sign of step 5's numerator for l = 0; for very long waves both speeds tend to the same constant −βΛ²
     (no dispersion). · **plain** Long waves carry energy west, short waves east.
  9. **did** Add a uniform eastward current U · **tex** $c_x=U-\dfrac{\beta}{k^2+l^2+f_0^2/c^2}$ (13.120); $c_x=0\ \Rightarrow\ \lambda=2\pi\sqrt{U/\beta}$ · **why**
     Seen from the ground the frequency is ω + Uk (Chapter 7's $\omega_0=\omega+\mathbf U\cdot\mathbf K$ (7.9)), so the
     phase speed shifts by U; the wavelength is for F = 0, l = 0. · **plain** A westerly
     flow can hold a Rossby wave still.
- **Result:** $\omega=-\dfrac{\beta k}{k^2+l^2+f_0^2/c^2}$; $c_x=\dfrac\omega k<0$ always; (ours) $c_{gx}=\dfrac{\beta(k^2-l^2-f_0^2/c^2)}{(k^2+l^2+f_0^2/c^2)^2}$, $c_{gy}=\dfrac{2\beta kl}{(k^2+l^2+f_0^2/c^2)^2}$ — "phase always
  west; energy west for waves longer than 2πΛ, east for shorter ones".
- **Check:** units — β/k² = (1/(m s))·m² = m/s ✓. Limits — c → ∞: ω = −βk/K², and $c_{gx}=\beta(k^2-l^2)/K^4$; long-wave limit
  −βΛ². Numbers (first baroclinic mode, 35° N, Λ = 43.15 km) — 150 km wave: c_x = −0.82 cm/s, c_gx = +0.43 cm/s; 1000 km
  wave: −3.25 and −2.81 cm/s; ω_max = 4.05 × 10⁻⁷ s⁻¹ (180 days). Stationary wave — U = 17 m/s, 35° N: 5983 km. Code —
  complex-step derivatives of `GFD.rossby_omega`.
- **What it means:** the ocean adjusts to a change of wind by long Rossby waves crossing it westward (years at
  mid-latitudes, months in the tropics); atmospheric wave trains spread *downstream*, eastward; mountains anchor
  stationary waves in the westerlies. **Fails when:** within a few degrees of the equator; in strongly sheared currents
  (D24); for amplitudes large enough to advect themselves.
- **Traps:** with ω > 0 the zonal wavenumber is negative (trap T8); the "maximum" phase speed is a maximum of magnitude;
  the group-velocity arrows on the circles point inward; slip #12 — a round β and β at a stated latitude give speeds a
  fifth apart.

### D24 · The Rayleigh–Kuo criterion (13.124) from (13.122) — ★★, 8 steps, in C15 (notebook)
- **Goal:** find a condition without which a zonal current U(y) on a β-plane cannot be unstable. The book says Chapter
  11's analysis "carries over" and does not show it. · **Start:** $\Big(\dfrac{\partial}{\partial t}+\mathbf u\cdot\nabla\Big)(\zeta+f)=0$ (13.122) — "with constant depth,
  absolute vorticity is conserved following the horizontal flow". · **Plan:** • split the flow into a zonal current plus a
  small perturbation • linearise • introduce a stream function and normal modes • multiply by the complex conjugate,
  integrate, and take the imaginary part. · **Tools:** linear stability and normal modes with a complex phase speed (Ch.
  11, P214, reminders); integration by parts (P218a); real and imaginary parts of a complex equation (P260); necessary
  vs sufficient (P255). · **Assumptions:** inviscid; barotropic (no depth variation); β-plane; walls at y₁, y₂ (or decay)
  where the perturbation vanishes.
- **Steps:**
  1. **did** Decompose the flow · **tex** $u=U(y)+u',\quad v=v',\quad\zeta=\bar\zeta+\zeta',\quad\bar\zeta=-\dfrac{dU}{dy}$ · **why** A steady zonal current satisfies
     the Start exactly (nothing varies along x and v = 0); primes mark a small disturbance. · **plain** Basic flow plus
     perturbation.
  2. **did** Linearise · **tex** $\dfrac{\partial\zeta'}{\partial t}+U\dfrac{\partial\zeta'}{\partial x}+v'\dfrac{d}{dy}\big(\bar\zeta+f\big)=0$ · **why** Insert step 1 and drop products of two primed
     quantities; $\bar\zeta+f$ depends on y only, so only v′ advects it. · **plain** The perturbation vorticity is carried by
     the current and created by moving fluid across the background vorticity gradient.
  3. **did** Introduce the stream function · **tex** $\dfrac{\partial}{\partial t}\big(\nabla^2\psi\big)+U\dfrac{\partial}{\partial x}\big(\nabla^2\psi\big)+\Big(\beta-\dfrac{d^2U}{dy^2}\Big)\dfrac{\partial\psi}{\partial x}=0$ (13.123) · **why** With
     $u'=-\partial\psi/\partial y$, $v'=\partial\psi/\partial x$: $\zeta'=\nabla^2\psi$; and $d(\bar\zeta+f)/dy=-U''+\beta$. · **plain** One equation for one unknown.
  4. **did** Insert a normal mode · **tex** $(U-c)\Big[\dfrac{d^2}{dy^2}-k^2\Big]\hat\psi+\Big[\beta-\dfrac{d^2U}{dy^2}\Big]\hat\psi=0,\qquad\psi=\hat\psi(y)\,e^{ik(x-ct)}$ · **why** ∂/∂t → −ikc, ∂/∂x → ik;
     divide by ik. A complex $c=c_r+ic_i$ with $c_i>0$ means growth like $e^{kc_it}$. · **plain** Rayleigh's equation of
     Chapter 11 with β added to the vorticity gradient.
  5. **did** Divide by U − c, multiply by the conjugate, integrate · **tex**
     $\displaystyle\int_{y_1}^{y_2}\hat\psi^*\Big(\dfrac{d^2\hat\psi}{dy^2}-k^2\hat\psi\Big)dy+\displaystyle\int_{y_1}^{y_2}\dfrac{\beta-U''}{U-c}\,\lvert\hat\psi\rvert^2\,dy=0$ · **why** If $c_i\neq0$ then U − c is never zero and the division is
     allowed; multiplying by $\hat\psi^*$ turns $\hat\psi$ into the real quantity $\lvert\hat\psi\rvert^2$. · **plain** A weighted average of the
     equation over the channel.
  6. **did** Integrate the first term by parts · **tex** $-\displaystyle\int_{y_1}^{y_2}\Big(\Big\lvert\dfrac{d\hat\psi}{dy}\Big\rvert^2+k^2\lvert\hat\psi\rvert^2\Big)dy+\displaystyle\int_{y_1}^{y_2}\dfrac{\beta-U''}{U-c}\,\lvert\hat\psi\rvert^2\,dy=0$ · **why**
     Integration by parts gives $\int\hat\psi^*\hat\psi''dy=[\hat\psi^*\hat\psi']-\int\lvert\hat\psi'\rvert^2dy$, and the boundary term vanishes because $\hat\psi=0$ on the walls. · **plain**
     The first integral is real and negative.
  7. **did** Take the imaginary part · **tex** $c_i\displaystyle\int_{y_1}^{y_2}\dfrac{\beta-U''}{\lvert U-c\rvert^2}\,\lvert\hat\psi\rvert^2\,dy=0$ · **why** Write
     $\dfrac1{U-c}=\dfrac{U-c_r+ic_i}{\lvert U-c\rvert^2}$; the first integral of step 6 has no imaginary part. · **plain** Either the
     mode does not grow, or this integral is zero.
  8. **did** Draw the conclusion · **tex** $c_i\neq0\ \Rightarrow\ \dfrac{d}{dy}\big(\bar\zeta+f\big)=\beta-\dfrac{d^2U}{dy^2}$ changes sign in $(y_1,y_2)$ (13.124) · **why** The
     weight $\lvert\hat\psi\rvert^2/\lvert U-c\rvert^2$ is positive, so the integral can vanish only if β − U″ is positive in some places and
     negative in others. · **plain** Instability needs the gradient of absolute vorticity to change sign.
- **Result:** a necessary condition for barotropic instability is that $\dfrac{d}{dy}(\bar\zeta+f)=\beta-\dfrac{d^2U}{dy^2}$ changes sign within the
  flow — "Rayleigh's inflection-point criterion with β".
- **Check:** limit — β = 0: Rayleigh's criterion of Chapter 11 (U″ must change sign). Number — our easterly jet
  $U=-U_0\,\mathrm{sech}^2(y/L)$ with U₀ = 14 m/s at 12° N: the largest U″ is 2U₀/L², so the criterion is met for L < 1118 km and
  not beyond; a westerly jet of the same shape needs L < 646 km. Code — `GFD.rayleigh_kuo_criterion`.
- **What it means:** β stabilises: a current must be sharp enough for its own curvature U″ to beat β somewhere. The
  quantity that must change sign is the cross-stream gradient of the basic state's potential vorticity — the same
  gradient that Rossby waves propagate on. **Fails when (as a test):** it is necessary, not sufficient: a flow that
  passes may still be stable.
- **Traps:** necessary, not sufficient; the β term comes from v′ d(ζ̄ + f)/dy; the stream-function sign differs from
  Chapter 11's (trap T14) — the eigenvalue c is unaffected; the division in step 5 needs c_i ≠ 0.

### D25 · The Eady basic state (13.128) and perturbation equation (13.136) — ★★★, 13 steps, in C16 (notebook · `eady_instability`)
- **Goal:** derive the equation governing small quasi-geostrophic disturbances on a wind that increases linearly with
  height in a uniformly stratified, rotating layer. · **Start:** the f-plane, hydrostatic, inviscid set (13.125):
  $\dfrac{\partial u}{\partial t}+u\dfrac{\partial u}{\partial x}+v\dfrac{\partial u}{\partial y}-fv=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial x}$, $\dfrac{\partial v}{\partial t}+u\dfrac{\partial v}{\partial x}+v\dfrac{\partial v}{\partial y}+fu=-\dfrac1{\rho_0}\dfrac{\partial p}{\partial y}$, $0=-\dfrac{\partial p}{\partial z}-\rho g$, $\dfrac{\partial u}{\partial x}+\dfrac{\partial v}{\partial y}+\dfrac{\partial w}{\partial z}=0$,
  $\dfrac{\partial\rho}{\partial t}+u\dfrac{\partial\rho}{\partial x}+v\dfrac{\partial\rho}{\partial y}+w\dfrac{\partial\rho}{\partial z}=0$ — "momentum, hydrostatics, mass, density; no friction, f constant". · **Plan:** • split into a
  basic state and a perturbation • find the basic wind from the thermal wind • linearise the vorticity equation • express
  ζ′ and w′ through the pressure perturbation • substitute. · **Tools:** ordering in a small parameter (P328);
  cross-differentiation (P177); thermal wind (C03, D04); linearisation. · **Assumptions:** quasi-geostrophic (Ro ≪ 1);
  f-plane; uniform N; hydrostatic; constant ∂ρ̄/∂y; vertical advection of momentum neglected (trap T15).
- **Steps:**
  1. **did** Decompose · **tex** $u=U(z)+u',\ \ v=v',\ \ w=w',\ \ \rho=\bar\rho(y,z)+\rho',\ \ p=\bar p(y,z)+p'$ (13.126) · **why** A zonal wind that
     varies only with height, in balance with a density field that varies northward and upward; primes are small. ·
     **plain** Basic state plus perturbation.
  2. **did** Write the balance of the basic state · **tex** $fU=-\dfrac1{\rho_0}\dfrac{\partial\bar p}{\partial y},\qquad0=-\dfrac{\partial\bar p}{\partial z}-\bar\rho g$ (13.127) · **why** Insert the unprimed
     fields into the Start: with no perturbation only the geostrophic and hydrostatic balances remain. · **plain** The
     basic wind is geostrophic and hydrostatic.
  3. **did** Eliminate the basic pressure · **tex** $\dfrac{dU}{dz}=\dfrac{g}{f\rho_0}\dfrac{\partial\bar\rho}{\partial y}$ (13.128), so $U=\dfrac{U_0z}{H}$ · **why** ∂/∂z of the first and ∂/∂y of
     the second share the mixed derivative of p̄ (the thermal-wind move of D04); a constant density gradient and U(0) = 0
     give a linear wind. · **plain** Denser toward the pole means a westerly wind growing with height.
  4. **did** Form the vorticity equation of the total flow · **tex** $\dfrac{\partial\zeta}{\partial t}+u\dfrac{\partial\zeta}{\partial x}+v\dfrac{\partial\zeta}{\partial y}-(\zeta+f)\dfrac{\partial w}{\partial z}=0$ (13.129) · **why**
     Cross-differentiate the two momentum equations as in D17 steps 1–6 (now β = 0) and replace the horizontal
     divergence by −∂w/∂z from continuity. · **plain** Vorticity is advected and created by vertical stretching.
  5. **did** Linearise · **tex** $\dfrac{\partial\zeta'}{\partial t}+U\dfrac{\partial\zeta'}{\partial x}-f\dfrac{\partial w'}{\partial z}=0$ (13.130) · **why** The basic flow has no vertical vorticity
     (U depends on z only); drop products of primed quantities and ζ′ against f. · **plain** Stretching of the planetary
     vorticity is the only source.
  6. **did** Use geostrophic perturbation velocities · **tex** $u'\simeq-\dfrac1{\rho_0f}\dfrac{\partial p'}{\partial y},\qquad v'\simeq\dfrac1{\rho_0f}\dfrac{\partial p'}{\partial x}$ (13.131) · **why** To lowest
     order in the Rossby number the perturbation is geostrophic too (D03). · **plain** The disturbance winds follow the
     disturbance isobars.
  7. **did** Compute the perturbation vorticity · **tex** $\zeta'=\dfrac{\partial v'}{\partial x}-\dfrac{\partial u'}{\partial y}=\dfrac1{\rho_0f}\nabla_H^2p'$ (13.132) · **why** Differentiate step 6;
     ρ₀ and f are constants. · **plain** A low is a cyclonic vortex.
  8. **did** Linearise the density equation · **tex** $\dfrac{\partial\rho'}{\partial t}+U\dfrac{\partial\rho'}{\partial x}+v'\dfrac{\partial\bar\rho}{\partial y}-\dfrac{\rho_0N^2w'}{g}=0$ (13.133) · **why** Keep terms with one
     primed factor; write $\partial\bar\rho/\partial z=-\rho_0N^2/g$. The term $v'\,\partial\bar\rho/\partial y$ is where the basic state's tilt enters. · **plain**
     Density changes by advection along the flow, by north–south motion across the basic gradient, and by vertical
     motion.
  9. **did** Use hydrostatic balance of the perturbation · **tex** $0=-\dfrac{\partial p'}{\partial z}-\rho'g\ \Rightarrow\ \rho'=-\dfrac1g\dfrac{\partial p'}{\partial z}$ (13.134) · **why** Subtract
     the basic-state balance of step 2 from the hydrostatic equation of the total flow. · **plain** Density anomalies
     are vertical gradients of the pressure anomaly.
  10. **did** Solve the density equation for w′ · **tex** $w'=-\dfrac1{\rho_0N^2}\Big[\Big(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\Big)\dfrac{\partial p'}{\partial z}-\dfrac{dU}{dz}\dfrac{\partial p'}{\partial x}\Big]$ (13.135) · **why** Insert
      step 9, step 6 for v′, and step 3 in the form $\partial\bar\rho/\partial y=(f\rho_0/g)\,dU/dz$; multiply by $g/(\rho_0N^2)$. · **plain** Rising
      motion is diagnosed from the pressure field.
  11. **did** Differentiate w′ in z · **tex**
      $\dfrac{\partial w'}{\partial z}=-\dfrac1{\rho_0N^2}\Big[\Big(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\Big)\dfrac{\partial^2p'}{\partial z^2}+\dfrac{dU}{dz}\dfrac{\partial^2p'}{\partial x\,\partial z}-\dfrac{dU}{dz}\dfrac{\partial^2p'}{\partial x\,\partial z}-\dfrac{d^2U}{dz^2}\dfrac{\partial p'}{\partial x}\Big]$ · **why** Product rule: the z-derivative also acts on U(z)
      inside the advection operator and on dU/dz in the last term; N is uniform. · **plain** Four terms, two of them
      equal and opposite.
  12. **did** Cancel and use the linear wind · **tex** $\dfrac{\partial w'}{\partial z}=-\dfrac1{\rho_0N^2}\Big(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\Big)\dfrac{\partial^2p'}{\partial z^2}$ · **why** The two shear terms cancel
      exactly; the last vanishes because $d^2U/dz^2=0$ for $U=U_0z/H$. · **plain** The stretching depends only on the curvature
      of the pressure in z.
  13. **did** Substitute into the vorticity equation · **tex** $\Big(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\Big)\Big[\nabla_H^2p'+\dfrac{f^2}{N^2}\dfrac{\partial^2p'}{\partial z^2}\Big]=0$ (13.136) · **why** Insert
      steps 7 and 12 into step 5 and multiply by ρ₀f. · **plain** The bracket — the perturbation's potential vorticity —
      is carried along by the basic wind.
- **Result:** $\Big(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\Big)\Big[\nabla_H^2p'+\dfrac{f^2}{N^2}\dfrac{\partial^2p'}{\partial z^2}\Big]=0$ — "quasi-geostrophic perturbations on an eastward current U(z)".
- **Check:** units — ∇²p′ and (f²/N²) ∂²p′/∂z² are both Pa/m² ✓. Limit — U = 0: the bracket is steady. Structure — the
  bracket is ρ₀f times (relative vorticity + stretching), the three-dimensional cousin of D22's bracket.
- **sympy check** (`check_src`; executed as written — every assertion passes and the final line is printed):
  ```python
  import sympy as sp                                            # symbolic algebra
  x, y, z, t = sp.symbols("x y z t", real=True)                 # east, north, up, time
  f, N, rho0, g, U0, H = sp.symbols("f N rho_0 g U_0 H", positive=True)   # constants of the Eady problem
  p = sp.Function("p")(x, y, z, t)                              # perturbation pressure p'
  U = U0 * z / H                                                # step 3: the basic flow, linear in z
  adv = lambda F: sp.diff(F, t) + U * sp.diff(F, x)             # d/dt + U d/dx (following the basic flow)
  drho_dy = f * rho0 / g * sp.diff(U, z)                        # step 3: thermal wind solved for d(rho-bar)/dy
  vp = sp.diff(p, x) / (rho0 * f)                               # step 6: geostrophic v'
  up = -sp.diff(p, y) / (rho0 * f)                              # step 6: geostrophic u'
  zeta = sp.diff(vp, x) - sp.diff(up, y)                        # step 7: perturbation vorticity
  assert sp.simplify(zeta - (sp.diff(p, x, 2) + sp.diff(p, y, 2)) / (rho0 * f)) == 0   # (13.132)
  rhop = -sp.diff(p, z) / g                                     # step 9: hydrostatic perturbation
  wsym = sp.symbols("w")                                        # unknown w' in the density equation
  dens = adv(rhop) + vp * drho_dy - rho0 * N**2 * wsym / g      # step 8: (13.133)
  wp = sp.solve(dens, wsym)[0]                                  # step 10: solve it for w'
  book_w = -(adv(sp.diff(p, z)) - sp.diff(U, z) * sp.diff(p, x)) / (rho0 * N**2)   # (13.135)
  assert sp.simplify(wp - book_w) == 0
  vort = adv(zeta) - f * sp.diff(wp, z)                         # step 5: (13.130) with zeta' and w' put in
  lap = sp.diff(p, x, 2) + sp.diff(p, y, 2)                     # horizontal Laplacian of p'
  book = adv(lap + f**2 / N**2 * sp.diff(p, z, 2))              # (13.136)
  assert sp.simplify(sp.expand(vort * rho0 * f - book)) == 0    # steps 11-13: the two shear terms have cancelled
  print("D25 checks passed: (13.132), (13.135) and the Eady equation (13.136)")
  ```
- **What it means:** in the interior nothing can create perturbation potential vorticity, so a normal mode has none there
  — the whole instability lives in what happens at the two lids (D26). **Fails when:** β matters (the Charney problem,
  beyond the book), N or the shear vary with height (interior potential-vorticity gradients appear), or Ro is not small.
- **Traps:** the density equation keeps $v'\,\partial\bar\rho/\partial y$ — that is where the shear enters; when ∂/∂z acts on w′ the two shear
  terms cancel; vertical advection of momentum is neglected silently (trap T15); α in the next block is a scaled
  wavenumber, not the thermal expansion coefficient (trap T13).


### D26 · The Eady boundary conditions, the 2 × 2 system and the phase speed (13.141) — ★★★, 12 steps, in C16 (notebook · `eady_instability`)
- **Goal:** solve the Eady equation between two rigid lids and find the phase speed of its waves — and with it, when
  they grow. The book says "after some straightforward algebra". · **Start:** $\Big(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\Big)\Big[\nabla_H^2p'+\dfrac{f^2}{N^2}\dfrac{\partial^2p'}{\partial z^2}\Big]=0$ (13.136) with
  $U=U_0z/H$ and $w'=0$ at z = 0 and z = H — "the Result of D25; no flow through the lids". · **Plan:** • insert a wave • solve
  the vertical structure with hyperbolic functions centred at mid-depth • write w′ = 0 on each lid in terms of the
  pressure • demand a non-trivial solution of the resulting 2 × 2 system • reduce the determinant to a quadratic in c. ·
  **Tools:** normal modes (Ch. 11, reminder); non-trivial solutions of a homogeneous 2 × 2 system (P330); hyperbolic
  identities and coth (P331); the quadratic formula with a possibly negative discriminant (P159, reminder). ·
  **Assumptions:** rigid horizontal lids; unbounded in x and y; U ≠ c in the interior (step 2).
- **Steps:**
  1. **did** Insert a wave · **tex** $p'=\hat p(z)\,e^{i(kx+ly-\omega t)}$ (13.137) · **why** The coefficients depend on z only, so the solution may
     be taken as a plane wave in x, y, t with an unknown vertical structure. · **plain** A wave whose strength varies with
     height.
  2. **did** Substitute into the Start · **tex** $\dfrac{d^2\hat p}{dz^2}-\alpha^2\hat p=0$ (13.138) · **why** The operator gives the factor $ik(U-c)$ with $c=\omega/k$,
     non-zero in the interior; the bracket gives $-(k^2+l^2)\hat p+(f^2/N^2)\hat p''$. · **plain** Between the lids the perturbation's
     potential vorticity is zero.
  3. **did** Name the scaled wavenumber · **tex** $\alpha^2\equiv\dfrac{N^2}{f^2}\big(k^2+l^2\big)$ (13.139) · **why** A definition; 1/α is the height over which a
     disturbance of horizontal wavenumber K is felt, and αH = K × (NH/f). · **plain** A wave's vertical reach is its
     horizontal scale times f/N.
  4. **did** Write the solution about mid-depth · **tex** $\hat p=A\cosh\alpha\Big(z-\dfrac H2\Big)+B\sinh\alpha\Big(z-\dfrac H2\Big)$ (13.140) · **why** cosh and sinh are two
     independent solutions; centring them at H/2 makes the two lids mirror images and shortens the algebra. · **plain**
     An even part and an odd part about the middle of the layer.
  5. **did** Express w′ = 0 through the pressure · **tex** $\Big(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\Big)\dfrac{\partial p'}{\partial z}-\dfrac{U_0}{H}\dfrac{\partial p'}{\partial x}=0$ at $z=0$ and $z=H$ · **why** Set the
     bracket of D25's $w'=-\dfrac1{\rho_0N^2}\Big[\Big(\dfrac{\partial}{\partial t}+U\dfrac{\partial}{\partial x}\Big)\dfrac{\partial p'}{\partial z}-\dfrac{dU}{dz}\dfrac{\partial p'}{\partial x}\Big]$ (13.135) to zero; U = 0 at the lower lid, U₀ at the
     upper. · **plain** On each lid, advection of the temperature anomaly balances advection of the basic temperature.
  6. **did** Insert the wave · **tex** $(U-c)\dfrac{d\hat p}{dz}-\dfrac{U_0}{H}\hat p=0$ at $z=0\ (U=0)$ and $z=H\ (U=U_0)$ · **why** ∂/∂t → −ikc and ∂/∂x → ik;
     divide by ik. · **plain** Two conditions on the vertical structure.
  7. **did** Evaluate the structure on the lids · **tex** $\hat p(0)=A\cosh X-B\sinh X,\quad\hat p'(0)=\alpha(-A\sinh X+B\cosh X),\quad X\equiv\dfrac{\alpha H}2$ · **why** At
     z = 0 the argument is −X; cosh is even and sinh odd. At z = H the argument is +X, so the same expressions hold
     with +B sinh X and +A sinh X. · **plain** The even and odd parts add on one lid and subtract on the other.
  8. **did** Write the two conditions · **tex** $A\Big[\alpha c\sinh X-\dfrac{U_0}H\cosh X\Big]+B\Big[-\alpha c\cosh X+\dfrac{U_0}H\sinh X\Big]=0$;
     $A\Big[\alpha(U_0-c)\sinh X-\dfrac{U_0}H\cosh X\Big]+B\Big[\alpha(U_0-c)\cosh X-\dfrac{U_0}H\sinh X\Big]=0$ · **why** Insert step 7 into step 6 and collect the
     terms with A and with B. These are the book's pair. · **plain** Two homogeneous equations for A and B.
  9. **did** Set the determinant to zero and expand · **tex** $\alpha^2c\,(U_0-c)\sinh\alpha H-\alpha\dfrac{U_0^2}{H}\cosh\alpha H+\dfrac{U_0^2}{H^2}\sinh\alpha H=0$ · **why** A
     non-trivial (A, B) needs a vanishing determinant; after multiplying out, use $2\sinh X\cosh X=\sinh2X$ and
     $\cosh^2X+\sinh^2X=\cosh2X$ with 2X = αH. · **plain** One equation for the phase speed.
  10. **did** Divide by α² sinh αH · **tex** $c^2-U_0c+\dfrac{U_0^2}{\alpha H}\coth\alpha H-\dfrac{U_0^2}{\alpha^2H^2}=0$ · **why** sinh αH > 0 for αH > 0; coth = cosh/sinh. ·
      **plain** A quadratic in c.
  11. **did** Solve the quadratic · **tex** $c=\dfrac{U_0}2\pm\dfrac{U_0}{\alpha H}\sqrt{\dfrac{\alpha^2H^2}4-\alpha H\coth\alpha H+1}$ · **why** Quadratic formula; take the factor
      $U_0/(\alpha H)$ out of the root. If the radicand is negative the two roots are complex conjugates. · **plain** Both waves
      move, on average, at the mid-level wind speed.
  12. **did** Factor the radicand · **tex** $c=\dfrac{U_0}2\pm\dfrac{U_0}{\alpha H}\sqrt{\Big(\dfrac{\alpha H}2-\tanh\dfrac{\alpha H}2\Big)\Big(\dfrac{\alpha H}2-\coth\dfrac{\alpha H}2\Big)}$ (13.141) · **why** Multiply the two brackets out:
      $(X-\tanh X)(X-\coth X)=X^2-X(\tanh X+\coth X)+1$; then use the identity $\tanh X+\coth X=2\coth2X$ of the primer. · **plain** The sign of a product of
      two simple factors decides stability.
- **Result:** $c=\dfrac{U_0}2\pm\dfrac{U_0}{\alpha H}\sqrt{\Big(\dfrac{\alpha H}2-\tanh\dfrac{\alpha H}2\Big)\Big(\dfrac{\alpha H}2-\coth\dfrac{\alpha H}2\Big)}$ — "the Eady phase speed: complex, hence growth, when the second
  factor is negative".
- **Check:** units — c in m/s ✓ (αH dimensionless). Limits — αH → ∞: both factors → αH/2, c → U₀/2 ± U₀/2, i.e. c → 0 and c →
  U₀: two separate edge waves, each riding with the wind at its own lid. Numbers — αH = 1: c = (0.5 ± 0.2511 i) U₀; αH = 3:
  c = 0.3384 U₀ and 0.6616 U₀. Independent route — a Chebyshev eigen-solve (`ch13.eady_numeric_eigs`).
- **sympy check** (`check_src`; executed as written — every assertion passes; it prints c = (0.5 + 0.215819 i) U₀ at αH = 1.4):
  ```python
  import sympy as sp                                            # symbolic algebra
  z = sp.symbols("z", real=True)                                # height, 0 <= z <= H
  al, H, U0 = sp.symbols("alpha H U_0", positive=True)          # scaled wavenumber, depth, top speed
  c, A, B = sp.symbols("c A B")                                 # phase speed (may be complex) and the two constants
  phat = A * sp.cosh(al * (z - H / 2)) + B * sp.sinh(al * (z - H / 2))   # step 4: (13.140)
  bc = lambda zz: ((U0 * zz / H - c) * sp.diff(phat, z) - U0 / H * phat).subs(z, zz)   # steps 5-6: w' = 0 at z = zz
  rows = [sp.expand(bc(0)), sp.expand(bc(H))]                   # step 7: the two boundary conditions
  M = sp.Matrix([[r.coeff(A), r.coeff(B)] for r in rows])       # step 8: coefficients of A and B
  X = al * H / 2                                                # the half-argument alpha H / 2
  det = sp.simplify(M.det().rewrite(sp.exp))                    # step 9: the determinant (as exponentials)
  quad = c**2 - U0 * c + U0**2 / (al * H) * sp.coth(al * H) - U0**2 / (al * H)**2   # step 10: the quadratic
  ratio = sp.simplify((det / quad.rewrite(sp.exp)))             # they must differ by a factor without c
  assert c not in ratio.free_symbols
  rad = (X - sp.tanh(X)) * (X - sp.coth(X))                     # step 12: the radicand of (13.141)
  c_plus = U0 / 2 + U0 / (al * H) * sp.sqrt(rad)                # (13.141), upper sign
  val = quad.subs(c, c_plus).subs({al: sp.Rational(7, 5), H: 1, U0: 1})   # a wavenumber inside the unstable band
  assert abs(sp.N(val, 30)) < 1e-25                             # steps 11-12: (13.141) solves the quadratic
  val2 = quad.subs(c, c_plus).subs({al: 3, H: 2, U0: 5})        # and one beyond the cut-off
  assert abs(sp.N(val2, 30)) < 1e-25
  print("D26 checks passed: determinant -> quadratic in c -> (13.141); c at alpha*H = 1.4:",
        sp.N(c_plus.subs({al: sp.Rational(7, 5), H: 1, U0: 1}), 6))
  ```
- **What it means:** the first factor is always positive (x > tanh x); the second is negative for long waves (x < coth x).
  So every wave longer than a cut-off is unstable: a sheared, stratified, rotating flow has no stable range at long
  wavelengths. Whenever c is complex, one of the conjugate pair grows. **Fails when:** the lids are not rigid (a
  tropopause, a free surface), β is included, or the shear is not uniform.
- **Traps:** the lid condition is w′ = 0 written in terms of the pressure — at z = H the extra term with U₀ appears;
  solving about mid-depth is what makes the determinant a quadratic in c − U₀/2; tanh and coth must not be swapped;
  growth needs the product under the root to be negative.

### D27 · The Eady cut-off, the unstable band (13.142), the fastest wave and the maximum growth rate — ★★, 6 steps, in C16 (notebook)
- **Goal:** turn the phase-speed formula into the numbers a climate scientist uses: which wavelengths grow, which grows
  fastest, and how fast in days. The book says "it can be shown" and stops at the wavelength; the growth rate in
  physical units is **ours, computed** (no citation). · **Start:** $c=\dfrac{U_0}2\pm\dfrac{U_0}{\alpha H}\sqrt{\Big(\dfrac{\alpha H}2-\tanh\dfrac{\alpha H}2\Big)\Big(\dfrac{\alpha H}2-\coth\dfrac{\alpha H}2\Big)}$ (13.141) —
  "the Result of D26". · **Plan:** • find when the radicand is negative • solve the marginal condition numerically •
  write the growth rate k c_i and scale it • maximise it numerically • put units back. · **Tools:** graphical roots (P316);
  `brentq` (P108, reminder); `minimize_scalar` (P170, reminder). · **Assumptions:** l = 0 for the wavelengths and the
  maximum; x ≡ αH/2.
- **Steps:**
  1. **did** Find the sign of each factor · **tex** $x-\tanh x>0\ \text{for all}\ x>0;\qquad x-\coth x<0\iff x<\coth x$ · **why** tanh x < x
     for every positive x, while coth x falls from +∞ toward 1 and so crosses the rising line x exactly once. · **plain**
     Only the second factor can change sign.
  2. **did** Solve the marginal condition · **tex** $\dfrac{\alpha_cH}2=\coth\dfrac{\alpha_cH}2\ \Rightarrow\ \alpha_cH=2.3994$ · **why** The radicand is zero where x = coth x;
     this transcendental equation is solved by `brentq` (x_c = 1.19968). Ours, computed to five digits. · **plain**
     Waves with αH below 2.3994 are unstable.
  3. **did** Translate into wavelength · **tex** $\dfrac{HN}{f}<\dfrac{\alpha_cH}{k}$ (13.142; the book prints $\alpha_cH$ as a rounded number), i.e.
     $\lambda>\dfrac{2\pi}{\alpha_cH}\Lambda=2.6187\,\Lambda,\quad\Lambda\equiv\dfrac{HN}f$ · **why** For l = 0, α = Nk/f, so αH = kΛ; αH < α_cH is the printed inequality, and
     λ = 2π/k. · **plain** Only waves longer than 2.62 Eady radii grow.
  4. **did** Write the growth rate · **tex** $\sigma=kc_i=\dfrac{fU_0}{NH}\,G(x),\qquad G(x)=\sqrt{(x-\tanh x)(\coth x-x)}$ · **why** A wave
     $e^{ik(x-ct)}$ grows like $e^{kc_it}$; from the Start $c_i=(U_0/\alpha H)\,G$, and for l = 0, k/α = f/N. · **plain** The growth
     rate is a universal curve G times fU₀/(NH).
  5. **did** Maximise G · **tex** $G_{max}=0.30982\ \text{at}\ x=0.80306,\ \text{i.e.}\ \alpha H=1.6061,\quad\lambda=\dfrac{2\pi}{1.6061}\Lambda=3.9120\,\Lambda$ · **why** G is zero
     at x = 0 and at x_c and positive between; a bounded scalar maximiser finds the peak. Ours, computed. · **plain**
     The fastest-growing wave is 3.91 Eady radii long.
  6. **did** Put the units back · **tex** $\sigma_{max}=0.30982\,\dfrac{f}{N}\dfrac{dU}{dz},\qquad\text{e-folding time}=\dfrac1{\sigma_{max}}$ · **why** U₀/H is the shear dU/dz;
     the result no longer contains the depth H. · **plain** Storms grow faster where the shear is strong, the
     stratification weak and the latitude high.
- **Result:** unstable for $\alpha H<\alpha_cH=2.3994$ (wavelengths above 2.6187 Λ); fastest at αH = 1.6061 (3.9120 Λ) with (ours)
  $\sigma_{max}=0.30982\,\dfrac fN\dfrac{dU}{dz}$ — "the Eady growth rate".
- **Check:** units — (s⁻¹/s⁻¹)(s⁻¹) = s⁻¹ ✓. Limits — x → 0: G ≈ x/√3 → 0 (long waves grow slowly, although c_i is largest
  there); x → x_c: G → 0. Numbers — atmosphere (35° N, N = 1.1 × 10⁻² s⁻¹, H = 9 km, U₀ = 27 m/s): Λ = 1183 km, e-folding 1.64
  days, fastest wavelength 4630 km, cut-off 3099 km; ocean (N = 5 × 10⁻³ s⁻¹, H = 1 km, U₀ = 0.1 m/s): 60 km, 22 days, 234 km.
- **What it means:** the atmosphere's storms and the ocean's mesoscale eddies are the same instability at two Rossby
  radii; "Eady growth rate" maps are this formula evaluated from analysed winds and temperatures. **Fails when:** the
  real profile departs from uniform shear and N; β and friction lower the growth rate and remove the long-wave
  instability in practice.
- **Traps:** the growth rate is k c_i, not c_i; this Λ = NH/f has no π (trap T11); the constants are ours, computed —
  quote them to five digits and never as the book's rounded values; l ≠ 0 reduces the growth rate by the factor k/K.

### D28 · Fjørtoft's ratios (13.145) — ★, 5 steps, in C17 (notebook)
- **Goal:** if nonlinear interactions move the energy at one wavenumber to one smaller and one larger wavenumber while
  conserving both energy and enstrophy, find how it is shared. · **Start:** $S_0=S_1+S_2$ and $K_0^2S_0=K_1^2S_1+K_2^2S_2$ with
  $K_1<K_0<K_2$ — "energy conserved; enstrophy (K² × energy) conserved", the three-mode form of
  $\dfrac d{dt}\displaystyle\int_0^\infty S(K)\,dK=0$ and $\dfrac d{dt}\displaystyle\int_0^\infty K^2S(K)\,dK=0$ (13.143)–(13.144). · **Plan:** • treat the two sums as two linear equations for S₁ and
  S₂ • eliminate one unknown • form the ratios. · **Tools:** two linear equations in two unknowns (P57, reminder);
  difference of two squares. · **Assumptions:** all the energy goes to exactly two wavenumbers; inviscid; two-dimensional.
- **Steps:**
  1. **did** Write the two constraints as a system · **tex** $S_1+S_2=S_0,\qquad K_1^2S_1+K_2^2S_2=K_0^2S_0$ · **why** S₀ and the three wavenumbers
     are given; S₁ and S₂ are the unknowns. Two equations, two unknowns: no freedom is left. · **plain** The split is
     completely determined by the two conservation laws.
  2. **did** Eliminate S₁ · **tex** $\big(K_2^2-K_1^2\big)S_2=\big(K_0^2-K_1^2\big)S_0$ · **why** Multiply the first equation by K₁² and subtract it from
     the second. · **plain** The share that goes to the small scale.
  3. **did** Eliminate S₂ · **tex** $\big(K_2^2-K_1^2\big)S_1=\big(K_2^2-K_0^2\big)S_0$ · **why** Multiply the first by K₂² and subtract the second. ·
     **plain** The share that goes to the large scale.
  4. **did** Divide · **tex** $\dfrac{S_1}{S_2}=\dfrac{K_2^2-K_0^2}{K_0^2-K_1^2}=\dfrac{K_2-K_0}{K_0-K_1}\,\dfrac{K_2+K_0}{K_1+K_0}$ (13.145), first member · **why** Step 3 over step 2; then
     factor each difference of squares. · **plain** The energy ratio is larger than the ratio of the wavenumber gaps,
     because K₂ + K₀ > K₁ + K₀.
  5. **did** Multiply by K₁²/K₂² · **tex** $\dfrac{K_1^2S_1}{K_2^2S_2}=\dfrac{K_1^2}{K_2^2}\,\dfrac{K_2^2-K_0^2}{K_0^2-K_1^2}$ (13.145), second member · **why** The enstrophy at a wavenumber is K²
     times its energy. · **plain** The enstrophy ratio: for comparable gaps it is smaller than 1.
- **Result:** $\dfrac{S_1}{S_2}=\dfrac{K_2-K_0}{K_0-K_1}\,\dfrac{K_2+K_0}{K_1+K_0}$ and $\dfrac{K_1^2S_1}{K_2^2S_2}=\dfrac{K_1^2}{K_2^2}\,\dfrac{K_2^2-K_0^2}{K_0^2-K_1^2}$ — "more energy to the larger scale, more enstrophy to
  the smaller".
- **Check:** units — both ratios are pure numbers ✓. Number — K₁ = K₀/3, K₂ = 3K₀/2: S₁/S₂ = 1.40625, enstrophy ratio
  0.06944; S₁ = 0.5844 S₀, S₂ = 0.4156 S₀. Limit — K₁ → K₀: all the energy stays at K₀ (S₂ → 0). Code — `np.linalg.solve`.
- **What it means:** in nearly two-dimensional (geostrophic) turbulence energy cascades to large scales and enstrophy to
  small scales — the reverse of Chapter 12's three-dimensional cascade; eddies merge and grow. **Fails when:** the flow
  is three-dimensional (vortex stretching destroys enstrophy conservation) or strongly dissipative.
- **Traps:** all the energy is assumed to go to just two wavenumbers; the second ratio is the inverse statement, not an
  independent one; the spectrum here is one-sided with no factor ½ (trap T16).

### D29 · The enstrophy spectrum K²S(K), the −3 range and the Rhines length — ★★, 7 steps, in C17 (notebook)
- **Goal:** justify three statements the book makes in passing: that enstrophy is distributed as K²S(K); that the
  enstrophy-cascading range has a spectrum proportional to K⁻³; and that eddies stop growing at a length √(u/β). All
  three are presented as **ours, computed** — no citation. · **Start:** $\overline{u^2}=\displaystyle\int_0^\infty S(K)\,dK$ and $\overline{\zeta^2}=\displaystyle\int_0^\infty K^2S(K)\,dK$ (the book's
  definitions in §13.18), and the barotropic vorticity equation $\Big(\dfrac{\partial}{\partial t}+\mathbf u\cdot\nabla\Big)(\zeta+f)=0$ (13.122). · **Plan:** • get the K²
  weighting from the vorticity of one Fourier mode • list the dimensions • find the exponents of each inertial range by
  dimensional analysis • compare two terms of the vorticity equation for the Rhines length. · **Tools:** derivative of a
  Fourier mode (P142, reminder); the Π theorem (Ch. 1, reminder); dominant balance (P198, reminder); enstrophy (P333). ·
  **Assumptions:** isotropic, two-dimensional, non-divergent; an inertial range in which only one flux matters.
- **Steps:**
  1. **did** Take the vorticity of one Fourier mode · **tex** $(u,v)=(\hat u,\hat v)\,e^{i(kx+ly)}\ \Rightarrow\ \hat\zeta=i\big(k\hat v-l\hat u\big)$ · **why** By definition $\zeta=\partial v/\partial x-\partial u/\partial y$,
     and each derivative of a Fourier mode is a multiplication by ik or il. · **plain** In wavenumber space taking a curl
     is algebra.
  2. **did** Use non-divergence · **tex** $k\hat u+l\hat v=0\ \Rightarrow\ \lvert\hat\zeta\rvert^2=\big(k^2+l^2\big)\big(\lvert\hat u\rvert^2+\lvert\hat v\rvert^2\big)=K^2\big(\lvert\hat u\rvert^2+\lvert\hat v\rvert^2\big)$ · **why** A non-divergent
     velocity is perpendicular to its wavevector, so the cross product in step 1 has the full length K times the
     speed. · **plain** Each mode's enstrophy is K² times its energy: the enstrophy spectrum is K²S(K).
  3. **did** List the dimensions · **tex** $[S]=\mathrm{m^3\,s^{-2}},\quad[K]=\mathrm{m^{-1}},\quad[\varepsilon]=\mathrm{m^2\,s^{-3}},\quad[\alpha]=\mathrm{s^{-3}}$ · **why** S dK is a velocity squared;
     ε is energy per mass per time; α (the enstrophy flux in this section) is vorticity squared per time. · **plain**
     The enstrophy flux contains no length at all.
  4. **did** Find the exponents of the enstrophy range · **tex** $S=C_Z\,\alpha^a K^b:\quad-2=-3a,\quad3=-b\ \Rightarrow\ S\propto\alpha^{2/3}K^{-3}$ · **why** If the
     spectrum there depends only on α and K, both sides must have the same powers of seconds and metres (one
     dimensionless group). · **plain** Toward small scales the spectrum falls as K⁻³.
  5. **did** Do the same for the energy range · **tex** $S=C_E\,\varepsilon^a K^b:\quad-2=-3a,\quad3=2a-b\ \Rightarrow\ S\propto\varepsilon^{2/3}K^{-5/3}$ · **why** The same
     argument with ε in place of α; it is Kolmogorov's argument of Chapter 12, although here the energy flux is toward
     large scales. · **plain** Toward large scales the spectrum has the familiar −5/3 slope.
  6. **did** Compare two terms of the vorticity equation · **tex** $\mathbf u\cdot\nabla\zeta\sim\dfrac{u^2}{l^2},\qquad\beta v\sim\beta u$ · **why** For eddies of speed u
     and size l, ζ ~ u/l and its gradient ~ u/l²; the term βv comes from $\mathbf u\cdot\nabla f=\beta v$. · **plain** Eddies advect each
     other's vorticity; β makes them radiate Rossby waves.
  7. **did** Equate them · **tex** $\dfrac{u^2}{l^2}=\beta u\ \Rightarrow\ l\sim\sqrt{\dfrac u\beta}$ · **why** Small eddies are dominated by advection (turbulence),
     large ones by β (waves); the crossover is where the two terms are equal. · **plain** The inverse cascade stalls at
     this size, and the flow organises into east–west jets of about this width.
- **Result:** enstrophy spectrum $K^2S(K)$; $S\propto\alpha^{2/3}K^{-3}$ in the enstrophy range and $S\propto\varepsilon^{2/3}K^{-5/3}$ in the energy range;
  $l\sim\sqrt{u/\beta}$ — "two inertial ranges and a stopping scale".
- **Check:** units — $\alpha^{2/3}K^{-3}$ = s⁻² m³ ✓; $\varepsilon^{2/3}K^{-5/3}$ = m^{4/3} s⁻² m^{5/3} = m³ s⁻² ✓; √(u/β) = √((m/s)(m s)) = m ✓. Code —
  `DIM.pi_groups` returns exponents (2/3, −3). Number — 35° N: √(u/β) = 833 km for u = 13 m/s and 65 km for u = 0.08 m/s.
- **What it means:** energy injected at the scale of baroclinic eddies moves up-scale until β turns eddies into waves:
  the alternating jets of the giant planets and of the ocean are of this width. **Fails when:** there is no clean
  inertial range (our 64² demo has none: the K⁻³ range is not resolved there and no slope is quoted from it), or dissipation and forcing overlap the
  range.
- **Traps:** α is the enstrophy flux [s⁻³] in this section (trap T13); the exponents are exact rationals from the Π
  theorem, the constants C_E, C_Z are not given by it; nothing here is a benchmark — it is dimensional reasoning.


---

## Closing check — design gate (bullets only, so that the ledger parser of Part E is not disturbed)

Counted by a scratch script kept outside the repository (it parses this file the way `tools/nbkit.py`, `tools/coverage_check.py`
and `tools/embed_check.py` do), on the file as saved on 2026-10-07.

- **CORE blocks:** 17 of 17 (`#### C01 —` … `#### C17 —`, each with one `nb.core` row), every one with plain words, the idea,
  a worked example with easy numbers, a code row with *explain*, a from-scratch `nb.check_agree` row, at least one visual
  with see / read / change, and a "What would change if…" link. Blocks without an explainer of their own (C01, C04, C07,
  C09, C13, C14, C17) each have at least one Python visual of their own (listed in A.19).
- **Sections:** 18 of 18 `nb.section` rows (13.1–13.18); §13.1–§13.3 open the C01 story before `nb.core("C01", …)`; §13.16
  continues block C15 through the re-entry rule.
- **Curation ids:** 218 of 218 placed (17 C · 168 N, each opened in bold · 31 `nb.recap` rows · 2 S in one `nb.pointer` row).
- **Derivations:** 29 of 29 `### Dnn ·` blocks in Part F and 29 `nb.derivation` rows in Part A, each in the CORE block the
  curation names; step counts equal the curation's (6, 9, 6, 7, 7, 14, 8, 9, 7, 9, 13, 10, 7, 8, 8, 11, 12, 7, 11, 7, 7, 10,
  9, 8, 13, 12, 6, 5, 7 = 253 steps); every step has did / tex / why / plain with a *why* of 6–35 words; every block has
  Goal, Start, Plan, Tools, Assumptions, Result, Check, What it means, Fails when, Traps; no step uses "it can be shown",
  "similarly", "clearly" or "as in the book".
- **★★★ checks:** the nine blocks D06, D11, D12, D16, D17, D19, D22, D25, D26 carry a `check_src` listing; all nine were
  extracted from this file and executed (sympy, each under 3 s): every assertion passes.
- **Explainers:** E1–E10 and B1, each with title, summary, CORE ids, reference explainer, meta, physics (JS ↔ Python), views,
  controls with the meaning of their extreme values, presets, status, readouts, depth features (≥ 2 beyond Explain,
  Code and Derivation), Explain sections ending in *Reading the current setting* with one text per regime, Derivation
  tab, Code, walkthrough (6, 6, 7, 6, 6, 7, 6, 7, 7, 7 steps; 5 for B1; none over 45 words), equations, ≥ 3 check
  questions, parity rows (≥ 4 each for E1–E10, 3 for B1; at least one southern-hemisphere row in E1, E2, E3, E4, E7, E8), and a fit plan for 360×640 and for
  844×345. Every ★★★ derivation of a block that has an explainer is in it (D06 in E3, D11 in E5, D12 in E6, D16 in E8,
  D22 in E9, D25 and D26 in E10); D17 (C13) and D19 (C14) are notebook-only, as the curation rules (B1 would carry D19).
- **Primers:** 29 new (P307–P335), each with an `nb.primer` row in Part A and a "primer (in Cxx)" row in Part E whose
  Concept is the primer term verbatim.
- **Ledger:** 231 rows, none without "Explained by"; every C / R id in that column exists in the curation; no primer is
  placed after its first use.
- **Part C:** 158 function rows (GFD 99 plus one row of constants · VM 11 · SW 20 including the `ShallowWater` class ·
  ch13 28), plus one row naming the 12 sympy engines and one naming the 8 `fig_*` helpers — 178 callables in all; every `ch13.` / `GFD.` / `VM.` / `SW.` name used in Parts A, B, D, E, F is
  in Part C; every Part C row is called by a storyboard row, mirrored by an explainer, or listed in C.7 as test-only. No
  name or signature was changed after the first save (19:39 EDT); later edits to Part C touched wording only (the
  word "value" in front of Returns cells that began with a symbol; `SW.relative_vorticity` added to the C.7 list).
- **Equation mentions:** 0 number-only mentions in the header, Part A, Part B and Part F (scratch scanner: every sentence,
  table cell or step field that carries a book number also carries the equation in TeX). Not counted, by design:
  heading lines (`#### Cnn —`, `### Dnn ·`), `viz:equations` meta lists, the `ref '…'` label of an equation card that shows
  its TeX, code listings, and the "Book eq." column of Part C. Numbers stand only beside the equation printed under
  them; the form $\omega^2=f^2\sin^2\theta+N^2\cos^2\theta$ is labelled "unnumbered, p. 669" everywhere.
- **Sentences and cells starting with a symbol:** 0 found by the scratch scanner (Greek letters and single-letter
  variables at the start of a sentence, a table cell outside a symbol column, or a why / plain field).
- **Control characters:** 0. **Public check:** `tools/check_public.py analysis/ch13_design.md` → OK.
- **Slips and traps:** 13 slip rows and 17 trap rows in the header; twelve slips have a callout in the house form where
  they are taught (#13 is deliberately not asserted); each trap has a callout row.
- **Cross-check against the implementer's code as it stood at 20:45 EDT (read-only, for information):** all Part C names
  exist; 21 functions carry extra trailing optional arguments (compatible); of 79 spot evaluations of this design's
  *expect* numbers, 77 agree; one differs — `shallow_water_omega` returns a slow root 0.18 % away from `numpy.roots` and
  from the trigonometric form (external mode, 3100 km, 35° N: −8.904 × 10⁻⁶ against −8.888 × 10⁻⁶ s⁻¹; its three roots do
  not sum to zero) — reported to the orchestrator, not edited; and one check of mine was too strict (the density flux of
  the Eady mode is zero, not negative, exactly on the lids).
- **Correction pass of 2026-10-07 (after the implementer's measurements; each claim re-run by the designer before
  editing):** weight-1 orthogonality of the vertical modes holds for both lids — the surface term belongs to the energy
  relation (Part C.2, C08 rows 7 and 20, D11 step 13, Check and Traps, E5); `orthogonality_matrix(modes, kind="psi" or
  "energy", …)`, `modal_amplitudes(c_n, p_n_t, p_n=None, rho0=1.0, g=G0)`, `book_slips()` keys, `traps()`,
  `rayleigh_kuo_eigs(..., **kw)` with β as its own keyword, and the NaN return of `shallow_water_omega` where the
  discriminant is negative are now in the contract; the WKB error falls as 1/(H_N m), not as its square; the Rayleigh–Kuo
  growth rates 0.1191, 0.0744, 0.0151 (β = 0, 0.3, 0.6; none at 0.7) replace the hypothesis row; the cached turbulence
  runs are 64² and every statement about them is qualitative (no slope is quoted); rigid-lid `Modes` start at the first
  baroclinic mode and `project` drops the depth mean; the westward step flow is a labelled remark, not a slip; E8's bump
  presets stay a JS convolution labelled "computed in the explainer; the step case is checked against Python". The
  earlier remark above about a slow-root discrepancy in `shallow_water_omega` is obsolete: re-run, it returns
  −8.88828 × 10⁻⁶ s⁻¹. Item 5 of the orchestrator's list (energy drift of the nonlinear model, Poincaré errors of `run`)
  needed no edit: this design states no value for either.
- **Still a hypothesis in its row:** the three `ekman_solve` transports with $K_v(z)$ (C05 row 9) and the tilt of the
  fastest Eady mode (C16 row 16) — not computed by the designer.
