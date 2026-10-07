# Chapter 12 verification — Turbulence                       2026-10-07 · code at `c21490a` (HEAD `5e1db2f` differs only in `analysis/ch12_design.md`)

## Verdict: PASS — 117 tests in `tests/test_ch12.py`, 0 failing (122 s); full suite 1737 passed (42 min 44 s); 42 of 42 planted mutants killed

Every script ran, every test passes, every computable CORE row has ≥ 2 independent evidence levels (one of V1/V2/V3/V5 in
each), every one of the 215 names of design Part C is called by at least one test (a self-audit test enforces it), the seven
planted wrong variants of the book-slip table all fail their checks, and nothing is `unverified` without an Open item. No
physics defect was found in `fluidpy/`; the ten flagged items are settled below (one design expectation is wrong, six
docstring statements need correcting — listed, none of them changes a number a builder uses).

## Environment
python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 · Windows 11.
`tests/test_ch12.py`: **117 tests, 122 s** (114 in 50 s with `-m "not slow"`). After the full-suite run four illustrative inputs of the test file (a scalar half-width, a fitted spreading rate, one B of the (12.92) loop) were changed so that no number coincides with an entry of the book's table; the file was re-run alone: 117 passed in 103 s. Full suite: **1737 passed in 2564.16s (0:42:44)**, 0 failed (`.venv/Scripts/python.exe -m pytest -q -p no:cacheprovider`; 1620 before this phase + 117).

## The ten flagged items — rulings

| # | Item | Ruling | Evidence |
|---|---|---|---|
| 1 | `WT.fit_log_law` returns κ = 0.383, B = 3.99 on a Spalding profile at δ⁺ = 5000 over 30 < y⁺ < 0.15 δ⁺ | **Function right, design expectation wrong.** The fit recovers (κ, B) of an exact log profile to 1e-10. Spalding's curve reaches the logarithm from below (U⁺ − log law = −0.66 at y⁺ = 30, −0.16 at 100, −0.03 at 300; exactly ln(1 + c/y⁺)/κ with c = O(U⁺³)), so a window that starts at y⁺ = 30 is biased low — and the same happens with real data: on the Lee & Moser DNS at Re_τ = 5186 the window from 30 gives κ = 0.400, B = 4.92, the window from 3√Re_τ = 216 gives 0.3845, 4.31 and from 350 gives 0.3837, 4.28, inside the published κ = 0.384 ± 0.004. **What the notebook should use:** (a) the DNS subset in `reference/ch12` with `window=(350.0, 0.15)` → *expect* κ = 0.384 ± 0.004 (B ≈ 4.28); or (b) the composite profile at δ⁺ = 10⁵ with `window=(1000.0, 0.15)` → 0.4099, 5.001; at δ⁺ = 5000 the best the composite can do is `window=(300.0, 0.15)` → 0.405, 4.80 (state it as "0.405 ± 0.005, 4.8 ± 0.1", not 0.41 ± 0.01, 5.0 ± 0.2). | `test_fit_log_law_V1_exact_recovery_and_the_bias_of_a_composite_profile`, `test_log_law_V5_lee_moser_von_karman_constant_from_the_dns_profile` |
| 2 | `rans_eddy_viscosity_residual(printed=True)` reads the printed $\partial P/\partial x_j$ of (12.97) as "j summed" | **Accept, with a wording correction.** Page re-read: (12.97) prints $\frac{\partial U_i}{\partial t}+U_j\frac{\partial U_i}{\partial x_j}=-\frac1\rho\frac{\partial P}{\partial x_j}+\frac{\partial}{\partial x_j}\big([\nu+\nu_T](\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i})-\frac23\bar e\,\delta_{ij}\big)$ — the pressure term carries a lone j, so the equation is ill-formed (free-index mismatch); the corrected $\partial P/\partial x_i$ is what substituting (12.94) into (12.30) gives, and the code's corrected residual equals the (12.30) residual with the stress of (12.94) to 2e-5 on a manufactured field. The planted variant realises the ill-formed term as Σ_j ∂P/∂x_j; it is off by exactly (Σ_j ∂P/∂x_j − ∂P/∂x_i)/ρ and fails. Docstring should say "ill-formed; realised here as a sum" — the summation convention itself does not sum an index that occurs once. | `test_slip06_rans_eddy_viscosity_V1_corrected_matches_12_30_and_printed_differs` |
| 3 | `dispersion_regime` inclusive boundaries; an E10 parity row sits on t = 0.3 Λ_t | **Accept; pinned.** t = 0.3 Λ_t → "ballistic", t = 3 Λ_t → "diffusive", for Λ_t = 0.1, 1, 7, 10, 123.4, with both neighbours at ±1e-9; `dispersion_regime(3.0, 10.0) == "ballistic"`. The JS mirror must compare `t/Λ <= 0.3` on the **ratio** (as the Python does), not `t <= 0.3*Λ`. Mutant M25 (exclusive boundaries) is killed. | `test_dispersion_regime_V7_inclusive_boundaries_and_both_sides` |
| 4 | `surface_layer_wind(unstable="businger_dyer")` without ψ_m(z₀/L) by default | **Accept both.** With `psi_at_z0=True` the function equals $\frac{u_*}\kappa\int_{z_0}^{z}\phi_m\,dz/z$, φ_m = (1 − 16ζ)^{−1/4}, by quadrature to 1e-10 and U(z₀) = 0 exactly. The default equals the common form $\frac{u_*}\kappa[\ln\frac z{z_0}-\psi_m(z/L_M)]$ with ψ_m obtained independently as ∫₀^ζ(1 − φ_m)/ζ′dζ′ to 1e-10; the two differ by exactly (u_*/κ)ψ_m(z₀/L_M) > 0, smaller than 4.2 (u_*/κ) z₀/∣L_M∣ (0.004 m/s for L_M = −20 m, z₀ = 0.03 m), and U(z₀) = −(u_*/κ)ψ_m(z₀/L_M) ≠ 0. | `test_surface_layer_wind_V5_businger_dyer_integral_with_and_without_the_z0_term` |
| 5 | `surface_layer_regime` nominal crossover ∣z/L_M∣ = 1; "forced convection" for either sign | **Accept as documented** (the book's meaning: z ≪ ∣L_M∣). Both sides of the boundary, both signs, ±∞ and the `crossover` argument are pinned; a test shows the same string for a stable and an unstable layer while `surface_layer_state(...)["verdict"]` says "stable"/"unstable" — E9 must take stability from `verdict`, never from `layer`. | `test_surface_layer_regime_V7_boundaries_both_signs_and_neutral` |
| 6 | `channel_energy_budget` / `_at`: magnitudes flat, signs nested; closure, identity, table | **Accept.** Flat: `pressure_work + transport − mean_dissipation − production` = 0 to 1e-12 at every height; nested `mean` sums to 0 and equals minus the flat magnitudes; `turbulence` closes by residual (labelled model). Integral identity work = dissipation + production: residual 6.8e-5 → 5.2e-7 for n = 200 → 1600 (order 2), and exactly (1e-10) by an independent quadrature of my own closed forms; ∫U⁺/Re_τ dy⁺ = U_bulk⁺. `channel_energy_budget_at` equals my closed forms to 1e-11 and U⁺ by quadrature to 1e-9. The JSON equals a fresh `explainer_tables()` and `write_reference_tables` reproduces the file byte for byte. **One note for the explainer builders:** the stored y⁺ is itself rounded to 6 digits, so evaluating the function at the *stored* height reproduces steep or cancelling columns (`transport`, `uv_plus` near the centre, `mean_dissipation` far out) only to ≈ 3e-5 relative — parity rows that read a table value need rtol 1e-4, or should call the function at a round y⁺. | `test_channel_energy_budget_V4_*`, `_V3_*`, `test_channel_energy_budget_at_V1_*`, `test_explainer_tables_V1_*` |
| 7 | No public DNS in `reference/ch12/`, no `SOURCES.md` | **Done.** `reference/ch12/make_refs.py` downloads the five Lee & Moser mean-profile files (HTTP 200, 2026-10-07) and writes a thinned, un-interpolated subset (281 rows) + `benchmarks.json` + `SOURCES.md` (source, URL, terms, date). Used as V5 for `viscous_sublayer`, `log_law_indicator`, `fit_log_law`, `log_law`, `spalding_uplus`, `stress_partition`, `channel_mixing_length`. | five `*_V5_*lee_moser*` tests |
| 8a | `model_spectrum` c_L = 6.78 "unverified" | **Verified two ways.** (i) Numerically: c_L is the value for which ∫E dK = (ε̄L)^{2/3}, i.e. L = ē^{3/2}/ε̄; the root is 6.415, 6.745, 6.779 at L/η = 5.6e3, 1.8e5, 5.6e6 → 6.78 in the high-Reynolds limit (C = 1.5, p₀ = 2). (ii) Open sources quoting Pope's p₀ = 2, c_L = 6.78 (arXiv:1705.04917; ATOMIX wiki), read 2026-10-07. Docstring should drop "unverified". | `test_model_spectrum_V4_dissipation_and_energy_integrals` |
| 8b | `plane_jet_stress_profile`: C₃G = −½F∫F, opposite to analysis row 121 | **Code and design are right; analysis row 121 has the wrong sign.** Page re-read: (12.63) prints $\{\delta U'_{CL}/U_{CL}\}F^2-\{\delta U'_{CL}/U_{CL}+\delta'\}F'\int_0^\xi F\,d\xi=\{\Psi/U_{CL}^2\}G'$; with −½ and +½ the left side is d(−½F∫F)/dξ (sympy, any F). Independently of every similarity formula, integrating $U\,\partial U/\partial x+V\,\partial U/\partial y=-\partial\overline{uv}/\partial y$ (12.61) outward from the axis with V from continuity reproduces the coded −mean(uv) to 2e-5 of its maximum: negative for y > 0, where ∂U/∂y < 0, so ν_T > 0. Mutant M15 (the row-121 sign) is killed. | `test_plane_jet_V2_derivation_D15_*`, `test_plane_jet_stress_V1_sign_from_the_momentum_equation_by_quadrature` |
| 9 | Lapse-rate convention: Kundu's Γ = dT/dz computed, meteorological Γ = −dT/dz alongside | **Both verdicts agree** for dT/dz = −6.5, 0, +10, −12 K/km and Γ_a itself: the Kundu text carries `>`/`<`/`=` and the meteorological text the mirrored comparison with all numbers negated; both equal `core.stratification.lapse_rate_stability` verbatim; Ri is identical whichever convention is asked for; Ri from (in-situ, Γ_a) equals Ri from dθ/dz to 1e-12 and N² equals ch01's `brunt_vaisala_sq_from_lapse`; the no-Γ_a mutant calls the standard atmosphere unstable (and M21 is killed). **Note for E9:** the strings do not begin with the verdict — an unstable layer reads "stable ⇔ dT/dz > Γa: −12.0 < −9.8 K/km" (criterion, then the comparison that actually holds). Show the bare word from `["verdict"]` beside it. | `test_gradient_richardson_V1_in_situ_and_potential_routes_and_both_conventions`, `test_surface_layer_state_V1_*` |
| 10 | `stratified_tke_budget(dUdz=…)` added after the implementer's full-suite run | **Covered and the full suite re-run.** One height with `dUdz`, a profile with and without it (they agree to the (Δz/z)²/3 differencing error of a log profile), Rf = `flux_richardson` = z/L_M, and the one-height call without `dUdz` raises. Full suite: **1737 passed in 2564.16s (0:42:44)**, 0 failed (`.venv/Scripts/python.exe -m pytest -q -p no:cacheprovider`; 1620 before this phase + 117). | `test_stratified_tke_budget_V1_terms_and_optional_shear_argument` |

## Validation table (CORE rows; every computable row has ≥ 2 independent levels incl. one of V1/V2/V3/V5)

| ID · concept (equation) | fluidpy target | Evidence | Numbers | Label |
|---|---|---|---|---|
| C01 · ensemble average and its rules, $\overline{u^m}=\frac1N\sum_nu^m$ (12.1), (12.2), (12.4)–(12.9) | `TS.ensemble_average`, `moment`, `central_moment`, `time_average`, `check_averaging_rules`, `product_average_split`, `ch12.time_average_exp_cos` | V1 · V2 · V3 · stat | hand loops 1e-13; six linear rules < 1e-11, product rule = covariance; Example 12.1 window integral by sympy; sliding average order 2.00; Gaussian moments within 5 s.e. (N = 2e5); scatter ∝ N^−0.53 | analytic |
| C02 · $R_{11}(\tau)=\overline{u_1(t)u_1(t+\tau)}$ (12.17), $\Lambda_t=\int_0^\infty r_{11}d\tau$ (12.18), (12.19) | `TS.autocorrelation`, `cross_correlation`, `integral_scale`, `taylor_microscale`, `correlation_spectrum_pair` | V1 · V2 · stat | FFT = direct 1e-10; cosine within its geometric-sum bound; Λ_t, λ_t of three closed pairs ≤ 1e-5; D02, D04 by sympy; OU record error ∝ T^−0.55 (16 seeds) | analytic |
| C03 · $S_e(\omega)=\frac1{2\pi}\int R_{11}e^{-i\omega\tau}d\tau$ (12.20)–(12.22) | `TS.spectrum_from_correlation`, `correlation_from_spectrum`, `spectrum_variance`, `periodogram` | V2 · V1 · V3 · V4 · stat | transforms of e^{−τ/τc}, e^{−τ²/τc²} by sympy, ∫S = σ², S(0) = σ²Λ/π; Filon rule within its interpolation bound, order 2.00; Parseval 1e-12; = scipy `welch`/4π within 2 %; Lorentzian within 5 s.e. per bin | symbolic |
| C04 · Reynolds-averaged momentum (12.30), $\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}$ | `ch12.rans_sympy`, `rans_momentum_residual`, `mean_stress_tensor`, `TS.reynolds_stress`, `displaced_parcel_uv` | V2 · V3 · V1 · stat | finite four-member ensemble: mean(member equation) − (12.30) ≡ 0, dropping the stress ≠ 0; engine checks 5/5; manufactured 3-D field order 2.00; parcels within 5 s.e. | symbolic |
| C05 · $\bar\varepsilon=30\nu\overline{u^2}/\lambda_f^2=15\nu\overline{u^2}/\lambda_g^2$ (12.43), (12.40)–(12.41) | `gradient_moments_isotropic_sympy`, `isotropic_correlation_tensor`, `transverse_from_longitudinal`, `isotropic_scales`, `dissipation_isotropic`, `dissipation_rate` | V2 · V1 · V4 · stat | D07: g = f + (r/2)f′ solved from ∂R_ij/∂r_j = 0, Λ_g = Λ_f/2, λ_g² = λ_f²/2; D08: moments (2, 4, −1), full double sum = 30; three routes agree 2e-5; ε̄ = νΣK²∣û∣² 1e-10 on a 32³ field; ratio 2 (3-D) and 3 (2-D) within 5 s.e. | symbolic |
| C06 · energy budgets (12.46), (12.47), production $-\overline{u_iu_j}\,\partial U_i/\partial x_j$ | `mean_energy_budget_sympy`, `tke_budget_sympy`, `mean_energy_budget`, `tke_budget`, `channel_energy_budget(_at)` | V2 · V4 · V3 · V1 | D09, D10 by the finite ensemble (≡ 0; wrong sign / factor / triple-correlation ½ ≠ 0); ½ trace of (12.35) = (12.47); Poiseuille work = dissipation 1e-5, residual = its exact truncation term; channel identity 1e-10 by quadrature; order 1.99 | symbolic |
| C07 · $\eta=(\nu^3/\bar\varepsilon)^{1/4}$, $u_K=(\nu\bar\varepsilon)^{1/4}$ (12.50), (12.48)–(12.52) | `kolmogorov_scales`, `kolmogorov_exponents`, `scale_separation`, `dissipation_outer_scaling`, `scale_table`, `cascade_tiers` | V2 · V1 | exponents (¾, −¼), (¼, ¼), (½, −½) by a rational solve and by ch01's `solve_exponents`; pint units; ηu_K/ν = 1 to 1e-14; slopes −¾, −½, +½ to 1e-12 | symbolic |
| C08 · $S_{11}=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ (12.54, corrected), (12.53), (12.55) | `inertial_spectrum_1d/3d`, `kolmogorov_constants`, `one_dimensional_from_3d`, `fit_inertial_range`, `model_spectrum`, `kolmogorov_normalize_spectrum` | V2 · V1 · V4 · V5 | D12: (⅔, −5/3) and (5/4, ¼); 18/55 by sympy and by quadrature 1e-8; printed +5/3 fails units and slope; 2ν∫K²E dK = ε̄ 1e-8; collapse 1e-12; C₁ = 0.491 vs 0.53 ± 0.055 (Sreenivasan 1995) | symbolic |
| C09 · plane jet $U=C_5(J_s/\rho)^{1/2}x^{-1/2}F(y/x)$ (12.66), (12.62), (12.63), (12.68) | `plane_jet_*`, `plane_jet_similarity_sympy`, `general_similarity_check`, `free_shear_exponents`, `free_shear_flow` | V2 · V4 · V1 | D13–D16 by sympy (two explicit profile pairs, generic δ, U_CL, Ψ); J_s constant 1e-10 at four stations, scalar flux 1e-9; stress from the momentum equation by quadrature 2e-5; seven exponent systems exact | symbolic |
| C10 · $U^+=f(y^+)$ (12.80), $u_*^2=\tau_0/\rho$ (12.81), $U^+=y^+$ (12.82), $\bar\tau=\tau_0(1-2y/h)$ | `WT.wall_units`, `friction_velocity`, `viscous_sublayer`, `channel_total_stress`, `channel_pressure_gradient`, `law_of_the_wall_groups` | V2 · V1 · V5 | Π groups exact; D17 by sympy; round trips 1e-14; DNS: U⁺ = y⁺ − y⁺²/(2Re_τ) within 5e-4 for y⁺ < 1 at five Re_τ, a few per cent by y⁺ = 5 | benchmark |
| C11 · $U^+=\frac1\kappa\ln y^++B$ (12.88), (12.84)–(12.89), (12.93) | `WT.log_law`, `log_law_defect`, `fit_log_law`, `log_law_indicator`, `friction_law_from_overlap`, `rough_wall_log_law`, `overlap_matching_sympy` | V2 · V1 · V5 | D19 by `dsolve`; exact recovery 1e-10; DNS Re_τ = 5186: κ = 0.3837 (window 350–0.15δ⁺), 0.3845 (from 3√Re_τ), indicator 0.3834 — all inside 0.384 ± 0.004 | benchmark |
| C12 · $\overline{u_iu_j}=\frac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$ (12.94), (12.98)–(12.101) | `eddy_viscosity_stress`, `mixing_length_wall_profile`, `mixing_length_intercept`, `shear_flow_eddy_viscosity_solve`, `channel_mixing_length` | V2 · V1 · V3 · V5 (approximate) | D21: root, limits, closed form, B = (ln 4κ − 1)/κ by sympy; van Driest profile by independent quadrature 1e-8; trapezoid order 2.00; BVP solver order 2.00 and = closed-form channel to 8e-5; vs DNS ∣ΔU⁺∣ ≤ 0.68, C_f within 1.1 % (Re_τ ≥ 550) | converged |
| C13 · k–ε: (12.103), $\nu_T=C_\mu\bar e^2/\bar\varepsilon$ (12.104), (12.105) | `k_epsilon_decay`, `k_epsilon_rhs`, `k_epsilon_eddy_viscosity`, `k_epsilon_loglayer_kappa`, `K_EPSILON_CONSTANTS` | V2 · V1 · V5 | D23: n = 1/(C_ε2 − 1) and κ² = √C_μ(C_ε2 − C_ε1)σ_ε solved by sympy; closed form = `solve_ivp` = RK4 to 1e-7; constants = published set | symbolic |
| C14 · $\mathrm{Rf}=\frac{-g\alpha\overline{wT'}}{-\overline{uw}\,dU/dz}$ (12.107), (12.108), (12.109) | `flux_richardson`, `gradient_richardson_thermal`, `turbulence_regime`, `stratified_tke_budget`, `turbulent_prandtl` | V1 · V2 · V7 | hand value 1e-14, signs, ±inf/NaN at zero shear; D24: Ri = (ν_T/κ_T)Rf by sympy and numerically 1e-13; two temperature routes 1e-12; regime boundaries | analytic |
| C15 · $L_M=-u_*^3/(\kappa\alpha g\overline{wT'})$ (12.110), $\mathrm{Rf}=z/L_M$ (12.111), log-linear profile | `monin_obukhov_length`, `monin_obukhov_from_fluxes`, `flux_richardson_surface_layer`, `WT.surface_layer_wind`, `dimensionless_shear`, `surface_layer_state` | V1 · V2 · V5 · V7 | Rf from (12.107) with the log-law gradient = z/L_M 1e-13; profile by sympy integration; φ_m = (κz/u_*)dU/dz 1e-6; Businger–Dyer integral by quadrature 1e-10 | analytic |
| C16 · $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\,t\int_0^t(1-\tau/t)r_\alpha d\tau$ (12.119), (12.120)–(12.129) | `taylor_dispersion*`, `taylor_dispersion_rate`, `eddy_diffusivity_*`, `langevin_particles`, `random_walk`, `dispersion_regime`, `smoke_plume_width` | V2 · V1 · V4 · stat · V7 | D26–D28 by sympy (generic r; exponential; Gaussian); (12.118) = (12.119) = closed forms ≤ 1e-8; 2e4 Langevin particles within 5 s.e. at 40 times, mean(Xu) = D_T; mean(R_n²) = nL² within 5 s.e. in 1, 2, 3 dimensions | symbolic |

**Counts per validation label (CORE rows):** symbolic 9 · analytic 4 · benchmark 2 · converged 1 · conserved 0 as a primary label
(V4 is a second level in C03, C05, C06, C08, C09, C16) · qualitative 0 · unverified 0.
**Tests by primary level:** V1 45 · V2 30 (23 named `_V2_derivation`) · V3 6 · V4 9 · V5 8 · V6 2 (private) · V7 8 · stat 9.
Coded NOTE items and helpers: all 215 contract names are exercised (`test_contract_V1_every_part_c_name_is_called_by_a_test`);
labelled **qualitative**: `jet_tke_budget`, `k_epsilon_channel`, the free-shear amplitude constants; labelled **approximate**
against DNS: `spalding_uplus`, `stress_partition`, `channel_mixing_length`, `mixing_length_intercept(A_plus=26)`.

## Derivations (curation §4b) — every ★★ / ★★★ row re-derived, independently of the chapter's sympy engines

| D | Result | Test | How |
|---|---|---|---|
| D01 ★ | rules (12.4)–(12.9), $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$ | `test_averaging_rules_V2_derivation_D01_*` | three-member symbolic ensemble |
| D02 ★★ | Schwartz inequality (12.16) | `test_schwartz_V2_derivation_D02_*` | minimum of the quadratic in λ; discriminant |
| D04 ★★ | $\lambda_t^2=-2/r''(0)=2\overline{u^2}/\overline{(du/dt)^2}$ (12.19) | `test_taylor_microscale_V2_derivation_D04_*` | random-phase sum; osculating parabola by series |
| D05 ★★ | pair (12.20)–(12.21), (12.22), $S_e(0)=\overline{u^2}\Lambda_t/\pi$ | `test_spectrum_pair_V2_derivation_D05_*` | sympy integrals, two correlations |
| D06 ★★ | (12.27), (12.28), (12.30) | `test_rans_V2_derivation_D06_*` | finite ensemble of polynomial solenoidal fields, exact averages |
| D07 ★★★ | (12.40) → (12.41), g = f + (r/2)f′, Λ_g = Λ_f/2, λ_g = λ_f/√2 | `test_isotropic_tensor_V2_derivation_D07_*` | divergence of the general form, solve for g, three f |
| D08 ★★★ | moments 2 : 4 : −1, factor 30 (12.42)–(12.43) | `test_isotropic_dissipation_V2_derivation_D08_*` | −∂²R_ij/∂r_k∂r_l at 0, full double sum |
| D09 ★★ | mean-flow energy budget (12.46) | `test_mean_energy_V2_derivation_D09_*` | U_i × mean(member equation) |
| D10 ★★★ | turbulent energy budget (12.47) | `test_tke_budget_V2_derivation_D10_*` | mean(u_i × member equation); intermediate line mean(u_i·mean equation) = 0 |
| D12 ★★ | (12.53), (12.54) with −5/3, (12.55), 18/55 | `test_inertial_spectrum_V2_derivation_D12_*` | exponent solves; isotropic integral |
| D13 ★★ | (12.61) → (12.62) | `test_plane_jet_V2_derivation_D13_*` | flux form identity |
| D14 ★★★ | similarity equation (12.63) | `test_plane_jet_V2_derivation_D14_*` | generic δ, U_CL, Ψ; steps for ∂U/∂x and V checked line by line |
| D15 ★★ | γ = −½, C₃G = −½F∫F, V̇ ∝ x^{1/2} (12.64)–(12.68) | `test_plane_jet_V2_derivation_D15_*` | solve; derivative identity for any F |
| D16 ★★ | exponents of the seven flows | `test_free_shear_exponents_V2_derivation_D16_*` | seven linear systems |
| D17 ★★ | $\bar\tau=\tau_0(1-2y/h)$, (12.90), (12.91) | `test_channel_stress_V2_derivation_D17_*` | integrate, symmetry, plug balance |
| D19 ★★ | (12.85)–(12.89), friction law | `test_log_law_V2_derivation_D19_*` | chain rule, `dsolve`, sum |
| D21 ★★ | mixing-length root, limits, closed form, intercept | `test_mixing_length_V2_derivation_D21_*` | solve, limits, derivative |
| D23 ★★ | n = 1/(C_ε2 − 1), κ² = √C_μ(C_ε2 − C_ε1)σ_ε | `test_k_epsilon_V2_derivation_D23_*` | substitution into the pair and into (12.105) |
| D26 ★★ (+ D27, D28 ★) | (12.117)–(12.119), limits, D_T | `test_taylor_dispersion_V2_derivation_D26_*` | Leibniz rule on a generic r; two closed cases |
| D18, D20, D24, D25 ★ | wall groups; rough wall; Ri = Pr_T Rf; Rf = z/L_M and the log-linear profile | `test_wall_units_V2_*`, `test_rough_wall_V1_derivation_D20_*`, `test_richardson_V2_derivation_D24_*`, `test_monin_obukhov_V1_*_D25` | — |
| (12.112), slip #15 | temperature-variance budget | `test_temperature_variance_V2_derivation_budget_and_printed_factor` | finite ensemble; the printed κ∂mean(T′²)/∂z does not close |

Intermediate lines of design Part F checked explicitly: D14 steps 4, 9 and 13 (∂U/∂x, V in similarity form, the divided
equation), D10's use of mean(u_i × mean equation) = 0, D07's divergence r_i(rF′ + 4F + G′/r) and D15 step 8 (the exact
derivative −½(FI)′) — all correct. The other derivations were re-derived from their starting equations to their results
without comparing every written line (Part F was being edited while this phase ran). D03, D11, D22 (★) are covered by
V1/V2 tests of their results (`R_uv(τ) = R_vu(−τ)`; exponent solves; units of C_μē²/ε̄).

## Book slips — planted wrong variants (all must fail; all do)

| Slip | Wrong variant | What the test shows |
|---|---|---|
| #1 (12.54) | `inertial_spectrum_1d(printed=True)` | slope +5/3, units not m³/s², 10⁵ × too large at k₁ = 100 |
| #3 (12.70) | the integrand "at y = 0" | ρ∫ȲU dy = Ṁ_s at every x (1e-9); the y = 0 value falls as 1/x |
| #5 (12.74) | exponential family in `general_similarity_check` | c₂ ≡ 0 (sympy and numerically), momentum-flux exponent −1.4 ≠ 0; power family passes |
| #6 (12.97) | `rans_eddy_viscosity_residual(printed=True)` | off by (Σ_j∂P/∂x_j − ∂P/∂x_i)/ρ, > 0.1 |
| #12 (12.129) | `eddy_diffusivity_asymptote(printed=True)`; `dispersion_regime` | 100 × the exact D_T at t = 0.01 Λ_t; regimes named ballistic/diffusive on the right sides |
| #13 Fig. 12.27 | `smoke_plume_width` slopes | d ln Z/d ln x = 1 near the source, ½ far away (1e-3) |
| #15 (12.112) | `temperature_variance_sympy()["printed_check"]` | False; and my own ensemble shows the missing ½ |

## Functions used by the notebook and explainers (design Part C)
All 215 names (C.1: 39 · C.2: 41 + `LOG_LAW_CONSTANTS` · C.3: the rest) exist on `ch12`, `ch12` re-exports every public
function of `TS` and `WT` by identity, each docstring carries Book/Assumptions/Validation, and each name occurs in a test
body (three `test_contract_*` tests). All 47 `py:` parity expressions of the storyboards (E1–E10, B1) evaluate, eleven of them
are pinned against closed forms written out in the test. Drawing helpers of `scripts/ch12_drawings.py`: one smoke test each.
`rans_momentum_residual` (contract-only): manufactured 3-D field, order 2.00.

## Convergence studies

| Scheme | Grids | Observed order | Design |
|---|---|---|---|
| `TS.time_average` (trapezoid window) | 201 … 1601 samples | 2.00 | 2 |
| `TS.spectrum_from_correlation` (Filon) | lag step h, 2h | 2.00 (exp.), 2.00 (Gauss.) | 2 |
| `rans_momentum_residual` (central differences) | h = 0.08 … 0.01 | 2.00 | 2 |
| `rans_2d_residual` | h = 0.08 … 0.01 | 2.00 | 2 |
| `channel_mixing_length` U_cl⁺ (cumulative trapezoid, stretched grid) | n = 200 … 1600 | 1.99 | 2 |
| `mixing_length_wall_profile(method="trapezoid")` | 201 … 1601 | 2.00 | 2 |
| `shear_flow_eddy_viscosity_solve` (conservative FD, smooth ν_T) | n = 21 … 161 | 2.00 | 2 |
| `WT.boundary_layer_stress_from_profile` (Blasius) | (21, 201) → (41, 401) | error ratio between 3 and 5 | 2 |
| `mean_energy_budget` transport integral (Poiseuille) | 401 → 1601 | 2.0 ± 0.25 | 2 |
| statistical: scatter of an N-member mean | N = 25 … 1600 | −0.53 | −½ |
| statistical: rms error of R(τ) from an OU record | T = 500 … 8000 s, 16 seeds | −0.55 | −½ |
| `frozen_field_probe` error vs u_rms/U₀ | 0.002 … 0.008 | 0.99 | 1 |

## Conservation / invariant residuals

| Invariant | Residual |
|---|---|
| Parseval, `TS.periodogram` (Σ w S = variance) | < 1e-12 |
| ∇·u of `synthetic_solenoidal_field` (2-D, 3-D) | < 1e-11 u_rms/Δx; commutation of mean and divergence < 1e-11 |
| ε̄ = ν Σ K²∣û∣² (32³ field) | 1e-10 relative |
| 2ν∫K²E dK = ε̄ (`model_spectrum`, no L) | < 1e-8 |
| plane-jet momentum flux ρ∫U²dy = J_s (x = 0.5, 1, 2, 8) | < 1e-10 |
| slot-fluid flux ρ∫ȲU dy = Ṁ_s | < 1e-9 |
| free-shear invariants of six flows (x = 1, 2, 4) | < 1e-9 |
| channel: work = dissipation + production (independent quadrature) | 1e-10; function's `identity_residual` 5.2e-7 at n = 1600 |
| plume: ∫c dz = Q/U | < 1e-9 |
| d mean(X²)/dt = 2 mean(Xu) (Langevin particles) | within sampling noise (< 0.5 % mean difference) |

## Benchmarks used (`reference/ch12/SOURCES.md`, all read 2026-10-07)

| Value | Ours | Source |
|---|---|---|
| κ = 0.384 ± 0.004 (channel DNS, Re_τ = 5186) | 0.3837 (fit 350 < y⁺ < 0.15 δ⁺), 0.3845 (from 3√Re_τ), 0.3834 (indicator) | Lee & Moser, J. Fluid Mech. 774, 395 (2015), doi:10.1017/jfm.2015.268, abstract (arXiv:1410.7809) |
| DNS mean profiles U⁺(y⁺), dU⁺/dy⁺ at Re_τ = 182, 543, 1001, 1995, 5186 | sublayer within 5e-4 (y⁺ < 1); Spalding(0.41, 5.0) within 0.77–0.85 wall units for y/δ < 0.15 (**approximate**); mixing-length channel within 0.54–0.68 wall units and C_f within +0.5 … +1.1 % for Re_τ ≥ 550 (**approximate**; 1.16 and +5.5 % at Re_τ = 182) | same; files at turbulence.oden.utexas.edu/channel2015/data |
| one-dimensional Kolmogorov constant 0.53, s.d. 0.055 | 18·1.5/55 = 0.491 (inside one s.d.; 7 % below the mean) | Sreenivasan, Phys. Fluids 7, 2778 (1995) via the ATOMIX wiki |
| C_μ, C_ε1, C_ε2, σ_k, σ_ε = 0.09, 1.44, 1.92, 1.0, 1.3 | identical | Launder & Sharma (1974) (OpenFOAM, SimScale documentation) |
| Prandtl's law 1/√f = 2.0 log₁₀(Re√f) − 0.8 | residual 1e-10; log-law-derived constants 1.986 and −1.020, f higher by 5–9 % for Re_d = 10⁴–10⁷ | Prandtl (1935) via McKeon et al., J. Fluid Mech. 538, 429 (2005) |
| φ_m = (1 − 16ζ)^{−1/4}, ζ < 0 | identical; ψ_m proved to be its integral | AMS Glossary, "Businger–Dyer relationship" |
| Pope's c_L = 6.78, p₀ = 2 | root of ∫E dK = (ε̄L)^{2/3}: 6.779 at L/η = 5.6e6 | Pope (2000) via arXiv:1705.04917 and the ATOMIX wiki |

## Numbers from the text (book vs ours) — private values redacted to relative differences
Two `_V6_` tests read `tests/book_values_ch12.json` (git-ignored; they skip when it is absent). Reproduced: the stirred-vessel
Kolmogorov length (< 0.5 %); the two-sided spectral constant from the three-dimensional one (< 2 %, the book rounds); all
velocity and scalar exponents of the free-shear table (exact); the fuel-jet example — stoichiometric fractions (< 0.5 %),
distance (< 0.5 %), centreline speed (< 2 %, two-digit book value); the flat-plate example — Re_x (< 1 %), θ, δ*, δ₉₉ (< 0.5 %);
(12.92) for the book's three (κ, B) pairs (< 2 %); the pipe-law constants of the exercise (1e-5); the convection example
(< 6 %, one-digit book values); the k–ε constants (exact); Rf_cr and the log-linear coefficient (exact defaults).

## Figures reproduced with our code (`outputs/ch12/verify/`, local; `tests/ch12_verify_figures.py`)

| Figure | PNG | Visual verdict |
|---|---|---|
| 12.18 analogue | `fig_12_18_law_of_the_wall.png` | DNS points at three Re_τ fall on U⁺ = y⁺ below y⁺ ≈ 5, bend through the buffer layer and run parallel to the log line; the indicator peaks near 5.9 at y⁺ ≈ 10 and sits on a plateau at 1/0.384 for 300 < y⁺ < 800 only at Re_τ = 5186; Spalding lies slightly below the data in the buffer layer. |
| 12.12 analogue | `fig_12_12_spectrum_collapse.png` | Three model spectra collapse on one roll-off for k₁η ≳ 0.03, follow the −5/3 line over a range that lengthens with L/η, and flatten at their own outer scale. |
| 12.13 analogue | `fig_12_13_plane_jet_similarity.png` | Profiles widen and slow downstream, collapse exactly in ξ = y/x; the stress is odd, zero on the axis, negative for ξ > 0 (peak ∣·∣ ≈ 0.022 U_CL² near ξ ≈ 0.07). |
| 12.10/12.15 analogue | `fig_12_10_channel_energy_budgets.png` | Near the wall transport (+1) balances direct dissipation (−1); production rises to 0.245 at y⁺ ≈ 10 (below the bound ¼) and is mirrored by the residual sink; everything decays outward. |
| 12.5 analogue | `fig_12_5_correlation_and_spectrum.png` | Measured r(τ) of an OU record lies on e^{−τ/τc} (Λ_t = 0.498 vs 0.5); the periodogram follows the Lorentzian and lifts above it only near the Nyquist frequency (aliasing of the sampled record). |
| 12.21 analogue | `fig_12_21_surface_layer_wind.png` | Stable profile bends to more wind aloft, unstable (Businger–Dyer) to less, all three meet at z₀; the book's log-linear form in strongly unstable air turns back to U < 0 above z ≈ 4 m — outside its range, as documented. |
| 12.25–12.27 analogue | `fig_12_25_taylor_dispersion.png` | Particle mean(X²) lies on Taylor's curve, slope 2 then 1 with the crossover near t ≈ 2Λ_t; the plume envelope is a wedge near the source and a parabola far away. |

## Deviations & justifications (every `# DEVIATION` in the code)
1. (12.97): ∂P/∂x_i coded, the printed ∂P/∂x_j kept as `printed=True` — slip #6, item 2 above.
2. (12.54): exponent −5/3 coded, +5/3 kept as `printed=True` — slip #1.
3. `plane_jet_stress_profile`: C₃G = −½F∫F against analysis row 121 — item 8b above (row 121 is wrong).
4. `mixing_length_stress`: l_T²∣dU/dy∣dU/dy instead of the printed (dU/dy)² — keeps the sign of the stress; odd in the shear (tested).
5. `temperature_variance_budget`: κ∂(½mean(T′²))/∂z — slip #15, proved by the finite ensemble.
6. `eddy_diffusivity_asymptote`: condition t ≫ Λ_t for (12.129) — slip #12.
7. `TS.periodogram`: segment averages and tapers are rescaled so that Σ w S equals the variance of the samples used — a
   normalisation, stated in the docstring; Parseval is evidence only for one boxcar segment (that is the case the V4 test uses).

## Open items (verbatim for the orchestrator)
1. **Design expectation to change (lesson-designer / notebook-builder):** the `fit_log_law` *expect* "κ = 0.41 ± 0.01,
   B = 5.0 ± 0.2 from a composite profile at δ⁺ = 5000 over 30 < y⁺ < 0.15 δ⁺" cannot be met by a correct fit (it returns
   0.383, 3.99). Use the Lee & Moser subset with `window=(350.0, 0.15)` (expect 0.384 ± 0.004), or the composite profile at
   δ⁺ = 10⁵ with `window=(1000.0, 0.15)` (0.4099, 5.001), or δ⁺ = 5000 with `window=(300.0, 0.15)` (0.405, 4.80).
2. **Analysis row 121** writes C₃G = +½F∫F; (12.63) gives −½F∫F. Correct the analysis note (code and design are right).
3. **Docstring corrections for the implementer (no number changes):** (a) `model_spectrum`: drop "unverified" for c_L = 6.78
   (verified, see item 8a); (b) `skin_friction_zpg`: the three laws agree within 10 % up to Re_x = 10⁸ and within 12 % at 10⁹,
   not "within 10 % to 10⁹"; (c) `stress_partition`: equal viscous and Reynolds shares at y⁺ ≈ 9.9 (Spalding, κ = 0.41,
   B = 5.0), not ≈ 11; (d) `mixing_length_intercept`: the damped value 5.277 is 5.5 % above B = 5.0 and 1.5 % above 5.2 — say
   which intercept "within 5 %" refers to, and label it approximate (with κ = 0.384 it gives 5.10 against the DNS fit 4.28);
   (e) `explainer_tables`: rows reproduce the function "to the 6 digits stored" only at the un-rounded heights — at the
   stored y⁺ to ≈ 3e-5; (f) `rans_eddy_viscosity_residual`: "ill-formed, realised as a sum" (item 2); (g)
   `WT.dimensionless_shear`: the docstring attributes the stable coefficient 4.7 to the AMS Glossary; today's search summary
   of that entry gives 1 + 5 z/L (4.7 is Businger et al. 1971) — the page itself returned HTTP 403, so this needs a human look.
4. **Explainer builders:** E10 — compare `t/Λ_t <= 0.3` on the ratio (item 3); E9 — take the stability word from
   `["verdict"]`, the strings `verdict_kundu`/`verdict_met` start with the criterion "stable ⇔ …" even for an unstable layer,
   and `layer == "forced convection"` says nothing about stability (items 5, 9); E5 — table parity needs rtol 1e-4 (item 6).
5. **qualitative (by design):** `jet_tke_budget` (only the production term follows from the similarity solution);
   `k_epsilon_channel` (optional demo outside the contract: converges, κ_model = 0.433; no grid study, no DNS comparison);
   the amplitude constants and half-widths of the free-shear table — no public table confirmed, `FREE_SHEAR_CONSTANTS` is
   empty and every example passes labelled illustrative constants (exponents are exact and tested).
6. **approximate (models against DNS, bands stated in the tests):** `spalding_uplus` (≤ 0.85 wall units, worst in the buffer
   layer), `stress_partition` (viscous share within 0.08 of the DNS dU⁺/dy⁺ for y/δ < 0.2), `channel_mixing_length`
   (∣ΔU⁺∣ ≤ 0.68, C_f ≤ +1.1 % for Re_τ ≥ 550; +5.5 % at Re_τ = 182, where no logarithmic layer exists).
7. **Data terms:** the Lee & Moser files carry no licence statement, only the request to cite the paper; a 281-row subset
   is committed with that citation. If the orchestrator prefers not to redistribute, delete the CSV and let
   `make_refs.py` fetch it (the five V5 tests then skip with "run reference/ch12/make_refs.py").
8. `skin_friction_zpg(law="power_fifth")`: the coefficient 0.074 is still unverified from ch09 (carried open item); only its
   algebra is tested here.
9. Secondary sources: the k–ε constants, Prandtl's law and the Businger–Dyer form were confirmed through search summaries of
   the cited pages (two primary pages returned HTTP 403 to the fetch tool); the Sreenivasan constant through the ATOMIX
   summary of the paper.

## Discrimination (planted wrong variants on a scratch copy of `fluidpy/`; the 114 default tests)
Baseline on the scratch copy: 114 passed, 3 deselected. Each mutant is one changed line; the run stops at the first failing
test (`-x`), so the test named is the first that notices, not the only one. **42 of 42 killed, 0 survived.**

| # | Wrong variant | Result · first failing test |
|---|---|---|
| M01 | central_moment keeps the mean (12.11) | killed · `test_ensemble_average_V1_matches_hand_loops_and_moment_identities` |
| M02 | time_average divides by 2Δt (12.2) | killed · `test_time_average_V1_sliding_window_matches_closed_form_and_limits` |
| M03 | spectrum_variance forgets the even half (12.22) | killed · `test_periodogram_V4_parseval_and_welch_cross_check` |
| M04 | taylor_microscale λ² = −2/a instead of −1/a (12.19) | killed · `test_integral_and_micro_scales_V1_closed_form_pairs` |
| M05 | reynolds_stress sign (12.30) | killed · `test_reynolds_stress_V1_covariance_sign_rotation_and_mean_stress` |
| M06 | mean stress μS̄ instead of 2μS̄ (12.30) | killed · `test_reynolds_stress_V1_covariance_sign_rotation_and_mean_stress` |
| M07 | RANS residual: stress divergence not divided by ρ0 (12.30) | killed · `test_slip06_rans_eddy_viscosity_V1_corrected_matches_12_30_and_printed_differs` |
| M08 | isotropic tensor with r f′ (2-D form) instead of (r/2) f′ (12.41) | killed · `test_isotropic_relations_V1_gaussian_f_three_routes_to_the_dissipation` |
| M09 | dissipation ν instead of ν/2 (12.42) | killed · `test_dissipation_rate_V1_strain_rotation_and_V4_spectral_identity_on_a_3d_field` |
| M10 | ε̄ = 15νu²/λ_f² (λ_f for λ_g) (12.43) | killed · `test_isotropic_relations_V1_gaussian_f_three_routes_to_the_dissipation` |
| M11 | production tensor with the transposed gradient (12.35) | killed · `test_reynolds_stress_V1_covariance_sign_rotation_and_mean_stress` |
| M12 | sign of the exchange term in the mean budget (12.46) | killed · `test_channel_energy_budget_V1_generic_budget_functions_agree_with_the_model_channel` |
| M13 | η = (ν³ε̄)^{1/4} (12.50) | killed · `test_kolmogorov_scales_V2_exponents_by_dimensional_analysis_D11` |
| M14 | two-sided spectrum not halved (12.55) | killed · `test_inertial_spectrum_V2_derivation_D12_one_group_and_the_18_over_55` |
| M15 | jet stress profile with the sign of analysis row 121 (12.63) | killed · `test_plane_jet_V2_derivation_D15_exponents_stress_profile_and_volume_flux` |
| M16 | volume flux with ∫F² instead of ∫F (12.68) | killed · `test_plane_jet_V4_momentum_flux_invariant_and_growing_volume_flux` |
| M17 | eddy-viscosity hypothesis sign (12.94) | killed · `test_slip06_rans_eddy_viscosity_V1_corrected_matches_12_30_and_printed_differs` |
| M18 | ν_T = C_μ ē^{3/2}/ε̄ (12.104) | killed · `test_k_epsilon_V2_derivation_D23_decay_exponent_and_log_layer_constant` |
| M19 | C_ε1 and C_ε2 swapped (12.105) | killed · `test_k_epsilon_decay_V1_closed_form_ivp_and_rk4_from_scratch` |
| M20 | flux Richardson sign (12.107) | killed · `test_flux_richardson_V1_hand_value_signs_and_degenerate_inputs` |
| M21 | Γ_a added instead of subtracted (12.108) | killed · `test_gradient_richardson_V1_in_situ_and_potential_routes_and_both_conventions` |
| M22 | Monin–Obukhov length without the minus (12.110) | killed · `test_stratified_tke_budget_V1_terms_and_optional_shear_argument` |
| M23 | dispersion rate without the factor 2 (12.117) | killed · `test_taylor_dispersion_V1_closed_forms_two_integral_forms_and_limits` |
| M24 | exponential dispersion without the 2 (12.119) | killed · `test_parity_V1_every_design_parity_expression_evaluates` |
| M25 | dispersion_regime with exclusive boundaries (item 3) | killed · `test_parity_V1_every_design_parity_expression_evaluates` |
| M26 | log law with log10 (12.88) | killed · `test_parity_V1_every_design_parity_expression_evaluates` |
| M27 | channel dP/dx = −τ0/h (half height) (12.90) | killed · `test_channel_stress_V2_derivation_D17_linear_total_stress_and_pressure_gradients` |
| M28 | total stress τ0(1 − y/h) (12.77) | killed · `test_channel_stress_V2_derivation_D17_linear_total_stress_and_pressure_gradients` |
| M29 | rough-wall law with log10 (12.93) | killed · `test_rough_wall_V1_derivation_D20_log_law_with_roughness_length` |
| M30 | Spalding bracket without the cubic term | killed · `test_fit_log_law_V1_exact_recovery_and_the_bias_of_a_composite_profile` |
| M31 | channel mixing length: τ⁺ missing inside the root (12.100) | killed · `test_channel_energy_budget_V4_pointwise_closure_and_integral_identity` |
| M32 | log-linear φ_m = 1 − βζ | killed · `test_monin_obukhov_V1_definition_sign_and_surface_layer_richardson_number_D25` |
| M33 | η/L ∼ Re^{−1/2} (12.51) | killed · `test_parity_V1_every_design_parity_expression_evaluates` |
| M34 | C1 = (55/18) C | killed · `test_inertial_spectrum_V2_derivation_D12_one_group_and_the_18_over_55` |
| M35 | Businger–Dyer ψ_m without the π/2 | killed · `test_surface_layer_wind_V5_businger_dyer_integral_with_and_without_the_z0_term` |
| M36 | eddy diffusivity (12.129) coded with the printed condition | killed · `test_slip12_eddy_diffusivity_V7_long_time_constant_and_printed_condition_fails` |
| M37 | parcel uv without the minus sign (N50) | killed · `test_parity_V1_every_design_parity_expression_evaluates` |
| M38 | Langevin position variance: wrong cross term | killed · `test_langevin_particles_stat_dispersion_matches_taylor_formula_within_5_standard_errors` |
| M39 | pipe dP/dx = −2τ0/d (12.91) | killed · `test_channel_stress_V2_derivation_D17_linear_total_stress_and_pressure_gradients` |
| M40 | Nagib–Chauhan coefficient 0.1663 → 0.1636 (12.92) | killed · `test_wall_correlations_V1_nagib_chauhan_zpg_fits_and_pipe_friction` |
| M41 | eddy diffusivity for the exponential correlation without Λ_t (12.127) | killed · `test_parity_V1_every_design_parity_expression_evaluates` |
| M42 | van Driest damping exp(−y⁺/A⁺) → exp(−2y⁺/A⁺) in the wall profile | killed · `test_explainer_tables_V1_json_reproduces_the_functions` |

Several mutants are first caught by the pinned closed forms of the parity test because it runs early; two of them were
re-run against later tests only to confirm an independent kill: M24 also fails `test_taylor_dispersion_V1_*`,
`test_langevin_particles_stat_*` and `test_slip13_*`; M42 (first caught by the cached table) also fails
`test_mixing_length_V1_van_driest_profile_by_independent_quadrature`.
