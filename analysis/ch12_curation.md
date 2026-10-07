# Chapter 12 — Turbulence: curation
(from `analysis/ch12.md` — 243 inventory rows (#1–#243), 129 equation numbers (12.1)–(12.129), 46 derivations in §2b (a-D1…a-D46),
17 load-bearing clusters (§3) and 15 suspected slips (§9); `book.yaml` ch12 (13 sections, 5 viz seeds) and `policy.tier_a_full_treatment: [12, 18]`,
`coverage: exhaustive`. Curated 2026-10-06, concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (P01–P279 +
P218a; **ch12 primers start at P280**), `notation.md`, `viz_patterns.md`, `analysis/ch11_curation.md` (format), the `interactive-viz` skill §4–§5.)

Counts: A 16 · B 210 · C 17 · RECAP 8 · SKIP 2 · derivations written out 28 (★ 10 ★★ 14 ★★★ 4) · demoted to statements 16

Tier words (parsed by `tools/nbkit.py`): **CORE 16 · NOTE 217 · RECAP 8 · SKIP 2** = 243 rows; **DERIVATION 28**.
Depth: A 16 (all CORE) · B 210 (204 NOTE + 6 RECAP) · C 17 (13 NOTE + 2 RECAP: R02, R07 + 2 SKIP: S01, S02).
**Reconciliation with analysis §2:** the 243 rows `#1…#243` each appear exactly once in §2 below (last column). A = 16 of 243 rows (6.6 %).

**Decisions that shape this chapter.**
1. **Sixteen A items from the analyst's seventeen clusters** (6.6 % of the rows — the chapter is wide, so the spine is kept thin).
   Merges: (a) the closure problem and the Reynolds-stress transport equation (12.35) (cluster 5) are B items of the RANS block
   **C04** — (12.35) is stated term by term and checked by sympy, and its only use downstream, the turbulent kinetic-energy budget,
   is derived directly in D10; (b) the isotropic correlations f, g (cluster 6) and the isotropic dissipation (cluster 7) are one
   block **C05** whose A row is (12.43); (c) the outer estimate ε̄ ~ ΔU³/L, the cascade and the Kolmogorov scales are one block **C07**;
   (d) Table 12.1 and the exponents of the other free shear flows are B items of the plane-jet block **C09**, with their derivation
   written out (D16) because the book prints only the table; (e) linear total stress, wall units and the sublayer are one block
   **C10**; the defect law, the overlap matching and the rough-wall law are inside the log-law block **C11**; (f) eddy viscosity and
   the mixing length are one block **C12** (the book-yaml seed treats them as one idea); (g) the limits, the random walk and the
   eddy diffusivity (12.127)–(12.129) are inside the Taylor block **C16**. §12.11 keeps **two** A blocks (C14 Richardson numbers,
   C15 Monin–Obukhov length and the surface-layer profile) because the reader is a climate-dynamics student. Next candidates if a
   reviewer asks for one more A (all B now): the closure problem (N60), Table 12.1 (N123), the eddy diffusivity (N213).
2. **SEEN rows are RECAP** (eight, R01–R08, inventory order): the Boussinesq set in flux form (ch04), the stress sign convention
   (ch02), the pressure-gradient / wall-stress balances (12.90)–(12.91) (ch08 D07, (8.8)), the stratified set-up with potential
   temperature and the two lapse-rate conventions (ch01, ch07, ch11), the gradient Richardson number (12.108) (ch11), the Gaussian
   spreading of a line vortex and σ² = 2Dt (ch01, ch03, ch08). Each reuses the earlier chapter's function; what the row adds
   (the turbulent mean, the in-situ-gradient form of Ri) is stated.
3. **SKIP is minimal** — S01 (exercises), S02 (literature). The analyst's 29 SKIP-class rows hold ideas (definitions, prose, steps)
   and are B or C items here. Exercise results the lesson uses are B items or D rows (listed in S01's pointer); exercise inputs
   and printed answers stay private.
4. **Derivations (D01–D28).** 30 of the 46 analysis derivations are written out as 28 D rows (three pairs merged: a-D8 + a-D9 in D06,
   a-D13 + a-D14 in D07, a-D31 + a-D32 in D19; a-D40 split between D24 and D25); none is added beyond the analyst's list
   (a-D39, the k–ε consistency checks, is "ours" there and is D23). Written out because **the book never writes them**: D02
   (Schwartz inequality), D04 (Taylor microscale, Ex. 12.9), D07 ((12.40)–(12.41), Ex. 12.17–12.18, ★★★), D08 ((12.43), Ex. 12.19,
   ★★★), D09 ((12.46), Ex. 12.15), D10 ((12.47), ★★★), D16 (exponents of Table 12.1, Ex. 12.25–12.28), D23 (k–ε consistency).
   Sixteen are demoted to statements (§4c) — among them four the book never writes but the lesson does not need step by step:
   a-D7 (periodogram, shown numerically three ways), a-D11 (mean scalar equation: the D06 moves again), a-D12 ((12.35): the D10
   moves with the second index kept free, checked by sympy), a-D34 (pipe friction law: checked against quadrature). The four ★★★
   are D07, D08 (C05 has no explainer — notebook only), D10 (shown in `turbulent_energy_budget`) and D14 (shown in
   `turbulent_jet_similarity`).
5. **Lapse-rate and temperature convention (analysis §9 C1; project rule).** Code computes with Kundu's Γ ≡ dT/dz (Γ_a ≈ −9.8 K/km,
   stable when dT/dz > Γ_a); every block that touches stability (N82, C14, C15, E9) shows the meteorological Γ ≡ −dT/dz
   (Γ_d ≈ +9.8 K/km, stable when Γ < Γ_d) beside it — two-row table, both inequalities in every legend, badge and slider trace.
   T̄ and T′ in the buoyancy terms of (12.30), (12.47), (12.106)–(12.112) are **potential** temperature; `gradient_richardson_thermal`
   takes the in-situ gradient plus Γ_a and a test asserts it equals the dθ/dz route. Heat-flux sign: ⟨wT′⟩ > 0 is upward
   (daytime, unstable) ⇒ Rf < 0 and L_M < 0, the same sign as the meteorological Obukhov length.
6. **A conventions block opens the chapter** (before C01; not a primer and not an inventory row): the overloaded symbols of analysis
   §9 C3 (κ von Kármán vs thermal diffusivity; k conductivity, wavenumber and the "k" of k–ε, which the book writes ē; e;
   λ / Λ; η; f, g, F, G; R_ij; δ; θ; α; τ; σ; S; N; L; U₀; h, d), the normalisations of §9 C4 (two-sided spectra in angular
   frequency with 1/2π in the forward transform; ⟨u²⟩ = one component; un-normalised skewness and kurtosis; which Reynolds number;
   y from the wall, h = full channel height) and the silent switches of §9 C5 (Boussinesq → constant density in §12.8–12.10 →
   Boussinesq; ensemble averages in the theory, time averages in every measurement). Each block repeats the row it needs.
7. **Printed slips (analysis §9 #1–#15) become named callouts "slip #k"** ("the page prints X, the correct form is Y", table in §8);
   `ch12.book_slips()` carries them. Planted wrong variants a test must fail: slip #1 (+5/3), #3 (source integral at y = 0),
   #5 (exponential similarity family), #6 (∂P/∂x_j), #12 (condition of (12.129)), #13 (swapped plume regimes).
8. **Climate hooks named where they are real:** C01 (30-year normals; eddy fluxes as deviations from a time or zonal mean),
   C03 (spectra of atmospheric and oceanic records; Ch. 13 wave and eddy spectra), C04 (eddy momentum and heat fluxes ⟨u′v′⟩,
   ⟨v′T′⟩ of the general circulation), C07 (millimetre Kolmogorov scale against kilometre grids — why sub-grid mixing is
   parameterised), C08 (−5/3 in the atmosphere and ocean; the −3 range of geostrophic turbulence in Ch. 13), C09 (plumes, wakes of
   islands and mountains), C11 (roughness length z₀, neutral drag coefficient, logarithmic wind profile), C12 (eddy viscosity of
   the Ekman layer, Ch. 13), C13 (turbulent-kinetic-energy closures of ocean and boundary-layer schemes), C14 (Ri-dependent mixing,
   nocturnal stable layer vs daytime convective layer), C15 (Monin–Obukhov similarity in bulk surface-flux formulas), C16
   (Richardson's 4/3 law, tracer and pollutant dispersion, eddy diffusivities).
9. **Book values stay private** (`tests/book_values_ch12.json`): Table 12.1's constants and half-widths, the per-flow (κ, B) pairs,
   the wake strength, the worked-example inputs and answers (Examples 12.1–12.3), the stirred-vessel and convection numbers, the
   fully-turbulent thresholds and all exercise data. Published values are cited to their authors (Sreenivasan 1995; Lee & Moser
   2015; Launder & Sharma 1974; Prandtl's friction law via McKeon et al. 2005; AMS Glossary for Businger–Dyer) and stored in
   `reference/ch12/`. κ and B are required keywords — no silent defaults. The fitted coefficients of the zero-pressure-gradient
   formulas are displayed equations of the book; if the orchestrator prefers to treat them as book-quoted numbers they move to the
   private JSON with Monkewitz et al. (2007) cited instead (a one-line decision, flagged in §8).

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
Order = book order, which is also the dependency order (analysis §3: A → B → C → D → E → F → G → H → I). Prerequisite notes:
C02 and C03 are needed again by C05 (same functions with r for τ and k₁ for ω) and by C16 (Λ_t); the turbulent kinetic-energy
budget C06 is placed before the cascade C07 because ε̄ ~ production is its steady state; C10 → C11 → C12 → C15 is one chain
(u_* → log law → mixing length → surface layer); N06 (Richardson's 4/3 law) is named in §12.2 and stated at the end of C16 where
both of its ingredients (C08, D_T) exist. RECAP rows are reminders inside the block named as their A parent. The conventions block
(decision 6) precedes C01.

**§12.1 Introduction** — no A item of its own: N01–N04 open the C01 block (the five marks of turbulence as a checklist, a solenoidal synthetic field against white noise, scope).

**§12.2 Historical Notes** — no A item of its own: N05 (a one-line timeline that links each name to its block) opens C01; N06 (Richardson's 4/3 law) is named there and stated at the end of C16.

**§12.3 Nomenclature and Statistics for Turbulent Flow**
- **C01** Ensemble averages and moments (12.1), the time average (12.2), the commutation rules (12.4)–(12.9): once you see an average as a plain sum over realizations, it commutes with every linear operation (sums, constants, derivatives, integrals) and fails on exactly one thing — a product — so ⟨ũṽ⟩ = Ū V̄ + ⟨uv⟩, and that leftover covariance is the whole subject of the chapter.
  - B: N01, N02, N07, N08, N09, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N20, N21, N22, N23
  - C: N03, N04, N05   ·   RECAP: –

**§12.4 Correlations and Spectra**
- **C02** Autocorrelation in lag form (12.17) and the integral time scale Λ_t (12.18): once the statistics do not depend on the time origin, the correlation is a function of the lag only, even in τ, equal to 1 at zero lag and decaying; its area Λ_t is the signal's memory, a record of length Δt holds about Δt/Λ_t independent samples, and the curvature at the origin gives the Taylor microscale.
  - B: N24, N25, N26, N27, N28, N29, N30, N31, N32, N33, N34
  - C: –   ·   RECAP: –
- **C03** The energy spectrum as the Fourier transform of the autocorrelation (12.20)–(12.22): once you transform the even function R₁₁(τ), you get a real, even S_e(ω) whose integral is the variance — a long memory means a narrow spectrum and S_e(0) = ⟨u²⟩Λ_t/π; with Taylor's frozen-field swap t ↔ x/U₀ the same pair gives wavenumber spectra.
  - B: N35, N36, N37, N38, N39, N40
  - C: –   ·   RECAP: –

**§12.5 Averaged Equations of Motion**
- **C04** The Reynolds-averaged momentum equation (12.30) and the Reynolds stress: once you write every field as mean + fluctuation and average the Boussinesq equations with the rules of C01, every term keeps its form except the nonlinear one, which leaves ∂⟨u_iu_j⟩/∂x_j — a stress carried by the fluctuations (negative ⟨uv⟩ in a positive shear, by the displaced-parcel argument) for which there is no equation: the closure problem.
  - B: N42, N43, N44, N45, N46, N47, N48, N49, N50, N51, N52, N53, N54, N56, N57, N58, N59, N60
  - C: N41, N55, N61   ·   RECAP: R01, R02

**§12.6 Homogeneous Isotropic Turbulence**
- **C05** Isotropic correlations and the dissipation (12.43): once the statistics have no preferred position or direction, the two-point tensor can only be built from δ_ij and r_ir_j, incompressibility ties the transverse function to the longitudinal one (g = f + r f′/2), and the nine-term dissipation collapses to ε̄ = 30ν⟨u²⟩/λ_f² = 15ν⟨u²⟩/λ_g² — one curve, one measurable gradient.
  - B: N62, N63, N64, N65, N66, N67, N68, N69, N70, N71, N72, N73, N74, N75, N76, N77
  - C: –   ·   RECAP: –

**§12.7 Turbulent Energy Cascade and Spectrum**
- **C06** The turbulent kinetic-energy budget (12.47) with the mean-flow budget (12.46): once you multiply the mean equation by U_i and the fluctuation equation by u_i and average, the same term −⟨u_iu_j⟩∂U_i/∂x_j appears in both with opposite signs — energy leaves the mean flow as shear production, is moved by transport terms, is added or removed by buoyancy and ends as viscous dissipation ε̄.
  - B: N78, N79, N80, N81, N82, N83
  - C: –   ·   RECAP: –
- **C07** Kolmogorov scales (12.50) and the scale separation (12.51): once the large eddies set the energy supply ε̄ ~ (ΔU)³/L regardless of viscosity, and the smallest eddies can know only ν and ε̄, dimensional analysis leaves one length η = (ν³/ε̄)^{1/4} and one velocity u_K = (νε̄)^{1/4} (their Reynolds number is 1), and η/L ~ Re_L^{−3/4} — the cascade gets longer, not weaker, as Re grows.
  - B: N84, N85, N86, N87, N88, N89, N90, N91, N92
  - C: –   ·   RECAP: –
- **C08** Kolmogorov's −5/3 law (12.54) and the universal spectrum (12.53): once the inertial range can depend on neither L nor ν, the only combination of ε̄ and k₁ with the units of S₁₁ is ε̄^{2/3}k₁^{−5/3}; in Kolmogorov scaling all high-Re spectra fall on one curve, with an inertial range that widens as Re_L^{3/4}.
  - B: N93, N94, N95, N96, N97
  - C: N98, N99   ·   RECAP: –

**§12.8 Free Turbulent Shear Flows**
- **C09** The self-similar plane jet (12.66): once the thin-layer equations give an invariant (the momentum flux J_s) and the profiles are assumed self-preserving, the similarity equation forces δ ∝ x and the invariant forces U_CL ∝ x^{−1/2}, so the volume flux grows as x^{1/2} (entrainment); the same two-line recipe with a different invariant gives every row of Table 12.1.
  - B: N100, N101, N102, N103, N104, N105, N106, N107, N108, N109, N110, N111, N112, N113, N114, N115, N116, N117, N118, N119, N120, N121, N122, N123, N124, N125, N126, N128, N130
  - C: N127, N129   ·   RECAP: –

**§12.9 Wall-Bounded Turbulent Shear Flows**
- **C10** The law of the wall (12.80) in wall units: once the mean momentum balance shows the total stress is linear across a channel and nearly constant near the wall, the wall stress alone sets a velocity u_* = (τ₀/ρ)^{1/2} and a length ν/u_*, so near any smooth wall U⁺ = f(y⁺), with U⁺ = y⁺ in the viscous sublayer where the stress is all viscous.
  - B: N131, N132, N133, N134, N135, N136, N137, N138, N139, N140, N141, N142
  - C: –   ·   RECAP: R03, R04
- **C11** The logarithmic law (12.88) by overlap matching: once the inner law (viscosity, no δ) and the defect law (δ, no viscosity) must both hold where y⁺ ≫ 1 and y/δ ≪ 1, the gradient y dU/dy can depend on neither variable — it is the constant u_*/κ — and integrating gives the logarithm in both forms; on a rough wall the same logarithm is anchored at the roughness length y₀.
  - B: N143, N144, N145, N146, N147, N148, N149, N150, N151, N152, N153, N154, N156, N157, N158, N159, N160, N161
  - C: N155   ·   RECAP: –

**§12.10 Turbulence Modeling**
- **C12** Eddy viscosity (12.94) and the mixing length (12.99)–(12.101): once the Reynolds stress is modelled as a viscosity ν_T ~ l_T u_T that belongs to the flow, the parcel argument of C04 gives −⟨uv⟩ = l_T²(dU/dy)², and with l_T = κy the constant-stress layer integrates from U⁺ = y⁺ to a logarithm — the model contains one length and returns C11.
  - B: N162, N163, N164, N165, N166, N167, N168, N169, N170, N171, N172, N173, N174, N181
  - C: –   ·   RECAP: –
- **C13** The k–ε model (12.103)–(12.105): once the velocity scale is taken from ē and the length scale from ē^{3/2}/ε̄, the eddy viscosity is C_μē²/ε̄ and two modelled transport equations close the mean flow; its constants are not free — C_ε2 is fixed by the decay of grid turbulence and the log layer requires κ² = √C_μ(C_ε2 − C_ε1)σ_ε.
  - B: N175, N176, N177, N178, N180
  - C: N179   ·   RECAP: –

**§12.11 Turbulence in a Stratified Medium**
- **C14** The flux and gradient Richardson numbers (12.107)–(12.109): once the energy budget of a horizontally uniform stratified shear flow is reduced to production, buoyancy and dissipation, their ratio Rf = buoyant destruction / shear production says whether turbulence can live (Rf ≳ ¼: it cannot), and the eddy coefficients turn it into the measurable Ri = (ν_T/κ_T)Rf — with temperature meaning potential temperature in both conventions of the lapse rate.
  - B: N182, N183, N184, N185
  - C: –   ·   RECAP: R05, R06
- **C15** The Monin–Obukhov length (12.110) and the stratified surface layer: once the log-layer values dU/dz = u_*/κz and −⟨uw⟩ = u_*² are put into Rf, the flux Richardson number is simply z/L_M — below ∣L_M∣ shear rules and the profile is logarithmic, above it buoyancy rules, and the wind profile bends as U = (u_*/κ)[ln(z/z₀) + 5z/L_M].
  - B: N186, N187, N188, N189, N190, N191, N192, N193, N194, N195
  - C: –   ·   RECAP: –

**§12.12 Taylor's Theory of Turbulent Dispersion**
- **C16** Taylor's dispersion formula (12.119): once the displacement is written as the time integral of the Lagrangian velocity, the mean-square spread depends only on the velocity variance and its autocorrelation — ballistic X_rms = u_rms t while the particle remembers its velocity (t ≪ Λ_t), a random walk X_rms = u_rms(2Λ_t t)^{1/2} once it has forgotten it, and an eddy diffusivity D_T that grows from ⟨u²⟩t to the constant ⟨u²⟩Λ_t.
  - B: N06, N196, N197, N198, N199, N200, N201, N202, N203, N204, N205, N206, N207, N208, N209, N210, N211, N212, N213, N214, N215, N216
  - C: N217   ·   RECAP: R07, R08

**§12.13 Concluding Remarks** — no A item of its own: N217 closes the C16 block with pointers to Ch. 13.

**End matter** — S01 (Exercises 12.1–12.38), S02 (literature cited).

## 2. Chapter map (depth) — every inventory row exactly once
Rows are in analysis §2 order (inventory row number in the last column). IDs: `C` = A/CORE (teaching order), `N` = B or C NOTE (inventory
order), `R` = RECAP, `S` = SKIP. No ASCII pipe inside cells (∣x∣ is written with U+2223). "slip #k" = analysis §9 suspected slip k;
"(§4c a-Dk)" = the analysis derivation is given as a statement; ⟨ ⟩ stands for the book's over-bar.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) | Row |
|---|---|---|---|---|---|---|---|
| N01 | Def: five marks of turbulence (fluctuations, nonlinearity, 3-D vorticity, dissipation, diffusivity) and the working definition; random waves are not turbulence | 12.1 | B | NOTE | C01 | stated as the opening checklist of C01 (a five-row table, each mark tied to the block that quantifies it); laminar vs turbulent (ch08), vortex stretching (ch05) and critical parameters (ch11) recalled in one line each | #1 |
| N02 | Thm: a turbulent velocity field obeys the conservation laws, a random vector field need not (∇·u = 0 vs white noise) | 12.1 | B | NOTE | C01 | stated with `TS.synthetic_solenoidal_field` vs `TS.white_noise_field` and `ch12.divergence_rms` (two panels, one number each); the synthetic field is labelled kinematic (no cascade) | #2 |
| N03 | Fig. 12.1 turbulent boundary layer: largest eddy l ≈ layer thickness δ, ragged turbulent / irrotational interface | 12.1 | C | NOTE | C01 | named with a sketch from the synthetic field; pointer: outer scale L in C07, entrainment in C09, δ in C10 | #3 |
| N04 | Scope: incompressible, no Coriolis force, 3-D fluctuations; nearly 2-D geostrophic turbulence deferred | 12.1 | C | NOTE | C01 | named; pointer: Ch. 13 (2-D turbulence sends energy upscale — the opposite of C07's cascade) | #4 |
| N05 | History: Reynolds (1883), Taylor (1915, 1921, 1935–36), Prandtl and von Kármán, Richardson, Kolmogorov and Obukhov (1941) | 12.2 | C | NOTE | C01 | named as a one-line timeline, each name linked to its block (C04, C16, C05, C12, C11, C07, C08, C15) | #5 |
| N06 | Thm: Richardson's four-thirds law K ∝ l^{4/3} (dimensionally K ~ ε̄^{1/3} l^{4/3}), consistent with the −5/3 spectrum | 12.2 | B | NOTE | C16 | named in the §12.2 timeline and stated where it is used, at the end of C16: `richardson_diffusivity`, exponents by `core.dimensional.solve_exponents`, one number for a 1 km patch | #6 |
| N07 | Def: realization, ensemble, ensemble average (over-bar), expected value ⟨ ⟩ as N → ∞; Reynolds decomposition first named | 12.3 | B | NOTE | C01 | stated with `TS.make_ensemble` (Ornstein–Uhlenbeck fluctuations on a prescribed mean) and `ensemble_average`; ch11's Lorenz sensitivity recalled as the reason for statistics instead of trajectories | #7 |
| C01 | Eq. (12.1) the m-th moment as an ensemble average ⟨u^m⟩ = lim_{N→∞} (1/N) Σ_n u(x, t: n)^m — with the time average (12.2), the commutation rules (12.4)–(12.9) and "the mean of a product is not the product of the means" | 12.3 | A | CORE | – | load-bearing: every equation of the chapter is an average; the rules (12.4)–(12.9) are the only moves Reynolds averaging uses, and the one rule that fails (products) is where the Reynolds stress comes from | #8 |
| N08 | Def: stationary in time; Eq. (12.2) time average over a window Δt | 12.3 | B | NOTE | C01 | stated with `TS.time_average` (centred sliding window; cumulative-sum primer); the window slider of E1 | #9 |
| N09 | Def: homogeneous (stationary in space); Eq. (12.3) volume average | 12.3 | B | NOTE | C01 | stated with `TS.volume_average` on the synthetic field; one number | #10 |
| N10 | Fig. 12.2 two sample records: stationary (constant mean) and non-stationary (decaying mean) | 12.3 | B | NOTE | C01 | stated with two `make_ensemble` members side by side | #11 |
| N11 | Which average is used in the field (time or space, the weather cannot be repeated); the window must satisfy t_c ≪ Δt ≪ drift time | 12.3 | B | NOTE | C01 | stated with the ergodicity primer and the climate hook (30-year normals; eddy fluxes as deviations from a time or zonal mean, Ch. 13); quantified by N23 and N32 | #12 |
| N12 | Eq. (12.4) averaging commutes with addition | 12.3 | B | NOTE | C01 | derived in D01 (A chain); residual from `TS.check_averaging_rules` | #13 |
| N13 | Eq. (12.5) averaging commutes with multiplication by a constant | 12.3 | B | NOTE | C01 | derived in D01; `check_averaging_rules` (entry "constant") | #14 |
| N14 | Eq. (12.6) averaging commutes with time differentiation | 12.3 | B | NOTE | C01 | derived in D01 (the chain the book writes out); ⚠️ exact for ensemble averages, only approximate for a finite time window; used by C16 in (12.115) | #15 |
| N15 | Eq. (12.7) averaging commutes with time integration | 12.3 | B | NOTE | C01 | derived in D01; used by C16 in (12.116) | #16 |
| N16 | Eq. (12.8) averaging commutes with spatial differentiation | 12.3 | B | NOTE | C01 | derived in D01; the rule behind mean continuity (D06) | #17 |
| N17 | Eq. (12.9) averaging commutes with spatial integration | 12.3 | B | NOTE | C01 | derived in D01; `check_averaging_rules` (entry "integral dx") | #18 |
| N18 | Rel: the average of an average is the average; the average of a product is not the product of the averages; ours: ⟨ũṽ⟩ = Ū V̄ + ⟨uv⟩ | 12.3 | B | NOTE | C01 | derived in the last steps of D01 with a two-member counter-example; `TS.product_average_split`; the seed of the Reynolds stress (C04) | #19 |
| N19 | Eq. (12.10) the mean (first moment) | 12.3 | B | NOTE | C01 | stated as the m = 1 case of (12.1); `ensemble_average` | #20 |
| N20 | Fig. 12.3 ensemble averages of 2, 4, 8 members of a decaying signal vs the expected value and a sliding time average; ours: scatter of the N-member mean ∝ N^{−1/2} | 12.3 | B | NOTE | C01 | stated with animation A1 and `TS.standard_error` (slope −½ on log–log axes; standard-error primer); the N slider of E1 | #21 |
| N21 | Eq. (12.11) central moments | 12.3 | B | NOTE | C01 | stated with `TS.central_moment` | #22 |
| N22 | Def: variance, skewness, kurtosis (the book's are the raw 2nd–4th central moments), standard deviation, rms; ours (Ex. 12.1): the moment identities | 12.3 | B | NOTE | C01 | stated with `TS.statistics(normalized=)` — both conventions in a two-row table with the Gaussian reference values 0 and 3 (histogram primer) | #23 |
| N23 | Ex. 12.1 time average of a decaying mean plus a cosine: the window multiplies the two parts by sinh(Δt/2τ)/(Δt/2τ) and sin(ωΔt/2)/(ωΔt/2); good window 1 ≪ ωΔt ≪ ωτ | 12.3 | B | NOTE | C01 | stated (§4c a-D2) with `time_average_exp_cos`, a sympy check of the window integral, slider figure F1 and E1; our own A, B, τ, ω | #24 |
| N24 | Eq. (12.12) two-point, two-time correlation tensor R_ij | 12.4 | B | NOTE | C02 | stated with `TS.correlation` (ensemble of pairs); ⚠️ R_ij is no longer the rotation tensor of ch02–ch03 | #25 |
| N25 | Def: uncorrelated, weakly / strongly correlated, anticorrelated; cross-correlation (i ≠ j) vs autocorrelation | 12.4 | B | NOTE | C02 | stated as a four-row table with scatter thumbnails from `TS.correlated_pair` | #26 |
| N26 | Eq. (12.13) autocorrelation function R_11 | 12.4 | B | NOTE | C02 | stated as the i = j case of (12.12) | #27 |
| N27 | Eq. (12.14) cross-correlation coefficient r_12 | 12.4 | B | NOTE | C02 | stated with `TS.correlation_coefficient`, cross-checked by `np.corrcoef` | #28 |
| N28 | Eq. (12.15) autocorrelation coefficient r_11 (equals 1 at zero separation) | 12.4 | B | NOTE | C02 | stated; the normalised curve every later figure plots | #29 |
| N29 | Eq. (12.16) Schwartz inequality ⇒ every correlation coefficient lies in [−1, 1] | 12.4 | B | NOTE | C02 | derived in D02 (the book states it without proof); ⚠️ trap: the bound is needed on ∣⟨uv⟩∣, the page prints no absolute value | #30 |
| C02 | Eq. (12.17) stationary process: the correlation depends on the lag only, R_11(τ) = ⟨u_1(t)u_1(t + τ)⟩ = R_11(−τ) — with the integral time scale Λ_t = ∫₀^∞ r_11 dτ (12.18) | 12.4 | A | CORE | – | load-bearing: the lag-only correlation and its area Λ_t (the memory time) decide how long to average (C01), set the zero-frequency spectrum (C03) and the scales of C05, and are the whole input of Taylor's dispersion theory (C16) | #31 |
| N30 | Fig. 12.4 three series (u, a delayed copy, v) and their auto- and cross-correlations: peak at the aligning lag | 12.4 | B | NOTE | C02 | stated with animation A2 (`TS.cross_correlation` of a delayed copy: the peak sits at the delay) and R_uv(τ) = R_vu(−τ) | #32 |
| N31 | Eq. (12.18) integral time scale Λ_t (the memory) | 12.4 | B | NOTE | C02 | stated with `TS.integral_scale` (to the first zero and over the whole record) and the equal-area rectangle; exact Λ_t = τ_c for r = e^{−τ/τ_c} | #33 |
| N32 | Def: correlation time t_c (first zero of r_11); number of independent samples in a record N ≈ Δt/t_c | 12.4 | B | NOTE | C02 | stated with `correlation_time`, `effective_samples`; closes the loop to N20 (the N of N^{−1/2} is the number of independent samples) | #34 |
| N33 | Eq. (12.19) temporal Taylor microscale λ_t² ≡ −2/r_11″(0) (osculating parabola, Ex. 12.9) | 12.4 | B | NOTE | C02 | derived in D04 (the book gives only the definition); `TS.taylor_microscale`; ⚠️ an Ornstein–Uhlenbeck signal has a cusp at τ = 0, so a smooth-spectrum signal is used | #35 |
| N34 | Fig. 12.5 autocorrelation coefficient with the equal-area rectangle (height 1, width Λ_t) and t_c | 12.4 | B | NOTE | C02 | stated: our figure with the rectangle, the first zero and the osculating parabola of N33 | #36 |
| C03 | Eq. (12.20) the energy spectrum S_e(ω) = (1/2π) ∫ R_11(τ) e^{−iωτ} dτ — the Fourier-transform pair (12.20)–(12.21) and the variance decomposition ⟨u_1²⟩ = ∫ S_e dω (12.22) | 12.4 | A | CORE | – | load-bearing: the spectrum sorts the variance by scale; (12.45), Kolmogorov's law (C08) and the temperature spectra are this pair with k for ω, and `TS.periodogram` becomes the core spectrum routine Ch. 13 reuses | #37 |
| N35 | Eq. (12.21) inverse transform R_11(τ) = ∫ S_e(ω) e^{+iωτ} dω | 12.4 | B | NOTE | C03 | derived in D05; `TS.correlation_from_spectrum` (round trip) | #38 |
| N36 | Eq. (12.22) the spectrum distributes the variance over frequency | 12.4 | B | NOTE | C03 | derived in D05 (τ = 0 in (12.21)); `TS.spectrum_variance` with the two-sided / one-sided factor 2 | #39 |
| N37 | Rel: zero-frequency value S_e(0) = ⟨u_1²⟩Λ_t/π | 12.4 | B | NOTE | C03 | derived in the last steps of D05; `integral_scale_from_spectrum` | #40 |
| N38 | Thm (Ex. 12.8): the same spectrum from a finite record — the windowed periodogram S_e = lim (1/2πT)∣∫u e^{−iωt}dt∣² | 12.4 | B | NOTE | C03 | stated (§4c a-D7) with `TS.periodogram` in the book's normalisation, a Parseval check and the three-routes figure (correlation → transform, periodogram, `scipy.signal.welch`); FFT-normalisation and Welch primers | #41 |
| N39 | Eq. (12.23) spatial correlation tensor R_ij(r) of a homogeneous field | 12.4 | B | NOTE | C03 | stated with `TS.spatial_correlation` on the synthetic field | #42 |
| N40 | Thm: Taylor's frozen-turbulence hypothesis x = U₀t; ∂u/∂x ≈ −(1/U)∂u/∂t (Ex. 12.11); ours: k_1 = ω/U₀, S_11(k_1) = U₀S_e(ω) | 12.4 | B | NOTE | C03 | stated with `TS.taylor_frozen`, `frequency_to_wavenumber_spectrum` (density change-of-variable primer) and `frozen_field_probe` (error grows with u_rms/U₀) | #43 |
| N41 | Why the mean is what matters (flight time ≫ fluctuation time; a range of scales too wide to resolve) | 12.5 | C | NOTE | C04 | named as the motivation of C04; pointer: scale separation in C07; the book's numbers private | #44 |
| N42 | Eq. (12.24) Reynolds decomposition of every field: ũ_i = U_i + u_i, p̃ = P + p, ρ̃ = ρ̄ + ρ′, T̃ = T̄ + T′ | 12.5 | B | NOTE | C04 | stated with `TS.reynolds_decompose` and the notation row of the conventions table (tilde = total, capital / over-bar = mean, lower case / prime = fluctuation) | #45 |
| N43 | Eq. (12.25) the means are expected values | 12.5 | B | NOTE | C04 | stated; asserted by `reynolds_decompose` | #46 |
| N44 | Eq. (12.26) fluctuations have zero mean | 12.5 | B | NOTE | C04 | stated; residual at round-off returned by `reynolds_decompose` | #47 |
| R01 | Rel: the Boussinesq set in flux form | 12.5 | B | RECAP | C04 | ch04 C13: (4.10), (4.86), (4.89) (D27–D28), `core.navier_stokes`; the flux-form move (adding ũ_i ∂ũ_j/∂x_j = 0) is step 1 of D06 | #48 |
| N45 | Eq. (12.27) mean continuity ∂U_i/∂x_i = 0 | 12.5 | B | NOTE | C04 | derived in D06 (A chain); `mean_divergence` | #49 |
| N46 | Eq. (12.28) the fluctuation is divergence-free too | 12.5 | B | NOTE | C04 | derived in D06 (total minus mean) | #50 |
| N47 | Eq. (12.29) the decomposition substituted into the momentum equation | 12.5 | B | NOTE | C04 | a middle line of D06 | #51 |
| N48 | Rel: term-by-term averages; only the product term leaves something new, ∂⟨u_iu_j⟩/∂x_j | 12.5 | B | NOTE | C04 | derived in D06 (which cross terms vanish and why); numeric check with `product_average_split` | #52 |
| C04 | Eq. (12.30) Reynolds-averaged momentum equation ∂U_i/∂t + U_j ∂U_i/∂x_j = −g[1 − α(T̄ − T₀)]δ_i3 + (1/ρ₀) ∂τ̄_ij/∂x_j with τ̄_ij = −Pδ_ij + 2μS̄_ij − ρ₀⟨u_iu_j⟩ — the Reynolds stress and the closure problem | 12.5 | A | CORE | – | load-bearing: the mean-flow equation of every later section (jets, walls, models, stratified flow) and of Ch. 13's eddy-viscosity Ekman layers; its one new term −ρ₀⟨u_iu_j⟩ is what the rest of the chapter measures, scales or models | #53 |
| N49 | Def: Reynolds stress tensor −ρ₀⟨u_iu_j⟩: symmetric, six components, normal and shear stresses; ≫ viscous stress except very near a wall | 12.5 | B | NOTE | C04 | stated with `TS.reynolds_stress` and `anisotropy_tensor` (covariance-matrix primer; ch02 principal axes) | #54 |
| N50 | Fig. 12.6 why ⟨uv⟩ < 0 when dU/dy > 0: a parcel displaced upward keeps its lower mean speed | 12.5 | B | NOTE | C04 | stated with `displaced_parcel_uv` (Monte-Carlo parcels), animation A3 and E3; the argument the mixing length reuses (C12) | #55 |
| N51 | Rel: Reynolds stress = mean momentum flux carried by the fluctuations, ρ₀⟨(U + u)v⟩ = ρ₀⟨uv⟩ | 12.5 | B | NOTE | C04 | stated in three lines; the flux is counted directly across a plane in the parcel script | #56 |
| R02 | Fig. 12.7 positive directions of the Reynolds shear stresses on an element | 12.5 | C | RECAP | C04 | ch02 §2.4 stress sign convention (Fig. 2.4, stress cube); drawn inside the parcel figure | #57 |
| N52 | Eq. (12.31) mean temperature equation with the divergence of ⟨u_jT′⟩ | 12.5 | B | NOTE | C04 | stated (§4c a-D10): the D06 moves applied to (4.89); `rans_sympy` (heat part); ⚠️ κ is the thermal diffusivity here (`kappa_th`) | #58 |
| N53 | Eq. (12.32) the same in terms of the mean heat flux Q_j = −k ∂T̄/∂x_j + ρ₀C_p⟨u_jT′⟩ | 12.5 | B | NOTE | C04 | stated with `mean_heat_flux`; Fourier's law (ch01 (1.2)) recalled | #59 |
| N54 | Def: turbulent heat flux ρ₀C_p⟨u_jT′⟩; upward over a heated surface | 12.5 | B | NOTE | C04 | stated with `turbulent_heat_flux` and the W/m² ↔ K m/s conversion (primer); its sign feeds C14 and C15 | #60 |
| N55 | Def: passive (conserved) scalar; mixture density from the volume fraction | 12.5 | C | NOTE | C04 | named with `mixture_density`; pointer: jet dilution N117–N121 and Example 12.2 (N126) | #61 |
| N56 | Eq. (12.33) instantaneous scalar conservation | 12.5 | B | NOTE | C04 | stated (flux form, Fick's law of ch01 (1.1) recalled) | #62 |
| N57 | Eq. (12.34) mean scalar equation with the turbulent flux ⟨u_jY′⟩ | 12.5 | B | NOTE | C04 | stated (§4c a-D11; Ex. 12.12): the D06 moves at constant ρ_m; `mean_scalar_flux` | #63 |
| N58 | Def: Reynolds averaging, RANS equations = (12.27), (12.30), (12.32), (12.34); constant-density form | 12.5 | B | NOTE | C04 | stated with the constant-density equation written out — the form used from §12.8 to §12.10 (P = deviation from hydrostatic; ⚠️ a silent switch of assumptions) | #64 |
| N59 | Eq. (12.35) transport equation for the Reynolds stress ⟨u_iu_j⟩ | 12.5 | B | NOTE | C04 | stated (§4c a-D12) term by term (advection, turbulent transport, production, pressure–velocity, dissipation, viscous diffusion, buoyancy); `reynolds_stress_budget_sympy` checks one component and that half its trace is (12.47) (D10) | #65 |
| N60 | Def: closure problem — first moments need second, second need third; three responses: RANS models, DNS, LES | 12.5 | B | NOTE | C04 | stated with `closure_count` (4 equations for 10 unknowns at first order; the triple correlations at second) and the three responses (models → C12, C13; DNS cost → C07; LES → ch10 pointer) | #66 |
| N61 | With Reynolds stresses the mean momentum equation has no Bernoulli integral | 12.5 | C | NOTE | C04 | named; pointer: ch04 `which_bernoulli` | #67 |
| N62 | Def: homogeneous and isotropic turbulence; grid turbulence as its approximate realisation | 12.6 | B | NOTE | C05 | stated as the opening of C05 | #68 |
| N63 | Fig. 12.8 scatter of (u, v) samples: round cloud ⇒ ⟨uv⟩ = 0; cloud along v = −u ⇒ ⟨uv⟩ < 0 | 12.6 | B | NOTE | C05 | stated with `TS.correlated_pair` and the principal axes of the covariance (covariance-ellipse primer) | #69 |
| N64 | Eq. (12.36) homogeneity and isotropy: no gradients of statistics, equal normal stresses, equal same-direction gradient moments | 12.6 | B | NOTE | C05 | stated with `isotropy_report` on a 3-D synthetic field (FAST: 48³, cached) | #70 |
| N65 | Eq. (12.37) equal cross-direction gradient moments | 12.6 | B | NOTE | C05 | stated; `isotropy_report` | #71 |
| N66 | Eq. (12.38) longitudinal and transverse correlation coefficients f(r), g(r) | 12.6 | B | NOTE | C05 | stated with `TS.longitudinal_transverse_correlation`; the figure of f and g from a synthetic field | #72 |
| N67 | Eq. (12.39) integral scales Λ_f, Λ_g and Taylor microscales λ_f, λ_g | 12.6 | B | NOTE | C05 | stated; the same `integral_scale` and `taylor_microscale` as C02 with r for τ | #73 |
| N68 | Eq. (12.40) the most general isotropic two-point tensor R_ij = F(r) r_i r_j + G(r) δ_ij | 12.6 | B | NOTE | C05 | derived in D07 (Ex. 12.17; the book quotes F and G); `isotropic_correlation_tensor` | #74 |
| N69 | Eq. (12.41) incompressible isotropic tensor in terms of f alone | 12.6 | B | NOTE | C05 | derived in D07 (Ex. 12.18); `isotropic_tensor_divergence_sympy` (must vanish); slip #14 (the exercise cites (12.39) for (12.40)) | #75 |
| N70 | Rel: g = f + (r/2) f′; Λ_g = Λ_f/2; λ_g = λ_f/√2 | 12.6 | B | NOTE | C05 | derived in the last steps of D07; `transverse_from_longitudinal`, `isotropic_scales`; ⚠️ three-dimensional results (a 2-D field obeys g = d(rf)/dr) | #76 |
| N71 | Fig. 12.9 geometry of the longitudinal and transverse correlations | 12.6 | B | NOTE | C05 | stated with arrows drawn on the synthetic field above the f, g curves | #77 |
| N72 | Def: turbulent kinetic energy per unit mass from the trace, R_ii(0) = 2ē | 12.6 | B | NOTE | C05 | stated with `TS.turbulent_kinetic_energy`; ⚠️ ⟨u²⟩ in §12.6 is one component, ē = (3/2)⟨u²⟩; `e` is not ch01's internal energy | #78 |
| N73 | Eq. (12.42) mean dissipation rate of the fluctuations ε̄ = (ν/2)⟨(∂u_i/∂x_j + ∂u_j/∂x_i)²⟩ | 12.6 | B | NOTE | C05 | stated as the average of ch04's (4.58); `dissipation_rate`; the starting line of D08 | #79 |
| C05 | Eq. (12.43) isotropic dissipation ε̄ = 6ν{⟨(∂u_1/∂x_1)²⟩ + ⟨(∂u_1/∂x_2)²⟩ + ⟨(∂u_1/∂x_2)(∂u_2/∂x_1)⟩} = −15ν⟨u²⟩f″(0) = 30ν⟨u²⟩/λ_f² = 15ν⟨u²⟩/λ_g² — built on the isotropic correlations f, g (12.38)–(12.41) | 12.6 | A | CORE | – | key method: symmetry collapses a nine-component two-point tensor to one function f(r) and the dissipation to one measurable gradient; it fixes what the Taylor microscale means and feeds R_λ, (12.52) and the model spectra of C08 | #80 |
| N74 | Rel: gradient moments from the correlation tensor, ⟨∂_k u_i ∂_l u_j⟩ = −∂²R_ij/∂r_k∂r_l at r = 0 (Ex. 12.19) | 12.6 | B | NOTE | C05 | derived in D08 (the three moments are 2 : 4 : −1 in units of ⟨u²⟩/λ_f²); `gradient_moments_isotropic_sympy` | #81 |
| N75 | Eq. (12.44) Taylor-scale Reynolds number R_λ | 12.6 | B | NOTE | C05 | stated with `taylor_reynolds_number` (which microscale: a factor √2); the book's fully-turbulent threshold private | #82 |
| N76 | Eq. (12.45) one-dimensional wavenumber spectrum S_11(k_1) | 12.6 | B | NOTE | C05 | stated: C03's pair with (r_1, k_1) for (τ, ω); units m³/s² (spectral-density units primer); the input of C08 | #83 |
| N77 | The measurement recipe: time record → frozen-turbulence conversion → R_11 → S_11 | 12.6 | B | NOTE | C05 | stated as a four-step pipeline cell (`taylor_frozen` → `autocorrelation` → `spectrum_from_correlation`) against `periodogram` | #84 |
| N78 | Eq. (12.46) kinetic-energy budget of the mean flow, Ē = ½U_i² | 12.7 | B | NOTE | C06 | derived in D09 (Ex. 12.15; the book gives the result with labelled terms); `mean_energy_budget` | #85 |
| N79 | Rel: meaning of the terms — divergence terms only move energy; ⟨u_iu_j⟩∂U_i/∂x_j = ⟨u_iu_j⟩S̄_ij | 12.7 | B | NOTE | C06 | stated: Gauss's theorem (ch02) for the transport terms, the symmetric contraction of ch02 §2.10 for the exchange term | #86 |
| N80 | Rel: direct viscous dissipation of the mean flow is ~ 1/Re of the shear production | 12.7 | B | NOTE | C06 | stated (§4c a-D18) with `mean_to_turbulent_dissipation_ratio` | #87 |
| C06 | Eq. (12.47) turbulent kinetic-energy budget ∂ē/∂t + U_j ∂ē/∂x_j = transport − 2ν⟨S′_ijS′_ij⟩ − ⟨u_iu_j⟩ ∂U_i/∂x_j + gα⟨u_3T′⟩ — with the mean-flow budget (12.46) and the shear-production term they share | 12.7 | A | CORE | – | load-bearing: the production −⟨u_iu_j⟩∂U_i/∂x_j leaves the mean flow and enters the turbulence; the budget starts the cascade argument (C07), is the equation k–ε models (C13) and is what the Richardson numbers divide (C14); the turbulent twin of ch11's (11.88) | #88 |
| N81 | Def: shear production −⟨u_iu_j⟩∂U_i/∂x_j (opposite signs in (12.46) and (12.47)); dissipation ε̄ = 2ν⟨S′_ijS′_ij⟩ > 0; buoyant production / destruction gα⟨u_3T′⟩ | 12.7 | B | NOTE | C06 | stated with `shear_production`, `buoyant_production`; the paired bars of E5; ch11 `disturbance_energy_budget` recalled | #89 |
| N82 | Fig. 12.10 unstable layer: upward heat flux, mean potential energy falls, turbulent energy rises | 12.7 | B | NOTE | C06 | stated with `mixing_potential_energy_change`; ⚠️ T here is potential temperature — both lapse-rate conventions in the caption | #90 |
| N83 | Rel: isotropic turbulence has no shear production | 12.7 | B | NOTE | C06 | stated (§4c a-D18) in one line (⟨u_iu_j⟩ ∝ δ_ij contracts to ∇·U = 0); a test of `shear_production` | #91 |
| N84 | Fig. 12.11 outer scales: L = cross-stream extent of the velocity difference ΔU | 12.7 | B | NOTE | C07 | stated with a schematic for a free layer and a duct | #92 |
| N85 | Eq. (12.48) energy input rate to the large eddies ~ (ΔU)³/L | 12.7 | B | NOTE | C07 | derived in D11 (A chain); `dissipation_outer_scaling` | #93 |
| N86 | Eq. (12.49) stationary turbulence: ε̄ ~ (ΔU)³/L, independent of ν | 12.7 | B | NOTE | C07 | derived in D11; the viscosity sets where, not how much, energy is dissipated | #94 |
| N87 | Def: the (Richardson) cascade — each tier strained by the next larger one, inviscid while u′l′/ν ≫ 1; ours: u′(l′) ~ (ε̄l′)^{1/3} | 12.7 | B | NOTE | C07 | stated with `cascade_tiers` (sizes, velocities, turnover times, eddy Reynolds numbers down to η) and animation A4; ch05 vortex stretching recalled | #95 |
| C07 | Eq. (12.50) Kolmogorov scales η = (ν³/ε̄)^{1/4}, u_K = (νε̄)^{1/4} — with the outer estimate ε̄ ~ (ΔU)³/L (12.48)–(12.49) and the scale separation η/L ~ Re_L^{−3/4} (12.51) | 12.7 | A | CORE | – | load-bearing: the two ends of the cascade and a separation of scales that grows with the Reynolds number — the reason turbulence must be parameterised in every ocean and atmosphere model (millimetres against kilometres) | #96 |
| N88 | Rel: the Kolmogorov-scale Reynolds number ηu_K/ν is one | 12.7 | B | NOTE | C07 | stated; asserted in `kolmogorov_scales` | #97 |
| N89 | Eq. (12.51) η/L ~ Re_L^{−3/4} | 12.7 | B | NOTE | C07 | derived in D11; `scale_separation`, `dns_grid_points` (ours: (L/η)³ ∝ Re^{9/4}) | #98 |
| N90 | Value: η of order millimetres in the ocean and atmosphere; a stirred-vessel estimate | 12.7 | B | NOTE | C07 | stated with our own cases from `scale_table` (atmospheric boundary layer, ocean thermocline, a kitchen mixer); the book's worked number private | #99 |
| N91 | Eq. (12.52) Taylor microscale between η and L: λ_T/L ∝ Re_L^{−1/2} | 12.7 | B | NOTE | C07 | stated (§4c a-D20) with `scale_separation` (the factor 15 or 30 of (12.43) exposed); why λ_T is not the size of the dissipating eddies | #100 |
| N92 | Rel: R_λ ~ Re_L^{1/2}; ordering η < λ_T < Λ < L at high Reynolds number | 12.7 | B | NOTE | C07 | stated with `scale_ordering`; the book's thresholds private | #101 |
| N93 | Eq. (12.53) universal (Kolmogorov-scaled) form of the spectrum at high wavenumber, S_11/(u_K²η) = Φ(k_1η) | 12.7 | B | NOTE | C08 | derived in D12 (A chain); `kolmogorov_normalize_spectrum`, Π groups from `core.dimensional.pi_groups` | #102 |
| C08 | Eq. (12.54) Kolmogorov's inertial-subrange law S_11(k_1) = C_1 ε̄^{2/3} k_1^{−5/3} for 2π/L ≪ k_1 ≪ 2π/η (⚠️ slip #1: the page prints k_1^{+5/3}) — with the universal form (12.53) | 12.7 | A | CORE | – | load-bearing: the chapter's signature result; the same one-group dimensional argument gives the temperature spectra (12.113)–(12.114) and Richardson's 4/3 law, and Ch. 13's geostrophic-turbulence spectra are read against it | #103 |
| N94 | Eq. (12.55) the double-sided normalisation behind the quoted one-dimensional constant | 12.7 | B | NOTE | C08 | stated with the one- vs two-sided factor 2 in `inertial_spectrum_1d` | #104 |
| N95 | Rel: three-dimensional spectrum S(K) = C ε̄^{2/3} K^{−5/3}; ours: one-sided C_1 = (18/55) C | 12.7 | B | NOTE | C08 | stated with `kolmogorov_constants` (the factor checked by `one_dimensional_from_3d` on a pure power law); Sreenivasan (1995) cited for the measured constant | #105 |
| N96 | Fig. 12.12 measured one-dimensional spectra collapse in Kolmogorov scaling; −5/3 line; roll-off near the dissipation scale | 12.7 | B | NOTE | C08 | stated with our model spectra at several R_λ in the same scaling (slider figure F3, E4); the book's read-off numbers private | #106 |
| N97 | Rel: a model spectrum for the whole range (Pao's exponential roll-off; Pope's energy-range factor) | 12.7 | B | NOTE | C08 | stated with `model_spectrum` (forms cited) and its two integral checks (∫E dK = ē, 2ν∫K²E dK = ε̄); log-grid quadrature primer | #107 |
| N98 | Value: first confirmation of the −5/3 law in a tidal channel at very large Reynolds number | 12.7 | C | NOTE | C08 | named; pointer: our Fig. 12.23 analogue (N195); the number private | #108 |
| N99 | Universality of the small scales motivates closure models and large-eddy simulation | 12.7 | C | NOTE | C08 | named; pointer: C12, C13 and ch10 (LES) | #109 |
| N100 | Fig. 12.13 three generic free shear flows — jet, wake, shear layer (a plume is a buoyant jet); Def: free vs wall-bounded shear flow | 12.8 | B | NOTE | C09 | stated with mean profiles from `free_shear_flow` at several stations; laminar jet and wake (ch09, ch04) recalled | #110 |
| N101 | Def: entrainment; self-preservation (profiles at different x collapse in local scales) | 12.8 | B | NOTE | C09 | stated as the two opening definitions of C09; the similarity collapse of ch08 §8.4 and ch09 recalled | #111 |
| N102 | Eq. (12.56) self-preserving form of the mean velocity U = U_CL(x) F(y/δ(x)) | 12.8 | B | NOTE | C09 | stated as the starting ansatz of D14; `plane_jet_mean_velocity` | #112 |
| N103 | Eq. (12.57) self-preserving form of the Reynolds shear stress −⟨uv⟩ = Ψ(x) G(y/δ(x)) | 12.8 | B | NOTE | C09 | stated with D14; `plane_jet_reynolds_stress`; ⚠️ F, G are not the tensor functions of (12.40) | #113 |
| N104 | Eq. (12.58) two-dimensional mean continuity | 12.8 | B | NOTE | C09 | stated; `plane_jet_cross_velocity` | #114 |
| N105 | Eq. (12.59) stream-wise mean momentum | 12.8 | B | NOTE | C09 | stated as the starting line of D13; `rans_2d_residual` | #115 |
| N106 | Eq. (12.60) cross-stream mean momentum | 12.8 | B | NOTE | C09 | stated; `rans_2d_residual` | #116 |
| N107 | Eq. (12.61) thin-layer form at high Reynolds number with no imposed pressure gradient | 12.8 | B | NOTE | C09 | derived in D13 (A chain; sizes of kept and dropped terms from `thin_shear_layer_terms`) | #117 |
| N108 | Rel: conservative form and its cross-stream integral | 12.8 | B | NOTE | C09 | derived in D13 | #118 |
| N109 | Eq. (12.62) the momentum flux per unit span J_s is the same at every station | 12.8 | B | NOTE | C09 | derived in D13; `jet_momentum_flux_per_span` flat in x (the invariant bars of E6); laminar J of ch09 recalled | #119 |
| N110 | Rel: V eliminated by continuity; the similarity forms substituted | 12.8 | B | NOTE | C09 | the first moves of D14 | #120 |
| N111 | Eq. (12.63) similarity equation of the plane jet (three coefficient brackets) | 12.8 | B | NOTE | C09 | derived in D14 (★★★; the book writes "somewhat tedious"); `plane_jet_similarity_sympy`, `plane_jet_stress_profile` | #121 |
| N112 | Eq. (12.64) simplest similarity condition: each coefficient is a constant | 12.8 | B | NOTE | C09 | derived in D15 | #122 |
| N113 | Rel: linear growth δ = (C_2 − C_1)(x − x_o) and the virtual origin | 12.8 | B | NOTE | C09 | derived in D15; `virtual_origin_fit` | #123 |
| N114 | Eq. (12.65) the invariant fixes the decay exponent: 2γ + 1 = 0 | 12.8 | B | NOTE | C09 | derived in D15 | #124 |
| C09 | Eq. (12.66) far field of the plane turbulent jet U(x, y) = C_5 (J_s/ρ)^{1/2} x^{−1/2} F(y/x) — from the invariant J_s (12.62) and the similarity equation (12.63); δ ∝ x, volume flux ∝ x^{1/2} (12.68) | 12.8 | A | CORE | – | key method: an invariant plus self-preservation gives the exponents without solving for the turbulence — the recipe behind every row of Table 12.1 (jets, wakes, plumes, shear layers) and behind entrainment; plumes return in Ch. 13 | #125 |
| N115 | Eq. (12.67) far-field Reynolds shear stress | 12.8 | B | NOTE | C09 | stated; ⚠️ slip #4: the undetermined constants are C_3 and C_5 (the text names C_3, C_4) | #126 |
| N116 | Eq. (12.68) the volume flux per unit span grows as x^{1/2} (entrainment) | 12.8 | B | NOTE | C09 | derived in the last steps of D15; `plane_jet_volume_flux`, `plane_jet_entrainment_velocity` | #127 |
| N117 | Eq. (12.69) similarity form of the mean mass fraction | 12.8 | B | NOTE | C09 | stated; `plane_jet_mass_fraction` | #128 |
| N118 | Eq. (12.70) conservation of slot fluid | 12.8 | B | NOTE | C09 | stated (§4c a-D25); `scalar_flux_per_span`; ⚠️ slip #3: the source integral is evaluated at x = 0 | #129 |
| N119 | Eq. (12.71) far-field mean mass fraction ∝ x^{−1/2} | 12.8 | B | NOTE | C09 | stated (§4c a-D25) | #130 |
| N120 | Eq. (12.72) uniform slot exit: velocity in nozzle variables | 12.8 | B | NOTE | C09 | stated; `free_shear_centerline` | #131 |
| N121 | Eq. (12.73) uniform slot exit: mass fraction in nozzle variables | 12.8 | B | NOTE | C09 | stated; `free_shear_centerline` (scalar part) | #132 |
| N122 | Rel: Gaussian fit of the profiles, specified by the half-width ξ_{1/2} | 12.8 | B | NOTE | C09 | stated with `gaussian_profile`, `profile_integrals` (Gaussian integral P194) | #133 |
| N123 | Table 12.1 self-similar far-field laws for seven flows (planar / round jet, planar / round plume, shear layer, planar / round wake) | 12.8 | B | NOTE | C09 | stated as our own table of exact exponents (`free_shear_exponents`) and amplitude scales (`free_shear_flow`); the book's constants and half-widths stay private, examples pass labelled illustrative constants; the flow modes of E6 | #134 |
| N124 | Rel (ours, Ex. 12.25–12.28): the exponents of Table 12.1 from each flow's invariant and growth law | 12.8 | B | NOTE | C09 | derived in D16 (the book prints only the table); `local_reynolds_number_exponent` | #135 |
| N125 | Eq. (12.74) more general similarity: the three coefficients need only share their x-dependence | 12.8 | B | NOTE | C09 | stated (§4c a-D26) with `general_similarity_check`; ⚠️ slip #5: the quoted exponential family zeroes the middle coefficient and breaks (12.62) | #136 |
| N126 | Ex. 12.2 fuel jet into air: distance to the stoichiometric mixture on the axis | 12.8 | B | NOTE | C09 | stated with `stoichiometric_mass_fraction`, `round_jet_distance_for_mass_fraction` for our own gas, nozzle and speed; the book's inputs and answers private | #137 |
| N127 | Fig. 12.14 profiles across a plane jet of ē, the normal stresses and the shear stress | 12.8 | C | NOTE | C09 | named with a qualitative model panel (shear stress from `plane_jet_stress_profile`); pointer: N128 and C06 | #138 |
| N128 | Eq. (12.75) turbulent kinetic-energy budget of the jet in the thin-layer approximation | 12.8 | B | NOTE | C09 | stated (§4c a-D28) with `jet_tke_budget` (labelled qualitative); ⚠️ slip #10: the triple correlation is ½⟨u_i²u_j⟩ throughout | #139 |
| N129 | Fig. 12.15 terms of the jet's energy budget across the layer | 12.8 | C | NOTE | C09 | named; pointer: the three balances (axis, mid-layer, edge) read from `jet_tke_budget` | #140 |
| N130 | Rel (ours): a constant eddy viscosity across the jet reproduces the laminar sech² shape with turbulent exponents | 12.8 | B | NOTE | C09 | stated with `plane_jet_eddy_viscosity_profile` against the Gaussian; laminar ghost from `core.jets`; forward pointer to C12 | #141 |
| N131 | Def: wall-bounded turbulence has two length scales (viscous l_ν near the wall, thickness δ) and no Reynolds-number independence on smooth walls | 12.9 | B | NOTE | C10 | stated as the opening of C10; ch09's drag curve recalled | #142 |
| N132 | Fig. 12.16 mean turbulent profiles are blunter, with a larger wall slope, than laminar ones | 12.9 | B | NOTE | C10 | stated with `WT.composite_profile` against the parabola (`core.laminar.channel_flow`) and Blasius | #143 |
| N133 | Eq. (12.76) mean momentum in fully developed channel flow; total stress τ̄ = μ ∂U/∂y − ρ₀⟨uv⟩ | 12.9 | B | NOTE | C10 | derived in D17 (A chain); `WT.total_stress`, `channel_momentum_residual` | #144 |
| N134 | Rel: wall-normal balance integrated from the wall | 12.9 | B | NOTE | C10 | a step of D17 | #145 |
| N135 | Eq. (12.77) the stream-wise pressure gradient is the same at every y | 12.9 | B | NOTE | C10 | a step of D17 | #146 |
| N136 | Thm: the total stress varies linearly across a channel (and a pipe); ours: τ̄(y) = τ₀(1 − 2y/h) | 12.9 | B | NOTE | C10 | derived in D17; `WT.channel_total_stress`, `stress_partition`; ⚠️ h is the full height | #147 |
| N137 | Fig. 12.17 total stress across a channel and across a zero-pressure-gradient boundary layer (constant-stress layer near the wall) | 12.9 | B | NOTE | C10 | stated with `stress_partition` (viscous vs Reynolds part) and the boundary-layer stress from a composite profile | #148 |
| N138 | Eq. (12.78) stream-wise mean momentum of a flat-plate boundary layer | 12.9 | B | NOTE | C10 | stated; `boundary_layer_stress_from_profile`; laminar form (9.9) recalled | #149 |
| N139 | Def: inner layer, outer layer, overlap region | 12.9 | B | NOTE | C10 | stated with `WT.layer_name`; matching primer | #150 |
| N140 | Eq. (12.79) inner-layer dependence U = U(ρ, τ₀, ν, y) | 12.9 | B | NOTE | C10 | the starting line of D18; `law_of_the_wall_groups` | #151 |
| C10 | Eq. (12.80) law of the wall U⁺ ≡ U/u_* = f(y⁺), y⁺ = yu_*/ν — with the friction velocity u_*² = τ₀/ρ (12.81), the linear total stress (12.76)–(12.77) and the viscous sublayer U⁺ = y⁺ (12.82) | 12.9 | A | CORE | – | load-bearing: wall units collapse every smooth-wall flow near the wall; u_* is the velocity scale of the log law (C11), the mixing length (C12), the Monin–Obukhov length (C15) and Ch. 13's surface and Ekman layers | #152 |
| N141 | Eq. (12.81) friction velocity u_* (and the viscous wall unit l_ν = ν/u_*) | 12.9 | B | NOTE | C10 | derived in D18; `WT.friction_velocity`, `viscous_length`, `friction_reynolds_number` | #153 |
| N142 | Eq. (12.82) viscous sublayer U⁺ = y⁺ (to about y⁺ ≈ 5) | 12.9 | B | NOTE | C10 | derived in the last steps of D18; `WT.viscous_sublayer` | #154 |
| N143 | Eq. (12.83) outer-layer dependence U = U(ρ, τ₀, δ, y) | 12.9 | B | NOTE | C11 | the starting line of D19; `defect_law_groups` | #155 |
| N144 | Eq. (12.84) velocity defect law (U_∞ − U)/u_* = F(y/δ) | 12.9 | B | NOTE | C11 | derived in D19; `WT.velocity_defect` | #156 |
| N145 | Eq. (12.85) velocity gradient from the inner law | 12.9 | B | NOTE | C11 | a step of D19 | #157 |
| N146 | Eq. (12.86) velocity gradient from the outer law | 12.9 | B | NOTE | C11 | a step of D19 | #158 |
| N147 | Eq. (12.87) matching the gradients in the overlap: −ξ dF/dξ = y⁺ df/dy⁺ = 1/κ | 12.9 | B | NOTE | C11 | derived in D19 (the separation argument spelled out); `log_law_indicator` | #159 |
| C11 | Eq. (12.88) logarithmic law U⁺ = (1/κ) ln y⁺ + B — by matching the law of the wall to the defect law (12.84) in the overlap (12.85)–(12.87); rough-wall form U⁺ = (1/κ) ln(y/y₀) (12.93) | 12.9 | A | CORE | – | load-bearing: the most used result of wall turbulence — wall functions, skin-friction laws, the roughness length and the neutral wind profile of the atmospheric surface layer are all this logarithm; Ch. 13's bulk drag law starts here | #160 |
| N148 | Eq. (12.89) logarithmic law in outer (defect) form | 12.9 | B | NOTE | C11 | derived in D19; `log_law_defect`, `friction_law_from_overlap` (ours: U_∞⁺ = (1/κ) ln δ⁺ + A + B) | #161 |
| N149 | Fig. 12.18 measured boundary-layer profiles in wall units with the layer names: viscous sublayer, buffer layer, logarithmic layer, wake region | 12.9 | B | NOTE | C11 | stated with `composite_profile` at three δ⁺ and public channel DNS points (Lee & Moser 2015); slider figure F5 and E7 | #162 |
| N150 | Rel: zero-pressure-gradient boundary-layer fits for θ, δ*, δ_99 and C_f | 12.9 | B | NOTE | C11 | stated with `WT.zpg_boundary_layer`; ours: C_f = 2/(U_∞⁺)² | #163 |
| N151 | Rel: two other skin-friction correlations (Schultz-Grunow; White) | 12.9 | B | NOTE | C11 | stated with `skin_friction_zpg(law=)`; log₁₀ vs ln flagged | #164 |
| N152 | Ex. 12.3 boundary-layer thicknesses on a wing at landing speed | 12.9 | B | NOTE | C11 | stated with our own chord and speed, the laminar (Blasius) values alongside; the book's inputs and answers private | #165 |
| N153 | Rel: Spalding's single implicit formula for the whole inner layer | 12.9 | B | NOTE | C11 | stated with `spalding_yplus`, `spalding_uplus` (`brentq`, P108) | #166 |
| N154 | Rel: Coles' outer profile — log law plus a wake term W(y/δ) | 12.9 | B | NOTE | C11 | stated with `coles_wake`, `composite_profile`; the wake strength is passed explicitly (the book's value private) | #167 |
| N155 | Open questions (power-law overlap, stress-gradient layers); the sublayer is universal, the wake is not | 12.9 | C | NOTE | C11 | named; pointer: constants by flow type (N156) | #168 |
| N156 | Value: log-law constants (κ, B) by flow type (channel, pipe, boundary layer) | 12.9 | B | NOTE | C11 | stated with `WT.LOG_LAW_CONSTANTS` (cited public values only; the book's pairs private); κ and B are required keywords in every function | #169 |
| R03 | Eq. (12.90) channel: dP/dx = −2τ₀/h | 12.9 | B | RECAP | C10 | ch08 D07 (force balance on a slug); new for the turbulent mean: the last step of D17; `WT.channel_pressure_gradient`; slip #7 (the proof is Exercise 12.32, not 12.31) | #170 |
| R04 | Eq. (12.91) pipe: dP/dx = −4τ₀/d | 12.9 | B | RECAP | C10 | ch08 (8.8) `core.laminar.pipe_wall_stress` (parity test); `WT.pipe_pressure_gradient` | #171 |
| N157 | Eq. (12.92) empirical link between κ and B across flows and pressure gradients | 12.9 | B | NOTE | C11 | stated with `nagib_chauhan_kappa`, `nagib_chauhan_B` | #172 |
| N158 | Fig. 12.19 smooth vs rough wall; Def: hydrodynamically smooth / rough | 12.9 | B | NOTE | C11 | stated with `log_law` and `rough_wall_log_law` on one plot | #173 |
| N159 | Eq. (12.93) rough-wall logarithmic law; y₀ = roughness length | 12.9 | B | NOTE | C11 | derived in D20; `rough_wall_log_law`, `friction_velocity_from_wind`, `drag_coefficient_neutral`; climate hook: z₀ and the neutral bulk drag coefficient | #174 |
| N160 | Rel (Ex. 12.34): smooth-pipe friction law from the log law; Prandtl's law; laminar 64/Re | 12.9 | B | NOTE | C11 | stated (§4c a-D34) with `pipe_bulk_velocity_loglaw`, `pipe_friction_factor_turbulent` against ch08's `pipe_friction_factor` | #175 |
| N161 | Rel (Ex. 12.29): laminar vs turbulent flat-plate skin friction over six decades of Re_x | 12.9 | B | NOTE | C11 | stated as one log–log figure (`blasius_skin_friction` vs `skin_friction_zpg`) | #176 |
| N162 | Def: eddy viscosity ν_T and eddy diffusivities κ_T, κ_mT — properties of the flow, not of the fluid; turbulent-viscosity and gradient-diffusion hypotheses | 12.10 | B | NOTE | C12 | stated as the opening of C12, with the molecular analogues of ch01 and ch04 recalled | #177 |
| C12 | Eq. (12.94) turbulent-viscosity hypothesis ⟨u_iu_j⟩ = (2/3)ē δ_ij − ν_T(∂U_i/∂x_j + ∂U_j/∂x_i) — with ν_T ~ l_T u_T (12.98) and the mixing-length model l_T = κy that returns the logarithmic law (12.99)–(12.101) | 12.10 | A | CORE | – | load-bearing: the closure used by every ocean and atmosphere model (Ch. 13's Ekman layers need ν_T); the mixing length is the smallest model that reproduces the log law and shows what a closure constant is | #178 |
| N163 | Eq. (12.95) gradient diffusion of heat | 12.10 | B | NOTE | C12 | stated with `gradient_diffusion_flux` | #179 |
| N164 | Eq. (12.96) gradient diffusion of a passive scalar | 12.10 | B | NOTE | C12 | stated; `gradient_diffusion_flux` | #180 |
| N165 | Eq. (12.97) RANS momentum with an effective viscosity ν + ν_T | 12.10 | B | NOTE | C12 | stated (§4c a-D35) with `rans_eddy_viscosity_residual`; ⚠️ slip #6: the pressure gradient carries the free index i (∂P/∂x_i) | #181 |
| N166 | Why the molecular analogy is imperfect: eddies are as large as the gradient length | 12.10 | B | NOTE | C12 | stated with ch01's Knudsen number as the contrast (here l_T/L = O(1)) | #182 |
| N167 | Eq. (12.98) eddy coefficients scale as a turbulent length times a turbulent velocity | 12.10 | B | NOTE | C12 | stated with `eddy_diffusivity_estimate`; closed again by C16's D_T = ⟨u²⟩Λ_t | #183 |
| N168 | Eq. (12.99) unidirectional mean shear flow with an eddy viscosity | 12.10 | B | NOTE | C12 | the starting line of D21; `shear_flow_eddy_viscosity_solve` (Picard-iteration primer) | #184 |
| N169 | Fig. 12.20 and Rel: the mixing-length estimate −⟨uv⟩ = l_T²(dU/dy)², l_T = κy at a wall | 12.10 | B | NOTE | C12 | derived in D21 (the scaling chain, with the displaced parcel of N50); `mixing_length_stress` (sign-preserving), `mixing_length_eddy_viscosity` | #185 |
| N170 | Eq. (12.100) mixing-length model near a wall | 12.10 | B | NOTE | C12 | a step of D21 | #186 |
| N171 | Rel: first integral (constant total stress); ours: its exact solution dU⁺/dy⁺ = 2/(1 + √(1 + 4κ²y⁺²)) | 12.10 | B | NOTE | C12 | derived in D21; `mixing_length_wall_profile` | #187 |
| N172 | Eq. (12.101) outside the sublayer the model returns the logarithmic law | 12.10 | B | NOTE | C12 | derived in D21; the intercept B the model implies with and without van Driest damping is shown, not hidden | #188 |
| N173 | Eq. (12.102) vertical acceleration of a buoyant fluctuation | 12.10 | B | NOTE | C12 | stated (§4c a-D37) with `convective_velocity_scale` | #189 |
| N174 | Rel: free-fall velocity over the layer depth and the resulting eddy diffusivity | 12.10 | B | NOTE | C12 | stated with `convective_eddy_diffusivity` for our own layer (ratio to the molecular value); the book's numbers private | #190 |
| N175 | Rel: one-equation model ingredients u_T = c√ē, ε̄ = C_ε ē^{3/2}/l_T, gradient transport of ē | 12.10 | B | NOTE | C13 | stated with `one_equation_closure`; slip #8 (the viscous transport term carries u_i, as in (12.47)) | #191 |
| C13 | Eq. (12.103) the modelled turbulent-kinetic-energy equation — the k–ε model with ν_T = C_μ ē²/ε̄ (12.104) and the dissipation equation (12.105) | 12.10 | A | CORE | – | key method: the standard two-equation closure of engineering CFD and the template of the turbulent-kinetic-energy closures in ocean and atmosphere models; its constants are tied to two measurable facts (the decay exponent and the von Kármán constant) | #192 |
| N176 | Eq. (12.104) k–ε eddy viscosity ν_T = C_μ ē²/ε̄ | 12.10 | B | NOTE | C13 | derived in D22; `k_epsilon_eddy_viscosity`, `k_epsilon_length_scale` | #193 |
| N177 | Eq. (12.105) modelled dissipation equation | 12.10 | B | NOTE | C13 | stated term by term (a modelled equation, built by analogy, not derived); `k_epsilon_rhs` | #194 |
| N178 | Value: the five standard model constants C_μ, C_ε1, C_ε2, σ_e, σ_ε | 12.10 | B | NOTE | C13 | stated with `K_EPSILON_CONSTANTS` (Launder & Sharma 1974, public); slip #9 (five constants, six printed entries) | #195 |
| N179 | The closed set; wall functions and their limits; sensitivity to inlet values; Reynolds-stress closures | 12.10 | C | NOTE | C13 | named; pointer: the log law used as a boundary condition (C11) and ch10's RANS note | #196 |
| N180 | Rel (ours): two consistency checks — decay of homogeneous turbulence ē ∝ (t + t₀)^{−n}, n = 1/(C_ε2 − 1), and the log layer κ² = √C_μ (C_ε2 − C_ε1) σ_ε | 12.10 | B | NOTE | C13 | derived in D23 (the book never writes it); `k_epsilon_decay`, `k_epsilon_loglayer_kappa` | #197 |
| N181 | Rel (ours): a whole turbulent channel from the mixing length (van Driest damping) plus the linear total stress | 12.10 | B | NOTE | C12 | stated with `channel_mixing_length` against public DNS (labelled approximate); the stage of E8 and the data behind E5 | #198 |
| R05 | Def: set-up for stratified turbulence — stability set by the density gradient in excess of the adiabatic one, so "temperature" means potential temperature | 12.11 | B | RECAP | C14 | ch01 §1.10 (1.29)–(1.32) `core.stratification`, ch07 (7.127), ch11 §11.7; opens C14 with the two-convention table (Kundu Γ = dT/dz, Γ_a ≈ −9.8 K/km; meteorology Γ = −dT/dz, +9.8 K/km) and the ch01 stability badge | #199 |
| N182 | Eq. (12.106) turbulent kinetic-energy budget for a horizontally uniform shear flow | 12.11 | B | NOTE | C14 | derived in D24 (A chain); `stratified_tke_budget`; ⚠️ slip #10: the triple correlation written as ⟨ew⟩ | #200 |
| C14 | Eq. (12.107) flux Richardson number Rf = −gα⟨wT′⟩ / (−⟨uw⟩ dU/dz) = buoyant destruction / shear production — with the gradient Richardson number Ri = N²/(dU/dz)² (12.108) and Ri = (ν_T/κ_T) Rf (12.109) | 12.11 | A | CORE | – | load-bearing for the reader's field: the number that says whether stratification can switch turbulence off; the basis of Ri-dependent mixing schemes in ocean and atmosphere models and of Ch. 13's stratified boundary layers | #201 |
| N183 | Value: turbulence stops being self-supporting near Rf_cr ≈ ¼; large −Rf means convection dominates | 12.11 | B | NOTE | C14 | stated with `turbulence_regime`; ⚠️ two different "one quarter" statements: the theorem Ri > ¼ of ch11 (linear, sufficient for stability) vs the observation Rf_cr ≈ ¼ (Ri = Pr_T/4) | #202 |
| R06 | Eq. (12.108) gradient Richardson number Ri ≡ N²/(dU/dz)² = αg(dT̄/dz)/(dU/dz)² | 12.11 | B | RECAP | C14 | ch11 (11.66)–(11.67) `gradient_richardson`, N² from ch01 and ch07; new here: `gradient_richardson_thermal` takes the in-situ gradient and Γ_a (Kundu sign) and returns the verdict in both conventions | #203 |
| N184 | Eq. (12.109) Ri = (ν_T/κ_T) Rf | 12.11 | B | NOTE | C14 | derived in the last steps of D24; `turbulent_prandtl`, `flux_from_gradient_richardson` | #204 |
| N185 | Def: turbulent Prandtl number ν_T/κ_T — > 1 when stable, small when unstable, ≈ 1 when neutral | 12.11 | B | NOTE | C14 | stated (internal waves carry momentum but not heat, ch07); Reynolds analogy named | #205 |
| C15 | Eq. (12.110) Monin–Obukhov length L_M ≡ −u_*³ / (κ α g ⟨wT′⟩) — with Rf = z/L_M in the logarithmic layer (12.111) and the log-linear wind profile U = (u_*/κ)[ln(z/z₀) + 5 z/L_M] | 12.11 | A | CORE | – | load-bearing for the reader's field: the height at which buoyancy takes over from shear; the scaling length of every surface-flux (bulk) formula in weather and climate models | #206 |
| N186 | Eq. (12.111) in the logarithmic layer Rf = z/L_M | 12.11 | B | NOTE | C15 | derived in D25; `flux_richardson_surface_layer` | #207 |
| N187 | Fig. 12.21 unstable atmosphere: forced convection below ∣L_M∣, free convection above | 12.11 | B | NOTE | C15 | stated with `surface_layer_regime` and a two-zone sketch | #208 |
| N188 | Rel: log-linear wind profile of the stratified surface layer | 12.11 | B | NOTE | C15 | derived in the last steps of D25 (ours: from the dimensionless shear φ_m = 1 + βz/L_M); `WT.surface_layer_wind`; the book's coefficient 5 beside Businger–Dyer's 4.7 and the unstable-side form | #209 |
| N189 | Fig. 12.22 ln z against U: neutral straight, stable more shear, unstable less | 12.11 | B | NOTE | C15 | stated with slider figure F7 and E9 | #210 |
| N190 | Eq. (12.112) budget of temperature variance; ε̄_T = κ⟨(∂T′/∂x_j)²⟩ | 12.11 | B | NOTE | C15 | stated (§4c a-D41) with `temperature_variance_budget`, `temperature_variance_sympy`; slip #15 (factor ½ in the molecular transport) | #211 |
| N191 | Def: dissipation rate of temperature variance; temperature spectrum S_T(K) normalised to the variance | 12.11 | B | NOTE | C15 | stated; `TS.spectrum_variance` reused | #212 |
| N192 | Eq. (12.113) Obukhov–Corrsin spectrum S_T ∝ ε̄_T ε̄^{−1/3} K^{−5/3} | 12.11 | B | NOTE | C15 | stated (§4c a-D42) with `scalar_spectrum`; exponents from `core.dimensional` with temperature as an extra dimension | #213 |
| N193 | Rel: Batchelor scale η_T = η(κ/ν)^{1/2} for ν/κ ≫ 1 | 12.11 | B | NOTE | C15 | stated with `batchelor_scale` (sea water: our number) | #214 |
| N194 | Eq. (12.114) Batchelor's spectrum S_T ∝ K^{−1} in the viscous-convective subrange | 12.11 | B | NOTE | C15 | stated (§4c a-D42); second branch of `scalar_spectrum` | #215 |
| N195 | Fig. 12.23 measured velocity and temperature spectra in a tidal channel (−5/3, then −1 for temperature) | 12.11 | B | NOTE | C15 | stated with `model_spectrum` and `scalar_spectrum` for sea water (our parameters) | #216 |
| N196 | Fig. 12.24 sample particle paths from a point source; Def: Lagrangian position X(a, t); stationary homogeneous turbulence, zero mean velocity | 12.12 | B | NOTE | C16 | stated with `langevin_particles` (Ornstein–Uhlenbeck / Langevin primer) and animation A5; Lagrangian description (ch03) recalled | #217 |
| N197 | Eq. (12.115) rate of change of the mean-square displacement | 12.12 | B | NOTE | C16 | derived in D26 (A chain); `dispersion_rate_from_particles` | #218 |
| N198 | Eq. (12.116) the same in terms of the Lagrangian velocity | 12.12 | B | NOTE | C16 | derived in D26 (both sides evaluated by `dispersion_rate_from_particles`) | #219 |
| N199 | Def: Lagrangian autocorrelation coefficient r_α(τ) | 12.12 | B | NOTE | C16 | stated; `TS.autocorrelation` on particle velocities; ⚠️ α is a no-sum index here, not the expansion coefficient | #220 |
| N200 | Eq. (12.117) dispersion rate = 2⟨u_α²⟩ times the running integral of the correlation | 12.12 | B | NOTE | C16 | derived in D26; ⚠️ trap: r_α(t′ − t) = r_α(t − t′) only because r is even; `taylor_dispersion_rate` | #221 |
| N201 | Eq. (12.118) mean-square displacement as a double integral | 12.12 | B | NOTE | C16 | derived in D26; `taylor_dispersion(form="double")` | #222 |
| N202 | Rel: integration by parts of the double integral | 12.12 | B | NOTE | C16 | a step of D26 (P218a); the two forms compared numerically | #223 |
| C16 | Eq. (12.119) Taylor's formula ⟨X_α²⟩(t) = 2⟨u_α²⟩ t ∫₀ᵗ (1 − τ/t) r_α(τ) dτ — with the ballistic and diffusive limits (12.120)–(12.123) and the eddy diffusivity D_T (12.127)–(12.129) | 12.12 | A | CORE | – | load-bearing: shows from first principles why a cloud spreads as t at first and as √t later, and why an eddy diffusivity is not a constant — the basis of tracer and pollutant dispersion and of the eddy coefficients of C12 and Ch. 13 | #224 |
| N203 | Fig. 12.25 "small t" and "large t" marked on the correlation curve | 12.12 | B | NOTE | C16 | stated: C02's correlation figure with two markers | #225 |
| N204 | Eq. (12.120) short-time (ballistic) limit ⟨X_α²⟩ ≃ ⟨u_α²⟩t² | 12.12 | B | NOTE | C16 | derived in D27 | #226 |
| N205 | Eq. (12.121) rms displacement grows linearly at first (t ≪ Λ_t) | 12.12 | B | NOTE | C16 | derived in D27; `dispersion_regime` | #227 |
| N206 | Eq. (12.122) long-time (diffusive) limit ⟨X_α²⟩ = 2⟨u_α²⟩Λ_t t | 12.12 | B | NOTE | C16 | derived in D27 (with the constant offset −2⟨u²⟩∫₀^∞ τ r dτ the book drops) | #228 |
| N207 | Eq. (12.123) rms displacement grows as √t later (t ≫ Λ_t) | 12.12 | B | NOTE | C16 | derived in D27; slip #11 (the text cites "(11.119)" for (12.119)) | #229 |
| N208 | Rel (ours): closed form for an exponential correlation, ⟨X²⟩ = 2⟨u²⟩Λ_t²[t/Λ_t − 1 + e^{−t/Λ_t}]; Gaussian correlation (Ex. 12.38) | 12.12 | B | NOTE | C16 | derived in the last steps of D27 (exponential case); the Gaussian case stated; `taylor_dispersion_exponential`, `taylor_dispersion_gaussian`; matched by the Langevin particles | #230 |
| N209 | Eq. (12.124) one step of a random walk with uncorrelated directions | 12.12 | B | NOTE | C16 | stated (§4c a-D45) with `random_walk` | #231 |
| N210 | Fig. 12.26 one realisation of a random walk; Rel: the induction ⟨R_n²⟩ = nL² | 12.12 | B | NOTE | C16 | stated with walker paths and the ensemble rms (induction primer) | #232 |
| N211 | Eq. (12.125) rms distance after n uncorrelated steps (R_n)_rms = L√n | 12.12 | B | NOTE | C16 | stated; asserted from `random_walk` within 5 standard errors; n = t/Δt and L = u_rmsΔt with Δt ≈ 2Λ_t tie it to (12.123) | #233 |
| N212 | Fig. 12.27 time-averaged smoke plume in a uniform wind: width ∝ x near the source, ∝ √x far away (t = x/U) | 12.12 | B | NOTE | C16 | stated with `smoke_plume_width`, `plume_concentration`; ⚠️ slip #13: the printed caption swaps the two regimes; the plume mode of E10 | #234 |
| R07 | Rel: the molecular yardstick — a suddenly introduced line vortex spreads as a Gaussian of standard deviation √(2νt) | 12.12 | C | RECAP | C16 | ch08 Exercise 8.26 (spin-up), ch03 / ch05 Lamb–Oseen, ch01 Gaussian spreading; ⚠️ σ² = 2νt here against the core radius σ² = 4νt of ch03 | #235 |
| R08 | Eq. (12.126) a diffusivity from the growth of a variance, ν = ½ dσ²/dt | 12.12 | B | RECAP | C16 | ch01 σ² = 2Dt, `core.diffusion.gaussian_spreading`; `diffusivity_from_variance`; the starting line of D28 | #236 |
| N213 | Eq. (12.127) effective (eddy) diffusivity D_T ≡ ½ d⟨X_α²⟩/dt = ⟨u_α²⟩∫₀ᵗ r_α dτ | 12.12 | B | NOTE | C16 | derived in D28; `eddy_diffusivity_taylor` (exponential case D_T = ⟨u²⟩Λ_t(1 − e^{−t/Λ_t})) | #237 |
| N214 | Eq. (12.128) at short times D_T ≅ ⟨u_α²⟩t | 12.12 | B | NOTE | C16 | derived in D28 | #238 |
| N215 | Eq. (12.129) at long times D_T ≅ ⟨u_α²⟩Λ_t | 12.12 | B | NOTE | C16 | derived in D28; ⚠️ slip #12: the page prints t ≪ Λ_t for both limits — the constant value holds for t ≫ Λ_t | #239 |
| N216 | Turbulent diffusion is not molecular diffusion with a bigger constant: a growing patch is spread by an ever wider range of eddies | 12.12 | B | NOTE | C16 | stated with a constant-D Gaussian ghost of equal late-time width and Richardson's law (N06); ch10's numerical vs eddy diffusivity recalled | #240 |
| N217 | Concluding remarks (research continues; symposia) | 12.13 | C | NOTE | C16 | named in the closing cell; pointer: Ch. 13 (rotating, stratified, geostrophic turbulence) and the DNS / LES literature | #241 |
| S01 | Exercises 12.1–12.38 | Ex. | C | SKIP | – | pointer: one line; never reproduced. Results the lesson uses are B items or D rows: 12.1 (N22), 12.6 (D05), 12.8 (N38), 12.9 (D04), 12.11 (N40), 12.12 (N57), 12.15 (D09), 12.16 (N59), 12.17–12.19 (D07, D08; slip #14), 12.25–12.28 (D16), 12.29 (N161), 12.32 (D17), 12.34 (N160), 12.38 (N208) | #242 |
| S02 | Literature cited and supplemental reading | – | C | SKIP | – | pointer: one bibliography line; public sources cited where used (Sreenivasan 1995; Lee & Moser 2015; Launder & Sharma 1974; Prandtl's friction law via McKeon et al. 2005; AMS Glossary for Businger–Dyer; Pao 1965; Spalding 1961; Coles 1956; Taylor 1921) | #243 |

## 3. Section coverage

| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 12.1 | Introduction | – (covered by B/C items opening the C01 block) | N01, N02 | N03, N04 | – | – |
| 12.2 | Historical Notes | – (covered by a C item opening C01 and a B item closing C16) | N06 | N05 | – | – |
| 12.3 | Nomenclature and Statistics for Turbulent Flow | C01 | N07, N08, N09, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N20, N21, N22, N23 | – | – | – |
| 12.4 | Correlations and Spectra | C02, C03 | N24, N25, N26, N27, N28, N29, N30, N31, N32, N33, N34, N35, N36, N37, N38, N39, N40 | – | – | – |
| 12.5 | Averaged Equations of Motion | C04 | N42, N43, N44, N45, N46, N47, N48, N49, N50, N51, N52, N53, N54, N56, N57, N58, N59, N60 | N41, N55, N61 | R01, R02 | – |
| 12.6 | Homogeneous Isotropic Turbulence | C05 | N62, N63, N64, N65, N66, N67, N68, N69, N70, N71, N72, N73, N74, N75, N76, N77 | – | – | – |
| 12.7 | Turbulent Energy Cascade and Spectrum | C06, C07, C08 | N78, N79, N80, N81, N82, N83, N84, N85, N86, N87, N88, N89, N90, N91, N92, N93, N94, N95, N96, N97 | N98, N99 | – | – |
| 12.8 | Free Turbulent Shear Flows | C09 | N100, N101, N102, N103, N104, N105, N106, N107, N108, N109, N110, N111, N112, N113, N114, N115, N116, N117, N118, N119, N120, N121, N122, N123, N124, N125, N126, N128, N130 | N127, N129 | – | – |
| 12.9 | Wall-Bounded Turbulent Shear Flows | C10, C11 | N131, N132, N133, N134, N135, N136, N137, N138, N139, N140, N141, N142, N143, N144, N145, N146, N147, N148, N149, N150, N151, N152, N153, N154, N156, N157, N158, N159, N160, N161 | N155 | R03, R04 | – |
| 12.10 | Turbulence Modeling | C12, C13 | N162, N163, N164, N165, N166, N167, N168, N169, N170, N171, N172, N173, N174, N175, N176, N177, N178, N180, N181 | N179 | – | – |
| 12.11 | Turbulence in a Stratified Medium | C14, C15 | N182, N183, N184, N185, N186, N187, N188, N189, N190, N191, N192, N193, N194, N195 | – | R05, R06 | – |
| 12.12 | Taylor's Theory of Turbulent Dispersion | C16 | N196, N197, N198, N199, N200, N201, N202, N203, N204, N205, N206, N207, N208, N209, N210, N211, N212, N213, N214, N215, N216 | – | R07, R08 | – |
| 12.13 | Concluding Remarks | – (covered by a C item closing C16) | – | N217 | – | S01, S02 |

Notes. Rows are listed under their own book section (the § column of the chapter map), wherever their A parent sits: N06 (§12.2)
is taught at the end of C16, N181 (§12.10, the mixing-length channel) inside C12. The exercise row S01 and the bibliography row
S02 are listed with §12.13. Every section has at least one item; §12.1, §12.2 and §12.13 have no A item and are covered by B/C
items inside the C01 and C16 blocks.

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Numbering continues at **P280** (the designer assigns numbers). Already primed or taught and only *reminded* (one sentence, no new
primer): P10 seeded generator, P11 Gaussian components, P12 Poisson scatter 1/√N, P13 log–log slopes and `np.polyfit`, P15
`assert np.allclose`, P16 animate, P17 slider_figure, P18 show_viz, P22 `np.gradient`, P25 partial derivative, P26/P98 Taylor
series, P31/P94 `solve_ivp`, P37/P203 trapezoid and Simpson, P38 product rule, P40/P117 sympy, P45 Euler's formula, P48
inequalities under a sign change (the lapse-rate move), P49 chain rule, P60 `fractions.Fraction`, P62 `np.einsum`, P83 iterated
integrals, P84 fundamental theorem of calculus, P87 `quad`, P106 substitution, P108 `brentq`, P109 differentiation under the
integral, P119/P120 isotropic tensors and the deviatoric part, P130 order-of-magnitude scaling, P133 scaled variables, P142
Fourier modes and the FFT, P144 improper integrals, P159 quadratic formula, P167 separation argument, P168 sinh/cosh, P171
sum-to-product, P172 Fourier integral, P181 `np.fft.rfft`, P188 two-length scaling, P189 Leibniz rule, P190
`cumulative_trapezoid`, P192 Picard iteration, P194 Gaussian integral, P195 erf, P197 exponent matching, P206 chain rule with a
moving similarity variable, P218a integration by parts, P252 caching, P255 necessary vs sufficient, P261 even and odd functions;
ch01 C69 Buckingham Π (`core.dimensional.pi_groups`, `solve_exponents`), ch01 potential temperature and the two lapse-rate
conventions (E4 badge), ch02 summation convention, δ_ij, trace, symmetric–antisymmetric contraction and principal axes, ch03
Lagrangian description, ch04 flux form and the Boussinesq set, ch09 thin-layer scaling and thicknesses, ch10 `FD.stretched_grid`,
ch11 (11.88) energy budget.

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| Probability vocabulary: random variable, probability density, histogram (`ax.hist`, `np.histogram`), Gaussian reference values | C01 (N07, N22), C05 (N63) | new; the book assumes it; P11 and P12 used a Gaussian and a Poisson law without the vocabulary |
| Standard error of a sample mean σ/√N and the "within 5 standard errors" test of a Monte-Carlo estimate | C01 (N20), C04 (N50), C16 (N211) | extends P12 from counts to means; every statistical assertion of the chapter uses it |
| Ergodicity: for a stationary process a long time average equals the ensemble average | C01 (N11), C02, C03 | the book uses it silently (theory in ensemble averages, data in time averages) |
| Ornstein–Uhlenbeck process and the Langevin equation; the exact update u_{n+1} = u_n e^{−Δt/τ_c} + σ(1 − e^{−2Δt/τ_c})^{1/2} ξ_n | C01 (`make_ensemble`), C02, C16 (`langevin_particles`) | new modelling tool: a signal with a chosen variance and memory, the test bed for every estimator |
| Sliding-window averages with `np.cumsum` (and `np.convolve`); what happens at the record's ends | C01 (N08, N23) | new Python idiom |
| The sinc function and `np.sinc` (numpy's is sin(πx)/(πx)); zeros at whole periods | C01 (N23) | new; the window factor of Example 12.1 |
| Covariance, the covariance matrix and its ellipse (bivariate Gaussian; sampling a pair with a chosen correlation) | C02 (N25, N27), C04 (N49), C05 (N63) | new; the Reynolds stress is a covariance matrix; ch02's principal axes are recalled |
| A quadratic that is never negative has discriminant ≤ 0 | C02 (D02) | the one move of the Schwartz proof; P159 gives the formula, not this use |
| Correlation by FFT with zero padding against the direct sum (`np.correlate`); biased vs unbiased estimates; noisy tails | C02 (`TS.autocorrelation`) | new numerical tool |
| Fourier-transform pair: where the 2π sits, angular vs cyclic frequency, real + even ⇒ a cosine transform | C03 (D05), C05 (N76) | extends P172 (a packet's Fourier integral) to a transform pair with a fixed normalisation |
| The discrete Fourier transform as a Riemann sum: scaling `np.fft` output to a spectral density, Parseval, one- vs two-sided | C03 (N38), C08 (N94) | extends P181; the commonest source of factors 2 and 2π |
| Segment averaging (Welch), windows and leakage; `scipy.signal.welch` as a cross-check | C03 (N38) | new; why one raw periodogram is as noisy as its mean |
| Change of variable in a density, S(k)dk = S(ω)dω, and the units of a spectral density (per unit ω or k) | C03 (N40), C05 (N76), C08 (D12) | new; needed before S₁₁ can enter a dimensional analysis |
| A sympy averaging operator (linear; kills a single fluctuation; keeps products of fluctuations) | C04 (`rans_sympy`), C06 (`tke_budget_sympy`) | new Python tool; makes the Reynolds rules executable |
| Counting independent components of a symmetric tensor (6 of 9; 10 triple correlations) | C04 (N49, N60) | new; the arithmetic of the closure problem |
| Derivatives of functions of r = ∣r∣: ∂r/∂r_j = r_j/r, ∂r_i/∂r_j = δ_ij, δ_jj = 3 | C05 (D07, D08) | new index-calculus moves; ch02 gives δ_ij but not these |
| Random-phase synthetic fields with a prescribed spectrum (`np.fft.ifftn`; solenoidal by a stream function or a projection) — kinematic, not dynamic | C01 (N02), C05 (N64, N66) | new Python tool; extends P142 |
| Quadrature on a logarithmic grid (K = e^s, dK = K ds) | C08 (N97), C05 | new; spectra span decades |
| Semi-log axes: a logarithm is a straight line on `semilogx`; slope per decade = 2.303/κ | C10, C11 (N149), C15 (N189) | P13 covers log–log only |
| Inner, outer and overlap limits: two descriptions that must agree where both hold (matching) | C10 (N139), C11 (D19) | new method; P167's "function of one variable = function of another ⇒ constant" is recalled inside it |
| Composite profile: inner law + outer correction, valid across both layers | C11 (N153, N154) | new |
| A nonlinear diffusion problem by Picard iteration on the eddy viscosity (freeze ν_T, solve the linear problem, update) | C12 (N168, N181) | extends P192 from time stepping to a steady boundary-value problem |
| Dividing two ODEs to eliminate time (dε̄/dē) and recognising a power-law solution | C13 (D23) | new move |
| Kinematic vs dynamic fluxes: H = ρc_p⟨wT′⟩ [W/m²], τ₀ = ρu_*² [Pa]; α = 1/T for a perfect gas | C04 (N54), C14, C15 | new unit bridge between the book's symbols and observed fluxes |
| Stability parameter ζ = z/L and integrating a flux–profile relation (κz/u_*) dU/dz = φ_m(ζ) | C15 (D25, N188) | new; the form every surface-layer scheme uses |
| A double integral over a triangle: iterated form, and trading it for a single weighted integral | C16 (D26) | extends P83 and P218a |
| Proof by induction | C16 (N210) | new logic tool |

Glosses (one sentence where used, no primer): Wiener–Khinchin theorem (the name of (12.20)); dissipation anomaly; local isotropy;
virtual origin; van Driest damping; wall function; Reynolds analogy; Businger–Dyer functions; virtual temperature; Batchelor
and Obukhov–Corrsin ranges; mole vs mass fraction and stoichiometry (inside N126); log₁₀ vs ln (inside N151, N160).

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Analysis §2b numbers appear as (a-DNN) in the Result column. Every D row stands for one `nb.derivation` block. "book never writes
it" marks rule (b): the result belongs to a B item of the named A block and the book states it without derivation.
| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | commutation rules (12.4)–(12.9), ⟨Ū⟩ = Ū, and ⟨ũṽ⟩ = Ū V̄ + ⟨uv⟩ ≠ Ū V̄ from the definition (12.1) (a-D1) | C01 | ★ | 8 | finite sums and linear operators, limit N → ∞ taken last, expanding a product of two sums | the rules hold because a finite sum commutes with any linear operation — a product is not linear; (12.5) with m = 0 gives ⟨A⟩ = A; commuting with ∂/∂t is exact for ensemble averages only (a sliding time window commutes only approximately); the same N for both variables | notebook · `reynolds_averaging_window` |
| D02 | Schwartz inequality (12.16) and −1 ≤ r ≤ 1 (a-D3; book never writes it) | C02 | ★★ | 6 | ⟨(u + λv)²⟩ ≥ 0 for every real λ, discriminant of a quadratic (primer) | the quadratic in λ has at most one real root, so b² − 4ac ≤ 0 — not ≥; the bound is on ∣⟨uv⟩∣ (the page prints no absolute value); equality only when v ∝ u | notebook |
| D03 | stationary lag forms (12.17), R_11(τ) = R_11(−τ), and R_ij(−τ) = R_ji(τ) for the cross-correlation (a-D4) | C02 | ★ | 5 | shift of the time origin (substitution, P106), stationarity | the shift is allowed only because the statistics do not depend on the origin; the cross-correlation is not even — its mirror swaps the indices | notebook · `reynolds_averaging_window` · `correlation_and_spectrum` |
| D04 | Taylor microscale (12.19) as the intercept of the osculating parabola; ours: λ_t² = 2⟨u²⟩/⟨(du/dt)²⟩ (a-D5; Ex. 12.9, book never writes it) | C02 | ★★ | 7 | Taylor series to second order (P26), even functions (P261), stationarity (d/dt of ⟨u²⟩ = 0) | no linear term because r is even and smooth; r″(0) < 0, hence the minus sign in the definition; an Ornstein–Uhlenbeck signal has a cusp (r ≈ 1 − ∣τ∣/τ_c), so its microscale does not exist | notebook · `correlation_and_spectrum` |
| D05 | the Fourier pair (12.20)–(12.21), S_e real and even, the variance integral (12.22) and S_e(0) = ⟨u_1²⟩Λ_t/π (a-D6) | C03 | ★★ | 9 | complex exponential and Euler's formula (P45), Fourier-transform pair (primer), even functions (P261), improper integrals (P144) | the 1/2π sits in the forward transform and ω is angular; the sine part vanishes because R is even, so S is real; halving the range doubles the integral for S_e(0); two-sided vs one-sided is a factor 2 | notebook · `correlation_and_spectrum` |
| D06 | mean continuity (12.27), ∂u_i/∂x_i = 0 (12.28) and the Reynolds-averaged momentum equation (12.30) from the Boussinesq set via (12.29) (a-D8 + a-D9) | C04 | ★★ | 13 | Reynolds decomposition (12.24)–(12.26), the rules of D01, product rule (P38), summation convention (ch02), ch04's (4.40) for the viscous term | start from the flux form (add ũ_i ∂ũ_j/∂x_j = 0); ⟨U_iu_j⟩ = U_i⟨u_j⟩ = 0 but ⟨u_iu_j⟩ ≠ 0; return ∂(U_iU_j)/∂x_j to U_j∂U_i/∂x_j with (12.27); the new term moves to the right with a minus sign and is read as a stress; T̄ is potential temperature | notebook · `reynolds_stress_parcels` |
| D07 | isotropic two-point tensor (12.40) with F, G; incompressible form (12.41); g = f + (r/2)f′, Λ_g = Λ_f/2, λ_g = λ_f/√2 (a-D13 + a-D14; Ex. 12.17–12.18, book never writes it) | C05 | ★★★ | 14 | isotropic tensors (P119), derivatives of functions of ∣r∣ (primer), ∂R_ij/∂r_j = 0 from (12.28), integration by parts (P218a), Taylor series (P26), sympy (P40) | only δ_ij and r_ir_j are available (no ε_ijk r_k term by reflection symmetry); take r along e₁ to read F and G; δ_jj = 3 and r_jr_j = r² in the divergence; g = (1/2r) d(r²f)/dr makes the Λ integral one line (boundary term needs f to decay faster than 1/r); results are three-dimensional | notebook |
| D08 | isotropic dissipation (12.42) → (12.43): the three gradient moments 2 : 4 : −1 and ε̄ = 30ν⟨u²⟩/λ_f² = 15ν⟨u²⟩/λ_g² = 15ν⟨(∂u_1/∂x_1)²⟩ (a-D15; Ex. 12.19, book never writes it) | C05 | ★★★ | 13 | expanding a squared sum over two indices, symmetry counting with (12.36)–(12.37), ⟨∂_ku_i ∂_lu_j⟩ = −∂²R_ij/∂r_k∂r_l at r = 0 (homogeneity), Taylor expansion f ≈ 1 − r²/λ_f² (P26), sympy (P40) | the square gives 9 + 9 terms that group into three different moments, not equal ones; the minus sign in the gradient-moment identity comes from differentiating at the two points; λ_f vs λ_g is a factor 2 in ε̄; ⟨u²⟩ is one component | notebook |
| D09 | kinetic-energy budget of the mean flow (12.46) from U_i × (12.30) (a-D16; Ex. 12.15, book never writes it) | C06 | ★★ | 9 | product rule to make divergences (P38), mean continuity (12.27), symmetric contraction (ch02 §2.10) | U_i∂τ̄_ij/∂x_j = ∂(U_iτ̄_ij)/∂x_j − τ̄_ij∂U_i/∂x_j — the second piece is the exchange and the dissipation; +⟨u_iu_j⟩∂U_i/∂x_j is a loss for the mean flow because ⟨uv⟩ dU/dy < 0; the gravity term uses ρ̄ = ρ₀[1 − α(T̄ − T₀)] | notebook · `turbulent_energy_budget` |
| D10 | turbulent kinetic-energy budget (12.47): fluctuation equation = total − mean, multiply by u_i, average (the i = j half-trace of (12.35)) (a-D17; book never writes it) | C06 | ★★★ | 15 | the rules of D01, product rule on triple products (P38), ∂u_i/∂x_i = 0 (12.28), index notation (ch02), sympy averaging operator (primer) | the fluctuation equation keeps +∂⟨u_iu_j⟩/∂x_j, which averages away after multiplying by u_i; u_iu_j∂U_i/∂x_j gives the production with a minus sign; ⟨u_i∂p/∂x_i⟩ = ∂⟨pu_i⟩/∂x_i needs continuity; νu_i∇²u_i splits into a transport divergence minus a positive dissipation (the book's 2ν⟨S′S′⟩ form differs from ν⟨(∂_ju_i)²⟩ by a divergence); slip #2 (the page labels the left side Ē) and slip #10 (use ½⟨u_i²u_j⟩) | notebook · `turbulent_energy_budget` |
| D11 | ε̄ ~ (ΔU)³/L (12.48)–(12.49), Kolmogorov scales (12.50), ηu_K/ν = 1 and η/L ~ Re_L^{−3/4} (12.51) (a-D19) | C07 | ★ | 9 | steady state of the budget (C06), order-of-magnitude scaling (P130), exponent solve (ch01 Π theorem, `solve_exponents`) | ε̄ is per unit mass (m²/s³); two equations for the two exponents of ν^a ε̄^b; ν drops out of ε̄ but not of η; ours: u_K/ΔU = Re^{−1/4}, τ_η ΔU/L = Re^{−1/2} | notebook · `energy_cascade_spectrum` |
| D12 | universal form (12.53), the inertial-range law (12.54) with exponent −5/3, and the two-sided normalisation (12.55) (a-D21; slip #1 corrected) | C08 | ★★ | 9 | units of a spectral density (primer), Π theorem with (S_11, ε̄, ν, k_1) (ch01), log–log slope (P13) | S_11 has units m³/s² (velocity² per wavenumber); with ν dropped there is one group, so S_11 ε̄^{−2/3} k_1^{5/3} = const — the exponent of k_1 in S_11 is −5/3 (the page prints +5/3); the two forms of (12.53) are the same by (12.50); the quoted constant is for the two-sided spectrum | notebook · `energy_cascade_spectrum` |
| D13 | thin-layer jet equations (12.61) and the momentum-flux invariant (12.62) from (12.58)–(12.60) (a-D22) | C09 | ★★ | 10 | two-length scaling (P188, ch09), continuity × U added to momentum, integration across the layer, boundary values at ±∞ | the y-equation gives P + ρ⟨v²⟩ = P_∞, which is why ∂P/∂x is dropped; the boundary term [VU + ⟨uv⟩] vanishes because U and ⟨uv⟩ → 0 although V does not; ρ replaces ρ_s downstream (dilution) | notebook · `turbulent_jet_similarity` |
| D14 | the similarity equation (12.63) from (12.61) with the self-preserving forms (12.56)–(12.57) (a-D23; the book writes "somewhat tedious") | C09 | ★★★ | 14 | chain rule with ξ = y/δ(x) (P206), Leibniz rule (P189), integration by parts in ξ (P218a), sympy (P40) | ∂U/∂x = U′_CL F − U_CL F′ξδ′/δ (the similarity variable moves with x); V = −∫₀^y ∂U/∂x dy needs ∫ξF′dξ = ξF − ∫F dξ; the two ξFF′ terms cancel; divide by U_CL²/δ to expose the three brackets | notebook · `turbulent_jet_similarity` |
| D15 | δ ∝ x, U_CL ∝ x^{−1/2} from (12.64)–(12.65), the far field (12.66)–(12.67) and the volume flux ∝ x^{1/2} (12.68) (a-D24) | C09 | ★★ | 10 | "a function of x equal to a function of ξ is a constant" (P167), separable ODE, exponent matching (P197), Gaussian integral (P194) | with δ = x, δU′_CL/U_CL = C_1 integrates to a power law x^{C_1}; the invariant demands 2γ + 1 = 0; C_5 = C_4(ρ/J_s)^{1/2} (slip #4); ours: the entrainment velocity is ½ dV̇/dx on each side and the local Reynolds number grows as x^{1/2} | notebook · `turbulent_jet_similarity` |
| D16 | exponents of Table 12.1 for the round jet, plane and round wake, plane and round plume and shear layer from each flow's invariant and growth law (a-D27; Ex. 12.25–12.28, book never writes it) | C09 | ★★ | 12 | D13–D15 as the template, exponent matching (P197), exact rationals (`fractions.Fraction`, P60) | the invariant changes with the flow (momentum flux; momentum deficit = drag; buoyancy flux); wakes linearise advection to U_∞∂/∂x, so dδ/dx ~ ΔU/U_∞ instead of a constant; round flows integrate over 2πr dr; the round jet's local Reynolds number is constant, the round wake's falls | notebook · `turbulent_jet_similarity` |
| D17 | linear total stress τ̄ = τ₀(1 − 2y/h) from (12.76)–(12.77) and the pressure-gradient balances (12.90)–(12.91) (a-D29) | C10 | ★★ | 9 | fully developed flow (∂/∂x of statistics = 0, V = 0), integration of a constant, symmetry about the centreline, control-volume balance (ch08 D07) | τ̄ depends on y only and dP/dx on x only, so each is a constant; τ̄(h) = −τ̄(0) by symmetry gives (12.90) with h the full height; the wall-normal equation moves ρ⟨v²⟩ into the pressure but not into dP/dx; slip #7 | notebook · `law_of_the_wall` |
| D18 | law of the wall (12.80) with u_* (12.81) by the Π theorem, and the viscous sublayer U⁺ = y⁺ (12.82) (a-D30) | C10 | ★ | 7 | Π theorem (ch01): 5 variables, 3 dimensions ⇒ 2 groups, repeating variables (ρ, τ₀, ν); integration of a constant | δ and U_∞ are deliberately left out of the inner list; the stress is constant across the sublayer only because the layer is thin (D17); U⁺ = y⁺ is 5 % off by y⁺ ≈ 5–7 | notebook · `law_of_the_wall` |
| D19 | defect law (12.84), overlap matching (12.85)–(12.87), the logarithmic law in inner and outer form (12.88)–(12.89); ours: the friction law U_∞⁺ = (1/κ) ln δ⁺ + A + B (a-D31 + a-D32) | C11 | ★★ | 12 | Π theorem (ch01), chain rule (P49), matching in the overlap (primer), separation argument (P167), ∫dy/y = ln y | the defect is used because the outer flow knows U only up to a shift; multiply both gradients by y/u_* before equating; the constant is 1/κ on both sides with opposite signs for the defect; ln, not log₁₀; the overlap exists only when δ⁺ ≫ 1 | notebook · `law_of_the_wall` |
| D20 | rough-wall law (12.93) from (12.88) with U = 0 at y = y₀; ours: equivalent B = −(1/κ) ln y₀⁺ and the neutral drag coefficient C_D = [κ/ln(z/z₀)]² (a-D33) | C11 | ★ | 5 | logarithm rules, definition of a drag coefficient | y₀ is not the size of the roughness elements; viscosity has left the problem, so the friction coefficient no longer depends on Re; y ≤ y₀ is outside the law | notebook · `law_of_the_wall` |
| D21 | mixing length: −⟨uv⟩ = l_T²(dU/dy)², the wall model (12.100), its first integral, the exact slope dU⁺/dy⁺ = 2/(1 + √(1 + 4κ²y⁺²)) and the logarithmic limit (12.101) (a-D36) | C12 | ★★ | 11 | (12.94) for simple shear, the parcel estimate of C04, quadratic formula (P159), limits y⁺ → 0 and y⁺ → ∞, quadrature (P190) | keep the sign: l_T²∣dU/dy∣dU/dy; the integration constant is τ₀/ρ = u_*² evaluated at the wall; take the positive root; the plain l_T = κy gives the right slope 1/κ but far too small an intercept B — wall damping is what fixes it | notebook · `mixing_length_closure` |
| D22 | modelled energy equation (12.103) and the k–ε eddy viscosity (12.104) from (12.47) and (12.98) (a-D38) | C13 | ★ | 6 | gradient-diffusion model for the transport divergence, substitution, units check | the three transport terms are lumped into one gradient flux (ν_T/σ_e)∂ē/∂x_j; l_T = ē^{3/2}/ε̄ and u_T = √ē give ν_T ∝ ē²/ε̄; the production keeps its exact form | notebook |
| D23 | consistency of the k–ε constants: decay ē ∝ (t + t₀)^{−n} with n = 1/(C_ε2 − 1), and the log layer κ² = √C_μ (C_ε2 − C_ε1) σ_ε, ē = u_*²/√C_μ (a-D39; ours, book never writes it) | C13 | ★★ | 12 | dividing two ODEs (primer), separable ODE and power laws, production = dissipation in the log layer, ε̄ = u_*³/(κy) from C11 | homogeneous decay removes transport and production; integrate dε̄/dē = C_ε2 ε̄/ē first; in the log layer ē is uniform so only the diffusion of ε̄ survives; the standard constants give κ ≈ 0.43, not exactly 0.41 | notebook |
| D24 | the reduced budget (12.106) from (12.47), the flux Richardson number (12.107) and Ri = (ν_T/κ_T) Rf (12.109) (a-D40, first part) | C14 | ★ | 8 | horizontal homogeneity (∂/∂x, ∂/∂y of statistics = 0), the closures (12.94)–(12.95), N² from the potential-temperature gradient (ch01) | signs: upward heat flux ⟨wT′⟩ > 0 makes Rf < 0 (unstable); −⟨uw⟩ = ν_T dU/dz and ⟨wT′⟩ = −κ_T dT̄/dz; T̄ is potential temperature — with the in-situ gradient subtract Γ_a (Kundu) or compare with Γ_d (meteorology); Rf_cr ≈ ¼ is an observation, Ri > ¼ of ch11 a theorem | notebook · `stratified_surface_layer` |
| D25 | Rf = z/L_M (12.111) from (12.107), (12.110) and the log-layer values, and the log-linear wind profile by integrating φ_m = 1 + βz/L_M (a-D40, second part; profile step ours) | C15 | ★ | 7 | log law of C11 (dU/dz = u_*/κz, −⟨uw⟩ = u_*²), substitution, ∫(1/z + β/L_M)dz | L_M carries a minus sign so that it is positive when stable (downward heat flux); κ here is von Kármán's constant, not the diffusivity; the profile uses z₀ as the lower limit; the linear correction is the stable-side form (unstable side: Businger–Dyer) | notebook · `stratified_surface_layer` |
| D26 | Taylor's formula: (12.115)–(12.117), the double integral (12.118) and its integration by parts to (12.119) (a-D43) | C16 | ★★ | 11 | the rules (12.6)–(12.7) of D01, chain rule (P49), substitution τ = t − t′ (P106), even functions (P261), integration by parts (P218a) | the average passes inside the time integral by (12.7); r(t′ − t) = r(t − t′) only because r is even; the boundary term of the integration by parts gives the t∫r dτ piece; α is a no-sum index | notebook · `taylor_dispersion` |
| D27 | the ballistic limit (12.120)–(12.121), the diffusive limit (12.122)–(12.123) with its constant offset, and the closed form for r = e^{−τ/Λ_t} (a-D44) | C16 | ★ | 8 | limits of an integral (r ≈ 1; τ/t → 0), improper integrals (P144), exponential integrals, sympy (P40) | the long-time form drops −2⟨u²⟩∫₀^∞ τ r dτ, a constant, so ⟨X²⟩ approaches a straight line that does not pass through the origin; Λ_t here is the Lagrangian integral scale; slip #11 | notebook · `taylor_dispersion` |
| D28 | eddy diffusivity: ν = ½ dσ²/dt (12.126) from Gaussian spreading, D_T (12.127) as half the slope of (12.119), and its limits (12.128)–(12.129) (a-D46; slip #12 corrected) | C16 | ★ | 7 | variance of a Gaussian (ch01 σ² = 2Dt), differentiation under the integral (P109), limits | σ² = 2νt per coordinate (ch03's core radius uses 4νt); D_T is not a constant — it grows as ⟨u²⟩t and saturates at ⟨u²⟩Λ_t only for t ≫ Λ_t (the page prints ≪ for both); D_T = u_T l_T with l_T = u_rms Λ_t recovers (12.98) | notebook · `taylor_dispersion` |

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C01** | a-D2 Example 12.1 window factors | ⟨u⟩(t) = [sinh(Δt/2τ)/(Δt/2τ)] A e^{−t/τ} + [sin(ωΔt/2)/(ωΔt/2)] B cos ωt; good window 1 ≪ ωΔt ≪ ωτ | N23 | a B item whose integral the book writes out; the code cell does the window integral with sympy and the slider figure shows both factors |
| **C03** | a-D7 periodogram form of the spectrum (★★★) | S_e(ω) = lim_{T→∞} (1/2πT)∣∫_{−T/2}^{T/2} u e^{−iωt} dt∣² (Ex. 12.8) | N38 | the book never writes it, but the lesson needs the result, not the double-integral proof: the triangular weight (1 − ∣τ∣/T) → 1 is said in one sentence and the equality is shown numerically three ways (Parseval to round-off) |
| **C04** | a-D10 mean temperature and heat-flux form | (12.31), (12.32) | N52, N53 | the D06 moves applied to (4.89), which the book writes; `rans_sympy` prints the steps |
| **C04** | a-D11 mean scalar equation | (12.34) (Ex. 12.12) | N57 | the D06 moves a third time at constant ρ_m; one sentence on where the constant-density assumption enters; checked by `rans_sympy` |
| **C04** | a-D12 Reynolds-stress transport (★★★) | (12.35) (Ex. 12.16) | N59 | the lesson needs its trace, which D10 derives in full from the same fluctuation equation; (12.35) is D10 with the second index kept free (multiply the i-equation by u_j, the j-equation by u_i, add) — stated term by term and verified by `reynolds_stress_budget_sympy` (one component, and trace/2 = (12.47)) |
| **C06** | a-D18 no shear production in isotropic turbulence; mean viscous dissipation ~ 1/Re of production | ⟨u_iu_j⟩∂U_i/∂x_j = ⟨u_1²⟩∂U_i/∂x_i = 0; ratio ~ ν/(UL) = 1/Re | N80, N83 | two one-line arguments the book writes |
| **C07** | a-D20 Taylor-microscale scaling | (12.52) λ_T/L ∝ Re_L^{−1/2}; R_λ ~ Re_L^{1/2}; ours: λ/η ~ Re^{1/4} | N91, N92 | a proportionality chain the book writes; numbers from `scale_separation` instead |
| **C09** | a-D25 scalar field of the plane jet | (12.69)–(12.73): Ȳ = C_6 (Ṁ_s/√(ρJ_s)) x^{−1/2} H(y/x) | N117–N121 | algebra the book writes; the second invariant (slot-fluid flux) is shown flat in x by `scalar_flux_per_span` (slip #3) |
| **C09** | a-D26 general similarity condition | (12.74) with the power family δ ~ x^m, U_CL ~ x^n, Ψ ~ x^{2n+m−1} | N125 | a remark of the book; the sympy check of both families is shown (slip #5: only m + 2n = 0 survives the invariant) |
| **C09** | a-D28 jet energy budget | (12.75) | N128 | a reduction of (12.47) the book writes; the terms are plotted from a labelled model (slip #10) |
| **C11** | a-D34 pipe friction law from the log law | U_av = u_*[(1/κ) ln(au_*/ν) + B − 3/(2κ)]; f^{−1/2} = 2.0 log₁₀(Re_d f^{1/2}) − 0.8 (Ex. 12.34) | N160 | the book never writes it, but it is an application outside the A chain: the area integral of the log law is one integration by parts, named, and checked numerically (`pipe_bulk_velocity_loglaw` vs quadrature; derived constants vs Prandtl's) |
| **C12** | a-D35 effective-viscosity form of RANS | (12.97) with ∂P/∂x_i (slip #6) | N165 | one substitution of (12.94) into (12.30); the −⅔ē δ_ij term is read as an extra pressure |
| **C12** | a-D37 convective velocity and eddy diffusivity | (12.102): w ~ (gLΔT/T)^{1/2}, κ_T ~ wL | N173, N174 | a scaling chain the book writes; our own layer's numbers instead |
| **C15** | a-D41 temperature-variance budget | (12.112), with κ∂(½⟨T′²⟩)/∂z in the molecular transport (slip #15) | N190 | the recipe of D10 with T′ for u_i, said in words by the book; `temperature_variance_sympy` shows the factor ½ |
| **C15** | a-D42 scalar spectra and the Batchelor scale | (12.113) S_T ∝ ε̄_T ε̄^{−1/3} K^{−5/3}; η_T = η(κ/ν)^{1/2}; (12.114) S_T ∝ K^{−1} | N192, N193, N194 | the D12 argument with temperature as an extra dimension; exponents returned by `core.dimensional.solve_exponents` in the cell |
| **C16** | a-D45 random walk | (12.124)–(12.125): ⟨R_n²⟩ = nL², (R_n)_rms = L√n | N209, N210, N211 | the book writes every step; stated with the induction primer and confirmed by simulation (within 5 standard errors) |

## 5. Interactive explainers (5–10 + backup)
Ten explainers; every A item except C05 (isotropic relations — static figures and two notebook ★★★ derivations serve it better) and
C13 (k–ε — the backup) is attached to at least one. Every window shows its equations with numbers, has an Explain tab ("Explanation &
interpretation", numbered sections with the reader's numbers and a "Reading the current setting" paragraph), a synced Code tab, a
Derivation tab for its D rows, a 4–8-step walkthrough and ≥ 3 check questions. **Colour code across the chapter:** mean = purple
(accent), fluctuation and turbulent energy = teal, Reynolds stress / shear production = orange, viscous stress / dissipation = rose,
buoyancy and temperature = blue, scalar = amber, laminar or reference ghosts = muted. Reference explainers (`interactive-viz` §5):
FDV = `forced_damped_vibrations.html`, AFE = `angular_frequency_explorer_1.html`, APS = `amplitude_phase_second_order_II_3.html`.
**JS physics stays light:** seeded Ornstein–Uhlenbeck signals and particles with the exact update (`Viz.rng`), closed forms, one
Newton solve (Spalding) and 1-D quadratures; no table is needed except E5's channel profile (ours, cached, parity-checked at ≥ 2
points). Two sign conventions (E9): verdict word → fluidpy criterion text in the chosen convention → the other convention's bare
relation, pinned by exact-text selftest rows (ch01 E4 pattern).

### E1 · reynolds_averaging_window
- A: C01, C02 (also shows N07, N08, N18, N20, N23 Example 12.1, N32 independent samples, N11) · **Confusion removed:** "what is 'the mean' of a signal that never repeats, and how long must I average?" — the ensemble mean is a sum over realizations; one record can replace it only if the window is long against the memory Λ_t and short against the drift, and a window too short lets the fluctuation leak into the "mean".
- **Why interactive:** dragging the window Δt and the number of members N shows the two estimates of the mean converge and fail in different ways (scatter ∝ N^{−1/2}; the sinc and sinh-ratio factors of Example 12.1), and the product panel shows ⟨ũṽ⟩ − ŪV̄ refusing to vanish — a static Fig. 12.3 shows three values of N, not the trade-off.
- **Stage:** (1) the record: grey members, the ensemble mean of N members (purple), the sliding time average of one member (teal) with the window drawn as a band, the true mean as a ghost; (2) error of each estimate against N and against Δt (log–log, the −½ line; the window factors with the current point); (3) (hidePortrait) the mean of the product split into ŪV̄ and ⟨uv⟩ as two bars.
- **Controls:** members N · window Δt · memory τ_c · signal (chips: stationary, decaying mean, mean + wave) · transport t (new members appear).
- **Equations:** (12.1) $\langle u^m\rangle=\lim_{N\to\infty}\frac1N\sum_n u(\mathbf x,t{:}n)^m$; (12.2) $\overline{u}=\frac1{\Delta t}\int_{t-\Delta t/2}^{t+\Delta t/2}u\,dt$; (12.6) $\overline{\partial u/\partial t}=\partial\bar u/\partial t$; the Example 12.1 factors $\frac{\sinh(\Delta t/2\tau)}{\Delta t/2\tau}$, $\frac{\sin(\omega\Delta t/2)}{\omega\Delta t/2}$; $\overline{\tilde u\tilde v}=\bar u\bar v+\overline{uv}$; (12.17) $R_{11}(\tau)=\overline{u_1(t)u_1(t+\tau)}$.
- **Mirrors** `TS.make_ensemble`, `TS.ensemble_average`, `TS.time_average`, `TS.standard_error`, `ch12.time_average_exp_cos`, `TS.product_average_split`, `TS.effective_samples`.
- derivations: D01, D03 · depth features: explain (the two window factors with your Δt, τ, ω; the standard error with your N and the number of independent samples Δt/Λ_t; the product split; the reading), code, **linked views**, **transport**, **presets** (too short a window, one whole period, good window, N = 2 / 8 / 64), **status** ("window 0.4 periods: 76 % of the wave leaks into the mean" / "good window: 1 ≪ ωΔt ≪ ωτ"), **terms** (ŪV̄ + ⟨uv⟩ = ⟨ũṽ⟩) · follows: FDV (signal + toggled curves on one time slider, Explain panel) · **aha:** a mean is only as good as the number of *independent* samples in it — N members, or Δt/Λ_t memory times.

### E2 · correlation_and_spectrum
- A: C02, C03 (also shows N30 delay peak, N31, N33 Taylor microscale, N37 S_e(0), N38 periodogram, N40 frozen field) · **Confusion removed:** "what do the autocorrelation and the spectrum each tell me, and why are they the same information?" — r(τ) says how long the signal remembers itself; S_e(ω) says which frequencies carry the variance; they are a Fourier pair, so a long memory is a narrow spectrum.
- **Why interactive:** stretching the memory time narrows the spectrum in front of the reader while the area under S stays equal to the variance and S(0) follows Λ_t; sliding the lag shows the product u(t)u(t + τ) being averaged — the reciprocity is a motion, not a pair of static plots.
- **Stage:** (1) the signal and its copy shifted by the lag τ, their product shaded; (2) r(τ) with the current lag, the equal-area rectangle Λ_t, the first zero t_c and the osculating parabola λ_t; (3) S_e(ω) on log axes: exact transform (bold), periodogram of the record (faint), S(0) marker, shaded area = variance.
- **Controls:** correlation shape (modes: exponential, Gaussian, damped cosine) · memory time τ_c · oscillation frequency ω₀ · lag τ · axis (frequency / wavenumber with probe speed U₀).
- **Equations:** (12.17); (12.18) $\Lambda_t=\int_0^\infty r_{11}\,d\tau$; (12.19) $\lambda_t^2=-2/r_{11}''(0)$; (12.20) $S_e(\omega)=\frac1{2\pi}\int R_{11}(\tau)e^{-i\omega\tau}d\tau$; (12.21) $R_{11}(\tau)=\int S_e(\omega)e^{+i\omega\tau}d\omega$; (12.22) $\overline{u_1^2}=\int S_e\,d\omega$; $S_e(0)=\overline{u_1^2}\Lambda_t/\pi$.
- **Mirrors** `TS.correlation_spectrum_pair` (§8), `TS.integral_scale`, `TS.taylor_microscale`, `TS.spectrum_from_correlation`, `TS.spectrum_variance`, `TS.frequency_to_wavenumber_spectrum`.
- derivations: D03, D04, D05 · depth features: explain (Λ_t, t_c, λ_t and S(0) with your settings; the variance by integrating S; the half-width of S times Λ_t; the reading), code, **linked views**, **modes**, **presets** (short memory, long memory, a hidden wave, Gaussian pair), **inspector** (click the spectrum: the cosine-transform arithmetic at that ω), **status** ("memory 0.5 s ⇒ spectrum flat to ω ≈ 2 rad/s") · follows: APS (linked windows with a numbered live derivation) · **aha:** stretch the memory and the spectrum squeezes toward zero frequency — its height there is the integral scale and its area is always the variance.

### E3 · reynolds_stress_parcels
- A: C04 (also shows N42, N48, N49, N50 Fig. 12.6, N51 momentum flux, N63 scatter plots, N60 closure count) · **Confusion removed:** "how can fluctuations that average to zero push on the mean flow, and why is ⟨uv⟩ negative when the shear is positive?" — a parcel carried upward arrives with the lower speed of where it came from (u < 0 with v > 0), one carried downward arrives fast (u > 0 with v < 0): every exchange gives uv < 0, a net downward flux of x-momentum.
- **Why interactive:** the reader releases parcels in a shear and watches each land in the (u, v) scatter; turning the shear down, to zero and negative tilts the cloud, and the correlation coefficient and the stress bar follow — the sign argument of Fig. 12.6 becomes a count.
- **Stage:** (1) mean profile U(y) with parcels hopping between levels, each coloured by the sign of its uv; (2) the (u, v) scatter with its covariance ellipse and principal axes, quadrant counts; (3) (hidePortrait) the stress bars: viscous μ dU/dy, Reynolds −ρ⟨uv⟩, total.
- **Controls:** shear dU/dy (both signs) · displacement l_rms · v_rms · how well a parcel keeps its momentum (correlation 0–1) · transport (parcels accumulate).
- **Equations:** (12.24) $\tilde u_i=U_i+u_i$; (12.30) with $\bar\tau_{ij}=-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}$; $\rho_0\overline{(U+u)v}=\rho_0\overline{uv}$; the parcel estimate $u\approx-\ell\,dU/dy$, $\overline{uv}\approx-\overline{v\ell}\,dU/dy$; (12.14) $r_{12}$.
- **Mirrors** `ch12.displaced_parcel_uv`, `ch12.parcel_uv_expected` (§8), `TS.reynolds_stress`, `TS.correlation_coefficient`, `ch12.mean_stress_tensor`.
- derivations: D06 · depth features: explain (⟨uv⟩ from your shear, l_rms, v_rms with the sampling error; r_uv; the stress in Pa for air and for water; the eddy viscosity −⟨uv⟩/(dU/dy) it implies; the reading), code, **linked views**, **transport**, **presets** (strong shear, no shear — isotropic cloud, reversed shear, parcels that forget), **terms** (viscous + Reynolds = total), **inspector** (click a parcel: its ℓ, u, v and uv), end-of-run card · follows: `forward_noising_lab` (click a sample to see its arithmetic) + FDV (Explain) · **aha:** no single fluctuation pushes the flow — the push is the *correlation* between going up and being slow.

### E4 · energy_cascade_spectrum
- A: C07, C08 (also shows N85–N87 cascade tiers, N89, N91 Taylor microscale, N92, N94–N97, N90 real cases, slip #1) · **Confusion removed:** "if viscosity is what dissipates the energy, why does the dissipation rate not depend on viscosity — and what does a higher Reynolds number change?" — the large eddies fix the supply ΔU³/L; viscosity only decides how small the eddies must get before they can dissipate it, so a higher Re makes the cascade longer.
- **Why interactive:** dragging Re_L slides η away from L on the log axis, adds tiers to the cascade and widens the −5/3 range of the spectrum while its level stays put; changing ν at fixed ΔU and L moves the right end only — the reader sees what depends on what.
- **Stage:** (1) the cascade ladder on a log length axis from L to η, tiers with their velocities and turnover times, λ_T marked; (2) the model spectrum in Kolmogorov scaling or in physical units (toggle), the −5/3 line, inertial range shaded, the printed +5/3 as a ghost; (3) (hidePortrait) a table of real cases with the current row highlighted.
- **Controls:** ΔU · L · ν (chips: air, water) · view (physical / Kolmogorov-scaled) · case presets.
- **Equations:** (12.49) $\bar\varepsilon\sim(\Delta U)^3/L$; (12.50) $\eta=(\nu^3/\bar\varepsilon)^{1/4}$, $u_K=(\nu\bar\varepsilon)^{1/4}$; (12.51) $\eta/L\sim\mathrm{Re}_L^{-3/4}$; (12.52) $\lambda_T/L\propto\mathrm{Re}_L^{-1/2}$; (12.53) $S_{11}/(u_K^2\eta)=\Phi(k_1\eta)$; (12.54) $S_{11}=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$.
- **Mirrors** `ch12.dissipation_outer_scaling`, `ch12.kolmogorov_scales`, `ch12.scale_separation`, `ch12.cascade_tiers`, `ch12.model_spectrum`, `ch12.inertial_spectrum_1d`, `ch12.dns_grid_points`.
- derivations: D11, D12 · depth features: explain (ε̄, η, u_K, τ_η, λ_T with your numbers; decades of inertial range; grid points for a DNS; the exponent solve; the reading), code, **linked views**, **presets** (kitchen mixer, wind tunnel, atmospheric boundary layer, ocean thermocline), **status** ("Re_L = 10⁷: η/L = 5.6 × 10⁻⁶ — 5.25 decades below L"), **inspector** (click the spectrum: the units arithmetic of ε̄^{2/3}k₁^{−5/3}), **notes** (a model grid of 1 km against η) · follows: `overfitting_curves` (a marker and verdict on a two-slider figure) + AFE (table with the current row, presets) · **aha:** viscosity does not set how much energy is dissipated, only how far down the ladder the eddies must go to do it.

### E5 · turbulent_energy_budget
- A: C06 (also shows N78 mean-flow budget, N80, N81, N136 linear stress, N181 mixing-length channel, N149 production peak in the buffer layer; ch11 (11.88) recalled) · **Confusion removed:** "where does the energy of the turbulence come from and where does it go?" — the term −⟨uv⟩ dU/dy is a loss in the mean-flow budget and the same amount is a gain in the turbulent budget; it peaks near the wall, and what is produced is dissipated locally or carried away by transport.
- **Why interactive:** clicking a height in a channel shows the two budgets' bars side by side with the exchange term mirrored; sliding Re_τ moves the production peak (always near y⁺ ≈ 12 in wall units) and shrinks the direct viscous loss of the mean flow — the "same term, opposite sign" is seen, not read.
- **Stage:** (1) channel profiles U⁺(y⁺), −⟨uv⟩⁺ and production with a cursor; (2) paired term bars at the cursor: mean flow (pressure work, viscous dissipation, loss to turbulence) and turbulence (production, dissipation + transport as modelled, labelled); (3) (hidePortrait) the integrated budget across the channel: work by the pressure gradient = direct dissipation + production.
- **Controls:** Re_τ · cursor height y⁺ · wall damping on/off · axis (wall units / outer units).
- **Equations:** (12.46) with $-2\nu\bar S_{ij}\bar S_{ij}+\overline{u_iu_j}\,\partial U_i/\partial x_j$; (12.47) with $-\overline{u_iu_j}\,\partial U_i/\partial x_j-\bar\varepsilon$; for $U(y)$: production $=-\overline{uv}\,dU/dy$; the ratio $2\nu\bar S_{ij}\bar S_{ij}/\text{production}\sim1/\mathrm{Re}$.
- **Mirrors** `ch12.channel_energy_budget` (§8), `ch12.mean_energy_budget`, `ch12.tke_budget`, `ch12.shear_production`, `ch12.channel_mixing_length`.
- derivations: D09, D10 (★★★, all 15 steps) · depth features: explain (each term at your y⁺ with numbers; where production peaks and why (viscous stress = Reynolds stress there); the integrated balance; what is modelled and what is exact; the reading), code, **linked views**, **terms** (two bar groups that share one mirrored term), **inspector** (click a height), **presets** (sublayer, buffer-layer peak, log layer, centreline), **status** ("y⁺ ≈ 12: production at its maximum, ¼ in wall units, where viscous and Reynolds stress are equal") · follows: `fid_formula_lab` (clickable terms and bars) + FDV (Explain) · **aha:** production is one term seen from two sides — the mean flow's loss is the turbulence's income, and it is largest where viscous and turbulent stresses are equal.

### E6 · turbulent_jet_similarity
- A: C09 (also shows N101 entrainment, N109 invariant, N116, N122, N123–N124 Table 12.1 flows, N130 eddy-viscosity sech², N125 slip #5) · **Confusion removed:** "how can anyone predict how a turbulent jet spreads without solving for the turbulence?" — the momentum flux cannot change and the profiles keep their shape; those two facts alone force δ ∝ x and U_CL ∝ x^{−1/2}, and the missing fluid is drawn in from the sides.
- **Why interactive:** toggling raw ↔ rescaled makes the profiles at five stations collapse onto one F(ξ); a wrong-exponent slider breaks either the collapse or the invariant bar; switching the flow (round jet, wake, plume) changes the invariant and the exponents follow — the argument is checked by breaking it.
- **Stage:** (1) the jet in the x–y plane: mean-velocity colour, half-width lines, entrainment arrows, the laminar jet (x^{2/3}, x^{−1/3}) as a ghost; (2) profiles at stations, raw or rescaled; (3) invariant bars along x: momentum flux (flat), volume flux (rising), scalar flux (flat).
- **Controls:** flow (modes: plane jet, round jet, plane wake, plane plume) · view raw / rescaled · trial decay exponent n · station x · scalar on/off.
- **Equations:** (12.56) $U=U_{CL}(x)F(y/\delta)$; (12.62) $J_s=\rho\int U^2dy$; (12.63) the three brackets; (12.64); (12.65) $2\gamma+1=0$; (12.66) $U=C_5(J_s/\rho)^{1/2}x^{-1/2}F(y/x)$; (12.68) $\dot V\propto x^{1/2}$; (12.71) $\bar Y\propto x^{-1/2}$.
- **Mirrors** `ch12.plane_jet_mean_velocity`, `ch12.jet_momentum_flux_per_span`, `ch12.plane_jet_volume_flux`, `ch12.free_shear_exponents`, `ch12.free_shear_flow`, `ch12.general_similarity_check`, `core.jets.free_jet` (ghost).
- derivations: D13, D14 (★★★, all 14 steps), D15, D16 · depth features: explain (the two exponent equations for your flow; J_s, V̇ and U_CL at your station; the entrainment velocity; the local Reynolds number trend; the reading), code, **linked views**, **modes** (flows), **presets** (correct exponents, wrong decay, laminar jet), **status** ("n = −0.4: momentum flux grows as x^{0.2} — not allowed"), **terms** (the three invariant bars) · follows: AFE (modes and linked views) with the ch09 `free_jet_similarity` stage (raw ↔ rescaled, invariant bars) · **aha:** conservation fixes one exponent and shape-preservation the other — no turbulence model is needed to know how a jet spreads.

### E7 · law_of_the_wall
- A: C10, C11 (also shows N136 linear stress, N137, N139 layers, N142 sublayer, N144 defect law, N147 indicator, N149, N153 Spalding, N154 Coles, N159 rough wall, N156) · **Confusion removed:** "why do profiles from different flows and Reynolds numbers fall on one curve near a wall, and where does the logarithm come from?" — near the wall only τ₀ and ν matter (wall units), far from it only τ₀ and δ; where both descriptions hold, y dU/dy can depend on neither, so it is a constant.
- **Why interactive:** sliding Re_τ stretches the log region while the inner part of the curve stays fixed in wall units and the wake peels off at a different y⁺; switching to outer units makes the opposite collapse; the indicator y⁺dU⁺/dy⁺ shows the plateau 1/κ appear only when δ⁺ is large — one static Fig. 12.18 cannot show both collapses.
- **Stage:** (1) U⁺ against y⁺ on a semi-log axis with the sublayer line, the log line and layer bands (sublayer, buffer, log, wake), public DNS points as dots; (2) the stress partition: viscous and Reynolds parts of the linear total stress; (3) (hidePortrait) the indicator y⁺dU⁺/dy⁺ with the 1/κ level.
- **Controls:** Re_τ (log) · κ · B · scaling (inner / outer) · wall (smooth / rough with y₀⁺) · wake strength (optional).
- **Equations:** (12.76)–(12.77) ⇒ $\bar\tau=\tau_0(1-2y/h)$; (12.80) $U^+=f(y^+)$; (12.81) $u_*^2=\tau_0/\rho$; (12.82) $U^+=y^+$; (12.84) $(U_\infty-U)/u_*=F(y/\delta)$; (12.87) $y^+df/dy^+=1/\kappa$; (12.88) $U^+=\frac1\kappa\ln y^++B$; (12.93) $U^+=\frac1\kappa\ln(y/y_0)$.
- **Mirrors** `WT.wall_units`, `WT.viscous_sublayer`, `WT.log_law`, `WT.spalding_uplus`, `WT.composite_profile`, `WT.stress_partition`, `WT.layer_name`, `WT.log_law_indicator`, `WT.rough_wall_log_law`.
- derivations: D17, D18, D19, D20 · depth features: explain (u_*, l_ν and δ⁺ from a wall stress in air or water; U⁺ at your y⁺ from each law; where sublayer and log lines cross; the extent of the log layer in decades; C_f from the friction law; the reading), code, **linked views**, **presets** (Re_τ = 180, 1000, 5200, atmospheric surface layer 10⁶, rough wall), **status** ("y⁺ = 12: buffer layer — neither law holds"), **inspector** (click the profile: the wall-unit arithmetic), **modes** (inner / outer scaling) · follows: APS (crosshair readouts across linked windows, numbered live derivation) · **aha:** the logarithm is what is left when the answer may depend on neither the viscous length nor the outer length.

### E8 · mixing_length_closure
- A: C12 (also shows N162, N167, N169–N172, N181 channel, N132 laminar comparison, N176 ν_T of k–ε named) · **Confusion removed:** "what does 'modelling the Reynolds stress' mean in practice, and how can one constant give the whole profile?" — replace −⟨uv⟩ by l_T²(dU/dy)² with l_T = κy, keep the exact linear total stress, and the mean profile follows by one integration; κ sets the log slope, the wall damping sets the intercept.
- **Why interactive:** sliding κ rotates the log line and sliding the damping constant A⁺ shifts it, with the implied intercept B read off live; switching the model off returns the laminar parabola at the same pressure gradient — the reader sees exactly what the closure adds and which constant controls what.
- **Stage:** (1) the channel profile U(y) from the model against the laminar parabola at the same pressure gradient (ghost) and public DNS points; (2) U⁺(y⁺) semi-log with the slope 1/κ and the intercept B marked; (3) (hidePortrait) eddy viscosity ν_T/ν and mixing length across the channel.
- **Controls:** κ · A⁺ (damping; 0 = none) · Re_τ · model on/off (laminar) · outer cap of l_T (optional).
- **Equations:** (12.94) $\overline{u_iu_j}=\tfrac23\bar e\delta_{ij}-\nu_T(\partial U_i/\partial x_j+\partial U_j/\partial x_i)$; (12.98) $\nu_T\sim l_Tu_T$; $-\overline{uv}=l_T^2(dU/dy)^2$; (12.100); $\nu\,dU/dy+\kappa^2y^2(dU/dy)^2=\tau_0/\rho$; $dU^+/dy^+=2/(1+\sqrt{1+4\kappa^2y^{+2}})$; (12.101) $U/u_*\cong\frac1\kappa\ln y+\mathrm{const}$.
- **Mirrors** `ch12.mixing_length_wall_profile`, `ch12.mixing_length_intercept` (§8), `ch12.channel_mixing_length`, `ch12.mixing_length_eddy_viscosity`, `core.laminar.channel_flow`.
- derivations: D21 · depth features: explain (the quadratic at your y⁺ with numbers and its positive root; ν_T/ν there; the intercept B implied by your κ and A⁺ against the measured one; bulk velocity and C_f against laminar at the same pressure gradient; the reading), code, **linked views**, **presets** (no damping — intercept too low, van Driest 26, laminar, Re_τ = 5200), **status** ("κ = 0.41, A⁺ = 26: B = 5.3 — matches smooth-wall data" / "no damping: B = −1.2"), **inspector** (click a y⁺: the quadratic's arithmetic) · follows: `stride_padding_playground` (classic presets, formula with numbers, badges) + FDV (Explain) · **aha:** a closure is a guess for one length; κ turns the profile into a logarithm and the near-wall damping decides where that logarithm sits.

### E9 · stratified_surface_layer
- A: C14, C15 (also shows N182, N183 Rf_cr, R05 conventions, R06 Ri with the in-situ gradient, N184, N185, N186–N189, N159 z₀, N54 heat-flux units) · **Confusion removed:** "when does stratification kill turbulence, what is the difference between Rf and Ri, and what does the Monin–Obukhov length measure?" — Rf compares what buoyancy removes with what shear supplies; in the surface layer it equals z/L_M, so L_M is the height where buoyancy catches up with shear; Ri is its measurable cousin, larger by the turbulent Prandtl number.
- **Why interactive:** dragging the surface heat flux through zero flips the sign of L_M, bends the wind profile to the other side of the neutral logarithm and moves the Rf = Rf_cr height up and down the column, while the stability badge states the regime in both lapse-rate conventions — three numbers and a profile that only make sense together.
- **Stage:** (1) the surface layer: wind profile U(z) on a semi-log height axis with the neutral logarithm as a ghost, the height ∣L_M∣ marked, forced / free convection or stable zones shaded; (2) T(z) and θ(z) side by side with the adiabat (both conventions in the legend); (3) (hidePortrait) budget bars at the cursor height: shear production, buoyancy, dissipation, with Rf and Ri.
- **Controls:** friction velocity u_* (or wind at 10 m) · surface heat flux H (both signs) · roughness length z₀ (chips: sea, grass, forest) · cursor height z · convention toggle (Kundu / meteorology).
- **Equations:** (12.106); (12.107) $\mathrm{Rf}=\dfrac{-g\alpha\overline{wT'}}{-\overline{uw}\,dU/dz}$; (12.108) $\mathrm{Ri}=N^2/(dU/dz)^2$; (12.109) $\mathrm{Ri}=(\nu_T/\kappa_T)\mathrm{Rf}$; (12.110) $L_M=-u_*^3/(\kappa\alpha g\overline{wT'})$; (12.111) $\mathrm{Rf}=z/L_M$; $U=\dfrac{u_*}\kappa\Big[\ln\dfrac z{z_0}+5\dfrac z{L_M}\Big]$; $N^2=g\alpha(dT/dz-\Gamma_a)$ with $\Gamma_a\approx-9.8$ K/km (Kundu) $=g\alpha(\Gamma_d-\Gamma)$ (meteorology).
- **Mirrors** `ch12.surface_layer_state` (§8), `ch12.monin_obukhov_from_fluxes`, `ch12.flux_richardson_surface_layer`, `ch12.gradient_richardson_thermal`, `ch12.turbulence_regime`, `WT.surface_layer_wind`, `core.stratification.lapse_rate_stability`.
- derivations: D24, D25 · depth features: explain (⟨wT′⟩ from H; L_M with every factor; Rf and Ri at your height; the height where Rf = ¼; the wind at 10 m against neutral; the same lapse rate in both conventions; the reading), code, **linked views**, **presets** (sunny afternoon, neutral overcast, clear calm night, strong wind night, sea surface), **status** (verdict → criterion text in the chosen convention → the other convention's relation), **terms** (production, buoyancy, dissipation), **notes** (what a climate model's surface scheme does with these numbers) · follows: FDV (Explain panel with regime-dependent reading) with the ch01 `parcel_stability` badge · **aha:** L_M is a height — below it the wind makes the turbulence, above it the heat flux makes (or kills) it.

### E10 · taylor_dispersion
- A: C16 (also shows N196 particle paths, N199, N208 exponential closed form, N211 random walk, N212 smoke plume, N213–N215 D_T, N216, N06 Richardson's law named, slips #12 and #13) · **Confusion removed:** "why does a puff spread in proportion to time at first and to the square root of time later, and why is 'the' eddy diffusivity not a constant?" — while a particle remembers its velocity it flies straight (X ∝ t); after a few memory times its steps are independent and it random-walks (X ∝ √t); D_T is the running integral of the velocity correlation, so it grows before it saturates.
- **Why interactive:** the particles spread on one clock with the ⟨X²⟩ curve on log axes changing slope from 2 to 1 as t passes Λ_t; dragging Λ_t moves the bend, and a constant-diffusivity ghost of the same late-time width shows how wrong Fickian diffusion is near the source; the plume mode maps t to x/U.
- **Stage:** (1) a cloud of Langevin particles from a point source with the ±X_rms envelope (ballistic and diffusive asymptotes dashed), or the time-averaged plume behind a chimney; (2) ⟨X²⟩(t) on log–log axes: particles (dots), Taylor's formula (bold), the two limits, the constant-D ghost; (3) (hidePortrait) r(τ) with the running area and D_T(t) rising to ⟨u²⟩Λ_t.
- **Controls:** u_rms · memory Λ_t · system (modes: puff in time / plume in a wind U) · correlation shape (exponential / Gaussian) · transport t.
- **Equations:** (12.117) $\frac{d}{dt}\overline{X_\alpha^2}=2\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$; (12.119) $\overline{X_\alpha^2}=2\overline{u_\alpha^2}\,t\int_0^t(1-\tau/t)r_\alpha\,d\tau$; (12.121) $X_{rms}=u_{rms}t$; (12.123) $X_{rms}=u_{rms}\sqrt{2\Lambda_tt}$; (12.125) $(R_n)_{rms}=L\sqrt n$; (12.127) $D_T=\overline{u_\alpha^2}\int_0^tr_\alpha\,d\tau$; (12.128) $D_T\cong\overline{u_\alpha^2}t$; (12.129) $D_T\cong\overline{u_\alpha^2}\Lambda_t$ for $t\gg\Lambda_t$.
- **Mirrors** `ch12.taylor_dispersion_exponential`, `ch12.taylor_dispersion_gaussian`, `ch12.eddy_diffusivity_exponential` (§8), `ch12.langevin_particles`, `ch12.dispersion_regime`, `ch12.smoke_plume_width`.
- derivations: D26, D27, D28 · depth features: explain (⟨X²⟩ at your t from the closed form with numbers; both limits and how far off each is; D_T now against its final value; the plume width at your distance; the reading), code, **linked views**, **transport** (play / step / scrub, end-of-run card with the measured slopes), **modes** (puff / plume), **presets** (short memory, long memory, t = Λ_t, chimney in a 5 m/s wind), **status** ("t = 0.3 Λ_t: ballistic — D_T still growing"), **terms** (D_T reached / still to come) · follows: AFE (one clock, modes, end-of-run summary, regime notes) · **aha:** turbulent spreading is a random walk whose step is set by how long a particle remembers its velocity — until then it is not diffusion at all.

### B1 · k_epsilon_decay (backup)
- A: C13 (also shows N175, N176–N178, N180) · **Confusion removed:** "where do the five k–ε constants come from?" — two of them are pinned by experiments anyone can picture: C_ε2 by how fast grid turbulence decays (n = 1/(C_ε2 − 1)), and the combination √C_μ(C_ε2 − C_ε1)σ_ε by the von Kármán constant of the log layer.
- Stage: ē(t) and ε̄(t) on log–log axes from the model with the decay exponent read off · ν_T(t) and the length scale ē^{3/2}/ε̄ · the implied κ as a dial against 0.41. Controls: C_ε2, C_ε1, C_μ, σ_ε, initial ē and ε̄. Equations (12.103)–(12.105) with their terms named. Mirrors `ch12.k_epsilon_decay`, `ch12.k_epsilon_loglayer_kappa`, `ch12.k_epsilon_eddy_viscosity`. Derivations D22, D23 · depth: linked views, transport, presets (standard constants, C_ε2 = 2.0, measured decay n = 1.3), status, terms (production, dissipation, transport of ē). Follows FDV. Built only if an explainer above fails review.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
Animations (5, `animate` + `show_animation`; ≤ 90 frames, FAST → 40; dpi 80):
| # | A | What moves | Why | Player |
|---|---|---|---|---|
| A1 | C01 (N20) | members of a decaying ensemble appear one by one while the running ensemble mean settles on the expected value; a sliding time average of one member beside it | what "N → ∞" means, and that the scatter shrinks as N^{−1/2} | frames (stop at N = 2, 4, 8, 64) |
| A2 | C02 (N30) | a copy of the signal slides past the original; the product is shaded and the correlation curve is drawn point by point, peaking at the delay | the lag correlation as an operation, not a formula | video |
| A3 | C04 (N50) | parcels exchanged across a mean shear, each dropping a dot in the (u, v) plane; the running ⟨uv⟩ converges to −⟨vℓ⟩ dU/dy | the sign of the Reynolds shear stress built up from single events | video |
| A4 | C07, C08 (N87) | the cascade ladder on a log axis as Re_L rises from 10³ to 10⁸: tiers are added at the small end, the spectrum's −5/3 range widens | scale separation grows with Re; the large-scale end does not move | frames (one frame per decade of Re) |
| A5 | C16 (N196) | Langevin particles leaving a point source with the ballistic (∝ t) and diffusive (∝ √t) envelopes; ⟨X²⟩ on log–log axes beside it | the change of regime at t ≈ Λ_t on one clock | video |

Interactive figures (plotly `slider_figure`, precomputed; ≤ 40 steps):
| # | A | What the slider controls | Why |
|---|---|---|---|
| F1 | C01 (N23) | window Δt: the sliding average of A e^{−t/τ} + B cos ωt against the true mean, with the two window factors | seed 1 on the published page: too short, whole period, too long |
| F2 | C02, C03 | memory time τ_c: r(τ) and S_e(ω) side by side, Λ_t rectangle and S(0) marked | the Fourier reciprocity |
| F3 | C07, C08 (N96) | Re_L: model spectra in Kolmogorov scaling with the −5/3 line and the inertial range shaded | seed 2: the range widens as Re^{3/4} |
| F4 | C09 | station x: raw profiles U(x, y) and the rescaled U/U_CL vs y/x for the plane jet (turbulent exponents; laminar as a second trace) | collapse and the two growth laws |
| F5 | C10, C11 (N149) | Re_τ: composite U⁺(y⁺) on a semi-log axis with layer bands and DNS points | seed 3: inner collapse, outer departure |
| F6 | C12 | damping constant A⁺ (and a κ dropdown): mixing-length U⁺(y⁺) with the implied intercept B | seed 4: what each closure constant does |
| F7 | C14, C15 (N189) | L_M (from stable through neutral to unstable): wind profile ln z vs U with Rf(z) and the verdict in both lapse-rate conventions in the trace names | the stratified surface layer on the published page |
| F8 | C16 | Λ_t: X_rms(t) and D_T(t) from Taylor's formula with both limits and a constant-D ghost | seed 5: t then √t |
A live `ipywidgets` cell (kernel only) pairs with F3 (free ΔU, L, ν) and one with F7 (free u_*, H, z₀). Static figures accompany every
A block — our analogues of Figs. 12.2–12.6, 12.8–12.13, 12.16–12.27, the three-routes spectrum figure, the f / g figure, the scale
table, the exponent table for Table 12.1, the k–ε decay (closed form vs `solve_ivp`), all generated by our code.

## 7. From-scratch moments
Each shows a hand-written version beside the tested function, followed by `assert np.allclose(...)` (statistical comparisons use the
same seed, so they agree to round-off).
| A | § | Hand-written | Compared with |
|---|---|---|---|
| **C01** | 12.3 | a loop over members accumulating Σu and Σu², then mean and variance; the product split ⟨ũṽ⟩ − ŪV̄ by hand | `TS.ensemble_average`, `TS.central_moment`, `TS.product_average_split` |
| **C02** | 12.4 | the direct-sum autocorrelation (a loop over lags) and Λ_t by `np.trapezoid` to the first zero | `TS.autocorrelation(method="fft")`, `TS.integral_scale` |
| **C03** | 12.4 | the cosine transform of R(τ) by `np.trapezoid`, and a periodogram from `np.fft.rfft` scaled by hand (Parseval asserted) | `TS.spectrum_from_correlation`, `TS.periodogram` |
| **C04** | 12.5 | the 2 × 2 Reynolds stress as a covariance computed with explicit sums; ⟨uv⟩ of the parcels against −⟨vℓ⟩ dU/dy | `TS.reynolds_stress`, `ch12.displaced_parcel_uv` |
| **C05** | 12.6 | ε̄ three ways for a Gaussian f: 30ν⟨u²⟩/λ_f², 15ν⟨u²⟩/λ_g² and −15ν⟨u²⟩f″(0) by finite differences | `ch12.dissipation_isotropic`, `ch12.isotropic_scales` |
| **C06** | 12.7 | production −⟨uv⟩ dU/dy and its integral by `np.gradient` and `np.trapezoid` on the model channel | `ch12.shear_production`, `ch12.tke_budget` |
| **C07** | 12.7 | the 2 × 2 exponent solve for η = ν^a ε̄^b with `np.linalg.solve`, then η for an atmospheric case | `ch12.kolmogorov_scales`, `core.dimensional.solve_exponents` |
| **C08** | 12.7 | the slope of a band of the model spectrum by `np.polyfit` on the logs (−5/3 within the stated band) | `ch12.fit_inertial_range`, `ch12.inertial_spectrum_1d` |
| **C09** | 12.8 | the momentum flux ρ∫U²dy by `np.trapezoid` at three stations (equal) and the volume flux (growing as x^{1/2}) | `ch12.jet_momentum_flux_per_span`, `ch12.plane_jet_volume_flux` |
| **C10** | 12.9 | u_*, l_ν, y⁺ and U⁺ from τ₀, ρ, ν in four lines for air over a plate | `WT.wall_units`, `WT.friction_velocity` |
| **C11** | 12.9 | κ and B from a straight-line fit of U⁺ against ln y⁺ (`np.polyfit`) on a composite profile | `WT.fit_log_law` |
| **C12** | 12.10 | the positive root of ν dU/dy + κ²y²(dU/dy)² = u_*² point by point, then `cumulative_trapezoid` | `ch12.mixing_length_wall_profile` |
| **C13** | 12.10 | an RK4 loop for dē/dt = −ε̄, dε̄/dt = −C_ε2 ε̄²/ē against the power law | `ch12.k_epsilon_decay` |
| **C14** | 12.11 | Rf and Ri by hand from a heat flux, a stress and two gradients (Kundu sign), then the same Ri from dθ/dz | `ch12.flux_richardson`, `ch12.gradient_richardson_thermal` |
| **C15** | 12.11 | L_M from u_* and H in one line; the wind at 10 m from the log-linear formula | `ch12.monin_obukhov_from_fluxes`, `WT.surface_layer_wind` |
| **C16** | 12.12 | a random-walk loop (cumulative sum of unit steps) giving ⟨R_n²⟩ = nL²; Taylor's integral by `np.trapezoid` for e^{−τ/Λ_t} | `ch12.random_walk`, `ch12.taylor_dispersion`, `ch12.taylor_dispersion_exponential` |
§12.1, §12.2 and §12.13 have no computable A item (covered inside C01 and C16). Every other section has at least one from-scratch moment.

## 8. Notes for the implementer
Functions and helpers that the A items' figures and the explainers need beyond analysis §4 (B/C functions optional):
1. `TS.correlation_spectrum_pair(kind, sigma, tau_c, omega0=0.0)` → callables r(τ), S_e(ω) and the exact Λ_t, λ_t, S(0) for kind ∈ {"exponential", "gaussian", "damped_cosine"} (book normalisation) — E2's closed forms, F2, and V1 targets for `autocorrelation` and `periodogram`. Also `TS.smooth_signal(n, dt, spectrum, seed)` (random-phase signal with a prescribed smooth spectrum) so that the Taylor microscale exists (D04).
2. `ch12.parcel_uv_expected(dUdy, l_rms, v_rms, correlation=1.0)` → −correlation · v_rms · l_rms · dU/dy, the expectation of `displaced_parcel_uv` (E3 parity rows; 5-standard-error test).
3. `ch12.channel_energy_budget(Re_tau, kappa, A_plus, n=400)` → y⁺, U⁺, −⟨uv⟩⁺, the mean-flow terms (pressure work, viscous dissipation, loss to turbulence) and the turbulence terms (production; dissipation + transport as the residual, labelled model) with their integrals across the channel; integral identity work = dissipation + production asserted. Cached table for E5.
4. `ch12.inertial_range_decades(Re_L)` = log₁₀(L/η) and `ch12.scale_table(case)` with our own cases (E4 presets and table).
5. `ch12.wrong_exponent_fluxes(x, n, m)` → momentum and volume flux of a profile family U_CL ∝ x^n, δ ∝ x^m (E6's wrong-exponent slider; flat only when m + 2n = 0); `ch12.free_shear_profile(flow, xi)` (a Gaussian in the right similarity variable, constants passed explicitly).
6. `ch12.mixing_length_intercept(kappa, A_plus=None)` → the additive constant B implied by the model (limit of U⁺ − (1/κ) ln y⁺), with and without van Driest damping (E8 status, F6); `WT.log_law_crossing(kappa, B)` → y⁺ where U⁺ = y⁺ meets the log law (E7).
7. `ch12.surface_layer_state(u_star, H, T, z0, z, rho, cp, kappa, beta=5.0, Pr_T=1.0, convention="kundu")` → dict: ⟨wT′⟩, L_M, Rf(z), Ri(z), φ_m, U(z), regime, the height where Rf = Rf_cr, and the verdict text in both conventions (one entry point for E9 and F7); `ch12.gradient_richardson_surface_layer(z, L_M, Pr_T, beta)`.
8. `ch12.eddy_diffusivity_exponential(t, u2, Lambda_t)` = ⟨u²⟩Λ_t(1 − e^{−t/Λ_t}) and `ch12.dispersion_local_slope(t, Lambda_t)` (d ln⟨X²⟩/d ln t, from 2 to 1) for E10's status and end-of-run card.
9. `ch12.book_slips()` → the table below (id, where, printed, corrected, how a test tells them apart), and `ch12.conventions()` → the symbol / normalisation table of decision 6 as a pandas frame for the opening block.
10. `ch12.free_shear_exponents(flow, return_equations=True)` → the two linear equations (invariant, growth law) as text plus their exact solution (D16's table and E6's Explain tab).
11. `TS.periodogram` returns angular frequency (or wavenumber) and a two-sided density by default; `one_sided=True` doubles it. Every spectrum function states its normalisation in the docstring and a Parseval test pins it.

Printed slips to carry as callouts (analysis §9; the item each one sits in):
| slip | where | printed | taught as | sits in |
|---|---|---|---|---|
| #1 | (12.54) | exponent of k₁ is +5/3 | −5/3 (units and the text's own name for the law) | C08, D12, E4 |
| #2 | (12.47) | the label under the left side names Ē | the turbulent ē | C06, D10 |
| #3 | (12.70) | source integral at y = 0 | at x = 0, as in (12.62) | N118 |
| #4 | after (12.67) | constants called C₃, C₄; later C₅, C₆ | C₃ and C₅; the tabulated pair multiplies (12.72) and (12.73) | N115, D15 |
| #5 | (12.74) | exponential family δ ~ e^{ax}, U_CL ~ e^{−ax} offered | its middle coefficient is zero and it breaks (12.62); the power family with m + 2n = 0 is the valid one | N125 |
| #6 | (12.97) | ∂P/∂x_j | ∂P/∂x_i (free index) | N165 |
| #7 | p. 590 | proof cited as Exercise 12.31 | Exercise 12.32 | R03, D17 |
| #8 | p. 595 | 2ν⟨u_jS′_ij⟩ in the modelled transport | 2ν⟨u_iS′_ij⟩ as in (12.47) | N175 |
| #9 | p. 595 | "five" constants, six printed entries | five: C_μ, C_ε1, C_ε2, σ_e, σ_ε | N178 |
| #10 | (12.75), (12.106) | triple correlation as ½⟨ev⟩ and ⟨ew⟩ | ½⟨u_i²u_j⟩ throughout | N128, N182, D10 |
| #11 | p. 604 | "(11.119)" | (12.119) | N207, D27 |
| #12 | (12.129) | condition t ≪ Λ_t | t ≫ Λ_t | N215, D28, E10 |
| #13 | Fig. 12.27 caption | √x near the source, x far away | linear near, √x far (as the text and (12.121), (12.123) say) | N212, E10 |
| #14 | Exercise 12.18a | cites (12.39) for R_ij | (12.40) | N69, D07 |
| #15 | (12.112) | molecular transport κ∂⟨T′²⟩/∂z | κ∂(½⟨T′²⟩)/∂z | N190 |
Traps that are not slips: (12.16) without an absolute value (N29, D02); r_α(t′ − t) in (12.117) (N200, D26); the three different gradient
moments of (12.43) (D08); ν set as an italic v in (12.109).

Open decisions for the orchestrator: (a) the fitted coefficients of the zero-pressure-gradient formulas (N150) — displayed equations
(default) or private book numbers; (b) `FREE_SHEAR_CONSTANTS` has no default set until the verifier confirms a public table (examples pass
labelled illustrative constants; exponents are exact and public); (c) re-export `ch11.gradient_richardson` from `core.stratification`.

## 9. Implementation guidance (conventions, runtime budget, caching)
- **Conventions in code.** `kappa` = von Kármán constant, `kappa_th` = thermal diffusivity; `e` = turbulent kinetic energy (the "k" of
  k–ε), `k1` and `K` = wavenumbers, `k_th` = conductivity; `u2` = one-component variance; spectra two-sided in angular frequency or
  wavenumber unless `one_sided=True`; `h` = full channel height, `delta` = half-height or boundary-layer thickness; `dTdz` in
  Kundu's sign with `Gamma_a` negative, verdict strings through `core.stratification.lapse_rate_stability` in both conventions;
  heat flux `wT` [K m/s] and `H` [W/m²] positive upward. κ, B, the wake strength and the free-shear constants are required keywords.
- **Seeds and tolerances.** Every sampled quantity uses `np.random.default_rng(seed)` with the seed in the docstring; statistical
  assertions use 5 standard errors; synthetic fields are labelled kinematic (no cascade) in every caption.
- **Runtime budget (< 5 min on Colab CPU).** Everything is closed-form or 1-D except: the 3-D synthetic field (FAST 48³, full 64³,
  < 2 s, cached), 10⁴ Langevin particles × 10³ steps (FAST 2000 × 400), the nine sympy engines (each < 10 s — cache their results
  to `outputs/ch12/cache/` with P252 and show the printed steps), the mixing-length channel at five Re_τ (cached), and the optional
  `k_epsilon_channel` (not needed by any A figure: skip in the notebook unless cached). Animations ≤ 90 frames (FAST 40), slider
  figures ≤ 40 steps × ≤ 4 traces. Expected total ≈ 150–200 s.
- **Reference data.** `reference/ch12/` holds the few columns of the Lee & Moser (2015) channel profiles the figures use, with
  `SOURCES.md`; the verifier reads the numbers from the files, not from the analysis note.
- **Explainer parity.** Each explainer's `selftest()` mirrors the fluidpy functions named above at ≥ 2 points; stochastic stages use
  the closed-form expectation for parity (E1 window factors and standard error, E3 `parcel_uv_expected`, E10 the exponential closed
  form), never a sampled value.
