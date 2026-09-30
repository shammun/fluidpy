# Chapter 9 reference data — sources

Produced by `reference/ch09/make_refs.py` → `benchmarks.json`. Public sources only; every page or PDF was fetched and
read on 2026-09-30 (PDF text extracted with PyMuPDF where the fetch returned a binary). Book-quoted numbers and printed forms
are **not** here (they live in the git-ignored `tests/book_values_ch09.json`). Label use in `tests/test_ch09.py`: a published
*number* our code reproduces from its own implementation is V5; a published closed form reproduced identically is a V1 "form
cross-check" (not V5); secondary (encyclopaedic) sources get the looser tolerance and say so.

| value / form | what it is | source | URL | how obtained | verified | label |
|---|---|---|---|---|---|---|
| `belden_falkner_skan` | Falkner–Skan wall shear κ = f″(0) (their normalisation f‴ + ff″ + β(1 − f′²) = 0, η = y√((m+1)U/2νx), m = β/(2 − β)): β = 0.5 → 0.927680039836653; β = 0 → 0.469599988361; β = −0.12 → 0.28176052424; β = −0.198837735 → 0 (separation); reversed branch β = −0.12 → −0.1429351943576, β = −0.02 → −0.065168585542904 | E. R. Belden, Z. A. Dickman, S. J. Weinstein, A. D. Archibee, E. Burroughs, N. S. Barlow, "Asymptotic Approximant for the Falkner–Skan Boundary-Layer equation", arXiv:1907.09912 (2019), Fig. 2 caption and Eqs. (1)–(2) | https://arxiv.org/abs/1907.09912 | PDF downloaded, text extracted; RK4 shooting with domain-length convergence quoted by the authors | 2026-09-30 | V5 |
| `blasius_wikipedia` | Blasius f″(0) = 0.332057336215196 (small-η coefficient); δ* ≈ 1.72, θ ≈ 0.665 (rounded); the page also prints 0.332043934904293 for the boundary-condition conversion (a different, less converged number on the same page) | Wikipedia, "Blasius boundary layer" | https://en.wikipedia.org/wiki/Blasius_boundary_layer | page read; digit strings quoted | 2026-09-30 | V5 (constant), rounded thicknesses V1-level only |
| `hiemenz_weidman_turner` | Hiemenz stagnation-point flow F‴ + FF″ − F′² + 1 = 0: F″(0) = 1.232588 (σ = 0 of a stretching-plate expansion) | P. Weidman & M. R. Turner, "Stagnation-point flows with stretching surfaces: A unified formulation and new results" (Univ. Colorado / Univ. Surrey preprint), §5.1 Eqs. (5.3), (5.5) | https://personalpages.surrey.ac.uk/m.turner/Crane_Revision_FINAL.pdf | PDF downloaded, text extracted, §5.1 read | 2026-09-30 | V5 (7 significant figures) |
| `hiemenz_delta_star_wikipedia` | Hiemenz δ* = 0.6479 δ (δ = √(ν/a)) | Wikipedia, "Stagnation point flow" | https://en.wikipedia.org/wiki/Stagnation_point_flow | page read; only this number is printed | 2026-09-30 | V5 (4 s.f.) |
| `bickley_jet_wikipedia` | Bickley plane jet: u = 0.4543 (M²/νρ²x)^{1/3} sech²ξ, ξ = 0.2752 (M/ν²ρ)^{1/3} y/x^{2/3}, Q = 2ρ∫₀^∞u dy = 3.3019 (Mνρ²x)^{1/3}, v = 0.5503 (Mν/ρx²)^{1/3}(2ξ sech²ξ − tanh ξ); M = momentum flux per unit span (book's J) | Wikipedia, "Bickley jet" (citing Bickley, Phil. Mag. 23 (1937) 727–731) | https://en.wikipedia.org/wiki/Bickley_jet | page read; formulas quoted | 2026-09-30 | V5 (4–5 s.f.) |
| `thwaites_agrawal` | Thwaites' fit L(m) ≈ 0.45 + 6m with m = −λ, θ² = 0.45ν/U_e⁶ ∫U_e⁵ dr, laminar separation at m ≈ 0.09 | R. Agrawal, S. T. Bose, K. P. Griffin, P. Moin, "An extension of Thwaites' method for turbulent boundary layers", arXiv:2310.16337, Eqs. (2.3)–(2.5) | https://arxiv.org/abs/2310.16337 | PDF downloaded, text extracted, §2 read | 2026-09-30 | V5 (form) |
| `skin_friction_wikipedia` | laminar c_f = 0.664 Re_x^{−1/2}; turbulent local c_f = 0.0576 Re_x^{−1/5} (Prandtl one-seventh-power law) ⇒ plate-mean 0.0720 Re_L^{−1/5} (× 5/4) | Wikipedia, "Skin friction drag" | https://en.wikipedia.org/wiki/Skin_friction_drag | page read | 2026-09-30 | V5 secondary (3 %) |
| `karman_ratio` | staggered point-vortex street stable only at b/a = cosh⁻¹(√2)/π ≈ 0.281 (Kármán 1911–12) | Horváth et al., J. Geophys. Res. Atmos. 125 (2020), doi:10.1029/2019JD032121; A Stroll down Kármán Street, Resonance 10(8), ias.ac.in; arXiv:1807.00203 | https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JD032121 | search-result text quoted (page itself not opened) | 2026-09-30 | V5 corroboration (the exact value is V1/V2) |
| `morrison_sphere_drag` | sphere C_D correlation — reused from ch04 | F. A. Morrison (2016), Michigan Tech — `reference/ch04/SOURCES.md` | https://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdf | reused (verified 2026-09-23) | 2026-09-23 | V5 (correlation) |

## Looked for and not found (recorded, not used)

| what | outcome |
|---|---|
| Glauert (1956) wall-jet numerical constants (J. Fluid Mech. 1, 625) | not fetched; no independent published number for the wall-jet constants is asserted. The wall jet is verified by V2 (sympy), V3 (two routes: IVP vs the implicit relation), V4 (invariant), V7 (scaling f_∞³ = 72 f″(0)). |
| Lienhard (1966) cylinder C_D(Re) table (Fig. 9.21) | PDF not accessed; `cylinder_cd_schematic` stays a qualitative schematic (no numbers asserted). |
| Roshko (1954) Strouhal numbers | secondary sources only; `strouhal_of_re` plateau tested as a band, not a number. |
| Howarth (1938) retarded-flow and cylinder separation numbers (0.12 L; 109.6°) | not confirmed; nothing asserted (Open item). |
| Primary Schlichting/White coefficient of the turbulent flat-plate law | not accessed; only the Wikipedia secondary value (0.0576 local, 0.0720 mean) — the code's 0.074 is asserted to 3 % of it. |

## Measured against the sources (our code, 2026-09-30)

| quantity | published | ours | relative difference |
|---|---|---|---|
| Blasius f″(0) (Wikipedia small-η coefficient) | 0.332057336215196 | 0.3320573362151941 (Töpfer IVP) | 5.7 × 10⁻¹⁵ |
| Blasius κ = f″(0)√2 (Belden normalisation) | 0.469599988361 | 0.469599988361036 | 7.7 × 10⁻¹⁴ (their 12 printed digits) |
| FS β = 0.5 (n = 1/3) κ | 0.927680039836653 | 0.927680039836652 | 1.1 × 10⁻¹⁵ |
| FS β = −0.12 (n = −0.0566) κ | 0.28176052424 | 0.281760524240411 | 1.5 × 10⁻¹² |
| FS separation β | −0.198837735 | −0.198837735047 | 2.4 × 10⁻¹⁰ (9 printed digits) |
| FS reversed branch β = −0.12 κ | −0.1429351943576 | −0.1429351943568 | 5.5 × 10⁻¹² |
| FS reversed branch β = −0.02 κ | −0.065168585542904 | −0.0651685855429459 | 6.4 × 10⁻¹³ |
| Hiemenz F″(0) | 1.232588 | 1.2325876568 | 2.8 × 10⁻⁷ (7 printed digits) |
| Hiemenz δ*/δ | 0.6479 | 0.6479004744 | 7.3 × 10⁻⁷ |
| Bickley u₀ coefficient | 0.4543 | 0.454280 | 4.4 × 10⁻⁵ (rounding of 4 printed digits) |
| Bickley ξ coefficient | 0.2752 | 0.275161 | 1.4 × 10⁻⁴ (rounding of 4 printed digits) |
| Bickley Q coefficient | 3.3019 | 3.301927 | 8.3 × 10⁻⁶ |
| Bickley v coefficient (√6/3 · C^{−1/3}) | 0.5503 | 0.550321 | 3.9 × 10⁻⁵ |
| Kármán ratio b/a | 0.281 | 0.280550 | 1.6 × 10⁻³ (3 printed digits) |
| plate-mean turbulent coefficient | 0.0720 (= 5/4 × 0.0576) | 0.074 (code default) | +2.8 % (secondary source; inside the 3 % band) |
