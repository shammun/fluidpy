# Chapter 8 reference data — sources

Produced by `reference/ch08/make_refs.py` → `benchmarks.json`. Public sources only; every page was fetched and read on
the date given. Book-quoted numbers and printed closed forms are **not** here (they live in the git-ignored
`tests/book_values_ch08.json`). Label use in `tests/test_ch08.py`: published *numbers* our code reproduces from its own
implementation are V5; a published closed form reproduced identically is a V1 "form cross-check" (not V5).

| value / form | what it is | source | URL | how obtained | verified | label |
|---|---|---|---|---|---|---|
| `san_andres_slider` | inclined slider bearing: w = (6μUL²B/h₂²)W(α), W(α) = [ln α + 2(1 − α)/(1 + α)]/(1 − α)², α = h₁/h₂ (inlet/exit); α_opt = 2.1889, W(α_opt) = 0.0267; P_max(α) = (α − 1)/(4α(1 + α)), maximal (0.043) at α = 2.414 | L. San Andrés, *Modern Lubrication*, Notes 2 Appendix "One-dimensional fluid film bearings", Texas A&M (© 2009, rev. Aug 2012), pp. 5–6 | https://rotorlab.tamu.edu/me626/Notes_pdf/Notes02_App_1D_bearings.pdf | PDF downloaded, text extracted with PyMuPDF, Eqs. (7)–(8b) and the MATHCAD optimum lines read | 2026-09-28 | V5 |
| `huppert_planar_current` | planar viscous gravity current h_t − β(h³h_x)_x = 0, ∫₀^{x_N} h dx = A: h = (3/10)^{1/3}η_N^{2/3}(A²/β)^{1/5}t^{−1/5}(1 − y²)^{1/3}, x_N = η_N(βA³)^{1/5}t^{1/5}, η_N = [(1/5)(3/10)^{1/3}π^{1/2}Γ(1/3)/Γ(5/6)]^{−3/5} = 1.411… | T. V. Ball & H. E. Huppert, "Similarity solutions and viscous gravity current adjustment times" (preprint), Appendix (a), Eqs. (A1)–(A7); citation of record H. E. Huppert, J. Fluid Mech. 121, 43–58 (1982) | https://warwick.ac.uk/fac/sci/maths/people/staff/tball/publications/similarity_preprint.pdf | PDF downloaded, text extracted with PyMuPDF, Appendix (a) read (the preprint's similarity variable prints t^{1/5} where t^{−1/5} is meant; (A4)–(A6) are consistent and are what we use) | 2026-09-28 | V5 |
| `oseen_drag` | Oseen F = 6πμau(1 + (3/8)Re), Re = ρua/μ (**radius**); C_d = (12/Re)(1 + (3/8)Re); Proudman–Pearson F = 6πμaU(1 + (3/8)Re + (9/40)Re² ln Re + O(Re²)) | Wikipedia, "Oseen equations" | https://en.wikipedia.org/wiki/Oseen_equations | page read; formulas quoted | 2026-09-28 | V1 form (exact rationals of the theory) |
| `stokes_law` | F = 6πμRv; settling v = (2/9)(ρ_p − ρ_f)gR²/μ; p(r, θ) = −(3μRu/2) cos θ/r² | Wikipedia, "Stokes' law" | https://en.wikipedia.org/wiki/Stokes%27_law | page read; formulas quoted | 2026-09-28 | V1 form |
| `stokes_second_problem` | u = Ue^{−√(ω/2ν)y} cos(ωt − √(ω/2ν)y); penetration depth δ = √(2ν/ω) | Wikipedia, "Stokes problem" | https://en.wikipedia.org/wiki/Stokes_problem | page read; formula quoted | 2026-09-28 | V1 form |
| `hagen_poiseuille` | Δp = 8μLQ/(πR⁴); u = (G/4μ)(R² − r²); Darcy f = 64/Re, Re = ρvd/μ (mean velocity, diameter) | Wikipedia, "Hagen–Poiseuille equation" | https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_equation | page read; formulas quoted | 2026-09-28 | V1 form |
| `taylor_couette` | v_θ = Ar + B/r, A = Ω₁(μ − η²)/(1 − η²), B = Ω₁R₁²(1 − μ)/(1 − η²), μ = Ω₂/Ω₁, η = R₁/R₂; Rayleigh: stable iff (rv_θ)² increases outward (μ > η² for co-rotation) | Wikipedia, "Taylor–Couette flow" | https://en.wikipedia.org/wiki/Taylor%E2%80%93Couette_flow | page read; formulas quoted | 2026-09-28 | V1 form |
| `elementary_charge` | e = 1.602 176 634 × 10⁻¹⁹ C (exact) | NIST, CODATA 2022 recommended values | https://physics.nist.gov/cgi-bin/cuu/Value?e | page read | 2026-09-28 | V5 |
| `morrison_sphere_drag` | sphere C_D correlation (24/Re for Re < 2) — reused | F. A. Morrison (2016), Michigan Tech — `reference/ch04/SOURCES.md` | https://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdf | reused (verified 2026-09-23) | 2026-09-23 | V5 (correlation) |
| `ussa1976_sutherland` | Sutherland β = 1.458 × 10⁻⁶ kg/(m s K^{1/2}), S = 110.4 K — reused | NASA-TM-X-74335 (USSA-1976) — `reference/ch01/` | (see `reference/ch01/SOURCES.md`) | reused (verified 2026-09-12) | 2026-09-12 | V5 (input property) |

Not a benchmark (context only): the Royal Society Phil. Trans. A 378 (2020) 20190521 Stokes-layer paper (search result
only); the pipe transition "Re ≈ 2300" (no primary page read — the book's 2000–3000 band stays V6).

## Measured against the sources (our code, 2026-09-28)

| quantity | published | ours | comment |
|---|---|---|---|
| optimum inlet/exit ratio K_opt | 2.1889 | 2.18870 (`slider_optimum_taper`) | −9 × 10⁻⁵ relative (printed 5 figures agree to 2 × 10⁻⁴ absolute) |
| maximum dimensionless load W(K_opt) | 0.0267 | 0.026707 | agrees to the printed 3 figures |
| San Andrés W(K) vs our exact load | identity | 1e-12 on a sweep K = 1.05 … 6 | same closed form after K = 1 + α |
| peak-pressure optimum K, P_max | 2.414, 0.043 | 2.4142 (= 1 + √2), 0.04289 | agrees to printed precision |
| η_N (planar viscous current) | 1.411… | 1.411245 | agrees to printed precision; Γ-function formula reproduced to 1e-14 |
| Oseen / Proudman–Pearson with Re_a = Re/2 | radius forms | identical to 1e-14 | the book's 3/16 on Re = 2aU/ν is Wikipedia's 3/8 on Re_a |
| synthetic Millikan (seed 0, 40 drops, 1 % speed noise) | e (exact) | relative error 5.1 × 10⁻⁴ | all 40 integers recovered |
| Morrison correlation between Stokes and Oseen, 0.1 ≤ Re ≤ 5 | book's claim | holds at 40 log-spaced points | correlation, not data |
