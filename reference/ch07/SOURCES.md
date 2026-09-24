# Chapter 7 reference data — sources

Produced by `reference/ch07/make_refs.py` → `benchmarks.json`. Public sources only; every page was fetched and read on
the date given. Book-quoted numbers and printed closed forms are **not** here (they live in the git-ignored
`tests/book_values_ch07.json`). Label use in `tests/test_ch07.py`: published *numbers* our code reproduces from its own
implementation are V5; a published closed form reproduced identically is a V1 "form cross-check" (not V5).

| value / form | what it is | source | URL | how obtained | verified | label |
|---|---|---|---|---|---|---|
| `fenton_mckee_1990` | explicit wavelength L = L₀{tanh[(2π√(d/g)/T)^{3/2}]}^{2/3} (Eq. 21); max wavelength error 1.7 % ("always better than 1.7 %"); optimal exponent ν = 1.49 (rounded to 3/2) | J. D. Fenton & W. D. McKee, "On calculating the lengths of water waves", Coastal Engineering 14, 499–513 (1990), doi:10.1016/0378-3839(90)90032-R | https://johndfenton.com/Papers/Fenton90c+McKee-On-calculating-the-lengths-of-water-waves.pdf | PDF fetched, text extracted with pymupdf; abstract and p. 507 read | 2026-09-24 | V5 |
| `fenton_2006_note` | the same formula "accurate to within 1.5 %"; Guo's formula with exponent 5/2 "maximum error of about 0.7 %, about half the error of (2)"; both exact in both limits | J. D. Fenton, "A note on two approximations to the linear dispersion relation for surface gravity water waves" (Oct. 26, 2006) | https://johndfenton.com/Papers/Dispersion-Relation.pdf | PDF fetched, text extracted | 2026-09-24 | V5 (loose statements) |
| `guo_2002` | kh = x²(1 − e^{−x^β})^{−1/β}, β = 2.4908; maximum relative error 0.75 % | J. Guo, "Simple and explicit solution of wave dispersion equation", Coastal Engineering 45, 71–74 (2002), doi:10.1016/S0378-3839(02)00039-X | https://www.sciencedirect.com/science/article/abs/pii/S037838390200039X | **search-confirmed**: the abstract's 0.75 % and β = 2.4908 read from search-engine summaries; ScienceDirect, ResearchGate and academia.edu returned HTTP 403 | 2026-09-24 | V5 (search-confirmed) |
| `capillary_wave` | ω² = \|k\|[((ρ − ρ′)/(ρ + ρ′))g + σk²/(ρ + ρ′)]; air–water c_m = 0.23 m/s at λ_m = 1.7 cm; λ_m = 2π√(σ/((ρ − ρ′)g)) | Wikipedia, "Capillary wave" | https://en.wikipedia.org/wiki/Capillary_wave | page read; formulas and numbers quoted (2 significant figures) | 2026-09-24 | V5 (numbers) + V1 form |
| `iapws_surface_tension` | σ(20 °C) = 72.74 mN/m (input property for the capillary numbers; `ch01.surface_tension_water`) | IAPWS R1-76(2014) | https://iapws.org/public/documents/CH-L9/Surf-H2O-2014.pdf | reused from ch01 (V5 there) | 2026-09-24 | input |
| `stokes_limit` | limiting deep-water Stokes wave H/λ = 0.1410633 ± 4·10⁻⁷, crest angle 120°; third-order c = (1 + ½(ka)²)√(g/k) | HandWiki, "Physics:Stokes wave" (the page attributes the number to Schwartz & Fenton); crest angle 2π/3 also in Dyachenko, Lushnikov & Korotkevich, Stud. Appl. Math. 137, 419–472 (2016), arXiv:1507.02784 | https://handwiki.org/wiki/Physics:Stokes_wave ; https://arxiv.org/abs/1507.02784 | pages read; number, angle and formula quoted | 2026-09-24 | V5 (constant) + V1 form (speed) |
| `stokes_drift` | deep water ū_S ≈ ωka²e^{2kz}; "at z = −λ/4 it is about 4 % of its value at the mean free surface" | Wikipedia, "Stokes drift" | https://en.wikipedia.org/wiki/Stokes_drift | page read; formula and sentence quoted (the finite-depth form is not on the page) | 2026-09-24 | V5 (4 %) + V1 form |
| `belanger` | y₂/y₁ = ½(√(1 + 8Fr₁²) − 1); ΔE = (y₂ − y₁)³/(4y₁y₂) | Wikipedia, "Hydraulic jumps in rectangular channels" (Bélanger 1828) | https://en.wikipedia.org/wiki/Hydraulic_jumps_in_rectangular_channels | page read; formulas quoted | 2026-09-24 | V1 form |
| `cnoidal` | dimensional KdV; cnoidal η = η₂ + H cn²((x − ct)/Δ \| m) with η₂, Δ, c, λ = 2ΔK(m); solitary limit; Ursell U = Hλ²/h³ | Wikipedia, "Cnoidal wave" | https://en.wikipedia.org/wiki/Cnoidal_wave | page read; formulas quoted | 2026-09-24 | V1 form |
| `st_andrews_cross` | four beams at arccos(ω/N) from the vertical from an oscillating body | D. E. Mowbray & B. S. H. Rarity, J. Fluid Mech. 28, 1–16 (1967) | https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/theoretical-and-experimental-investigation-of-the-phase-configuration-of-internal-waves-of-small-amplitude-in-a-density-stratified-liquid/BD4159D2420DA60E9640520CD06B68AB | citation confirmed by search | 2026-09-24 | pedigree |

## Measured against the sources (our code, 2026-09-24)

| quantity | published | ours | comment |
|---|---|---|---|
| Fenton–McKee max wavelength error | 1.7 % (1990 primary) · 1.5 % (2006 note) | 1.658 % (in λ), 1.631 % (in kd) | the primary 1990 bound holds and rounds to 1.7 %; the 2006 note's "within 1.5 %" is not reproduced in λ or kd (it is ≈ the *period* error, 1.49 %) |
| Fenton–McKee optimal exponent ν | 1.49 | 1.483 (min–max in λ), 1.486 (in kd) | 0.5 % below the printed value (grid/measure difference; within the 1 % benchmark tolerance) |
| Guo max error, β = 2.4908 | 0.75 % | 0.753 % (in kh) | reproduced to printed precision |
| Guo max error, β = 5/2 (`guo_kh`) | "about 0.7 %" (2006 note) | 0.789 % (in kh) | loose statement; the ratio to Fenton–McKee is 0.48 ("about half" ✓) |
| capillary minimum (σ from IAPWS at 20 °C, ρ = 998.2 kg/m³, g = 9.81) | 0.23 m/s at 1.7 cm | 0.2312 m/s at 1.712 cm | agree to the printed 2 significant figures |
| deep Stokes drift at z = −λ/4 | about 4 % | e^{−π} = 4.32 % | agrees to the printed precision |

Attribution note: `fluidpy.ch07_gravity_waves.STOKES_LIMIT_STEEPNESS`'s docstring (and analysis §8) credit the number to
Dyachenko, Lushnikov & Korotkevich (2016) "as reported by HandWiki"; the HandWiki page itself attributes it to Schwartz &
Fenton, and the Dyachenko et al. abstract gives the 120° crest angle but no steepness value. The value is right; the
citation text should be corrected (open item in `reports/ch07_verification.md`).
