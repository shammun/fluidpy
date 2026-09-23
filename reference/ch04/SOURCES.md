# reference/ch04 — sources of the public benchmark data

`benchmarks.json` is written by `make_refs.py`. Chapter 4 (Conservation Laws) needs no datasets: the evidence in
`tests/test_ch04.py` is overwhelmingly Tier 1 (exact solutions of the Navier–Stokes equations, sympy re-derivations,
conservation budgets, convergence orders). The public sources below pin a few **published numbers** and **forms** the
chapter's examples touch. Every row was read from the cited page by the math-verifier on the date given (page fetched
and its text read). No value comes from the textbook; book-quoted numbers stay in the git-ignored
`tests/book_values_ch04.json`.

| value / key in `benchmarks.json` | what it is | source | URL | how obtained | verified |
|---|---|---|---|---|---|
| `morrison_sphere_drag` → C_D = 24/Re + 2.6(Re/5)/(1 + (Re/5)^1.52) + 0.411(Re/2.63e5)^−7.94/(1 + (Re/2.63e5)^−8.00) + 0.25(Re/1e6)/(1 + Re/1e6); valid to Re = 1e6; plateau ≈ 0.14; 24/Re for Re < 2 | explicit sphere-drag correlation (our Fig. 4.21 curve) | F. A. Morrison, "Data Correlation for Drag Coefficient for Sphere", Michigan Tech (2016), Eq. (1); *An Introduction to Fluid Mechanics*, CUP 2013, Fig. 8.13 | https://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdf | PDF downloaded, text extracted with PyMuPDF, formula and remarks read | 2026-09-23 |
| `wgs84` → a = 6 378 137.0 m, 1/f = 298.257223563, b ≈ 6 356 752.314245 m, ω = 72.92115e-6 rad/s | reference ellipsoid and Earth rotation rate (Coriolis, centrifugal, oblateness 2(a − b)) | NGA WGS 84 via Wikipedia, "World Geodetic System" | https://en.wikipedia.org/wiki/World_Geodetic_System | page fetched, defining-parameter table read | 2026-09-23 |
| `vena_contracta_sharp_orifice` → C_c = 0.611 | contraction coefficient of a sharp orifice | Wikipedia, "Vena contracta" ("The typical value may be taken as 0.611 for a sharp orifice") | https://en.wikipedia.org/wiki/Vena_contracta | page fetched, sentence read | 2026-09-23 |
| `belanger_jump` → h₂/h₁ = (√(1 + 8Fr₁²) − 1)/2, Fr₁ = v₁/√(gh₁); a moving surge is a stationary jump in its own frame | form cross-check of Example 4.3's bore speed | Wikipedia, "Hydraulic jump" | https://en.wikipedia.org/wiki/Hydraulic_jump | page fetched, equation and the moving-surge sentence read | 2026-09-23 |
| `capillary_length_water_20C` → λ_c = √(γ/Δρg); 2.71 mm (water–air, 20 °C) | capillary length (the δ of Example 4.7) | Wikipedia, "Capillary length" | https://en.wikipedia.org/wiki/Capillary_length | page fetched; the page does not state its σ or g, so our comparison uses IAPWS σ (ch01) and a 3-digit rounding band | 2026-09-23 |
| `tsiolkovsky` → Δv = v_e ln(m₀/m_f) | Example 4.4 with g = F_S = 0 (form) | Wikipedia, "Tsiolkovsky rocket equation" | https://en.wikipedia.org/wiki/Tsiolkovsky_rocket_equation | page fetched, equation read | 2026-09-23 |
| `prandtl` → air 0.70–0.73 (250–1000 K); water 10.6 (280 K), 5.9 (300 K); monatomic 2/3; Eucken 4γ/(9γ − 5) | Prandtl numbers (4.116) | Wikipedia, "Prandtl number" | https://en.wikipedia.org/wiki/Prandtl_number | page fetched, values read | 2026-09-23 |
| `added_mass_sphere` → (2/3)πr³ρ = ½ displaced mass | added mass of a sphere (C12's accelerating sphere) — form cross-check | Wikipedia, "Added mass" | https://en.wikipedia.org/wiki/Added_mass | page fetched, formula read | 2026-09-23 |
| `taylor_green` → u = U₀ sin kx cos ky F, v = −U₀ cos kx sin ky F, F = e^{−2νk²t}, p = (ρU₀²/4)(cos 2kx + cos 2ky)F² | exact NS solution used as a test field — form cross-check | Wikipedia, "Taylor–Green vortex" | https://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex | page fetched, formulas read | 2026-09-23 |

Reused from `reference/ch01/` (verified there 2026-09-12): USSA-1976 sea-level speed of sound 340.29 m/s
(`ussa1976_table1.csv`, NASA-TM-X-74335), IAPWS R1-76(2014) surface tension of water (`iapws_sigma.csv`).

## Labels
Only published **numbers** are labelled V5 in the report (Morrison C_D values of its own formula, WGS-84 Ω²a and
2(a − b), C_c = 0.611, λ_c = 2.71 mm, Pr ranges and 2/3, USSA c). The Wikipedia **forms** (Bélanger, Tsiolkovsky,
added mass, Taylor–Green) are recorded as V1 "form cross-checks (not V5)", as in ch03.

## Rejected or not used as V5
| candidate | why |
|---|---|
| Wikipedia "Volume viscosity" | does not state μ_v = λ + ⅔μ in the book's form (analysis §8); (4.32)–(4.37) are checked by V2 only |
| Engineering ToolBox / nuclear-power.com water Pr ≈ 7.0 at 20 °C | search results only, not fetched; the water check uses Wikipedia's 300 K value instead |
| Any number from Examples 4.1–4.8 or the chapter's text | book-derived; private `tests/book_values_ch04.json` (V6), never in `reference/` |
