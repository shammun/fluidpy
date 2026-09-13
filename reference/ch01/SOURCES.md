# reference/ch01 — sources of the public benchmark data

All files in this folder are written by `make_refs.py`. Every value was read from the public source below by the
math-verifier on the date given (PDFs downloaded and their text extracted with pymupdf, or the scanned page rendered
and read; HTML tables parsed). No value comes from the textbook. Tolerances used by `tests/test_ch01.py` are stated in
`reports/ch01_verification.md`.

| value / file | what it is | source | DOI / URL | how obtained | verified |
|---|---|---|---|---|---|
| `constants.json` k_B = 1.380649e-23 J/K | Boltzmann constant (exact) | CODATA (NIST CUU) | https://physics.nist.gov/cgi-bin/cuu/Value?k | page fetched | 2026-09-12 |
| `constants.json` N_A = 6.02214076e23 1/mol | Avogadro constant (exact) | CODATA (NIST CUU) | https://physics.nist.gov/cgi-bin/cuu/Value?na | page fetched | 2026-09-12 |
| `constants.json` R = 8.314462618… J/(mol K) | molar gas constant (exact, 2022 CODATA) | CODATA (NIST CUU) | https://physics.nist.gov/cgi-bin/cuu/Value?r | page fetched | 2026-09-12 |
| `ussa1976_constants.json` g0, P0, T0, β = 1.458e-6, γ = 1.40 | adopted constants, Table 2 | U.S. Standard Atmosphere 1976, NASA-TM-X-74335 | https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf | scanned Table 2 (doc p. 2, PDF p. 18) rendered and read | 2026-09-12 |
| `ussa1976_constants.json` S = 110.4 K, r0 = 6356.766 km | Sutherland constant, effective Earth radius | same, Table 2 **corrected by the NASA errata sheet** bound into the same PDF (p. 242: Table 2 prints S = 110 K and r0 ×1000 too large) | same | errata text extracted | 2026-09-12 |
| `ussa1976_constants.json` R* = 8314.32 J/(kmol K) | gas constant of the standard | same, doc p. 3 (PDF p. 19) (Table 2's 8.31432e-3 is the misprint named in the errata) | same | page rendered and read | 2026-09-12 |
| `ussa1976_constants.json` composition | Table 3: M_i and F_i of 10 dry-air species (gives M0 = 28.9644 kg/kmol) | same, doc p. 3 | same | page rendered and read | 2026-09-12 |
| `ussa1976_constants.json` layers | Table 4: base geopotential heights and gradients (−6.5, 0, +1.0, +2.8, 0, −2.8, −2.0 K/km′) | same, doc p. 3 | same | page rendered and read | 2026-09-12 |
| `ussa1976_table1.csv` | T, p, ρ, c, g at geometric Z = 0…50 km (SI Table 1) | Public Domain Aeronautical Software, "Tables of the U.S. Standard Atmosphere, 1976" (computed from NASA-TM-X-74335, Part 4) | https://www.pdas.com/bigtables.html | HTML downloaded, table cells parsed by script | 2026-09-12 |
| `ussa1976_table2.csv` | μ, ν, H_p, n, mean particle speed, mean free path, M (SI Table 2) | same | same | same | 2026-09-12 |
| `iapws_sigma.csv`, `benchmarks.json → iapws_sigma_equation` | σ(t) of water, Table 1 (experimental, uncertainty, calculated columns); σ = Bτ^μ(1 + bτ), B = 235.8 mN/m, b = −0.625, μ = 1.256, T_c = 647.096 K | IAPWS R1-76(2014) Revised Release on Surface Tension of Ordinary Water Substance | https://iapws.org/public/documents/CH-L9/Surf-H2O-2014.pdf | PDF text extracted (the text layer drops the minus of b; the calculated column fixes the sign) | 2026-09-12 |
| `benchmarks.json → iapws_viscosity_check` | μ(298.15 K, 998 kg/m³) = 889.735100 µPa s; μ(373.15 K, 1000 kg/m³) = 307.883622 µPa s | IAPWS R12-08, Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance, Table 4 (computer-program verification) | https://iapws.org/technical-guidance/release/viscosity.download | PDF text extracted | 2026-09-13 |
| `benchmarks.json → jennings` | mean free path of air 67.3 nm at 300 K, 1 atm (Jennings 1988); formula λ = √(π/8) μ/(0.4987445 √(ρP)); MD value 38.5 ± 1 nm | Tsalikis, Mavrantzas & Pratsinis, "A new equation for the mean free path of air", Aerosol Sci. Technol. 58(8) (2024), CC BY 4.0, quoting Jennings, J. Aerosol Sci. 19, 159 (1988) | doi:10.1080/02786826.2024.2333859 ; doi:10.3929/ethz-b-000669211 | open-access PDF text extracted | 2026-09-13 |
| `benchmarks.json → taylor_blast` | S(1.4)^−5 = 0.856 (K in E = KρD⁵/t² for a spherical blast, γ = 1.4); S(1.4) = 1.032 | J. S. Díaz, "Explosion analysis from images: Trinity and Beirut", quoting G. I. Taylor, Proc. R. Soc. Lond. A 201, 159 (1950) | https://arxiv.org/abs/2009.05674 | PDF text extracted (p. 3) | 2026-09-13 |
| `benchmarks.json → ams_lapse_rate` | dry-adiabatic lapse rate g/c_pd ≈ 9.8 °C/km (rate of decrease, i.e. the meteorology convention) | American Meteorological Society, Glossary of Meteorology, "adiabatic lapse rate" | https://glossary.ametsoc.org/wiki/Dry-adiabatic_lapse_rate | HTML downloaded with a browser user agent (WebFetch got HTTP 403) | 2026-09-13 |
| `benchmarks.json → unesco_eos80_check` | EOS-80 check values: ρ(S, t68, p) and V at S = 0, 35; t68 = 5, 25 °C; p = 0, 10000 dbar (only the p = 0 rows are used: `seawater_density_eos80` is the one-atmosphere equation) | N. P. Fofonoff & R. C. Millard Jr., *Algorithms for computation of fundamental properties of seawater*, Unesco Technical Papers in Marine Science 44 (1983), p. 19 (PDF p. 23), quoting Unesco Report 38, p. 191 | https://darchive.mblwhoilibrary.org/bitstreams/f77d18e9-e756-58eb-b042-a8870de55e3b/download (same document: https://www.jodc.go.jp/info/ioc_doc/UNESCO_tech/059832eb.pdf ; ERIC ED261897) | scanned PDF downloaded, OCR text extracted and the page image rendered and read; t is IPTS-68 (t90 = t68/1.00024) | 2026-09-13 |

## Rejected or not used as V5
| candidate | why |
|---|---|
| Engineering ToolBox c_p of air ≈ 1005 J/(kg K) | search-snippet only, not a primary source; the code derives C_p from γ = 1.40 (USSA Table 2, verified above) and R, so no independent c_p benchmark is needed |
| Sutherland constants from CFD-Online / Fluent snippets | a snippet quoted β = 1.458e-5 (a typo); replaced by the primary NASA-TM-X-74335 Table 2 + errata |
| Gill (1982) / wiki snippet of the UNESCO check value | superseded 2026-09-13 by the primary Unesco Tech. Pap. 44 table above (V5) |
| Kell (1975) water density table, maximum at 3.98 °C | primary paper not fetched; `water_density` gets V7/V1 self-consistency evidence only |
| PDAS Table 2 mean free path 6.6332e-8 m at sea level | USSA uses a hard-sphere collision-diameter formula, not the Jennings formula implemented in `mean_free_path_jennings`; not comparable |
