# Chapter 12 — reference data and where every number comes from

Nothing in this folder comes from the textbook. Book-quoted numbers used as test evidence live only in the git-ignored
`tests/book_values_ch12.json`. Access date of every source below: **2026-10-07** (verifier), unless a row says otherwise.

## Files

| file | what it is | source | URL | how obtained | licence / terms | verified |
|---|---|---|---|---|---|---|
| `lee_moser_2015_channel_mean.csv` | mean velocity of turbulent channel flow from direct numerical simulation, wall units: y/δ, y⁺, U⁺, dU⁺/dy⁺ at Re_τ = 182.088, 543.496, 1000.512, 1994.756, 5185.897 (281 rows: 43 + 47 + 48 + 51 + 92) | M. Lee & R. D. Moser, "Direct numerical simulation of turbulent channel flow up to Re_τ ≈ 5200", J. Fluid Mech. 774, 395–415 (2015), doi:10.1017/jfm.2015.268, arXiv:1410.7809 | https://turbulence.oden.utexas.edu/channel2015/data/LM_Channel_{0180,0550,1000,2000,5200}_mean_prof.dat | downloaded by `make_refs.py` (HTTP 200 for all five files), then **thinned, not interpolated**: the wall row, ~72 rows log-spaced in y⁺, and for Re_τ = 5200 every third row in 300 ≤ y⁺ ≤ 800; 8 significant digits | no licence statement in the files; their header asks that the paper above be cited when the data are used. A small subset is redistributed here with that citation; the full files stay in the git-ignored `data/online/ch12/` | 2026-10-07 |
| `benchmarks.json` | the scalar published values of the table below, with their citations, plus the header numbers (Re_τ, u_τ, ν, row counts) of the five DNS files | see below | see below | typed from the sources below by `make_refs.py` | cited values | 2026-10-07 |
| `explainer_tables.json` | **ours**: the mixing-length channel (Re_τ = 180, 550, 1000, 5200; κ = 0.41, A⁺ = 26) and the van Driest wall profile, for the explainers | computed by `fluidpy.ch12_turbulence.explainer_tables` (`scripts/ch12_tables.py`) | — | deterministic; `tests/test_ch12.py` recomputes rows and the whole file | ours | 2026-10-07 |

## Scalar values (`benchmarks.json`)

| value | what it is | source | URL | how obtained | verified |
|---|---|---|---|---|---|
| κ = 0.384 ± 0.004 at Re_τ = 5186 | von Kármán constant of the logarithmic region of the channel DNS | Lee & Moser (2015), abstract | https://arxiv.org/abs/1410.7809 | page fetched; the abstract states the value and its uncertainty; no additive constant is given there | 2026-10-07 |
| C₁ = 0.53, standard deviation 0.055, 95 % interval of the mean ± 0.03 | one-dimensional (longitudinal, one-sided) Kolmogorov constant | K. R. Sreenivasan, Phys. Fluids 7, 2778–2784 (1995), as summarised by the ATOMIX wiki "Spectra in the inertial subrange" | https://atomix.app.uib.no/Spectra_in_the_inertial_subrange | page fetched (secondary summary of the primary paper) | 2026-10-07 |
| C₁ = 18C/55 ≈ 27/55 (so C = 1.5) | isotropic relation between the one- and three-dimensional constants | same ATOMIX page; the factor 18/55 is also re-derived by sympy in `tests/test_ch12.py` | same | page fetched; factor proved symbolically | 2026-10-07 |
| C_μ = 0.09, C_ε1 = 1.44, C_ε2 = 1.92, σ_k = 1.0, σ_ε = 1.3 | constants of the standard k–ε model | Launder & Sharma (1974); Launder & Spalding (1974) | https://www.openfoam.com/documentation/guides/latest/api/classFoam_1_1RASModels_1_1LaunderSharmaKE.html · https://www.simscale.com/docs/content/simulation/model/turbulenceModel/kEpsilon.html | web search summary quoting both pages (cfd-online's wiki page returned HTTP 403 to the fetch tool) | 2026-10-07 |
| 1/√f = 2.0 log₁₀(Re√f) − 0.8 | Prandtl's friction law for smooth pipes | Prandtl (1935), as quoted by McKeon, Zagarola & Smits, J. Fluid Mech. 538, 429 (2005) | https://authors.library.caltech.edu/records/yn6nj-pkr96 | web search summary | 2026-10-07 |
| φ_m = (1 − 16 z/L)^{−1/4} for z/L < 0 | unstable flux–profile (Businger–Dyer) relation | AMS Glossary of Meteorology, "Businger–Dyer relationship" | https://glossary.ametsoc.org/wiki/businger-dyer-relationship/ | web search summary (the page itself returned HTTP 403 to the fetch tool). The same summary gives the stable form as 1 + 5 z/L (4.7 is the coefficient of Businger et al. 1971). The integrated ψ_m is **not** taken from a citation: the tests prove it is ∫₀^ζ (1 − φ_m)/ζ′ dζ′ by quadrature | 2026-10-07 |
| p₀ = 2, c_L = 6.78 | energy-range factor of Pope's model spectrum | S. B. Pope, *Turbulent Flows* (2000), §6.5, as quoted in arXiv:1705.04917 and the ATOMIX wiki "Pope Model Shear Spectrum" | https://arxiv.org/abs/1705.04917 · https://atomix.app.uib.no/Pope_Model_Shear_Spectrum | web search summary; **and** proved numerically in the tests: c_L = 6.78 is the value for which ∫E dK = (ε̄L)^{2/3}, i.e. L = ē^{3/2}/ε̄, at high Reynolds number (root found: 6.779 at L/η = 5.6 × 10⁶) | 2026-10-07 |

## Not used as benchmarks (and why)

- Constants of the free-shear-flow table (jets, wakes, plumes, shear layer): no public table was confirmed; `FREE_SHEAR_CONSTANTS`
  is empty, examples pass labelled illustrative constants, the book's set stays private (V6 test only).
- Per-flow (κ, B) pairs of Nagib & Chauhan (2008): only the pipe value was confirmed from an open page by the analyst; the
  correlation κB = 1.6[exp(0.1663B) − 1] is tested as an identity, the book's pairs only in the private V6 test.
- Model results (mixing-length channel, Spalding profile, k–ε channel) are compared with the DNS under the label
  "approximate"/"qualitative": the DNS is the benchmark, the models are not.
