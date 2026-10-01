# reference/ch11 — sources

Published benchmarks used as V5 evidence for chapter 11 (Instability), and the files we compute ourselves. Nothing in this
folder is taken from the textbook; the book's own (rounded) numbers live only in the git-ignored `tests/book_values_ch11.json`.
Regenerate everything with `.venv/Scripts/python.exe reference/ch11/make_refs.py` (heavy pieces cached in `outputs/ch11/cache`).

## Published benchmarks (`benchmarks.json`)

| value/table | what it is | source | DOI/URL | how obtained | verified |
|---|---|---|---|---|---|
| Ra_c = 1707.762, K_c = 3.117 | rigid–rigid Rayleigh–Bénard onset | Chandrasekhar (1961), *Hydrodynamic and Hydromagnetic Stability*, quoted in the Nek5000 examples and arXiv:nlin/0302057 | https://www.mcs.anl.gov/~fischer/nek5000/examples.pdf ; https://arxiv.org/pdf/nlin/0302057 | search snippets (analyst) | 2026-09-30 |
| Ra_c = 27π⁴/4, K_c = π/√2 | stress-free onset | analytic (Rayleigh 1916) | — | computed | — |
| Ra_c = 1100.65, K_c = 2.682 | rigid bottom, free top | Chandrasekhar (1961) | primary not fetched | our Chebyshev value 1100.6496 at 2.6823 agrees | 2026-09-30 |
| Ra = 17610.39, K = 5.365 | gravest odd mode, rigid–rigid | Chandrasekhar (1961) p. 39 | primary not fetched | our Chebyshev and determinant values agree (17610.394 at 5.3647) | 2026-09-30 |
| Ta_c → 1707.76 (μ → 1) | narrow-gap Taylor = rigid Bénard | arXiv:2601.14806 (small-gap Couette–Taylor) | https://arxiv.org/html/2601.14806 | read online (analyst) | 2026-09-30 |
| Re_c = 5772.22, k_c = 1.02056, c_r = 0.264 | plane Poiseuille critical point | Orszag, J. Fluid Mech. 50, 689 (1971); confirmed in Hack & Zaki, CTR Res. Briefs 2019, and Kachuma & Sobey, Oxford NA-07/21 | https://web.stanford.edu/group/ctr/ResBriefs/2019/29_Hack.pdf ; https://www.cs.ox.ac.uk/files/725/NA-07-21.pdf | read online (analyst); the Oxford report's 1.0255 is a transcription slip (our c_i there is −4.5e-6) | 2026-09-30 |
| c = 0.23752649 + 0.00373967i | plane Poiseuille, Re = 10⁴, k = 1 | Orszag (1971), as reproduced in later papers; eigentools documentation (max c_i = 3.740e-3) | https://eigentools.readthedocs.io/en/latest/notebooks/Orr%20Somerfeld%20pseudospectra.html | read online (analyst) | 2026-09-30 |
| Re_δ*,c = 519.2, αδ* = 0.303, ω = 0.120 | Blasius (parallel flow) critical point | Thomas, via Gallagher, Griffiths & Stephen, Phys. Fluids 28, 074107 (2016), Table II; Jordinson (1970): 520 | https://pure-oai.bham.ac.uk/ws/files/29212736/Griffiths_The_effect_of_non_Newtonian_viscosity_on_the_stability_of_the_Blasius_boundary_layer.pdf | read online (analyst) | 2026-09-30 |
| Re_c ≈ 4.0 at k ≈ 0.2 | Bickley jet, sinuous mode | Tatsumi & Kakutani, J. Fluid Mech. 4, 261 (1958), via later papers | https://www.sciencedirect.com/science/article/abs/pii/S0142727X16306518 | search summary; approximate | 2026-09-30 |
| k_max ≈ 0.4446, (kc_i)_max ≈ 0.1897 | tanh shear layer, most amplified inviscid mode | Michalke, J. Fluid Mech. 19, 543 (1964) | digits not verified online | label "approximate"; ours 0.44492, 0.18970 | 2026-09-30 |
| k = 1, c = 0 (φ = sech y); Bickley k = 2 / 1, c = 2/3 | inviscid neutral modes | analytic (Drazin & Reid 1981) | — | sympy/numeric residual 0 | — |
| kh = 1.278465 | piecewise-linear shear layer neutral wavenumber | root of (kh − 1)² = e^{−2kh} (Rayleigh 1880) | — | computed | — |
| J = k(1 − k) | neutral curve of U = tanh z, N² = J sech²z | exact neutral mode \|tanh z\|^{1−k} sech^k z, c = 0 (checked numerically here; Drazin & Reid 1981) | — | computed | — |
| r_H = 24.7368 (Pr = 10, b = 8/3) | Lorenz convection states lose stability | Lorenz (1963); Wikipedia "Lorenz system" | https://en.wikipedia.org/wiki/Lorenz_system | read online (analyst) | 2026-09-30 |
| δ = 4.669201609, A₂ = 1 + √6, A₃ = 3.5440903, A₄ = 3.5644073, A_∞ = 3.5699456 | Feigenbaum constant, logistic period doublings | Feigenbaum (1978); Wikipedia "Feigenbaum constants" | https://en.wikipedia.org/wiki/Feigenbaum_constants | read online (analyst) | 2026-09-30 |

## Our computed tables (labelled "ours" inside each file)

| file | what | written by |
|---|---|---|
| `benard_neutral_curves.csv` | marginal Ra(K): rigid–rigid, free–free, rigid–free, odd mode (K = 0.5 … 10) | `ch11.benard_neutral_table` |
| `taylor_critical.csv` | narrow-gap Ta_c(μ), k_c(μ) vs (11.54) | `ch11.write_reference_tables` |
| `tg_growth_map.csv` | kc_i(k, J) for U = tanh z, N² = J sech²z | `ch11.tg_growth_map` |
| `rayleigh_spectra.json` | leading unstable Rayleigh eigenvalue c(k) per profile | `ch11.rayleigh_spectrum_table` |
| `explainer_tables.json` | the explainers' tables (≤ 4 s.f.) | `ch11.write_reference_tables` |
| `critical_points.json` | our critical points (Bénard, Taylor, Poiseuille, Blasius, Bickley, tanh, Lorenz, Feigenbaum ratios) | `ch11.write_reference_tables` |
| `os_neutral_*.csv`, `os_grid_*.csv`, `os_modes.json` | Orr–Sommerfeld neutral curves, (Re, k) grids with budgets, mode samples | `ch11.neutral_curve_tables` |
| `benchmarks.json` | the published values above | `make_refs.py` |
