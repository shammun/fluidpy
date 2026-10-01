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

### Verifier re-check (math-verifier, 2026-09-30)
Each value used as V5 evidence in `tests/test_ch11.py` was re-read today in the source itself (not a search snippet) where
the source is reachable:

| value | re-read in | what the source prints | status |
|---|---|---|---|
| Re_c = 5772.22; α_c = 1.02056 ± 0.00001; c(10⁴, 1) = 0.23752649 + 0.00373967i; c_r(Re_c) = 0.26400 | Orszag (1971), J. Fluid Mech. 50, 689 — PDF copy at http://www.damtp.cam.ac.uk/user/tong/fluids/Orszag.pdf (abstract, §4 text, Tables 2 and 5) | exactly these digits | verified |
| Ra_c = 1707.76, K_c = 3.117 (rigid–rigid); 1100.65, 2.682 (rigid–free); 657.511, 2.2214 (free–free) | Nek5000 examples PDF (https://www.mcs.anl.gov/~fischer/nek5000/examples.pdf), table "Critical Rayleigh number for 3 types of boundaries" (column Ra_c, k_c, attributed to Chandrasekhar); arXiv:nlin/0302057 text (1707.76, 3.117) | exactly these digits | verified (secondary source quoting Chandrasekhar 1961) |
| Blasius Re_δ*,c = 519.2, α_c = 0.303, ω_c = 0.120 (Thomas); Gallagher et al. own n = 1 row 519.12, 0.3022, 0.1198 | Gallagher, Griffiths & Stephen, Phys. Fluids 28, 074107 (2016), Table II (PDF above) | exactly these digits | verified |
| δ = 4.669201609…, period-doubling a_n = 3, 3.4494897, 3.5440903, 3.5644073, 3.5687594 | Wikipedia "Feigenbaum constants" (logistic-map table) | exactly these digits; the accumulation point only as "≈ 3.5699…" | verified (A_∞ to 7 digits not used in tests) |
| r_H ≈ 24.74 (Pr = 10, b = 8/3), formula σ(σ + β + 3)/(σ − β − 1) | Wikipedia "Lorenz system" | formula and 24.74 | verified |
| Bickley Re_c ≈ 4.0 at k ≈ 0.2 | Tatsumi & Kakutani (1958), via the abstract of Int. J. Heat Fluid Flow paper (sciencedirect S0142727X16306518) | "critical Reynolds number … 4.0 at a non-dimensional wavenumber 0.2" | verified as an *approximate* value: our Re_c = 4.017 agrees to 0.4 %, our k_c = 0.1728 differs by 14 % (confirmed by an independent half-domain solve in the tests; the 1958 analysis was approximate) |
| tanh most amplified k ≈ 0.4446, k c_i ≈ 0.1897 | Michalke (1964) — primary not reachable; 0.4446 appears in secondary snippets | — | approximate (tested at 1e-3, our 0.44492, 0.18970) |
| odd mode Ra = 17610.39 at K = 5.365 | Chandrasekhar (1961) p. 39 — not reachable online today | — | not verified online: the tests use two independent routes (Chebyshev parity filter and the odd sin/sinh determinant agree to 1e-8) and the private V6 book value instead |

### Verifier re-check of our computed tables (math-verifier, loop 2, 2026-10-01; replaces the loop-1 notes)
- `tg_growth_map.csv` and the `tg_map` key of `explainer_tables.json` (regenerated in loop 2 with the box rule
  y_max = max(30, 12/k) and the tan-map scale rule s = min(0.5, max(0.035, 0.025/k)), N = 100): identical to each other
  entry for entry; a full live recomputation (620 points) rounds to the CSV in all 620 entries (largest relative rounding
  difference 4.1e-4). 335 positive entries = exactly the 335 grid points with J < k(1 − k); none on or above the exact
  neutral curve. (Loop 1 had 299: the 36 missing ones were unstable points at k ≥ 0.65 reported as 0 — closed.)
- Independent check of the values: shooting on the unbounded layer (no box, no Chebyshev grid) agrees with the default
  `tg_growth` to ≤ 1.7e-6 absolute at 35 points, including the points next to the neutral curve at k = 0.65 … 0.97.
- Stated limitation, measured: off the grid a weak mode (k c_i < 0.004) within 0.0002 (k = 0.05) … 0.006 (k = 0.95) of the
  neutral curve is still reported as 0 (CSV header, docstrings) — confirmed by bisection in 13 columns.
- `rayleigh_spectra.json` / the `rayleigh` key (regenerated in loop 2 with the box max(profile box, 12/k)): identical to
  each other; a full live recomputation rounds to the file in every entry; the k = 0.1 entries now equal the
  box-converged values (shear layer c_i 0.8364, jet 0.0928 + 0.215i; the fixed boxes gave 0.836 and 0.09273 + 0.2149i).
- Unchanged against commit `8ae4f15`: the keys `benard`, `taylor`, `lorenz_sweep`, `note` of `explainer_tables.json`,
  `benchmarks.json`, `critical_points.json`, `os_modes.json`.

### Verifier addition after the review (math-verifier, 2026-10-01) — Taylor numbers away from μ = 1 (review Should-fix 7)
| value/table | what it is | source | DOI/URL | how obtained | verified |
|---|---|---|---|---|---|
| Ta_c = 1707.76 [1 − 0.00761 ((1 − μ)/(1 + μ))²] as μ → 1, with Ta = −2AΩ₁d⁴(1 + μ)/ν² (= ½(1 + μ) × the Ta = −4AΩ₁d⁴/ν² of (11.52)) | narrow-gap critical Taylor number for co-rotating cylinders, leading order in 1 − μ | Wikipedia "Taylor–Couette flow", section "Taylor's criterion" (revision 1373953105, 2026-09-08) | https://en.wikipedia.org/wiki/Taylor%E2%80%93Couette_flow | wikitext read with `action=raw` (formula and the definition of Ta typed from it) | 2026-10-01 |

- Tertiary source and an asymptotic formula: used at 1e-5 for μ ≥ 0.5 and 1e-3 for 0 ≤ μ < 0.5 (the general rule is 1 %).
  Ours × ½(1 + μ) against the formula: μ = 1: 1.0e-6 (the source rounds 1707.762); 0.75: 1.2e-6; 0.5: 3.0e-6; 0.25: 1.7e-5;
  0: 1.1e-4 (1694.950 vs 1694.764 — the next order in 1 − μ). The coefficient measured from our Ta_c at μ = 0.9 and 0.8 is
  0.007603 and 0.007602 (the source prints 0.00761: 1e-3 relative in the coefficient, < 1e-7 in Ta_c there).
- **Not found: a counter-rotating value (μ < 0, e.g. Chandrasekhar's narrow-gap table at μ = −1).** Tried today: Chandrasekhar
  (1961) is not readable online (archive.org copy is borrow-only); arXiv:2601.14806 (only μ = 1: T_c ≈ 1708, α_c ≈ 3.117);
  arXiv:2608.10951 (T_c(μ) only as a figure; states μ_c ≈ −0.8 below which the first instability is non-axisymmetric);
  the Royal Society narrow-gap review (rsta 381, 20220134) and MDPI Fluids 6, 306 return HTTP 403; Scholarpedia's
  certificate has expired; arXiv:physics/0502069 has no such table. No number for μ < 0 is cited anywhere in the tests;
  our μ = −0.5 and μ = −1 values stay supported by the Galerkin route and N-convergence only.

### Verifier re-check of the tables regenerated in post-review loop 2 (math-verifier, 2026-10-01)
No published value changed: `benchmarks.json` is byte-identical to what `make_refs.py --no-tables` writes, and the Blasius
benchmarks (Thomas 519.2 / 0.303 / 0.120; Gallagher et al. n = 1 row 519.12; Jordinson 520) are untouched. What changed is
**ours**, and each file was compared with a live recomputation and with a route that shares no code with the solver that
wrote it. This list supersedes the line "Unchanged against commit `8ae4f15` … `critical_points.json`, `os_modes.json`" above.

| file | what changed | live recomputation | independent route |
|---|---|---|---|
| `rayleigh_spectra.json`, `explainer_tables.json → rayleigh` | Rayleigh eigenvalues now from `core.stability.rayleigh_eigs_contour` (collocation on a complex path). Six entries moved: jet k = 1.8 → 0.6324 + 0.02432i and k = 1.9 → 0.6497 + 0.01163i (were 0 — false zeros), jet k = 1.6, 1.7, tanh layer k = 0.9, sin k = 0.8 in the 4th figure. Only the sub-keys `jet`, `shear_layer`, `sin`, `shear_layer_walls` of the explainer copy differ from HEAD | all 56 non-zero entries round to the file; the explainer copy equals the file | real-axis shooting written in the test file: ≤ 1.0e-10 on all 56; `rayleigh_shoot`: ≤ 1.05e-10; the eigenvalues satisfy the real-axis integral identities of (11.83)–(11.84) to 7e-13; exact neutral wavenumbers (2, 1, 1, √3/2) approached linearly |
| `os_neutral_bickley.csv` | wavelength-scaled box max(40, 12/k); k scanned to 2.0. New picture: one unstable band whose lower edge leaves k ≥ 0.02 above Re = 7.2 (NaN = unstable down to k = 0.02, edge not resolved), a stable gap from Re ≈ 17.5, upper branch finite to Re = 1000 (→ 2) | 7 of 40 rows recomputed: identical at 4 s.f., same NaN pattern | compound-matrix shooting (exact free-stream start, no box): neutral k agree to ≤ 2e-7 at 7 neutral points; eigenvalues to ≤ 7e-8 at 14 (Re, k) points (3.5e-7 at Re = 1000) |
| `os_neutral_bickley_longwave.csv` (**new**) | the long-wave unstable band below the stable gap (Re = 19.3, 22.2, 25.6: upper edge 0.04472, 0.03091, 0.0222; lower edge below k = 0.02, unknown), with the `same_mode` flag | the three finite rows recomputed: identical | compound matrix: k_long_upper to 8e-8; sign of c_i on both sides of the gap at Re = 17.6 |
| `os_grid_bickley.csv`, `os_grid_tanh.csv`, `os_grid_blasius.csv` | least-damped *discrete* mode in the wavelength-scaled box (600 rows each, none without a mode; 93 / 237 / 425 rows lie below the edge c_i = −k/Re of the continuous spectrum) | 13 rows recomputed: c, P, Λ, E equal at 4 s.f. | compound matrix: ≤ 6e-9 (Bickley, tanh), ≤ 1.1e-7 (Blasius), 1.5e-6 at Re = 6000 (the tables' N = 80; N = 120 agrees to 1e-10) |
| `os_neutral_blasius.csv` | box max(20, 12/k) δ*: lower branch up to +9.0 % (Re = 6000: 0.07248 → 0.07902); upper branch moved by ≤ 3.8e-4 | 3 rows recomputed: identical | compound matrix with a freshly integrated Blasius profile: lower branch at Re = 6000 = 0.0790171 |
| `os_modes.json` | only the Blasius preset (Re = 1000, k = 0.25): c 0.349803 + 0.01208i → 0.34978 + 0.01209i; the four other presets are identical to HEAD | live `ts_mode` rounds to the file | compound matrix: 1.1e-7 |
| `critical_points.json` | only `blasius`: Re_c 519.0765 → 519.0601180747, k_c 0.3037752 → 0.3037711, c_r, ω_c | default live solve equals the file to 1e-9 (slow test) | compound matrix: neutral Re at our k_c = 519.0601175 (1e-9 relative), minimum over k at 0.3037710 |

Not changed (byte-identical to HEAD): `benard_neutral_curves.csv`, `taylor_critical.csv`, `tg_growth_map.csv`,
`os_grid_poiseuille.csv`, `os_neutral_poiseuille.csv`, `os_neutral_tanh.csv`; in `explainer_tables.json` the keys `note`,
`benard`, `taylor`, `tg_map`, `lorenz_sweep`.

## Our computed tables (labelled "ours" inside each file)

| file | what | written by |
|---|---|---|
| `benard_neutral_curves.csv` | marginal Ra(K): rigid–rigid, free–free, rigid–free, odd mode (K = 0.5 … 10) | `ch11.benard_neutral_table` |
| `taylor_critical.csv` | narrow-gap Ta_c(μ), k_c(μ) vs (11.54) | `ch11.write_reference_tables` |
| `tg_growth_map.csv` | kc_i(k, J) for U = tanh z, N² = J sech²z | `ch11.tg_growth_map` |
| `rayleigh_spectra.json` | leading growing Rayleigh eigenvalue c(k) per profile (complex-path collocation; 0 = c_i ≤ 1e-4) | `ch11.rayleigh_spectrum_table` |
| `os_neutral_bickley_longwave.csv` | sinuous Bickley jet: the long-wave unstable band below the stable gap (k_long_lower, k_long_upper, c_r, same_mode) | `ch11.neutral_curve_tables` (`bickley_neutral_curve`) |
| `explainer_tables.json` | the explainers' tables (≤ 4 s.f.) | `ch11.write_reference_tables` |
| `critical_points.json` | our critical points (Bénard, Taylor, Poiseuille, Blasius, Bickley, tanh, Lorenz, Feigenbaum ratios) | `ch11.write_reference_tables` |
| `os_neutral_*.csv`, `os_grid_*.csv`, `os_modes.json` | Orr–Sommerfeld neutral curves, (Re, k) grids with budgets, mode samples | `ch11.neutral_curve_tables` |
| `benchmarks.json` | the published values above | `make_refs.py` |
