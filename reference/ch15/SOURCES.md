# reference/ch15 — sources

Everything in `benchmarks.json` is written by `make_refs.py`. Each entry was read first-hand on the date given. Nothing in
this folder comes from the textbook. The primary source is **NACA Report 1135** — a U.S. government report, in the public
domain — whose table rows were read from rendered page images and typed by hand, and whose equations are stored as
**forms** (a formula as the report writes it, under the report's own equation number): the tests evaluate our functions
against the form with our own inputs. The four Wikipedia pages are **secondary** sources, used for a second table and
for the friction and heating forms that the report does not contain. No cached model run is committed for this chapter.

| value / form | what it is | source | URL | how obtained | verified |
|---|---|---|---|---|---|
| `naca1135_table1_subsonic`: rows M = 0.50, 0.55, 0.60 (p/p_t, ρ/ρ_t, T/T_t, β, q/p_t, A/A*, V/a*) | isentropic flow, γ = 7/5, subsonic | Ames Research Staff, *Equations, Tables, and Charts for Compressible Flow*, NACA Report 1135 (1953), Table I, report p. 633 | https://www.nasa.gov/wp-content/uploads/2023/03/equations-tables-charts-compressibleflow-report-1135.pdf | PDF downloaded; page rendered at 4× and read as an image; typed by hand; then reproduced by our functions (20 of 21 entries to half a unit of the last digit; A/A* at M = 0.55 is printed 1.2550 where the formula gives 1.254948 — one unit, 2.4e-6 from the rounding boundary) | 2026-10-13 |
| `naca1135_table2_supersonic`: rows M = 1.60, 1.65, 1.70, 2.00, 2.05, 2.10 (the seven isentropic columns, ν, μ in degrees, M₂, p₂/p₁, ρ₂/ρ₁, T₂/T₁, p_t2/p_t1, p₁/p_t2) | isentropic flow, Prandtl–Meyer angle, Mach angle and normal-shock ratios, γ = 7/5 | same report, Table II, report p. 634 | same PDF | page rendered at 3.2× and read as an image; typed by hand; all 90 entries reproduced to half a unit of the last printed digit | 2026-10-13 |
| `naca1135_forms`: equations (105), (106), (111), (124), (127), (132), (137), (138), (139b), (142), (145), (146), (147), (167), (168), (169), (170), (171c); the weak-shock series in ξ = p₂/p₁; the series for the pressure in a small expansion; "attached for a wedge of semivertex angle less than sin⁻¹(1/γ) = 45.6° for γ = 7/5" | normal- and oblique-shock relations, Prandtl's relation for an oblique shock, the shock angle for sonic flow behind the shock and for maximum deflection, their near-sonic approximations (coefficients 0.2946 and 0.3208), weak-shock entropy series (coefficients 35/216 and 245/864 for γ = 7/5), the Prandtl–Meyer function | same report, pp. 620–626 | same PDF | pages rendered and read as images; equations transcribed by hand (θ in the report is our σ) | 2026-10-13 |
| `normal_shock_tables`: rows M₁ = 1.5, 2.0, 2.5, 3.0 (M₂, p₂/p₁, ρ₂/ρ₁, T₂/T₁, p₀₂/p₀₁) | normal-shock table, γ = 1.4 | Wikipedia, "Normal shock tables" (secondary) | https://en.wikipedia.org/wiki/Normal_shock_tables | page fetched and the four rows extracted; the M₁ = 2 row agrees with NACA Report 1135 | 2026-10-13 |
| `choked_flow`: p*/p₀ = (2/(γ + 1))^{γ/(γ−1)}, "p* = 0.528 p₀"; table of minimum upstream/downstream pressure ratios for eight gases; ṁ = C_d A √(γ ρ₀ P₀ (2/(γ + 1))^{(γ+1)/(γ−1)}) | critical pressure ratio and choked mass flow | Wikipedia, "Choked flow" (secondary) | https://en.wikipedia.org/wiki/Choked_flow | page fetched and read | 2026-10-13 |
| `fanno_flow`: dM²/M², 4fL*/D_h, the starred ratios, Δs/c_p | adiabatic constant-area flow with friction (ratios to the sonic state of the same process) | Wikipedia, "Fanno flow" (secondary) | https://en.wikipedia.org/wiki/Fanno_flow | page fetched and the formulas read | 2026-10-13 |
| `rayleigh_flow`: dM²/M², the starred ratios, "the maximum value of H occurs at M = 0.845" | frictionless constant-area flow with heat exchange | Wikipedia, "Rayleigh flow" (secondary) | https://en.wikipedia.org/wiki/Rayleigh_flow | page fetched and read (0.845 is stated without the formula 1/√γ; γ = 1.4 is implied) | 2026-10-13 |

## Read and deliberately not used as a benchmark
- NACA Report 1135, charts 2–4 (oblique-shock deflection, pressure and Mach number against shock angle): they would have
  to be digitised from a plot; the report's equations (132), (138), (167), (168) give the same information exactly, so the
  charts are not used.
- The analyst's row M = 1.90 of Table II (in `analysis/ch15.md` §8) was not re-read by the verifier and is therefore not
  stored; six other rows of the same page were.
- NASA Glenn, "Normal Shock Wave Equations" and "Isentropic Flow Equations" (read by the analyst): same forms as the
  report's; not stored a second time. Their stations are numbered 0 and 1 (trap T6).
- Course notes on linearised supersonic flow found by search (Ackeret's C_p = 2θ/√(M² − 1)): not read first-hand; the
  thin-airfoil tests rest on the limit of the exact shock–expansion result and on the report's expansion series instead.
- Taylor's blast-wave constant is reused from `reference/ch01/benchmarks.json` (see `reference/ch01/SOURCES.md`); the crude
  constant of this chapter is 9.4 % below it and is labelled qualitative.

## Not benchmarks
Our own numbers (δ_max = 22.97° at M₁ = 2, ν_max = 130.45°, the flux factor 0.6847, the nozzle thresholds) are computed,
not cited. The book's own numbers are in the git-ignored `tests/book_values_ch15.json` and are used only by the
`test_book_V6_*` tests.
