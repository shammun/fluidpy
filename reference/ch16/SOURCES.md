# reference/ch16 — sources

Everything in `benchmarks.json` is written by `make_refs.py`. Each entry was read first-hand on the date given. Nothing in
this folder comes from the textbook. The one primary source is Secomb & Pries (2013), whose equations and sentence were
read from the PDF's extracted text. The Wikipedia pages are **secondary** sources and are stored as **forms** (a formula
as the page writes it): the tests evaluate our functions against the form with our own inputs. The Imperial College
lecture notes give an independent derivation of the Dean equations from the Navier–Stokes equations; the tests solve
that form numerically (a Chebyshev–Fourier Newton solver in `tests/test_ch16.py`) and use no number from the notes.

| value / form | what it is | source | URL | how obtained | verified |
|---|---|---|---|---|---|
| `pries_secomb_pries_2013`: Eqs. (3)–(5) (coefficients 220, 1.3, 3.2, 2.44, 0.06, 0.645, 0.8, 0.075, 10⁻¹¹, 12, H = 0.45) and the sentence "the trend of apparent viscosity to decrease with decreasing diameter continues down to about 7 µm … critical minimum diameter, which is about 2.7 µm" | in vitro relative apparent viscosity of blood in glass tubes against diameter and discharge hematocrit (high shear) | T. W. Secomb & A. R. Pries, "Blood viscosity in microvessels: experiment and theory", C. R. Physique 14, 470–478 (2013), §2.1, doi:10.1016/j.crhy.2013.04.002 | https://comptes-rendus.academie-sciences.fr/physique/item/10.1016/j.crhy.2013.04.002.pdf | PDF fetched; text extracted with pymupdf; equations and sentence read and typed | 2026-10-14 |
| `moens_korteweg_wikipedia`: PWV = √(E_inc h/(2rρ)), r the radius | Moens–Korteweg pulse-wave speed (form) | Wikipedia, "Moens–Korteweg equation" (secondary) | https://en.wikipedia.org/wiki/Moens%E2%80%93Korteweg_equation | page fetched and read | 2026-10-14 |
| `windkessel_wikipedia`: I(t) = P/R + C dP/dt; P(t) = P(t_d)e^{−(t−t_d)/(RC)} | two-element Windkessel and its diastolic decay (forms) | Wikipedia, "Windkessel effect" (secondary) | https://en.wikipedia.org/wiki/Windkessel_effect | page fetched and read | 2026-10-14 |
| `murray_wikipedia`: P = (8μ/π)lQ²/r⁴ + πλlr²; Q ∝ r³; r³ = r₁³ + r₂³ (no branching angle on the page) | Murray's cost and cube law (forms) | Wikipedia, "Murray's law" (secondary; cites C. D. Murray, PNAS 12, 207 (1926)) | https://en.wikipedia.org/wiki/Murray%27s_law | page fetched and read | 2026-10-14 |
| `pulsatile_flow_wikipedia`: u, Q and τ_w of Womersley's solution with Λ_n = α n^{1/2} i^{3/2}, ∂p/∂x = Σ P′_n e^{inωt} | Womersley's oscillatory pipe flow (forms; their τ_w = μ∂u/∂r, opposite in sign to ours) | Wikipedia, "Pulsatile flow" (secondary; cites J. R. Womersley, J. Physiol. 127, 553 (1955)) | https://en.wikipedia.org/wiki/Pulsatile_flow | page fetched and read | 2026-10-14 |
| `entrance_length_wikipedia`: L_h = 0.0575 Re_D D (laminar), L_h = 1.359 D Re_D^{1/4} (turbulent) | hydrodynamic entrance length (forms) | Wikipedia, "Entrance length (fluid dynamics)" (secondary; Bergman & Incropera; Çengel) | https://en.wikipedia.org/wiki/Entrance_length_(fluid_dynamics) | page fetched and read | 2026-10-14 |
| `dean_equations_imperial_m4a33`: K(ψ_z v_x − ψ_x v_z) = 1 + ∇²v; K(ψ_z Ω_x − ψ_x Ω_z) = ∇²Ω − 2Kvv_z; Ω = −∇²ψ; ∇ψ = 0, v = 0 on the wall | the Dean equations derived from the Navier–Stokes equations for a slightly curved pipe (form; no numbers used) | Imperial College London, Department of Mathematics, course M4A33 lecture notes "Flow in curved pipes: the Dean equations" (author's page ~ajm8) | https://www.ma.ic.ac.uk/~ajm8/M4A33/Dean.pdf | PDF fetched; text extracted with pymupdf; eqs. (7)–(10) read; the vorticity equation re-derived by hand as the curl of the momentum equations | 2026-10-14 |

## Read and deliberately not used as a benchmark
- Wikipedia, "Womersley number" (read 2026-10-14): its human-vessel Womersley numbers (aorta 13.83 for a 25 mm vessel,
  artery 2.21, capillary 4.43 × 10⁻³) do not state the heart rate or which length (radius or diameter) they use, and the
  canine list assumes 2 Hz; only the order of magnitude is consistent with our `vessel_scales_table`, so no number is
  stored.
- Dean (1928) itself was not reached; the small-K flux coefficients usually quoted from it could not be read first-hand.
  The chapter's coefficients are therefore checked by two independent routes of our own (sympy series of the equations;
  spectral Newton solution of the Navier–Stokes form above), not against a cited number.
- Womersley (1955, 1957) and Morgan & Kiely (1954) for the viscous elastic-tube wave speed: not reached as text; the
  closed form (16.136) is checked against an independent 4 × 4 determinant of the un-eliminated equations instead.

## Not benchmarks
Our own numbers (the Starling-resistor choke 1.998989 × 10⁻⁵ m³/s at p₂ − p_e = −946.13 Pa for the design's vein, the
symmetric Murray angle 37.467°, c/c₀ at α = 8, the Dean K⁴ coefficient 0.0119312) are computed, not cited. The book's own
numbers are in the git-ignored `tests/book_values_ch16.json` and are used only by `test_book_V6_worked_numbers`.
