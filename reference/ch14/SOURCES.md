# reference/ch14 — sources

Everything in `benchmarks.json` is written by `make_refs.py`. Each entry was read first-hand on the date given. Nothing in
this folder comes from the textbook. Most entries are **forms** (a formula as the source writes it): the tests evaluate
our functions against the form with our own inputs. No cached model run is committed for this chapter.

| value / form | what it is | source | URL | how obtained | verified |
|---|---|---|---|---|---|
| `thin_airfoil_lift_slope`: c_l = 2πα; cambered c_l = c_l0 + 2πα | lift-curve slope of a thin section, per radian | Wikipedia, "Airfoil", section "Thin airfoil theory" (secondary; its own references, Abbott & von Doenhoff and Clancy, were not read) | https://en.wikipedia.org/wiki/Airfoil | read on the page | 2026-10-12 |
| `thin_airfoil_centre_of_pressure_symmetric`: 0.25 | centre of pressure of a thin symmetric section, "exactly one quarter of the chord behind the leading edge" | same page | https://en.wikipedia.org/wiki/Airfoil | read on the page | 2026-10-12 |
| `joukowsky_kutta_circulation`: Γ = 4πV∞R sin(α + sin⁻¹(μ_y/R)), R = √((1 − μ_x)² + μ_y²) | circulation fixed by the Kutta condition for z = ζ + 1/ζ and a circle centred at μ through ζ = 1 | Wikipedia, "Joukowsky transform" (secondary) | https://en.wikipedia.org/wiki/Joukowsky_transform | read on the page | 2026-10-12 |
| `induced_drag_coefficient`: C_D,i = C_L²/(π AR e); e = 1.0 for an elliptic distribution; AR = span²/area | induced drag coefficient | NASA Glenn Research Center, Beginner's Guide to Aeronautics, "Induced Drag Coefficient" (T. Benson); second source Wikipedia, "Lift-induced drag" | https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/induced.html · https://en.wikipedia.org/wiki/Lift-induced_drag | read on both pages | 2026-10-12 |
| `induced_drag_force`: D_i = L²/(½ρV²πb²) | induced drag from lift and span | Wikipedia, "Lift-induced drag" (secondary) | https://en.wikipedia.org/wiki/Lift-induced_drag | read on the page | 2026-10-12 |
| `induced_share_wide_body_cruise`: 0.37 | induced share of the cruise drag of a typical twin-engine wide-body ("approximately 37 %"; skin friction "almost 48 %") — evidence for trap T9 only, **not a test of our code** | Wikipedia, "Lift-induced drag" (secondary) | https://en.wikipedia.org/wiki/Lift-induced_drag | read on the page | 2026-10-12 |
| `wagner_jones_fit`: φ(s) = 1 − 0.165 e^{−0.0455 s} − 0.335 e^{−0.3 s}, s in semi-chords; φ(0) = 0.5; absolute error "less than 1 %" | R. T. Jones' approximation of Wagner's function | S. T. M. Dawson & S. L. Brunton, "Improved approximations to the Wagner function using sparse identification of nonlinear dynamics", arXiv:2104.15122v2 (2021), Eq. (1), Table 1 and the text under Fig. 2; second source BYU FLOW Lab, Aeroelasticity.jl documentation, "Wagner's Function" | https://arxiv.org/abs/2104.15122 · https://flow.byu.edu/Aeroelasticity.jl/dev/models/aerodynamics/wagner | PDF text extracted and read; BYU page fetched and read. **The primary report (R. T. Jones, NACA TN 667, 1938, as Dawson & Brunton cite it) was not read** | 2026-10-12 |
| `wagner_garrick_fit`: φ(s) = 1 − 2/(4 + s); absolute error "2 %" | Garrick's approximation of Wagner's function (an independent cross-check of the Jones fit) | Dawson & Brunton, arXiv:2104.15122v2, Eq. (2) and the text under Fig. 2 | https://arxiv.org/abs/2104.15122 | PDF text extracted and read | 2026-10-12 |

## Read and deliberately not used as a benchmark
- NASA Glenn, "Downwash Effects on Lift" (https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/downwash-effects-on-lift/,
  read 2026-10-12) writes Cl = Cl0/(1 + Cl0/(π AR)) with Cl0 a lift **coefficient**. Lifting-line theory for an untwisted
  elliptic wing gives that factor with the section lift **slope** K in the denominator, dC_L/dα = K/(1 + K/(πΛ)); the two
  agree only if "Cl0" is read as the slope. The tests verify our form analytically and with the solver instead.
- Wikipedia, "Lifting-line theory": its lift-slope factor carries the page's own "citation needed" flag (noted by the
  analyst); not used.
- A table of Glauert's induced-drag factor for rectangular wings was not found in a readable public source; the
  rectangular-wing span efficiency is ours by convergence only.
- Measured section data (Abbott & von Doenhoff, NACA reports) were not fetched; the stall curves stay qualitative.

## Not benchmarks
Standard public constants (2π, the low-Mach limit 0.3, 3√3/4) are computed, not cited. The book's own numbers are in the
git-ignored `tests/book_values_ch14.json` and are used only by the `test_book_V6_*` tests.
