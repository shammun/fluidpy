# reference/ch03 — sources of the public benchmark data

`benchmarks.json` is written by `make_refs.py`. Chapter 3 (Kinematics) needs no datasets: the evidence in
`tests/test_ch03.py` is overwhelmingly Tier 1 (analytic closed forms, sympy re-derivations, convergence orders,
conservation identities). The public sources below pin the **forms** of the chapter's vortex models and of the two
transport theorems, and one **published number** (the Lamb–Oseen radius of maximum speed). Every row was read from the
cited page by the math-verifier on the date given (page fetched and its text read). No value comes from the textbook;
book-quoted numbers stay in the git-ignored `tests/book_values_ch03.json`.

| value / key in `benchmarks.json` | what it is | source | URL | how obtained | verified |
|---|---|---|---|---|---|
| `swirl_lamb_oseen` → α = 1.256 (3 digits), v_θ = v_max(1 + 1/2α)(r_max/r)[1 − exp(−α r²/r_max²)] | Lamb–Oseen profile written with the radius of maximum speed: r_max² = α σ² (σ the Gaussian core radius of the book's (3.29)), so r_max/σ = √α; the normalisation (1 + 1/2α)(1 − e^{−α}) = 1 is exactly the book's condition e^x = 1 + 2x | J. R. Canivete Cuissa & O. Steiner, "An innovative and automated method for vortex identification. I. Description of the SWIRL algorithm", A&A (2022), Sect. 2.4, Eq. (10) | https://arxiv.org/abs/2210.05223 | PDF downloaded, text searched ("α = 1.256 establishes that r_max is the radius at which the rotational velocity v_θ is maximal") | 2026-09-23 |
| `rankine_vortex` → v_θ = (Γ/2π) r/a² (r ≤ a), (Γ/2π)/r (r > a); ω_z = 2Ω inside, 0 outside | the Rankine profile of the book's (3.28) (a = σ) | Wikipedia, "Rankine vortex" | https://en.wikipedia.org/wiki/Rankine_vortex | page fetched, formulas read | 2026-09-23 |
| `lamb_oseen_vortex` → v_θ = Γ/(2πr)(1 − e^{−r²/4νt}), ω_z = Γ/(4πνt) e^{−r²/4νt} | the Lamb–Oseen vortex; with σ² = 4νt it is the book's Gaussian vortex (3.29) (form only; the page gives no r_max) | Wikipedia, "Lamb–Oseen vortex" | https://en.wikipedia.org/wiki/Lamb%E2%80%93Oseen_vortex | page fetched, formulas read | 2026-09-23 |
| `reynolds_transport_theorem` → d/dt ∫_Ω(t) f dV = ∫ ∂f/∂t dV + ∮ (v_b·n) f dA, n outward, v_b the boundary velocity | the form and sign of the boundary term of (3.35); fixed region → d/dt passes inside | Wikipedia, "Reynolds transport theorem" | https://en.wikipedia.org/wiki/Reynolds_transport_theorem | page fetched, statement read | 2026-09-23 |
| `leibniz_integral_rule` → d/dx ∫_{a(x)}^{b(x)} f dt = f(b) b′ − f(a) a′ + ∫ ∂f/∂x dt | the signs of the two end terms of (3.30) | Wikipedia, "Leibniz integral rule" | https://en.wikipedia.org/wiki/Leibniz_integral_rule | page fetched, formula read | 2026-09-23 |

## Library / closed-form cross-checks used as evidence (not literature)
| check | what | how it enters the tests |
|---|---|---|
| `scipy.special.lambertw(z, k=-1)` | x* = −W₋₁(−e^{−1/2}/2) − ½ = 1.2564312086… | `gaussian_vortex_max_radius(method="lambertw")` vs brentq (1e-13) and vs a numerical maximisation of u_θ (`scipy.optimize.minimize_scalar`) |
| `scipy.linalg.expm`, `numpy.linalg.svd`, `numpy.linalg.det` | exact linear-flow maps, ellipse axes, Jacobi's formula | path lines of linear flows, strain ellipses, volume ratios |
| `scipy.integrate.solve_ivp` (independent of the stencils) | trajectories | material derivative = d/dt F along the path line |

## Rejected or not used as V5
| candidate | why |
|---|---|
| NASA NTRS 20140003974 (Ahmad, Proctor et al., wake-vortex model review) | not needed: the Lambert-W closed form (V1) already gives x* to 1e-13, and the SWIRL paper supplies the published 3-digit value |
| Any numeric table for Examples 3.1–3.2 or §3.5 | book-derived; private `tests/book_values_ch03.json` (V6), never in `reference/` |
