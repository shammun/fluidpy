# Chapter 6 reference data — sources

Produced by `reference/ch06/make_refs.py` → `benchmarks.json`. Public sources only; every page was fetched and read on
the date given. Book-quoted numbers and printed closed forms are **not** here (they live in the git-ignored
`tests/book_values_ch06.json`). Label use in `tests/test_ch06.py`: `added_mass_sphere` and `rayleigh_collapse` are
published *numbers* our code reproduces from its own derivation (V5); every other entry is a published closed form that
our code reproduces identically, i.e. a V1 "form cross-check" (not V5), and the Hess & Smith entry is the method's
pedigree only.

| value / form | what it is | source | URL | how obtained | verified | label |
|---|---|---|---|---|---|---|
| `added_mass_sphere` | M = (2/3)πr³ρ = ½ × displaced fluid mass; (m_p + ρV_p/2) dv/dt = ΣF | primary: H. Lamb, *Hydrodynamics*, 6th ed., Cambridge Univ. Press (1932), Ch. V §92 (section number as cited by the derivation review; the text was not fetched); secondary (fetched, values read): Wikipedia, "Added mass" | https://en.wikipedia.org/wiki/Added_mass | page read; formula and "half the volume … times the density" quoted | 2026-09-23 | V5 (coefficient ½) |
| `rayleigh_collapse` | empty-cavity filling time t ≈ 0.91468 R₀√(ρ/P∞) (Besant 1859, Rayleigh) | primary: Lord Rayleigh, "On the pressure developed in a liquid during the collapse of a spherical cavity", Phil. Mag. (6) 34(200), 94–98 (1917), doi:10.1080/14786440808635681 (citation confirmed by search, https://www.tandfonline.com/doi/abs/10.1080/14786440808635681); secondary (fetched, constant read): Wikipedia, "Rayleigh–Plesset equation" | https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Plesset_equation | page read; constant quoted (5 significant figures) | 2026-09-23 | V5 (0.91468) |
| `cylinder_form` | ϕ = Ur(1 + R²/r²)cos θ, V_r = U(1 − R²/r²)cos θ, V_θ = −U(1 + R²/r²)sin θ; surface C_p from +1 to −3; zero drag (d'Alembert) | Wikipedia, "Potential flow around a circular cylinder" | https://en.wikipedia.org/wiki/Potential_flow_around_a_circular_cylinder | page read; formulas and range quoted | 2026-09-23 | V1 form |
| `kutta_joukowski_form` | L′ = ρ∞V∞Γ, Γ taken on a clockwise contour; F̄ = (iρ/2)∮w′²dz | Wikipedia, "Kutta–Joukowski theorem" | https://en.wikipedia.org/wiki/Kutta%E2%80%93Joukowski_theorem | page read; formula and the clockwise-contour sentence quoted | 2026-09-23 | V1 form |
| `blasius_form` | F_x − iF_y = (iρ/2)∮(dw/dz)²dz; M = Re{−(ρ/2)∮z(dw/dz)²dz} | Wikipedia, "Blasius theorem" | https://en.wikipedia.org/wiki/Blasius_theorem | page read; formulas quoted | 2026-09-23 | V1 form |
| `joukowsky_form` | z = ζ + 1/ζ; unit circle → flat plate from −2 to +2; W = W̃/(1 − 1/ζ²) | Wikipedia, "Joukowsky transform" | https://en.wikipedia.org/wiki/Joukowsky_transform | page read; formulas quoted | 2026-09-23 | V1 form |
| `rankine_half_body_form` | ψ = Ur sin θ + mθ/2π; b = m/2πU; body r = b(π − θ)/sin θ | Wikipedia, "Rankine half body" | https://en.wikipedia.org/wiki/Rankine_half_body | page read; formulas quoted (the page gives no far-downstream width) | 2026-09-23 | V1 form |
| `sphere_form` | Φ = v(r + a³/2r²)cos θ; u = v(1 − a³/r³)cos θ r̂ − v(1 + a³/2r³)sin θ θ̂; cylindrical components; Stokes ψ convention u_z = (1/ϱ)∂Ψ/∂ϱ, u_ϱ = −(1/ϱ)∂Ψ/∂z | K. T. McDonald, "Pressure in Fluid Flow Past a Sphere", Princeton (March 19, 2015; updated March 25, 2015), eqs. (9), (10), (12), (14) | http://kirkmcd.princeton.edu/examples/sphereflow.pdf | PDF fetched, text extracted with pymupdf, equations read | 2026-09-23 | V1 form |
| `source_panel_method` | constant-strength source panels, zero normal velocity at the panel midpoints (method only) | J. L. Hess & A. M. O. Smith, Prog. Aerosp. Sci. 8, 1–138 (1967), doi:10.1016/0376-0421(67)90003-6 | https://www.sciencedirect.com/science/article/abs/pii/0376042167900036 | citation confirmed by search (ScienceDirect, ADS 1967PrAeS...8....1H) | 2026-09-23 | pedigree |

Sign conventions pinned by the sources: Wikipedia's Kutta–Joukowski page takes the circulation on a **clockwise**
contour, the book's convention in (6.36)–(6.40) and (6.62) (`Gamma_cw=` in `fluidpy`), so L′ = +ρUΓ_cw; McDonald's
Stokes stream function has the sign of the book's (6.75) (`core.potential.AxisymFlow`). Not verified and therefore not
used: a measured high-Reynolds-number cylinder pressure distribution (Fig. 6.10's dashed curve; candidates Roshko 1961,
Achenbach 1968) — the notebook's "real" band stays labelled *qualitative*.
