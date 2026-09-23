# Chapter 5 reference data — sources

Produced by `reference/ch05/make_refs.py` → `benchmarks.json`. Public sources only; every page was fetched and read on
the date given. Book-quoted numbers and printed closed forms are **not** here (they live in the git-ignored
`tests/book_values_ch05.json`). Label use in `tests/test_ch05.py`: published numbers or published results used as such
are V5; a page that only restates an analytic form we also derive is used as a V1 "form cross-check" (not V5).

| value / form | what it is | source | URL | how obtained | verified | label |
|---|---|---|---|---|---|---|
| `kelvin_ring_speed` | U = Γ/(4πR)[ln(8R/a) − 1/4], thin core of uniform vorticity (Kelvin) | Wikipedia, "Vortex ring" | https://en.wikipedia.org/wiki/Vortex_ring | page read; formula and "thin-core, uniform" qualifier quoted | 2026-09-23 | V5 |
| `point_vortex_pairs` | corotating pair Γ₁ = Γ₂ = 1 at ±l rotates at Ω₀ = 1/(4πl²); counter-rotating pair (1 at +l, −1 at −l) translates at V₀ = −i/(4πl) | C. García & S. V. Haziot, "Global bifurcation for corotating and counter-rotating vortex pairs", Commun. Math. Phys. (2023), doi:10.1007/s00220-023-04741-6 | https://arxiv.org/html/2204.11327 (§2.1) | HTML full text read; formulas quoted | 2026-09-23 | V5 |
| `burgers_vortex_form` | v_r = −α_w r, v_z = 2α_w z, v_θ = (Γ/2πr)[1 − e^{−α_w r²/2ν}], ω_z = (α_wΓ/2πν)e^{−α_w r²/2ν}; α_w = α/2 of the book's Exercise 5.12 form | Wikipedia, "Burgers vortex" | https://en.wikipedia.org/wiki/Burgers_vortex | page read; formulas quoted; parameter map derived by us | 2026-09-23 | V1 form |
| `hill_vortex_form` | ψ_in = −(3U/4)(1 − r²/a²)r² sin²θ, ψ_out = (U/2)(1 − a³/r³)r² sin²θ, ω_φ = −(15U/2a²) r sin θ | Wikipedia, "Hill's spherical vortex" | https://en.wikipedia.org/wiki/Hill%27s_spherical_vortex | page read; formulas quoted (the relation A = 15U/2a² in magnitude is ours) | 2026-09-23 | V1 form |
| `biot_savart_segment_form` | infinite line Γ/(2πr); segment Γ/(4πr)(cos A − cos B) | Wikipedia, "Biot–Savart law" (Aerodynamics) | https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law | page read; formulas quoted | 2026-09-23 | V1 form |
| `kelvin_bjerknes_form` | Kelvin: DΓ/Dt = 0 for a material curve in a barotropic ideal fluid with conservative body forces; Poincaré–Bjerknes: ∫(ω + 2Ω)·n dS conserved in a frame rotating at constant Ω | Wikipedia, "Kelvin's circulation theorem" | https://en.wikipedia.org/wiki/Kelvin%27s_circulation_theorem | page read; statements quoted | 2026-09-23 | V1 form |

Reused (already verified): WGS-84 Ω = 7.292115e-5 rad/s (`reference/ch04/benchmarks.json`, `core.rotating.OMEGA_EARTH`).

Sign conventions pinned by the sources: García & Haziot use counterclockwise-positive circulation (Ω₀ > 0 for a = 1),
the same as `fluidpy` (§9 of `analysis/ch05.md`); their counter-rotating pair with +1 on the **right** moves in −y, our
`opposite_pair` preset (+1 on the **left**) moves in +y — the mirror image, consistent.
