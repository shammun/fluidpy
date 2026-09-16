# reference/ch02 — sources of the public benchmark data

`benchmarks.json` is written by `make_refs.py`. Chapter 2 (Cartesian tensors) is pure mathematics: the evidence in
`tests/test_ch02.py` is overwhelmingly Tier 1 (analytic identities, sympy re-derivations, convergence orders). The
public sources below pin the **conventions** the chapter fixes for the whole project (index placement of Cauchy's
formula, active vs passive rotation matrices, the right-hand rule of Stokes' theorem) and two closed-form numbers.
Every row was read from the cited page by the math-verifier on the date given (page fetched and its text read).
No value comes from the textbook; book-quoted numbers stay in the git-ignored `tests/book_values_ch02.json`.

| value / key in `benchmarks.json` | what it is | source | URL | how obtained | verified |
|---|---|---|---|---|---|
| `divergence_theorem_sphere_example` → flux = 8π/3 | worked example: F = 2x i + y² j + z² k through the unit sphere; ∇·F = 2(1 + y + z); odd terms vanish over the ball; flux = 2·(4π/3) | Wikipedia, "Divergence theorem", section *Example* | https://en.wikipedia.org/wiki/Divergence_theorem | page fetched, example text read | 2026-09-16 |
| `levi_civita_identities` → ε_ijk ε_imn = δ_jm δ_kn − δ_jn δ_km; ε_jmn ε_imn = 2δ_ij; ε_ijk ε_ijk = 6; (a×b)_i = ε_ijk a_j b_k; det A = ε_ijk a_1i a_2j a_3k; pseudotensor | the product identities and the determinant formula in three dimensions | Wikipedia, "Levi-Civita symbol", section *Three dimensions* | https://en.wikipedia.org/wiki/Levi-Civita_symbol | page fetched, identities read | 2026-09-16 |
| `cauchy_stress_tensor` → T_j = σ_ij n_i; σ' = A σ Aᵀ; I₁, I₂, I₃ (components and principal values); 2-D principal stresses; τ_max = ½\|σ₁ − σ₃\| | Cauchy's stress theorem (first index contracted with n — the book's (2.15)), transformation rule in the *active* convention (A = Cᵀ), invariants, Mohr's formula | Wikipedia, "Cauchy stress tensor" | https://en.wikipedia.org/wiki/Cauchy_stress_tensor | page fetched, formulas read | 2026-09-16 |
| `rotation_matrix` → R(θ) = [[cos, −sin],[sin, cos]] active (counterclockwise); passive = Rᵀ; RᵀR = I, det R = +1; Rodrigues R = cos θ I + sin θ [u]ₓ + (1 − cos θ) u uᵀ | the active/passive (alibi/alias) distinction the book's C relies on: C for a frame rotated by +θ equals R(θ) and is applied as Rᵀ | Wikipedia, "Rotation matrix" | https://en.wikipedia.org/wiki/Rotation_matrix | page fetched, text read | 2026-09-16 |
| `stokes_theorem` → ∮ F·dr = ∬ (∇×F)·n dS, right-hand rule | statement and orientation rule of the Kelvin–Stokes theorem; the page has **no numeric example**, so the V1 case in the tests (solid-body rotation, circulation 2\|b\|πR²) is our own closed form | Wikipedia, "Stokes' theorem" | https://en.wikipedia.org/wiki/Stokes%27_theorem | page fetched, statement read | 2026-09-16 |

## Library cross-checks used as evidence (not literature)
| check | what | how it enters the tests |
|---|---|---|
| `scipy.spatial.transform.Rotation.from_rotvec(θ k).as_matrix()` | the *active* rotation matrix about axis k | `rotation_matrix_3d(k, θ)` must equal it entrywise (the book's passive C is the same matrix applied as Cᵀ) — `test_direction_cosines_V5_scipy_active_matrix` |
| `numpy.cross`, `numpy.linalg.eigh`, `numpy.linalg.eigvalsh`, `numpy.roots`, `numpy.linalg.det` | reference implementations | `cross`/`cross_einsum`, `principal_axes`, `characteristic_polynomial`, the ε determinant formula |

## Rejected or not used as V5
| candidate | why |
|---|---|
| Sommerfeld, *Mechanics of Deformable Bodies* (1964), p. 59 (the book's citation for the tetrahedron argument) | bibliographic only, not fetched; D05/D06 are re-derived symbolically instead (V2) |
| Any numeric table for Examples 2.1–2.6 | book-derived; private `tests/book_values_ch02.json` (V6), never in `reference/` |
