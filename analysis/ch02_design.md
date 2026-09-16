# Chapter 2 — Cartesian Tensors: lesson design
(from `analysis/ch02_curation.md` (A 16 · B 60 · C 18 = 94 rows; CORE 16 · NOTE 74 · RECAP 2 · SKIP 2; 15 derivations
written out; E1–E5 + backup B1) and `analysis/ch02.md` (§2b derivations, §4 rows F1–F26, §9 conventions and typos);
every equation re-read on the rendered pages p067–p089 (`tools/render_pages.py ch02 --eq 2.1 … 2.36`); 2026-09-16,
lesson-designer. The implementer works in parallel from analysis §4 + curation §8; Part C is the contract both sides
keep.)

**Binding conventions for every builder (analysis §9, curation decision 5).**
1. **Passive rotation.** $C_{ij} = \mathbf e_i\cdot\mathbf e'_j$ (row = old axis, column = new axis; the *columns* of C
   are the new unit vectors in old components). Components transform as $\mathbf x' = \mathbf C^{\rm T}\mathbf x$,
   $\mathbf x = \mathbf C\,\mathbf x'$, $\boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau\,\mathbf C$. The book's C for a
   frame turned by +θ *equals* Wikipedia's / scipy's active R(θ) but is applied as Rᵀ. Every C-returning function's
   docstring says "passive: columns = new axes"; the notebook has one ⚠️ callout (C03) and E1 a "what rotates?" toggle.
2. **Traction contracts the first index**: $f_i = \tau_{ji}n_j$, $\mathbf f = \mathbf n\cdot\boldsymbol\tau$; equals
   $\boldsymbol\tau\cdot\mathbf n$ only for symmetric τ (Ch. 4 proves symmetry; we never assume it in code: `traction`
   uses `einsum('ji,j->i')`). Divergence of a tensor contracts the **second** index, $(\nabla\cdot\boldsymbol\tau)_i =
   \partial\tau_{ij}/\partial x_j$. Velocity gradient `G[i, j] = ∂u_i/∂x_j`. Book double dot $\mathbf A:\mathbf B =
   A_{ij}B_{ji}$; Frobenius $A_{ij}B_{ij}$ — `double_dot(A, B, convention=)` with no silent default in prose.
3. **Stress sign**: tensile positive; on the face with outward normal $+\mathbf e_i$ the positive components point along
   $+\mathbf e_j$; on the $-\mathbf e_i$ face they reverse (so the traction there is $-\tau_{i\cdot}$).
4. **Stokes orientation**: choose the outside of A → n; t runs counterclockwise about n (right-hand rule: fingers along
   t, thumb along n); $\mathbf n_c \equiv \mathbf n\times\mathbf t$ is the in-surface normal to C pointing *into* A, so
   the book's $\mathbf t = \mathbf n_c\times\mathbf n$ holds and $(\mathbf n_c, \mathbf n, \mathbf t)$ is right-handed.
   `planar_loop` is counterclockwise about `normal` and satisfies $\oint(\mathbf x-\mathbf c)\times\mathbf t\,ds = +2A\,\mathbf n$
   (the analyst's F25 test). Flipping n flips both sides of (2.34).
5. **Grid layout for the whole project**: 3-D arrays `[k, j, i]` = (z, y, x) — x on the **last** axis, built with
   `np.meshgrid(z, y, x, indexing='ij')`; 2-D `[j, i]` = (y, x) with `indexing='xy'`; vector fields stack the component
   on axis 0 (`u[c, k, j, i]`). Stencils: second-order central inside, second-order one-sided at the edges (never
   `np.gradient`).
6. **Ex. 2.4's Γ is the off-diagonal element $S_{12}$** (analysis §9 typo 4): `example_2_4(Gamma)` takes Γ = S₁₂ and the
   notebook flags that for u₁(x₂) alone S₁₂ = ½ du₁/dx₂ (E3's simple-shear preset with du₁/dx₂ = 2 s⁻¹ gives S₁₂ = 1
   s⁻¹ = Ex. 2.4's Γ). Other typos taught corrected: Ex. 2.3's third term is $a x_3\mathbf e_3$; Ex. 2.6's second
   bracket holds $u_y$; (2.27)'s lower limit is $i = 1$; Ex. 2.2's "[?]" is ignored; Ex. 2.4's "E_ij" is S_ij.
7. **Colours (text, figures, explainers):** old frame / Cartesian `teal` · new (rotated) frame `orange` · normal stress
   `blue` · shear stress `rose` · symmetric part `teal` · antisymmetric part `orange` · flux in `blue` / out `orange` ·
   circulation `accent` purple · curl positive purple / negative green · ghosts and references `muted` grey.
8. **Depth.** A items are the full nine-part CORE block (Part A). B items are `nb.note` paragraphs (plain words, the
   equation, one number, one code line) inside the A block the curation names; C items one `nb.note` sentence with a
   pointer; R → `nb.recap`; S → `nb.pointer`. The 11 demoted derivations (curation §4c) are **stated, never derived**;
   the 15 D rows are written out in Part F and copied word for word into `nb.derivation` and the explainer Derivation
   tabs.
9. **Symbol collisions said once, where they bite** (analysis §9): Γ here = shear rate / S₁₂ (code `Gamma_shear`), not
   ch01's lapse rate nor Ch. 5's circulation (code `circulation`); ω = the vector of an antisymmetric tensor (vorticity
   in Ch. 3), not ch07's frequency; λ^k is an eigenvalue *label* (`lam[k]`), not a power; A = antisymmetric part *and*
   surface area *and* a generic tensor — the notebook writes 𝐀 for the tensor and $A$ (area) only inside integrals.
10. **Parsing.** Part E is the only part with table rows after the Part E heading; Part F has no line starting with `|`
    (`tools/coverage_check.py` reads only the ledger). Explainer headings are exactly `### E1 · rotation_of_axes` …
    `### B1 · index_machine` (`tools/embed_check.py`).

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---

## Part C — functions the builders will call (the implementer's contract)

**Status column.** **§4 F·n** = planned in `analysis/ch02.md` §4 with this signature (keep it). **§8** = added by the
curation's notes for the implementer. **NEW** = added by this design (not in analysis §4 nor curation §8) — flagged as
the phase asks. Every callable is reachable as `ch02.<name>` (the chapter module `fluidpy/ch02_cartesian_tensors.py`
re-exports `core.tensors`, `core.index_notation`, `core.grids`, `core.operators`, `core.integral_theorems`), which is
what `tools/shot.py` parity rows need: `py:` expressions are evaluated with **no builtins**, so rows index dicts / tuples
down to one number (`ch02.traction_2d([[0,1],[1,0]], 0.5236)["sigma_n"]`, `ch02.principal_axes([[0,1],[1,0]])[0][1]`).
Scalar in → float out (`core._util.as_scalar_if_0d`); lists accepted wherever arrays are (`np.asarray` first). Angles in
radians; `_deg` only at the interface. This chapter is pure mathematics: stresses in Pa when physical, strain rates in
1/s, lengths in m; docstrings cite § and Eq. and carry the validation label.

### C.0 Machinery (existing, `fluidpy/core/`)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST), `style.COLORS`, `style.savefig(fig, "ch02", name)` | FAST · palette · path | setup, every figure | WIP |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video"\|"frames")` | HTML | C02, C10/C11, C12, C15, C16 | WIP |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.animate_figure` · `interact.live(fn, **widgets)` | plotly Figure · widget | C05, C06, C09, C13, C14, C15 · C13 live | WIP |
| 0.4 | `embed.show_viz("ch02", slug)` | display | 5 explainer cells | WIP |
| 0.5 | `tools.convergence.observed_order(h, err) -> float` | slope | C09, C15, C16 | WIP |
| 0.6 | `core._util.as_scalar_if_0d`, `require_positive` | – | every function | WIP |

### C.1 `fluidpy/core/index_notation.py` (NEW module, §4 F2)
| # | Callable (signature) | Implements | Returns | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `expand_indices(term: str, dim: int = 3, symbols: dict \| None = None) -> sympy.Expr \| sympy.Array` | the summation convention: a repeated letter is summed 1..dim; `delta_ij`, `eps_ijk` evaluate to numbers; `u_i,i` comma → `Derivative(u_i(x1,x2,x3), x_i)`; free indices → an `Array` of expressions | expression / array | C01 (N07 N21 N23 N41), C03 (N13), C08 (N49 N50), §2.14 (N74), B1 | §4 F2 |
| 1.2 | `expand_indices_str(term: str, dim: int = 3) -> str` | same, rendered as a plain string `"a_1*b_1 + a_2*b_2 + a_3*b_3"` (the form B1 mirrors and the notebook prints) | str | C01, B1 | **NEW** (curation §8.10 asks for a string result) |
| 1.3 | `classify_indices(term: str) -> tuple[list[str], list[str]]` | (free, dummy) letters; raises if a letter appears 3× | tuple | C01 (N12), B1 | §4 F2 |
| 1.4 | `tensor_order(term: str) -> int` | number of free indices | int | C01 (N08) | §4 F2 |
| 1.5 | `rename_dummy(term: str, old: str, new: str) -> str` | renames a summed letter (raises for a free one) | str | C01 (N12) | §4 F2 |
| 1.6 | `comma_to_partial(term: str) -> str` | `"u_i,j"` → `"∂u_i/∂x_j"` text | str | §2.14 (N74) | §4 F2 |

### C.2 `fluidpy/core/tensors.py` (NEW module; §4 F1, F3–F6, F8–F10, F12–F14, F18–F21 + §8 additions)
| # | Callable (signature) | Implements (Eq.) | Units / returns | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `unit_vectors(n=3) -> ndarray (n,n)` (rows e_i) · `vector_from_components(comp, E) -> ndarray` (= comp @ E) | (2.1), (2.3) | any | C01 (N03), C02 (N10) | §4 F1 |
| 2.2 | `dot(a, b) -> float` · `inner(A, B) -> ndarray` · `outer(u, v)` · `tensor_product(A, B)` (= `np.multiply.outer`) · `contract(A, B, pattern: str)` (einsum string, e.g. `"ij,ki->kj"`) · `dot_tensor_vector(A, u, index: int)` (index 1 → A·u, 0 → Aᵀ·u) · `double_dot(A, B, convention: str = "book")` (`"book"` = A_ij B_ji, `"frobenius"` = A_ij B_ij; the WIP defaults to the book's — the notebook and explainers always pass it explicitly) · `trace(A) -> float` | (2.2), (2.9), (2.14), §2.5 text | any | C01, C06 (N29), C07 (N30–N33) | §4 F3 |
| 2.3 | `direction_cosines(E_old, E_new) -> C` (C_ij = e_i·e'_j) · `rotation_matrix_2d(theta) -> (2,2)` · `rotation_matrix_3d(axis, angle) -> (3,3)` (Rodrigues) · `rotation_angle(C) -> float` · `random_rotation(rng) -> C` (QR of a Gaussian matrix, det fixed to +1) | C02, N15, Ex. 2.1 | rad | C02, C03, C06, C07, C08, E1 | §4 F4 |
| 2.4 | `is_orthogonal(C, tol=1e-12) -> bool` · `is_proper_rotation(C, tol=1e-12) -> bool` · `orthogonality_residual(C) -> float` (max|CᵀC − I|) | D02 (N15) | – | C02, E1 status | §4 F5 |
| 2.5 | `transform_vector(u, C) -> u'` (= Cᵀu, (2.5)/(2.8)) · `inverse_transform_vector(up, C)` (= C u', (2.7)) · `transform_tensor(T, C)` (any order, generated einsum `'ia,jb,…,ij…->ab…'`, (2.12)/(2.13)) · `transforms_as_vector(component_fn, C, points=None, vector_params=(), rng=None) -> float` (residual; WIP signature) · `transforms_as_tensor(fn, C, points=None, order=2) -> float` | (2.5), (2.7), (2.8), (2.12), (2.13) | – | C03 (N16 N18), C06 (N27 N28 N29), E1 | §4 F6 |
| 2.6 | `cube_face_tractions(tau) -> dict[str, tuple[ndarray, ndarray]]` (face keys `"+1"`, `"-1"`, … `"-3"` (WIP) → (outward normal, traction ±τ[i, :])) · `stress_component_meaning(i, j) -> str` | §2.4 sign convention (C04, R01) | Pa | C04, E2 element | §4 F8 |
| 2.7 | `tetrahedron_face_areas(n, dA) -> ndarray(3)` (dA_i = n_i dA) | N35, D05 | m² | C05 | §4 F9 |
| 2.8 | `traction(tau, n) -> f` (`einsum('ji,j->i')`; accepts 2×2 and 3×3) · `normal_shear_stress(tau, n) -> (sigma_n, tau_s, shear_dir)` (calls the promoted `traction_components`) | (2.15) | Pa | C05, C13, E2 | §4 F10 (+§8.1 2×2, §8.8 promotion of `ch01.traction_components` into `core/tensors.py` with a re-export from `ch01_introduction`) |
| 2.9 | `mohr_circle_2d(tau) -> (center, radius)` (= (I₁/2, √(((τ₁₁−τ₂₂)/2)² + τ₁₂²))) | Mohr preview (N38, C13) | Pa | C06 slider, C13, E2 | §8.1 |
| 2.10 | `kronecker_delta(n=3)` · `levi_civita() -> (3,3,3)` · `permutation_sign(i, j, k) -> int` (0-based or 1-based accepted, documented) · `epsilon_delta_residual() -> float` (max over 81 cases of |Σ_k ε_ijk ε_klm − (δ_il δ_jm − δ_im δ_jl)|) · `triple_product(a, b, c) -> ndarray` (a × (b × c)) | (2.16)–(2.19) | – | C01 (N40 N41), C08 (N44), B1 | §4 F12 |
| 2.11 | `is_isotropic(T, rng=None, n_rotations=50, proper=True) -> tuple[bool, float]` | N42 | – | C08 | §4 F13 |
| 2.12 | `cross(u, v)` ((2.20) explicit) · `cross_einsum(u, v)` (`einsum('ijk,i,j->k', eps, u, v)`) · `angle_between(u, v) -> float` (arctan2(|u×v|, u·v)) | (2.20), (2.21), N45 | rad | C08 (N46–N50), C01 (N45) | §4 F14 |
| 2.13 | `is_symmetric(B, tol)`, `is_antisymmetric(B, tol)`, `independent_components(B) -> int`, `symmetric_part(B)`, `antisymmetric_part(B)`, `strain_rate_tensor(G)` (= symmetric part), `rotation_tensor(G)` (= G − Gᵀ = 2A: the book's R with vector ω = ∇×u; A = ½R has vector ½∇×u — review M1) | §2.10, D14 | any | C12 (N56), E3 | §4 F18 |
| 2.14 | `antisymmetric_from_vector(omega) -> R` (`-einsum('ijk,k->ij', eps, omega)`) · `vector_from_antisymmetric(R) -> omega` (`-0.5*einsum('ijk,ij->k', eps, R)`) | (2.26), (2.27), D15 | any | C12 (N57 N58), E3 | §4 F19 |
| 2.15 | `symmetric_double_contraction(tau, B) -> (P, P_S, P_A)` (Frobenius sums; P_A = 0 for symmetric τ) | (2.28), (2.29) | any | C12 (N59 N60) | §4 F20 |
| 2.16 | `principal_axes(tau) -> (lam, B)` (λ ascending; B columns = unit eigenvectors; det B = +1; **raises** for non-symmetric input; stable for repeated λ) · `characteristic_polynomial(tau) -> ndarray(4)` (coefficients of λ³ − I₁λ² + I₂λ − I₃) · `diagonalize(tau) -> (C, tau_prime)` (C = B; τ' = CᵀτC) · `invariants(A) -> (I1, I2, I3)` (I₂ = ½(I₁² − A_ij A_ji)) · `normal_stress_bounds(tau, rng=None, n=1000) -> (min_normal, max_normal, max_shear)` (Monte-Carlo unit normals) | §2.11 facts (D17), Exercise 2.9 (D18) | Pa | C07, C13, E2 | §4 F21 (+§8.6) |

### C.3 `fluidpy/core/grids.py` (NEW module, §4 F15)
| # | Callable (signature) | Implements | Used by | Status |
|---|---|---|---|---|
| 3.1 | `grid(bounds, n) -> Grid` with fields `x, y, z` (1-D), `X, Y, Z` (3-D, layout `[k, j, i]`), `h` (per-axis spacing, tuple) | the project layout (convention 5) | C09–C11, C14, C15, tests | §4 F15 |
| 3.2 | `grid2d(bounds, n) -> Grid2D` (`x, y, X, Y, h`; `[j, i]`, `indexing='xy'`) | 2-D figures and E4/E5 parity | C09 figures, C10/C11 heatmaps, C16 | §4 F15 (+§8.5) |

### C.4 `fluidpy/core/operators.py` (NEW module, §4 F16)
| # | Callable (signature) | Implements (Eq.) | Used by | Status |
|---|---|---|---|---|
| 4.1 | `partial(f, direction: int, h: float, bc: str = "onesided") -> ndarray` (second-order central inside, second-order one-sided at the edges; `bc="periodic"` option) — **`direction` is the coordinate direction (0 = x, 1 = y, 2 = z)**; the function maps it to the array axis of the `[k, j, i]` layout itself (WIP signature `partial(f, direction, h, bc, ndim_space=None)`) | ∂/∂x_i on the grid (R02) | C09 from-scratch, every operator | §4 F16 (WIP ✓) |
| 4.2 | `gradient(phi, h) -> ndarray (ndim, …)` (component axis 0, ordered (x, y[, z])) · `directional_derivative(grad_phi, n) -> ndarray` (∇φ·n) | (2.22), ∂φ/∂n = ∇φ·n | C09 | §4 F16 |
| 4.3 | `divergence(u, h) -> ndarray` · `vector_gradient(u, h) -> G` with `G[i, j, …] = ∂u_i/∂x_j` · `tensor_divergence(T, h, index: int = 1) -> ndarray` (contracts index 1 = the second index by default) | (2.23), N52 | C10, E4 heatmap | §4 F16 |
| 4.4 | `curl(u, h) -> ndarray (3, …)` (`einsum('ijk,kj...->i...', eps, G)`) · `curl_components(u, h) -> tuple` ((2.25) written out) · `is_solenoidal(u, h, tol=1e-8) -> bool` · `is_irrotational(u, h, tol=1e-8) -> bool` · `laplacian(phi, h)` (bonus, ch04+) | (2.24), (2.25), N54 | C11, C16, E5 | §4 F16 |

### C.5 `fluidpy/core/integral_theorems.py` (NEW module, §4 F23–F26 + §8.3)
| # | Callable (signature) | Implements (Eq.) | Used by | Status |
|---|---|---|---|---|
| 5.1 | `volume_integral_box(f_fn, bounds, n) -> float` (composite midpoint, 3-D) · `gauss_gradient_box(Q_fn, bounds, n) -> (lhs, rhs)` (scalar Q: both sides are 3-vectors) · `divergence_theorem_box(Q_fn, bounds, n, div_fn=None, rule="midpoint") -> (lhs, rhs)` · `flux_through_box_faces(Q_fn, bounds, n, rule) -> dict` (six faces with Fig. 2.4's letters EADH … DCGH — the D22 animation reads its values) · `flux_through_sphere(Q_fn, R, n_theta, n_phi) -> float` · `divergence_theorem_sphere(Q_fn, R, n) -> (lhs, rhs)` | (2.30), N64, benchmark 8π/3 | C14, E4 | §4 F23 |
| 5.2 | `divergence_theorem_rect2d(Q_fn, bounds, n, div_fn=None) -> (lhs, rhs, face_fluxes: dict)` · `flux_through_faces(Q_fn, bounds, n) -> dict[str, float]` (face keys `"+x" "-x" "+y" "-y"` (WIP), signed n·Q × length; 2-D = unit depth) · `integral_divergence_2d(Q_fn, x0, h, n_side=8) -> float` | (2.30), (2.32) in 2-D for the page figures and E4 parity | C14 figure, C15, E4 | §8.3 |
| 5.3 | `divergence_theorem_tiled(Q_fn, bounds, tiles: int, n) -> (sum_of_tile_boundary_fluxes, outer_flux, interior_cancellation)` (2-D; splits the box into tiles × tiles sub-boxes and calls `flux_through_faces` on each; the third number is Σ over interior faces, expected 0 to round-off) | D25 step 8 on screen (interior faces cancel) | C14 figure, E4 parity | **NEW — not in the WIP yet** (≈ 15 lines on top of `flux_through_faces`) |
| 5.4 | `integral_gradient(Q_fn, x0, h, n_face=8) -> ndarray` · `integral_divergence(Q_fn, x0, h, n_face=8) -> float` · `integral_curl(Q_fn, x0, h, n_face=8) -> ndarray` (six-face midpoint sums on a cube of side h divided by h³; error O(h²)) | (2.31), (2.32), (2.33) | C15 (N66 N67), E4 | §4 F24 |
| 5.5 | `boundary_tangent(n_c, n) -> t` (= n_c × n, unit) · `planar_loop(center, normal=None, radius, n) -> Loop` (fields `points, tangents, ds, length, area, center, normal, dim`; counterclockwise about `normal`; **a 2-D `center` with `normal=None` gives a 2-D loop for the plane `VectorField`s, a 3-D `center` a 3-D loop** — WIP) · `rectangle_loop(center, normal=None, a, b, n) -> Loop` · `planar_disc(center, normal=None, radius, nr, ntheta) -> Surface` | N69, Fig. 2.10 (convention 4) | C16, E5 | §4 F25 (+§8.3) |
| 5.6 | `circulation(u_fn, loop) -> float` (midpoint sum Σ u·t ds; `u_fn` returns 2 or 3 components) · `curl_flux(u_fn, surface, curl_fn=None, h_curl=1e-3) -> float` · `stokes_theorem_check(u_fn, loop, surface, curl_fn=None) -> StokesCheck(lhs, rhs, hypothesis_ok, note)` (a dataclass that iterates as (lhs, rhs); `hypothesis_ok=False` with a `note` when the field's `singular_at` lies inside the loop — the gap lhs − rhs is then reported, never asserted; parity rows use `.hypothesis_ok`) · `integral_curl_component(u_fn, x0, n, h, n_side=8) -> float` (square loop of side h ⊥ n; (2.35)) | (2.34), (2.35), Ex. 2.6 | C16 (N71 N73), C11, E5 | §4 F26 (+§8.4 `singular_at`) |

### C.6 `fluidpy/ch02_cartesian_tensors.py` (chapter module; re-exports all of C.1–C.5)
| # | Callable (signature) | Implements | Used by | Status |
|---|---|---|---|---|
| 6.1 | `polar_components(u1, u2, theta) -> (u_r, u_theta)` · `cartesian_from_polar(u_r, u_theta, theta)` · `polar_components_deg(u1, u2, theta_deg)` | Ex. 2.1 = (2.5) with j ∈ {r, θ}; C = [[cos θ, −sin θ],[sin θ, cos θ]] | C03 (N19), E1 polar mode | §4 F7 |
| 6.2 | `shear_flow_stress(a) -> (2,2)` · `example_2_2(a, phi) -> dict(n, f, magnitude, angle_rad, sigma_n, tau_s, tau_rot)` (angle by arctan2 mod 2π; τ_rot = CᵀτC with C = rotation_matrix_2d(φ)) | Ex. 2.2 (N38) | C05, C06, E2 preset | §4 F11 |
| 6.3 | `traction_2d(tau, phi) -> dict(n, f, sigma_n, tau_s)` (any symmetric or not 2×2 τ; f = n·τ first index) · `stress_vs_angle(tau, phi_array) -> dict(sigma_n, tau_s, f1, f2)` · `principal_angle_2d(S) -> float` (½ arctan2(2S₁₂, S₁₁ − S₂₂)) | (2.15) in 2-D; Mohr; Ex. 2.4 angle | C05 slider, C06 slider, C13, E2 | §8.1, §4 F22 |
| 6.4 | `example_2_4(Gamma) -> dict(S, lam, B, C, S_prime, angle_rad)` (Γ = S₁₂; λ ordered (Γ, −Γ) as the book, b¹ first; C = 45° rotation with det +1; docstring states the ½ du₁/dx₂ factor) | Ex. 2.4 (N62) | C13, E2/E3 presets | §4 F22 (+§8.6) |
| 6.5 | `radial_field(a, dim=3) -> VectorField` · `solid_body_rotation_field(b) -> VectorField` (3-D) · `shear_field(Gamma) -> VectorField` (2-D, u₁ = Γ x₂) · `irrotational_vortex_field(K) -> VectorField` (2-D, u_θ = K/r; `singular_at = (0.0, 0.0)`) · `potential_field(phi_expr=None, dim=2) -> VectorField` (default φ = x₁² − x₂²; parity rows call `ch02.potential_field()`) · `point_source_field(m) -> VectorField` (2-D, `singular_at`) · `smooth_test_field(dim=3) -> VectorField` · `exact_div_curl(field_expr: sympy Matrix) -> (div, curl)` · `example_2_3(a, b, n, L) -> dict` (Ex. 2.3 numerically and symbolically, N55) — a `VectorField` (`core/fields.py`) is callable on `(X, Y)` or `(X, Y, Z)` arrays **matching its `dim`**, returns components stacked on axis 0, and carries `div_fn`, `curl_fn`, `grad_fn`, sympy `exprs` and `singular_at`; 2-D fields need 2-D loops/points (see 5.5) | Ex. 2.3 (N55), E4/E5 presets, convergence tests | C10, C11, C14–C16, E4, E5 | §4 F17, §8.4 |
| 6.6 | `linear_flow_map(G, t) -> ndarray` (= `scipy.linalg.expm(G t)`) · `deform_square(G, t, n_side=10, half_width=1.0, boundary_only=False) -> ndarray (2, N)` (tracers x(t) = e^{Gt} x₀ of a square) · `velocity_gradient_preset(name, Gamma=1.0, dim=2) -> ndarray` (WIP names `"simple_shear"` [[0, Γ],[0, 0]], `"solid_body_rotation"` [[0, −Γ],[Γ, 0]], `"pure_strain"` [[Γ, 0],[0, −Γ]], `"uniaxial_extension"` [[Γ, 0],[0, 0]], `"irrotational_strain"` [[0, Γ],[Γ, 0]]; `dim=3` pads with zeros) · `material_line_angle(G, t, theta0) -> float` | C12 animation, E3 | C12, E3 | §8.2 |
| 6.7 | `operator_convergence(op="gradient"\|"divergence"\|"curl", ns=(8, 16, 32, 64), bc="onesided", L=1.0) -> dict(n, h, err, order)` and `integral_definition_convergence(kind="divergence"\|"curl"\|"gradient", x0=(0.3, −0.2, 0.5), hs=(0.4, 0.2, 0.1, 0.05), field=None, n_face=8) -> dict(h, err, values, exact, order)` — the WIP's two helpers (this design's `convergence_study` is these two); the notebook calls the three operator studies once in C09 and the three integral studies once in C15, stores the dicts and reuses them in C15/C16 (curation §8.11) | V3 evidence made visible | C09, C15, C16 | WIP ✓ (`ch02.operator_convergence`, `ch02.integral_definition_convergence`) |
| 6.8 | `EXAMPLE_TENSOR = [[2, 1, 0],[1, 3, 1],[0, 1, 4]]` (symmetric, invariants 9, 24, 18) · `PRESSURE_SHEAR(p, a) = [[-p, a, 0],[a, -p, 0],[0, 0, -p]]` | the running numbers of C04, C07, C13 | notebook | notebook-local constants defined in the first code cell of C04 / C07 (no library change) |

### C.7 `scripts/` (drawing helpers, no physics; imported as `from scripts.ch02_drawings import …`)
`ch02_fig2_1_position_vector.py` (plotly arrow + dashed projections), `ch02_fig2_2_rotated_axes.py` (two frames, one
arrow, both component sets; `frame_rotation_frames(theta_list)` for the C02 animation), `ch02_fig2_3_polar_components.py`,
`ch02_fig2_4_stress_cube.py` (`stress_cube_figure(tau, show_hidden=False)` → go.Figure with ≤ 18 `go.Cone`),
`ch02_fig2_5_tetrahedron.py` (`tetrahedron_figure(tau, n)`), `ch02_fig2_6_shear_traction.py` (2-D element + f arrow at φ),
`ch02_fig2_7_gradient.py`, `ch02_fig2_8_principal_axes.py` (deforming square + eigen-axes + Mohr circle),
`ch02_fig2_9_gauss.py` (`gauss_box_figure(Q_fn)` faces coloured by n·Q), `ch02_fig2_10_stokes.py` (cap + n, n_c, t),
`ch02_matrix_product_highlight.py` (`matrix_product_figure(i, j)` — (2.11) row/column boxes), `ch02_examples.py`
(Examples 2.1–2.6 numerically and symbolically), `ch02_convergence.py` (wraps 6.7 and plots). Figures → `outputs/ch02/`.

**Count:** C.1 6 · C.2 16 rows (≈ 60 callables) · C.3 2 · C.4 4 rows · C.5 6 rows · C.6 8 rows · C.7 scripts.
**Beyond analysis §4 and curation §8 (NEW):** `expand_indices_str` (1.2 — already in the WIP), `divergence_theorem_tiled` (5.3 —
**the one function still missing from the WIP**); the convergence helpers (6.7) and `flux_through_box_faces`, `example_2_3` exist in the
WIP under those names and are adopted here. Everything else is F1–F26 or curation §8. **WIP cross-check (probe of 2026-09-16):** all
other Part C names exist with compatible signatures; the design was aligned to the WIP's `partial(f, direction, …)`, face keys
`'+x' …`, `cube_face_tractions` keys `'+1' …`, `VectorField` objects with `dim`, `StokesCheck`, and the preset names.

**Contract details the notebook and explainers rely on (please keep):** `traction` never symmetrises; `principal_axes`
returns λ ascending but `example_2_4` returns the book's order (Γ, −Γ); `rotation_matrix_2d(theta)` returns the
*passive* C (so `transform_vector(x, C)` = Cᵀx); `flux_through_faces` returns signed fluxes summing to `rhs`;
`stokes_theorem_check` never asserts when `hypothesis_ok` is False; `grid.h` is a tuple `(hx, hy, hz)` even when equal (WIP ✓: node grids, h = (b − a)/(n − 1));
`gradient` component order is (x, y, z) regardless of the array layout `[k, j, i]`; `expand_indices_str` uses `*` and
` + ` with single spaces (B1 compares the string).

---

## Part A — notebook storyboard (`notebooks/build_ch02.py` → `notebooks/ch02_cartesian_tensors.ipynb`)

**Section check (book order; §2.6 is merged with §2.4 because Cauchy's formula is taught before (2.12) — curation
decision 2; §2.7 and §2.8 are merged because the cross product is ε's first job). Cell numbers are estimates for the
builder's planning (± 3).**
- §2.1 → C01 · N01–N08, N12, N21–N24, N40, N41, N45 — cells 6–36
- §2.2 → C02 (N09, N10, N14, N15, N17 · D02, D03), C03 (N11, N13, N16, N18, N19, N20 · D01 · **E1**) — cells 37–78
- §2.3 → no A item; pointer paragraph tagged → C01, C03 + one code cell — cells 79–81
- §2.4, 2.6 → C04 (N25, R01, N26), C05 (N34–N39 · D05 · **E2**), C06 (N27–N29 · D06) — cells 82–128
- §2.5 → C07 (N30–N33 · D18) — cells 129–146
- §2.7, 2.8 → C08 (N42–N44, N46–N50 · D09) — cells 147–172
- §2.9 → R02, C09 (N51), C10 (N52), C11 (N53–N55 · D12) — cells 173–216
- §2.10 → C12 (N56–N61 · D14, D15 · **E3**) — cells 217–241
- §2.11 → C13 (N62, N63 · D17 ★★★) — cells 242–264
- §2.12 → C14 (N64, N65 · D25 · **E4**), C15 (N66–N68, N72 · D21, D22) — cells 265–302
- §2.13 → C16 (N69–N71, N73 · D26 · **E5**) — cells 303–326
- §2.14 → N74 (tagged → C10) — cells 327–330
- end → S01, S02, summary — cells 331–334

Every CORE block below has the nine parts in order (question · idea · primers · maths with D rows · tiny example · code
+ "What does the code above do?" · from-scratch check · ≥ 1 visual · notes + "what would change"). *expect:* gives the
numbers the executed cell must print (computed by hand in this file; the builder re-checks them against the WIP
`fluidpy`). *explain:* is the numbered "What does the code above do?" list. Markdown drafts are in our words; equation
LaTeX is transcribed from the page images.

### A.0 Front matter
1. `nb.title(big_idea="A vector or a tensor is one physical thing; its components are what a particular set of axes
   sees. This chapter builds the index machinery (a repeated index is a sum), the rotation rule that defines vectors and
   tensors, the two tensors fluid mechanics runs on (stress and velocity gradient), and the two theorems that turn
   'inside a volume' into 'on its boundary' (Gauss and Stokes).", roadmap=["§2.1–2.3 the summation convention, the
   direction-cosine matrix and how components rotate", "§2.4–2.6 the stress tensor, Cauchy's traction f = n·τ, and the
   rule that makes τ a tensor", "§2.5, 2.7, 2.8 contraction and invariants; δ, ε and the cross product", "§2.9 gradient,
   divergence, curl — on a grid you will reuse in every later chapter", "§2.10–2.11 symmetric + antisymmetric parts;
   principal axes", "§2.12–2.14 Gauss and Stokes; comma notation"], prerequisites=["vectors and matrices (Ch. 1
   primers P53–P58 for determinants, rank, null space)", "partial derivatives (Ch. 1 P25) and first-order Taylor
   expansion (P26)", "the definite integral (P27) and the trapezoid rule (P37)", "stress as force per area (Ch. 1 §1.3,
   P05) and the free-body diagram (P09)"])`
2. `nb.explainer_index([("rotation_of_axes", "Which rotates — the arrow or the ruler?", "C02 C03 C06: passive rotation,
   x' = Cᵀx, τ' = CᵀτC"), ("cauchy_traction_principal_axes", "Cut the point any way you like: what pushes on the cut?",
   "C05 C13 C04: f = n·τ, Mohr's circle, principal axes"), ("strain_vs_rotation_split", "Is simple shear a rotation?",
   "C12: G = S + A, the hidden vector of A"), ("gauss_flux_box", "What leaks out of a box?", "C14 C15: Gauss' theorem
   → divergence as outflux per volume"), ("stokes_circulation_loop", "How much does the flow go round a loop?", "C16
   C11: Stokes' theorem → curl as circulation per area")])`
3. `nb.setup()` — identical to ch01; the setup cell also runs `from fluidpy import ch02_cartesian_tensors as ch02`,
   `import numpy as np, sympy as sp, matplotlib.pyplot as plt, plotly.graph_objects as go`.
4. `nb.md` — "**Notation used in this notebook.** Indices $i, j, k, l, m, n$ run over 1, 2, 3 (Python: 0, 1, 2 — we
   say which whenever a printed index appears). $\equiv$ means 'is defined as'. Boldface $\mathbf u$ is the object;
   $u_i$ its components. Colours: teal = the original axes, orange = rotated axes, blue = normal stress, rose = shear."
   (gloss for the ledger rows "index letters and ≡".)
5. `nb.md` — a two-column table "Where this chapter is used later": (2.12) → Ch. 3 strain rate, Ch. 4 stress;
   (2.15) → Ch. 4 §4.3 wall forces, Ch. 8–9 drag; ε → Ch. 4 Coriolis 2Ω × u, Ch. 5 vorticity; (2.30) → Ch. 4 §4.2
   continuity; (2.34) → Ch. 5 Kelvin, Ch. 6/14 lift. Climate hook: the Coriolis term and the geostrophic balance of Ch.
   13 are one ε_ijk and one ∇p away.

### A.1 §2.1 Scalars, Vectors, Tensors, Notation — C01 (+N01–N08, N12, N21–N24, N40, N41, N45)
1. `nb.section("2.1", "Scalars, Vectors, Tensors, Notation", intro="**What is this section about?** Three kinds of
   quantity (one number, three numbers, nine numbers), the index notation that writes all of them as symbols with
   subscripts, and the one convention — a repeated index means 'add over 1, 2, 3' — that lets the whole book drop its
   Σ signs.")`
2. `nb.note` — **N01 [B]** "**A scalar (zero-order tensor)** is one number at each point, the same whatever axes you
   choose: pressure $p$, temperature $T$, density $\rho$ — the fields of Ch. 1. Zero free indices." · **N02 [B]** "**A
   vector (first-order tensor)** has a size and a direction: position $\mathbf x$, velocity $\mathbf u$, gravity
   $\mathbf g$. Its three components $u_1, u_2, u_3$ depend on the axes; the arrow does not. §2.2 sharpens this into a
   test (C03)." · **N06 [B]** "**A second-order tensor** has one component for each *pair* of directions, $3\times3 =
   9$: the stress $\tau_{ij}$ (which face, which force direction — C04). Two free indices." · **N24 [C]** "Hierarchy 1
   → 3 → 9: Newtonian fluid mechanics needs orders 0–2 (one fourth-order exception, the viscosity tensor of Ch. 4 §4.5,
   N28)." (four `nb.note` calls.)
3. `nb.core("C01", "The summation convention: a repeated index is a sum", question="Why does the book never write a
   Σ sign — and how do you read $a_i b_i$, $A_{ik}B_{kj}$ or $C_{im}C_{jn}\tau_{ij}$ without one?")`
4. `nb.md` — **The problem in plain words:** "Every equation of fluid mechanics is really three or nine equations at
   once (one per direction, one per pair). Written out, the momentum equation fills a page; in index notation it is one
   line. The price is a reading rule you must run in your head: **if a letter appears twice in a term, add that term up
   for the letter = 1, 2, 3.** Once the rule is automatic, $\tau_{ji}n_j$ *is* a matrix–vector product and
   $\partial u_i/\partial x_i$ *is* the divergence."
5. `nb.md` — **The idea:**
   ```
   a_i b_i     →  i is repeated  →  a_1 b_1 + a_2 b_2 + a_3 b_3           (one number:   no free index)
   A_ik B_kj   →  k is repeated  →  A_i1 B_1j + A_i2 B_2j + A_i3 B_3j     (nine numbers: free i, j)
   C_im C_jn τ_ij → i, j repeated → 9 terms for each (m, n)               (nine numbers: free m, n)
   ```
   "**Free index = shape of the answer; repeated (dummy) index = a loop.** A dummy letter can be renamed at will
   ($a_i b_i = a_k b_k$); a free letter must match on both sides of an equation."
6. `nb.note` — **N03 [B]** (2.1): "$\mathbf x = \mathbf e_1 x_1 + \mathbf e_2 x_2 + \mathbf e_3 x_3$ — the position
   vector is its components times the unit vectors; the subscript on $\mathbf e$ names an *axis*, not a component."
   equation `\mathbf x = \mathbf e_1 x_1 + \mathbf e_2 x_2 + \mathbf e_3 x_3`, ref "2.1". · **N04 [C]** "Column form
   $\mathbf x = [x_1\ x_2\ x_3]^{\rm T}$, transpose swaps rows and columns; the book also writes triplets $(x_1, x_2,
   x_3)$ — in code `x = np.array([1., 2., 3.])`, `x.T`." (inside the same note.)
7. `nb.primer("np.einsum index strings", "The code form of the summation convention: `np.einsum('i,i', a, b)` sums over
   the repeated letter i; `'ik,kj->ij'` sums over k and keeps i, j; letters after `->` are the free indices (the shape
   of the answer). Nothing else is needed to write any formula of this chapter.", code="import numpy as np\na, b =
   np.array([1., 2., 3.]), np.array([4., 5., 6.])\nprint(np.einsum('i,i', a, b))            # 32.0 = 4 + 10 + 18\nA =
   np.arange(9.).reshape(3, 3)\nprint(np.einsum('ik,kj->ij', A, A)[0, 1])   # P_12 = A_1k A_k2 (Python row 0, col 1)")`
8. `nb.primer("matrix multiplication, transpose and identity", "Row-by-column: $(AB)_{ij} = \sum_k A_{ik}B_{kj}$ — walk
   along row i of A and down column j of B, multiply pairs, add. In Python `A @ B`; `A.T` swaps rows and columns
   ($A^{\rm T}_{ij} = A_{ji}$); `np.eye(3)` is the identity (1 on the diagonal, 0 elsewhere) — it leaves any matrix
   unchanged.", code="A = np.array([[1., 2.], [3., 4.]]); B = np.array([[0., 1.], [1., 0.]])\nprint(A @ B)          #
   [[2, 1], [4, 3]]: row i of A · column j of B\nprint(A.T)            # [[1, 3], [2, 4]]\nprint(np.eye(2) @ A)  # A
   again")`
9. `nb.md` — **The maths, step by step (three equations expanded by hand):** "1. Dot product (2.2): $\mathbf a\cdot
   \mathbf b = a_1b_1 + a_2b_2 + a_3b_3 = \sum_{i=1}^{3} a_ib_i \equiv a_ib_i$ — one repeated letter, one sum, a scalar.
   2. Matrix product (2.9): $P_{ij} = \sum_{k=1}^{3} A_{ik}B_{kj} \equiv A_{ik}B_{kj}$ — k summed, i and j free, nine
   numbers. 3. Kronecker substitution (2.17): $\delta_{ij}u_j = \delta_{i1}u_1 + \delta_{i2}u_2 + \delta_{i3}u_3 = u_i$
   — only the term with j = i survives, so δ *replaces* its summed index."
10. `nb.note` — **N07 [B]** (2.2) with equation `\mathbf a\cdot\mathbf b = a_1 b_1 + a_2 b_2 + a_3 b_3 = \sum_{i=1}^{3}
    a_i b_i \equiv a_i b_i`, ref "2.2": "the ≡ is the summation convention at work." · **N45 [C]** "Also $\mathbf u\cdot
    \mathbf v = uv\cos\theta$ (the projection of one on the other) and $\mathbf u\cdot\mathbf v$ is the trace of the
    nine products $u_iv_j$ — the contraction idea of §2.5 (C07); `angle_between(u, v)`."
11. `nb.note` — **N40 [B]** (2.16): "**The Kronecker delta** $\delta_{ij} = 1$ if $i = j$, 0 otherwise — the identity
    matrix in index clothes (`kronecker_delta()` = `np.eye(3)`); we need it now because the rotation matrix's
    orthogonality (C02) is written with it." equation `\delta_{ij} = \begin{cases} 1 & i = j \\ 0 & i \neq j\end{cases}`,
    ref "2.16". · **N41 [B]** (2.17): "**δ substitutes its summed index**: $\delta_{ij}u_j = u_i$, $\delta_{ij}A_{jk} =
    A_{ik}$, and — the classic slip — $\delta_{ii} = 3$ (a summed pair), not 1." equation `\delta_{ij} u_j = u_i`, ref
    "2.17".
12. `nb.note` — **N21 [B]** (2.9): "**Matrix product as an index sum**: $P_{ij} = A_{ik}B_{kj}$ sums the *adjacent*
    index k; `inner(A, B)` is `np.einsum('ik,kj->ij')` and equals `A @ B`." equation `P_{ij} = \sum_{k=1}^{3} A_{ik}
    B_{kj} \equiv A_{ik} B_{kj}`, ref "2.9". · **N22 [C]** "(2.10) $\mathbf P = \mathbf A\cdot\mathbf B$: a single dot
    means one index summed — the same dot that later reads $\mathbf n\cdot\boldsymbol\tau$ (C05) and
    $\nabla\cdot\mathbf u$ (C10)." · **N23 [B]** (2.11): "Explicitly, $P_{12} = A_{11}B_{12} + A_{12}B_{22} + A_{13}
    B_{32}$: row 1 of A against column 2 of B — the figure below boxes them."
13. `nb.note` — **N08 [B]** "**Boldface for meaning, indices for manipulation.** $\mathbf a\mathbf b$ is ambiguous for
    tensors ($\mathbf a\mathbf b \neq \mathbf b\mathbf a$) and hides the order; indices never do: **order = number of
    free indices**. `tensor_order('A_ij B_kl')` → 4." · **N12 [B]** "**Free vs dummy**: $x_iC_{ij} = x_kC_{kj}$
    (rename the summed letter freely); the free j must appear on both sides. A letter appearing three times is an
    error." (§2.2 material, taught here.)
14. `nb.worked_example("(1, 2, 3)·(4, 5, 6) and δ_ij u_j for i = 2", "1. $a_ib_i$ with i = 1, 2, 3: $1\times4 +
    2\times5 + 3\times6 = 4 + 10 + 18 = 32$. 2. $\delta_{2j}u_j = \delta_{21}u_1 + \delta_{22}u_2 + \delta_{23}u_3 = 0 +
    u_2 + 0 = u_2$. 3. $P_{12}$ for $A = B = \begin{bmatrix}1&2&3\\4&5&6\\7&8&9\end{bmatrix}$: $1\cdot2 + 2\cdot5 +
    3\cdot8 = 2 + 10 + 24 = 36$.")`
15. `nb.code[explain]` — *code:* `a, b = np.array([1., 2., 3.]), np.array([4., 5., 6.])`; `print(ch02.expand_indices_str
    ("a_i b_i"))` · `print(ch02.dot(a, b))` · `print(ch02.expand_indices_str("delta_ij u_j"))` · `A = np.arange(1.,
    10.).reshape(3, 3)`; `print(ch02.inner(A, A)[0, 1], (A @ A)[0, 1])` · `print(ch02.classify_indices("x_i C_ij"),
    ch02.tensor_order("A_ij B_kl"))` · `print(ch02.rename_dummy("x_i C_ij", "i", "k"))`. *expect:* `a_1*b_1 + a_2*b_2 +
    a_3*b_3` · 32.0 · `u_i` (or the three-line array u_1, u_2, u_3) · 36.0 36.0 · (['j'], ['i']) 4 · `x_k C_kj`.
    *explain:* 1. the expander parses subscripts, finds the repeated letter and writes the hidden sum; 2. `dot` is
    `np.einsum('i,i')`; 3. δ collapses the sum; 4. `inner` is (2.9) and equals `@`; 5. free/dummy classification and the
    order; 6. renaming a dummy changes nothing.
16. `nb.check_agree` — **from scratch (curation §7):** `mine = 0.0; for i in range(3): mine += a[i]*b[i]` → `assert
    np.allclose(mine, ch02.dot(a, b))`; `P = np.zeros((3, 3)); for i in range(3): for j in range(3): for k in range(3):
    P[i, j] += A[i, k]*A[k, j]` → `assert np.allclose(P, ch02.inner(A, A))`; `assert np.allclose(P, np.einsum('ik,kj->ij',
    A, A))`. Markdown: "The triple loop *is* the summation convention; einsum only writes it faster."
17. `nb.figure` — (2.11) row/column highlight: `from scripts.ch02_matrix_product_highlight import
    matrix_product_figure; fig = matrix_product_figure(i=1, j=2, A=A, B=A)` (row 1 of A teal, column 2 of B orange, the
    three products and their sum 36 written under the cell P₁₂). *see:* "three boxed pairs, one sum." *read:* "the
    highlighted row and column are the two copies of the dummy k walking together." *change:* "pick i = 3, j = 1: row
    3 and column 1 light up, P₃₁ = 7·1 + 8·4 + 9·7 = 102."
18a. `nb.primer("plotly 3-D arrows, lines and meshes", "Ch. 1 used `go.Surface` and `go.Scatter3d` points. Arrows are
    `go.Cone(x, y, z, u, v, w)` (position and direction), edges are `go.Scatter3d(mode='lines')`, faces are `go.Mesh3d`
    (vertices + triangles). All stay rotatable on the published page.", code="import plotly.graph_objects as go\nfig =
    go.Figure(go.Cone(x=[0], y=[0], z=[0], u=[1], v=[0], w=[0], sizemode='absolute', sizeref=0.5))\nfig.update_layout
    (height=250, scene_aspectmode='cube'); fig.show()   # one arrow along x")` — placed here because the next cell is the
    chapter's first 3-D figure (also used by C02, C04, C05, C14, C16).
18. `nb.plotly` — **N05 [B]** Fig. 2.1 as our own rotatable 3-D figure: `from scripts.ch02_fig2_1_position_vector import
    position_vector_figure; position_vector_figure(x=[1, 2, 3])` — arrow OP, dashed projections onto the axes and the
    x₁x₂ plane, e₁ e₂ e₃ as short unit arrows. *see/read/change* in a following `nb.md`: "rotate the view: the numbers
    (1, 2, 3) are the shadows of one arrow on three rulers — §2.2 changes the rulers."
19. `nb.md` — **What would change if…** "…you meet $\varepsilon_{ijk}u_iv_j$ (§2.7)? The same rule: i and j are summed
    (nine terms), k is free — a vector. And $\partial u_i/\partial x_i$ (§2.9) is a sum of three derivatives, a scalar.
    Every later formula is read this way; try `expand_indices_str` on it when in doubt."

### A.2 §2.2 Rotation of Axes: Formal Definition of a Vector — C02 (+N09, N10, N14, N15, N17 · D02, D03), C03 (+N11, N13, N16, N18, N19, N20 · D01 · E1)
1. `nb.section("2.2", "Rotation of Axes: Formal Definition of a Vector", intro="**What is this section about?** Keep
   the arrow, turn the rulers. The nine cosines between old and new axes form a matrix C; the components of *any*
   vector change by the same rule x' = Cᵀx — and that rule becomes the definition of a vector.")`
2. `nb.core("C02", "The direction-cosine matrix C_ij = e_i · e'_j", question="Two sets of axes share an origin. What is
   the one table of numbers that relates them — and what does a single entry $C_{ij}$ measure?")`
3. `nb.md` — **The problem in plain words:** "A weather model uses axes east/north/up; a wind sensor on a tilted mast
   reports its own three components; a stress calculation wants axes normal and tangent to a wall. Same wind, same
   stress, different numbers. We need the dictionary between two frames — and it turns out to be nine cosines."
4. `nb.md` — **The idea:**
   ```
   old axes  e_1, e_2, e_3  (teal)        new axes  e'_1, e'_2, e'_3  (orange), same origin O
   C_ij = e_i · e'_j = cos(angle between old axis i and new axis j)
   row i  = "how much of old axis i lies along each new axis"
   column j = "the new axis e'_j written in old components"   ← the whole matrix, one picture
   ```
   "For a turn by θ about $\mathbf e_3$: $\mathbf e'_1 = (\cos\theta, \sin\theta, 0)$, $\mathbf e'_2 = (-\sin\theta,
   \cos\theta, 0)$, $\mathbf e'_3 = \mathbf e_3$ — these are the columns of C."
5. `nb.note` — **N09 [B]** "**The rotated frame** O1'2'3' shares the origin; its unit vectors are $\mathbf e'_j$; the
   vector $\mathbf x$ is *the same arrow* with new components $x'_j$." · **N17 [C]** "Fig. 2.2 is our own 3-D drawing
   below (two frames, one arrow) and the animation after it." · **N10 [B]** (2.3): "$\mathbf x = x'_1\mathbf e'_1 +
   x'_2\mathbf e'_2 + x'_3\mathbf e'_3$ — the same vector spelled in the primed basis; in code
   `vector_from_components(xp, E_new)` rebuilds the very same array as `vector_from_components(x, E_old)`." equation
   `\mathbf x = x'_1 \mathbf e'_1 + x'_2 \mathbf e'_2 + x'_3 \mathbf e'_3`, ref "2.3".
6. `nb.primer("orthonormal basis, projection and completeness", "A basis is orthonormal when each vector has length 1
   and any two are perpendicular: $\mathbf e_i\cdot\mathbf e_j = \delta_{ij}$. The component of a vector along a unit
   vector is the dot product (its projection). **Completeness**: a vector is the sum of its projections, $\mathbf v =
   \sum_j (\mathbf v\cdot\mathbf e'_j)\,\mathbf e'_j$ — nothing is lost because three perpendicular unit vectors span
   3-D space.", code="E = ch02.rotation_matrix_3d([0, 0, 1], np.deg2rad(30)).T   # rows = the new unit vectors\nprint
   (np.round(E @ E.T, 12))                       # identity: orthonormal\nv = np.array([1., 2., 3.])\nprint(sum((v @
   E[j]) * E[j] for j in range(3)))          # [1 2 3]: the sum of the projections is v")`
7. `nb.primer("cosines of angles between unit vectors", "For unit vectors the dot product *is* the cosine of the angle
   between them. Two identities we use: $\cos(\pi/2 - \theta) = \sin\theta$ and $\cos(\theta + \pi/2) = -\sin\theta$.
   Python trig works in radians: `np.deg2rad(30)` = 0.5236.", code="th = np.deg2rad(30)\nprint(np.cos(np.pi/2 - th),
   np.sin(th))       # 0.5 0.5\nprint(np.cos(th + np.pi/2), -np.sin(th))     # -0.5 -0.5")`
8. `nb.primer("np.linalg.norm and np.linalg.qr", "`np.linalg.norm(v)` is the length $\sqrt{v_iv_i}$. `np.linalg.qr(M)`
   factors a matrix into an orthogonal Q and a triangular R — we use Q of a random matrix as a random rotation (after
   fixing det Q = +1).", code="rng = np.random.default_rng(0)\nQ, _ = np.linalg.qr(rng.normal(size=(3, 3)))\nprint
   (np.round(Q.T @ Q, 12))                     # identity: Q is orthogonal\nprint(np.linalg.norm(Q @ [3., 4., 0.]))
   # 5.0: a rotation keeps lengths")`
9. `nb.md` — **The maths, step by step:** "1. Define $C_{ij} \equiv \mathbf e_i\cdot\mathbf e'_j$ (nine numbers;
   dimensionless). 2. Column j of C is $\mathbf e'_j$ in old components: $(\mathbf e'_j)_i = \mathbf e_i\cdot\mathbf e'_j
   = C_{ij}$. 3. Row i of C is $\mathbf e_i$ in new components. 4. Because both bases are orthonormal, C is
   **orthogonal**: $C_{ij}C_{kj} = \delta_{ik}$ and $C_{ji}C_{jk} = \delta_{ik}$ — matrix form $\mathbf C\mathbf C^{\rm T}
   = \mathbf C^{\rm T}\mathbf C = \mathbf I$ (D02 below). 5. Hence $\mathbf C^{-1} = \mathbf C^{\rm T}$ and $\det\mathbf C
   = +1$ for a rotation (−1 would be a mirror)."
10. `nb.derivation("D02", "Why C is orthogonal: C Cᵀ = Cᵀ C = I and det C = +1", ref="Exercise 2.8", …)` — Part F · D02
    (8 steps).
11. `nb.note` — **N15 [B]** "**Orthogonality of C** — the result of D02 — with the code check `is_orthogonal(C)`,
    `orthogonality_residual(C)` ≈ 1e-16 and `is_proper_rotation(C)` (det +1)." equation `C_{ij} C_{kj} = C_{ji} C_{jk}
    = \delta_{ik}, \qquad \det \mathbf C = +1`, ref "Exercise 2.8".
12. `nb.derivation("D03", "The inverse transformation x_j = x'_i C_ji", ref="2.7", …)` — Part F · D03 (4 steps). (It
    uses (2.5), stated in the *start* line and derived in C03's D01 just below — the notebook says "(2.5) is derived in
    the next block; here we only need that it holds".)
13. `nb.note` — **N14 [B]** (2.7): "**Back to the old components**: $x_j = x'_iC_{ji}$ — now the *second* index of C is
    summed (matrix form $\mathbf x = \mathbf C\mathbf x'$). Round trip in code: `inverse_transform_vector
    (transform_vector(x, C), C)` returns x to 1e-16." equation `x_j = \sum_{i=1}^{3} x'_i C_{ji} \equiv x'_i C_{ji}`,
    ref "2.7".
14. `nb.worked_example("a 30° turn about e₃", "1. $\cos 30° = 0.866$, $\sin 30° = 0.5$. 2. New axes in old components:
    $\mathbf e'_1 = (0.866, 0.5, 0)$, $\mathbf e'_2 = (-0.5, 0.866, 0)$, $\mathbf e'_3 = (0, 0, 1)$. 3. Columns → $C =
    \begin{bmatrix}0.866 & -0.5 & 0\\ 0.5 & 0.866 & 0\\ 0 & 0 & 1\end{bmatrix}$. 4. Check one entry of $C^{\rm T}C$:
    row 1 · row 1 of Cᵀ = $0.866^2 + 0.5^2 = 0.75 + 0.25 = 1$; row 1 · row 2 = $0.866\cdot(-0.5) + 0.5\cdot0.866 = 0$
    ✓. 5. $\det C = \cos^2 30° + \sin^2 30° = 1$ ✓ (a rotation, not a mirror).")`
15. `nb.code[explain]` — *code:* `C = ch02.rotation_matrix_3d([0, 0, 1], np.deg2rad(30))`; `print(np.round(C, 3))`;
    `E_old = ch02.unit_vectors()`; `E_new = C.T` (rows = new axes); `C2 = ch02.direction_cosines(E_old, E_new)`;
    `assert np.allclose(C, C2)`; `print(ch02.orthogonality_residual(C), ch02.is_proper_rotation(C))`; `rng =
    np.random.default_rng(0); R = ch02.random_rotation(rng); print(ch02.is_orthogonal(R), np.round(np.linalg.det(R),
    12))`. *expect:* the 3×3 above · 0.0 (≤ 2e-16) True · True 1.0. *explain:* 1. Rodrigues' formula builds the passive
    C (columns = new axes) for a turn by 30° about e₃; 2. the same C from nine dot products; 3. orthogonality residual
    and det; 4. a random rotation for later invariance tests.
16. `nb.check_agree` — **from scratch:** `C_mine = np.zeros((3, 3)); for i in range(3): for j in range(3): C_mine[i, j]
    = E_old[i] @ E_new[j]` → `assert np.allclose(C_mine, ch02.direction_cosines(E_old, E_new))`; `assert np.allclose
    (C_mine[:2, :2], ch02.rotation_matrix_2d(np.deg2rad(30)))`.
17. `nb.md` — ⚠️ **Common confusion (passive vs active):** "Wikipedia's rotation matrix $R(\theta)$ rotates a *vector*
    counterclockwise; it has exactly the entries of our C for a frame turned by +θ. But the book uses C *passively*:
    the arrow stays, the axes turn, so the new components are $\mathbf x' = \mathbf C^{\rm T}\mathbf x$ — the transpose.
    Turning the axes by +θ looks, from the axes' point of view, like turning the arrow by −θ. `scipy.spatial.transform
    .Rotation.as_matrix()` is active too. Code: every C-returning function says 'passive'; E1's toggle shows both."
18. `nb.plotly` — **N09/N17** 3-D figure: `from scripts.ch02_fig2_2_rotated_axes import rotated_axes_figure;
    rotated_axes_figure(C, x=[1, 2, 3])` — teal old frame, orange new frame, one arrow, dashed projections in both
    frames, the two component sets in the legend. Following `nb.md` see/read/change: "the arrow is drawn once; two
    sets of dashed shadows fall from it."
19. `nb.animation` — curation §6 row 1 (`player="frames"`, 10 frames at 0°, 10°, … 90°; FAST 6 frames): left panel the
    x₁x₂ plane with the fixed arrow (1, 2), teal axes fixed, orange axes turning; right panel two bar groups, $x_i$
    (teal, constant) and $x'_j$ (orange, moving); a text box with the live 2×2 C and "CᵀC − I = 0". *code intent:*
    `ch02.rotation_matrix_2d(th)`, `ch02.transform_vector([1, 2], C2)`. *see:* "the arrow never moves, the orange
    rulers do; the orange bars change, the teal ones do not." *read:* "at 90° the roles swap: $x'_1 = x_2$, $x'_2 =
    -x_1$." *change:* "turn by −30° instead: C becomes its transpose."
20. `nb.md` — **What would change if…** "…the second frame were *mirrored* ($\mathbf e'_3 = -\mathbf e_3$)? C is still
    orthogonal but det C = −1: not a rotation. §2.7 will show that ε changes sign under such a frame — the book's
    'rotation' always means det C = +1."
21. `nb.core("C03", "How components transform: x'_j = x_i C_ij — the formal definition of a vector", question="If the
    arrow does not change but its numbers do, exactly how do the numbers change — and which triples of numbers deserve
    the name 'vector'?")`
22. `nb.md` — **The problem in plain words:** "Any three numbers can be stacked into a column. Are (temperature,
    pressure, density) a vector? Are the three squares $(x_1^2, x_2^2, x_3^2)$? Physics says no: turn the axes and
    those triples do not turn like an arrow. The book's answer is operational — **a vector is anything whose components
    change under a rotation exactly as the position vector's do.**"
23. `nb.md` — **The idea:** "The same arrow casts two shadows: onto the old axes and onto the new. Each new shadow is a
    weighted sum of the old ones, the weights being cosines: $x'_j = x_1C_{1j} + x_2C_{2j} + x_3C_{3j}$. Projection does
    the work: dot $\mathbf x$ with $\mathbf e'_j$ in both spellings of $\mathbf x$."
24. `nb.note` — **N11 [B]** (2.4): "**The first projection**: $\mathbf x\cdot\mathbf e'_1 = x_1\,\mathbf e_1\cdot\mathbf
    e'_1 + x_2\,\mathbf e_2\cdot\mathbf e'_1 + x_3\,\mathbf e_3\cdot\mathbf e'_1 = x'_1$. For θ = 30° and x = (1, 2, 0):
    $1\times0.866 + 2\times0.5 + 0 = 1.866$." equation `\mathbf x\cdot\mathbf e'_1 = x_1\,\mathbf e_1\cdot\mathbf e'_1 +
    x_2\,\mathbf e_2\cdot\mathbf e'_1 + x_3\,\mathbf e_3\cdot\mathbf e'_1 = x'_1`, ref "2.4".
25. `nb.derivation("D01", "The transformation rule for components", ref="2.5", …)` — Part F · D01 (6 steps).
26. `nb.note` — **N13 [C]** (2.6): "(2.5) with other letters, $x'_i = x_kC_{ki}$, says the same nine things —
    `expand_indices` of both prints identical sums." · **N16 [B]** "**Matrix forms**: $\mathbf x' = \mathbf C^{\rm T}
    \mathbf x$ (2.5: the first index of C summed = rows of Cᵀ) and $\mathbf x = \mathbf C\mathbf x'$ (2.7). In code
    `transform_vector(x, C)` is `C.T @ x` — asserted below." equation `\mathbf x' = \mathbf C^{\rm T}\cdot\mathbf x,
    \qquad \mathbf x = \mathbf C\cdot\mathbf x'`, ref "§2.3". · **N18 [B]** (2.8): "**The formal definition**: $\mathbf
    u$ is a Cartesian vector if and only if $u'_j = u_iC_{ij}$ for every rotation. Pass: $\mathbf x$, $a\mathbf x$,
    $\mathbf b\times\mathbf x$, $\nabla\phi$. Fail: $(x_1^2, x_2^2, x_3^2)$ and $(|\mathbf x|, 0, 0)$ — the residual table
    below shows the failures are not small." equation `u'_j = \sum_{i=1}^{3} u_i C_{ij} \equiv u_i C_{ij}`, ref "2.8".
27. `nb.worked_example("x = (1, 2) seen from axes turned by 30°", "1. $C = \begin{bmatrix}0.866 & -0.5\\ 0.5 &
    0.866\end{bmatrix}$ (columns = new axes). 2. $x'_1 = x_1C_{11} + x_2C_{21} = 1\times0.866 + 2\times0.5 = 1.866$.
    3. $x'_2 = x_1C_{12} + x_2C_{22} = 1\times(-0.5) + 2\times0.866 = 1.232$. 4. Length check: $1^2 + 2^2 = 5$ and
    $1.866^2 + 1.232^2 = 3.482 + 1.518 = 5.000$ ✓ — a rotation cannot change a length. 5. Back with (2.7): $x_1 =
    x'_1C_{11} + x'_2C_{12} = 1.866\times0.866 + 1.232\times(-0.5) = 1.616 - 0.616 = 1.000$ ✓.")`
28. `nb.note` — **N19 [B]** Ex. 2.1 (stated, curation §4c D04): "**Polar components** are (2.5) with j ∈ {r, θ}: $u_r =
    u_1\cos\theta + u_2\sin\theta$, $u_\theta = -u_1\sin\theta + u_2\cos\theta$, i.e. $C = \begin{bmatrix}\cos\theta &
    -\sin\theta\\ \sin\theta & \cos\theta\end{bmatrix}$ — the 2-D matrix of the tiny example. With u = (1, 2), θ = 30°:
    $u_r = 1.866$, $u_\theta = 1.232$; the point's polar frame is our rotated frame." · **N20 [C]** "Fig. 2.3 is the
    'polar' mode of the explainer below (our drawing)."
29. `nb.code[explain]` — *code:* `C2 = ch02.rotation_matrix_2d(np.deg2rad(30))`; `x = np.array([1., 2.])`; `xp =
    ch02.transform_vector(x, C2)`; `print(np.round(xp, 3), np.round(C2.T @ x, 3))`; `print(np.round(ch02.
    inverse_transform_vector(xp, C2), 12))`; `print(np.round(ch02.polar_components(1., 2., np.deg2rad(30)), 3))`;
    `print(np.linalg.norm(x), np.linalg.norm(xp))`. *expect:* [1.866 1.232] [1.866 1.232] · [1. 2.] · (1.866, 1.232) ·
    2.236 2.236. *explain:* 1. the passive 2-D C; 2. (2.5) as `transform_vector` and as `C.T @ x` — identical; 3. (2.7)
    round trip; 4. Ex. 2.1 is the same call; 5. lengths agree.
30. `nb.check_agree` — **from scratch:** `xp_mine = np.array([sum(x[i]*C2[i, j] for i in range(2)) for j in range(2)])`
    → `assert np.allclose(xp_mine, ch02.transform_vector(x, C2))`.
31. `nb.code[explain]` — the **N18 residual table**: `rng = np.random.default_rng(1); R = ch02.random_rotation(rng);
    pts = rng.normal(size=(20, 3))`; candidates as a dict of lambdas: `"x": lambda x: x`, `"3x": lambda x: 3*x`,
    `"b×x": lambda x: np.cross([0, 0, 1], x)`, `"x_i²": lambda x: x**2`, `"(|x|,0,0)": lambda x: np.array([np.linalg.
    norm(x), 0, 0])`; `for name, fn in ...: print(f"{name:10s} residual {ch02.transforms_as_vector(fn, R, pts):.2e}")`.
    *expect:* first three ≤ 1e-14; `x_i²` and `(|x|,0,0)` of order 1 (e.g. 2.7, 1.9). *explain:* 1. one random rotation
    and 20 random points; 2. for each candidate, evaluate the triple in the old frame and transform it with (2.8), and
    evaluate the same *formula* in the new frame's coordinates; 3. residual = largest difference — zero only for true
    vectors. `nb.figure` — bar chart of the five residuals on a log axis (floor 1e-16). *see:* "three bars at the floor,
    two at order one." *read:* "a vector's formula gives the same arrow whichever axes you compute it in." *change:*
    "replace `x**2` by `x**3`: fails too; `x * (x @ x)` passes (a scalar times a vector)."
32. `nb.explainer("rotation_of_axes", heading="Which rotates — the arrow or the ruler?", why="One object, two sets of
    numbers: only motion shows that the arrow stays put while C, x'_j and (later) τ'_mn move together — and the toggle
    settles the passive/active confusion once.", tries=["Drag θ to 90° and read C: the columns are the new axes (0,1)
    and (−1,0).", "Flip 'what rotates' to the arrow: the same C now moves the vector the other way (x' = Cx).", "Polar
    mode: move the point round its circle and watch u_r, u_θ — Ex. 2.1 live.", "Tensor mode, pure shear: find the angle
    where τ'₁₂ = 0 (45°) — that is §2.11's principal axis."])`
33. `nb.md` — **What would change if…** "…you applied the rule *twice*, once per index? You would get the transformation
    of a 3×3 array — and that is precisely the definition of a second-order tensor, (2.12), which we reach in §2.4 after
    meeting the stress tensor and Cauchy's formula."

### A.3 §2.3 Multiplication of Matrices — no A item (paragraphs tagged → C01, C03)
1. `nb.section("2.3", "Multiplication of Matrices", intro="**What is this section about?** The matrix product is the
   summation convention with one adjacent index — we already used it: (2.9)–(2.11) in C01 and the matrix forms x' =
   Cᵀx, x = Cx' in C03. This short section only collects them.")`
2. `nb.md` — "📝 (2.9) $P_{ij} = A_{ik}B_{kj}$ and its boxes (2.11) were expanded in **C01** (N21–N23); the rule 'a
   single dot sums one index' is (2.10) (N22). (2.6) written with the summed index adjacent, $x'_i = C^{\rm T}_{ik}x_k$,
   gives the matrix form of **C03** (N16): $\mathbf x' = \mathbf C^{\rm T}\cdot\mathbf x$, and (2.7) is $\mathbf x =
   \mathbf C\cdot\mathbf x'$."
3. `nb.code[explain]` — `C = ch02.rotation_matrix_3d([0, 0, 1], np.deg2rad(30)); x = np.array([1., 2., 3.])`;
   `assert np.allclose(ch02.transform_vector(x, C), C.T @ x)`; `assert np.allclose(ch02.inner(C.T, C), np.eye(3))`;
   `print(ch02.expand_indices_str("C_ki x_k"))` (the (2.6) form: `C_1i*x_1 + C_2i*x_2 + C_3i*x_3`). *explain:* 1–2.
   the two matrix forms of C03; 3. CᵀC = I as `inner`; 4. the adjacent-index form.

### A.4 §2.4, 2.6 Second-Order Tensors and the Force on a Surface — C04 (+N25, R01, N26), C05 (+N34–N39 · D05 · E2), C06 (+N27–N29 · D06)
1. `nb.section("2.4, 2.6", "Second-Order Tensors and the Force on a Surface", intro="**What is this section about?**
   The stress tensor: nine numbers that tell the force per area on *any* plane through a point. We meet its sign
   convention (C04), derive Cauchy's formula f = n·τ from a tiny tetrahedron (C05, the book's §2.6), and only then
   prove the rule that makes τ a tensor, (2.12) (C06). The book states (2.12) first and cites a tetrahedron argument
   it never shows; deriving (2.15) first lets us *prove* (2.12) in eight moves.")`
2. `nb.core("C04", "The stress tensor τ_ij and its sign convention", question="Nine stress components sit on a cube.
   Which index is the face, which is the force — and which way do positive arrows point?")`
3. `nb.md` — **The problem in plain words:** "Ch. 1 split stress into 'normal' (push/pull) and 'shear' (slide) on one
   surface. But at a point in a fluid there are infinitely many surfaces. Pick three perpendicular ones — the faces of
   a tiny cube — and write down the force per area on each: three faces × three force directions = nine numbers. That
   array is the stress tensor. Ch. 4 will balance forces on this cube to get the equation of motion, so every sign
   here matters."
4. `nb.md` — **The idea (rule table):**
   ```
   τ_ij :  i = the face (its outward normal is ±e_i)      j = the direction of the force per area
   face with outward normal +e_i : positive τ_ij points along +e_j   (tensile normal stress is positive)
   face with outward normal −e_i : the arrows reverse (traction −τ_i·)   ← Newton III at a point (R01)
   i = j : normal stress (τ_11, τ_22, τ_33)      i ≠ j : shear stress (τ_12 = force along 2 on the face ⊥ 1)
   ```
5. `nb.recap("R01", "Opposite faces carry equal and opposite stresses", "Ch. 1 §1.7 balanced pressure on the two faces
   of a box (D05, P09): at a point the two sides of one plane push on each other with equal and opposite force — Newton's
   third law. So the −e₂ face of Fig. 2.4 carries exactly −(τ₂₁, τ₂₂, τ₂₃). Across a *finite* cube the values differ by
   a Taylor step, and that difference is Ch. 4's net force.", where="Ch. 1 §1.7 C20")`
6. `nb.md` — one-line reminder: "3-D arrows, edges and faces are `go.Cone`, `go.Scatter3d(mode='lines')` and
   `go.Mesh3d` — primed in C01 (A.1 #18a)."
7. `nb.md` — **The maths, step by step:** "1. On the face whose outward normal is $+\mathbf e_2$ (ABCD in Fig. 2.4),
   the force per area is the vector $(\tau_{21}, \tau_{22}, \tau_{23})$ [Pa]: the first index names the face, the
   second the component. 2. On the opposite face ($-\mathbf e_2$) it is $-(\tau_{21}, \tau_{22}, \tau_{23})$ (R01).
   3. Collect the three '+' faces as rows: $\boldsymbol\tau = \begin{bmatrix}\tau_{11}&\tau_{12}&\tau_{13}\\ \tau_{21}
   &\tau_{22}&\tau_{23}\\ \tau_{31}&\tau_{32}&\tau_{33}\end{bmatrix}$ — **row i = traction on the +e_i face**. 4.
   Pressure will enter Ch. 4 as $\tau_{ij} = -p\,\delta_{ij} + (\text{viscous})$: a pressure pushes *inward*, so its
   normal stress is negative."
8. `nb.note` — **N26 [C]** "The 3×3 array above is what the code carries everywhere (`tau[i, j]`, Python indices 0–2);
   it is a *tensor* only because it obeys (2.12) — C06." · **N25 [B]** "Fig. 2.4 is our rotatable cube below: nine
   arrows on the three visible faces, a toggle for the hidden faces (reversed arrows)."
9. `nb.worked_example("a pressure p with one shear a", "τ = [[−p, a, 0],[a, −p, 0],[0, 0, −p]] with p = 3 Pa, a = 1 Pa.
   1. Face +e₁: traction row 1 = (−3, 1, 0) Pa — pushed *into* the cube by 3 Pa (compression, negative normal stress)
   and dragged along +e₂ by 1 Pa. 2. Face −e₁: (+3, −1, 0) Pa (reversed). 3. Face +e₃: (0, 0, −3): pure pressure, no
   shear. 4. τ₁₂ = τ₂₁ = 1 Pa: the shear on the face ⊥ 1 along 2 equals the shear on the face ⊥ 2 along 1 — Ch. 4 shows
   this symmetry is forced by angular momentum; here we simply chose it.")`
10. `nb.code[explain]` — *code:* `p, a = 3.0, 1.0; tau = np.array([[-p, a, 0], [a, -p, 0], [0, 0, -p]])` (the notebook-local `PRESSURE_SHEAR`); `faces = ch02.cube_face_tractions
    (tau)`; `for name, (n, f) in faces.items(): print(f"{name}: n = {n}, traction = {f} Pa")`; `total = sum(f for _, f
    in faces.values()); print("sum over the six faces:", total)`; `print(ch02.stress_component_meaning(2, 3))`. *expect:*
    `+1: n = [1 0 0], traction = [-3. 1. 0.]`, `-1: … [3. -1. -0.]`, … (WIP keys `'+1'` … `'-3'`); sum `[0. 0. 0.]`; "τ_23: on the face whose
    outward normal is +x2, the force per unit area along x3 …". *explain:* 1. the running stress state; 2. six faces →
    (outward normal, traction) by the sign rule (row i, sign of the normal); 3. the six tractions cancel — stresses *at a
    point* carry no net force (R01); 4. the sentence the rule table says, generated.
11. `nb.check_agree` — **from scratch:** `mine = {}; for i in range(3): e = np.eye(3)[i]; mine[f"+{i+1}"] = (e, tau[i]);
    mine[f"-{i+1}"] = (-e, -tau[i])` → `for k in mine: assert np.allclose(mine[k][1], faces[k][1])`.
12. `nb.plotly` — `from scripts.ch02_fig2_4_stress_cube import stress_cube_figure; stress_cube_figure(tau,
    show_hidden=False)` (9 cones; normal stresses blue, shears rose; face letters A–H as in the book's naming since Ex.
    2.5 uses them). `nb.md` see/read/change: *see:* "on each visible face one arrow points straight out or in (normal)
    and two lie in the face (shear)." *read:* "the face's index is the first subscript, the arrow's direction the
    second; negative normal stress points *into* the cube." *change:* "set a = 0: only the three inward pressure
    arrows remain, equal on every face — an isotropic stress, Ch. 1's pressure."
13. `nb.md` — **What would change if…** "…you cut the point with a *slanted* plane? None of the nine arrows sits on it —
    yet the force per area there is fixed by them. That is Cauchy's formula, next."
14. `nb.core("C05", "Cauchy's traction formula f = n·τ", question="Given the nine stresses on the coordinate faces,
    what is the force per unit area on a plane with an arbitrary unit normal $\mathbf n$?")`
15. `nb.md` — **The problem in plain words:** "The drag on a wing, the friction on a river bed, the pressure on a dam:
    each is the force per area on a surface whose normal points wherever the surface happens to face. The nine $\tau_{ij}$
    were measured on three special planes. Cauchy's formula says the answer on *any* plane is a matrix–vector product,
    $f_i = \tau_{ji}n_j$, and every wall force of the book is $\oint \mathbf n\cdot\boldsymbol\tau\,dA$."
16. `nb.md` — **The idea:** "Slice a corner off the cube: a tiny tetrahedron with three coordinate faces (outward
    normals $-\mathbf e_1, -\mathbf e_2, -\mathbf e_3$) and one slanted face with normal $\mathbf n$. Balance forces on
    it. Face forces scale like (size)², volume forces (weight, inertia) like (size)³ — shrink it and only the faces
    matter. The slanted face's traction must then balance the three coordinate faces'."
17. `nb.primer("limits and orders of smallness", "When an element of size h shrinks, a quantity ∝ h² (an area) shrinks
    more slowly than one ∝ h³ (a volume); their ratio h³/h² = h → 0. So in a balance of face terms (∝ h²) and volume
    terms (∝ h³), the volume terms drop out in the limit — Ch. 1 §1.7 used the same idea once, glossed. Here it is the
    hinge of the derivation.", code="for h in [1e-1, 1e-2, 1e-3]:\n    print(h, h**3 / h**2)     # the volume/area
    ratio is h itself → 0")`
18. `nb.primer("vector area of a closed surface", "Give every patch of a closed surface the vector n dA (outward normal
    times area). Their sum is zero: the surface has no net 'direction' (think of a closed box: each face's vector is
    cancelled by the opposite face's). For the tetrahedron this gives $\mathbf n\,dA = \mathbf e_1 dA_1 + \mathbf e_2 dA_2
    + \mathbf e_3 dA_3$, i.e. $dA_i = n_i\,dA$ — the coordinate face i is the shadow of the slanted face.", code="n =
    np.array([0.6, 0.0, 0.8]); dA = 2.0\nprint(ch02.tetrahedron_face_areas(n, dA))      # [1.2 0.  1.6] = n * dA\nprint
    (n * dA - ch02.tetrahedron_face_areas(n, dA))   # closed-surface sum: zero")`
19. `nb.note` — **N34 [B]** "**A surface element is a vector**: $d\mathbf A = \mathbf n\,dA$ — size dA, direction the
    outward normal. Fluxes (C14) and forces (here) are dot products with it." · **N35 [B]** "**The tetrahedron balance**
    in the 1-direction: $f_1\,dA = \tau_{11}dA_1 + \tau_{21}dA_2 + \tau_{31}dA_3$ with $dA_i = n_i\,dA$ — the book's two
    lines, which are steps 5–7 of the derivation below (`tetrahedron_face_areas`)." equation `f_1\,dA = \tau_{11}\,dA_1
    + \tau_{21}\,dA_2 + \tau_{31}\,dA_3, \qquad dA_i = n_i\,dA`, ref "§2.6".
20. `nb.derivation("D05", "Cauchy's traction formula from a shrinking tetrahedron", ref="2.15", …)` — Part F · D05
    (9 steps).
21. `nb.note` — **N37 [C]** "(2.15) $\mathbf f = \mathbf n\cdot\boldsymbol\tau$ is to τ what $u_n = \mathbf u\cdot\mathbf n$
    is to $\mathbf u$ — contract with the normal — except that $\mathbf f$ is a vector." · ⚠️ callout: "**Which index?**
    $f_i = \tau_{ji}n_j$ contracts the *first* index (row = face). $\boldsymbol\tau\cdot\mathbf n$ contracts the second
    and equals $\mathbf n\cdot\boldsymbol\tau$ only when τ is symmetric — true for the stress (Ch. 4), false for a
    general tensor; the code never assumes it."
22. `nb.worked_example("Ex. 2.2 with a = 1 Pa, φ = 30°", "τ = [[0, a],[a, 0]] (pure shear in a channel, x₁ along the
    flow). 1. $\mathbf n = (\cos30°, \sin30°) = (0.866, 0.5)$. 2. $f_1 = \tau_{11}n_1 + \tau_{21}n_2 = 0 + 1\times0.5 =
    0.5$ Pa. 3. $f_2 = \tau_{12}n_1 + \tau_{22}n_2 = 1\times0.866 + 0 = 0.866$ Pa. 4. $|\mathbf f| = \sqrt{0.25 + 0.75}
    = 1 = |a|$. 5. Direction: $\theta = \arctan2(0.866, 0.5) = 60°$; for a = −1 Pa both components flip → 240°. 6.
    Normal part $\sigma_n = \mathbf f\cdot\mathbf n = 0.5\times0.866 + 0.866\times0.5 = 0.866 = a\sin60°$; shear part
    $\tau_s = |\mathbf f - \sigma_n\mathbf n| = 0.5 = a\cos60°$ — in general $\sigma_n = a\sin2\phi$, $\tau_s =
    a\cos2\phi$ (stated; the Mohr-circle preview of C13).")`
23. `nb.primer("np.arctan2", "`np.arctan2(y, x)` returns the angle of the point (x, y) in the correct quadrant (−π,
    π]; `np.arctan(y/x)` cannot tell (−1, −1) from (1, 1). Add 2π to negative results when you want 0–360°.",
    code="print(np.rad2deg(np.arctan2(0.866, 0.5)))              # 60.0\nprint(np.rad2deg(np.arctan2(-0.866, -0.5)) %
    360)       # 240.0")`
24. `nb.note` — **N38 [B]** Ex. 2.2 (stated, curation §4c D07): the numbers of the tiny example plus the (2.12) route
    "$\tau'_{11} = \sqrt3 a/2 = 0.866a$, $\tau'_{12} = a/2$ — the same two numbers as σ_n and τ_s, because the rotated
    frame's 1'-axis *is* n; repeated inside C06." · **N39 [C]** "Fig. 2.6 (channel, element at 30°) is the 'channel
    shear' preset of the explainer and the slider figure below."
25. `nb.code[explain]` — *code:* `tau2 = ch02.shear_flow_stress(1.0)`; `n = np.array([np.cos(np.deg2rad(30)), np.sin
    (np.deg2rad(30))])`; `f = ch02.traction(tau2, n); print(np.round(f, 4))`; `res = ch02.example_2_2(1.0, np.deg2rad
    (30)); print(np.round(res["magnitude"], 4), np.round(np.rad2deg(res["angle_rad"]), 1), np.round(res["sigma_n"], 4),
    np.round(res["tau_s"], 4))`; `res_neg = ch02.example_2_2(-1.0, np.deg2rad(30)); print(np.round(np.rad2deg(res_neg
    ["angle_rad"]), 1))`; 3-D: `f3 = ch02.traction(tau, [0.6, 0, 0.8]); print(f3)` (tau from C04: (−3·0.6 + 0, 1·0.6,
    −3·0.8) = (−1.8, 0.6, −2.4)). *expect:* [0.5 0.866] · 1.0 60.0 0.866 0.5 · 240.0 · [-1.8 0.6 -2.4]. *explain:* 1.
    the 2×2 shear stress; 2. the unit normal at 30°; 3. (2.15) by einsum 'ji,j->i'; 4. Ex. 2.2 packaged: magnitude,
    angle, normal and shear parts; 5. a < 0 flips the direction; 6. the same formula on the 3-D state of C04.
26. `nb.check_agree` — **from scratch and the wrong way round:** `f_mine = np.array([sum(tau2[j, i]*n[j] for j in
    range(2)) for i in range(2)])` → `assert np.allclose(f_mine, ch02.traction(tau2, n))`; then `tau_ns = np.array
    ([[0., 2.], [0., 0.]])` (non-symmetric): `print(ch02.traction(tau_ns, n), tau_ns @ n)` → `[0. 1.732]` vs `[1. 0.]`
    — "different: the index order is not a formality; only a symmetric τ lets you forget it."
27. `nb.plotly` — Fig. 2.5 tetrahedron: `from scripts.ch02_fig2_5_tetrahedron import tetrahedron_figure;
    tetrahedron_figure(tau, n=[0.6, 0, 0.8])` (coordinate faces with their reversed stress arrows, slanted face with f
    = (−1.8, 0.6, −2.4) as one bold arrow, n as a thin arrow). *see:* "three shadows and one slanted face; the bold
    arrow is the balance of the three faces' arrows." *read:* "f is not along n: it has a normal part (blue) and a
    shear part (rose)." *change:* "n = e₃ gives f = row 3 = (0, 0, −3): Cauchy's formula returns the coordinate face."
28. `nb.plotly` — curation §6 slider (C05): `slider_figure(lambda phi: {...}, "φ", np.linspace(0, 180, 37), unit="°",
    …)` drawing the 2-D element, n, f, σ_n n (blue) and the shear part (rose) for a = +1 and a = −1 (two panels or two
    traces) with the trace name carrying σ_n = a sin 2φ and τ_s = a cos 2φ (`ch02.stress_vs_angle`). *see:* "f turns
    twice as fast as n." *read:* "at 45° f is along n (pure normal stress a, no shear); at 0° and 90° it is pure shear."
    *change:* "a < 0 mirrors every arrow through the origin."
29. `nb.explainer("cauchy_traction_principal_axes", heading="Cut the point any way you like: what pushes on the cut?",
    why="The traction is a function of the cut direction; only dragging n and watching f, σ_n, τ_s and the Mohr point
    move together makes (2.15) a picture — and lets you *find* the shear-free planes before §2.11 proves they are
    eigenvectors.", tries=["Channel-shear preset: drag φ to 45° — the shear vanishes and σ_n = a.", "Hydrostatic
    preset: f stays along n whatever φ; every plane is principal.", "Toggle 'contract second index' on the
    non-symmetric demo: f changes — the index order matters.", "Read the Explain tab at φ = 30°: f₁ = 0·0.866 +
    1·0.5 = 0.5 Pa, term by term."])`
30. `nb.md` — **What would change if…** "…a second observer used rotated axes? Both $\mathbf f$ and $\mathbf n$ are
    vectors (2.8), and (2.15) must hold for both observers — that forces the nine $\tau_{ij}$ to transform in one
    particular way. That way is the definition of a tensor: C06."
31. `nb.core("C06", "The transformation rule of a second-order tensor: τ' = Cᵀ τ C", question="How do the nine stress
    components change when the axes turn — and why does *that* rule, not the array itself, define a tensor?")`
32. `nb.md` — **The problem in plain words:** "The same stress state must give the same force on the same plane
    whichever axes you compute in. Two observers, two 3×3 arrays. The only way they can agree about every plane is if
    the arrays are related by one C per index: $\tau'_{mn} = C_{im}C_{jn}\tau_{ij}$. Ch. 3 will need it for the strain
    rate, Ch. 4 for the stress, Ch. 12 for the Reynolds stress — it is the working definition of 'tensor'."
33. `nb.md` — **The idea:** "(2.8) once per index. A vector: $u'_n = u_iC_{in}$. A tensor has two slots, each behaves
    like a vector: $\tau'_{mn} = C_{im}\,C_{jn}\,\tau_{ij}$. Matrix form: $\boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol
    \tau\,\mathbf C$ — Cᵀ on the left for the first index, C on the right for the second."
34. `nb.derivation("D06", "The tensor transformation rule from Cauchy's formula", ref="2.12", …)` — Part F · D06
    (8 steps). Preceded by one sentence: "The book states (2.12) and points to a tetrahedron argument (Sommerfeld); we
    have the tetrahedron result (2.15) in hand, so we derive (2.12) from it."
35. `nb.note` — **N27 [B]** "**Tensor ≠ matrix.** Any 3×3 array is a matrix; its entries are the components of a tensor
    only if they obey (2.12) in every frame. The code below declares an array 'the same in every frame' and measures the
    (2.12) residual: not small." · **N28 [B]** (2.13): "**Fourth order**: one C per index, four of them, $3^4 = 81$
    components — `transform_tensor(A4, C)`; the viscosity tensor of Ch. 4 §4.5 is the one we will meet." equation
    `A'_{mnpq} = C_{im} C_{jn} C_{kp} C_{lq} A_{ijkl}`, ref "2.13". · **N29 [B]** "**Examples that pass**: the stress
    $\tau_{ij}$, the velocity gradient $\partial u_i/\partial x_j$ (N52), and the outer product $u_iv_j$ of two vectors
    (Exercise 2.10) — asserted below with `outer`."
36. `nb.worked_example("Ex. 2.2's rotated stress, a = 1 Pa, θ = 30°", "1. $C = \begin{bmatrix}0.866 & -0.5\\ 0.5 &
    0.866\end{bmatrix}$. 2. $\tau'_{11} = C_{i1}C_{j1}\tau_{ij} = C_{11}C_{21}\tau_{12} + C_{21}C_{11}\tau_{21} = 2\times
    0.866\times0.5\times1 = 0.866$ Pa. 3. $\tau'_{12} = C_{11}C_{22}\tau_{12} + C_{21}C_{12}\tau_{21} = 0.866^2 - 0.5^2 =
    0.75 - 0.25 = 0.5$ Pa. 4. $\tau'_{22} = 2C_{12}C_{22}\tau_{12} = -2\times0.5\times0.866 = -0.866$ Pa. 5. So $\tau'
    = \begin{bmatrix}0.866 & 0.5\\ 0.5 & -0.866\end{bmatrix}$: the 1'-axis is n, so $\tau'_{11} = \sigma_n$ and
    $\tau'_{12} = \tau_s$ of C05 ✓; the trace is still 0 (C07).")`
37. `nb.code[explain]` — *code:* `C2 = ch02.rotation_matrix_2d(np.deg2rad(30))`; `tp = ch02.transform_tensor(tau2,
    C2); print(np.round(tp, 4))`; `print(np.round(C2.T @ tau2 @ C2, 4))`; `print(np.round(ch02.example_2_2(1.0,
    np.deg2rad(30))["tau_rot"], 4))`; `u, v = np.array([1., 2., 3.]), np.array([0., 1., -1.]); R = ch02.random_rotation
    (np.random.default_rng(2))`; `print(ch02.transforms_as_tensor(lambda x: ch02.outer(u, v), R, pts))` — hmm, `outer`
    of constant vectors is frame-dependent; use position-dependent: `lambda x: ch02.outer(x, np.cross([0,0,1], x))`
    → residual ≤ 1e-14; `print(ch02.transforms_as_tensor(lambda x: np.diag([1., 2., 3.]), R, pts))` → order 1 (N27);
    `A4 = np.random.default_rng(3).normal(size=(3, 3, 3, 3)); print(ch02.transform_tensor(A4, R).shape)`. *expect:*
    [[0.866 0.5],[0.5 −0.866]] three times · ~1e-15 · ~1 · (3, 3, 3, 3). *explain:* 1. (2.12) as `transform_tensor`
    (einsum 'im,jn,ij->mn'); 2. the matrix form CᵀτC — identical; 3. Ex. 2.2's route; 4. an outer product of two
    vector fields passes the tensor test; 5. a fixed array "the same in every frame" fails; 6. (2.13): four C's, 81
    numbers.
38. `nb.check_agree` — **from scratch:** `tp_mine = np.zeros((2, 2)); for m in range(2): for n_ in range(2): for i in
    range(2): for j in range(2): tp_mine[m, n_] += C2[i, m]*C2[j, n_]*tau2[i, j]` → `assert np.allclose(tp_mine,
    ch02.transform_tensor(tau2, C2))`; `assert np.allclose(tp_mine, C2.T @ tau2 @ C2)`.
39. `nb.plotly` — curation §6 slider (C06): θ from 0 to 180° (37 steps): left panel $\tau'_{11}, \tau'_{12}, \tau'_{22}$
    vs θ for the pure-shear (a = 1) and uniaxial ([[1, 0],[0, 0]]) presets with a moving dot; right panel Mohr's circle
    ($\sigma_n$ vs $\tau_s$, centre I₁/2, radius from `mohr_circle_2d`) with the current point. *see:* "three cosine
    curves of period 180°; the dot runs round a circle." *read:* "τ'₁₂ crosses zero at 45° for pure shear — a
    principal axis (C13); the circle's radius is the largest shear any plane can carry." *change:* "add a pressure
    (−p on the diagonal): the curves shift down by p, the circle slides left; nothing else changes."
40. `nb.md` — **What would change if…** "…you contracted the two indices of τ' with each other? $\tau'_{mm} = C_{im}
    C_{jm}\tau_{ij} = \delta_{ij}\tau_{ij} = \tau_{ii}$: the trace is the same in every frame. That is the next
    section — contraction makes invariants."

### A.5 §2.5 Contraction and Multiplication — C07 (+N30–N33 · D18)
1. `nb.section("2.5", "Contraction and Multiplication", intro="**What is this section about?** Setting two indices
   equal and summing (contraction) lowers the order by two; multiplying raises it. Contracted quantities with no free
   index are the same in every frame — the invariants of a tensor — and the einsum patterns here are used in every
   later notebook.")`
2. `nb.core("C07", "Contraction: the trace, the three invariants and the double dot", question="Which combinations of
   the nine numbers $\tau_{ij}$ are the same for every observer — and why exactly those?")`
3. `nb.md` — **The problem in plain words:** "The pressure in a flowing fluid is defined in Ch. 4 as $-\tau_{ii}/3$;
   the rate of volume expansion in Ch. 3 is $S_{ii}$; turbulence models in Ch. 12 are built from invariants of the
   Reynolds stress. Each is a number that must not depend on how you set up your axes. Contraction is the recipe that
   produces such numbers."
4. `nb.md` — **The idea:**
   ```
   A_ij  →  set j = i and sum  →  A_ii = A_11 + A_22 + A_33        (a chain of indices closing on itself)
   under a rotation each index brings one C; a closed pair  C_im C_jm  is  δ_ij  (D02)  →  the C's vanish
   A_ij A_ji  = a closed chain of two  →  invariant, and it is the one that enters I₂;    A_ij A_ij = tr(A Aᵀ) is also invariant* — but it is NOT the λ-coefficient of the characteristic polynomial
   ```
   "(*both scalars survive a rotation — CᵀAC leaves tr(A Aᵀ) unchanged too. The distinction is that I₂ = ½(I₁² − A_ij A_ji) needs the closed chain A_ij A_ji; the two coincide when A is symmetric.)"
5. `nb.primer("Vieta's formulas", "For a cubic with roots $\lambda^1, \lambda^2, \lambda^3$: $(\lambda - \lambda^1)
   (\lambda - \lambda^2)(\lambda - \lambda^3) = \lambda^3 - (\sum\lambda^k)\lambda^2 + (\sum_{k<l}\lambda^k\lambda^l)
   \lambda - \lambda^1\lambda^2\lambda^3$ — the coefficients are the sum, the pair-sums and the product of the roots.",
   code="import sympy as sp\nl, a, b, c = sp.symbols('lambda a b c')\nprint(sp.expand((l - a)*(l - b)*(l - c)))   #
   lambda**3 - (a+b+c) lambda**2 + (ab+ac+bc) lambda - abc")`
6. `nb.md` — **The maths, step by step:** "1. Contraction: $\sum_j A_{jj} \equiv A_{jj} = A_{11} + A_{22} + A_{33}$, the
   trace. 2. Invariance: apply (2.12) and set n = m — D18 below shows every closed index chain survives the rotation
   unchanged. 3. The three invariants: $I_1 = A_{ii}$, $I_2 = \tfrac12(I_1^2 - A_{ij}A_{ji})$, $I_3 = \det\mathbf A$;
   they are the coefficients of $\det(\mathbf A - \lambda\boldsymbol\delta) = 0$ (Exercise 2.9), so the eigenvalues of
   C13 are invariants too. 4. Products raise the order (N30); contractions of products (2.14) are matrix products in
   disguise (N31)."
7. `nb.derivation("D18", "Why the trace, I₂ and the determinant do not depend on the axes", ref="Exercise 2.9", …)` —
   Part F · D18 (8 steps).
8. `nb.note` — **N30 [C]** "**Multiplication raises the order**: $P_{ijkl} = A_{ij}B_{kl}$ has four free indices
   (`tensor_product(A, B)` = `np.multiply.outer`); it transforms by (2.13)." · **N31 [B]** (2.14) as a 4-row table
   "index string → matrix form": `A_ij B_ki` → (B·A)_kj · `A_ij B_ik` → (Aᵀ·B)_jk · `A_ij B_kj` → (A·Bᵀ)_ik · `A_ij B_jk`
   → (A·B)_ik, each asserted with `contract(A, B, pattern)`; "rearrange until the summed index is adjacent, then it is
   a matrix product." equation `A_{ij}B_{ki} = (\mathbf B\cdot\mathbf A)_{kj}, \quad A_{ij}B_{ik} = (\mathbf A^{\rm T}
   \cdot\mathbf B)_{jk}, \quad A_{ij}B_{kj} = (\mathbf A\cdot\mathbf B^{\rm T})_{ik}, \quad A_{ij}B_{jk} = (\mathbf A
   \cdot\mathbf B)_{ik}`, ref "2.14". · **N32 [B]** "**Tensor · vector, two ways**: $A_{ij}u_j = (\mathbf A\cdot\mathbf
   u)_i$ vs $A_{ij}u_i = (\mathbf A^{\rm T}\cdot\mathbf u)_j$ — `dot_tensor_vector(A, u, index=1 or 0)`; exactly the
   distinction of C05 ($\mathbf n\cdot\boldsymbol\tau$ vs $\boldsymbol\tau\cdot\mathbf n$)." · **N33 [B]** ⚠️ "**Double
   dot, two conventions**: the book's $\mathbf A:\mathbf B = A_{ij}B_{ji}$ (transpose pairing) vs the common Frobenius
   $A_{ij}B_{ij}$ (Wikipedia, most continuum texts). For A = [[1, 2],[3, 4]], B = [[0, 1],[1, 0]]: book 2·1 + 3·1 = 5,
   Frobenius 2·1 + 3·1 = 5 too — symmetric B hides the difference; take B = [[0, 1],[0, 0]]: book A₁₂B₂₁ + A₂₁B₁₂ = 3,
   Frobenius A₁₂B₁₂ = 2. They agree whenever one operand is symmetric (Ch. 4's dissipation $\tau_{ij}\partial u_i/
   \partial x_j$). `double_dot(A, B, convention=…)` — say which."
9. `nb.worked_example("invariants of A = [[2, 1, 0],[1, 3, 1],[0, 1, 4]]", "1. $I_1 = 2 + 3 + 4 = 9$. 2. $A_{ij}A_{ji}$
   = sum of squares (A symmetric) = $4 + 9 + 16 + 2(1 + 0 + 1) = 33$, so $I_2 = \tfrac12(81 - 33) = 24$. 3. $I_3 =
   2(3\cdot4 - 1\cdot1) - 1(1\cdot4 - 0) + 0 = 22 - 4 = 18$. 4. Characteristic polynomial $\lambda^3 - 9\lambda^2 +
   24\lambda - 18 = 0$; its roots multiply to 18 and add to 9 (Vieta) — C13 finds them numerically: 1.27, 3.00, 4.73
   (sum 9.00 ✓).")`
10. `nb.code[explain]` — *code:* `A = np.array([[2., 1, 0], [1, 3, 1], [0, 1, 4]])` (the notebook-local `EXAMPLE_TENSOR`, defined here and reused in C13); `print(ch02.trace(A), ch02.invariants
    (A))`; `print(ch02.characteristic_polynomial(A))`; `R = ch02.random_rotation(np.random.default_rng(4)); Ap =
    ch02.transform_tensor(A, R); print(np.round(Ap, 3)); print(np.round(ch02.invariants(Ap), 10))`; `print(ch02.
    double_dot(A, Ap, convention="book"), ch02.double_dot(A, Ap, convention="frobenius"))` (equal: both symmetric);
    `B = np.array([[0., 1.], [0., 0.]]); A2 = np.array([[1., 2.], [3., 4.]]); print(ch02.double_dot(A2, B, "book"),
    ch02.double_dot(A2, B, "frobenius"))`. *expect:* 9.0 (9.0, 24.0, 18.0) · [1, −9, 24, −18] · a scattered 3×3 ·
    (9.0, 24.0, 18.0) · two equal numbers · 3.0 2.0. *explain:* 1–2. the invariants and the polynomial; 3. rotate: the
    entries change, the three invariants do not; 4–5. the two double-dot conventions, equal for symmetric operands,
    different otherwise.
11. `nb.check_agree` — **from scratch:** `I1 = sum(A[i, i] for i in range(3))`; `I2 = 0.5*(I1**2 - sum(A[i, j]*A[j, i]
    for i in range(3) for j in range(3)))`; `I3 = A[0,0]*(A[1,1]*A[2,2] - A[1,2]*A[2,1]) - A[0,1]*(A[1,0]*A[2,2] -
    A[1,2]*A[2,0]) + A[0,2]*(A[1,0]*A[2,1] - A[1,1]*A[2,0])` (cofactors, P53) → `assert np.allclose((I1, I2, I3),
    ch02.invariants(A))`; also the four (2.14) patterns: `for pat, mat in [("ij,ki->kj", B3 @ A), …]: assert
    np.allclose(ch02.contract(A, B3, pat), mat)`.
12. `nb.figure` — invariants flat over 50 random rotations: x = rotation number, three flat lines I₁, I₂, I₃ (purple
    shades) and the nine entries $A'_{ij}$ as grey scattered dots; `rng = np.random.default_rng(5)`; loop `R =
    ch02.random_rotation(rng)`. *see:* "grey dots wander over [−2, 5]; three lines do not move." *read:* "each dot is
    an observer's opinion about one component; the lines are facts." *change:* "use a non-symmetric A and plot
    $A_{ij}A_{ij}$ too: it wanders — not a closed chain."
13. `nb.md` — **What would change if…** "…you contracted a *third-order* object? $\varepsilon_{ijk}$ contracted with two
    vectors is the cross product — the alternating tensor of the next section."

### A.6 §2.7, 2.8 Kronecker Delta, Alternating Tensor, and the Cross Product — C08 (+N42–N44, N46–N50 · D09)
1. `nb.section("2.7, 2.8", "Kronecker Delta and Alternating Tensor; Vector, Dot, and Cross Products", intro="**What is
   this section about?** Two constant tensors that look the same to every observer: δ (met in C01) and the alternating
   tensor ε_ijk. ε turns the cross product, the curl and the rotation tensor into index sums, and one identity — the
   epsilon–delta relation — closes every vector identity you will ever need.")`
2. `nb.core("C08", "The alternating tensor ε_ijk", question="How can 27 numbers, all 0 or ±1, encode 'right-handed
   perpendicular' — and give the cross product, the curl and the Coriolis term as index sums?")`
3. `nb.md` — **The problem in plain words:** "The Coriolis acceleration $2\boldsymbol\Omega\times\mathbf u$, the vorticity
   $\nabla\times\mathbf u$, the torque $\mathbf r\times\mathbf F$ — fluid mechanics is full of cross products, and the
   component formula (2.20) is a mess to manipulate. ε_ijk is the bookkeeping device that turns 'cross' into 'sum',
   and (2.19) is the one identity that turns products of two ε's back into δ's."
4. `nb.md` — **The idea:**
   ```
   ε_ijk = +1  for ijk = 123, 231, 312   (cyclic: walk round 1→2→3→1)
         = −1  for ijk = 321, 213, 132   (anticyclic: walk backwards)
         =  0  whenever two indices agree                 → 6 nonzero entries out of 27
   swap two indices → sign flips;  move one index two places (a cyclic shift) → sign kept
   ```
5. `nb.primer("permutations, cyclic order and parity", "A permutation of (1, 2, 3) is a reordering. It is *even* if it
   takes an even number of pairwise swaps to reach from 123 (123, 231, 312 — the cyclic shifts) and *odd* otherwise
   (132, 213, 321). ε_ijk is +1 on even, −1 on odd permutations, 0 if not a permutation.", code="from itertools import
   permutations\nfor p in permutations((1, 2, 3)):\n    print(p, ch02.permutation_sign(*p))   # (1,2,3) 1, (1,3,2) -1,
   (2,1,3) -1, (2,3,1) 1, (3,1,2) 1, (3,2,1) -1")`
6. `nb.primer("numpy arrays with three axes and np.transpose", "`eps[i, j, k]` indexes a 3×3×3 array (shape (3, 3, 3));
   `np.transpose(eps, (1, 2, 0))` reorders the axes so that `new[i, j, k] = eps[j, k, i]` — the code form of 'move an
   index'.", code="eps = ch02.levi_civita()\nprint(eps.shape, eps[0, 1, 2], eps[0, 2, 1], eps[0, 0, 1])   # (3,3,3) 1
   -1 0  (Python 0-based = ε_123, ε_132, ε_112)\nprint(np.array_equal(np.transpose(eps, (1, 2, 0)), eps))   # True:
   two-place move keeps ε")`
7. `nb.primer("right-hand rule and orientation", "Point the fingers of your right hand along the first vector, curl
   them toward the second: the thumb gives the direction of the cross product. e₁ × e₂ = e₃ fixes a right-handed
   frame; the same rule orients the boundary of a surface in §2.13 (thumb along n, fingers along t).", code="print
   (np.cross([1, 0, 0], [0, 1, 0]))   # [0 0 1] = e3: right-handed")`
8. `nb.md` — **The maths, step by step:** "1. Definition (2.18) as in the idea box. 2. **Index moves** (N43): swapping
   two indices flips the sign, $\varepsilon_{ijk} = -\varepsilon_{ikj}$; moving one index two places keeps it,
   $\varepsilon_{ijk} = \varepsilon_{jki} = \varepsilon_{kij}$. 3. **Cross product** (2.21): $(\mathbf u\times\mathbf
   v)_k = \varepsilon_{ijk}u_iv_j$ — nine terms, two survive for each k; k = 1: $\varepsilon_{231}u_2v_3 + \varepsilon_
   {321}u_3v_2 = u_2v_3 - u_3v_2$ ✓ (2.20). 4. **Epsilon–delta** (2.19): $\varepsilon_{ijk}\varepsilon_{klm} = \delta_{il}
   \delta_{jm} - \delta_{im}\delta_{jl}$; contracted once $\varepsilon_{pqi}\varepsilon_{pqj} = 2\delta_{ij}$, twice
   $\varepsilon_{pqr}\varepsilon_{pqr} = 6$ — D09. 5. **Isotropy** (N42): δ and ε (up to a factor) are the only
   constant tensors of order 2 and 3 unchanged by every proper rotation."
9. `nb.note` — **N46 [B]** "**Cross product, the geometric definition**: $\mathbf w = \mathbf u\times\mathbf v$ has
   length $uv\sin\theta$, is perpendicular to both, and $(\mathbf u, \mathbf v, \mathbf w)$ is right-handed; $\mathbf u
   \times\mathbf v = -\mathbf v\times\mathbf u$; $\mathbf e_1\times\mathbf e_2 = \mathbf e_3$." · **N47 [B]** (2.20)
   stated (curation §4c D10: 'because $\mathbf e_i\times\mathbf e_j = \varepsilon_{ijk}\mathbf e_k$ and the product
   distributes'): equation `\mathbf u\times\mathbf v = (u_2v_3 - u_3v_2)\mathbf e_1 + (u_3v_1 - u_1v_3)\mathbf e_2 +
   (u_1v_2 - u_2v_1)\mathbf e_3`, ref "2.20"; `cross(u, v)` asserted equal to `np.cross`. · **N48 [C]** "Determinant
   form: $\mathbf u\times\mathbf v = \det[\mathbf e_1\ \mathbf e_2\ \mathbf e_3;\ u_1\ u_2\ u_3;\ v_1\ v_2\ v_3]$ —
   one sympy line reproduces (2.20)." · **N49 [B]** (2.21): equation `(\mathbf u\times\mathbf v)_k = \sum_{i=1}^{3}
   \sum_{j=1}^{3}\varepsilon_{ijk}u_iv_j \equiv \varepsilon_{ijk}u_iv_j = \varepsilon_{kij}u_iv_j`, ref "2.21";
   `cross_einsum` = `np.einsum('ijk,i,j->k')`; `expand_indices_str("eps_ijk u_i v_j")` prints (2.20). · **N50 [C]**
   "k = 1 check: only (i, j) = (2, 3) and (3, 2) survive — visible in the printed expansion."
10. `nb.note` — **N43 [B]** "**Index moves on ε** — with `np.transpose` assertions: `(1, 2, 0)` and `(2, 0, 1)` return
    ε; `(0, 2, 1)` returns −ε." · **N42 [B]** "**Isotropic tensors**: $\delta'_{ij} = C_{im}C_{jn}\delta_{mn} =
    C_{im}C_{jm} = \delta_{ij}$ (D02) in every frame; ε is isotropic under proper rotations (det C = +1) and flips sign
    under a reflection — a *pseudo*tensor, which is why the book says 'rotation'. `is_isotropic(delta)` True,
    `is_isotropic(eps, proper=True)` True, `is_isotropic(eps, proper=False)` False, a random tensor False. Pointer:
    the isotropic Newtonian stress model of Ch. 4 §4.5 is built from δ alone." · **N44 [B]** (2.19) with equation
    `\sum_{k=1}^{3}\varepsilon_{ijk}\varepsilon_{klm} \equiv \varepsilon_{ijk}\varepsilon_{klm} = \delta_{il}\delta_{jm}
    - \delta_{im}\delta_{jl}`, ref "2.19": "the book says 'verify by choosing values'; D09 below proves it, and
    `epsilon_delta_residual()` checks all 81 cases; the vector triple product $\mathbf a\times(\mathbf b\times\mathbf c)
    = (\mathbf a\cdot\mathbf c)\mathbf b - (\mathbf a\cdot\mathbf b)\mathbf c$ is its first payoff."
11. `nb.derivation("D09", "The epsilon–delta relation, its contractions and the triple product", ref="2.19", …)` —
    Part F · D09 (11 steps).
12. `nb.worked_example("one case of (2.19) and one cross product", "1. Take i = 1, j = 2, l = 1, m = 2: left side
    $\sum_k\varepsilon_{12k}\varepsilon_{k12} = \varepsilon_{123}\varepsilon_{312} = 1\times1 = 1$ (k = 3 is the only
    term); right side $\delta_{11}\delta_{22} - \delta_{12}\delta_{21} = 1 - 0 = 1$ ✓. 2. Swap l and m (l = 2, m = 1):
    left $\varepsilon_{123}\varepsilon_{321} = -1$, right $\delta_{12}\delta_{21} - \delta_{11}\delta_{22} = -1$ ✓. 3.
    Cross product (1, 2, 3) × (4, 5, 6): $k = 1$: $2\cdot6 - 3\cdot5 = -3$; $k = 2$: $3\cdot4 - 1\cdot6 = 6$; $k = 3$:
    $1\cdot5 - 2\cdot4 = -3$ → (−3, 6, −3); check ⊥: $(1, 2, 3)\cdot(-3, 6, -3) = -3 + 12 - 9 = 0$ ✓.")`
13. `nb.code[explain]` — *code:* `eps = ch02.levi_civita(); print(eps[0, 1, 2], eps[2, 1, 0], np.count_nonzero(eps))`;
    `print(ch02.epsilon_delta_residual())`; `print(np.einsum('pqi,pqj->ij', eps, eps)); print(np.einsum('pqr,pqr',
    eps, eps))`; `u, v = np.array([1., 2., 3.]), np.array([4., 5., 6.]); print(ch02.cross(u, v), ch02.cross_einsum(u,
    v), np.cross(u, v))`; `print(ch02.expand_indices_str("eps_ijk u_i v_j"))`; `a, b, c = rng.normal(size=(3, 3));
    print(np.allclose(ch02.triple_product(a, b, c), (a @ c)*b - (a @ b)*c))`; `print(ch02.is_isotropic(np.eye(3),
    rng)[0], ch02.is_isotropic(eps, rng, proper=True)[0], ch02.is_isotropic(eps, rng, proper=False)[0], ch02.
    is_isotropic(A, rng)[0])`. *expect:* 1 −1 6 · 0.0 · 2·I₃ · 6 · [−3 6 −3] ×3 · the three-component expansion ·
    True · True True False False. *explain:* 1. six nonzero entries; 2. all 81 cases of (2.19) hold exactly; 3. the
    two contractions 2δ and 6; 4. three cross products agree; 5. the expansion shows the k = 1 check; 6. the triple
    product identity; 7. isotropy: δ always, ε only for proper rotations, a random tensor never.
14. `nb.check_agree` — **from scratch:** `eps_mine = np.zeros((3, 3, 3)); for p in permutations(range(3)): eps_mine[p]
    = ch02.permutation_sign(*[q + 1 for q in p])` → `assert np.array_equal(eps_mine, eps)`; 81-case check as a
    quadruple loop `for i, j, l, m in product(range(3), repeat=4): lhs = sum(eps[i,j,k]*eps[k,l,m] for k in range(3));
    rhs = (i==l)*(j==m) - (i==m)*(j==l); assert lhs == rhs`.
15. `nb.figure` — ε as an unfolded 3×9 coloured grid (three 3×3 slices k = 1, 2, 3; +1 purple, −1 green, 0 white) with
    the six nonzero cells labelled by their index triple. *see:* "six coloured cells in 27, a diagonal pattern shifting
    by one per slice." *read:* "each slice k holds the 2-D antisymmetric pattern of the remaining two indices — the
    matrix (2.26) of C12 in disguise." *change:* "relabel 1 ↔ 2 everywhere (a mirror): every colour flips."
16. `nb.md` — **What would change if…** "…the second vector were the operator $\partial/\partial x_j$? $\varepsilon_{ijk}
    \partial u_k/\partial x_j$ — the curl (C11). And $-\varepsilon_{ijk}\omega_k$ packs a vector into an antisymmetric
    matrix (C12). ε is the bridge between vectors and antisymmetric tensors."

### A.7 §2.9 Gradient, Divergence, and Curl — R02, C09 (+N51), C10 (+N52), C11 (+N53–N55 · D12)
1. `nb.section("2.9", "Gradient, Divergence, and Curl", intro="**What is this section about?** The del operator acting
   three ways: on a scalar it gives the direction of steepest climb (gradient), on a vector either how much the field
   spreads (divergence) or how much it spins (curl). We compute all three on a grid you will meet again in every
   chapter, and check the stencils converge at second order.")`
2. `nb.recap("R02", "The del operator ∇", "Ch. 1 primed the partial derivative (P25) and ∇ as the vector that 'points
   uphill', and finite differences as its numerical form (P21, P22). New here: the index form $\nabla = \mathbf e_i\,
   \partial/\partial x_i$ (2.22) — an operator with one free index, so it raises the order of whatever it acts on by
   one, or lowers it by one when contracted.", where="Ch. 1 §1.5 primers P21, P22, P25")`
3. `nb.core("C09", "The gradient ∇φ: perpendicular to the level sets, along the steepest climb", question="On a
   contour map of a scalar φ, which way is steepest, how steep is it, and how steep is it in some *other* direction
   n?")`
4. `nb.md` — **The problem in plain words:** "Pressure maps drive the wind: air accelerates from high to low pressure,
   fastest where the isobars are closest — down $-\nabla p$. A temperature map drives heat flux down $-\nabla T$ (Ch. 1
   Fourier). The gradient turns a scalar map into a vector field: direction of fastest increase, magnitude the rate."
5. `nb.md` — **The idea:**
   ```
   contours of φ (level sets)  ──  ∇φ is perpendicular to them  (moving along a contour changes nothing)
   |∇φ| = the slope in that steepest direction;  in any other unit direction n:  ∂φ/∂n = ∇φ · n = |∇φ| cos(angle)
   ```
6. `nb.primer("level sets and the directional derivative", "A level set (contour line in 2-D, surface in 3-D) is where
   φ takes one value. The directional derivative ∂φ/∂n is the rate of change of φ per metre when you walk along the
   unit vector n; by the chain rule (P49) it equals ∇φ·n, so it is zero along a contour (n ⊥ ∇φ) and largest along
   ∇φ.", code="phi = lambda x, y: x**2 + y**2\ng = np.array([2*1.0, 2*2.0])                  # ∇φ at (1, 2) = (2, 4)\nfor n
   in ([1, 0], [0, 1], g/np.linalg.norm(g), [-2, 1]/np.sqrt(5)):\n    print(np.round(g @ np.asarray(n), 3))     # 2.0
   4.0 4.472 0.0 (along the contour)")`
7. `nb.primer("np.meshgrid and the project grid layout", "`np.meshgrid` turns 1-D coordinate arrays into arrays that
   hold the coordinate of every grid point. **Project convention (every later chapter):** 3-D arrays are indexed `[k,
   j, i]` = (z, y, x) — x on the last axis — built with `np.meshgrid(z, y, x, indexing='ij')`; 2-D arrays are `[j, i]`
   = (y, x) from `indexing='xy'`. Vector fields stack the component on axis 0: `u[c, k, j, i]`. `core.grids.grid` and
   `grid2d` do this once so nobody transposes a curl by accident.", code="g2 = ch02.grid2d(((-1, 1), (-1, 1)), 5)\nprint
   (g2.X.shape, g2.X[0], g2.Y[:, 0])   # (5, 5); x varies along the last axis, y along the first\ng3 = ch02.grid(((0,
   1), (0, 2), (0, 3)), 4)\nprint(g3.X.shape, g3.h)              # (4, 4, 4) as (z, y, x); spacings (hx, hy, hz)")`
8. `nb.primer("numpy broadcasting", "An operation between a full grid array and a scalar (or a lower-dimensional
   array) is applied at every grid point at once; `X**2 + Y**2` on meshgrid arrays evaluates φ everywhere without a
   loop. Slices like `phi[:, 2:] - phi[:, :-2]` subtract neighbours along one axis for the whole grid.", code="X, Y =
   g2.X, g2.Y\nphi = X**2 + Y**2                      # φ at all 25 points\nprint((phi[:, 2:] - phi[:, :-2]).shape)   #
   (5, 3): central differences along x, interior only")`
9. `nb.primer("plt.contour, plt.quiver and plt.streamplot", "`ax.contour(X, Y, phi, levels)` draws level curves;
   `ax.quiver(X, Y, U, V)` draws an arrow at each point; `ax.streamplot(x, y, U, V)` draws curves tangent to a vector
   field (1-D x, y arrays; the field on an `[j, i]` grid — our layout).", code="fig, ax = plt.subplots(figsize=(3,
   3))\nax.contour(X, Y, phi, levels=6); ax.quiver(X, Y, 2*X, 2*Y, scale=40); ax.set_aspect('equal'); plt.show()")`
10. `nb.md` — **The maths, step by step:** "1. (2.22) $\nabla = \mathbf e_i\,\partial/\partial x_i$; on a scalar,
    $(\nabla\phi)_i = \partial\phi/\partial x_i$ — a vector (it passes (2.8) by the chain rule; the residual table of
    C03 said so). 2. Along a contour φ is constant, so its rate of change along the tangent t is $\nabla\phi\cdot\mathbf
    t = 0$: ∇φ ⊥ level sets. 3. In direction n: $\partial\phi/\partial n = \nabla\phi\cdot\mathbf n = |\nabla\phi|\cos
    \alpha$ — maximal ($|\nabla\phi|$) for n ∥ ∇φ (Fig. 2.7). 4. On the grid: the central difference $(\phi_{i+1} -
    \phi_{i-1})/2h$ with error ∝ h² (P21); at the edges a one-sided second-order stencil $(-3\phi_i + 4\phi_{i+1} -
    \phi_{i+2})/2h$ keeps the order."
11. `nb.worked_example("φ = x² + y² at (1, 2)", "1. $\partial\phi/\partial x = 2x = 2$, $\partial\phi/\partial y = 2y =
    4$: ∇φ = (2, 4), $|\nabla\phi| = \sqrt{20} = 4.472$. 2. Along $\mathbf n = (1, 0)$: $\partial\phi/\partial n = 2$.
    3. Along the contour direction $(-2, 1)/\sqrt5$: $(-4 + 4)/\sqrt5 = 0$. 4. Along ∇φ itself: 4.472 — the steepest.
    5. Central difference with h = 0.1 at x = 1: $(1.1^2 - 0.9^2)/0.2 = (1.21 - 0.81)/0.2 = 2.000$ — exact, because a
    quadratic has no third derivative.")`
12. `nb.code[explain]` — *code:* `g = ch02.grid2d(((-2, 2), (-2, 2)), 41)` (h = 0.1); `phi = g.X**2 + g.Y**2`; `grad
    = ch02.gradient(phi, g.h)` (shape (2, 41, 41)); `j, i = 30, 30` (y = 1.0? — pick indices for (1, 2): x index of 1.0
    is 30, y index of 2.0 is 40 → edge; use point (1, 1): i = j = 30); `print(grad[:, 30, 30])`; `print(ch02.
    directional_derivative(grad, [1, 0])[30, 30], ch02.directional_derivative(grad, np.array([-1, 1])/np.sqrt(2))[30,
    30])`; `ns = (8, 16, 32) if FAST else (8, 16, 32, 64); study = {op: ch02.operator_convergence(op, ns) for op in ("gradient",
    "divergence", "curl")}; print({op: round(s["order"], 2) for op, s in study.items()})`.
    *expect:* [2. 2.] (exact) · 2.0 0.0 · order ≈ 2.0 (each of gradient, divergence, curl). *explain:* 1. a 41×41 grid
    on [−2, 2]²; 2. φ at every point by broadcasting; 3. `gradient` = `partial` along each axis, components stacked on
    axis 0 as (x, y); 4. the value at (1, 1); 5. directional derivatives; 6. the cached convergence study (reused in
    C15, C16): errors of all three operators on a sin/cos field vs the exact sympy derivatives, slope 2 on log–log
    axes.
13. `nb.check_agree` — **from scratch (curation §7):** `dphidx = np.zeros_like(phi); dphidx[:, 1:-1] = (phi[:, 2:] -
    phi[:, :-2])/(2*g.h[0]); dphidx[:, 0] = (-3*phi[:, 0] + 4*phi[:, 1] - phi[:, 2])/(2*g.h[0]); dphidx[:, -1] =
    (3*phi[:, -1] - 4*phi[:, -2] + phi[:, -3])/(2*g.h[0])` → `assert np.allclose(dphidx, ch02.partial(phi, 0, g.h[0]))` (direction 0 = x; the library maps it to array axis 1
    of the `[j, i]` layout) and `assert np.allclose(dphidx, grad[0])`. Markdown: "x is the *last* array axis in the
    `[j, i]` layout — the one place to get this wrong, so `partial` takes the coordinate direction, not the axis."
14. `nb.figure` — **N51 [B]** Fig. 2.7 as ours: contours of $\phi = x^2 + y^2/4$ (ellipses), the ∇φ quiver (purple),
    at the probe (1, 1) the gradient arrow bold and a direction n at 30° (orange) with the label "∂φ/∂n = ∇φ·n =
    2.0·0.866 + 0.5·0.5 = 1.98". *see:* "arrows cross every contour at right angles and lengthen where contours crowd."
    *read:* "|∇φ| is the local contour density; n's projection on ∇φ is the slope you would feel walking along n."
    *change:* "rotate n to lie along a contour: ∂φ/∂n → 0."
15. `nb.plotly` — curation §6 slider (C09): direction angle of n, 36 steps (0–350°); contour plot fixed, the n arrow
    turning at the probe, the trace name "∂φ/∂n = 1.98" updating (`directional_derivative`). *see/read/change:*
    "maximum when n ∥ ∇φ, zero along the contour, negative downhill."
16. `nb.figure` — the convergence plot from `study`: log–log error vs h for gradient / divergence / curl with a dashed
    slope-2 reference (`observed_order` printed in the legend). *see:* "three straight lines of slope 2." *read:*
    "halving h quarters the error — the stencil is second order also at the edges (one-sided stencil)." *change:*
    "replace the edge stencil by first order: the slope of the edge error drops to 1 (the ch01 `np.gradient` lesson)."
    Gloss (ledger "observed order of convergence"): "the slope of log(error) vs log(h) — Ch. 1 P13 read slopes on
    log–log axes; `tools.convergence.observed_order` fits it."
17. `nb.md` — **What would change if…** "…∇ acted on a vector instead? Two choices: contract the operator's index with
    the vector's (divergence, C10) or cross them with ε (curl, C11)."
18. `nb.core("C10", "The divergence ∇·u = ∂u_i/∂x_i", question="How much does a velocity field *spread out* at a
    point — and why is that one number the most used quantity in the book?")`
19. `nb.md` — **The problem in plain words:** "Mass conservation (Ch. 4): if more fluid leaves a tiny box than enters,
    the density inside must drop. 'More leaves than enters' per unit volume is the divergence. For water and slow air
    ∇·u = 0 — the incompressibility condition that shapes every chapter after this one."
20. `nb.md` — **The idea:**
    ```
    ∂u_1/∂x_1 : does u_1 grow along x_1?  (stretching along 1)        ∇·u = sum of the three stretching rates
    a box with more leaving than entering  →  ∇·u > 0  (source);  the same in and out  →  ∇·u = 0  (solenoidal)
    ```
21. `nb.md` — **The maths, step by step:** "1. Contract ∇'s index with u's: $\nabla\cdot\mathbf u = \partial u_i/\partial
    x_i = \partial u_1/\partial x_1 + \partial u_2/\partial x_2 + \partial u_3/\partial x_3$ (2.23) — a scalar. 2. It is
    the trace of the velocity gradient $G_{ij} = \partial u_i/\partial x_j$ (N52): $\nabla\cdot\mathbf u = G_{ii}$, hence
    an invariant (C07). 3. Ex. 2.3: $\mathbf u = a\mathbf x$: $\partial(ax_i)/\partial x_i = a\,\delta_{ii} = 3a$ (the
    index route: $\partial x_i/\partial x_j = \delta_{ij}$ — gloss: each coordinate depends only on itself); the book's
    third term reads $ax_3\mathbf e_3$ (its 'ax₂e₂' twice is a misprint). $\mathbf u = \mathbf b\times\mathbf x$:
    $\partial(\varepsilon_{lmk}b_lx_m)/\partial x_k = \varepsilon_{lmk}b_l\delta_{mk} = \varepsilon_{lkk}b_l = 0$."
22. `nb.note` — **N52 [B]** ⚠️ "**Two index-order conventions to fix once.** The *gradient of a vector* raises the
    order: `vector_gradient(u)[i, j] = ∂u_i/∂x_j` (Kundu's Ch. 3 usage, $S_{ij} = \tfrac12(\partial u_i/\partial x_j +
    \partial u_j/\partial x_i)$). The *divergence of a tensor* lowers it and contracts the **second** index, $(\nabla
    \cdot\boldsymbol\tau)_i = \partial\tau_{ij}/\partial x_j$ — the opposite index from Cauchy's $f_i = \tau_{ji}n_j$;
    the two agree only for symmetric τ (Ch. 4's momentum equation uses $\partial\tau_{ij}/\partial x_j$).
    `tensor_divergence(T, h, index=1)` exposes the choice." equation `(\nabla\cdot\boldsymbol\tau)_i = \sum_{j=1}^{3}
    \frac{\partial\tau_{ij}}{\partial x_j} \equiv \frac{\partial\tau_{ij}}{\partial x_j}`, ref "§2.9". (N74 is placed in
    §2.14 with a back-reference to this block.)
23. `nb.worked_example("Ex. 2.3 with a = 1: u = x", "1. $u = (x_1, x_2, x_3)$. 2. $\partial u_1/\partial x_1 = 1$,
    likewise for 2 and 3. 3. $\nabla\cdot\mathbf u = 1 + 1 + 1 = 3 = 3a$ ✓ — every direction stretches at unit rate.
    4. On the grid with h = 0.1: central differences of a linear field are exact, so the code prints 3.000 at every
    point.")`
24. `nb.code[explain]` — *code:* `g3 = ch02.grid(((-1, 1),)*3, 24 if FAST else 48)`; `u_rad = ch02.radial_field(1.0)
    (g3.X, g3.Y, g3.Z)` (shape (3, n, n, n)); `div = ch02.divergence(u_rad, g3.h); print(div.min(), div.max())`; `u_rot
    = ch02.solid_body_rotation_field([0, 0, 1.0])(g3.X, g3.Y, g3.Z); print(np.abs(ch02.divergence(u_rot, g3.h)).max(),
    ch02.is_solenoidal(u_rot, g3.h))`; `G = ch02.vector_gradient(u_rad, g3.h); print(G[:, :, 5, 5, 5])`; `print
    (np.allclose(np.einsum('ii...', G), div))`; sympy twin: `X = sp.Matrix(sp.symbols('x1 x2 x3')); print(ch02.
    exact_div_curl(sp.Symbol('a')*X))`. *expect:* 3.0 3.0 (to 1e-12) · ~1e-13 True · the identity matrix · True ·
    (3a, [0, 0, 0]). *explain:* 1. a 3-D grid; 2. the radial field evaluated by broadcasting; 3. its divergence — 3
    everywhere; 4. the rotation field is solenoidal; 5. the velocity gradient of u = x is δ_ij; 6. divergence = trace
    of G; 7. sympy agrees symbolically.
25. `nb.check_agree` — **from scratch:** `div_mine = ch02.partial(u_rad[0], 0, g3.h[0]) + ch02.partial(u_rad[1], 1,
    g3.h[1]) + ch02.partial(u_rad[2], 2, g3.h[2])` → `assert np.allclose(div_mine, div)`. Markdown: "component c is
    differentiated along coordinate direction c; inside `partial` direction 0 (x) is array axis 2 — the layout rule in
    action."
26. `nb.figure` — heatmaps (2-D slices z = 0) of ∇·u for u = a x (uniform 3, orange) and for u = b × x (0, white),
    with the quiver of each field on top. *see:* "spreading arrows and a uniform orange square; circling arrows and a
    blank square." *read:* "divergence measures spreading, not speed — the rotating field is fast at the rim yet has
    zero divergence." *change:* "u = (x, −y, 0): stretch in x, squeeze in y, divergence 0 — solenoidal without
    rotating."
27. `nb.animation` — curation §6 row 2 (`player="video"`, 40 frames; FAST 24): tracer particles in a x (a ring of dots
    expanding from the origin) and in b × x (the ring turning rigidly) side by side, titles printing "∇·u = 3, ∇×u = 0"
    and "∇·u = 0, ∇×u = 2b" (uses `linear_flow_map` for exact positions x(t) = e^{Gt}x₀ with G = aI and G = [b×]).
    *see:* "one ring grows, the other spins." *read:* "divergence is the rate the ring's area grows (2a per unit area
    in 2-D); curl is twice the ring's spin rate." *change:* "a < 0: the ring shrinks — a sink."
28. `nb.md` — **What would change if…** "…you crossed ∇ with u instead of dotting? The curl — spin instead of spread."
29. `nb.core("C11", "The curl ∇×u = ε_ijk ∂u_k/∂x_j", question="How much does a velocity field *spin* at a point — and
    why does a perfectly straight shear flow have curl?")`
30. `nb.md` — **The problem in plain words:** "Vorticity — the curl of the velocity — is the subject of Ch. 5 and half
    of Ch. 13: hurricanes, ocean gyres and the planetary vorticity $2\boldsymbol\Omega$ are all $\nabla\times\mathbf u$.
    A paddle wheel dropped into the flow turns at half the curl. It turns in a whirlpool — and also in a straight
    river whose speed varies across the channel."
31. `nb.md` — **The idea:**
    ```
    (∇×u)_3 = ∂u_2/∂x_1 − ∂u_1/∂x_2 :  does u_2 grow to the right?  (+)   does u_1 grow upward?  (−)
    a paddle wheel spins if the flow on one side is faster than on the other   →  shear flow has curl
    ε_ijk ∂_j u_k = "cross ∇ with u" with the operator on the left: ∂ acts on u_k, never the other way round
    ```
32. `nb.md` — **The maths, step by step:** "1. (2.21) with $u_i \to \partial/\partial x_j$: $(\nabla\times\mathbf u)_i =
    \varepsilon_{ijk}\,\partial u_k/\partial x_j$ (2.24), operator ordering: ∂ acts on u. 2. D12 below expands it into
    the three components (2.25). 3. Ex. 2.3: $\mathbf u = a\mathbf x$: $\varepsilon_{ijk}\partial(ax_k)/\partial x_j =
    a\varepsilon_{ijk}\delta_{jk} = a\varepsilon_{ijj} = 0$ — irrotational. $\mathbf u = \mathbf b\times\mathbf x$:
    $\varepsilon_{ijk}\partial_j(\varepsilon_{lmk}b_lx_m) = \varepsilon_{ijk}\varepsilon_{lmk}b_l\delta_{mj} = \varepsilon_
    {ijk}\varepsilon_{ljk}b_l = 2\delta_{il}b_l = 2b_i$ (the 2δ contraction of D09): **the curl of a solid-body rotation
    is twice its angular velocity** — Ch. 3's vorticity = 2Ω."
33. `nb.derivation("D12", "The three components of the curl", ref="2.25", …)` — Part F · D12 (5 steps).
34. `nb.note` — **N53 [B]** (2.25) stated as D12's result with `curl_components` asserted equal to `curl`: equation
    `(\nabla\times\mathbf u)_1 = \frac{\partial u_3}{\partial x_2} - \frac{\partial u_2}{\partial x_3}, \quad (\nabla
    \times\mathbf u)_2 = \frac{\partial u_1}{\partial x_3} - \frac{\partial u_3}{\partial x_1}, \quad (\nabla\times
    \mathbf u)_3 = \frac{\partial u_2}{\partial x_1} - \frac{\partial u_1}{\partial x_2}`, ref "2.25". · **N54 [B]**
    "**Solenoidal** (∇·u = 0, 'no magnetic monopoles' in the word's origin) and **irrotational** (∇×u = 0). Incompressible
    flow (Ch. 4) is solenoidal; potential flow (Ch. 6) is both; Ch. 3 explains why 'irrotational' means no local spin.
    `is_solenoidal`, `is_irrotational`." · **N55 [B]** Ex. 2.3 stated (curation §4c D13): "$a\mathbf x$: divergence 3a,
    curl 0; $\mathbf b\times\mathbf x$: divergence 0, curl 2b — the worked numbers of C10 and here, by hand above and
    by sympy below; pointer: vorticity = 2Ω (Ch. 3)."
35. `nb.worked_example("Ex. 2.3 with b = e₃: u = (−x₂, x₁, 0)", "1. $(\nabla\times\mathbf u)_3 = \partial u_2/\partial
    x_1 - \partial u_1/\partial x_2 = 1 - (-1) = 2$. 2. $(\nabla\times\mathbf u)_1 = \partial u_3/\partial x_2 -
    \partial u_2/\partial x_3 = 0 - 0 = 0$, likewise component 2. 3. $\nabla\times\mathbf u = (0, 0, 2) = 2\mathbf b$ ✓.
    4. A straight shear flow $u_1 = \Gamma x_2$: $(\nabla\times\mathbf u)_3 = 0 - \Gamma = -\Gamma$ — curl without a
    single curved streamline (the paddle wheel turns clockwise for Γ > 0).")`
36. `nb.code[explain]` — *code:* `w = ch02.curl(u_rot, g3.h); print(w[:, 5, 5, 5], np.abs(w[2] - 2).max())`; `print
    (np.allclose(np.stack(ch02.curl_components(u_rot, g3.h)), w))`; `print(ch02.is_irrotational(u_rad, g3.h), ch02.
    is_irrotational(u_rot, g3.h))`; `print(ch02.exact_div_curl(sp.Matrix([0, 0, 1]).cross(X)))`; shear: `u_sh =
    ch02.shear_field(1.0)(g3.X, g3.Y, g3.Z); print(ch02.curl(u_sh, g3.h)[2, 5, 5, 5])`. *expect:* [0 0 2] ~1e-13 · True
    · True False · (0, [0, 0, 2]) · −1.0. *explain:* 1. the curl of the rotation field: (0, 0, 2b) everywhere; 2. the
    three explicit components agree with the einsum form; 3. the radial field is irrotational, the rotating one is not;
    4. sympy twin; 5. the straight shear flow has curl −Γ.
37. `nb.check_agree` — **from scratch:** `dx = lambda f, c: ch02.partial(f, c, g3.h[c])` (coordinate direction c); `w_mine =
    np.stack([dx(u_rot[2], 1) - dx(u_rot[1], 2), dx(u_rot[0], 2) - dx(u_rot[2], 0), dx(u_rot[1], 0) - dx(u_rot[0],
    1)])` → `assert np.allclose(w_mine, w)`.
38. `nb.figure` — heatmaps of $(\nabla\times\mathbf u)_3$ (purple/green) for b × x (uniform 2), for the shear flow
    (uniform −Γ) and for a x (0), each with its quiver and a small paddle-wheel glyph whose arrow marks the spin sense.
    *see:* "the straight shear flow is as purple/green as the whirlpool." *read:* "curl is the *difference* of speeds
    across the wheel, not curvature of streamlines." *change:* "the irrotational vortex $u_\theta = K/r$ (E5): curl 0
    everywhere except the centre — a whirlpool with no local spin."
39. `nb.explainer("stokes_circulation_loop", heading="How much does the flow go round a loop?", why="Curl and
    circulation are two views of one thing; dragging a loop over a shear flow (curl without curves) and over an
    irrotational vortex (curves without curl) removes the confusion no static figure can.", tries=["Shear preset: the
    streamlines are straight, yet the wheel turns and Γ ≠ 0.", "Irrotational-vortex preset: shrink the loop away from
    the centre — Γ/A → 0; enclose the centre — Γ = 2πK whatever the size (the ⚠️ status).", "Flip the orientation: both
    sides of (2.34) change sign together.", "Rectangle mode: read the four side sums of Ex. 2.6 in the Explain tab."])`
    — embedded here for C11's aha (C16 refers back to it; the embed happens once).
40. `nb.md` — **What would change if…** "…you split the velocity gradient $G_{ij}$ into its symmetric and antisymmetric
    halves? The antisymmetric half *is* the curl, packed into a matrix — C12."

### A.8 §2.10 Symmetric and Antisymmetric Tensors — C12 (+N56–N61 · D14, D15 · E3)
1. `nb.section("2.10", "Symmetric and Antisymmetric Tensors", intro="**What is this section about?** Every second-order
   tensor splits uniquely into a symmetric part (6 numbers) and an antisymmetric part (3 numbers). The antisymmetric
   part is a vector in disguise, and a symmetric tensor cannot see it — the two facts Ch. 3 (strain vs rotation) and
   Ch. 4 (dissipation) are built on.")`
2. `nb.core("C12", "Every tensor is a symmetric plus an antisymmetric part — and the antisymmetric part is a vector",
   question="A fluid element both stretches and spins. How does the velocity gradient split those two motions apart —
   and why does the spin hide inside three numbers?")`
3. `nb.md` — **The problem in plain words:** "Ch. 3 will take the velocity gradient $G_{ij} = \partial u_i/\partial x_j$
   (nine numbers) and ask: how fast is the fluid element being deformed, and how fast is it turning? The answer is a
   split: $G = S + A$, S symmetric (the strain rate — deformation), A antisymmetric (the rotation rate — spin). Ch. 4
   then shows viscous friction costs energy only through S. So the split is not algebra for its own sake; it separates
   the motion that heats the fluid from the motion that does not."
4. `nb.md` — **The idea:**
   ```
   B = ½(B + Bᵀ) + ½(B − Bᵀ) = S + A       S_ij = S_ji (6 free numbers)     A_ij = −A_ji, A_ii = 0 (3 free numbers)
   A =  [  0   −ω₃   ω₂ ]                  three numbers → one vector ω ;  A·x = ω × x  (A rotates x about ω)
        [  ω₃   0   −ω₁ ]
        [ −ω₂   ω₁   0  ]
   simple shear u₁ = Γx₂ :  G = [[0, Γ],[0, 0]] = ½[[0, Γ],[Γ, 0]] + ½[[0, Γ],[−Γ, 0]]  = stretch at 45° + spin at Γ/2
   ```
5. `nb.primer("scipy.linalg.expm", "For a linear velocity field u = G·x the trajectory is $\mathbf x(t) = e^{\mathbf
   Gt}\mathbf x_0$, where the matrix exponential $e^{\mathbf Gt} = \mathbf I + \mathbf Gt + (\mathbf Gt)^2/2 + \dots$.
   `scipy.linalg.expm(G*t)` computes it; for a pure rotation matrix it returns the rotation by angle |ω|t.", code="from
   scipy.linalg import expm\nG_rot = np.array([[0., -1.], [1., 0.]])          # u = (−y, x): solid-body rotation,
   ω = 1\nprint(np.round(expm(G_rot*np.pi/2), 12))          # rotation by 90°: [[0, −1],[1, 0]]")`
6. `nb.md` — **The maths, step by step:** "1. Add and subtract ½B_ji: $B_{ij} = \tfrac12(B_{ij} + B_{ji}) + \tfrac12
   (B_{ij} - B_{ji}) = S_{ij} + A_{ij}$ — D14 shows the split is unique and frame-independent. 2. Counting: a symmetric
   3×3 has 3 diagonal + 3 off-diagonal = 6 independent entries; an antisymmetric one has zero diagonal and 3 (N56).
   3. Pack a vector into A: $R_{ij} = -\varepsilon_{ijk}\omega_k$ (2.26)–(2.27) and unpack it: $\omega_k = -\tfrac12
   \varepsilon_{ijk}R_{ij}$ — D15, which also shows $\mathbf R\cdot\mathbf x = \boldsymbol\omega\times\mathbf x$. 4. A
   symmetric τ against any B: $\tau_{ij}B_{ij} = \tau_{ij}S_{ij}$, because $\tau_{ij}A_{ij} = 0$ (N59–N61)."
7. `nb.note` — **N56 [B]** "**Definitions**: symmetric $B_{ij} = B_{ji}$ (6 independent), antisymmetric $B_{ij} =
   -B_{ji}$ (zero diagonal, 3 independent); `independent_components` → 6 / 3 / 9."
8. `nb.derivation("D14", "The unique split into symmetric and antisymmetric parts", ref="§2.10", …)` — Part F · D14
   (6 steps).
9. `nb.note` — **N57 [B]** (2.26): "**The antisymmetric tensor of a vector ω**" equation `\mathbf R = \begin{bmatrix}0 &
   -\omega_3 & \omega_2\\ \omega_3 & 0 & -\omega_1\\ -\omega_2 & \omega_1 & 0\end{bmatrix}`, ref "2.26";
   `antisymmetric_from_vector(omega)`; pointer: R is Ch. 3's rotation tensor of the vorticity. · **N58 [B]** (2.27):
   equation `R_{ij} = -\varepsilon_{ijk}\omega_k, \qquad \omega_k = -\tfrac12\,\varepsilon_{ijk}R_{ij}`, ref "2.27" —
   "the book prints the lower limit of the first sum as 'i−1'; it is i = 1. Round trip in code to 1e-16; and $\mathbf R
   \cdot\mathbf x = \boldsymbol\omega\times\mathbf x$ (D15)."
10. `nb.derivation("D15", "The vector hidden in an antisymmetric tensor", ref="2.26, 2.27", …)` — Part F · D15 (7 steps).
11. `nb.note` — **N59 [B]** (2.28): equation `P = \tau_{kl}B_{kl} = \tau_{kl}(S_{kl} + A_{kl}) = \tau_{ij}S_{ij} +
    \tau_{ij}A_{ij}`, ref "2.28" — "`symmetric_double_contraction(τ, B)` returns (P, P_S, P_A)." · **N60 [B]** (2.29)
    stated (curation §4c D16): "swapping A's indices and renaming dummies gives $P = \tau_{kl}S_{kl} - \tau_{kl}A_{kl}$;
    comparing with (2.28), $X = \tau_{ij}A_{ij}$ satisfies $X = -X$, so $X = 0$." equation `P = \tau_{kl}S_{kl} -
    \tau_{kl}A_{kl}`, ref "2.29". · **N61 [C]** "Hence $\tau_{ij}B_{ij} = \tau_{ij}S_{ij} = \tfrac12\tau_{ij}(B_{ij} +
    B_{ji})$ — like the integral of even × odd over a symmetric interval vanishing; Ch. 4 §4.5: the dissipation
    $\tau_{ij}\partial u_i/\partial x_j$ sees only the strain rate."
12. `nb.worked_example("simple shear u₁ = Γx₂ with Γ = 2 s⁻¹", "1. $G = \begin{bmatrix}0 & 2\\ 0 & 0\end{bmatrix}$
    s⁻¹ ($G_{12} = \partial u_1/\partial x_2$). 2. $S = \tfrac12(G + G^{\rm T}) = \begin{bmatrix}0 & 1\\ 1 &
    0\end{bmatrix}$, $A = \tfrac12(G - G^{\rm T}) = \begin{bmatrix}0 & 1\\ -1 & 0\end{bmatrix}$; S + A = G ✓. 3. The
    vector of A (3-D, $A_{12} = -a_3$): $a_3 = -A_{12} = -1$ s⁻¹ — clockwise spin at 1 rad/s, i.e. $\mathbf a = \tfrac12\nabla\times\mathbf u = \tfrac12\boldsymbol\omega$
    (C11 gave curl = ω₃ = −Γ = −2; the book's R = 2A has vector ω itself) ✓. 4. Check $\mathbf A\cdot\mathbf x = \mathbf a\times\mathbf x$ at x = (1, 0, 0):
    A·x = (0, −1, 0); a × x = (0, 0, −1) × (1, 0, 0) = (0·0 − (−1)·0, (−1)·1 − 0·0, 0) = (0, −1, 0) ✓. 5. S:S = 2, A:A
    (Frobenius) = 2, S:A = 0 — equal parts stretch and spin, orthogonal to each other.")`
13. `nb.code[explain]` — *code:* `G = ch02.velocity_gradient_preset("simple_shear", Gamma=2.0, dim=3)` ([[0, 2, 0],[0, 0, 0],[0,
    0, 0]]); `S, A_ = ch02.symmetric_part(G), ch02.antisymmetric_part(G); print(S); print(A_)`; `assert np.allclose(S +
    A_, G)`; `omega = ch02.vector_from_antisymmetric(A_); print(omega)`; `print(np.allclose(ch02.antisymmetric_from_
    vector(omega), A_))`; `x = np.array([1., 0., 0.]); print(A_ @ x, np.cross(omega, x))`; `print(ch02.independent_
    components(S), ch02.independent_components(A_), ch02.independent_components(G))`; `print(ch02.symmetric_double_
    contraction(S, G))`; cross-check with C11: `Gnum = ch02.vector_gradient(u_sh, g3.h)[:, :, 5, 5, 5]; print(ch02.
    vector_from_antisymmetric(ch02.antisymmetric_part(Gnum)), 0.5*ch02.curl(u_sh, g3.h)[:, 5, 5, 5])` (Γ = 1 there:
    both (0, 0, −0.5)). *expect:* S, A as in the example (3-D padded) · [0 0 −1] · True · [0 −1 0] [0 −1 0] · 6 3 9 ·
    (2.0, 2.0, 0.0) · [0 0 −0.5] [0 0 −0.5]. *explain:* 1. the preset G; 2–3. the two parts and their sum; 4. the vector
    of A by (2.27); 5. the inverse map; 6. A·x = ω × x; 7. component counts; 8. (2.28)–(2.29): P = P_S, P_A = 0; 9.
    the vector of A of a *numerical* velocity gradient is half the curl — the sign convention pinned (G[i, j] =
    ∂u_i/∂x_j).
14. `nb.check_agree` — **from scratch:** `S_mine, A_mine = 0.5*(G + G.T), 0.5*(G - G.T)`; `om_mine = np.array([A_mine
    [2, 1], A_mine[0, 2], A_mine[1, 0]])` (reading (2.26): ω₁ = R₃₂, ω₂ = R₁₃, ω₃ = R₂₁) → `assert np.allclose(S_mine,
    S)`, `assert np.allclose(om_mine, omega)`.
15. `nb.animation` — curation §6 row 3 (`player="video"`, 48 frames, FAST 24): three panels, a 10×10 square of tracers
    under G, S and A from the same start for the simple-shear preset (Γ = 2, t up to 0.8 s): under G the square leans;
    under S it stretches along +45° and squeezes along −45° (eigen-axes drawn, N62/N63 preview); under A it turns
    rigidly clockwise at 1 rad/s (a paddle wheel glyph). `deform_square(G, t)` etc. *see:* "lean = stretch + spin."
    *read:* "the S panel keeps the axes at 45° fixed and changes lengths; the A panel keeps lengths and changes angles."
    *change:* "solid-body preset: the S panel does nothing at all (S = 0)."
16. `nb.explainer("strain_vs_rotation_split", heading="Is simple shear a rotation?", why="The decomposition is a
    statement about motion: only watching the same square under G, S and A on one clock shows that shear is half
    stretch, half spin — and the term bars make τ:A = 0 a thing you see.", tries=["Simple-shear preset: play; the
    stretch factor e^{Γt/2} and the turned angle Γt/2 appear on the end card.", "Set G₁₂ = G₂₁: A vanishes, the wheel
    stops.", "Solid-body preset: S = 0, the square keeps its shape.", "Watch the S:A bar stay at zero while you drag any
    slider."])`
17. `nb.md` — **What would change if…** "…you looked for the directions in which S *only* stretches, with no shear?
    Those are its eigenvectors — C13."

### A.9 §2.11 Eigenvalues and Eigenvectors of a Symmetric Tensor — C13 (+N62, N63 · D17 ★★★)
1. `nb.section("2.11", "Eigenvalues and Eigenvectors of a Symmetric Tensor", intro="**What is this section about?**
   For a real symmetric tensor there are three perpendicular directions on which it acts as a pure stretch; in those
   axes it is diagonal, and its three eigenvalues bound the normal stress (or strain rate) on every plane. The book
   lists the facts; we prove them (D17) and then read Ex. 2.4 by hand.")`
2. `nb.core("C13", "Principal axes of a symmetric tensor", question="Which planes through a point carry *no* shear —
   and why are the normal stresses on them the largest and smallest of all?")`
3. `nb.md` — **The problem in plain words:** "A material fails on the plane where the shear is largest; a fluid
   element stretches fastest along one line and shrinks along another. E2 let you *find* the shear-free planes by
   dragging n. The eigenvalue problem finds them exactly, for any symmetric tensor, and tells you the extreme normal
   stresses without searching."
4. `nb.md` — **The idea:**
   ```
   ask: is there a plane whose traction is purely normal?   τ·b = λ b   (f parallel to n = b)
   → (τ − λδ)·b = 0 has a nonzero b only if det(τ − λδ) = 0 → cubic in λ → three roots λ¹ λ² λ³ (real!)
   → three perpendicular b's ;  axes along them:  τ' = diag(λ¹, λ², λ³) ;  n·τ·n ∈ [λ_min, λ_max] for every n
   ```
5. `nb.primer("eigenvalues and eigenvectors", "$\mathbf A\cdot\mathbf b = \lambda\mathbf b$: the matrix only *scales*
   the vector b by λ. Nonzero b exist exactly when $\det(\mathbf A - \lambda\mathbf I) = 0$ (P53: a zero determinant
   means dependent columns), a polynomial in λ of degree n. `np.linalg.eigh` solves the symmetric case (returns λ
   ascending and unit eigenvectors as columns); `np.roots` finds polynomial roots.", code="A = np.array([[2., 1.], [1.,
   2.]])\nlam, B = np.linalg.eigh(A)\nprint(lam, np.round(B, 4))                 # [1. 3.], columns (1,−1)/√2, (1,1)/√2\n
   print(np.allclose(A @ B[:, 1], lam[1]*B[:, 1]))  # A·b = λb\nprint(np.roots([1, -4, 3]))                  # the
   characteristic polynomial λ² − 4λ + 3")`
6. `nb.primer("complex conjugate", "For $z = a + ib$ the conjugate is $\bar z = a - ib$; $z\bar z = a^2 + b^2 = |z|^2
   \ge 0$, and $z = \bar z$ exactly when z is real. Conjugating a product conjugates each factor; a real number is its
   own conjugate. (Ch. 1 P45 introduced i² = −1.)", code="z = 3 + 4j\nprint(z.conjugate(), (z*z.conjugate()).real,
   abs(z)**2)   # (3-4j) 25.0 25.0")`
7. `nb.primer("quadratic form and the Rayleigh quotient", "$n_i\tau_{ij}n_j = \mathbf n\cdot\boldsymbol\tau\cdot\mathbf n$
   is a quadratic form — for a unit n it is the *normal* component of the traction on the plane ⊥ n. Its values over
   all unit n range between the smallest and largest eigenvalue (the Rayleigh quotient bound, proved in D17). Repeated
   eigenvalues: any orthonormal pair in the eigenplane works (Gram–Schmidt picks one; `eigh` does it for you) — a
   gloss.", code="tau = np.array([[0., 1.], [1., 0.]])\nfor deg in (0, 30, 45, 90):\n    n = np.array([np.cos(np.deg2rad
   (deg)), np.sin(np.deg2rad(deg))])\n    print(deg, np.round(n @ tau @ n, 3))   # 0.0 0.866 1.0 0.0 — never outside
   [−1, 1]")`
8. `nb.md` — **The maths, step by step:** "1. Eigen-equation $\tau_{ij}b_j = \lambda b_i$; nonzero b needs $\det|\tau_{ij}
   - \lambda\delta_{ij}| = 0$, the cubic $\lambda^3 - I_1\lambda^2 + I_2\lambda - I_3 = 0$ (C07, D18). 2. Facts (D17):
   (1) the three λ are real; (2) eigenvectors of distinct λ are orthogonal; (3) with $C = [\mathbf b^1\,\mathbf b^2\,
   \mathbf b^3]$ (columns = new axes, exactly C02's C), $\boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau\mathbf C =
   \mathrm{diag}(\lambda^1, \lambda^2, \lambda^3)$; (4) for every unit n, $\lambda_{\min} \le n_i\tau_{ij}n_j \le
   \lambda_{\max}$ and the shear on any plane is at most $(\lambda_{\max} - \lambda_{\min})/2$. 3. ⚠️ The book's fact
   (4) says the λ bound 'the elements τ_ij'; the sharp and useful statement is about the *normal stress on every
   plane* and the shear bound (Mohr's circle, Ch. 4). 4. 2-D shortcut: principal angle $\tfrac12\arctan2(2S_{12}, S_{11}
   - S_{22})$, principal values $\tfrac12(\tau_{11} + \tau_{22}) \pm \sqrt{(\tfrac12(\tau_{11} - \tau_{22}))^2 +
   \tau_{12}^2}$."
9. `nb.derivation("D17", "Real eigenvalues, orthogonal axes, diagonal form and the bounds — for a real symmetric
   tensor", ref="§2.11 facts (1)–(4)", …, check_src=…)` — Part F · D17 (★★★, 15 steps + sympy check cell).
10. `nb.note` — **N62 [B]** Ex. 2.4 stated (curation §4c D19) and traced in the tiny example below; ⚠️ "The book writes
    '2S₁₂ = du₁/dx₂ = Γ' *and* S = [[0, Γ],[Γ, 0]]. Those disagree by a factor 2 (for u₁(x₂) alone, S₁₂ = ½ du₁/dx₂).
    Every later line of the example uses S₁₂ = Γ, so we take **Γ ≡ S₁₂** (`example_2_4(Gamma)`); with du₁/dx₂ = 2 s⁻¹
    as in C12, Γ = 1 s⁻¹. Also 'det|E_ij − λδ_ij|' means S_ij." Interpretation: "stretching at rate Γ along b¹ (45°),
    compression at −Γ along b² (135°) — the S panel of C12's animation; Ch. 3 §3.4 makes this the definition of strain
    rate." · **N63 [B]** "Fig. 2.8 (the 45° axes) is our slider figure below."
11. `nb.worked_example("Ex. 2.4 by hand with Γ = 1 s⁻¹", "1. $S = \begin{bmatrix}0 & 1\\ 1 & 0\end{bmatrix}$. 2.
    $\det\begin{bmatrix}-\lambda & 1\\ 1 & -\lambda\end{bmatrix} = \lambda^2 - 1 = 0 \Rightarrow \lambda^1 = 1,
    \lambda^2 = -1$ (Vieta: sum 0 = trace ✓, product −1 = det ✓). 3. Eigenvector for λ = 1: first row of (S − λI)·b = 0
    gives $-b_1 + b_2 = 0 \Rightarrow b_1 = b_2$; normalise $b_1^2 + b_2^2 = 1 \Rightarrow \mathbf b^1 = (1, 1)/\sqrt2$.
    4. For λ = −1: $b_1 + b_2 = 0 \Rightarrow \mathbf b^2 = (-1, 1)/\sqrt2$ (sign chosen so that $\det[\mathbf b^1\,
    \mathbf b^2] = +1$: a 45° rotation, not a mirror). 5. $C = \tfrac1{\sqrt2}\begin{bmatrix}1 & -1\\ 1 & 1\end{bmatrix}$
    — C02's 2-D matrix with θ = 45° ✓. 6. $S'_{12} = C_{i1}C_{j2}S_{ij} = C_{11}C_{22}S_{12} + C_{21}C_{12}S_{21} =
    \tfrac12 - \tfrac12 = 0$; $S'_{11} = 2C_{11}C_{21}S_{12} = 1$; $S'_{22} = 2C_{12}C_{22}S_{12} = -1$: $S' = \mathrm
    {diag}(1, -1)$ ✓. 7. Bounds: on the 30° plane $n\cdot S\cdot n = \sin60° = 0.866 \in [-1, 1]$ ✓; max shear
    $(1 - (-1))/2 = 1$, reached at φ = 0° and 90° (the original axes, where S has only shear) ✓.")`
12. `nb.code[explain]` — *code:* `res = ch02.example_2_4(1.0); print(res["lam"], np.round(res["B"], 4), np.round(np.
    rad2deg(res["angle_rad"]), 1)); print(np.round(res["S_prime"], 12))`; `lam, B = ch02.principal_axes(A) (EXAMPLE_
    TENSOR); print(np.round(lam, 4), np.round(np.linalg.det(B), 12))`; `C, tp = ch02.diagonalize(A); print(np.round(tp,
    10))`; `print(np.round(np.roots(ch02.characteristic_polynomial(A)), 4))`; `print(ch02.normal_stress_bounds(A,
    np.random.default_rng(6)))`; `print(ch02.principal_angle_2d([[0, 1], [1, 0]]))`; `try: ch02.principal_axes([[0, 1],
    [0, 0]]) except ValueError as e: print(e)`. *expect:* [1. −1.] [[0.7071 −0.7071],[0.7071 0.7071]] 45.0 · diag(1, −1)
    · [1.2679 3. 4.7321] 1.0 · diag(1.2679, 3, 4.7321) · the same three roots · (≈1.27, ≈4.73, ≈1.73) (Monte-Carlo:
    within [λ_min, λ_max], max shear ≤ (4.7321 − 1.2679)/2 = 1.732) · 0.7854 · "principal_axes needs a symmetric
    tensor". *explain:* 1. Ex. 2.4 packaged (book order (Γ, −Γ), b¹ first); 2. the 3×3 example's principal values by
    `eigh`, det C = +1; 3. τ' diagonal; 4. the same λ from the characteristic polynomial (D18); 5. 1000 random planes
    never leave the bounds; 6. the 2-D angle; 7. non-symmetric input is refused — facts (1)–(4) do not hold for it.
13. `nb.check_agree` — **from scratch (curation §7):** for `S2 = [[0, 1],[1, 0]]`: `tr, det = S2.trace(), np.linalg.det
    (S2); lam_mine = np.array([(tr + np.sqrt(tr**2 - 4*det))/2, (tr - np.sqrt(tr**2 - 4*det))/2])`; `b1 = np.array
    ([-(S2[0, 0] - lam_mine[0]) ... ])` — simpler: `b1 = np.array([S2[0, 1], lam_mine[0] - S2[0, 0]]); b1 /= np.linalg.
    norm(b1); b2 = np.array([-b1[1], b1[0]])` (perpendicular, det +1) → `assert np.allclose(lam_mine, res["lam"])`,
    `assert np.allclose(np.abs(b1 @ res["B"][:, 0]), 1)`, `assert np.allclose(np.c_[b1, b2].T @ S2 @ np.c_[b1, b2],
    np.diag(lam_mine))`.
14. `nb.plotly` — curation §6 slider (C13): Γ from 0.2 to 3 (30 steps): left, the square deformed by S for a fixed
    time with the eigen-axes at 45° and arrows ±Γ; right, Mohr's circle of S (centre 0, radius Γ) with the point for
    the 30° plane. *see:* "the axes never move; the arrows and the circle scale." *read:* "the principal directions
    depend on the *shape* of S, the eigenvalues on its size." *change:* "add S₁₁ = 1: the angle drops below 45°
    (`principal_angle_2d`) and the circle's centre moves right."
15. `nb.live` — S₁₁, S₂₂, S₁₂ sliders (`live(fn, S11=(-2, 2, 0.1), S22=…, S12=…)`) drawing the element, eigen-axes,
    Mohr's circle and printing λ and the angle (paired with the slider figure above for the page).
16. `nb.md` — **What would change if…** "…τ were not symmetric? Complex eigenvalues and skewed eigenvectors are
    possible and the bounds fail — which is why `principal_axes` refuses, and why Ch. 4's proof that the stress *is*
    symmetric matters."

### A.10 §2.12 Gauss' Theorem — C14 (+N64, N65 · D25 · E4), C15 (+N66–N68, N72 · D21, D22)
1. `nb.section("2.12", "Gauss' Theorem", intro="**What is this section about?** A derivative integrated over a volume
   equals the field itself integrated over the boundary with the outward normal — for a field of any order. Read
   backwards on a tiny box it *defines* the divergence as outflux per unit volume, and the Cartesian formula (2.23)
   falls out face by face.")`
2. `nb.core("C14", "Gauss' theorem: derivative inside = normal outside", question="Why does adding up ∂Q/∂x_i over a
   whole volume only depend on what Q does on the surface?")`
3. `nb.md` — **The problem in plain words:** "Ch. 4 will write conservation laws for a blob of fluid: the mass inside
   changes only by what flows through its skin; the momentum inside changes by the forces on its skin (C05!). Turning
   'inside' into 'on the skin' — and back — is Gauss' theorem. Without it there is no differential form of any
   conservation law."
4. `nb.md` — **The idea:**
   ```
   1-D:  ∫_a^b f'(x) dx = f(b) − f(a)                 (fundamental theorem: the inside adds up to the two ends)
   3-D:  ∭_V ∂Q/∂x_i dV = ∯_A n_i Q dA                 (the ends become the boundary, with n_i = ±1 on two faces of a box)
   any V:  tile it with boxes — every interior face is counted twice with opposite n and cancels
   ```
5. `nb.primer("volume and surface integrals as midpoint sums", "$\iiint_V f\,dV$ adds f × (small volume) over a
   3-D grid of cells; the midpoint rule uses f at each cell centre (error ∝ h²). A surface integral $\iint_A g\,dA$
   does the same on a 2-D grid of the surface. Iterated integrals: integrate in x first, then y, then z.", code="n = 20;
   xs = (np.arange(n) + 0.5)/n                 # cell centres on [0, 1]\nX, Y, Z = np.meshgrid(xs, xs, xs, indexing='ij')\n
   print(np.sum(X**2) / n**3)                   # ∭ x² dV over the unit cube = 1/3 (0.3329 with n = 20)\nprint(np.sum
   (X[:, :, 0]**2) / n**2)             # ∬ x² dA over the unit square = 1/3")`
6. `nb.primer("fundamental theorem of calculus", "$\int_a^b f'(x)\,dx = f(b) - f(a)$: integrating a derivative gives
   the change of the function between the ends (Ch. 1 P27 defined the integral; this is its partner). It is the whole
   content of Gauss' theorem in one dimension.", code="f = lambda x: x**3\nxs = np.linspace(0, 2, 2001)\nprint(np.
   trapezoid(3*xs**2, xs), f(2) - f(0))          # 8.000 8")`
7. `nb.md` — **The maths, step by step:** "1. Statement (2.30): $\iiint_V \partial Q/\partial x_i\,dV = \iint_A n_i Q\,
   dA$ for Q of any order (scalar, vector, tensor — the other indices ride along), n outward. 2. For a vector, set Q →
   Q_i and sum on i (allowed: (2.30) is linear in Q; N64): $\iiint_V \nabla\cdot\mathbf Q\,dV = \iint_A \mathbf n\cdot
   \mathbf Q\,dA$ — the **divergence theorem**: volume integral of the divergence = net outflux. 3. Proof (D25): on a box
   by the fundamental theorem of calculus along one axis; on any V by tiling and cancelling interior faces."
8. `nb.derivation("D25", "Gauss' theorem, from the fundamental theorem of calculus", ref="2.30", …)` — Part F · D25
   (9 steps).
9. `nb.note` — **N64 [B]** equation `\iiint_V \frac{\partial Q_i}{\partial x_i}\,dV = \iint_A n_i Q_i\,dA \quad\text
   {or}\quad \iiint_V \nabla\cdot\mathbf Q\,dV = \iint_A \mathbf n\cdot\mathbf Q\,dA`, ref "§2.12" — "the outflux reading:
   $\mathbf n\cdot\mathbf Q\,dA$ is how much Q leaves through dA. Benchmark (Wikipedia's divergence-theorem example):
   $\mathbf F = (2x, y^2, z^2)$ through the unit sphere: $\nabla\cdot\mathbf F = 2 + 2y + 2z$, whose integral over the
   ball is $2\cdot\tfrac43\pi + 0 + 0 = 8\pi/3 = 8.378$; `divergence_theorem_sphere` gives both sides." · **N65 [B]**
   "Fig. 2.9 is our 3-D box below, faces coloured by n·Q."
10. `nb.worked_example("Q = (x, y, z) on the unit cube [0, 1]³", "1. Inside: $\nabla\cdot\mathbf Q = 1 + 1 + 1 = 3$;
    volume 1 → left side 3. 2. Outside, face by face: on x = 1, n = +e₁, n·Q = x = 1, area 1 → +1; on x = 0, n·Q = −x
    = 0 → 0; same for y and z → right side 1 + 0 + 1 + 0 + 1 + 0 = 3 ✓. 3. Q = (x², 0, 0): inside $\int_0^1 2x\,dx = 1$;
    outside +1 (x = 1) + 0 = 1 ✓ — only the two x-faces contribute, as D25 says.")`
11. `nb.code[explain]` — *code:* `Q = lambda x, y, z: np.stack([x, y, z])`; `lhs, rhs = ch02.divergence_theorem_box(Q,
    ((0, 1),)*3, 16); print(lhs, rhs)`; `Q2 = lambda x, y, z: np.stack([x**2, 0*y, 0*z]); print(ch02.divergence_
    theorem_box(Q2, ((0, 1),)*3, 16))`; `F = lambda x, y, z: np.stack([2*x, y**2, z**2]); lhs_s, rhs_s = ch02.
    divergence_theorem_sphere(F, 1.0, 24 if FAST else 48); print(lhs_s, rhs_s, 8*np.pi/3)`; `print(ch02.gauss_gradient_
    box(lambda x, y, z: x*y*z, ((0, 1),)*3, 16))`; 2-D: `Q2d = ch02.radial_field(1.0, dim=2); print(ch02.flux_through_faces(Q2d, ((0, 1), (0, 1)), 16))`; `print(ch02.divergence_theorem_tiled(Q2d, ((0, 1), (0, 1)), 4, 16))`. *expect:*
    3.0 3.0 · (1.0, 1.0) (exact: midpoint is exact for quadratics? — midpoint integrates x² with error h²/12·… the
    volume term ∫2x is exact; the face term ∫1 exact → 1.0 1.0) · 8.377 8.378 8.378 (quadrature error ~1e-3 at n = 24;
    ~1e-4 at 48) · two equal vectors (0.25, 0.25, 0.25) · {'+x': 1.0, '-x': 0.0, '+y': 1.0, '-y': 0.0} · (2.0,
    2.0, ~1e-15). *explain:* 1. the vector field as a function of the grid arrays; 2. both sides on the cube: 3 = 3;
    3. the x²-field: 1 = 1; 4. the sphere benchmark 8π/3; 5. the scalar form (2.30) gives a vector on both sides; 6.
    the 2-D face fluxes (right and top leak 1 each); 7. tiling 4×4: the tile boundaries sum to the outer flux and the
    interior faces cancel to round-off — D25 step 8 as a number.
12. `nb.check_agree` — **from scratch:** six-face midpoint sums on the unit cube for Q = (x, y, z): `n = 16; c =
    (np.arange(n) + 0.5)/n; A, B = np.meshgrid(c, c, indexing='ij'); flux = 0; for axis in range(3): for side, sgn in
    ((1.0, 1), (0.0, -1)): … flux += sgn * (component `axis` at that face).sum()/n**2`; volume: `X, Y, Z = np.meshgrid
    (c, c, c, indexing='ij'); div = 3*np.ones_like(X); vol = div.sum()/n**3` → `assert np.allclose((vol, flux),
    ch02.divergence_theorem_box(Q, ((0, 1),)*3, n))`.
13. `nb.plotly` — Fig. 2.9 as our 3-D box: `from scripts.ch02_fig2_9_gauss import gauss_box_figure; gauss_box_figure(F)`
    — faces coloured by n·Q (blue in / orange out) with the two totals in the title. *see:* "orange faces where Q
    leaves, blue where it enters." *read:* "the sum of the face colours (weighted by area) equals the integral of the
    divergence inside." *change:* "Q = b × x: every face is half blue, half orange, net zero — solenoidal."
14. `nb.plotly` — curation §6 slider (C14): mesh count n = 4 … 32 (8 steps): bars for ∭∇·Q dV and ∯n·Q dA for the
    sin/cos `smooth_test_field` on the unit cube and their gap on a log axis. *see:* "the two bars meet as the mesh
    refines; the gap falls by 4 per doubling." *read:* "both sides are approximations of one number." *change:* "a
    polynomial Q of degree ≤ 2: the gap is 1e-16 at every n."
15. `nb.explainer("gauss_flux_box", heading="What leaks out of a box?", why="The theorem has a left side that lives
    inside and a right side on the boundary; moving and resizing a box over a field, watching the face fluxes re-sum to
    the integrated divergence, tiling it to see interior faces cancel, and shrinking it until (1/V)∮ settles on ∇·Q is
    the whole (2.30) → (2.32) chain as an experiment.", tries=["Drag the box over the source: the total bar turns
    orange.", "Tile 4×4: interior arrows appear in opposite pairs and cancel — D25 step 8.", "Shrink h on the log slider
    and watch (1/V)∮ approach the dashed ∇·Q(x₀) with slope 2.", "b × x preset: every box, anywhere, has zero net
    flux."])`
16. `nb.md` — **What would change if…** "…the volume shrank to a point? Divide both sides by V and take the limit: the
    right side becomes 'outflux per unit volume' — a definition of the divergence that never mentions coordinates. C15."
17. `nb.core("C15", "Divergence as outflux per unit volume: the integral definitions", question="Can ∇·Q be defined
    without any coordinates at all — and how does the formula ∂Q_i/∂x_i come back out?")`
18. `nb.md` — **The problem in plain words:** "The Cartesian formula (2.23) is a recipe, not a meaning. The meaning is:
    put a tiny closed surface around the point, measure what leaks out, divide by the volume. That definition works in
    spherical coordinates on a planet as well as on a Cartesian grid, and it is exactly the box argument Ch. 4 uses to
    derive the continuity equation."
19. `nb.md` — **The idea:**
    ```
    (2.30) on a tiny box around x₀:  (1/V) ∯ n·Q dA  →  ∇·Q(x₀)  as V → 0          (2.32)
    face by face (Ex. 2.5):  Q on the +e₁ face ≈ Q + (Δx₁/2) ∂Q/∂x₁ ;  on the −e₁ face ≈ Q − (Δx₁/2) ∂Q/∂x₁
                             opposite faces: the Q's cancel, the slopes add  →  ∂Q₁/∂x₁ Δx₁Δx₂Δx₃ ;  divide by V
    ```
20. `nb.primer("mean-value theorem for integrals", "If f is continuous on a region V, then $\iiint_V f\,dV = f(\mathbf
    x^*)\,V$ for some point x* inside V: the integral equals the volume times some intermediate value of f. As V shrinks
    to a point x₀, x* is squeezed onto x₀.", code="xs = np.linspace(0.9, 1.1, 2001)              # a small interval
    around 1\nf = xs**2\nmean = np.trapezoid(f, xs) / 0.2\nprint(mean, np.sqrt(mean))                # 1.0033, and x* =
    1.0017 lies inside [0.9, 1.1]")`
21. `nb.md` — **The maths, step by step:** "1. (2.31) $\mathcal DQ = \lim_{V\to0}(1/V)\iint_A n_iQ\,dA$ — the generalised
    derivative (gradient of any order), from (2.30) by the mean-value theorem (D21). 2. Contract: (2.32) $\nabla\cdot
    \mathbf Q = \lim (1/V)\iint \mathbf n\cdot\mathbf Q\,dA$; cross: (2.33) $\nabla\times\mathbf Q = \lim (1/V)\iint
    \mathbf n\times\mathbf Q\,dA$. 3. Ex. 2.5 (D22): evaluate (2.32) on a box — Taylor to the six face centres, opposite
    faces cancel the zeroth order, divide by the volume, take the limit: (2.23) returns."
22. `nb.derivation("D21", "The integral definitions as small-volume limits of Gauss' theorem", ref="2.31–2.33", …)` —
    Part F · D21 (7 steps).
23. `nb.note` — **N66 [B]** (2.31) equation `\mathcal D Q = \lim_{V\to 0}\frac{1}{V}\iint_A n_i Q\,dA`, ref "2.31" —
    "`integral_gradient(Q_fn, x0, h)` converges to `gradient` with order 2 in h (printed below)." · **N67 [B]** (2.33)
    equation `\nabla\times\mathbf Q = \lim_{V\to 0}\frac{1}{V}\iint_A \mathbf n\times\mathbf Q\,dA`, ref "2.33" —
    "`integral_curl` → `curl`; Ch. 5 uses this picture for vorticity." · **N72 [C]** "Remark (§2.13): the integral
    definitions of ∇, ∇·, ∇× do not depend on the coordinate system (Exercises 2.16–2.18) — cylindrical and spherical
    forms in Ch. 3 / Appendix B." · **N68 [B]** "Ex. 2.5 is the derivation D22 below (its recipe animated after it);
    the same box argument gives the continuity equation in Ch. 4 §4.2."
24. `nb.derivation("D22", "Ex. 2.5: the Cartesian divergence recovered from the outflux definition", ref="2.32 → 2.23",
    …)` — Part F · D22 (9 steps).
25. `nb.worked_example("Q = (x², 0, 0) around x₀ = (1, 0, 0), box side h = 0.2", "1. Faces ⊥ x: at x = 1.1, n·Q = 1.21,
    area h² = 0.04 → +0.0484; at x = 0.9, n·Q = −0.81 → −0.0324. 2. The other four faces: Q₂ = Q₃ = 0 → 0. 3. Net
    outflux 0.0160; volume h³ = 0.008; ratio 2.000. 4. Exact ∇·Q = 2x = 2 at x₀ ✓ — exact even at finite h, because
    $(x+h/2)^2 - (x-h/2)^2 = 2xh$ has no h³ term. 5. For Q = (x³, 0, 0): ratio $3x^2 + h^2/4 = 3.01$ vs exact 3 — the
    error is $h^2/4$: order 2.")`
26. `nb.code[explain]` — *code:* `Qf = lambda x, y, z: np.stack([x**2, 0*y, 0*z])` (WIP callables take `(x, y, z)` arrays); `for h in (0.4, 0.2, 0.1):
    print(h, ch02.integral_divergence(Qf, [1, 0, 0], h))` · `Qc = lambda x, y, z: np.stack([x**3, 0*y, 0*z]); errs =
    [abs(ch02.integral_divergence(Qc, [1, 0, 0], h) - 3.0) for h in hs]; print(errs, observed_order(hs, errs))` ·
    `print(ch02.integral_gradient(lambda x, y, z: x*y*z, [1, 2, 3], 0.1))` (∇(xyz) = (yz, xz, xy) = (6, 3, 2)) ·
    `print(ch02.integral_curl(ch02.solid_body_rotation_field([0, 0, 1.0]), [0.3, 0.2, 0], 0.1))` (→ (0, 0, 2)) ·
    `hs = (0.4, 0.2, 0.1) if FAST else (0.4, 0.2, 0.1, 0.05); st = {k: ch02.integral_definition_convergence(k, hs=hs) for k in
    ("divergence", "gradient", "curl")}; print({k: round(s["order"], 2) for k, s in st.items()})`. *expect:* 2.0 2.0 2.0 · [0.04,
    0.01, 0.0025] 2.00 · [6. 3. 2.] · [0 0 2] · ≈ 2.0 for all three. *explain:* 1. (2.32) evaluated by six face sums on
    a cube; exact for a quadratic; 2. for a cubic the error is h²/4 — observed order 2; 3. (2.31) on a scalar gives the
    gradient; 4. (2.33) gives the curl; 5. the cached integral-definition study.
27. `nb.check_agree` — **from scratch (curation §7):** `def six_faces(Q, x0, h): tot = 0; for ax in range(3): for sgn
    in (+1, -1): p = np.array(x0, float); p[ax] += sgn*h/2; tot += sgn*Q(*p)[ax]*h**2 (midpoint of the face); return
    tot/h**3` → `assert np.allclose(six_faces(Qf, [1, 0, 0], 0.2), ch02.integral_divergence(Qf, [1, 0, 0], 0.2,
    n_face=1))` (one-point face rule = the hand version); and the three-h order: `assert abs(observed_order(hs, [abs
    (six_faces(Qc, …) - 3) for h in hs]) - 2) < 0.15`.
28. `nb.animation` — curation §6 row 4 (`player="frames"`, 12 frames): Ex. 2.5's bookkeeping: frame 1 the box and x₀;
    frames 2–7 the six face values appear one pair at a time (Q(x) ± ½Δx ∂Q/∂x written on each face); frames 8–9 the
    Q(x) terms of opposite faces fade (cancel) leaving the slope terms; frame 10 the sum ∂Q_i/∂x_i V; frames 11–12 the
    box shrinks and (1/V)∮ prints 2.000. *see:* "the constants die in pairs; the slopes survive." *read:* "the
    divergence is the part of the outflux that opposite faces do *not* cancel — the first-order term." *change:* "shift
    x₀ to (2, 0, 0): every face value doubles its slope term; (1/V)∮ → 4."
29. `nb.plotly` — curation §6 slider (C15): box side h (20 log steps from 1 to 0.01): (1/V)∮n·Q dA and (1/V)∮nQ dA
    (gradient) vs the point values, and the error on log axes for the cubic field. *see:* "the estimate glides onto the
    dashed exact line; the error line has slope 2." *read:* "(2.32) is a limit — the approximation gets better as h²."
    *change:* "linear Q: exact at every h."
30. `nb.md` — **What would change if…** "…you did the same with a *loop* instead of a closed surface? Circulation per
    unit area — the curl, by Stokes' theorem. §2.13."

### A.11 §2.13 Stokes' Theorem — C16 (+N69–N71, N73 · D26 · E5)
1. `nb.section("2.13", "Stokes' Theorem", intro="**What is this section about?** The circulation of a velocity field
   around a closed loop equals the flux of its curl through any surface the loop bounds. Read on a tiny loop it defines
   the normal curl as circulation per unit area — and Kelvin's theorem (Ch. 5), lift (Ch. 6, 14) and vortex tubes all
   live on it.")`
2. `nb.core("C16", "Stokes' theorem and circulation", question="Why does the total spin inside a loop equal the flow's
   tendency to run around the loop — and when does that fail?")`
3. `nb.md` — **The problem in plain words:** "Put a closed loop into a flow and add up the velocity component along the
   loop all the way round: that is the circulation. It is the quantity Kelvin's theorem conserves, the quantity that
   produces lift on a wing, the strength of a tornado. Stokes' theorem says it equals the curl integrated over the
   loop's interior — so circulation is 'total spin inside', not 'how curved the streamlines are'."
4. `nb.md` — **The idea:**
   ```
   one small rectangle:  ∮ u·t ds = (∂u₂/∂x₁ − ∂u₁/∂x₂) ΔA = (∇×u)·n ΔA        (Ex. 2.6 bookkeeping)
   tile the surface with rectangles: every interior edge is walked twice, opposite ways → cancels
   what is left:  ∮_C u·t ds = ∬_A (∇×u)·n dA        orientation: thumb along n, fingers along t (right hand)
   ```
5. `nb.primer("line integral of a vector field around a loop", "Parametrise the loop by arc length s (or an angle);
   at each point take the velocity component along the unit tangent, u·t, and add it up times the small length ds:
   $\oint_C \mathbf u\cdot\mathbf t\,ds$. Ch. 1 P35 integrated p dv along a path; this is the vector version. A circle:
   x = c + R(cos θ, sin θ), t = (−sin θ, cos θ), ds = R dθ.", code="th = (np.arange(200) + 0.5) * 2*np.pi/200        #
   midpoints round a unit circle\nx, y = np.cos(th), np.sin(th)\ntx, ty = -np.sin(th), np.cos(th)                   # unit
   tangent, counterclockwise\nu1, u2 = -y, x                                     # u = b × x with b = e3\nprint(np.sum
   ((u1*tx + u2*ty) * 2*np.pi/200))          # 6.2832 = 2π")`
6. `nb.md` — **The maths, step by step:** "1. Orientation (N69): pick the outside of A → n; t runs counterclockwise
   about n (right-hand rule); $\mathbf n_c = \mathbf n\times\mathbf t$ is the in-surface normal to C pointing into A, so
   $\mathbf t = \mathbf n_c\times\mathbf n$ and $(\mathbf n_c, \mathbf n, \mathbf t)$ is right-handed. 2. Statement (2.34)
   $\iint_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA = \oint_C \mathbf u\cdot\mathbf t\,ds$; the right side is the
   **circulation**. 3. Small loop (2.35): $\mathbf n\cdot(\nabla\times\mathbf u) = \lim_{A\to0}(1/A)\oint_C\mathbf u\cdot
   \mathbf t\,ds$ — the normal curl is circulation per unit area (N71). 4. Proof (D26): one rectangle (Ex. 2.6 run
   backwards), tiling, cancellation; corollary $\oint\nabla\phi\cdot\mathbf t\,ds = 0 \Rightarrow \nabla\times\nabla\phi
   = 0$."
7. `nb.note` — **N69 [B]** the orientation rule with a ⚠️: "the book writes $\mathbf t = \mathbf n_c\times\mathbf n$;
   this is right-handed only if $\mathbf n_c$ points along the surface *into* A (up the cap in Fig. 2.10). Our code fixes
   t by the right-hand rule about `normal` and checks $\oint(\mathbf x - \mathbf c)\times\mathbf t\,ds = 2A\,\mathbf n$
   (`planar_loop`, `boundary_tangent`)." · **N70 [B]** "Fig. 2.10 is our 3-D cap below with n, n_c, t at a boundary
   point." · **N71 [B]** (2.35) equation `\mathbf n\cdot(\nabla\times\mathbf u) = \lim_{A\to 0}\frac{1}{A}\oint_C
   \mathbf u\cdot\mathbf t\,ds`, ref "2.35" — "`integral_curl_component(u, x0, n, h)` converges to n·curl with order 2;
   pointer: Ch. 3 vorticity = 2 × angular velocity, Ch. 5." · **N73 [B]** Ex. 2.6 stated correctly (curation §4c D23):
   "for the rectangle Δy × Δz in the plane x = const with n = e_x, the sides at y ± Δy/2 run along ±e_z (integrand
   u_z), the sides at z ∓ Δz/2 run along ±e_y (integrand **u_y** — the book prints u_z there by mistake; its own limit
   $\partial u_z/\partial y - \partial u_y/\partial z$ confirms u_y); midpoint values × side lengths, divide by ΔyΔz,
   take the limit: (2.25)'s first component. The y- and z-components follow by cyclic relabelling. D26 runs the same
   bookkeeping in the x₃-plane."
8. `nb.derivation("D26", "Stokes' theorem for a planar surface: one rectangle, then tiling", ref="2.34", …)` — Part F ·
   D26 (10 steps).
9. `nb.worked_example("u = b × x with b = e₃ around the unit circle", "1. $\mathbf u = (-x_2, x_1, 0)$; on the circle
   $\mathbf x = (\cos\theta, \sin\theta, 0)$, $\mathbf t = (-\sin\theta, \cos\theta, 0)$: $\mathbf u\cdot\mathbf t =
   \sin^2\theta + \cos^2\theta = 1$. 2. Circulation $= \oint 1\,ds = 2\pi R = 2\pi = 6.283$. 3. Curl (C11) = 2b = (0, 0,
   2); flux through the disc $= 2\times\pi R^2 = 2\pi$ ✓. 4. Per unit area: $2\pi/\pi = 2 = (\nabla\times\mathbf u)_3$ ✓
   (2.35). 5. Reverse n → −e₃: t reverses, circulation −2π, flux −2π: both sides flip together.")`
10. `nb.code[explain]` — *code:* `u_fn = ch02.solid_body_rotation_field([0, 0, 1.0])`; `loop = ch02.planar_loop([0, 0,
    0], [0, 0, 1], 1.0, 200); disc = ch02.planar_disc([0, 0, 0], [0, 0, 1], 1.0, 20, 40)`; `print(ch02.circulation(u_fn,
    loop), 2*np.pi)`; `print(ch02.stokes_theorem_check(u_fn, loop, disc, curl_fn=lambda p: np.array([0, 0, 2.0])))`;
    `print(ch02.integral_curl_component(ch02.shear_field(1.0), [0, 0], None, 0.1))` (2-D field, 2-D point; → −1.0) · potential flow:
    `loop2 = ch02.planar_loop([0, 0], None, 1.0, 400)` (a 2-D loop for the 2-D fields); `print(ch02.circulation(ch02.
    potential_field(), loop2))` (default φ = x₁² − x₂²; → 0.0) · irrotational vortex: `uv = ch02.irrotational_vortex_
    field(1.0); off = ch02.planar_loop([2, 0], None, 0.5, 400); print(ch02.circulation(uv, loop2), 2*np.pi, ch02.
    circulation(uv, off))`; `chk = ch02.stokes_theorem_check(uv, loop2, ch02.planar_disc([0, 0], None, 1.0, 20, 40));
    print(chk.hypothesis_ok, chk.note)`. *expect:* 6.2832 6.2832
    · StokesCheck(lhs 6.283, rhs 6.283, hypothesis_ok True) · −1.0 · 0.0 · 6.2832 6.2832 ~1e-16 · False "the field is
    singular at a point inside A: Stokes' theorem does not apply" (probed WIP values). *explain:* 1.
    the loop and the disc (right-handed about e₃); 2. circulation = 2π; 3. both sides of (2.34); 4. Ex. 2.6's recipe on
    the shear flow: −Γ; 5. a gradient field has zero circulation (Exercise 2.20, the D26 corollary); 6. the irrotational
    vortex: 2πK around any loop enclosing the centre, 0 around one that does not; 7. the check flags the singular core
    — Stokes' hypothesis fails, not the theorem.
11. `nb.check_agree` — **from scratch (curation §7):** the midpoint circle sum of the primer cell (`Gamma_mine`) vs
    `ch02.circulation(u_fn, loop)`; a polar-grid sum of $(\nabla\times\mathbf u)_3 = 2$ over the disc: `r = (np.arange
    (20) + 0.5)/20; th = …; flux_mine = np.sum(2 * r[:, None] * (1/20) * (2*np.pi/40))` → `assert np.allclose
    (Gamma_mine, flux_mine, rtol=1e-3)`; `assert np.allclose(Gamma_mine, ch02.circulation(u_fn, loop), rtol=1e-6)`.
12. `nb.plotly` — Fig. 2.10 as ours: `from scripts.ch02_fig2_10_stokes import stokes_cap_figure; stokes_cap_figure()`
    — a hemispherical cap, n at a patch, and at a rim point n (horizontal, outward), n_c (up the cap), t (along the rim,
    counterclockwise seen from above). *see:* "three perpendicular unit vectors at the rim." *read:* "n_c × n = t —
    check with the right hand." *change:* "choose the inside as outside: n and t both reverse; n_c stays."
13. `nb.animation` — curation §6 row 5 (`player="frames"`, 10 frames): a square loop of side h shrinking about (0, 0)
    in the shear field u₁ = Γx₂ (h = 2 → 0.1), the running Γ/A printed against $(\nabla\times\mathbf u)_3 = -\Gamma$,
    the four side sums shown as arrows, a paddle wheel spinning at −Γ/2. *see:* "Γ/A is −Γ at every size here (linear
    field), the sides' contributions shrink together." *read:* "(2.35) is a statement about the ratio, not about the
    loop." *change:* "the sin/cos smooth field: Γ/A wanders at large h and settles at small h with error ∝ h²."
14. `nb.figure` — convergence of `integral_curl_component` (from the cached `st`): log–log error vs h, slope 2, next
    to the C15 curve. *see/read/change:* as C15's.
15. `nb.md` — pointer back to the embedded explainer (embedded once, in C11): "🎮 `stokes_circulation_loop` above: use
    the rectangle mode for Ex. 2.6's four sides and the irrotational-vortex preset for the failing hypothesis."
16. `nb.md` — **What would change if…** "…the loop enclosed a singular point (the vortex core)? The field is not
    differentiable inside, Stokes' hypothesis fails, and the circulation is 2πK however small the loop — the
    'irrotational vortex with circulation' of Ch. 5 and Ch. 6. Everywhere else, (2.34) holds to the last digit."

### A.12 §2.14 Comma Notation — N74 (tagged → C10)
1. `nb.section("2.14", "Comma Notation", intro="**What is this section about?** One more shorthand: a comma in the
   subscript means 'partial derivative with respect to the following index'.")`
2. `nb.note` — **N74 [B]** (2.36): equation `A_{,i} \equiv \partial A/\partial x_i, \qquad \nabla\cdot\mathbf u = u_{i,i},
   \qquad (\nabla\times\mathbf u)_i = \varepsilon_{ijk}u_{k,j}`, ref "2.36" — "the divergence and curl of C10/C11 in one
   line each. A comma index behaves like a tensor index (the chain rule, P49, gives $\partial/\partial x'_j = C_{ij}\,
   \partial/\partial x_i$ — the vector rule (2.8) for ∂), but only in Cartesian coordinates; the book adopts it in
   §5.6. Watch for the comma: $u_{i,j}$ is the velocity gradient of N52, $u_{ij}$ would be something else."
3. `nb.code[explain]` — `print(ch02.expand_indices_str("u_i,i"))` → `Derivative(u_1, x1) + Derivative(u_2, x2) +
   Derivative(u_3, x3)`; `print(ch02.comma_to_partial("u_i,j"))` → `∂u_i/∂x_j`; `print(ch02.expand_indices_str
   ("eps_ijk u_k,j"))` (the three components of (2.25)). *explain:* 1. the comma is parsed as ∂; 2. the text form; 3.
   the curl in comma notation expands to (2.25).

### A.13 End matter — S01, S02, summary
1. `nb.pointer("S01 — Exercises 2.1–2.20 are not reproduced (the text is the book's). The identities we rely on from
   them — (2.7) and orthogonality (2.2, 2.8), invariants (2.9), the tensor examples (2.10), isotropy (2.11), the dot
   and cross products (2.12–2.14), ε–δ (2.5, 2.7), ∇×∇φ = 0 and ∇·∇×u = 0 (2.19, 2.20) — are tested in
   `tests/test_ch02.py`.")`
2. `nb.pointer("S02 — Literature: Sommerfeld's tetrahedron argument for (2.12) is our D05 + D06; Aris and Prager are the
   classical tensor references. Not needed to continue.")`
3. `nb.summary(clicked=["C01 a repeated index is a loop; the free indices are the shape of the answer", "C02 C_ij =
   e_i·e'_j; the columns of C are the new axes; C Cᵀ = I, det C = +1", "C03 x' = Cᵀx defines a vector: any triple that
   fails the test is not one", "C04 τ_ij: first index the face, second the force; positive out of a + face", "C05 f =
   n·τ: the force per area on any plane, from a shrinking tetrahedron", "C06 τ' = CᵀτC — one C per index — is what
   'tensor' means", "C07 closed index chains (trace, I₂, det) are frame-independent", "C08 ε_ijk encodes right-handed
   perpendicularity; ε ε = δδ − δδ closes vector identities", "C09 ∇φ is perpendicular to the level sets and points up
   the steepest slope", "C10 ∇·u is the total stretching rate — what leaks out per unit volume", "C11 ∇×u is twice the
   local spin; a straight shear flow has it", "C12 G = S + A: stretch plus spin; A is the vector ½∇×u; symmetric tensors
   cannot see A", "C13 a symmetric tensor has three perpendicular shear-free axes; its eigenvalues bound every normal
   stress", "C14 Gauss: derivative inside = normal outside; interior faces cancel", "C15 divergence = outflux per unit
   volume; the Cartesian formula falls out of a box", "C16 Stokes: circulation round a loop = curl flux through it;
   curl = circulation per area"], feeds_forward=["Ch. 3: G = S + R from ∂u_i/∂x_j (C12), principal strain rates (C13),
   vorticity = 2 × angular velocity (C11, D15)", "Ch. 4: Cauchy's equation from ∮ n·τ dA and Gauss (C05, C14),
   continuity from the shrinking box (C15), dissipation τ:S (C12), the Newtonian stress from δ (C08)", "Ch. 5: vorticity
   identities with ε–δ (C08), Kelvin's circulation theorem (C16)", "Ch. 13: Coriolis 2Ω × u as ε_ijk (C08), rotating
   frames as C(t) (C02)"], left_out=["Exercises (S01) and the literature (S02)", "curvilinear coordinates: Ch. 3 and
   Appendix B (N72)"])`

### A.14 Placement table (ID → notebook section → host block → call) — all 94 curation rows
| ID | Depth | Section | Host block | Call (Part A row) |
|---|---|---|---|---|
| N01 | B | §2.1 | tagged → C01 (before the block) | A.1 #2 note |
| N02 | B | §2.1 | tagged → C01 | A.1 #2 note |
| N03 | B | §2.1 | C01 | A.1 #6 note |
| N04 | C | §2.1 | C01 | A.1 #6 note |
| N05 | B | §2.1 | C01 | A.1 #18 plotly |
| N06 | B | §2.1 | tagged → C01 | A.1 #2 note |
| C01 | A | §2.1 | CORE | A.1 #3 core |
| N07 | B | §2.1 | C01 | A.1 #10 note |
| N08 | B | §2.1 | C01 | A.1 #13 note |
| N09 | B | §2.2 | C02 | A.2 #5 note + #18 plotly |
| N10 | B | §2.2 | C02 | A.2 #5 note |
| N11 | B | §2.2 | C03 | A.2 #24 note |
| C02 | A | §2.2 | CORE | A.2 #2 core |
| C03 | A | §2.2 | CORE | A.2 #21 core |
| N12 | B | §2.2 (taught in §2.1) | C01 | A.1 #13 note |
| N13 | C | §2.2 | C03 | A.2 #26 note |
| N14 | B | §2.2 | C02 | A.2 #13 note (D03) |
| N15 | B | §2.2 | C02 | A.2 #11 note (D02) |
| N16 | B | §2.3 (taught in §2.2) | C03 | A.2 #26 note; §2.3 pointer A.3 #2 |
| N17 | C | §2.2 | C02 | A.2 #5 note + #18–19 |
| N18 | B | §2.2 | C03 | A.2 #26 note + #31 residual table |
| N19 | B | §2.2 | C03 | A.2 #28 note |
| N20 | C | §2.2 | C03 | A.2 #28 note (E1 polar mode) |
| N21 | B | §2.3 (taught in §2.1) | C01 | A.1 #12 note; §2.3 pointer A.3 #2 |
| N22 | C | §2.3 (taught in §2.1) | C01 | A.1 #12 note |
| N23 | B | §2.3 (taught in §2.1) | C01 | A.1 #12 note + #17 figure |
| N24 | C | §2.4 (taught in §2.1) | C01 | A.1 #2 note |
| C04 | A | §2.4 | CORE | A.4 #2 core |
| N25 | B | §2.4 | C04 | A.4 #8 note + #12 plotly |
| R01 | C | §2.4 | recap (in C04) | A.4 #5 recap |
| N26 | C | §2.4 | C04 | A.4 #8 note |
| C06 | A | §2.4 | CORE | A.4 #31 core |
| N27 | B | §2.4 | C06 | A.4 #35 note + #37 code |
| N28 | B | §2.4 | C06 | A.4 #35 note |
| N29 | B | §2.4 | C06 | A.4 #35 note |
| C07 | A | §2.5 | CORE | A.5 #2 core |
| N30 | C | §2.5 | C07 | A.5 #8 note |
| N31 | B | §2.5 | C07 | A.5 #8 note + #11 asserts |
| N32 | B | §2.5 | C07 | A.5 #8 note |
| N33 | B | §2.5 | C07 | A.5 #8 note + #10 code |
| N34 | B | §2.6 | C05 | A.4 #19 note |
| N35 | B | §2.6 | C05 | A.4 #19 note |
| C05 | A | §2.6 | CORE | A.4 #14 core |
| N36 | B | §2.6 | C05 | A.4 #27 plotly |
| N37 | C | §2.6 | C05 | A.4 #21 note |
| N38 | B | §2.6 | C05 | A.4 #22 example + #24 note |
| N39 | C | §2.6 | C05 | A.4 #24 note + #28 slider |
| N40 | B | §2.7 (taught in §2.1) | C01 | A.1 #11 note |
| N41 | B | §2.7 (taught in §2.1) | C01 | A.1 #11 note |
| N42 | B | §2.7 | C08 | A.6 #10 note |
| C08 | A | §2.7 | CORE | A.6 #2 core |
| N43 | B | §2.7 | C08 | A.6 #10 note |
| N44 | B | §2.7 | C08 | A.6 #10 note (D09) |
| N45 | C | §2.8 (taught in §2.1) | C01 | A.1 #10 note |
| N46 | B | §2.8 | C08 | A.6 #9 note |
| N47 | B | §2.8 | C08 | A.6 #9 note |
| N48 | C | §2.8 | C08 | A.6 #9 note |
| N49 | B | §2.8 | C08 | A.6 #9 note |
| N50 | C | §2.8 | C08 | A.6 #9 note |
| R02 | B | §2.9 | recap (opens §2.9, before C09) | A.7 #2 recap |
| C09 | A | §2.9 | CORE | A.7 #3 core |
| N51 | B | §2.9 | C09 | A.7 #14 figure |
| C10 | A | §2.9 | CORE | A.7 #18 core |
| N52 | B | §2.9 | C10 | A.7 #22 note |
| C11 | A | §2.9 | CORE | A.7 #29 core |
| N53 | B | §2.9 | C11 | A.7 #34 note (D12) |
| N54 | B | §2.9 | C11 | A.7 #34 note |
| N55 | B | §2.9 | C11 | A.7 #34 note + #23/#35 examples |
| N56 | B | §2.10 | C12 | A.8 #7 note |
| C12 | A | §2.10 | CORE | A.8 #2 core |
| N57 | B | §2.10 | C12 | A.8 #9 note |
| N58 | B | §2.10 | C12 | A.8 #9 note (D15) |
| N59 | B | §2.10 | C12 | A.8 #11 note |
| N60 | B | §2.10 | C12 | A.8 #11 note |
| N61 | C | §2.10 | C12 | A.8 #11 note |
| C13 | A | §2.11 | CORE | A.9 #2 core |
| N62 | B | §2.11 | C13 | A.9 #10 note + #11 example |
| N63 | B | §2.11 | C13 | A.9 #10 note + #14 slider |
| C14 | A | §2.12 | CORE | A.10 #2 core |
| N64 | B | §2.12 | C14 | A.10 #9 note |
| N65 | B | §2.12 | C14 | A.10 #9 note + #13 plotly |
| N66 | B | §2.12 | C15 | A.10 #23 note (D21) |
| C15 | A | §2.12 | CORE | A.10 #17 core |
| N67 | B | §2.12 | C15 | A.10 #23 note |
| N68 | B | §2.12 | C15 | A.10 #23 note (D22) + #28 animation |
| N69 | B | §2.13 | C16 | A.11 #7 note |
| N70 | B | §2.13 | C16 | A.11 #7 note + #12 plotly |
| C16 | A | §2.13 | CORE | A.11 #2 core |
| N71 | B | §2.13 | C16 | A.11 #7 note |
| N72 | C | §2.13 (taught in §2.12) | C15 | A.10 #23 note |
| N73 | B | §2.13 | C16 | A.11 #7 note |
| N74 | B | §2.14 | tagged → C10 (own section) | A.12 #2 note |
| S01 | C | Ex. | end | A.13 #1 pointer |
| S02 | C | Lit. | end | A.13 #2 pointer |

Counts: 16 CORE blocks (C01–C16) · 74 NOTE · 2 RECAP · 2 SKIP = 94 rows placed once. Explainer embeds: E1 in C03, E2 in
C05, E3 in C12, E4 in C14, E5 in C11 (C16 refers back) — five `nb.explainer` calls.

---

## Part B — explainer storyboards

Common to all six: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch02">`; tabs Walkthrough ·
Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function that
mirrors a `ch02` callable and is proved by `selftest()` parity rows (`py:` expressions use only `ch02.…`, `np.pi`, lists
and floats — no builtins); Explain is "Explanation & interpretation" in numbered sections with `Viz.work.step / line /
box / table / hint / interpret`; derivation steps are copied from Part F (same `did` titles, same count; phones shorten
*why* to its first sentence); colours as in the header (old frame teal, new orange, normal blue, shear rose, symmetric
teal, antisymmetric orange, flux in blue / out orange, circulation purple, curl + purple / − green, ghosts grey). Angles
are shown in degrees on controls and in radians inside every function (`Viz.fmt(deg)°`). Each fits 360×640 …
1920×1080 and the 1000×700 notebook frame with no scrolling (fit plan per explainer).

### E1 · rotation_of_axes
- **Title:** "Which rotates — the arrow or the ruler?" · **Summary:** "Turn the axes under a fixed vector (or a fixed
  stress state) and watch the direction-cosine matrix C, the new components x' = Cᵀx and the tensor components τ' =
  CᵀτC move together — one rule applied once per index." · **CORE:** C02, C03, C06 (also N15, N16, N18, N19, N20, N27,
  N38) · **Reference:** `amplitude_phase_second_order_II_3.html` (windows linked by one state, numbered live derivation
  in the explanation) + ch01 E5 `buckingham_pi_machine`'s matrix view for C.
- **meta:** `viz:sections 2.2 2.3 2.4` · `viz:equations 2.4 2.5 2.7 2.8 2.12` · `viz:fluidpy ch02.rotation_matrix_2d
  ch02.rotation_matrix_3d ch02.transform_vector ch02.inverse_transform_vector ch02.transform_tensor
  ch02.orthogonality_residual ch02.polar_components ch02.example_2_2 ch02.transforms_as_vector` · `viz:derivations
  D01 D02 D06`.
- **Physics (JS ↔ Python):** `rot2d(theta)` → [[c, −s],[s, c]] ↔ `ch02.rotation_matrix_2d(theta)` · `transformVector(x,
  C)` (= Cᵀx, loops over i) ↔ `ch02.transform_vector(x, C)` · `inverseTransform(xp, C)` ↔ `ch02.inverse_transform_vector`
  · `transformTensor(T, C)` (double loop C_im C_jn T_ij) ↔ `ch02.transform_tensor(T, C)` · `orthoResidual(C)` (max|CᵀC −
  I|) ↔ `ch02.orthogonality_residual(C)` · `detC(C)` · `polar(u1, u2, theta)` ↔ `ch02.polar_components` · `activeRotate(x,
  theta)` (= C x, the Wikipedia action) · `vectorTestResidual(fnName, theta)` for the four test triples ↔
  `ch02.transforms_as_vector(fn, C, points)` (same sample points, seeded) · `shearRotated(a, theta)` (τ' entries of Ex.
  2.2: a sin 2θ, a cos 2θ, −a sin 2θ) ↔ `ch02.example_2_2(a, theta)["tau_rot"]`.
- **Modes** (`mode` chips): **vector** (the arrow x) · **tensor** (a 2-D stress element τ with its three arrows) ·
  **polar** (Ex. 2.1: the point moves on a circle of radius r; the polar frame is the frame turned by the point's
  angle θ). **Toggle `what`**: "axes turn (passive, book)" / "arrow turns (active, Wikipedia R)".
- **Views:**
  1. `plane` "The plane" (`equal: true`) — teal axes 1, 2 fixed; orange axes 1', 2' turned by θ; the arrow x (bold
     black) with teal dashed shadows on the old axes and orange dashed shadows on the new; component labels x₁, x₂,
     x'₁, x'₂ at the shadow feet. Tensor mode: a square element with τ₁₁, τ₂₂ (blue) and τ₁₂ (rose) arrows on its faces
     in the old frame, and a faint rotated square with τ'₁₁, τ'₁₂, τ'₂₂ arrows in the new frame; polar mode: the point at
     angle θ on a dashed circle, e_r and e_θ (orange) at the point, u decomposed twice. Pointer: drag the arrow tip
     (vector mode) → sets `x1, x2`; drag anywhere on the orange axis → sets θ. Active toggle: the teal axes stay, the
     arrow itself turns by θ (ghost of the original arrow in grey).
  2. `matrix` "C = e_i · e'_j" — the 2×2 (3×3 in 3-D option) matrix as cells; hovering/clicking a cell lights row i
     (teal) and column j (orange) and draws the two unit vectors e_i, e'_j on view 1 with the angle arc; below it the
     live assembly line "x'₁ = x₁C₁₁ + x₂C₂₁ = 1×0.866 + 2×0.5 = **1.866**" for the selected j (chips j = 1, 2), and
     "CᵀC − I = 0.0e0 · det C = +1.000" (status colours). Tensor mode: the four-term assembly of τ'_mn for the selected
     (m, n). `hidePortrait: false` (this is the idea).
  3. `curves` "Components vs θ" (`hidePortrait: true`; its key number repeated in view 2's title) — x: θ −180…180°;
     vector mode: x'₁(θ) (orange), x'₂(θ) (orange dashed), the constant |x| (grey) — a moving dot at the current θ;
     tensor mode: τ'₁₁, τ'₂₂ (blue, blue dashed), τ'₁₂ (rose) vs θ with the special angles marked (θ_p where τ'₁₂ = 0,
     the bridge to E2); polar mode: u_r, u_θ vs θ.
- **Controls:** `theta` "Rotation of the axes $\theta$" −180…180 step 1 (default 30) ° help "positive = axes turned
  counterclockwise; the arrow never moves" · `x1` "$x_1$" −3…3 step 0.1 (default 1) · `x2` "$x_2$" −3…3 step 0.1
  (default 2) · `mode` chips vector / tensor / polar · `what` toggle "what rotates: axes ⟷ arrow" · tensor mode: `tpreset`
  chips "pure shear a" (Ex. 2.2, [[0, a],[a, 0]], a = 1 Pa) · "uniaxial" ([[1, 0],[0, 0]]) · "hydrostatic −p" (−I) ·
  "fixed array (not a tensor)" (the same numbers declared in every frame) · *optional* `mirror` toggle "mirror the 2'
  axis (det C = −1)" · *optional* `axis3d` toggle "3-D: rotate about (1, 1, 1)/√3" (matrix view shows 3×3 Rodrigues C;
  plane view shows the projection).
- **Depth features:** Explain + Code · **linked views** (3) · **presets** (θ = 0 "same axes", 45°, 90° "axes swapped",
  "polar frame at θ" (mode polar), "pure shear" (mode tensor), "−30°") · **status** ("✅ orthogonal: CᵀC − I = 2e-16, det C
  = +1" / "⚠️ det C = −1: a mirror, not a rotation" (mirror on) / "❌ not a tensor: (2.12) residual 0.83" (fixed-array
  preset) / "🔄 active: the arrow turned by θ, axes fixed — same C, applied as C not Cᵀ") · **inspector** (click a matrix
  cell: "C₁₂ = e₁·e'₂ = cos(angle between old 1 and new 2) = cos(120°) = **−0.500**" with the two unit vectors drawn) ·
  **modes** (vector / tensor / polar). **Readouts:** "x'₁, x'₂" · "|x| = |x'|" · "CᵀC − I" · tensor: "τ'₁₁ τ'₁₂ τ'₂₂".
- **Explain** ("Explanation & interpretation"):
  0. *What the windows show* — "**The plane**: teal axes are the original frame O12, orange axes the rotated frame
     O1'2'; the black arrow is one vector x (or, in tensor mode, one stress state). Dashed teal/orange lines are its
     shadows on each frame — the components. **C** is the table of cosines between old and new axes; the lit row and
     column show one entry being built. **Components vs θ** traces how the orange numbers change as you turn the
     axes."
  1. *The direction-cosine matrix with your angle* — "C_ij = e_i·e'_j: C₁₁ = cos θ = cos 30° = **0.866**, C₂₁ = e₂·e'₁ =
     cos(90° − θ) = sin θ = **0.500**, C₁₂ = e₁·e'₂ = cos(90° + θ) = −sin θ = **−0.500**, C₂₂ = cos θ = **0.866**. The
     columns (0.866, 0.5) and (−0.5, 0.866) are the new unit vectors written in old components — read them off the
     orange arrows."
  2. *The new components term by term (2.5)* — "x'₁ = x₁C₁₁ + x₂C₂₁ = 1×0.866 + 2×0.500 = **1.866**; x'₂ = x₁C₁₂ +
     x₂C₂₂ = 1×(−0.500) + 2×0.866 = **1.232**. This is x' = Cᵀx: the *first* index of C is summed." (boxed)
  3. *Nothing was lost: orthogonality* — "CᵀC − I: max entry **2e-16** (round-off); det C = cos²θ + sin²θ = **+1**
     (a rotation, not a mirror). Consequence: |x'|² = 1.866² + 1.232² = 3.482 + 1.518 = **5.000** = |x|² = 1² + 2² —
     a rotation keeps lengths (D02)."
  4. *Back again (2.7)* — "x₁ = x'₁C₁₁ + x'₂C₁₂ = 1.866×0.866 + 1.232×(−0.5) = 1.616 − 0.616 = **1.000**: now the
     *second* index of C is summed, x = Cx' (D03)."
  5. *Tensor mode — (2.12) entry by entry* (or hint "switch to tensor mode to see τ' = CᵀτC") — "τ'₁₁ = C_i1 C_j1 τ_ij =
     C₁₁C₂₁τ₁₂ + C₂₁C₁₁τ₂₁ = 2×0.866×0.5×1 = **0.866** Pa; τ'₁₂ = C₁₁C₂₂τ₁₂ + C₂₁C₁₂τ₂₁ = 0.75 − 0.25 = **0.500** Pa;
     τ'₂₂ = 2C₁₂C₂₂τ₁₂ = **−0.866** Pa. Two C's because τ has two indices. Trace τ'₁₁ + τ'₂₂ = 0 = τ₁₁ + τ₂₂: an
     invariant (C07)." For the fixed-array preset: "the array declared 'the same in every frame' differs from CᵀτC by
     **0.83** — it is a matrix, not a tensor (N27)."
  6. *Polar mode — Ex. 2.1* (or hint) — "u_r = u₁cos θ + u₂sin θ = **1.866**, u_θ = −u₁sin θ + u₂cos θ = **1.232**:
     exactly (2.5) with j ∈ {r, θ}; the polar frame is the Cartesian frame turned by the point's angle."
  7. *Passive or active?* — table with the current row highlighted: | what turns | formula | x' now | | axes (book) |
     x' = Cᵀx | (1.866, 1.232) | | arrow (Wikipedia R) | x' = Cx | (−0.134, 2.232) | — "Same matrix, opposite job.
     Turning the axes by +θ is, from the axes' point of view, turning the arrow by −θ."
  8. *Reading the current setting* — θ = 0: "the frames coincide, C = I, x' = x." · 0 < |θ| < 90°: "the orange
     shadows are a reshuffling of the teal ones; the arrow's length 2.236 is untouched." · θ = ±90°: "the axes have
     swapped roles: x'₁ = ±x₂, x'₂ = ∓x₁." · |θ| = 180°: "both components change sign — a half turn, still a rotation
     (det +1)." · tensor mode near θ_p = 45° (pure shear): "τ'₁₂ = a cos 2θ = **0.000**: in these axes the shear has
     vanished and the normal stresses are ±a — the principal axes of §2.11 (E2)." · mirror on: "det C = −1: lengths are
     kept but handedness is not; ε would change sign (N42)."
- **Derivation tab:**
  - **D01** (6 steps) `view: 'plane'`; goal `set {mode: 'vector', theta: 30, x1: 1, x2: 2}`; step 2 (dot with e'_k)
    `highlight ['view:matrix']` and sets the matrix column j = 1 lit; step 5 **live** "C₁₁ = e₁·e'₁ = 0.866, C₂₁ = e₂·e'₁
    = 0.500 → x'₁ = 1×0.866 + 2×0.500 = 1.866"; step 6 `watch: "the orange shadow x'₁ on the 1' axis has length 1.866"`;
    interpret `s => "With θ = … your components are (…, …): the same arrow, new shadows."`.
  - **D02** (8 steps) `view: 'matrix'`; goal `set {mode: 'vector'}`; step 2 highlights column 1 (e₁ expanded in the
    primed basis); step 4 **live** "C₁₁C₁₁ + C₁₂C₁₂ = 0.750 + 0.250 = 1.000 = δ₁₁; C₁₁C₂₁ + C₁₂C₂₂ = 0.433 − 0.433 =
    0 = δ₁₂"; step 7 live "det C = 0.866×0.866 − (−0.5)(0.5) = +1.000"; step 8 `set {mirror: true}` then back, `watch:
    "with the mirror on the status turns amber: det C = −1"`; interpret: "orthogonality is why |x'| = |x| = … in the
    readout".
  - **D06** (8 steps) `view: 'plane'`; goal `set {mode: 'tensor', tpreset: 'shear', theta: 30}`; step 2 highlights the
    traction arrow f on the element's 1' face; step 5 **live** "in the primed frame f'₁ = τ'₁₁n'₁ + τ'₂₁n'₂ with n' =
    (1, 0): f'₁ = τ'₁₁ = 0.866 Pa"; step 8 `set {theta: 45}` `watch: "at 45° the τ'₁₂ curve crosses zero"`; interpret:
    "τ' = CᵀτC with your θ gives τ'₁₁ = a sin 2θ = …".
- **Code:**
  ```python
  th = np.deg2rad(theta)                                  # θ = {{theta}}° → {{th}} rad
  C = ch02.rotation_matrix_2d(th)                         # passive: columns = new axes
  x = np.array([x1, x2])                                  # the fixed arrow ({{x1}}, {{x2}})
  xp = ch02.transform_vector(x, C)                        # (2.5): x' = Cᵀx = ({{xp1}}, {{xp2}})
  back = ch02.inverse_transform_vector(xp, C)             # (2.7): x = Cx' = ({{b1}}, {{b2}})
  res = ch02.orthogonality_residual(C)                    # max|CᵀC − I| = {{res}}
  tau = np.array([[0.0, a], [a, 0.0]])                    # Ex. 2.2's pure shear, a = {{a}} Pa
  tau_p = ch02.transform_tensor(tau, C)                   # (2.12): τ' = CᵀτC → τ'₁₁ = {{t11}}, τ'₁₂ = {{t12}}
  ur, ut = ch02.polar_components(x1, x2, th)              # Ex. 2.1: ({{ur}}, {{ut}})
  ```
- **Walkthrough (8 steps):** 1. "One arrow, two rulers" — "The black arrow is a wind vector. Teal axes read it as (1,
  2). Turn the rulers — the arrow does not care. What do the orange rulers read?" `set {mode: 'vector', theta: 30, x1:
  1, x2: 2}` highlight `view:plane` · 2. "Nine cosines (four here)" — "Every entry of C is the cosine between an old
  axis and a new one. Click a cell: the two unit vectors and their angle appear." `inspect: true`, highlight
  `view:matrix`, `eq: 'Cdef'` · 3. "Shadows add up" — "x'₁ = x₁C₁₁ + x₂C₂₁ = 1×0.866 + 2×0.5 = 1.866. Watch the first
  column light as the sum is built." `eq: 'x5'`, `derive: {id: 'D01', step: 5}` · 4. "Nothing lost" — "CᵀC = I, det C
  = +1: the arrow's length 2.236 is the same in both frames — the status line says so at every θ." `readouts:
  ['norm', 'resid']`, `derive: {id: 'D02', step: 4}` · 5. "Which one really turned?" — "Flip the toggle: now the arrow
  turns and the axes stay. Same C — but applied as C, not Cᵀ. This is the passive/active confusion, settled." `set
  {what: 'arrow'}`, `notes: true` · 6. "Twice for a tensor" — "Tensor mode: a stress state with pure shear a. Its new
  components need C twice, once per index: τ' = CᵀτC." `set {what: 'axes', mode: 'tensor', tpreset: 'shear', theta:
  30}`, `eq: 'tau12'`, `derive: {id: 'D06', step: 8}` · 7. "Find the quiet axes" — "Drag θ until τ'₁₂ = 0 (the rose
  curve crosses zero). At 45° the shear is gone and the normal stresses are ±a. §2.11 calls these principal axes."
  `controls: ['theta']`, `play: false` · 8. "Your turn: polar" — "Polar mode is Ex. 2.1. Predict u_r and u_θ for the
  point at 60°, then drag θ there and read them off." `set {mode: 'polar', theta: 60}`, `controls: ['theta', 'x1',
  'x2']`, `readouts: ['xp']`.
- **Equations:** `Cdef` "Direction cosines" ref '§2.2' `C_{ij} = \mathbf e_i\cdot\mathbf e'_j` live "C = [[0.866, −0.500],
  [0.500, 0.866]]" symbols [['C_{ij}', 'cosine between old axis i and new axis j', '–']] · `x4` "First projection" ref
  'Eq. (2.4)' `\mathbf x\cdot\mathbf e'_1 = x_1\,\mathbf e_1\cdot\mathbf e'_1 + x_2\,\mathbf e_2\cdot\mathbf e'_1 + x_3\,
  \mathbf e_3\cdot\mathbf e'_1 = x'_1` live "= 1×0.866 + 2×0.500 = 1.866" · `x5` "Components transform" ref 'Eq. (2.5)'
  `x'_j = x_i C_{ij} \quad(\mathbf x' = \mathbf C^{\rm T}\mathbf x)` live "(1.866, 1.232)" · `x7` "Inverse" ref 'Eq. (2.7)'
  `x_j = x'_i C_{ji} \quad(\mathbf x = \mathbf C\,\mathbf x')` live "(1.000, 2.000)" · `ortho` "Orthogonality" ref
  'Exercise 2.8' `C_{ij}C_{kj} = \delta_{ik},\quad \det\mathbf C = +1` live "residual 2e-16, det +1" · `tau12` "Tensor
  rule" ref 'Eq. (2.12)' `\tau'_{mn} = C_{im}C_{jn}\tau_{ij} \quad(\boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau
  \mathbf C)` live "τ'₁₁ = 0.866, τ'₁₂ = 0.500, τ'₂₂ = −0.866 Pa" note "pure shear: τ'₁₁ = a sin 2θ, τ'₁₂ = a cos 2θ".
- **Check yourself:** (1) "Set θ = 90°. What are x'₁, x'₂ for x = (1, 2)?" — "(2, −1): the new 1' axis is the old 2
  axis, and the new 2' axis is minus the old 1 axis; C = [[0, −1],[1, 0]]." `set {theta: 90}` · (2) "In tensor mode
  with pure shear, at which θ is τ'₁₁ largest, and how large?" — "45°: τ'₁₁ = a sin 90° = a; the shear τ'₁₂ = 0 there
  — a principal axis." `set {mode: 'tensor', tpreset: 'shear', theta: 45}` · (3) "Turn the mirror on. What is det C, and
  is |x'| still |x|?" — "det C = −1; lengths are still preserved (CᵀC = I holds) but the frame is left-handed; the book's
  'rotation' excludes this." `set {mirror: true}` · (4) "Flip 'what rotates' to the arrow with θ = 30°. Is x' the same
  as before?" — "No: (−0.134, 2.232) instead of (1.866, 1.232) — the active rotation by +θ equals the passive one by
  −θ." `set {what: 'arrow', theta: 30}`.
- **Selftest parity rows:** `{name: 'C(30°)[0][1]', js: rot2d(0.5235987755982988)[0][1], py: 'ch02.rotation_matrix_2d
  (0.5235987755982988)[0][1]', rtol: 1e-12}` · `{name: "x'₁", js: transformVector([1, 2], rot2d(0.5235987755982988))[0],
  py: 'ch02.transform_vector([1.0, 2.0], ch02.rotation_matrix_2d(0.5235987755982988))[0]', rtol: 1e-12}` · `{name:
  "x'₂", …[1], py: '…[1]'}` · `{name: 'round trip', js: inverseTransform(transformVector([1, 2], C), C)[1], py:
  'ch02.inverse_transform_vector(ch02.transform_vector([1.0, 2.0], ch02.rotation_matrix_2d(0.5235987755982988)),
  ch02.rotation_matrix_2d(0.5235987755982988))[1]', rtol: 1e-12}` · `{name: "τ'₁₁ shear", js: transformTensor([[0, 1],
  [1, 0]], rot2d(0.5235987755982988))[0][0], py: 'ch02.transform_tensor([[0.0, 1.0], [1.0, 0.0]], ch02.
  rotation_matrix_2d(0.5235987755982988))[0][0]', rtol: 1e-12}` · `{name: "τ'₁₂ shear", …[0][1]}` · `{name: 'u_r polar',
  js: polar(1, 2, 0.5235987755982988)[0], py: 'ch02.polar_components(1.0, 2.0, 0.5235987755982988)[0]', rtol: 1e-12}` ·
  `{name: 'ortho residual', js: orthoResidual(rot2d(0.7)), py: 'ch02.orthogonality_residual(ch02.rotation_matrix_2d
  (0.7))', atol: 1e-12}` · `{name: 'Ex. 2.2 τ_rot', js: shearRotated(1, 0.5235987755982988)[0], py: 'ch02.example_2_2
  (1.0, 0.5235987755982988)["tau_rot"][0][0]', rtol: 1e-12}` · invariants `{name: '|x| kept', js: norm(xp) − norm(x),
  expect: 0, atol: 1e-12}` · `{name: 'det C', js: detC(rot2d(1.2)), expect: 1, rtol: 1e-12}`.
- **Fit plan:** portrait phone: `plane` (top, square) + `matrix` (bottom) — the curves view hidden, its special angles
  named in the matrix view's title ("τ'₁₂ = 0 at 45°"); Explore shows theta, mode, what (+ x1, x2 optional at density
  3); Derivation tab keeps `plane` (D01, D06) or `matrix` (D02). Landscape / notebook frame: three views in a row
  (1.4 : 1 : 1). Tabs fall back to short labels on 768 px.

### E2 · cauchy_traction_principal_axes
- **Title:** "Cut the point any way you like: what pushes on the cut?" · **Summary:** "Drag the normal n of a cut
  through a stressed point: the traction f = n·τ, its normal and shear parts, the Mohr point and the σ_n(φ), τ_s(φ)
  curves move together; the two shear-free cuts are the eigenvectors, and their normal stresses bound everything." ·
  **CORE:** C05, C13, C04 (also N32, N38, N39, N62) · **Reference:** `forced_damped_vibrations.html` (system + graphs on
  one parameter; the explanation panel with boxed numbers and a regime-dependent reading).
- **meta:** `viz:sections 2.4 2.6 2.11` · `viz:equations 2.15 2.12` · `viz:fluidpy ch02.traction ch02.traction_2d
  ch02.normal_shear_stress ch02.principal_axes ch02.diagonalize ch02.invariants ch02.mohr_circle_2d
  ch02.principal_angle_2d ch02.normal_stress_bounds ch02.example_2_2 ch02.example_2_4` · `viz:derivations D05 D17`.
- **Physics (JS ↔ Python):** `traction(tau, n)` (f_i = Σ_j τ_ji n_j) ↔ `ch02.traction(tau, n)` · `traction2d(tau, phi)`
  → {n, f, sigma_n, tau_s} ↔ `ch02.traction_2d(tau, phi)` · `tractionSecondIndex(tau, n)` (τ·n, for the toggle) ·
  `invariants2d(tau)` → {I1, I2 = det} ↔ `ch02.invariants` (2-D) · `principal2d(tau)` → {lam: [λ_min, λ_max], B, angle}
  (quadratic formula; B det +1) ↔ `ch02.principal_axes(tau)` · `principalAngle(tau)` (½ atan2(2τ₁₂, τ₁₁ − τ₂₂)) ↔
  `ch02.principal_angle_2d` · `mohr(tau)` → {center: I1/2, radius} ↔ `ch02.mohr_circle_2d` · `maxShearAngle(tau)` =
  principal + 45° · `ex22(a, phi)` ↔ `ch02.example_2_2(a, phi)` · `ex24(Gamma)` ↔ `ch02.example_2_4(Gamma)` ·
  `stressVsAngle(tau, phis)` ↔ `ch02.stress_vs_angle`.
- **Views:**
  1. `element` "The stressed point and your cut" (`equal: true`) — a square element with its old-frame arrows (τ₁₁,
     τ₂₂ blue on the ± faces, τ₁₂ = τ₂₁ rose along the faces, C04's sign rule: arrows reversed on the − faces); a cut
     line through the centre at angle φ + 90° with its unit normal n (black), the traction arrow f = n·τ (bold black)
     decomposed into σ_n n (blue) and the shear part (rose) drawn as a right-angled pair; the two eigen-axes b¹, b² as
     faint grey ghosts through the centre with λ labels; when the non-symmetric toggle is on, a second dashed arrow
     τ·n in amber. Pointer: drag anywhere → sets φ to the pointer's angle; click the element face → inspector for that
     face's stresses.
  2. `curves` "σ_n(φ) and τ_s(φ)" — x: φ 0…180°, y: stress (Pa); σ_n (blue) and τ_s (rose, signed as f·s with s the
     counterclockwise tangent) curves with a moving dot at the current φ; dashed horizontal lines at λ_min and λ_max
     (grey, labelled); vertical ticks at the principal angles θ_p, θ_p + 90° ("τ_s = 0") and at θ_p ± 45° ("max
     shear"); title carries the Mohr numbers when view 3 is hidden: "centre 0.00, radius 1.00 Pa".
  3. `mohr` "Mohr's circle" (`hidePortrait: true`, `equal: true`) — axes σ_n (x) and τ_s (y); the circle with centre
     I₁/2 and radius (λ_max − λ_min)/2; the current point (σ_n, τ_s) as a dot, the diameter to its antipode (the
     perpendicular cut), the two principal points on the σ axis; the angle 2φ marked from the σ axis (double angle);
     a strain-rate mode relabels axes to 1/s. Pointer: click a point on the circle → sets φ.
- **Controls:** `t11` "$\tau_{11}$" −3…3 step 0.1 (default 0) Pa · `t22` "$\tau_{22}$" −3…3 (default 0) · `t12`
  "$\tau_{12} = \tau_{21}$" −3…3 (default 1) · `phi` "Cut normal angle $\phi$" 0…180 step 1 (default 30) ° help "angle
  of n from the 1-axis; the cut itself is perpendicular to n" · `mode` chips "stress (Pa)" / "strain rate (1/s), Ex.
  2.4" · `nonsym` toggle "non-symmetric demo τ₂₁ = 0 (contract first vs second index)" *optional* · `tet3d` toggle "3-D
  tetrahedron view" *optional* (replaces view 3 by a `Viz.three` tetrahedron with the same n; `hidePortrait`).
- **Depth features:** Explain + Code · **linked views** (3) · **presets** ("channel shear a > 0" = Ex. 2.2 {t11: 0, t22:
  0, t12: 1, phi: 30} · "channel shear a < 0" {t12: −1} · "uniaxial pull" {t11: 2, t22: 0, t12: 0} · "hydrostatic −p"
  {t11: −2, t22: −2, t12: 0} · "Couette water: τ₁₂ = μU/h = 1.0e-3 × 1 / 1e-2 = 0.1 Pa" {t12: 0.1} · "Ex. 2.4 strain rate
  Γ = 1" {mode: 'strain', t12: 1}) · **status** ("🎯 principal plane: τ_s = 0.000, σ_n = λ_max = 1.000 Pa" when |τ_s| <
  1e-3 · "🔪 maximum-shear plane: τ_s = ±(λ_max − λ_min)/2 = 1.000 Pa" within 0.5° · "🌊 hydrostatic: every plane is
  principal, f = −p n" when radius < 1e-9 · "⚠️ non-symmetric τ: n·τ = (…) ≠ τ·n = (…) — Ch. 4 shows the stress is
  symmetric" when nonsym · else "↗ f = (0.500, 0.866) Pa: σ_n = 0.866, τ_s = 0.500 — 50 % of the maximum shear") ·
  **term bars** ("f decomposed", Pa: f₁ = τ₁₁n₁ + τ₂₁n₂ as two stacked bars (0.000 + 0.500), f₂ = τ₁₂n₁ + τ₂₂n₂ (0.866 +
  0.000); a second group σ_n and τ_s with |f| as the total via √) · **inspector** (click the Mohr point or a curve
  point: "at φ = 30°: n = (0.866, 0.500); f₁ = 0×0.866 + 1×0.500 = 0.500; f₂ = 1×0.866 + 0×0.500 = 0.866; σ_n = f·n =
  0.433 + 0.433 = **0.866**; τ_s = f·s = 0.500×(−0.5) + 0.866×0.866 = **0.500** Pa").
- **Explain:**
  0. *What the windows show* — "**The element** is a point's stress state drawn as a square in its own axes (blue =
     normal stresses, rose = shear, C04's signs). The black line is your cut; **n** its normal; **f** the force per
     area the material on the other side exerts across the cut, split into a **blue** normal push/pull and a **rose**
     slide. Faint grey lines are the two directions with no slide at all. **The curves** show how the blue and rose
     parts change as the cut turns; **Mohr's circle** is the same information as one circle: every cut is a point on
     it."
  1. *Your normal* — "n = (cos φ, sin φ) = (cos 30°, sin 30°) = (**0.866**, **0.500**)."
  2. *Cauchy's formula term by term (2.15)* — "f_i = τ_ji n_j contracts the *first* index (row = face): f₁ = τ₁₁n₁ +
     τ₂₁n₂ = 0×0.866 + 1×0.500 = **0.500** Pa; f₂ = τ₁₂n₁ + τ₂₂n₂ = 1×0.866 + 0×0.500 = **0.866** Pa. |f| = √(0.25 +
     0.75) = **1.000** Pa, direction atan2(0.866, 0.5) = **60°** — Ex. 2.2's numbers (240° for a < 0)." (nonsym on: "τ·n
     would give (…): different — the order is not a formality.")
  3. *Normal and shear parts* — "σ_n = f·n = 0.500×0.866 + 0.866×0.500 = **0.866** Pa (pull, positive); τ_s = f·s with
     s = (−sin φ, cos φ): = **0.500** Pa. For pure shear these are a sin 2φ and a cos 2φ: the double angle is why the
     curves have period 180°."
  4. *The invariants of your τ* — "I₁ = τ₁₁ + τ₂₂ = **0.000** Pa, I₂ = det = τ₁₁τ₂₂ − τ₁₂² = **−1.000** Pa²; the
     eigenvalues solve λ² − I₁λ + I₂ = 0 (D18 in 2-D)."
  5. *Principal values and axes (D17)* — "λ = I₁/2 ± √((τ₁₁ − τ₂₂)²/4 + τ₁₂²) = 0 ± √(0 + 1) = **−1.000, +1.000** Pa;
     principal angle θ_p = ½ atan2(2τ₁₂, τ₁₁ − τ₂₂) = ½ atan2(2, 0) = **45°**; b¹ = (0.707, 0.707), b² = (−0.707,
     0.707) (det +1). Rotating the axes to b¹, b²: τ' = diag(1, −1) — Ex. 2.4 with Γ = 1."
  6. *The bounds (fact 4, sharp form)* — "For every cut, λ_min ≤ σ_n ≤ λ_max: −1.000 ≤ **0.866** ≤ 1.000 ✓. The shear
     never exceeds (λ_max − λ_min)/2 = **1.000** Pa; yours is 0.500 = **50 %** of it, reached at θ_p ± 45° = 0° and
     90°. Mohr: centre I₁/2 = 0.000, radius 1.000; your point sits at angle 2φ = 60° on the circle."
  7. *3-D tetrahedron* (or hint "turn on the 3-D view") — "the same n on a tetrahedron: the three coordinate faces
     carry −τ_j· (reversed arrows, C04), the slanted face carries f; face areas dA_i = n_i dA = (0.866, 0.500, 0) dA."
  8. *Reading the current setting* — principal: "no slide across this cut — the material here is pulled straight
     apart (or squeezed) by λ; a crack would open on this plane." · max shear: "the slide is as large as it can be —
     the plane a ductile material yields on; note it is 45° from the principal axes." · hydrostatic: "f = −p n on every
     cut: the stress state of a fluid at rest (Ch. 1); no direction is special, so there is no principal direction
     to find." · generic: "your cut carries both: … Pa of pull and … Pa of slide (… % of the maximum); turn φ by ±45°
     to trade one for the other." · strain-rate mode: "λ = ±Γ are the stretching and squeezing rates along the 45°
     axes — Ex. 2.4; the fluid element in E3's S panel does exactly this." · nonsym: "with τ₂₁ ≠ τ₁₂ the two
     contractions disagree; Ch. 4 proves the stress is symmetric, so for real stresses they never do."
- **Derivation tab:**
  - **D05** (9 steps) `view: 'element'`; goal `set {tet3d: true, phi: 30}` on big screens (phones keep `element` and a
    2-D triangle sketch drawn by the step's `set {sketch: 'tetra'}`); step 3 highlights the coordinate faces with their
    reversed arrows; step 6 **live** "dA₁ = n₁ dA = 0.866 dA, dA₂ = 0.500 dA"; step 8 `set {shrink: 0.3}` (the element
    drawn smaller) `watch: "the volume terms scale with h³, the face terms with h² — the element can be made as small
    as you like"`; step 9 live "f₁ = τ₁₁n₁ + τ₂₁n₂ = 0.500 Pa, f₂ = 0.866 Pa" `highlight ['readout:f']`; interpret "the
    black arrow on the cut is this formula evaluated with your τ and φ".
  - **D17** (15 steps, ★★★) `view: 'curves'`; goal `set {phi: 30}`; steps 1–5 (real eigenvalues) no set, step 5 live
    "discriminant (τ₁₁ − τ₂₂)²/4 + τ₁₂² = 1.000 ≥ 0: λ real" ; step 8 (orthogonal eigenvectors) `set {phi: 45}` `watch:
    "at φ = θ_p the rose part is zero: f is along n"` and highlight the grey ghost axes; step 11 (τ' diagonal) live "τ'
    = CᵀτC = diag(−1.000, 1.000)"; step 13 live "n·τ·n = Σ λ_k (n·b^k)² = 1×0.966² + (−1)×0.259² = 0.933 − 0.067 =
    0.866 Pa at φ = 30°"; step 14 `watch: "the blue curve never leaves the two dashed lines"`; step 15 live "max shear
    (λ_max − λ_min)/2 = 1.000 Pa at θ_p ± 45°"; interpret "with your τ the bounds are … ≤ σ_n ≤ …; your cut is at … %
    of the shear bound".
- **Code:**
  ```python
  tau = np.array([[t11, t12], [t12, t22]])          # Pa: [[{{t11}}, {{t12}}], [{{t12}}, {{t22}}]]
  phi = np.deg2rad(phi_deg)                         # cut normal at {{phi_deg}}°
  n = np.array([np.cos(phi), np.sin(phi)])          # n = ({{n1}}, {{n2}})
  f = ch02.traction(tau, n)                         # (2.15): f_i = τ_ji n_j = ({{f1}}, {{f2}}) Pa
  sig_n, tau_s, s_dir = ch02.normal_shear_stress(tau, n)   # {{sn}} Pa normal, {{ts}} Pa shear
  lam, B = ch02.principal_axes(tau)                 # λ = ({{l1}}, {{l2}}), columns of B = principal axes
  theta_p = ch02.principal_angle_2d(tau)            # {{thp}}° : the shear-free cut
  c, r = ch02.mohr_circle_2d(tau)                   # Mohr centre {{c}} Pa, radius {{r}} Pa = max shear
  C, tau_p = ch02.diagonalize(tau)                  # τ' = CᵀτC = diag(λ)
  ```
- **Walkthrough (8 steps):** 1. "Nine numbers, one point" — "This point in a channel flow carries a pure shear of 1
  Pa. Cut it along any plane: what pushes across the cut? Drag the normal and watch the black arrow." `set preset
  'channel shear a > 0'`, highlight `view:element` · 2. "Cauchy's formula" — "f₁ = τ₁₁n₁ + τ₂₁n₂ = 0 + 1×0.5 = 0.5 Pa;
  f₂ = 0.866 Pa. Row of τ = face; the first index is contracted with n." `eq: 'cauchy'`, `terms: true`, `derive: {id:
  'D05', step: 9}` · 3. "Split it" — "The blue part pulls straight across the cut (σ_n = 0.866 Pa); the rose part
  slides along it (τ_s = 0.5 Pa). Both change with φ — the curves show how." `inspect: true`, highlight `view:curves` ·
  4. "Find the quiet cuts" — "Drag φ to 45°: the rose part vanishes. Only a normal pull of 1 Pa remains. At 135° a pure
  squeeze of −1 Pa. No other cut is shear-free." `set {phi: 45}`, `controls: ['phi']` · 5. "Those are eigenvectors" —
  "A shear-free cut has f ∥ n: τ·b = λb. For a symmetric τ there are always two (three in 3-D), perpendicular, with
  real λ — proved step by step in the Derivation tab." `eq: 'eigen'`, `derive: {id: 'D17', step: 8}` · 6. "Bounds and
  Mohr" — "Every cut's (σ_n, τ_s) lies on one circle: centre (τ₁₁ + τ₂₂)/2, radius (λ_max − λ_min)/2. So σ_n never
  leaves [λ_min, λ_max] and the shear never exceeds the radius." `eq: 'mohr'`, `derive: {id: 'D17', step: 14}`,
  highlight `view:mohr` · 7. "A fluid at rest" — "Hydrostatic preset: the circle shrinks to a point. f = −p n on every
  cut; no direction is special. That is Ch. 1's pressure." `set preset 'hydrostatic −p'` · 8. "Your turn: Ex. 2.4" —
  "Strain-rate mode with Γ = 1: predict the angle of fastest stretching and its rate, then read λ and θ_p. Then add
  τ₁₁ = 1 and predict how θ_p moves." `set preset 'Ex. 2.4 strain rate Γ = 1'`, `controls: ['t11', 't12', 'phi']`,
  `readouts: ['lam', 'thp']`.
- **Equations:** `cauchy` "Cauchy's traction" ref 'Eq. (2.15)' `f_i = \tau_{ji}n_j,\qquad \mathbf f = \mathbf n\cdot
  \boldsymbol\tau` live "f = (0.500, 0.866) Pa" note "n·τ = τ·n only for symmetric τ (Ch. 4)" symbols [['f_i', 'force
  per unit area on the cut', 'Pa'], ['n_j', 'unit normal of the cut', '–'], ['\\tau_{ji}', 'stress: face j, direction
  i', 'Pa']] · `areas` "Tetrahedron faces" ref '§2.6' `dA_i = n_i\,dA` live "(0.866, 0.500) dA" · `split` "Normal and
  shear parts" ref '§2.6 (Ex. 2.2)' `\sigma_n = \mathbf f\cdot\mathbf n,\quad \tau_s = |\mathbf f - \sigma_n\mathbf n|`
  live "0.866, 0.500 Pa" · `eigen` "Eigen-equation" ref '§2.11' `\det|\tau_{ij} - \lambda\delta_{ij}| = 0,\quad (\tau_{ij}
  - \lambda\delta_{ij})b_j = 0` live "λ = −1.000, +1.000 Pa" · `diag` "Diagonal form" ref '§2.11 fact (3), Eq. (2.12)'
  `\boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau\,\mathbf C = \mathrm{diag}(\lambda^1, \lambda^2)` live "diag(−1.000,
  1.000)" · `mohr` "Mohr's circle (2-D)" ref 'Ch. 4 preview' `\lambda_{1,2} = \tfrac{\tau_{11}+\tau_{22}}{2} \pm \sqrt{
  \big(\tfrac{\tau_{11}-\tau_{22}}{2}\big)^2 + \tau_{12}^2}` live "0 ± 1.000" note "pure shear: σ_n = a sin 2φ, τ_s = a
  cos 2φ".
- **Check yourself:** (1) "Channel shear a = 1 Pa. At which two φ is the shear on the cut largest, and how large?" —
  "0° and 90° (the coordinate faces): τ_s = ±1 Pa = (λ_max − λ_min)/2; there σ_n = 0." `set {phi: 0}` · (2) "Uniaxial
  pull τ₁₁ = 2 Pa. Which cut carries the most shear and how much?" — "45°: τ_s = 1 Pa = half the pull; σ_n there is 1
  Pa too." `set preset 'uniaxial pull', {phi: 45}` · (3) "Set τ₁₁ = 1, τ₂₂ = −1, τ₁₂ = 0. Where are the principal
  axes?" — "θ_p = 0°: the given axes are already principal (τ is diagonal); λ = ±1 and the max-shear planes are at 45°."
  · (4) "Turn on the non-symmetric demo. Which arrow is the traction?" — "The black one (n·τ, first index). The amber
  τ·n differs; for a real stress Ch. 4 proves τ₁₂ = τ₂₁ and they coincide."
- **Selftest parity rows:** `{name: 'f₁ Ex. 2.2', js: traction([[0,1],[1,0]], [0.8660254037844387, 0.5])[0], py:
  'ch02.traction([[0.0, 1.0], [1.0, 0.0]], [0.8660254037844387, 0.5])[0]', rtol: 1e-12}` · `{name: 'σ_n', js:
  traction2d([[0,1],[1,0]], 0.5235987755982988).sigma_n, py: 'ch02.traction_2d([[0.0, 1.0], [1.0, 0.0]], 0.5235987755982988)
  ["sigma_n"]', rtol: 1e-12}` · `{name: 'τ_s', …tau_s…["tau_s"]}` · `{name: 'λ_max', js: principal2d([[0,1],[1,0]]).lam[1],
  py: 'ch02.principal_axes([[0.0, 1.0], [1.0, 0.0]])[0][1]', rtol: 1e-12}` · `{name: 'θ_p', js: principalAngle([[1,0.5],
  [0.5,-1]]), py: 'ch02.principal_angle_2d([[1.0, 0.5], [0.5, -1.0]])', rtol: 1e-12}` · `{name: 'Mohr radius', js:
  mohr([[1,0.5],[0.5,-1]]).radius, py: 'ch02.mohr_circle_2d([[1.0, 0.5], [0.5, -1.0]])[1]', rtol: 1e-12}` · `{name:
  'Ex. 2.2 angle', js: ex22(1, 0.5235987755982988).angle, py: 'ch02.example_2_2(1.0, 0.5235987755982988)["angle_rad"]',
  rtol: 1e-12}` · `{name: 'Ex. 2.4 angle', js: ex24(1).angle, py: 'ch02.example_2_4(1.0)["angle_rad"]', rtol: 1e-12}`
  · `{name: 'non-symmetric differs', js: traction([[0,2],[0,0]], n)[1] − tractionSecondIndex([[0,2],[0,0]], n)[1],
  expect: 1.732 − 0 … (compute: n·τ = (0, 1.732), τ·n = (1, 0)) → js value 1.732, expect: 1.7320508, rtol: 1e-6}` ·
  invariant `{name: 'σ_n within bounds (36 angles)', js: 1 if all inside else 0, expect: 1, atol: 0}`.
- **Fit plan:** portrait phone: `element` (top, square) + `curves` (bottom, with the Mohr numbers in its title); `mohr`
  hidden (Explain §6 carries its numbers); term bars behind a chip; Explore shows t12, phi, mode (+ t11, t22 optional).
  Derivation tab keeps `element` (D05) / `curves` (D17). Landscape: three views in a row (1.2 : 1 : 1); the optional 3-D
  view replaces `mohr` only when its toggle is on and the width ≥ 1000 px.

### E3 · strain_vs_rotation_split
- **Title:** "Is simple shear a rotation?" · **Summary:** "A square of fluid carried by u = G·x deforms and spins at
  once; the same square under S = ½(G + Gᵀ) only stretches, under A = ½(G − Gᵀ) only spins — three panels on one
  clock, with the hidden vector ω of A and the check that a symmetric tensor cannot see A." · **CORE:** C12 (also N52,
  N55, N56–N63) · **Reference:** `angular_frequency_explorer_1.html` (linked views on one clock, modes, presets, "Right
  now" notes, end-of-run summary).
- **meta:** `viz:sections 2.10 2.11` · `viz:equations 2.26 2.27 2.28 2.29` · `viz:fluidpy ch02.symmetric_part
  ch02.antisymmetric_part ch02.antisymmetric_from_vector ch02.vector_from_antisymmetric ch02.symmetric_double_contraction
  ch02.principal_axes ch02.linear_flow_map ch02.deform_square ch02.material_line_angle ch02.example_2_4` ·
  `viz:derivations D14 D15`.
- **Physics (JS ↔ Python):** `symPart(G)`, `antiPart(G)` ↔ `ch02.symmetric_part`, `ch02.antisymmetric_part` · `omegaOf(A)`
  (3-D embedding of the 2-D A: ω₃ = −A₁₂ = A₂₁) ↔ `ch02.vector_from_antisymmetric` · `antiFromOmega(w)` ↔
  `ch02.antisymmetric_from_vector` · `expm2(G, t)` (closed form for 2×2: eigen-decomposition or the series to 20 terms)
  ↔ `ch02.linear_flow_map(G, t)` · `deformSquare(G, t, n)` ↔ `ch02.deform_square` · `lineAngle(G, t, th0)` ↔
  `ch02.material_line_angle` · `principal2d(S)` ↔ `ch02.principal_axes` · `doubleDotFro(A, B)` and `symContraction(tau,
  B)` → {P, PS, PA} ↔ `ch02.symmetric_double_contraction` · `ex24(Gamma)` ↔ `ch02.example_2_4`.
- **Views:**
  1. `motion` "The same square under G, S and A" — three sub-panels side by side (one panel with chips G / S / A / all
     on phones): a 10×10 square of tracer dots (teal for S, orange for A, black for G) advected exactly by x(t) =
     e^{Mt}x₀ for M = G, S, A from the same start; faint grey streamlines of each linear field behind; on the S panel
     the principal axes of S as blue (stretch, λ > 0) and rose (squeeze, λ < 0) double arrows with rates; on the A
     panel a small paddle wheel turning at ω₃ (= −A₁₂) with the ω vector drawn out of the plane as ⊙/⊗; on the G panel
     the outline of the S-only and A-only shapes as ghosts. Pointer: click a tracer → inspector (its start, its
     position now under each of G, S, A).
  2. `matrix` "G = S + A, entry by entry" — three 2×2 boxes with signed bars for each entry, the identities S_ij =
     ½(G_ij + G_ji), A_ij = ½(G_ij − G_ji) written under the selected entry (chips); to the right the term bars "S:S",
     "A:A", "S:A" (Frobenius) — the (2.29) check S:A = 0 in every state; and "ω₃ = −A₁₂ = …".
  3. `angle` "Angle of a material line vs t" (`hidePortrait: true`; its end value in view 1's title) — x: t 0…t_end,
     y: angle (°) of the line that started at θ₀ = 0°: under S (teal: tends to the stretching axis), under A (orange:
     uniform rotation at ω₃), under G (black: the combination), with the analytic e^{Gt} ghost dashed; the current
     time as a vertical line.
- **Controls:** `g11` "$G_{11}$" −2…2 step 0.1 (default 0) 1/s · `g12` "$G_{12} = \partial u_1/\partial x_2$" −3…3 step
  0.1 (default 2) · `g21` "$G_{21}$" (default 0) · `g22` "$G_{22}$" (default 0) · `show` chips "G · S · A · all" (default
  all) · `t` transport 0…1.5 s (rate 0.5 s per real second, end 'hold' with an end card, loop optional) · *optional*
  `theta0` "material line start angle" 0…180 (default 0) · *optional* `preset3d` toggle "3-D G (matrix view only; ω
  drawn as a vector)".
- **Depth features:** Explain + Code · **linked views** (3) · **transport** (`param: 't'`, 0 → 1.5 s, rate 0.5,
  end 'hold', end card: "after 1.5 s: stretched ×e^{λt} = 4.48 along +45° and ×0.22 along −45° (S); turned −1.5 rad =
  −85.9° (A); the G shape is the shear x₁ → x₁ + 3x₂" — numbers for the default) · **presets** ("simple shear u₁ = Γx₂
  (Ex. 2.4)" {g12: 2, others 0} · "solid-body rotation b × x (Ex. 2.3)" {g12: −1, g21: 1} · "pure strain" {g11: 1, g22:
  −1} · "uniaxial extension" {g11: 1} · "irrotational shear (rotation removed)" {g12: 1, g21: 1} · "expansion a x"
  {g11: 1, g22: 1}) · **modes** (`show`: G / S / A / all) · **term bars** (S:S, A:A, S:A → total G:G; the S:A bar is
  pinned at 0 with a "(2.29)" tag) · **status** ("🔄 pure rotation: S = 0 — shape kept, spin ω₃ = 1.00 rad/s" when |S| <
  1e-9 · "↔ pure strain: A = 0 — no spin; stretch rates ±1.00 1/s" when |A| < 1e-9 · "🌀 simple shear: equal parts
  stretch (±1.00 1/s along ±45°) and spin (ω₃ = −1.00 rad/s) — Ex. 2.4" when G₁₁ = G₂₂ = 0 and exactly one of G₁₂,
  G₂₁ nonzero · "⤢ expansion: ∇·u = G_ii = 2.00 1/s" when S ∝ I · else "mixed: stretch rates … along …°, spin …
  rad/s") · **notes** (the "Right now" table: | motion | after t = 0.60 s | | S: stretch along b¹ | ×e^{λ₁t} = 1.82 | |
  S: squeeze along b² | ×0.55 | | A: turned | ω₃t = −0.60 rad = −34.4° | | G: shear displacement | Γt = 1.20 × x₂ |) ·
  **inspector** (click a tracer).
- **Explain:**
  0. *What the windows show* — "Three copies of one fluid square, released at t = 0 into three linear velocity fields
     that share one clock: **black** under the full gradient G, **teal** under its symmetric part S (the strain rate),
     **orange** under its antisymmetric part A (the rotation rate). Blue/rose double arrows are S's principal axes; the
     wheel spins with A. The matrix view shows the split entry by entry; the angle view tracks one material line."
  1. *Your velocity gradient* — "G = [[0, 2],[0, 0]] 1/s: u₁ = 2x₂, u₂ = 0 — a simple shear (Couette-like: faster
     above, slower below). G_ij = ∂u_i/∂x_j, first index the velocity component, second the coordinate."
  2. *The split (D14)* — "S_ij = ½(G_ij + G_ji): S₁₂ = ½(2 + 0) = **1.00**, S₁₁ = S₂₂ = 0; A_ij = ½(G_ij − G_ji): A₁₂ =
     ½(2 − 0) = **1.00**, A₂₁ = **−1.00**. Check: S + A = [[0, 1 + 1],[1 − 1, 0]] = G ✓. Unique: no other symmetric +
     antisymmetric pair adds to G."
  3. *The hidden vector (D15)* — "ω_k = −½ ε_ijk A_ij: ω₃ = −½(ε₁₂₃A₁₂ + ε₂₁₃A₂₁) = −½(1×1.00 + (−1)(−1.00)) = **−1.00**
     rad/s: a clockwise spin. Indeed A·x = ω × x: for x = (1, 0): A·x = (0, −1), ω × x = (0, 0, −1) × (1, 0, 0) = (0,
     −1, 0) ✓. ω is half the curl: (∇×u)₃ = ∂u₂/∂x₁ − ∂u₁/∂x₂ = 0 − 2 = −2 = 2ω₃ ✓ (C11)."
  4. *Principal axes and rates of S (Ex. 2.4, D17)* — "λ = ±√(S₁₁S₂₂… ) — for S = [[0, 1],[1, 0]]: λ = **+1.00, −1.00**
     1/s along θ_p = ½ atan2(2S₁₂, S₁₁ − S₂₂) = **45°** and 135°: stretching at 1.00 1/s along +45°, squeezing at 1.00
     1/s along −45°. (This is Ex. 2.4 with Γ = S₁₂ = 1; note du₁/dx₂ = 2 = 2S₁₂.)"
  5. *Where the square is at t* — live: "t = **0.60** s: under S a length along b¹ has grown by e^{λ₁t} = e^{0.60} =
     **1.82** and along b² shrunk to e^{−0.60} = **0.55** (area kept: 1.82 × 0.55 = 1.00, because S_ii = ∇·u = 0); under
     A it has turned by ω₃t = **−0.60 rad = −34.4°**; under G the top edge has slid by Γt·x₂ = 1.20 x₂ — a lean, not a
     turn."
  6. *A symmetric tensor cannot see A (2.28)–(2.29)* — "Take the probe tensor τ = S (symmetric). τ:S = S_ij S_ij =
     **2.00**, τ:A = S_ij A_ij = 1×1 + 1×(−1) = **0.00**. Swapping the indices of A flips its sign but not τ's, so the
     sum must equal its own negative: zero. This is why viscous dissipation (Ch. 4) sees only S."
  7. *Angle view* (or hint "turn the phone sideways / open the angle view") — "a material line starting along the
     1-axis: under A it turns uniformly at −57.3°/s; under S it swings toward the stretching axis (+45°) and stays;
     under G it does both — the black curve is the shear's actual line angle, atan(−Γt·…)."
  8. *Reading the current setting* — simple shear: "**half stretch, half spin.** The lean you see under G is a stretch
     along +45° plus a squeeze along −45° (teal) *and* a clockwise turn (orange). Neither alone looks like shear; their
     sum is exactly shear. A paddle wheel in this flow turns at Γ/2 even though every streamline is straight." · pure
     rotation: "S = 0: no deformation at all, the square is carried round rigidly at ω₃; ∇×u = 2ω₃ = … — Ex. 2.3." ·
     pure strain: "A = 0: no spin; the square stretches and squeezes along fixed axes; a wheel would not turn." ·
     irrotational shear: "G symmetric, so G = S: the same S as simple shear but no rotation — the stretch without the
     lean." · expansion: "S = aI: every direction stretches equally (∇·u = 2a in 2-D), the square grows and keeps its
     shape; A = 0."
- **Derivation tab:**
  - **D14** (6 steps) `view: 'matrix'`; goal `set {g12: 2, g11: 0, g21: 0, g22: 0, show: 'all'}`; step 1 highlights the
    entry chips; step 2 **live** "S₁₂ = ½(2.00 + 0.00) = 1.00 = S₂₁"; step 3 live "A₁₂ = ½(2.00 − 0.00) = 1.00, A₂₁ =
    −1.00, A₁₁ = A₂₂ = 0"; step 5 (uniqueness) `set {show: 'all'}` `watch: "both parts drawn: the teal square and the
    orange square together reproduce the black one"`; step 6 live "S transforms by (2.12): with C the 45° rotation, S'
    = diag(1.00, −1.00)"; interpret "your G splits into S = … and A = …".
  - **D15** (7 steps) `view: 'motion'`; goal `set {show: 'A'}`; step 2 **live** "R₁₂ = −ε₁₂₃ω₃ = −(−1.00) = 1.00 = A₁₂
    ✓"; step 4 live "ε_ijl A_ij = ε₁₂₃A₁₂ + ε₂₁₃A₂₁ = 1.00 + 1.00 = 2.00 = −2ω₃ → ω₃ = −1.00"; step 6–7 `watch: "the
    wheel turns clockwise: A·x = ω × x with ω pointing into the screen (⊗)"`; interpret "ω = ½∇×u = … rad/s: Ch. 3's
    vorticity is 2ω".
- **Code:**
  ```python
  G = np.array([[g11, g12], [g21, g22]])              # 1/s: [[{{g11}}, {{g12}}], [{{g21}}, {{g22}}]]
  S = ch02.symmetric_part(G)                          # ½(G + Gᵀ): S12 = {{s12}}
  A = ch02.antisymmetric_part(G)                      # ½(G − Gᵀ): A12 = {{a12}}
  assert np.allclose(S + A, G)                        # the split is exact
  A3 = np.pad(A, ((0, 1), (0, 1)))                    # embed in 3-D to read the vector
  omega = ch02.vector_from_antisymmetric(A3)          # (2.27): ω = (0, 0, {{w3}}) rad/s
  lam, B = ch02.principal_axes(S)                     # stretch rates {{l1}}, {{l2}} 1/s along B's columns
  x_t = ch02.deform_square(G, t)                      # tracers at t = {{t}} s: x(t) = expm(G t) x0
  P, P_S, P_A = ch02.symmetric_double_contraction(S, G)   # τ:B = {{P}}, τ:S = {{PS}}, τ:A = {{PA}} (= 0)
  ```
- **Walkthrough (7 steps):** 1. "A square in a shear flow" — "Faster fluid above, slower below: u₁ = 2x₂. Play. The
  black square leans. Is it stretching? Turning? Both?" `set preset 'simple shear'`, `set {show: 'G'}`, `play: true` ·
  2. "Split the gradient" — "Write G as ½(G + Gᵀ) + ½(G − Gᵀ). Both halves are simple: S = [[0, 1],[1, 0]], A = [[0,
  1],[−1, 0]]." `eq: 'split'`, `derive: {id: 'D14', step: 3}`, highlight `view:matrix` · 3. "S alone: stretch" — "Release
  the same square into u = S·x: it stretches along +45° and squeezes along −45° at 1 1/s each — no turning. Those axes
  are Ex. 2.4's principal axes." `set {show: 'S', t: 0}`, `play: true` · 4. "A alone: spin" — "Now u = A·x: the square
  turns rigidly, clockwise, at 1 rad/s. Three numbers of A are a vector: ω = (0, 0, −1). And A·x = ω × x." `set {show:
  'A', t: 0}`, `play: true`, `eq: 'omega'`, `derive: {id: 'D15', step: 6}` · 5. "Add them back" — "All three on one
  clock: lean = stretch + spin. Simple shear is half pure strain, half solid-body rotation — Ex. 2.4 and Ex. 2.3
  glued together." `set {show: 'all', t: 0}`, `play: true`, `notes: true` · 6. "A symmetric tensor is blind to A" —
  "The bars: S:S = 2, A:A = 2, but S:A = 0 — always. Drag any slider and watch the S:A bar stay at zero. Ch. 4's
  dissipation lives here." `terms: true`, `eq: 'blind'`, `controls: ['g12', 'g21']` · 7. "Your turn" — "Predict: set G₂₁
  = G₁₂ = 1. Does the wheel turn? Then G₂₁ = −G₁₂: does the square change shape?" `set {g12: 1, g21: 1, t: 0}`,
  `controls: ['g11', 'g12', 'g21', 'g22']`, `readouts: ['omega', 'lam']`.
- **Equations:** `split` "Symmetric + antisymmetric" ref '§2.10' `B_{ij} = \tfrac12(B_{ij} + B_{ji}) + \tfrac12(B_{ij} -
  B_{ji}) = S_{ij} + A_{ij}` live "S₁₂ = 1.00, A₁₂ = 1.00 1/s" · `R26` "The tensor of a vector" ref 'Eq. (2.26)' `\mathbf
  R = \begin{bmatrix}0 & -\omega_3 & \omega_2\\ \omega_3 & 0 & -\omega_1\\ -\omega_2 & \omega_1 & 0\end{bmatrix}` live
  "ω₃ = −1.00 → R₁₂ = 1.00" · `omega` "Two-way map" ref 'Eq. (2.27)' `R_{ij} = -\varepsilon_{ijk}\omega_k,\qquad \omega_k =
  -\tfrac12\varepsilon_{ijk}R_{ij}` live "ω = (0, 0, −1.00) rad/s; R·x = ω × x" note "the book's lower limit 'i−1' reads
  i = 1" · `strain` "Strain rate (Ex. 2.4)" ref '§2.11' `S_{ij} = \tfrac12\Big(\frac{\partial u_i}{\partial x_j} +
  \frac{\partial u_j}{\partial x_i}\Big)` live "S₁₂ = ½ du₁/dx₂ = 1.00 1/s" note "Γ ≡ S₁₂ here; the book's '2S₁₂ = Γ'
  line carries a stray 2" · `blind` "Symmetric × antisymmetric = 0" ref 'Eqs. (2.28), (2.29)' `\tau_{ij}B_{ij} = \tau_{ij}
  S_{ij},\qquad \tau_{ij}A_{ij} = 0` live "S:G = 2.00 = S:S; S:A = 0.00" · `expm` "Trajectories" ref 'our addition'
  `\mathbf x(t) = e^{\mathbf G t}\,\mathbf x_0` live "t = 0.60 s".
- **Check yourself:** (1) "Simple shear with G₁₂ = 2. At what rate does the paddle wheel turn, and which way?" — "ω₃ =
  −G₁₂/2 = −1 rad/s, clockwise: the fluid above moves right faster, dragging the top of the wheel right." `set preset
  'simple shear'` · (2) "Set G₁₂ = G₂₁ = 1. Does the S panel differ from the simple-shear case?" — "No: S is the same
  [[0, 1],[1, 0]]; only A changed (to zero). The stretch is identical, the spin is gone." `set {g12: 1, g21: 1}` · (3)
  "Solid-body preset. What is ∇·u, and what does the S panel show?" — "∇·u = G_ii = 0 and S = 0: nothing happens in the
  S panel; the square is carried round rigidly at ω₃ = 1 rad/s." `set preset 'solid-body rotation'` · (4) "Why is S:A
  always zero, even for random sliders?" — "Renaming i ↔ j leaves S_ij A_ij unchanged but A_ji = −A_ij, so the sum
  equals its own negative."
- **Selftest parity rows:** `{name: 'S12', js: symPart([[0,2],[0,0]])[0][1], py: 'ch02.symmetric_part([[0.0, 2.0], [0.0,
  0.0]])[0][1]', rtol: 1e-12}` · `{name: 'A21', js: antiPart([[0,2],[0,0]])[1][0], py: 'ch02.antisymmetric_part([[0.0,
  2.0], [0.0, 0.0]])[1][0]', rtol: 1e-12}` · `{name: 'ω₃', js: omegaOf(antiPart([[0,2],[0,0]])), py: 'ch02.
  vector_from_antisymmetric([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])[2]', rtol: 1e-12}` · `{name: 'R12 from
  ω', js: antiFromOmega([0,0,-1])[0][1], py: 'ch02.antisymmetric_from_vector([0.0, 0.0, -1.0])[0][1]', rtol: 1e-12}` ·
  `{name: 'expm shear', js: expm2([[0,2],[0,0]], 0.6)[0][1], py: 'ch02.linear_flow_map([[0.0, 2.0], [0.0, 0.0]], 0.6)
  [0][1]', rtol: 1e-10}` · `{name: 'expm rotation', js: expm2([[0,-1],[1,0]], 0.6)[0][0], py: 'ch02.linear_flow_map
  ([[0.0, -1.0], [1.0, 0.0]], 0.6)[0][0]', rtol: 1e-10}` · `{name: 'line angle', js: lineAngle([[0,2],[0,0]], 0.6, 0.3),
  py: 'ch02.material_line_angle([[0.0, 2.0], [0.0, 0.0]], 0.6, 0.3)', rtol: 1e-10}` · `{name: 'P_A = 0', js: symContraction
  (S, G).PA, py: 'ch02.symmetric_double_contraction([[0.0, 1.0], [1.0, 0.0]], [[0.0, 2.0], [0.0, 0.0]])[2]', atol:
  1e-12}` · `{name: 'Ex. 2.4 λ', js: ex24(1).lam[0], py: 'ch02.example_2_4(1.0)["lam"][0]', rtol: 1e-12}` · invariant
  `{name: 'area kept when S_ii = 0', js: det(expm2(S, 0.6)), expect: 1, rtol: 1e-9}`.
- **Fit plan:** portrait phone: `motion` as one panel with the G/S/A chips (top, square) + `matrix` (bottom); `angle`
  hidden (its end value in the motion title: "line angle now: S 38.2°, A −34.4°, G −29.5°"); the transport's step/speed
  hidden; end card via `Viz.card`. Landscape / notebook: three views (1.6 : 1 : 1) with `motion` split into three
  panels when its width ≥ 540 px.

### E4 · gauss_flux_box
- **Title:** "What leaks out of a box?" · **Summary:** "Move and resize a box over a 2-D field (unit depth): the four
  face fluxes n·Q add up to exactly the integrated divergence inside; tile the box and watch interior faces cancel;
  shrink it and watch (1/V)∮ n·Q dA settle on ∇·Q — Gauss' theorem and the outflux definition of divergence as one
  experiment." · **CORE:** C14, C15 (also C10, N54, N55, N64, N65, N66, N68) · **Reference:** `fid_formula_lab.html`
  (term-by-term bars summing to the total, clickable terms) with the `forced_damped_vibrations.html` explanation panel.
- **meta:** `viz:sections 2.9 2.12` · `viz:equations 2.23 2.30 2.31 2.32` · `viz:fluidpy ch02.divergence_theorem_rect2d
  ch02.flux_through_faces ch02.divergence_theorem_tiled ch02.integral_divergence_2d ch02.divergence ch02.radial_field
  ch02.solid_body_rotation_field ch02.point_source_field ch02.smooth_test_field` · `viz:derivations D25 D22 D21`.
- **Physics (JS ↔ Python):** field presets as functions `Q(x, y)` → [Q₁, Q₂] and their exact divergence `divQ(x, y)`
  (`radial`: (ax, ay), div 2a · `solid`: (−by, bx), div 0 · `xsq`: (x², 0), div 2x · `source`: m(x, y)/(2πr²), div 0
  away from 0 (flux 2π m… per unit depth m) · `smooth`: (sin x cos y, cos x sin y), div 2 cos x cos y · `corner_sink`:
  (−(x−1), −(y−1))·0.5, div −1) ↔ `ch02.radial_field(a)`, `ch02.solid_body_rotation_field`, `ch02.point_source_field(m)`,
  `ch02.smooth_test_field()` · `faceFluxes(Q, box, n)` → {left, right, bottom, top} (midpoint rule with n points per
  face; unit depth) ↔ `ch02.flux_through_faces(Q_fn, bounds, n)` · `boxTheorem(Q, box, n)` → {lhs (midpoint sum of
  divQ), rhs (Σ faces)} ↔ `ch02.divergence_theorem_rect2d` · `tiled(Q, box, tiles, n)` → {sumTiles, outer, interior}
  ↔ `ch02.divergence_theorem_tiled` · `integralDiv(Q, x0, h, n)` (= rhs/h²) ↔ `ch02.integral_divergence_2d(Q_fn, x0, h,
  n_side)` · `taylorFace(Q, x0, h, face)` (Q(x₀) ± (h/2) ∂Q/∂x by central differences — Ex. 2.5's face values).
- **Views:**
  1. `field` "The field and your box" (`equal: true`) — heatmap of ∇·Q (blue negative … white 0 … orange positive)
     over [−2, 2]², faint streamlines, the draggable box (centre `cx, cy`, side `h`) with arrows n·Q on each face
     (outward orange when positive, inward blue when negative, length ∝ |n·Q| at the face midpoint) and the face
     flux value labelled; tiles on: the box split into k×k sub-boxes with the interior face arrows drawn in opposite
     pairs (fading when "cancel" is played); Taylor toggle: on each face the value "Q₁ + (h/2)∂Q₁/∂x₁ = …" printed.
     Pointer: drag inside the box → move; drag a corner → resize (`h`); click a face → inspector.
  2. `bars` "Face fluxes add up" — term bars: left, bottom, right, top (signed, coloured by sign) stacking to the
     total ∮ n·Q dA (bold), next to a second bar ∬ ∇·Q dA (purple) of equal height; the gap printed as "|Δ| = 1.2e-15";
     tiles on: a third group "interior faces: +… − … = 0".
  3. `limit` "(1/V)∮ n·Q dA → ∇·Q(x₀)" (`hidePortrait: true`; the current ratio repeated in view 2's title) — x: box
     side h (log, 0.01…4), y: the ratio for the current centre (purple curve, computed for 30 sizes) with the dashed
     exact ∇·Q(x₀) (grey) and the current h as a dot; a second axis or a small inset: |ratio − exact| vs h on log–log
     with the slope-2 ghost; for linear/quadratic fields the curve is flat (exact) and the title says so.
- **Controls:** `preset` chips "a x (Ex. 2.3)" · "b × x (solenoidal)" · "(x², 0)" · "point source" · "smooth sin/cos" ·
  "corner sink" · `cx`, `cy` "Box centre" −2…2 step 0.05 (default 0.5, 0.5) m · `h` "Box side $h$" log 0.01…3 (default
  1) m · `tiles` chips "1 · 2×2 · 4×4" · `taylor` toggle "show Ex. 2.5's Taylor face values" *optional* · `nface`
  "points per face" 1…32 (default 8) *optional*.
- **Depth features:** Explain + Code · **linked views** (3) · **term bars** (four faces → total vs ∬∇·Q; interior
  faces group) · **presets** (six fields + "shrink to the point" {h: 0.05} + "tile it" {tiles: 4}) · **status** ("🟠
  net outflow 2.000 m³/s per m depth: a source inside (∇·Q = 2.00)" · "🔵 net inflow …: a sink inside" · "⚖️ solenoidal:
  ∮ n·Q dA = 0.0e0 — every box, anywhere" (|total| < 1e-12 and preset solenoidal) · "📐 h = 0.05: (1/V)∮ = 2.000 is
  within 0.0 % of ∇·Q(x₀) = 2.000" when h small · "⚠️ singular point inside the box: the divergence is a delta — the
  flux is m = … regardless of h (Gauss needs a smooth Q)" for the source preset with the origin inside) · **inspector**
  (click a face: "right face x = 1.000, n = (1, 0): n·Q at midpoints = Q₁ = 1.000 (constant along the face for a x);
  flux = 1.000 × h = 1.000 m³/s per m") · **transport-like scrub** on h (the `limit` view's dot moves; a "shrink" button
  animates h from 2 to 0.02 in 3 s).
- **Explain:**
  0. *What the windows show* — "**The field** is coloured by its divergence ∇·Q (orange where the field spreads out,
     blue where it converges, white where it is solenoidal); grey lines follow Q. The square is your box; each face
     carries an arrow for n·Q, the outward flux density. **The bars** add the four faces and compare with the
     divergence integrated inside; **the limit view** divides the flux by the box's area and shrinks the box."
  1. *The field at your box's centre* — "Q = a x with a = 1: Q(0.5, 0.5) = (**0.500**, **0.500**); ∇·Q = ∂Q₁/∂x₁ +
     ∂Q₂/∂x₂ = 1 + 1 = **2.000** (2-D; the 3-D field a x has 3a)."
  2. *Each face's flux (n·Q × length)* — "right (x = 1.000, n = +e₁): n·Q = Q₁ = 1.000 along the whole face, flux =
     1.000 × 1.000 = **+1.000**; left (x = 0, n = −e₁): −Q₁ = 0 → **0.000**; top (y = 1, n = +e₂): **+1.000**; bottom:
     **0.000**. Sum ∮ n·Q dA = **2.000** m³/s per metre of depth." (boxed)
  3. *The inside (2.30)* — "∬ ∇·Q dA = 2.000 × area 1.000 = **2.000**. Difference from the face sum: **0.0e0** — Gauss'
     theorem to round-off (for this linear field the midpoint rule is exact; for the smooth preset the gap shrinks like
     1/n²_face)."
  4. *Tiling (D25 step 8)* (or hint "set tiles to 4×4") — "16 tiles, each obeying (2.30); their 16 boundary sums add
     to **2.000**; the 24 interior faces appear twice with opposite normals and the same Q, contributing **0.0e0** in
     total; only the outer 16 face segments survive — the outer flux 2.000."
  5. *Ex. 2.5's face values (D22)* (or hint "toggle the Taylor face values") — "Q₁ on the right face = Q₁(x₀) + (h/2)
     ∂Q₁/∂x₁ = 0.500 + 0.500×1 = **1.000**; on the left = 0.500 − 0.500 = **0.000**. Right minus left: the 0.500's
     cancel, the slopes add to h ∂Q₁/∂x₁ = 1.000; times the face length h: **1.000** h² … divided by the area h²:
     ∂Q₁/∂x₁ = **1.000**. Same for y: total **2.000** = ∇·Q."
  6. *The limit (2.32), (D21)* — "(1/V)∮ n·Q dA at h = 1.000: **2.000**; exact ∇·Q(x₀) = 2.000; error **0.0e0**. For
     the (x², 0) preset the ratio is 2x₀ exactly at every h; for the smooth preset the error falls as h² (slope 2 on
     the log–log inset): the definition needs the limit."
  7. *At the current box* — live: h, centre, the four fluxes, total, ∬∇·Q, ratio.
  8. *Reading the current setting* — source inside: "more leaves than enters: **2.000** m³/s per metre of depth flows
     out of this box. Per unit area that is the divergence 2.000 — the fluid inside would have to be created (or, in
     Ch. 4, its density would have to fall) at this rate." · sink: "…flows *in*: convergence; the density would rise."
     · solenoidal: "in = out for *every* box: what enters through one face leaves through another — the incompressible
     flows of Ch. 4–9 look like this everywhere." · small box: "the box is small enough that (1/V)∮ has converged to the
     point value: this is the divergence *defined* without coordinates (2.32)." · singular source inside: "the field is
     not differentiable at the origin; the flux m = … is the same for every box enclosing it — Gauss' theorem needs a
     smooth Q inside V, and (2.32) gives a delta function, not a number."
- **Derivation tab:**
  - **D25** (9 steps) `view: 'field'`; goal `set {preset: 'radial', cx: 0.5, cy: 0.5, h: 1, tiles: 1}`; step 4 (n₁ = ±1 on
    two faces, 0 on four) highlights the left and right faces; step 5 **live** "∬ Q₁(right) dA − ∬ Q₁(left) dA = 1.000
    − 0.000 = 1.000 = ∬ ∂Q₁/∂x₁ dA"; step 7 `set {tiles: 4}`; step 8 `watch: "the interior arrows come in opposite
    pairs — their bar reads 0"` `highlight ['view:bars']`; step 9 `set {tiles: 1}` `watch: "only the outer faces are
    left: 2.000"`; interpret "for your box: inside 2.000 = outside 2.000".
  - **D22** (9 steps) `view: 'field'`; goal `set {preset: 'radial', taylor: true, h: 1}`; step 2 **live** "Q₁(right) =
    0.500 + 0.500 = 1.000, Q₁(left) = 0.500 − 0.500 = 0.000"; step 4 live "right − left: (1.000 − 0.000) × h = 1.000 — the
    0.500's cancelled"; step 6 highlights top/bottom; step 8 `set {h: 0.1}` `watch: "(1/V)∮ stays 2.000: exact for a
    linear field"`; step 9 live "∇·Q = ∂Q₁/∂x₁ + ∂Q₂/∂x₂ = 1 + 1 = 2.000"; interpret "with your field and box: … ".
  - **D21** (7 steps) `view: 'limit'`; goal `set {preset: 'smooth', cx: 0.5, cy: 0.5, h: 2}`; step 2 (mean-value
    theorem) live "∬ ∇·Q dA = 2 cos x* cos y* × 4.000 for some x* in the box"; step 3 live "(1/V)∮ = …"; step 4 `set
    {h: 0.05}` `watch: "the dot slides down the purple curve onto the dashed line"`; step 7 (cross version) hint "E5
    does the same with a loop"; interpret "at h = … the ratio is within … % of ∇·Q(x₀)".
- **Code:**
  ```python
  Q = ch02.radial_field(1.0, dim=2)                        # preset: {{preset}} (a plane VectorField)
  box = ((cx - h/2, cx + h/2), (cy - h/2, cy + h/2))       # centre ({{cx}}, {{cy}}), side {{h}} m
  fluxes = ch02.flux_through_faces(Q, box, n=8)            # {'-x': {{fl}}, '+x': {{fr}}, '-y': {{fb}}, '+y': {{ft}}}
  rhs = sum(fluxes.values())                               # ∮ n·Q dA = {{rhs}} m³/s per m depth
  lhs, rhs2, _ = ch02.divergence_theorem_rect2d(Q, box, n=8)   # ∬ ∇·Q dA = {{lhs}}  (2.30)
  ratio = ch02.integral_divergence_2d(Q, [cx, cy], h, n_side=8) # (1/V)∮ = {{ratio}} → ∇·Q(x0) = {{exact}}  (2.32)
  tiles = ch02.divergence_theorem_tiled(Q, box, tiles=4, n=8)  # interior faces cancel: {{interior}} (NEW 5.3)
  ```
- **Walkthrough (8 steps):** 1. "Put a box in a field" — "The field is Q = x: it spreads from the origin. Drag the box
  around. The arrows on its faces show how much of Q leaves through each face." `set preset 'a x', {cx: 0.5, cy: 0.5,
  h: 1}`, highlight `view:field` · 2. "Add the faces" — "Right +1, top +1, left 0, bottom 0: net outflow 2. Now the
  purple bar: the divergence integrated over the inside — also 2. Always." `terms: true`, `eq: 'gauss'`, highlight
  `view:bars` · 3. "Why always? Cut it up" — "Tile the box 4×4. Every interior face is shared by two tiles with
  opposite normals: their fluxes cancel. Only the outer faces remain — and each tile is a 1-D fundamental theorem."
  `set {tiles: 4}`, `derive: {id: 'D25', step: 8}` · 4. "Read the faces with Taylor" — "Ex. 2.5: on the right face Q₁ =
  Q₁(centre) + (h/2)∂Q₁/∂x₁, on the left Q₁(centre) − (h/2)∂Q₁/∂x₁. Subtract: the centre value cancels, the slope
  survives." `set {tiles: 1, taylor: true}`, `derive: {id: 'D22', step: 4}`, `eq: 'taylor'` · 5. "Shrink it" — "Press
  shrink: the flux falls like h², the area like h², their ratio holds at 2 = ∇·Q. That ratio *is* the divergence — no
  coordinates needed (2.32)." `set {taylor: false}`, `play: true` (the h scrub), `eq: 'limit'`, highlight `view:limit` ·
  6. "A field that never leaks" — "b × x: the flow goes round. Every box, anywhere, any size: in = out, ∇·Q = 0.
  Incompressible flow looks like this." `set preset 'b × x (solenoidal)'`, `controls: ['cx', 'cy', 'h']` · 7. "When it
  fails" — "Point source: the field is singular at the origin. Enclose it: the flux is m however small the box. Gauss
  needs a smooth field inside." `set preset 'point source', {cx: 0, cy: 0, h: 1}`, `notes: true` · 8. "Your turn" —
  "Smooth preset. Predict: does (1/V)∮ overshoot or undershoot ∇·Q(x₀) at h = 2? Then shrink and read the slope of
  the error." `set preset 'smooth sin/cos', {cx: 0.5, cy: 0.5, h: 2}`, `controls: ['h', 'cx', 'cy']`, `readouts:
  ['ratio', 'exact']`.
- **Equations:** `gauss` "Gauss' theorem" ref 'Eq. (2.30)' `\iiint_V \frac{\partial Q}{\partial x_i}\,dV = \iint_A n_i
  Q\,dA` live "2.000 = 2.000" symbols [['Q', 'a field of any order', 'any'], ['n_i', 'outward unit normal', '–']] · `divthm`
  "Divergence theorem" ref '§2.12' `\iiint_V \nabla\cdot\mathbf Q\,dV = \iint_A \mathbf n\cdot\mathbf Q\,dA` live "∬∇·Q
  dA = 2.000, ∮ n·Q dA = 2.000 (unit depth)" · `gen` "Generalised derivative" ref 'Eq. (2.31)' `\mathcal D Q = \lim_{V\to
  0}\frac1V\iint_A n_i Q\,dA` · `limit` "Divergence as outflux per volume" ref 'Eq. (2.32)' `\nabla\cdot\mathbf Q = \lim_
  {V\to 0}\frac1V\iint_A \mathbf n\cdot\mathbf Q\,dA` live "(1/V)∮ = 2.000 at h = 1.000; ∇·Q(x₀) = 2.000" · `cart`
  "Cartesian form" ref 'Eq. (2.23)' `\nabla\cdot\mathbf u = \frac{\partial u_i}{\partial x_i}` live "1 + 1 = 2.000" ·
  `taylor` "Ex. 2.5 face values" ref '§2.12' `[\mathbf Q]_{\rm right} = \mathbf Q(\mathbf x) + \frac{\Delta x_1}{2}\frac{
  \partial\mathbf Q}{\partial x_1} + \dots` live "1.000 = 0.500 + 0.500" · `sphere` "Benchmark" ref 'Wikipedia example'
  `\oint_{|\mathbf x| = 1} (2x, y^2, z^2)\cdot\mathbf n\,dA = \tfrac{8\pi}{3} = 8.378` (static; the notebook computes it).
- **Check yourself:** (1) "a x preset, box centred at (0, 0) with h = 1. What is the net flux?" — "2.000 again: the
  divergence is 2 everywhere, so any unit box leaks 2 — even though the field is zero at the centre." `set {cx: 0, cy:
  0, h: 1}` · (2) "(x², 0) preset: (1/V)∮ at centre x₀ = 1 for h = 1 and h = 0.1?" — "2.000 both: ∇·Q = 2x = 2, and the
  midpoint face values are exact for a quadratic." `set preset '(x², 0)', {cx: 1, cy: 0}` · (3) "Solenoidal preset:
  which face fluxes are nonzero for a box at (1, 0)?" — "All four, but they cancel in pairs (right/left, top/bottom
  after signs): the total is 0 to round-off." · (4) "Point source with the origin *outside* the box. Flux?" — "0: away
  from the origin the source field is solenoidal (∇·Q = 0), so Gauss applies and the net flux vanishes."
- **Selftest parity rows:** `{name: 'right flux a x', js: faceFluxes(radial, [[0,1],[0,1]], 8).right, py: 'ch02.
  flux_through_faces(ch02.radial_field(1.0, 2), ((0.0, 1.0), (0.0, 1.0)), 8)["+x"]', rtol: 1e-12}` (WIP keys '+x', '-x', '+y', '-y'; the JS may name them right/left/top/bottom but maps them) · `{name: 'rhs sum',
  js: boxTheorem(radial, box, 8).rhs, py: 'ch02.divergence_theorem_rect2d(ch02.radial_field(1.0, 2), ((0.0, 1.0), (0.0,
  1.0)), 8)[1]', rtol: 1e-12}` · `{name: 'lhs smooth', js: boxTheorem(smooth, [[0,1],[0,1]], 16).lhs, py: 'ch02.
  divergence_theorem_rect2d(ch02.smooth_test_field(2), ((0.0, 1.0), (0.0, 1.0)), 16)[0]', rtol: 1e-10}` · `{name:
  'ratio (x²,0)', js: integralDiv(xsq, [1, 0], 0.2, 8), py: 'ch02.integral_divergence_2d(ch02.shear_field(0.0) …' — use
  a named preset: 'ch02.integral_divergence_2d(ch02.radial_field(1.0, 2), [0.5, 0.5], 0.2, 8)', expect 2.0, rtol: 1e-12}`
  · `{name: 'tiled interior', js: tiled(radial, box, 4, 8).interior, py: 'ch02.divergence_theorem_tiled(ch02.
  radial_field(1.0, 2), ((0.0, 1.0), (0.0, 1.0)), 4, 8)[2]', atol: 1e-12}` (needs the NEW 5.3 function) · `{name: 'solenoidal total', js: boxTheorem
  (solid, [[0.5,1.5],[-0.5,0.5]], 8).rhs, py: 'ch02.divergence_theorem_rect2d(ch02.irrotational_vortex_field(1.0), ((0.5, 1.5), (-0.5, 0.5)), 8)[1]', atol: 1e-10}` (the 2-D solenoidal field the WIP offers is the vortex with its core outside the box; **implementer note:** a 2-D `solid_body_rotation_field(b, dim=2)` would make this row and the E4/E5 presets cleaner) · invariant `{name: 'Gauss gap (smooth, n = 32)', js: |lhs −
  rhs|, expect: 0, atol: 1e-4}`.
- **Fit plan:** portrait phone: `field` (top, square) + `bars` (bottom, its title carrying "(1/V)∮ = 2.000 → 2.000");
  `limit` hidden; Explore shows preset, h, tiles (+ cx, cy optional — the box is draggable); Derivation keeps `field`
  (D25, D22) / `bars` (D21 on phones, with the ratio in its title). Landscape / notebook: three views (1.5 : 1 : 1).

### E5 · stokes_circulation_loop
- **Title:** "How much does the flow go round a loop?" · **Summary:** "Drag and resize a loop over a 2-D velocity field:
  the circulation ∮ u·t ds equals the curl integrated inside; a straight shear flow has it, an irrotational vortex has
  it only when the loop encloses its core; shrink the loop and Γ/A → (∇×u)·n." · **CORE:** C16, C11 (also N53, N54,
  N55, N69, N71, N73) · **Reference:** `templates/viz_example_field.html` (heatmap + streamlines + tracers + probe) with
  the `forced_damped_vibrations.html` explanation panel; the paddle wheel from `angular_frequency_explorer_1.html`.
- **meta:** `viz:sections 2.9 2.13` · `viz:equations 2.24 2.25 2.34 2.35` · `viz:fluidpy ch02.circulation ch02.curl_flux
  ch02.stokes_theorem_check ch02.integral_curl_component ch02.planar_loop ch02.rectangle_loop ch02.boundary_tangent
  ch02.curl ch02.solid_body_rotation_field ch02.shear_field ch02.irrotational_vortex_field ch02.potential_field` ·
  `viz:derivations D26 D12`.
- **Physics (JS ↔ Python):** field presets `u(x, y)` → [u₁, u₂] with exact `curlz(x, y)`: `solid` (−by, bx; curl 2b) ·
  `shear` (Γy, 0; curl −Γ) · `vortex` (K(−y, x)/r²; curl 0, singular at 0, circulation 2πK round the core) · `potential`
  (∇(x² − y²) = (2x, −2y); curl 0) · `source` (m(x, y)/(2πr²); curl 0) · `smooth` (sin x cos y, −cos x sin y; curl
  2 sin x sin y) ↔ `ch02.solid_body_rotation_field` (3-D; the JS 2-D copy is (−by, bx)), `shear_field`, `irrotational_vortex_field`, `potential_field()`,
  `point_source_field`, `smooth_test_field(2)` · `circleLoop(c, R, n)` → {points, tangents, ds} (counterclockwise) ↔
  `ch02.planar_loop(center, [0,0,1], radius, n)` · `rectLoop(c, a, b, n)` ↔ `ch02.rectangle_loop` · `circulation(u,
  loop)` (Σ u·t ds) ↔ `ch02.circulation(u_fn, loop)` · `curlFlux(curlz, loop)` (polar/rectangular midpoint sum of the
  exact curl over the interior) ↔ `ch02.curl_flux(u_fn, surface, curl_fn=…)` · `curlComponent(u, x0, h, n)` (square
  loop of side h; (2.35)) ↔ `ch02.integral_curl_component(u_fn, x0, [0,0,1], h, n_side)` · `stokesCheck(u, loop)` →
  {lhs, rhs, gap, ok} ↔ `ch02.stokes_theorem_check` · `fourSides(u, c, a, b)` (Ex. 2.6's midpoint sums, one per side)
  · `curlz25(u, x0, h)` (central differences of (2.25)'s third component) ↔ `ch02.curl` on a small grid.
- **Views:**
  1. `field` "The flow and your loop" (`equal: true`) — heatmap of (∇×u)₃ (purple positive, green negative, white 0)
     over [−2, 2]², streamlines, ~150 tracer particles advected by RK4 when the transport runs; the loop (circle of
     radius R or rectangle a × b at centre c) drawn purple with u·t arrows along it (tangent arrows, length ∝ u·t,
     reversed when negative) and the orientation arrowhead; a paddle wheel at the centre turning at ½(∇×u)₃(c); a ⚠️
     marker on the singular core when a preset has one. Pointer: drag inside the loop → move; drag its edge → resize;
     click a loop point → inspector.
  2. `ut` "u·t along the loop" — x: arc length s from 0 to the perimeter (m); y: u·t (m/s), purple curve; the signed
     area under it shaded (purple above 0, green below) = the circulation; the current pointer position marked;
     rectangle mode: the four sides as four shaded segments labelled bottom / right / top / left with their sums (Ex.
     2.6's bookkeeping).
  3. `bars` "∮ u·t ds vs ∬ (∇×u)·n dA, and Γ/A vs size" (`hidePortrait: true`; numbers repeated in view 2's title) —
     left: two bars (circulation purple, curl flux grey-hatched purple) with the gap printed; right: Γ/A vs loop size
     (log x) for the current centre with the dashed (∇×u)₃(c) and the current size as a dot; the singular-core preset
     draws Γ/A ∝ 1/A rising instead.
- **Controls:** `preset` chips "b × x (Ex. 2.3)" · "shear u₁ = Γx₂" · "irrotational vortex u_θ = K/r" · "potential ∇φ" ·
  "point source" · "smooth sin/cos" · `cx`, `cy` "Loop centre" −2…2 step 0.05 (default 0, 0) m · `R` "Radius $R$ (or
  half-side)" log 0.05…1.8 (default 1) m · `shape` chips "circle / rectangle (Ex. 2.6)" · `flip` toggle "flip n (→ −e₃):
  reverse the orientation" · `t` transport 0…10 s for tracers and the wheel (rate 1, loop) *optional* · `aspect`
  "rectangle b/a" 0.25…4 (default 1) *optional*.
- **Depth features:** Explain + Code · **linked views** (3) · **presets** (six fields + "shrink the loop" {R: 0.05} +
  "enclose the core" {preset: 'vortex', cx: 0, cy: 0, R: 1} + "miss the core" {preset: 'vortex', cx: 1.5, cy: 0, R:
  0.4}) · **status** ("🌀 positive circulation Γ = 6.283 m²/s: net counterclockwise spin inside (curl 2.00)" · "↔ shear
  flow: straight streamlines, yet Γ = −3.142 — the wheel turns clockwise" (preset shear) · "0️⃣ irrotational: Γ → 0 as the
  loop shrinks (Γ/A = 0.000)" (curl-free preset, core outside) · "⚠️ singular core inside the loop: Stokes does not
  apply — Γ = 2πK = 6.283 whatever the size" (vortex with the core inside) · "↩ orientation flipped: both sides changed
  sign together") · **term bars** (rectangle mode: bottom, right, top, left → Γ; circle mode: quarter arcs) ·
  **inspector** (click the loop at s: "x = (0.707, 0.707): u = (−0.707, 0.707) m/s, t = (−0.707, 0.707), u·t = 1.000
  m/s, ds = 0.0157 m → contribution 0.0157 m²/s") · **transport** (tracers and the wheel; end 'loop').
- **Explain:**
  0. *What the windows show* — "**The flow** is coloured by the z-component of its curl (purple = counterclockwise
     spin, green = clockwise, white = none); grey streamlines and moving dots show the velocity itself. The **purple
     loop** is yours; the arrows along it are the velocity component *along* the loop, u·t. The **wheel** at the centre
     turns at half the local curl. **u·t along the loop** unrolls the loop into a graph; its shaded area is the
     circulation. **Bars** compare the circulation with the curl added up inside."
  1. *The field at the loop's centre and its curl (2.25)* — "u = b × x with b = 1: u(0, 0) = (0, 0); (∇×u)₃ = ∂u₂/∂x₁
     − ∂u₁/∂x₂ = 1 − (−1) = **2.000** 1/s — twice the angular velocity (Ex. 2.3). Central differences with h = 0.01 give
     2.000."
  2. *The circulation, added up round the loop (2.34, right side)* — "circle: t = (−sin θ, cos θ), u·t = b R (sin²θ +
     cos²θ) = 1.000 m/s at every point; ∮ u·t ds = 1.000 × 2πR = **6.283** m²/s." Rectangle mode instead: "bottom (t =
     +e₁, y = −0.5): Σ u₁ ds = Γ(−0.5)(1.0) = … ; right (t = +e₂): …; top (t = −e₁): …; left (t = −e₂): …; sum = **…**
     — Ex. 2.6's four sides."
  3. *The curl added up inside (2.34, left side)* — "∬ (∇×u)·n dA = 2.000 × πR² = 2.000 × 3.142 = **6.283** m²/s. Gap
     to the circulation: **1e-13** (a 400-point midpoint loop on a linear field is exact to round-off; for the smooth
     preset the gap falls as 1/n²)."
  4. *Circulation per unit area (2.35)* — "Γ/A = 6.283/3.142 = **2.000** 1/s = (∇×u)₃(centre). For the smooth preset
     the ratio only converges as R → 0 — follow the curve in the bars view; for the vortex with the core inside it
     grows like 1/A (no limit)."
  5. *Orientation* — "n = +e₃ (out of the screen): t runs counterclockwise (right hand: thumb along n, fingers along t);
     n_c = n × t points into the loop; flipping n reverses t and both sides of (2.34): Γ → **−6.283**, flux → −6.283."
  6. *The wheel* — "turns at ½(∇×u)₃ = **1.000** rad/s counterclockwise; in the shear preset it turns clockwise at Γ/2
     although no streamline is curved — curl is a difference of speeds across the wheel, not curvature."
  7. *At the current time* — live: t, the tracers' mean angular displacement, the wheel's angle, Γ, flux, Γ/A.
  8. *Reading the current setting* — solid body: "the whole plane spins like a disc; every loop, anywhere, has Γ = 2b
     × its area — the spin is uniform." · shear: "straight flow, faster above: a loop is pushed along on top and
     dragged back below, so it has net clockwise circulation −Γ × area; **curl without curves**." · irrotational, core
     outside: "the flow goes round the core, yet inside your loop there is no spin at all: Γ → 0 as the loop shrinks;
     the dots orbit but do not rotate about themselves." · vortex, core inside: "Γ = 2πK for every loop enclosing the
     core, big or small — the curl is all concentrated in the singular centre; Stokes' theorem needs a differentiable
     field inside the loop, and this one is not. Ch. 5's line vortex, Ch. 6's lift." · potential flow: "u = ∇φ: ∮∇φ·t ds
     = ∮dφ = 0 for every loop (D26 corollary, Exercise 2.20): ∇×∇φ = 0." · smooth: "the curl changes sign inside a big
     loop; positive and negative regions partly cancel; shrink to isolate one sign."
- **Derivation tab:**
  - **D26** (10 steps) `view: 'field'`; goal `set {preset: 'shear', shape: 'rect', cx: 0, cy: 0, R: 0.5, flip: false}`;
    step 1 highlights the four sides with their tangents; step 3 **live** "bottom + top: [u₁(y = −0.5) − u₁(y = +0.5)]
    Δx = (−0.5 − 0.5) × 1.0 = −1.000"; step 5 live "right + left: [u₂(x = +0.5) − u₂(x = −0.5)] Δy = 0"; step 6 live
    "Γ_rect = (∂u₂/∂x₁ − ∂u₁/∂x₂) ΔA = (0 − 1) × 1.0 = −1.000 ✓"; step 7 `set {tiles: 4}` (the rectangle's interior
    tiled in the field view with shared edges drawn twice in opposite directions) `watch: "every interior edge is
    walked twice, opposite ways"`; step 9 `set {shape: 'circle'}` `watch: "the staircase boundary of the tiles tends to
    the circle"`; step 10 `set {preset: 'potential'}` live "∮ ∇φ·t ds = 0.000"; interpret "for your loop: Γ = … = curl
    flux …".
  - **D12** (5 steps) `view: 'field'`; goal `set {preset: 'solid', R: 0.5}`; step 5 (i = 3) **live** "(∇×u)₃ = ∂u₂/∂x₁ −
    ∂u₁/∂x₂ = 1.000 − (−1.000) = 2.000 at the centre" `highlight ['readout:curl']`; interpret "the wheel turns at half
    of this".
- **Code:**
  ```python
  u = ch02.solid_body_rotation_field([0, 0, 1.0])            # preset: {{preset}}
  loop = ch02.planar_loop([cx, cy], None, R, n=400)          # 2-D loop, counterclockwise about +e3; R = {{R}} m
  Gamma = ch02.circulation(u, loop)                          # ∮ u·t ds = {{Gamma}} m²/s   (2.34, right)
  disc = ch02.planar_disc([cx, cy], None, R, 40, 80)
  flux = ch02.curl_flux(u, disc, curl_fn=lambda p: [0, 0, 2.0])   # ∬ (∇×u)·n dA = {{flux}}   (2.34, left)
  check = ch02.stokes_theorem_check(u, loop, disc)           # lhs − rhs = {{gap}}, hypothesis_ok = {{ok}}
  w3 = ch02.integral_curl_component(u, [cx, cy], None, h=0.1)  # (2.35): {{w3}} 1/s
  rect = ch02.rectangle_loop([cx, cy], None, 2*R, 2*R, n=100)  # Ex. 2.6's rectangle
  ```
- **Walkthrough (8 steps):** 1. "A loop in a whirlpool" — "The plane spins like a record. Put a loop anywhere and add
  up the velocity along it: that sum is the circulation. Play to see the dots orbit." `set preset 'b × x', {cx: 0, cy:
  0, R: 1}`, `play: true`, highlight `view:field` · 2. "Unroll the loop" — "Along this circle u·t = 1 m/s everywhere;
  the shaded area is 1 × 2π = 6.283 m²/s." `eq: 'circ'`, highlight `view:ut`, `inspect: true` · 3. "Stokes: it equals
  the curl inside" — "The curl is 2 everywhere (Ex. 2.3); 2 × π = 6.283. Circulation round the edge = total spin
  inside — the two bars match." `eq: 'stokes'`, `derive: {id: 'D12', step: 5}`, highlight `view:bars` · 4. "Straight
  flow, still spinning" — "Shear: every streamline is straight, yet the wheel turns and Γ ≠ 0. Faster on top, slower
  below — that difference is curl." `set preset 'shear u₁ = Γx₂'`, `notes: true` · 5. "One rectangle proves it" —
  "Rectangle mode is Ex. 2.6: bottom + top give −∂u₁/∂x₂ ΔA, right + left give +∂u₂/∂x₁ ΔA. Tile any surface with
  rectangles and the interior edges cancel." `set {shape: 'rect', R: 0.5}`, `terms: true`, `derive: {id: 'D26', step:
  6}` · 6. "Shrink it" — "Circulation per unit area → the normal curl (2.35). Drag R down and watch Γ/A hold at the
  point value." `set {shape: 'circle'}`, `controls: ['R']`, `eq: 'per_area'` · 7. "When Stokes fails" — "Irrotational
  vortex: curl zero everywhere — except at the core. Miss the core: Γ → 0. Enclose it: Γ = 2πK at any size. The field
  is not smooth inside, so the theorem does not apply." `set preset 'irrotational vortex', {cx: 0, cy: 0, R: 1}`,
  `readouts: ['Gamma', 'ratio']` · 8. "Your turn: potential flow" — "Predict the circulation of u = ∇φ round any loop.
  Then flip the orientation: what happens to both sides?" `set preset 'potential ∇φ'`, `controls: ['cx', 'cy', 'R',
  'flip']`.
- **Equations:** `curl24` "Curl, index form" ref 'Eq. (2.24)' `(\nabla\times\mathbf u)_i = \varepsilon_{ijk}\frac{\partial
  u_k}{\partial x_j}` · `curl25` "z-component" ref 'Eq. (2.25)' `(\nabla\times\mathbf u)_3 = \frac{\partial u_2}{\partial
  x_1} - \frac{\partial u_1}{\partial x_2}` live "= 1.000 − (−1.000) = 2.000 1/s" · `orient` "Orientation" ref '§2.13'
  `\mathbf t = \mathbf n_c\times\mathbf n,\quad \mathbf n_c = \mathbf n\times\mathbf t` note "t counterclockwise about n
  (right hand); n_c points along the surface into A" · `stokes` "Stokes' theorem" ref 'Eq. (2.34)' `\iint_A (\nabla\times
  \mathbf u)\cdot\mathbf n\,dA = \oint_C \mathbf u\cdot\mathbf t\,ds` live "6.283 = 6.283 m²/s" symbols [['\\Gamma =
  \\oint \\mathbf u\\cdot\\mathbf t\\,ds', 'circulation', 'm²/s']] · `per_area` "Curl as circulation per area" ref 'Eq.
  (2.35)' `\mathbf n\cdot(\nabla\times\mathbf u) = \lim_{A\to 0}\frac1A\oint_C \mathbf u\cdot\mathbf t\,ds` live "Γ/A =
  2.000 1/s" · `ex26` "Ex. 2.6, one rectangle" ref '§2.13' `\oint_{\rm rect}\mathbf u\cdot\mathbf t\,ds \approx [u_2(x_1
  + \tfrac{\Delta x_1}{2}) - u_2(x_1 - \tfrac{\Delta x_1}{2})]\Delta x_2 - [u_1(x_2 + \tfrac{\Delta x_2}{2}) - u_1(x_2 -
  \tfrac{\Delta x_2}{2})]\Delta x_1` live "0 − (1.000) = −1.000" note "the book's x-plane version has u_y in the second
  bracket" · `grad0` "Gradient fields" ref 'Exercise 2.20' `\oint_C \nabla\phi\cdot\mathbf t\,ds = 0 \Rightarrow \nabla
  \times\nabla\phi = 0`.
- **Check yourself:** (1) "b × x with R = 0.5 at centre (1, 1). Circulation?" — "2 × π × 0.25 = 1.571 m²/s: the curl is
  uniform, so only the area matters, not where the loop is." `set {cx: 1, cy: 1, R: 0.5}` · (2) "Shear flow Γ = 1, unit
  square loop. Sign and size of the circulation?" — "−1.000: top side pushed +x at speed 0.5 but walked in −x, bottom
  walked +x at speed −0.5; both negative; the vertical sides add nothing." `set preset 'shear', {shape: 'rect', R:
  0.5}` · (3) "Vortex K = 1: loop of radius 0.3 at (1.5, 0) vs radius 1 at the origin?" — "≈ 0 vs 6.283 = 2πK: the curl
  lives only in the core; Γ counts whether the core is inside." · (4) "Flip n. What changes?" — "t reverses, so u·t
  and Γ change sign — and so does (∇×u)·n; the theorem holds either way." `set {flip: true}`.
- **Selftest parity rows:** `{name: 'Γ solid R = 1', js: circulation(solid, circleLoop([0,0], 1, 400)), py: 'ch02.
  circulation(ch02.solid_body_rotation_field([0.0, 0.0, 1.0]), ch02.planar_loop([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0,
  400))', rtol: 1e-10}` · `{name: 'Γ shear rect', js: circulation(shear, rectLoop([0,0], 1, 1, 100)), py: 'ch02.
  circulation(ch02.shear_field(1.0), ch02.rectangle_loop([0.0, 0.0], None, 1.0, 1.0, 100))', rtol:
  1e-10}` (2-D field → 2-D loop; probed −1.0) · `{name: '(2.35) shear', js: curlComponent(shear, [0.3, 0.2], 0.1, 8), py: 'ch02.integral_curl_component
  (ch02.shear_field(1.0), [0.3, 0.2], None, 0.1, 8)', rtol: 1e-10}` (probed −1.0) · `{name: 'Γ vortex core inside',
  js: circulation(vortex, circleLoop([0,0], 0.7, 400)), py: 'ch02.circulation(ch02.irrotational_vortex_field(1.0),
  ch02.planar_loop([0.0, 0.0], None, 0.7, 400))', rtol: 1e-8}` (2-D loop) · `{name: 'Γ potential ≈ 0', js:
  circulation(potential, circleLoop([0.5,0.3], 0.8, 400)), py: 'ch02.circulation(ch02.potential_field(),
  ch02.planar_loop([0.5, 0.3], None, 0.8, 400))', atol: 1e-10}` (the WIP default φ = x₁² − x₂²; probed 0.0) · `{name: 'tangent', js: tangentOf([0,0,1], [1,0,0])[1], py: 'ch02.boundary_tangent([-1.0,
  0.0, 0.0], [0.0, 0.0, 1.0])[1]', rtol: 1e-12}` (n_c inward = −r̂ at the point (1, 0): t = n_c × n = (−1,0,0) × (0,0,1) =
  (0, 1, 0) ✓) · `{name: 'hypothesis flag', js: stokesCheck(vortex, circleLoop([0,0], 1, 400)).ok ? 1 : 0, py: '1 if
  … ' — no builtins: py: 'ch02.stokes_theorem_check(ch02.irrotational_vortex_field(1.0), ch02.planar_loop([0.0, 0.0],
  None, 1.0, 400), ch02.planar_disc([0.0, 0.0], None, 1.0, 20, 40)).hypothesis_ok', and the JS side yields `false`;
  shot.py converts booleans with float → 0.0; atol: 0}` (probed: False, lhs 1.5e-14, rhs 6.283) · invariant `{name: 'Stokes gap
  smooth', js: |Γ − flux| for smooth, R = 0.8, expect: 0, atol: 1e-6}` · `{name: 'flip changes sign', js: Γ(flip) + Γ,
  expect: 0, atol: 1e-12}`.
- **Fit plan:** portrait phone: `field` (top, square) + `ut` (bottom; title "Γ = 6.283 = flux 6.283 · Γ/A = 2.000");
  `bars` hidden; the transport hidden on the Derivation tab; Explore shows preset, R, shape, flip (+ cx, cy optional —
  the loop is draggable). Landscape / notebook: three views (1.5 : 1 : 1). Tracer count 150 → 80 below 500 px wide.

### B1 · index_machine
- **Title:** "What does a repeated index actually do?" · **Summary:** "Pick an index expression and watch the
  summation convention unfold it into its 3, 9 or 27 terms, each term lighting the matrix or ε-cube cells it multiplies,
  with your numbers plugged in." · **CORE:** C01, C08 (also N12, N21, N23, N41, N44, N49, N74) · **Reference:**
  `stride_padding_playground.html` (a formula with numbers plugged in, badges, click a cell) + ch01 E5's matrix view.
- **meta:** `viz:sections 2.1 2.3 2.7 2.8 2.14` · `viz:equations 2.2 2.9 2.17 2.18 2.19 2.21 2.12 2.36` · `viz:fluidpy
  ch02.expand_indices_str ch02.classify_indices ch02.tensor_order ch02.levi_civita ch02.epsilon_delta_residual
  ch02.cross_einsum ch02.double_dot ch02.dot ch02.inner` · `viz:derivations D09`.
- **Physics (JS ↔ Python):** `parse(term)` → operands with index letters ↔ `ch02.classify_indices` · `expandStr(term,
  dim)` → the sum as a string (same spacing as Python) ↔ `ch02.expand_indices_str(term, dim)` · `order(term)` ↔
  `ch02.tensor_order` · `evaluate(term, values)` (numeric value(s)) ↔ `ch02.dot`, `ch02.inner`, `ch02.cross_einsum`,
  `ch02.double_dot`, `ch02.transform_tensor` for the respective chips · `eps(i, j, k)` ↔ `ch02.levi_civita()` · `epsDelta
  (i, j, l, m)` residual ↔ `ch02.epsilon_delta_residual`.
- **Views:** 1. `expr` "The expression and its expansion" — the chip's expression with free indices in blue and dummy
  letters in orange; below, the expanded sum with the current term highlighted (transport steps through terms; 3, 9 or
  27 of them); badges "order 0/1/2", "free: j", "summed: i". 2. `cells` "Operand cells" — the operand matrices (2×2 or
  3×3) and, for ε terms, the unfolded 3×9 ε grid, with the cells of the current term lit (row/column highlight for
  A_ik B_kj as in (2.11)). 3. `value` "With your numbers" (`hidePortrait: true`) — the numeric result (scalar, vector as
  three bars, matrix as a 3×3 heat grid) and the running partial sum as the terms are stepped.
- **Controls:** `term` chips `a_i b_i` · `A_ik B_kj` · `delta_ij u_j` · `eps_ijk u_i v_j` · `eps_ijk eps_klm` · `C_im
  C_jn tau_ij` · `u_i,i` · `dim` toggle 2 / 3 · `dummy` select "rename the summed letter to k / m / p" · operand presets
  (a = (1, 2, 3), b = (4, 5, 6), A = 1…9, C = 30° rotation, τ = pure shear) · `step` transport over the terms (rate 1
  term/s, end 'hold').
- **Depth features:** Explain + Code · **transport** (term by term) · **presets** · **inspector** (click a term: its
  factors and product) · **modes** (dim 2 / 3) · **status** ("order 2: two free indices m, n → 9 numbers" / "scalar: no
  free index" / "❌ an index appears three times — not a valid term").
- **Explain:** 0. what the colours mean (blue free, orange summed); 1. "Your expression": free and dummy letters,
  order; 2. "The hidden sum": the expansion with dim terms per dummy letter (3 for one, 9 for two, 27 for three);
  3. "With your numbers": each term's value and the total (a·b = 4 + 10 + 18 = 32); 4. "Renaming": the same total with
  the dummy renamed; 5. "ε terms" (or hint): which of the 27 cells are nonzero and why only two survive per free index;
  6. "Reading the current setting": "a repeated index is a loop; the free indices are the shape of the answer" with
  the current shape.
- **Derivation tab:** **D09** (11 steps) `view: 'cells'`; goal `set {term: 'eps_ijk eps_klm', dim: 3}`; steps 4–7 step the
  ε grid through the cases with the (i, j, l, m) chips; step 8 live "ε_pqi ε_pqj summed: 2δ_ij → 2 on the diagonal";
  step 9 live "= 6"; step 11 `set {term: 'eps_ijk u_i v_j'}` `watch: "a × (b × c) with your vectors"`.
- **Code:** `print(ch02.expand_indices_str("a_i b_i"))  # {{expansion}}` · `ch02.dot(a, b)  # {{value}}` ·
  `ch02.classify_indices("x_i C_ij")  # {{free}}, {{dummy}}` · `ch02.tensor_order("C_im C_jn tau_ij")  # {{order}}` ·
  `ch02.epsilon_delta_residual()  # {{res}}` · `ch02.cross_einsum(u, v)  # {{cross}}`.
- **Walkthrough (6 steps):** 1. "Why no Σ?" (a_i b_i, step the three terms) · 2. "Free vs summed" (A_ik B_kj: k walks
  along a row and down a column, (2.11)) · 3. "δ eats an index" (delta_ij u_j → u_i) · 4. "ε: 27 cells, 6 alive"
  (eps_ijk u_i v_j, the k = 1 check) · 5. "Two C's for a tensor" (C_im C_jn tau_ij with the shear preset) · 6. "Your
  turn: comma" (u_i,i — predict the number of terms).
- **Equations:** (2.2), (2.9), (2.17), (2.18), (2.19), (2.21), (2.12), (2.36) with live expansions.
- **Check yourself:** (1) "How many terms does C_im C_jn τ_ij have for one (m, n) in 3-D?" — "9 (i and j each run
  1..3); 81 numbers computed in total for the 9 (m, n) pairs." · (2) "Is a_i b_i c_i a valid term?" — "No: i appears
  three times; the convention sums pairs only." · (3) "In 2-D, how many nonzero ε_ij terms are there?" — "Two: ε₁₂ = 1,
  ε₂₁ = −1 (the 2-D alternating symbol)."
- **Selftest parity rows:** `{name: 'expansion string', js: expandStr('a_i b_i', 3) === 'a_1*b_1 + a_2*b_2 + a_3*b_3'
  ? 1 : 0, expect: 1, atol: 0}` (the Python string is fixed in the contract) · `{name: 'a·b', js: evaluate('a_i b_i'),
  py: 'ch02.dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])', rtol: 1e-12}` · `{name: 'order', js: order('C_im C_jn tau_ij'), py:
  'ch02.tensor_order("C_im C_jn tau_ij")', atol: 0}` · `{name: 'ε–δ residual', js: epsDeltaResidual(), py: 'ch02.
  epsilon_delta_residual()', atol: 0}` · `{name: 'cross k=1', js: evaluate('eps_ijk u_i v_j')[0], py: 'ch02.cross_einsum
  ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])[0]', rtol: 1e-12}`.
- **Fit plan:** portrait: `expr` + `cells`; `value` hidden (its number in the `expr` badge).

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

This chapter is cheap: no ODE solves, no long time loops. The heavy items are the 3-D grids (48³ = 110 592 points ×
3 components × a handful of stencil passes), the two convergence studies, the sphere quadrature, the D17 sympy cell,
five matplotlib animations and six plotly slider figures. Estimates from the ch01 timings (≈ 0.12 s per animation frame
at dpi 80; a 3-D stencil pass on 48³ ≈ 15 ms; sympy 2×2 symbolic eigenvectors ≈ 1 s).

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly/pint imports | 8 s | 8 s |
| §2.1 C01 | sympy expansions (7 calls), matrix-product figure, 3-D position-vector plotly | 3 s | 3 s |
| §2.2 C02, C03 | 3-D two-frame plotly, frame-rotation animation (10 frames, frames player), 200-rotation orthogonality loop, residual table (20 points × 5 fns) | 5 s | 4 s (6 frames) |
| §2.3 | three asserts | < 1 s | < 1 s |
| §2.4, 2.6 C04–C06 | stress-cube plotly (9 cones), tetrahedron plotly, two slider figures (37 steps × ≤ 4 traces × ≤ 200 points each), the (2.12) residual test | 6 s | 5 s (25 steps) |
| §2.5 C07 | 50-rotation invariants figure, sympy Vieta line | 2 s | 2 s |
| §2.7, 2.8 C08 | 81-case quadruple loop (trivial), isotropy tests (4 × 50 rotations, order 2–3 arrays), ε grid figure | 2 s | 2 s |
| §2.9 C09–C11 | 41×41 2-D grid ops; 48³ 3-D grid: radial, rotation, shear fields, `divergence`, `curl`, `vector_gradient` (≈ 12 stencil passes ≈ 0.3 s); `operator_convergence` × 3 ops on n = (8, 16, 32, 64) (≈ 1.5 s, cached); heatmap figures; tracer video (40 frames); contour slider (36 steps) | 14 s | 8 s (24³ grid, levels [8, 16, 32], 24 frames, 18 slider steps) |
| §2.10 C12 | expm on 2×2 (µs), three-panel video (48 frames), 3-D vector-gradient slice from §2.9's cached arrays | 8 s | 4 s (24 frames) |
| §2.11 C13 | D17 sympy check (symbolic 2×2 eigenvectors + two simplifies ≈ 2 s), Monte-Carlo bounds (1000 normals), Γ slider (30 steps), live widget (not executed interactively) | 4 s | 3 s (300 normals, 15 steps) |
| §2.12 C14, C15 | `divergence_theorem_box` n = 16 (µs–ms), sphere n = 48 (48 × 96 × 3 faces ≈ 50 ms), mesh slider (8 steps, n up to 32 ≈ 0.2 s), `integral_definition_convergence` × 3 kinds on 4 levels (≈ 0.3 s, cached), face-bookkeeping frames player (12 frames), h slider (20 steps), coloured-box plotly | 7 s | 5 s (sphere n = 24, levels 3) |
| §2.13 C16 | loops of 200–400 points (ms), disc 20 × 40, shrinking-loop frames player (10 frames), cap plotly | 4 s | 3 s |
| §2.14 + end | sympy comma expansions | < 1 s | < 1 s |
| explainer cells (5 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 65 s** | **≈ 49 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: 3-D grids 48³ → 24³;
`operator_convergence` levels (8, 16, 32, 64) → (8, 16, 32) and `integral_definition_convergence` (0.4, 0.2, 0.1, 0.05) → (0.4,
0.2, 0.1); sphere quadrature 48 × 96 → 24 × 48; videos 40–48 frames → 24; frames players 10–12 frames → 6–8; slider figures
≤ 37 steps → ≤ 25, ≤ 4 traces × ≤ 400 points (each < 250 kB); Monte-Carlo normals 1000 → 300. **Cached arrays:** the 3-D
grid `g3` and the three fields `u_rad, u_rot, u_sh` are computed once in C10's first code cell and reused by C11 and C12;
the two convergence studies are computed once (C09 grid, C15 integral) and their figures reuse the dicts (C09, C15,
C16); the plotly 3-D figures use ≤ 18 cones and ≤ 200 mesh triangles (each < 300 kB). Outputs: 2 videos (< 2 MB each), 3
frame players (< 3 MB each), 6 sliders, 5 plotly 3-D figures — page well under 15 MB.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is
explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use (the primer term is the Concept text
before any parenthesis — the builder uses it verbatim in `nb.primer`); "gloss" = one sentence where used; "knowledge/
primers.md: <term> (ch01)" = a one-line reminder naming the ch01 primer. IDs are curation IDs (A items explain
themselves; B/C items are explained by their note, named in parentheses).

| Concept | First used in | Explained by |
|---|---|---|
| scalar, vector, second-order tensor (orders 0, 1, 2; 1 → 3 → 9 numbers) | §2.1 (before C01) | C01 (N01, N02, N06, N24 notes) |
| index letters i, j, k ∈ {1, 2, 3} (Python 0–2) and the ≡ symbol | §2.1 front matter | C01 (gloss in the notation cell, A.0 #4) |
| summation convention (repeated index = sum; free vs dummy index) | C01 | C01 |
| position vector in a basis, column and transpose notation | C01 | C01 (N03, N04 notes) |
| dot product as an index sum (2.2); u·v = uv cos θ | C01 | C01 (N07, N45 notes) |
| np.einsum index strings ('i,i', 'ik,kj->ij', 'ijk,i,j->k') | C01 | primer (in C01) |
| matrix multiplication, transpose and identity (`@`, `.T`, `np.eye`) | C01 | primer (in C01) |
| matrix product as an index sum (2.9)–(2.11); single dot = one index summed | C01 | C01 (N21, N22, N23 notes) |
| Kronecker delta δ_ij (2.16) and its substitution rule (2.17); δ_ii = 3 | C01 | C01 (N40, N41 notes) |
| boldface vs indicial notation; order = number of free indices | C01 | C01 (N08 note) |
| `for` loops as the explicit form of a sum (triple loop) | C01 from-scratch | C01 (comments in the check cell) |
| plotly 3-D arrows, lines and meshes (`go.Cone`, `go.Scatter3d` lines, `go.Mesh3d`) | C01 (N05 figure) | primer (in C01, A.1 #18a) |
| rotated frame with the same origin; primed components (2.3) | C02 | C02 (N09, N10 notes) |
| direction-cosine matrix C_ij = e_i·e'_j; columns = new axes | C02 | C02 |
| orthonormal basis, projection and completeness (Σ_j e'_j e'_jᵀ = I) | C02 (D02) | primer (in C02) |
| cosines of angles between unit vectors; cos(π/2 − θ) = sin θ; radians and `np.deg2rad` | C02 | primer (in C02) |
| np.linalg.norm and np.linalg.qr (random rotations) | C02 | primer (in C02) |
| orthogonality of C, C Cᵀ = Cᵀ C = I, det C = ±1 | C02 | C02 (N15 note, D02) |
| determinant of a product, det(AB) = det A det B; 3×3 by cofactors | C02 (D02), C07 | knowledge/primers.md: matrices, determinants and minors (ch01 P53) — reminder in D02's *why* |
| continuity argument for det C = +1 (a rotation connects to the identity) | C02 (D02 step 8) | C02 (gloss in D02 step 8's *why*) |
| Rodrigues rotation about an axis (`rotation_matrix_3d`) | C02 | C02 (code comment; the 2-D case is the tiny example) |
| passive (axes turn) vs active (vector turns) rotation | C02 | C02 (⚠️ callout A.2 #17) |
| inverse transformation (2.7): the summed slot of C moves | C02 | C02 (N14 note, D03) |
| projection of x on e'_1 (2.4) | C03 | C03 (N11 note) |
| bilinearity of the dot product (distributes over sums; scalars pull out) | C03 (D01 step 2) | C03 (gloss in D01's *why*) |
| transformation of components (2.5), (2.6), matrix forms x' = Cᵀx, x = Cx' | C03 | C03 (N13, N16 notes, D01) |
| formal definition of a vector (2.8); pass/fail residual test | C03 | C03 (N18 note) |
| Python `dict` of `lambda` functions as test candidates | C03 (#31) | knowledge/primers.md: functions as arguments and lambda (ch01 P29); Python dictionaries (ch01 P23) — reminders |
| polar components, Ex. 2.1 | C03 | C03 (N19 note) |
| stress tensor τ_ij: first index face, second index force; sign convention | C04 | C04 |
| stress as force per area, normal vs shear | C04 | knowledge/primers.md: stress (ch01 P05) — reminder |
| Newton's third law at a point: opposite faces carry opposite stresses | C04 | R01 (recap inside C04) |
| free-body diagram and Newton's second law | C05 (D05) | knowledge/primers.md: Newton's second law and momentum (ch01 P09) — reminder in D05's tools |
| surface element as a vector dA = n dA | C05 | C05 (N34 note) |
| tetrahedron geometry dA_i = n_i dA | C05 | C05 (N35 note); primer vector area (in C05) |
| vector area of a closed surface (Σ n dA = 0) | C05 (D05) | primer (in C05) |
| limits and orders of smallness (h² vs h³) | C05 (D05), C15 (D22) | primer (in C05) |
| Cauchy's traction formula (2.15); n·τ vs τ·n | C05 | C05 (N37 note, ⚠️ callout) |
| np.arctan2 (four-quadrant angle) | C05 | primer (in C05) |
| normal and shear parts of a traction, σ_n = f·n, τ_s = |f − σ_n n| | C05 | C05 (tiny example, N38 note); `normal_shear_stress` |
| double-angle forms σ_n = a sin 2φ, τ_s = a cos 2φ (Mohr preview) | C05 | C05 (N38 note, stated) |
| tensor transformation rule (2.12), τ' = CᵀτC; tensor ≠ matrix | C06 | C06 (N27 note, D06) |
| "true for every n ⇒ the coefficients agree" | C06 (D06 step 7) | C06 (gloss in D06's *why*) |
| fourth-order transformation (2.13) | C06 | C06 (N28 note) |
| outer product u_i v_j as a tensor | C06 | C06 (N29 note) |
| contraction; trace; invariants I₁, I₂, I₃; closed index chains | C07 | C07 (D18) |
| Vieta's formulas (cubic coefficients ↔ roots) | C07 (D18) | primer (in C07) |
| characteristic polynomial det(A − λδ) = 0 | C07 | C07 (D18 steps 6–7; primer eigenvalues in C13 for the meaning) |
| the four contractions (2.14) as matrix products; A·u vs Aᵀ·u | C07 | C07 (N31, N32 notes) |
| double dot A:B, book vs Frobenius convention | C07 | C07 (N33 note, ⚠️ callout) |
| tensor product raises the order (`np.multiply.outer`) | C07 | C07 (N30 note) |
| np.random.default_rng (seeded random rotations) | C07 | knowledge/primers.md: np.random.default_rng (ch01 P10) — reminder |
| alternating tensor ε_ijk (2.18); index moves | C08 | C08 (N43 note) |
| permutations, cyclic order and parity (`itertools.permutations`) | C08 | primer (in C08) |
| numpy arrays with three axes and np.transpose | C08 | primer (in C08) |
| right-hand rule and orientation | C08, C16 (N69) | primer (in C08) |
| cross product: geometric definition, components (2.20), determinant form, index form (2.21) | C08 | C08 (N46–N50 notes) |
| epsilon–delta relation (2.19), contractions 2δ and 6, vector triple product | C08 | C08 (N44 note, D09) |
| isotropic tensors; ε as a pseudotensor under reflections | C08 | C08 (N42 note) |
| `itertools.product` (quadruple loop over 81 cases) | C08 from-scratch | C08 (code comment; sibling of ch01 P55 `combinations`) |
| the del operator ∇ (2.22) in index form | C09 | R02 (Ch. 1 §1.5 P25) |
| partial derivative ∂/∂x_i | C09 | knowledge/primers.md: partial derivative (ch01 P25) — reminder in R02 |
| finite differences; central and one-sided second-order stencils | C09 | knowledge/primers.md: finite differences (ch01 P21), np.gradient (ch01 P22) — reminder in R02; the one-sided stencil written out in A.7 #10 |
| level sets and the directional derivative ∂φ/∂n = ∇φ·n | C09 | primer (in C09) |
| chain rule (∇φ ⊥ level sets; ∂/∂x'_j = C_ij ∂/∂x_i; dφ along a curve) | C09, §2.14 (N74), C16 (D26 step 10) | knowledge/primers.md: chain rule (ch01 P49) — reminder |
| gradient ∇φ: perpendicular to level sets, steepest climb | C09 | C09 |
| np.meshgrid and the project grid layout ([k, j, i] = (z, y, x); components on axis 0) | C09 | primer (in C09) |
| numpy broadcasting (grid arrays, neighbour slices) | C09 | primer (in C09) |
| plt.contour, plt.quiver and plt.streamplot | C09 | primer (in C09) |
| observed order of convergence (log–log slope; `tools.convergence.observed_order`) | C09, C15, C16 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder + gloss in A.7 #16 |
| divergence ∇·u (2.23); solenoidal fields | C10 | C10 (N54 note in C11) |
| velocity gradient G[i, j] = ∂u_i/∂x_j; divergence of a tensor (second index) | C10 | C10 (N52 note, ⚠️ callout) |
| ∂x_i/∂x_j = δ_ij; product rule for partial derivatives | C10 (Ex. 2.3 index route), C15 (D22) | C10 (gloss in A.7 #21); knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| curl ∇×u (2.24), components (2.25); irrotational fields; operator ordering | C11 | C11 (N53, N54 notes, D12) |
| curl of a solid-body rotation = 2b (Ex. 2.3) | C11 | C11 (N55 note) |
| paddle-wheel picture of curl (turns at half the curl) | C11 | C11 (idea box; proved by D15 + Ex. 2.3 in C12) |
| scipy.linalg.expm — the matrix exponential and linear trajectories x(t) = e^{Gt}x₀ | C10 animation, C12 | primer (in C12) — the C10 animation cell carries a one-line forward pointer to it |
| symmetric and antisymmetric tensors; component counts 6 and 3 | C12 | C12 (N56 note) |
| decomposition B = S + A, uniqueness | C12 | C12 (D14) |
| antisymmetric tensor of a vector (2.26), two-way map (2.27), R·x = ω × x | C12 | C12 (N57, N58 notes, D15) |
| τ_ij A_ij = 0 for symmetric τ (2.28), (2.29) | C12 | C12 (N59–N61 notes) |
| strain-rate tensor S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i) (Ex. 2.4 definition) | C12, C13 | C12 (tiny example) with the ⚠️ factor-2 note in C13 (N62) |
| eigenvalues and eigenvectors (A·b = λb, det(A − λI) = 0, `np.linalg.eigh`, `np.roots`) | C13 | primer (in C13) |
| complex conjugate and |z|² = z z̄ | C13 (D17) | primer (in C13) |
| quadratic form and the Rayleigh quotient (n·τ·n bounds); Gram–Schmidt gloss | C13 (D17) | primer (in C13) |
| square root of a negative number, i² = −1 | C13 (D17) | knowledge/primers.md: square root of a negative number (ch01 P45) — reminder |
| principal axes; diagonal form; bounds on normal and shear stress; Mohr's circle (2-D) | C13 | C13 (D17; N62 note) |
| quadratic formula for 2×2 eigenvalues | C13 from-scratch | C13 (tiny example step 2; primer eigenvalues) |
| live widgets (`live(...)`) | C13 | knowledge/primers.md: live widgets (ch01 P47) — reminder |
| Gauss' theorem (2.30); divergence theorem; outflux reading | C14 | C14 (N64 note, D25) |
| volume and surface integrals as midpoint sums; iterated integrals | C14 | primer (in C14) |
| fundamental theorem of calculus | C14 (D25) | primer (in C14) |
| definite integral; trapezoid rule (`np.trapezoid`) | C14 primer demo, C15 primer demo | knowledge/primers.md: definite integral (ch01 P27), trapezoid rule (ch01 P37) — reminders |
| tiling a volume; cancellation of interior faces | C14 (D25) | C14 (D25 steps 7–8; `divergence_theorem_tiled`) |
| sphere benchmark 8π/3 | C14 | C14 (N64 note; `reference/ch02/benchmarks.json`) |
| generalised derivative (2.31); integral definitions of divergence (2.32) and curl (2.33) | C15 | C15 (N66, N67 notes, D21) |
| mean-value theorem for integrals | C15 (D21) | primer (in C15) |
| first-order Taylor expansion (face values) | C15 (D22), C16 (D26) | knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder in D22's tools |
| midpoint rule for a face integral (exact for linear variation) | C15 (D22) | primer volume and surface integrals as midpoint sums (in C14) — reminder in D22 step 3's *why* |
| Ex. 2.5 face bookkeeping (EADH/FBCG …) | C15 | C15 (N68 note, D22) |
| coordinate-free remark (N72) | C15 | C15 (N72 note) |
| Stokes' theorem (2.34); circulation; normal curl as circulation per area (2.35) | C16 | C16 (N71 note, D26) |
| orientation rule t = n_c × n; n_c into the surface | C16 | C16 (N69 note, ⚠️ callout) |
| line integral of a vector field around a loop; parametrising a circle and a rectangle | C16 (D26) | primer (in C16) |
| line integral along a path (∫p dv) | C16 primer | knowledge/primers.md: line integral along a path (ch01 P35) — reminder inside the primer |
| Ex. 2.6 rectangle bookkeeping with the u_y correction | C16 | C16 (N73 note) |
| Green's theorem as the planar Stokes | C16 (D26) | C16 (D26 itself is the proof; named in D26's *what it means*) |
| singular fields and the failure of Stokes' hypothesis (irrotational vortex) | C16, E5 | C16 (A.11 #10 explain item 7, #16) |
| comma notation (2.36) | §2.14 | C10 (N74 note in its own section) |
| sympy: symbols, Matrix, diff, simplify, expand, eigenvects, cross | C01, C07, C10, C13 (D17 check) | knowledge/primers.md: sympy (ch01 P40), sympy Matrix and nullspace (ch01 P61) — reminders |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| tuple unpacking (`lam, B = …`) | C02 | knowledge/primers.md: tuple unpacking (ch01 P14) — reminder |
| slider_figure | C05 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| animate and show_animation | C02 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| show_viz | C03 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| np.linspace | C05 | knowledge/primers.md: np.linspace and np.logspace (ch01 P06) — reminder |
| matplotlib figures | C01 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| np.pad (embedding a 2×2 A in 3×3) | C12, E3 code | C12 (code comment: "pad with a zero row and column") |
| `np.stack` / `np.c_` (stacking components on axis 0; columns) | C09, C13 | C09 (code comment in the field cells: "components stacked on axis 0") |
| `np.count_nonzero`, `np.array_equal` | C08 | C08 (code comments) |
| `try/except ValueError` (principal_axes refusing non-symmetric input) | C13 | C13 (code comment: "the library raises on purpose") |
| climate hooks: Coriolis 2Ω × u, geostrophy, ∇p, planetary vorticity | front matter, C08, C09, C11 | front matter table (A.0 #5) and the block's plain-words paragraphs; taught in Ch. 4/13 |

Ledger rows: 123. Primers planned (`nb.primer`, new in ch02): 25 — np.einsum index strings · matrix multiplication,
transpose and identity · orthonormal basis, projection and completeness · cosines of angles between unit vectors ·
np.linalg.norm and np.linalg.qr · plotly 3-D arrows, lines and meshes · limits and orders of smallness · vector area of
a closed surface · np.arctan2 · Vieta's formulas · permutations, cyclic order and parity · numpy arrays with three axes
and np.transpose · right-hand rule and orientation · level sets and the directional derivative · np.meshgrid and the
project grid layout · numpy broadcasting · plt.contour, plt.quiver and plt.streamplot · scipy.linalg.expm · eigenvalues
and eigenvectors · complex conjugate · quadratic form and the Rayleigh quotient · volume and surface integrals as
midpoint sums · fundamental theorem of calculus · mean-value theorem for integrals · line integral of a vector field
around a loop. Reminders of ch01 primers: 20 (P01, P04, P05, P06, P09, P10, P13, P14, P15, P16, P17, P18, P21/P22, P23,
P25, P26, P27, P29, P35, P37, P38, P40, P45, P47, P49, P53, P61). Glosses: 9.

exist; the ledger line reads R01).

---

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live` and `set` are the explainer's). Every step is
one move; *why* names the rule and says why we make the move (≤ 35 words, ≥ 6). The book's own moves were read on the
rendered pages (p069–p070 for D01/D03, p075 for D05, p073 for D06, p078 for D09, p081 for D12, p082 for D14/D15,
p083–p084 for D17/D18, p085–p087 for D21/D22/D25, p088–p089 for D26); the gaps listed in `analysis/ch02.md` §2b are
filled. Colours: old frame teal, new frame orange, normal stress blue, shear rose, symmetric teal, antisymmetric orange,
flux in blue / out orange, circulation purple. (No line in this part starts with a table bar.)

### D01 · The transformation rule for components, Eq. (2.5) — ★, 6 steps, in C03 (notebook · `rotation_of_axes`)
- **Goal.** Find how the three numbers that describe a vector change when the axes are rotated — the rule that will
  *define* what a vector is.
- **Start.** $\mathbf x = x_i\mathbf e_i = x'_j\mathbf e'_j$ — *in words:* the same arrow written once in the old basis
  (2.1) and once in the new basis (2.3); only the numbers differ.
- **Plan.** (1) Dot both spellings with one new unit vector. (2) On the primed side orthonormality kills every term but
  one. (3) Name the cosines that remain: that is C. (4) Do it for a general new axis, not only e'_1.
- **Tools.** Dot product (2.2) (C01, N07) · bilinearity of the dot product (gloss in step 2) · orthonormality
  $\mathbf e'_i\cdot\mathbf e'_j = \delta_{ij}$ (primer orthonormal basis, C02) · Kronecker substitution (2.17) (C01,
  N41) · the direction cosines $C_{ij} = \mathbf e_i\cdot\mathbf e'_j$ (C02).
- **Assumptions.** Both bases orthonormal (step 3) · same origin, pure rotation — no translation (start line: the same
  arrow has the same tail).
- **Steps.**
  1. *did:* Write the arrow in both bases · *tex:* $x_i\,\mathbf e_i = x'_j\,\mathbf e'_j$ · *why:* Equations (2.1) and
     (2.3) describe one and the same arrow, so their right sides are equal. We put them side by side to compare
     coefficients. · *plain:* Two spellings, one vector.
  2. *did:* Dot both sides with $\mathbf e'_k$ · *tex:* $x_i\,(\mathbf e_i\cdot\mathbf e'_k) = x'_j\,(\mathbf e'_j\cdot
     \mathbf e'_k)$ · *why:* The dot product is bilinear: it distributes over the sums and lets the scalars $x_i$, $x'_j$
     pull out. We dot with a *new* unit vector to isolate one new component. · *plain:* Project both spellings onto
     the new axis k.
  3. *did:* Use orthonormality of the new basis · *tex:* $x_i\,(\mathbf e_i\cdot\mathbf e'_k) = x'_j\,\delta_{jk}$ ·
     *why:* The new unit vectors are perpendicular and of length 1, so $\mathbf e'_j\cdot\mathbf e'_k$ is 1 for j = k and
     0 otherwise — that is δ. This is the move the book leaves silent. · *plain:* On the primed side only the k-th
     term can survive.
  4. *did:* Let δ substitute its summed index · *tex:* $x_i\,(\mathbf e_i\cdot\mathbf e'_k) = x'_k$ · *why:*
     Kronecker substitution (2.17): $x'_j\delta_{jk} = x'_k$, because the sum over j has one nonzero term. We now
     have one new component alone on one side. · *plain:* The new component k is a weighted sum of the old ones.
  5. *did:* Name the weights as direction cosines · *tex:* $x'_k = x_i\,C_{ik},\qquad C_{ik} \equiv \mathbf e_i\cdot
     \mathbf e'_k$ · *why:* Each weight is the cosine of the angle between old axis i and new axis k (unit vectors);
     giving them a name turns the projection into a matrix. For k = 1 this line is (2.4). · *plain:* The weights are
     the nine cosines between old and new axes.
  6. *did:* Rename the free index k → j · *tex:* $x'_j = x_i\,C_{ij}$ · *why:* A free index may be renamed on both sides
     at once (N12); k = 1, 2, 3 are three equations and this one line holds them all — no "similarly" needed. · *plain:*
     The new components are the old ones multiplied into C, first index summed.
- **Result.** $x'_j = x_iC_{ij}$ (2.5), i.e. $\mathbf x' = \mathbf C^{\rm T}\mathbf x$ — *in words:* each new component is
  the old components weighted by the cosines in column j of C.
- **Check.** Units: both sides carry the unit of x (C is dimensionless) ✓. Limit θ = 0: C = I, x' = x ✓. Number: x =
  (1, 2), θ = 30°: $x'_1 = 1(0.866) + 2(0.5) = 1.866$, $x'_2 = 1(-0.5) + 2(0.866) = 1.232$; lengths 2.236 both ✓
  (`transform_vector`).
- **What it means.** The rule contains no property of x except its components: **anything** with three components that
  obey it is a Cartesian vector (2.8); anything that does not is a list of numbers. It fails for triples like $(x_1^2,
  x_2^2, x_3^2)$ and for components in a non-orthonormal basis (then δ is replaced by a metric — Appendix B).
- **Traps.** Writing $C_{ji}$ for $C_{ij}$ (row = old axis, column = new axis). Forgetting step 3, which is why only one
  primed term survives. Doing j = 1 and saying "similarly" — one indexed line does all three.

### D02 · Why C is orthogonal: C Cᵀ = Cᵀ C = I and det C = +1 — ★★, 8 steps, in C02 (notebook · `rotation_of_axes`)
- **Goal.** Show that the direction-cosine matrix undoes itself when transposed, so that going back to the old
  components is free (D03), lengths are preserved, and contracted index pairs of C become δ (D18). The book leaves this
  to Exercise 2.8.
- **Start.** $C_{ij} = \mathbf e_i\cdot\mathbf e'_j$ — *in words:* the definition of the direction cosines (C02).
- **Plan.** (1) Expand an old unit vector in the new basis (completeness). (2) Dot with another old unit vector: δ on
  one side, a product of two C's on the other. (3) Repeat with the roles of the bases swapped. (4) Take determinants
  for the sign.
- **Tools.** Completeness of an orthonormal basis, $\mathbf v = (\mathbf v\cdot\mathbf e'_j)\mathbf e'_j$ (primer, C02) ·
  orthonormality of both bases (same primer) · bilinearity of the dot product (D01 step 2) · matrix product as an
  index sum (2.9) (C01, N21) · det(AB) = det A det B (P53 reminder) · continuity argument for the sign (gloss in step
  8).
- **Assumptions.** Both bases orthonormal and complete in 3-D (steps 2, 4, 5) · "rotation" means reachable from the
  identity by turning (step 8).
- **Steps.**
  1. *did:* Write the product we want to evaluate · *tex:* $C_{ij}C_{kj} = (\mathbf e_i\cdot\mathbf e'_j)(\mathbf e_k\cdot
     \mathbf e'_j)$ · *why:* Just the definition of C inserted twice; j is summed. We choose this product because it is
     $(\mathbf C\mathbf C^{\rm T})_{ik}$. · *plain:* Row i of C dotted with row k of C.
  2. *did:* Expand $\mathbf e_i$ in the new basis · *tex:* $\mathbf e_i = (\mathbf e_i\cdot\mathbf e'_j)\,\mathbf e'_j =
     C_{ij}\,\mathbf e'_j$ · *why:* Completeness: any vector is the sum of its projections on an orthonormal basis. We
     expand an *old* axis in *new* axes because that produces exactly a row of C. · *plain:* An old axis is a sum of
     new axes weighted by its row of C.
  3. *did:* Dot with $\mathbf e_k$ · *tex:* $\mathbf e_i\cdot\mathbf e_k = C_{ij}\,(\mathbf e'_j\cdot\mathbf e_k) = C_{ij}
     C_{kj}$ · *why:* Bilinearity pulls $C_{ij}$ out; $\mathbf e'_j\cdot\mathbf e_k = \mathbf e_k\cdot\mathbf e'_j = C_{kj}$
     by the definition of C. The right side is now the product of step 1. · *plain:* The same product equals a dot
     product of two old axes.
  4. *did:* Use orthonormality of the old basis · *tex:* $C_{ij}C_{kj} = \delta_{ik}\qquad(\mathbf C\mathbf C^{\rm T} =
     \mathbf I)$ · *why:* $\mathbf e_i\cdot\mathbf e_k = \delta_{ik}$ for perpendicular unit vectors. This is the first
     orthogonality relation; in matrix form the summed second indices make $\mathbf C\mathbf C^{\rm T}$. · *plain:* The
     rows of C are perpendicular unit vectors.
  5. *did:* Swap the roles: expand a new axis in the old basis · *tex:* $\mathbf e'_i = (\mathbf e'_i\cdot\mathbf e_j)\,
     \mathbf e_j = C_{ji}\,\mathbf e_j$ · *why:* Completeness again, now of the old basis. The second relation does not
     follow from the first without this second argument — a common gap. · *plain:* A new axis is a sum of old axes
     weighted by its column of C.
  6. *did:* Dot with $\mathbf e'_k$ and use orthonormality of the new basis · *tex:* $C_{ji}C_{jk} = \delta_{ik}\qquad(
     \mathbf C^{\rm T}\mathbf C = \mathbf I)$ · *why:* Bilinearity and $\mathbf e_j\cdot\mathbf e'_k = C_{jk}$ on the right;
     $\mathbf e'_i\cdot\mathbf e'_k = \delta_{ik}$ on the left. Now the *first* indices are summed. · *plain:* The
     columns of C are perpendicular unit vectors too.
  7. *did:* Take the determinant of step 6 · *tex:* $(\det\mathbf C)^2 = 1\ \Rightarrow\ \det\mathbf C = \pm1$ · *why:*
     det(AB) = det A · det B and det Cᵀ = det C (P53), while det I = 1. We take determinants because the sign tells a
     rotation from a mirror. · *plain:* C cannot stretch volumes; it may or may not flip them.
  8. *did:* Fix the sign for a rotation · *tex:* $\det\mathbf C = +1$ · *why:* Turning the axes gradually from θ = 0
     changes det C continuously; it starts at det I = +1 and can never pass through 0 (it is ±1), so it stays +1. A
     mirror (−1) cannot be reached by turning. · *plain:* A genuine rotation keeps the handedness of the axes.
- **Result.** $C_{ij}C_{kj} = C_{ji}C_{jk} = \delta_{ik}$, $\det\mathbf C = +1$ — *in words:* Cᵀ is the inverse of C, and C
  is a proper rotation.
- **Check.** θ = 30° in 2-D: row 1 · row 1 = 0.75 + 0.25 = 1, row 1 · row 2 = 0.866(0.5) + (−0.5)(0.866) = 0 ✓; det =
  0.866² + 0.5² = 1 ✓ (`orthogonality_residual` ≈ 2e-16, `is_proper_rotation` True). Limit: C = I trivially ✓. 200
  random rotations from `random_rotation` all pass.
- **What it means.** A rotation is a relabelling of directions that loses nothing: lengths, angles and volumes are
  kept; the inverse costs a transpose (D03). Every time two C's share a summed index in a formula, they collapse to δ
  — the engine of D06 and D18. It fails for skewed or non-unit axes (then $\mathbf C^{-1} \neq \mathbf C^{\rm T}$).
- **Traps.** Proving CᵀC = I and assuming CCᵀ = I "obviously" — for square matrices it is true, but the basis argument
  needs steps 5–6 to say why. A det of −1 passes the orthogonality test and is still not a rotation.

### D03 · The inverse transformation, Eq. (2.7) — ★, 4 steps, in C02 (notebook)
- **Goal.** Get the old components back from the new ones — the book says "it can be shown (Exercise 2.2)".
- **Start.** $x'_j = x_iC_{ij}$ (2.5) — *in words:* the forward rule of D01 (derived in C03; here we only need that it
  holds).
- **Plan.** (1) Multiply by a C with the free index and sum. (2) Two C's with a shared summed index become δ (D02).
  (3) Let δ substitute. (4) Rename.
- **Tools.** Dummy renaming (N12) · orthogonality $C_{ij}C_{kj} = \delta_{ik}$ (D02, C02) · Kronecker substitution
  (2.17) (C01).
- **Assumptions.** Same as D01/D02: orthonormal bases (enters through D02 in step 2).
- **Steps.**
  1. *did:* Multiply both sides by $C_{kj}$ and sum on j · *tex:* $x'_j\,C_{kj} = x_i\,C_{ij}C_{kj}$ · *why:* Multiplying
     an identity by the same factor on both sides keeps it true; summing on the now-repeated j too. We choose $C_{kj}$
     so that the C-pair on the right can collapse. · *plain:* Contract the new components with a row of C.
  2. *did:* Collapse the C-pair with D02 · *tex:* $x'_j\,C_{kj} = x_i\,\delta_{ik}$ · *why:* Orthogonality $C_{ij}C_{kj} =
     \delta_{ik}$ (second indices summed) from D02 step 4. This is where "C is orthogonal" does its work. · *plain:*
     The two C's cancel each other.
  3. *did:* Let δ substitute its summed index · *tex:* $x_k = x'_j\,C_{kj}$ · *why:* Kronecker substitution (2.17):
     $x_i\delta_{ik} = x_k$. We now have an old component alone. · *plain:* An old component is a weighted sum of the
     new ones.
  4. *did:* Rename the letters to the book's (k → j, j → i) · *tex:* $x_j = x'_i\,C_{ji}$ · *why:* Renaming a free index on
     both sides and a dummy anywhere changes nothing (N12); this is (2.7). Note the summed index of C is now its
     *second* one. · *plain:* Going back uses C with its indices the other way round: x = C x'.
- **Result.** $x_j = x'_iC_{ji}$ (2.7), $\mathbf x = \mathbf C\,\mathbf x'$ — *in words:* the inverse of $\mathbf x' =
  \mathbf C^{\rm T}\mathbf x$ is multiplication by C itself, because $\mathbf C^{-1} = \mathbf C^{\rm T}$.
- **Check.** Round trip with the D01 numbers: $x_1 = 1.866(0.866) + 1.232(-0.5) = 1.616 - 0.616 = 1.000$ ✓, $x_2 =
  1.866(0.5) + 1.232(0.866) = 0.933 + 1.067 = 2.000$ ✓ (`inverse_transform_vector`). Limit θ = 0: identity ✓.
- **What it means.** Rotating back is as cheap as rotating forward — no matrix inversion. In (2.5) the first index of C
  is summed, in (2.7) the second: the two index placements are the whole difference between C and Cᵀ.
- **Traps.** Using (2.5) with the indices swapped and calling it the inverse without D02 — it happens to be right only
  because C is orthogonal.

### D05 · Cauchy's traction formula from a shrinking tetrahedron, Eq. (2.15) — ★★, 9 steps, in C05 (notebook · `cauchy_traction_principal_axes`)
- **Goal.** Find the force per unit area on a plane with *any* unit normal n from the nine stresses on the three
  coordinate planes.
- **Start.** $\sum\mathbf F = m\,\mathbf a$ on a tiny tetrahedron — *in words:* Newton's second law for the fluid inside
  a tetrahedron with three faces on the coordinate planes and one slanted face with outward normal n.
- **Plan.** (1) List the forces: the slanted face, the three coordinate faces, weight and inertia. (2) Write the
  coordinate faces' tractions with C04's sign rule. (3) Relate the face areas to n. (4) Shrink: volume terms vanish
  faster than face terms.
- **Tools.** Free-body diagram and Newton's second law (P09 reminder) · the stress sign convention on the −e_j faces
  (C04) · vector area of a closed surface, $dA_j = n_j\,dA$ (primer, C05) · orders of smallness h² vs h³ (primer, C05)
  · limits.
- **Assumptions.** Stresses continuous at the point (step 3: the face tractions are the point's τ) · body force and
  acceleration finite (step 4) · the tetrahedron shrinks to the point (steps 8–9).
- **Steps.**
  1. *did:* Write the force balance on the tetrahedron · *tex:* $\mathbf F_{\rm slant} + \mathbf F_1 + \mathbf F_2 +
     \mathbf F_3 + \mathbf F_{\rm body} = \rho\,dV\,\mathbf a$ · *why:* Newton's second law for the mass ρ dV inside; every
     force on the boundary and in the volume is listed (free-body diagram, P09). · *plain:* Four face forces plus the
     body force equal mass times acceleration.
  2. *did:* Write the slanted face's force · *tex:* $\mathbf F_{\rm slant} = \mathbf f\,dA$ · *why:* By definition $\mathbf
     f$ is the force per unit area on the plane with normal n, and dA is that face's area. This is the unknown we
     want. · *plain:* The slanted face carries the traction we are after.
  3. *did:* Write the coordinate faces' forces with the sign rule · *tex:* $(F_j)_i = -\,\tau_{ji}\,dA_j\quad(j = 1, 2,
     3\text{, no sum})$ · *why:* The face perpendicular to axis j has *outward* normal $-\mathbf e_j$; by C04's rule the
     traction on a − face is the reversed row j, so its i-component is $-\tau_{ji}$. · *plain:* Each coordinate face
     pushes with minus its row of τ.
  4. *did:* Write body force and inertia together · *tex:* $\mathbf F_{\rm body} - \rho\,dV\,\mathbf a = \rho(\mathbf g -
     \mathbf a)\,dV$ · *why:* Both are proportional to the mass, hence to the volume dV; we group them because they
     will share the same fate when the element shrinks. · *plain:* Weight and inertia scale with the volume.
  5. *did:* Assemble the i-component · *tex:* $f_i\,dA - \tau_{1i}dA_1 - \tau_{2i}dA_2 - \tau_{3i}dA_3 + \rho(g_i - a_i)
     dV = 0$ · *why:* Steps 2–4 substituted into step 1, component i. The book's first line is this one without the
     volume term (and for i = 1). · *plain:* Slanted face minus the three coordinate faces plus volume terms is zero.
  6. *did:* Relate the face areas to the normal · *tex:* $dA_j = n_j\,dA$ · *why:* The vector area of a closed surface
     vanishes: $\mathbf n\,dA - \mathbf e_1dA_1 - \mathbf e_2dA_2 - \mathbf e_3dA_3 = 0$ (primer); reading component j
     gives this. The coordinate face j is the shadow of the slanted face. · *plain:* Each coordinate face is the
     slanted face times a cosine.
  7. *did:* Substitute the areas and write the sum with an index · *tex:* $f_i\,dA - \tau_{ji}\,n_j\,dA + \rho(g_i -
     a_i)\,dV = 0$ · *why:* Step 6 into step 5; the three face terms share the pattern $\tau_{ji}n_j dA$ with j summed
     (summation convention, C01). · *plain:* The three faces together give τ contracted with n on its first index.
  8. *did:* Divide by dA and compare sizes · *tex:* $f_i - \tau_{ji}n_j + \rho(g_i - a_i)\,\frac{dV}{dA} = 0,\qquad
     \frac{dV}{dA}\propto h$ · *why:* dA ≠ 0. For a tetrahedron of size h, dA ∝ h² and dV ∝ h³, so dV/dA ∝ h (primer
     orders of smallness); ρ, g, a stay finite. · *plain:* The volume term is smaller than the face terms by the
     element's size.
  9. *did:* Let the element shrink to the point · *tex:* $f_i = \tau_{ji}\,n_j$ · *why:* As h → 0 the last term of step 8
     vanishes while the others do not depend on h (stresses continuous at the point). This is (2.15). · *plain:* The
     traction on any plane is τ contracted with the normal on its first index.
- **Result.** $f_i = \tau_{ji}n_j$, $\mathbf f = \mathbf n\cdot\boldsymbol\tau$ (2.15) — *in words:* the nine numbers on the
  coordinate planes determine the force per area on every plane through the point.
- **Check.** Units: Pa on both sides ✓. Limit n = e₃: $f_i = \tau_{3i}$, the traction on the coordinate face — the
  formula returns its inputs ✓. Number (Ex. 2.2, a = 1 Pa, φ = 30°): $f = (0\cdot0.866 + 1\cdot0.5,\ 1\cdot0.866 +
  0\cdot0.5) = (0.5, 0.866)$ Pa ✓ (`traction`, `example_2_2`). Pure pressure τ = −pδ: f = −p n for every n (Ch. 1's
  isotropic pressure) ✓.
- **What it means.** Stress is nine numbers, not a vector, precisely because the force per area depends on the cut;
  yet those nine suffice for every cut. Every wall force in the book is $\oint\mathbf n\cdot\boldsymbol\tau\,dA$. The
  formula contracts the *first* index; $\boldsymbol\tau\cdot\mathbf n$ equals it only when τ is symmetric (Ch. 4 proves
  that for the stress). It fails where stresses are discontinuous (a shock, an interface with surface tension — then
  the jump conditions of Ch. 4 apply).
- **Traps.** Forgetting the minus signs of step 3 (outward normal −e_j). Dropping the volume terms without saying why
  (they are O(h³) against O(h²)). Writing $\tau_{ij}n_j$ — that needs symmetry.

### D06 · The tensor transformation rule from Cauchy's formula, Eq. (2.12) — ★★, 8 steps, in C06 (notebook · `rotation_of_axes`)
- **Goal.** Show how the nine stress components must change under a rotation of axes, given only that force per area
  and the normal are vectors. The book states (2.12) and cites Sommerfeld's tetrahedron; we have the tetrahedron result
  (2.15) already, so we derive (2.12) from it.
- **Start.** $f_i = \tau_{ji}n_j$ in the old frame and $f'_n = \tau'_{mn}n'_m$ in the new frame — *in words:* Cauchy's
  formula holds for every observer, each with their own components of τ.
- **Plan.** (1) Transform f as a vector. (2) Insert Cauchy's formula in the old frame. (3) Express the old n through
  the new n'. (4) Compare with Cauchy's formula in the new frame, for every n'.
- **Tools.** f and n are vectors, (2.8) (C03) · Cauchy's formula (2.15) (C05) · the inverse transformation (2.7) (C02,
  D03) · dummy renaming (N12) · "true for every n' ⇒ the coefficients agree" (gloss in step 7).
- **Assumptions.** τ is *defined* in each frame by (2.15) (start and step 5) · C orthogonal (enters through (2.7) in
  step 3).
- **Steps.**
  1. *did:* Transform the traction as a vector · *tex:* $f'_n = f_i\,C_{in}$ · *why:* f is a physical force per area,
     hence a vector, so its components obey (2.8) with the first index of C summed. We start from f because both
     frames must agree on it. · *plain:* The new observer's traction components come from the old ones by C.
  2. *did:* Insert Cauchy's formula in the old frame · *tex:* $f'_n = \tau_{ji}\,n_j\,C_{in}$ · *why:* (2.15) gives $f_i$
     in terms of the old τ and the old n; substituting expresses the new traction through old quantities. · *plain:*
     The new traction, written with the old stresses and the old normal.
  3. *did:* Write the old normal through the new one · *tex:* $n_j = n'_m\,C_{jm}$ · *why:* n is a unit vector, so
     the inverse rule (2.7) applies (second index of C summed). We need n' because the new observer measures n'. ·
     *plain:* The old components of the normal from the new ones.
  4. *did:* Substitute the normal · *tex:* $f'_n = C_{jm}\,C_{in}\,\tau_{ji}\,n'_m$ · *why:* Step 3 into step 2; the
     factors are reordered (scalars commute) so that both C's sit in front. Everything on the right is now old τ, two
     C's and the new normal. · *plain:* The new traction equals two C's times the old τ times the new normal.
  5. *did:* Write Cauchy's formula in the new frame · *tex:* $f'_n = \tau'_{mn}\,n'_m$ · *why:* The new observer sees the
     same physics, so (2.15) holds with primed components everywhere — this *defines* the new observer's τ'. ·
     *plain:* The new observer's own stress tensor gives the same traction.
  6. *did:* Equate steps 4 and 5 · *tex:* $\big(\tau'_{mn} - C_{jm}C_{in}\tau_{ji}\big)\,n'_m = 0$ · *why:* Both are the
     same $f'_n$; subtracting and factoring the common $n'_m$ (a sum on m) collects the difference. The bracket does
     not depend on n'. · *plain:* A certain matrix, applied to any unit normal, gives zero.
  7. *did:* Conclude the bracket vanishes · *tex:* $\tau'_{mn} = C_{jm}\,C_{in}\,\tau_{ji}$ · *why:* Step 6 holds for every
     unit n'; choose n' = e'_1, e'_2, e'_3 in turn and each picks out one column of the bracket, which must be zero
     (compare coefficients). · *plain:* The new stress components are fixed by the old ones and C.
  8. *did:* Rename the dummies (j → i, i → j) · *tex:* $\tau'_{mn} = C_{im}\,C_{jn}\,\tau_{ij}\qquad(\boldsymbol\tau' =
     \mathbf C^{\rm T}\boldsymbol\tau\,\mathbf C)$ · *why:* Swapping the names of two summed letters changes nothing
     (N12); this is the book's (2.12). In matrix form $C_{im} = (\mathbf C^{\rm T})_{mi}$, so Cᵀ stands on the left and C
     on the right. · *plain:* One C per index — the vector rule applied twice.
- **Result.** $\tau'_{mn} = C_{im}C_{jn}\tau_{ij}$, $\boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau\,\mathbf C$ (2.12) —
  *in words:* a second-order tensor transforms with one direction-cosine matrix per index; whatever obeys this is a
  tensor.
- **Check.** Units: Pa ✓. Limit C = I: τ' = τ ✓. Number (Ex. 2.2, a = 1, θ = 30°): $\tau'_{11} = 2(0.866)(0.5) = 0.866$,
  $\tau'_{12} = 0.75 - 0.25 = 0.5$, $\tau'_{22} = -0.866$ ✓ (`transform_tensor` = `C.T @ tau @ C`). Trace 0 before and
  after (D18).
- **What it means.** The rule contains nothing about stress: any nine numbers that obey it are the components of a
  physical object, and any that do not are not (N27). Ch. 3's strain rate, Ch. 4's stress and Ch. 12's Reynolds stress
  all pass. The book's order (2.12 stated first, 2.15 derived later) is reversed here on purpose.
- **Traps.** Writing $\mathbf C\boldsymbol\tau\mathbf C^{\rm T}$ — that is the *active* convention; ours is passive. Forgetting
  the final rename and comparing $C_{jm}C_{in}\tau_{ji}$ with the book. Assuming symmetry of τ anywhere — not needed.

### D09 · The epsilon–delta relation, its contractions and the triple product, Eq. (2.19) — ★★, 11 steps, in C08 (notebook · backup `index_machine`)
- **Goal.** Prove the identity that turns a product of two alternating tensors into Kronecker deltas — the tool that
  closes every vector identity (curl of a curl, a × (b × c), the vorticity algebra of Ch. 5). The book says "verify by
  choosing some values".
- **Start.** $L_{ijlm} \equiv \varepsilon_{ijk}\,\varepsilon_{klm}$ (sum on k) — *in words:* the left side of (2.19) as a
  four-index object we will evaluate case by case.
- **Plan.** (1) Both sides vanish when i = j or l = m. (2) For i ≠ j only one k survives, and both ε's are nonzero only
  if {l, m} = {i, j}. (3) Evaluate the two surviving cases: +1 and −1. (4) Contract the result once, twice; apply it to
  a × (b × c).
- **Tools.** Definition of ε (2.18) and its index moves (C08, N43) · antisymmetry: swapping two indices flips the sign ·
  Kronecker substitution and δ_ii = 3 (C01, N41) · the index form of the cross product (2.21) (C08, N49) · permutation
  parity (primer, C08).
- **Assumptions.** Three dimensions (the count "exactly one k outside {i, j}" is 3-D; step 4).
- **Steps.**
  1. *did:* Fix the free indices and expand the sum · *tex:* $L_{ijlm} = \varepsilon_{ij1}\varepsilon_{1lm} + \varepsilon_
     {ij2}\varepsilon_{2lm} + \varepsilon_{ij3}\varepsilon_{3lm}$ · *why:* The summation convention on the repeated k
     (C01). We write the three terms out so that "which k survives" can be argued. · *plain:* Three products, one per
     value of the summed index.
  2. *did:* Dispose of i = j · *tex:* $i = j:\quad L_{iilm} = 0 = \delta_{il}\delta_{im} - \delta_{im}\delta_{il}$ · *why:*
     ε with two equal indices is zero (2.18), so every term of step 1 vanishes; the right side of (2.19) is a difference
     of two equal products. Both sides are antisymmetric in (i, j). · *plain:* Equal first indices: zero on both sides.
  3. *did:* Dispose of l = m the same way · *tex:* $l = m:\quad L_{ijll} = 0 = \delta_{il}\delta_{jl} - \delta_{il}\delta_
     {jl}$ · *why:* $\varepsilon_{kll} = 0$ kills the left side; the right side is again a difference of identical
     products. From now on i ≠ j and l ≠ m. · *plain:* Equal last indices: zero on both sides.
  4. *did:* Identify the only surviving k · *tex:* $i \neq j:\quad \varepsilon_{ijk} \neq 0 \iff k = k^*,\ \{i, j, k^*\}
     = \{1, 2, 3\}$ · *why:* In 3-D, with i and j distinct, exactly one value of k differs from both; for it ε is ±1, for
     the other two ε is 0 (2.18). The sum in step 1 has one term. · *plain:* Only the third index completing the set
     contributes.
  5. *did:* Require the second ε to be nonzero too · *tex:* $\varepsilon_{k^*lm} \neq 0 \iff \{l, m\} = \{i, j\}$ · *why:*
     l and m must be distinct from $k^*$ and from each other, and the only two values left are i and j. If {l, m} ≠
     {i, j}, L = 0 — and then the right side is 0 as well, since l or m equals $k^*$ and every δ there has a zero
     factor. · *plain:* Both ε's are alive only when the last pair is the first pair, in some order.
  6. *did:* Evaluate the case l = i, m = j · *tex:* $L_{ijij} = \varepsilon_{ijk^*}\,\varepsilon_{k^*ij} = \varepsilon_{ijk^*}
     ^2 = 1$ · *why:* Moving $k^*$ two places ($\varepsilon_{k^*ij} = \varepsilon_{ijk^*}$, N43) makes the two factors
     equal; a nonzero ε squared is 1. Right side: $\delta_{ii}\delta_{jj} - \delta_{ij}\delta_{ji} = 1\cdot1 - 0 = 1$ (here
     $\delta_{ii}$ with i fixed is 1, not 3). · *plain:* Same pair, same order: plus one on both sides.
  7. *did:* Evaluate the case l = j, m = i · *tex:* $L_{ijji} = \varepsilon_{ijk^*}\,\varepsilon_{k^*ji} = -\varepsilon_{ijk^*}
     ^2 = -1$ · *why:* $\varepsilon_{k^*ji} = -\varepsilon_{k^*ij}$ (one swap flips the sign), then as in step 6. Right
     side: $\delta_{ij}\delta_{ji} - \delta_{ii}\delta_{jj} = 0 - 1 = -1$. All 81 cases are now covered: (2.19) holds. ·
     *plain:* Same pair, swapped order: minus one on both sides.
  8. *did:* Contract (2.19) on j = m · *tex:* $\varepsilon_{ijk}\,\varepsilon_{klj} = \delta_{il}\delta_{jj} - \delta_{ij}
     \delta_{jl} = 3\delta_{il} - \delta_{il} = 2\delta_{il}$ · *why:* Set m = j and sum: now $\delta_{jj} = 3$ (a summed
     pair, C01 N41) and $\delta_{ij}\delta_{jl} = \delta_{il}$ by substitution. Two-place moves rewrite this as the book's
     $\varepsilon_{pqi}\varepsilon_{pqj} = 2\delta_{ij}$. · *plain:* One contraction leaves twice a delta.
  9. *did:* Contract again on i = l · *tex:* $\varepsilon_{pqr}\,\varepsilon_{pqr} = 2\delta_{rr} = 6$ · *why:* Set the
     remaining free pair equal and sum; $\delta_{rr} = 3$. This counts the nonzero entries of ε: six of them, each
     squared to 1. · *plain:* Fully contracted, ε with itself gives 6.
  10. *did:* Write a × (b × c) with two ε's · *tex:* $[\mathbf a\times(\mathbf b\times\mathbf c)]_m = \varepsilon_{mpq}\,a_p
      \,\varepsilon_{qij}\,b_i c_j$ · *why:* (2.21) twice: $(\mathbf b\times\mathbf c)_q = \varepsilon_{ijq}b_ic_j =
      \varepsilon_{qij}b_ic_j$ and $(\mathbf a\times\mathbf d)_m = \varepsilon_{pqm}a_pd_q = \varepsilon_{mpq}a_pd_q$ (two-
      place moves). The shared summed index q is adjacent in both ε's. · *plain:* A double cross product is two ε's
      contracted on one index.
  11. *did:* Apply (2.19) and let the deltas substitute · *tex:* $[\mathbf a\times(\mathbf b\times\mathbf c)]_m = (\mathbf a
      \cdot\mathbf c)\,b_m - (\mathbf a\cdot\mathbf b)\,c_m$ · *why:* $\varepsilon_{mpq}\varepsilon_{qij} = \delta_{mi}\delta_
      {pj} - \delta_{mj}\delta_{pi}$ by (2.19) with (i, j, k, l, m) → (m, p, q, i, j); then $\delta_{mi}b_i = b_m$,
      $\delta_{pj}a_pc_j = \mathbf a\cdot\mathbf c$, etc. · *plain:* The "BAC − CAB" rule falls out in one line.
- **Result.** $\varepsilon_{ijk}\varepsilon_{klm} = \delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}$ (2.19);
  $\varepsilon_{pqi}\varepsilon_{pqj} = 2\delta_{ij}$; $\varepsilon_{pqr}\varepsilon_{pqr} = 6$; $\mathbf a\times(\mathbf b
  \times\mathbf c) = (\mathbf a\cdot\mathbf c)\mathbf b - (\mathbf a\cdot\mathbf b)\mathbf c$ — *in words:* a product of two
  ε's sharing an index is a difference of two δ-pairs, "first with first, second with second, minus the crossed pair".
- **Check.** Dimensionless ✓. `epsilon_delta_residual()` = 0 over all 81 cases; `einsum('pqi,pqj')` = 2I, `einsum
  ('pqr,pqr')` = 6 ✓. Number: a = (1, 0, 0), b = (0, 1, 0), c = (1, 0, 0): b × c = (0, 0, −1), a × (b × c) = (0, 1, 0);
  right side (a·c)b − (a·b)c = 1·b − 0 = (0, 1, 0) ✓ (`triple_product`).
- **What it means.** Every identity of vector calculus that contains two crosses — ∇ × (∇ × u) = ∇(∇·u) − ∇²u, the
  Lagrange identity |u × v|² = u²v² − (u·v)², the curl of b × x in Ex. 2.3, the ω ↔ R map of D15 — is this one line
  plus index bookkeeping. It is a 3-D statement (in 2-D the alternating symbol has two indices and ε_ij ε_kl = δ_ik δ_jl
  − δ_il δ_jk).
- **Traps.** δ_ii = 3 when summed, 1 when the index is fixed (steps 6 vs 8). The sign when moving the summed k next
  to its partner. Treating the 81-case table as the proof — it is the check; the proof is the antisymmetry argument.

### D12 · The three components of the curl, Eq. (2.25) — ★, 5 steps, in C11 (notebook · `stokes_circulation_loop`)
- **Goal.** Turn the index form of the curl into the three explicit component formulas — and see why only two of the
  nine terms survive for each component.
- **Start.** $(\nabla\times\mathbf u)_i = \varepsilon_{ijk}\,\dfrac{\partial u_k}{\partial x_j}$ (2.24) — *in words:* the
  cross product (2.21) with the first vector replaced by the operator ∇ = e_j ∂/∂x_j, acting to the right.
- **Plan.** (1) Fix i = 1. (2) Keep only the (j, k) pairs where ε ≠ 0. (3) Insert ±1. (4) Relabel cyclically for i =
  2, 3.
- **Tools.** Enumeration of ε (2.18) (C08) · the cross product in index form (2.21) (C08, N49) · operator ordering: ∂
  acts on u_k (gloss in the start line) · partial derivative (P25 reminder).
- **Assumptions.** u differentiable (the partial derivatives exist; every step).
- **Steps.**
  1. *did:* Write (2.24) for i = 1 with both sums · *tex:* $(\nabla\times\mathbf u)_1 = \sum_{j=1}^{3}\sum_{k=1}^{3}
     \varepsilon_{1jk}\,\dfrac{\partial u_k}{\partial x_j}$ · *why:* The summation convention (C01) hides a double sum
     over j and k — nine terms. We write it out to see which terms are nonzero. · *plain:* Nine candidate terms for the
     first component.
  2. *did:* Drop the terms with a repeated index · *tex:* $(\nabla\times\mathbf u)_1 = \varepsilon_{123}\,\dfrac{\partial
     u_3}{\partial x_2} + \varepsilon_{132}\,\dfrac{\partial u_2}{\partial x_3}$ · *why:* ε vanishes whenever two indices
     agree (2.18); with i = 1 fixed, only (j, k) = (2, 3) and (3, 2) avoid a repeat. Seven terms die. · *plain:* Only
     the two "other" directions contribute.
  3. *did:* Insert the values of ε · *tex:* $(\nabla\times\mathbf u)_1 = \dfrac{\partial u_3}{\partial x_2} - \dfrac{\partial
     u_2}{\partial x_3}$ · *why:* 123 is cyclic, so $\varepsilon_{123} = +1$; 132 is anticyclic, so $\varepsilon_{132} = -1$
     (2.18). The minus sign comes from ε, not from the derivative. · *plain:* The first component is a difference of two
     cross-derivatives.
  4. *did:* Relabel cyclically for i = 2 · *tex:* $(\nabla\times\mathbf u)_2 = \dfrac{\partial u_1}{\partial x_3} - \dfrac{
     \partial u_3}{\partial x_1}$ · *why:* Replacing every index by its cyclic successor 1 → 2 → 3 → 1 leaves ε unchanged
     (a two-place move, N43), so step 3 maps to this line. The survivors are now (j, k) = (3, 1) and (1, 3). · *plain:*
     Second component: the same pattern one step round the cycle.
  5. *did:* Relabel once more for i = 3 and assemble · *tex:* $(\nabla\times\mathbf u)_3 = \dfrac{\partial u_2}{\partial x_1}
     - \dfrac{\partial u_1}{\partial x_2}$ · *why:* Another cyclic step; survivors (1, 2) and (2, 1). The three lines
     together are (2.25). Throughout, ∂/∂x_j acts on u_k — $u_k\,\partial_j$ would be an operator, not a number. ·
     *plain:* Third component: does u₂ grow along 1 more than u₁ grows along 2?
- **Result.** $(\nabla\times\mathbf u)_1 = \partial_2u_3 - \partial_3u_2$, $(\nabla\times\mathbf u)_2 = \partial_3u_1 -
  \partial_1u_3$, $(\nabla\times\mathbf u)_3 = \partial_1u_2 - \partial_2u_1$ (2.25) — *in words:* each component of the
  curl is the "cross-derivative" of the other two velocity components.
- **Check.** Units: (m/s)/m = 1/s ✓. Number (Ex. 2.3, u = b × x = (−x₂, x₁, 0)): $(\nabla\times\mathbf u)_3 = 1 - (-1)
  = 2 = 2b_3$, the other two 0 ✓ (`curl_components` = `curl`). Shear flow u₁ = Γx₂: $(\nabla\times\mathbf u)_3 = 0 -
  \Gamma = -\Gamma$ ✓.
- **What it means.** The curl measures the local spin: a paddle wheel turns at half of it (D15 + Ex. 2.3). It is
  nonzero for a straight shear flow and zero for the irrotational vortex away from its core — spin is about
  differences of speed across a point, not about curved paths. Requires differentiable u; at a vortex core the
  formula is undefined and the circulation (D26) takes over.
- **Traps.** Writing $u_k\partial_j$. Losing the minus from $\varepsilon_{132}$. Doing i = 1 and saying "similarly" —
  the cyclic relabelling is the reason, and it is said.

### D14 · The unique split into symmetric and antisymmetric parts — ★, 6 steps, in C12 (notebook · `strain_vs_rotation_split`)
- **Goal.** Show that any second-order tensor is the sum of a symmetric and an antisymmetric tensor in exactly one way,
  and that both parts are themselves tensors — the basis of Ch. 3's strain rate and rotation rate.
- **Start.** $B_{ij}$, any second-order tensor — *in words:* nine numbers that transform by (2.12).
- **Plan.** (1) Add and subtract half the transpose. (2) Check the symmetry of each half. (3) Suppose a second split and
  show it coincides. (4) Transform each half with (2.12).
- **Tools.** Transpose $B_{ji}$ (primer matrix multiplication, transpose and identity, C01) · the definitions of
  symmetric and antisymmetric (C12, N56) · linearity of (2.12) (C06) · dummy renaming (N12) · "both symmetric and
  antisymmetric ⇒ zero" (step 5).
- **Assumptions.** None beyond B being a tensor (step 6).
- **Steps.**
  1. *did:* Add and subtract $\tfrac12B_{ji}$ · *tex:* $B_{ij} = \tfrac12(B_{ij} + B_{ji}) + \tfrac12(B_{ij} - B_{ji})$ ·
     *why:* The right side is $\tfrac12B_{ij} + \tfrac12B_{ij}$ plus and minus the same $\tfrac12B_{ji}$ — an identity.
     We choose the transpose because swapping i and j is the operation whose behaviour we want to separate. · *plain:*
     Split B into "the part that ignores the swap" and "the part that flips".
  2. *did:* Name the first half and swap its indices · *tex:* $S_{ij} \equiv \tfrac12(B_{ij} + B_{ji}),\qquad S_{ji} = S_{ij}$ ·
     *why:* Swapping i ↔ j in the definition gives $\tfrac12(B_{ji} + B_{ij})$, the same sum in the other order. S is
     symmetric by construction. · *plain:* The first half does not notice the swap.
  3. *did:* Name the second half and swap its indices · *tex:* $A_{ij} \equiv \tfrac12(B_{ij} - B_{ji}),\qquad A_{ji} =
     -A_{ij},\quad A_{ii} = 0$ · *why:* Swapping gives $\tfrac12(B_{ji} - B_{ij}) = -A_{ij}$; for i = j the difference
     vanishes. A is antisymmetric with zero diagonal — three independent entries. · *plain:* The second half changes
     sign under the swap and has an empty diagonal.
  4. *did:* Suppose another split and subtract · *tex:* $B = \tilde S + \tilde A\ \Rightarrow\ \tilde S_{ij} - S_{ij} =
     A_{ij} - \tilde A_{ij} \equiv X_{ij}$ · *why:* If a second symmetric $\tilde S$ and antisymmetric $\tilde A$ also add
     to B, subtracting the two splits and moving terms gives one object X written two ways. This is the uniqueness
     half the book omits. · *plain:* Any other split differs from ours by one tensor X.
  5. *did:* Show X is zero · *tex:* $X_{ij} = X_{ji} = -X_{ij}\ \Rightarrow\ X_{ij} = 0$ · *why:* X is a difference of
     symmetric tensors (so symmetric) and a difference of antisymmetric ones (so antisymmetric); a number equal to its
     own negative is zero. Hence $\tilde S = S$, $\tilde A = A$. · *plain:* Only the zero tensor is both symmetric and
     antisymmetric — the split is unique.
  6. *did:* Transform S with (2.12) · *tex:* $S'_{mn} = \tfrac12(B'_{mn} + B'_{nm}) = C_{im}C_{jn}\,\tfrac12(B_{ij} + B_{ji})
     = C_{im}C_{jn}S_{ij}$ · *why:* (2.12) is linear, so it applies to each term; in the second term $B'_{nm} = C_{in}C_
     {jm}B_{ij}$, and renaming the dummies i ↔ j turns it into $C_{im}C_{jn}B_{ji}$. The same lines with a minus give
     A. · *plain:* The symmetric part of the transformed B is the transformed symmetric part — S and A are tensors.
- **Result.** $B_{ij} = S_{ij} + A_{ij}$ with $S = \tfrac12(\mathbf B + \mathbf B^{\rm T})$, $A = \tfrac12(\mathbf B -
  \mathbf B^{\rm T})$, unique, both tensors; 6 + 3 = 9 components — *in words:* every tensor is a stretch-like part plus
  a spin-like part, and every observer agrees on the split.
- **Check.** Units: those of B ✓. Limit: B symmetric ⇒ A = 0; B antisymmetric ⇒ S = 0 ✓. Number: B = [[0, 2],[0, 0]]:
  S = [[0, 1],[1, 0]], A = [[0, 1],[−1, 0]], S + A = B ✓ (`symmetric_part`, `antisymmetric_part`); `independent_components`
  6, 3, 9.
- **What it means.** For the velocity gradient $\partial u_i/\partial x_j$, S is the strain rate (how a fluid element
  deforms) and A the rotation rate (how it spins) — Ch. 3 §3.4; Ch. 4's viscous dissipation sees only S (N59–N61).
  Frame-independence (step 6) is what lets us compute the split in any convenient axes.
- **Traps.** Thinking the split depends on the axes. Forgetting the uniqueness half. Counting 3 + 3 = 6 for the
  antisymmetric part (the diagonal is zero: three numbers).

### D15 · The vector hidden in an antisymmetric tensor, Eqs. (2.26)–(2.27) — ★★, 7 steps, in C12 (notebook · `strain_vs_rotation_split`)
- **Goal.** Show that the three independent numbers of an antisymmetric tensor form a vector ω, give the two-way map,
  and show that R acting on x rotates it: R·x = ω × x. The book checks one entry; the inverse and the sign are ours.
- **Start.** $R_{ij} \equiv -\varepsilon_{ijk}\,\omega_k$ (2.27, first part) — *in words:* pack a vector into a 3×3 array
  with the alternating tensor.
- **Plan.** (1) Check R is antisymmetric and read off its entries. (2) Contract R with ε to get ω back (D09's 2δ). (3)
  Apply R to a vector and recognise a cross product.
- **Tools.** Definition and index moves of ε (C08, N43) · the 2δ contraction $\varepsilon_{ijl}\varepsilon_{ijk} =
  2\delta_{lk}$ (D09 step 8) · Kronecker substitution (C01) · the cross product in index form (2.21) (C08, N49).
- **Assumptions.** Three dimensions (ε has three indices; everywhere).
- **Steps.**
  1. *did:* Swap the indices of R · *tex:* $R_{ji} = -\varepsilon_{jik}\,\omega_k = +\varepsilon_{ijk}\,\omega_k = -R_{ij}$ ·
     *why:* Swapping two indices of ε flips its sign (one-place move, N43). So R is antisymmetric — the object we want
     to represent. · *plain:* The packed array is antisymmetric, as promised.
  2. *did:* Read off the entries · *tex:* $R_{12} = -\varepsilon_{123}\omega_3 = -\omega_3,\quad R_{13} = -\varepsilon_
     {132}\omega_2 = +\omega_2,\quad R_{23} = -\varepsilon_{231}\omega_1 = -\omega_1$ · *why:* For each (i, j) only the
     k different from both survives; ε_123 = ε_231 = +1, ε_132 = −1 (2.18). The diagonal is zero by step 1. With step
     1 for the lower triangle this is the matrix (2.26). · *plain:* Above the diagonal: −ω₃, +ω₂, −ω₁.
  3. *did:* Contract R with ε to invert · *tex:* $\varepsilon_{ijl}\,R_{ij} = -\varepsilon_{ijl}\,\varepsilon_{ijk}\,\omega_k$
     · *why:* Multiply the definition by $\varepsilon_{ijl}$ and sum on i, j. We do this because a product of two ε's
     sharing two indices collapses to a δ (D09). · *plain:* Contract both sides with the alternating tensor.
  4. *did:* Apply the 2δ contraction · *tex:* $\varepsilon_{ijl}\,R_{ij} = -2\,\delta_{lk}\,\omega_k = -2\,\omega_l$ ·
     *why:* $\varepsilon_{ijl}\varepsilon_{ijk} = 2\delta_{lk}$ (D09 step 8, the book's $\varepsilon_{pqi}\varepsilon_{pqj}
     = 2\delta_{ij}$), then δ substitutes. · *plain:* The contraction returns minus twice the vector.
  5. *did:* Solve for ω and rename · *tex:* $\omega_k = -\tfrac12\,\varepsilon_{ijk}\,R_{ij}$ · *why:* Divide by −2 and
     rename the free index l → k. This is the second half of (2.27); the sum runs over i = 1..3 and j = 1..3 (the book's
     lower limit "i−1" is a misprint for i = 1). · *plain:* Half the contraction of R with ε, sign reversed, gives ω
     back.
  6. *did:* Apply R to a vector x · *tex:* $R_{ij}\,x_j = -\varepsilon_{ijk}\,\omega_k\,x_j = +\varepsilon_{ikj}\,\omega_k\,x_j$
     · *why:* Insert the definition; then swap the last two indices of ε (one place, sign flips) to bring k before j —
     the order the cross-product formula needs. · *plain:* R times x is ε contracted with ω and x in the cross-product
     order.
  7. *did:* Recognise the cross product · *tex:* $R_{ij}\,x_j = (\boldsymbol\omega\times\mathbf x)_i$ · *why:* (2.21) in the
     form $(\mathbf u\times\mathbf v)_i = \varepsilon_{ikj}u_kv_j$ (free index first, then u's, then v's — a two-place move
     of $\varepsilon_{kji}$). With u = ω, v = x this is step 6. · *plain:* R rotates x about ω: R·x = ω × x.
- **Result.** $R_{ij} = -\varepsilon_{ijk}\omega_k$, $\omega_k = -\tfrac12\varepsilon_{ijk}R_{ij}$ (2.27), and $\mathbf R
  \cdot\mathbf x = \boldsymbol\omega\times\mathbf x$ — *in words:* an antisymmetric tensor *is* a vector, and acting with
  it means "cross with that vector".
- **Check.** Units: those of ω ✓. Round trip: ω = (1, 2, 3) → R → ω to 1e-16 (`antisymmetric_from_vector`,
  `vector_from_antisymmetric`). Number: ω = (0, 0, −1), x = (1, 0, 0): R = [[0, 1, 0],[−1, 0, 0],[0, 0, 0]], R·x = (0,
  −1, 0); ω × x = (0·0 − (−1)·0, (−1)·1 − 0·0, 0) = (0, −1, 0) ✓ (the sign-discrimination test of the analyst's F19).
- **What it means.** For the velocity gradient, the antisymmetric part A of Ch. 3 carries the vector $\tfrac12\nabla
  \times\mathbf u$ — half the vorticity (the book's R = 2A, (2.26)–(2.27), carries the vorticity ω itself) — and $\mathbf A\cdot\mathbf x = \tfrac12\boldsymbol\omega\times\mathbf x$ is the velocity
  of a solid-body rotation at angular velocity ω (Ex. 2.3 read backwards). That is *why* the antisymmetric part means
  "spin". The sign in (2.27) is chosen so that this works; with the opposite sign R·x = −ω × x.
- **Traps.** The sign in step 6 (a one-place move). The lower limit misprint. Forgetting that δ_lk ω_k = ω_l (not 3ω).

### D17 · Real eigenvalues, orthogonal axes, diagonal form and the bounds — for a real symmetric tensor — ★★★, 15 steps + sympy, in C13 (notebook · `cauchy_traction_principal_axes`)
- **Goal.** Prove the four facts §2.11 states without proof: for a real symmetric τ (1) the eigenvalues are real, (2)
  eigenvectors of distinct eigenvalues are perpendicular, (3) in the axes of the eigenvectors τ is diagonal, (4) the
  normal stress on every plane lies between the smallest and largest eigenvalue (and the shear is bounded by half
  their spread).
- **Start.** $\tau_{ij}\,b_j = \lambda\,b_i$ with $\tau_{ij} = \tau_{ji}$ real — *in words:* a direction b that τ only
  scales (no shear on the plane ⊥ b); λ and b may be complex until we prove otherwise.
- **Plan.** (1) Contract the eigen-equation with the conjugate eigenvector, conjugate, use symmetry: λ = λ̄. (2) Do the
  same with two eigenpairs: (λ¹ − λ²) b¹·b² = 0. (3) Put the unit eigenvectors as columns of C — exactly C02's C — and
  apply (2.12). (4) Expand any unit n in the eigenbasis and read the normal stress as a weighted mean of the λ's.
- **Tools.** Eigenvalues and eigenvectors (primer, C13) · complex conjugate and |z|² = z z̄ (primer, C13) · symmetry
  τ_ij = τ_ji (C12, N56) · dummy renaming (N12) · orthonormality and completeness (primer, C02) and D02 · the
  transformation rule (2.12) (C06) · quadratic form / Rayleigh quotient (primer, C13) · Gram–Schmidt for repeated λ
  (gloss in step 8).
- **Assumptions.** τ real (step 3) and symmetric (steps 4 and 7 — used twice, said both times) · the eigenvectors are
  normalised to unit length (steps 9, 12) · a right-handed choice of signs (step 9).
- **Steps.**
  1. *did:* Contract the eigen-equation with the conjugate eigenvector · *tex:* $\bar b_i\,\tau_{ij}\,b_j = \lambda\,\bar b_i
     b_i$ · *why:* Multiply both sides by $\bar b_i$ and sum on i; $\bar b_ib_i = \sum|b_i|^2 > 0$ for a nonzero b. We
     conjugate one factor so that a real number appears on the right. · *plain:* A complex "normal stress" of τ along b
     equals λ times a positive number.
  2. *did:* Take the complex conjugate of the whole line · *tex:* $b_i\,\bar\tau_{ij}\,\bar b_j = \bar\lambda\,b_i\bar b_i$ ·
     *why:* Conjugating a product conjugates each factor (primer complex conjugate); an equation stays true under
     conjugation. We do it to produce λ̄ next to λ. · *plain:* The mirror-image statement with λ̄.
  3. *did:* Use that τ is real · *tex:* $b_i\,\tau_{ij}\,\bar b_j = \bar\lambda\,b_i\bar b_i$ · *why:* $\bar\tau_{ij} =
     \tau_{ij}$ because the components are real numbers. Assumption "real" enters here. · *plain:* τ is unchanged by the
     mirror.
  4. *did:* Rename the dummies and use symmetry · *tex:* $b_i\,\tau_{ij}\,\bar b_j = \bar b_j\,\tau_{ji}\,b_i = \bar b_i\,
     \tau_{ij}\,b_j$ · *why:* Reorder the scalar factors, rename i ↔ j (both summed, N12), then $\tau_{ji} = \tau_{ij}$
     (symmetry, first use). The left side of step 3 has become the left side of step 1. · *plain:* Thanks to symmetry
     the two "normal stresses" are the same number.
  5. *did:* Subtract step 3 from step 1 · *tex:* $(\lambda - \bar\lambda)\,\bar b_ib_i = 0\ \Rightarrow\ \lambda = \bar
     \lambda$ · *why:* Equal left sides give equal right sides; $\bar b_ib_i > 0$ can be divided out, leaving λ = λ̄ —
     λ is real (fact 1). A real λ makes $(\tau - \lambda\delta)b = 0$ a real linear system, so b can be chosen real. ·
     *plain:* The eigenvalues of a real symmetric tensor are real numbers.
  6. *did:* Contract two eigen-equations crosswise · *tex:* $b^2_i\,\tau_{ij}\,b^1_j = \lambda^1\,b^2_ib^1_i,\qquad b^1_i\,
     \tau_{ij}\,b^2_j = \lambda^2\,b^1_ib^2_i$ · *why:* Take eigenpairs (λ¹, b¹) and (λ², b²) (superscripts are labels,
     not powers); contract the first equation with $b^2_i$ and the second with $b^1_i$. Same trick as step 1 with two
     vectors. · *plain:* Two mixed "normal stresses", each equal to an eigenvalue times b¹·b².
  7. *did:* Show the two left sides are equal · *tex:* $b^1_i\,\tau_{ij}\,b^2_j = b^2_j\,\tau_{ji}\,b^1_i = b^2_i\,\tau_{ij}\,
     b^1_j$ · *why:* Rename i ↔ j, then $\tau_{ji} = \tau_{ij}$ (symmetry, second use). Without symmetry this step fails
     and eigenvectors need not be perpendicular. · *plain:* τ between b¹ and b² does not care about the order.
  8. *did:* Subtract: orthogonality for distinct eigenvalues · *tex:* $(\lambda^1 - \lambda^2)\,\mathbf b^1\cdot\mathbf b^2
     = 0\ \Rightarrow\ \mathbf b^1\cdot\mathbf b^2 = 0\ \ (\lambda^1 \neq \lambda^2)$ · *why:* Equal left sides in step 6
     force the right sides equal. If λ¹ ≠ λ² the dot product must vanish (fact 2). If λ¹ = λ², any vector in the
     eigenplane is an eigenvector and an orthogonal pair can be chosen (Gram–Schmidt — `eigh` does this). · *plain:*
     Principal axes of different eigenvalues are perpendicular; equal ones let us choose them so.
  9. *did:* Build C from the unit eigenvectors · *tex:* $C_{ij} \equiv b^j_i = \mathbf e_i\cdot\mathbf b^j,\qquad \mathbf C^
     {\rm T}\mathbf C = \mathbf I$ · *why:* Take the new axes $\mathbf e'_j = \mathbf b^j$ (normalised, mutually
     perpendicular by step 8): then $C_{ij} = \mathbf e_i\cdot\mathbf e'_j$ is precisely C02's matrix, orthogonal by D02.
     Flip the sign of one b if det C = −1, to keep a rotation. · *plain:* The columns of C are the principal axes.
  10. *did:* Apply the transformation rule (2.12) · *tex:* $\tau'_{mn} = C_{im}\,C_{jn}\,\tau_{ij} = b^m_i\,\tau_{ij}\,b^n_j$ ·
      *why:* (2.12) for a second-order tensor (D06) with step 9's C inserted. The right side is "τ between b^m and b^n"
      — the quantity of step 6. · *plain:* The new components are τ sandwiched between two principal axes.
  11. *did:* Use the eigen-equation on b^n · *tex:* $\tau'_{mn} = \lambda^n\,b^m_ib^n_i = \lambda^n\,\delta_{mn}$ · *why:*
      $\tau_{ij}b^n_j = \lambda^nb^n_i$ (no sum on n), then orthonormality $\mathbf b^m\cdot\mathbf b^n = \delta_{mn}$
      (steps 8–9). So τ' is diagonal with the eigenvalues on the diagonal (fact 3). · *plain:* In the principal axes
      there is no shear at all — only λ¹, λ², λ³.
  12. *did:* Expand a unit normal in the eigenbasis · *tex:* $\mathbf n = c_k\,\mathbf b^k,\qquad c_k = \mathbf n\cdot
      \mathbf b^k,\qquad c_kc_k = 1$ · *why:* Completeness of the orthonormal set {b^k} (primer); $|\mathbf n|^2 = c_kc_k$
      by orthonormality. We expand n because τ is simple on the b's. · *plain:* Any cut direction is a mix of the three
      principal axes, with weights whose squares add to 1.
  13. *did:* Compute the normal stress on the plane ⊥ n · *tex:* $n_i\,\tau_{ij}\,n_j = c_kc_l\,\mathbf b^k\cdot\boldsymbol
      \tau\mathbf b^l = \lambda^k\,c_k^2$ · *why:* Insert step 12 twice; by steps 6–8 and 11, $\mathbf b^k\cdot\boldsymbol
      \tau\mathbf b^l = \lambda^l\delta_{kl}$, so only k = l survives. The normal stress is a weighted mean of the
      eigenvalues with weights $c_k^2$. · *plain:* The normal stress on any plane is an average of the principal
      stresses.
  14. *did:* Bound the weighted mean · *tex:* $\lambda_{\min} \le n_i\tau_{ij}n_j \le \lambda_{\max}$ · *why:* Replace every
      λ^k by λ_min (or λ_max) in $\lambda^kc_k^2$: the sum can only decrease (increase), and $c_kc_k = 1$ gives the bound
      (Rayleigh quotient, primer). Equality when n is the corresponding b^k (fact 4, sharp form). · *plain:* No cut
      carries more pull than λ_max or more squeeze than λ_min.
  15. *did:* Bound the shear too · *tex:* $\tau_s^2 = \lambda^k\lambda^k c_k^2 - (\lambda^kc_k^2)^2 \le \Big(\tfrac{\lambda_
      {\max} - \lambda_{\min}}{2}\Big)^2$ · *why:* $|\mathbf f|^2 = f_if_i = \lambda^k\lambda^kc_k^2$ by the same expansion,
      and $\tau_s^2 = |\mathbf f|^2 - \sigma_n^2$ is the *variance* of the values λ^k under the weights $c_k^2$; a variance
      of numbers confined to an interval of length L is at most (L/2)². · *plain:* The shear on any plane is at most
      half the spread of the eigenvalues (Mohr's circle radius).
- **Result.** For a real symmetric τ: real λ^k; orthonormal principal axes b^k; $\boldsymbol\tau' = \mathbf C^{\rm T}
  \boldsymbol\tau\mathbf C = \mathrm{diag}(\lambda^1, \lambda^2, \lambda^3)$ with $\mathbf C = [\mathbf b^1\ \mathbf b^2\
  \mathbf b^3]$; $\lambda_{\min} \le \mathbf n\cdot\boldsymbol\tau\cdot\mathbf n \le \lambda_{\max}$ and $\tau_s \le
  (\lambda_{\max} - \lambda_{\min})/2$ for every unit n — *in words:* a symmetric tensor is a pure stretch along three
  perpendicular axes, and those stretches bound what any plane can feel.
- **Check.** Units: those of τ ✓. Limit τ = −pδ: all λ = −p, every direction principal, τ_s = 0 ✓. Number (Ex. 2.4, Γ
  = 1): λ = ±1, b¹ = (1, 1)/√2, b² = (−1, 1)/√2, C = 45° rotation, S' = diag(1, −1); on the 30° plane n·S·n = sin 60° =
  0.866 ∈ [−1, 1], τ_s = cos 60° = 0.5 ≤ 1 ✓ (`principal_axes`, `diagonalize`, `normal_stress_bounds`). Sympy check
  cell below.
- **sympy check intent (`check_src`)** — re-runs the construction, not only the result, for a symbolic 2×2:
  ```python
  import sympy as sp                                            # symbolic algebra (Ch. 1 P40)
  p, q, r, th = sp.symbols('p q r theta', real=True)           # τ = [[p, q],[q, r]] real symmetric; θ = cut angle
  tau = sp.Matrix([[p, q], [q, r]])                             # the symmetric tensor
  lam = sp.symbols('lambda')
  poly = sp.expand((tau - lam*sp.eye(2)).det())                 # characteristic polynomial λ² − (p + r)λ + (pr − q²)
  disc = sp.expand(poly.coeff(lam, 1)**2 - 4*poly.coeff(lam, 0))   # its discriminant
  assert sp.simplify(disc - ((p - r)**2 + 4*q**2)) == 0         # step 5: (p−r)² + 4q² ≥ 0 → both roots real
  lam1, lam2 = [sp.simplify(s) for s in sp.solve(poly, lam)]    # the two eigenvalues
  b1 = sp.Matrix([q, lam1 - p]); b2 = sp.Matrix([q, lam2 - p])  # eigenvectors from the first row of (τ − λI)b = 0
  assert sp.simplify(b1.dot(b2)) == 0                           # step 8: orthogonal (uses p r − q² = λ1 λ2, Vieta)
  C = sp.Matrix.hstack(b1/sp.sqrt(b1.dot(b1)), b2/sp.sqrt(b2.dot(b2)))  # step 9: unit eigenvectors (sqrt(b·b), not .norm(): .norm() adds Abs() and blocks simplify)
  D = sp.simplify(C.T*tau*C)                                    # step 11: τ' = CᵀτC
  assert sp.simplify(D[0, 1]) == 0 and sp.simplify(D[1, 0]) == 0   # off-diagonals vanish
  assert sp.simplify(D[0, 0] - lam1) == 0 and sp.simplify(D[1, 1] - lam2) == 0   # diagonal = eigenvalues
  n = sp.Matrix([sp.cos(th), sp.sin(th)])                       # step 12: any unit normal
  c1, c2 = n.dot(C[:, 0]), n.dot(C[:, 1])                       # its weights on the principal axes
  assert sp.simplify((n.T*tau*n)[0] - (lam1*c1**2 + lam2*c2**2)) == 0   # step 13: σ_n = Σ λ_k c_k²
  assert sp.simplify(c1**2 + c2**2 - 1) == 0                    # weights add to 1 → the bounds of step 14 follow
  print("D17 verified symbolically: real λ, orthogonal b, CᵀτC diagonal, σ_n = Σ λ_k c_k²")
  ```
- **What it means.** Ch. 3: the strain rate has three perpendicular directions of pure stretching (no shearing), with
  rates λ^k; Ch. 4: principal stresses and Mohr's circle; Ch. 11: normal modes of a symmetric operator. The book's
  loose "λ bound the elements τ_ij" is true but the useful statement is about the normal stress on planes and the
  shear bound. Every fact fails for a non-symmetric tensor (complex λ, skew eigenvectors) — `principal_axes` refuses
  one.
- **Traps.** Symmetry is used twice (steps 4 and 7) — one use gives realness, the other orthogonality. Repeated
  eigenvalues: orthogonal eigenvectors exist but are not unique (test with projectors, not vectors). det C = −1 from an
  unlucky sign choice: flip one b. λ^k is a label, not a power.

### D18 · Why the trace, I₂ and the determinant do not depend on the axes — ★★, 8 steps, in C07 (notebook)
- **Goal.** Show that three combinations of the nine components of a tensor are the same for every observer, that they
  are the coefficients of the characteristic polynomial, and what they become in the principal frame. The book gives
  hints in Exercise 2.9.
- **Start.** $A'_{mn} = C_{im}\,C_{jn}\,A_{ij}$ (2.12) — *in words:* how any second-order tensor's components change
  under a rotation (D06).
- **Plan.** (1) Contract (2.12) on m = n: a C-pair collapses to δ. (2) Do the same for the chain A_ij A_ji. (3) Take
  determinants. (4) Expand det(A − λδ) and identify the coefficients; evaluate them in the principal frame.
- **Tools.** The transformation rule (2.12) (C06) · orthogonality as "a contracted C-pair is δ", $C_{im}C_{jm} = \delta_
  {ij}$ (D02, C02) · Kronecker substitution (C01) · det(AB) = det A det B (P53 reminder) · Vieta's formulas (primer, C07)
  · the principal frame (D17 fact 3, forward pointer; step 8 only needs "τ' = diag(λ)").
- **Assumptions.** A any second-order tensor (steps 1–7); for step 8 A symmetric (so that a principal frame exists).
- **Steps.**
  1. *did:* Contract (2.12) on its free indices · *tex:* $A'_{mm} = C_{im}\,C_{jm}\,A_{ij}$ · *why:* Set n = m and sum
     (contraction, C07). The two C's now share their summed second index m. We contract because a fully contracted
     quantity has no free index — a candidate scalar. · *plain:* The trace in the new frame, written with the old
     components.
  2. *did:* Collapse the C-pair · *tex:* $A'_{mm} = \delta_{ij}\,A_{ij} = A_{ii}$ · *why:* $C_{im}C_{jm} = \delta_{ij}$ (D02,
     second indices summed), then δ substitutes. The trace is the same in both frames: $I_1 = A_{ii}$ is invariant. ·
     *plain:* Every observer measures the same trace.
  3. *did:* Transform the two-link chain · *tex:* $A'_{mn}A'_{nm} = C_{im}C_{jn}A_{ij}\;C_{kn}C_{lm}A_{kl}$ · *why:* Apply
     (2.12) to each factor with fresh dummy letters (i, j) and (k, l) — reusing a letter would wrongly sum it. We choose
     $A_{ij}A_{ji}$ because its indices close into a loop. · *plain:* The chain in the new frame, with four C's.
  4. *did:* Pair the C's and collapse · *tex:* $A'_{mn}A'_{nm} = (C_{im}C_{lm})(C_{jn}C_{kn})\,A_{ij}A_{kl} = \delta_{il}
     \delta_{jk}A_{ij}A_{kl} = A_{ij}A_{ji}$ · *why:* Scalars commute, so group the C's sharing m and those sharing n;
     each pair is a δ (D02); the δ's substitute l → i and k → j. Hence $I_2 = \tfrac12(I_1^2 - A_{ij}A_{ji})$ is
     invariant. · *plain:* Any closed chain of indices survives a rotation unchanged.
  5. *did:* Take the determinant of the matrix form · *tex:* $\det\mathbf A' = \det(\mathbf C^{\rm T}\mathbf A\,\mathbf C) =
     (\det\mathbf C)^2\det\mathbf A = \det\mathbf A$ · *why:* det of a product is the product of dets (P53); det Cᵀ = det
     C and (det C)² = 1 (D02 step 7). So $I_3 = \det\mathbf A$ is invariant. · *plain:* The determinant is the third
     frame-independent number.
  6. *did:* Expand the characteristic determinant in λ · *tex:* $\det(\mathbf A - \lambda\boldsymbol\delta) = -\lambda^3 +
     A_{ii}\,\lambda^2 - M\,\lambda + \det\mathbf A$ · *why:* Expanding a 3×3 determinant whose diagonal entries are
     $A_{ii} - \lambda$ (cofactors, P53): the λ³ and λ² coefficients come from the diagonal product, the constant term
     is det A, and M is the sum of the three principal 2×2 minors. · *plain:* A cubic whose coefficients are built
     from A alone.
  7. *did:* Identify M with I₂ · *tex:* $M = (A_{11}A_{22} - A_{12}A_{21}) + (A_{22}A_{33} - A_{23}A_{32}) + (A_{11}A_{33} -
     A_{13}A_{31}) = \tfrac12(A_{ii}A_{jj} - A_{ij}A_{ji})$ · *why:* Expand $A_{ii}A_{jj} = (A_{11} + A_{22} + A_{33})^2$ and
     subtract $A_{ij}A_{ji}$ (diagonal squares and the pairs $A_{12}A_{21}$ …): the diagonal squares cancel, each cross
     product appears twice, hence the ½. So the cubic is $\lambda^3 - I_1\lambda^2 + I_2\lambda - I_3 = 0$ (sign
     flipped). · *plain:* The middle coefficient is exactly the second invariant.
  8. *did:* Evaluate in the principal frame · *tex:* $I_1 = \lambda^1 + \lambda^2 + \lambda^3,\quad I_2 = \lambda^1\lambda^2
     + \lambda^2\lambda^3 + \lambda^3\lambda^1,\quad I_3 = \lambda^1\lambda^2\lambda^3$ · *why:* For symmetric A the
     principal frame has A' = diag(λ) (D17 fact 3); the invariants may be computed in any frame, so use this one:
     trace, sum of 2×2 minors, determinant of a diagonal matrix — Vieta's formulas for the cubic's roots. · *plain:*
     The invariants are the sum, the pair-sums and the product of the eigenvalues — so the eigenvalues are
     frame-independent too.
- **Result.** $I_1 = A_{ii}$, $I_2 = \tfrac12(I_1^2 - A_{ij}A_{ji})$, $I_3 = \det\mathbf A$ are unchanged by (2.12);
  $\det(\mathbf A - \lambda\boldsymbol\delta) = -(\lambda^3 - I_1\lambda^2 + I_2\lambda - I_3)$; in the principal frame they
  are the symmetric functions of the λ's — *in words:* three numbers, and hence the eigenvalues, belong to the tensor,
  not to the observer.
- **Check.** Units: A, A², A³ respectively ✓. Limit A = δ: I = (3, 3, 1), polynomial $(1 - \lambda)^3$ ✓. Number: A =
  [[2, 1, 0],[1, 3, 1],[0, 1, 4]]: I₁ = 9, $A_{ij}A_{ji} = 29 + 4 = 33$, I₂ = ½(81 − 33) = 24, I₃ = 18; roots 1.268,
  3.000, 4.732 sum to 9.000 and multiply to 18.00 ✓ (`invariants`, `characteristic_polynomial`, `principal_axes`).
- **What it means.** The pressure $-\tau_{ii}/3$ (Ch. 4), the volume strain rate $S_{ii} = \nabla\cdot\mathbf u$ (Ch. 3),
  and the turbulence invariants of Ch. 12 are physical because they are contractions. $A_{ij}A_{ij} = \mathrm{tr}(A A^T)$ is
  also frame-independent, but it is not the coefficient of λ in the cubic — I₂ needs the closed chain $A_{ij}A_{ji}$ (equal to it only for symmetric A).
- **Traps.** Reusing a dummy letter in step 3. Using $A_{ij}A_{ij}$ for I₂. The sign pattern of the cubic. Reading λ^k
  as a power.

### D21 · The integral definitions as small-volume limits of Gauss' theorem, Eqs. (2.31)–(2.33) — ★★, 7 steps, in C15 (notebook · `gauss_flux_box`)
- **Goal.** Turn Gauss' theorem, read on a shrinking volume, into coordinate-free definitions of the gradient, the
  divergence and the curl — the book calls (2.31) "the limiting form" and leaves the limit implicit.
- **Start.** $\displaystyle\iiint_V \frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA$ (2.30) on a small volume V
  around x₀ — *in words:* Gauss' theorem (D25) for any Q, applied to a small region.
- **Plan.** (1) Replace the volume integral by "value at some interior point × volume" (mean-value theorem). (2)
  Divide by V. (3) Shrink V to x₀. (4) Contract with e_i for the divergence, cross with ε for the curl.
- **Tools.** Gauss' theorem (2.30) (C14, D25) · mean-value theorem for integrals (primer, C15) · limits and continuity
  (primer orders of smallness, C05) · the index form of ∇·, (2.23) (C10) · the index form of ∇×, (2.24) (C11).
- **Assumptions.** ∂Q/∂x_i continuous near x₀ (steps 2 and 4) · V shrinks to x₀ in every direction (step 4; then the
  shape does not matter).
- **Steps.**
  1. *did:* Apply Gauss' theorem to a small V around x₀ · *tex:* $\displaystyle\iiint_V \frac{\partial Q}{\partial x_i}\,dV
     = \iint_A n_iQ\,dA$ · *why:* (2.30) holds for any volume with a piecewise smooth boundary; we pick a small one
     around the point where we want a derivative. · *plain:* The inside sum of the derivative equals the boundary sum
     of Q.
  2. *did:* Replace the volume integral by an interior value times V · *tex:* $\displaystyle\iiint_V \frac{\partial Q}{
     \partial x_i}\,dV = \frac{\partial Q}{\partial x_i}(\mathbf x^*)\;V$ · *why:* Mean-value theorem for integrals: a
     continuous integrand over V equals its value at some point x* in V times the volume (primer). We do this to turn
     an integral into a single value. · *plain:* The integral is the volume times the derivative somewhere inside.
  3. *did:* Divide by V · *tex:* $\displaystyle\frac{\partial Q}{\partial x_i}(\mathbf x^*) = \frac1V\iint_A n_iQ\,dA$ ·
     *why:* V > 0. The right side is now "boundary flux of Q per unit volume", a ratio that stays finite as V shrinks. ·
     *plain:* The derivative at an interior point equals the flux per unit volume.
  4. *did:* Shrink V to the point · *tex:* $\displaystyle\frac{\partial Q}{\partial x_i}(\mathbf x_0) = \lim_{V\to0}\frac1V
     \iint_A n_iQ\,dA \equiv \mathcal D Q$ · *why:* x* is trapped in V, so x* → x₀; continuity of ∂Q/∂x_i lets the left
     side tend to its value at x₀. The left side never knew the shape of V, so neither does the limit. This is (2.31). ·
     *plain:* The derivative of Q is the limit of the boundary flux per unit volume — no coordinates needed.
  5. *did:* Specialise to a vector Q and contract · *tex:* $\displaystyle\frac{\partial Q_i}{\partial x_i} = \lim_{V\to0}
     \frac1V\iint_A n_iQ_i\,dA$ · *why:* Put Q → Q_j in step 4 (Gauss holds component by component, the theorem being
     linear in Q), then set j = i and sum — contraction on both sides. · *plain:* Adding the three diagonal derivatives
     adds the three flux components.
  6. *did:* Recognise the divergence and the outflux · *tex:* $\displaystyle\nabla\cdot\mathbf Q = \lim_{V\to0}\frac1V
     \iint_A \mathbf n\cdot\mathbf Q\,dA$ · *why:* $\partial Q_i/\partial x_i$ is (2.23) and $n_iQ_i = \mathbf n\cdot\mathbf
     Q$ is the outward flux density. This is (2.32). · *plain:* Divergence = what leaks out per unit volume.
  7. *did:* Cross instead of dot · *tex:* $\displaystyle(\nabla\times\mathbf Q)_k = \varepsilon_{kij}\frac{\partial Q_j}{
     \partial x_i} = \lim_{V\to0}\frac1V\iint_A (\mathbf n\times\mathbf Q)_k\,dA$ · *why:* Multiply step 4 (with Q → Q_j)
     by $\varepsilon_{kij}$ and sum on i, j — a fixed linear combination passes through the limit; $\varepsilon_{kij}
     \partial_iQ_j$ is (2.24) and $\varepsilon_{kij}n_iQ_j = (\mathbf n\times\mathbf Q)_k$ by (2.21). This is (2.33). ·
     *plain:* Curl = the boundary's "n × Q" per unit volume.
- **Result.** $\mathcal DQ = \lim (1/V)\iint n_iQ\,dA$ (2.31), $\nabla\cdot\mathbf Q = \lim (1/V)\iint \mathbf n\cdot\mathbf
  Q\,dA$ (2.32), $\nabla\times\mathbf Q = \lim (1/V)\iint \mathbf n\times\mathbf Q\,dA$ (2.33) — *in words:* gradient,
  divergence and curl are boundary sums per unit volume in the limit of a vanishing volume.
- **Check.** Units: (2.32) — [Q]·m²/m³ = [Q]/m on both sides ✓. Limit: Q constant: the boundary integral of n vanishes
  (vector area of a closed surface, C05 primer), so all three derivatives are 0 ✓. Number: Q = (x², 0, 0), x₀ = (1, 0,
  0), cube h = 0.2: (1/V)∮ = (1.21 − 0.81)(0.04)/0.008 = 2.000 = 2x₀ exactly ✓; Q = (x³, 0, 0): 3.01 vs 3 at h = 0.2,
  error h²/4, order 2 (`integral_divergence`, `integral_definition_convergence("divergence")`).
- **What it means.** The operators of §2.9 are properties of the field, not of Cartesian axes: the same recipe in
  spherical or cylindrical coordinates gives the curvilinear formulas of Appendix B (N72), and Ch. 4 derives continuity
  from exactly (2.32) on a fluid box. The limit exists only where Q is smooth; at a point source (E4) the flux is fixed
  while V → 0 and (1/V)∮ diverges — a delta function, not a derivative.
- **Traps.** Forgetting the 1/V. Assuming shape-independence without continuity. Reading (2.32) and (2.33) as new
  theorems — they are (2.31) contracted or crossed.

### D22 · Ex. 2.5: the Cartesian divergence recovered from the outflux definition — ★★, 9 steps, in C15 (notebook · `gauss_flux_box`)
- **Goal.** Evaluate (2.32) on a small box and watch the formula $\partial Q_i/\partial x_i$ (2.23) come back — face
  by face — so that the coordinate formula is seen as a consequence of "outflux per volume".
- **Start.** $\displaystyle\nabla\cdot\mathbf Q = \lim_{V\to0}\frac1V\iint_A \mathbf n\cdot\mathbf Q\,dA$ (2.32) on a box
  $\Delta x_1\Delta x_2\Delta x_3$ centred at x — *in words:* the outflux definition, applied to a small rectangular
  box with faces perpendicular to the axes (Fig. 2.4's letters: EADH/FBCG ⊥ x₁, ABCD/EFGH ⊥ x₂, ABFE/DCGH ⊥ x₃).
- **Plan.** (1) Take the two faces ⊥ x₁; write their outward normals. (2) Taylor-expand Q from the centre to each face
  centre. (3) Integrate over each face (midpoint rule). (4) Add the pair: constants cancel, slopes add. (5) Same for
  the other two pairs; add; divide by the volume; take the limit.
- **Tools.** First-order Taylor expansion (P26 reminder) · midpoint rule: a face integral of a linear function equals
  the centre value × area (primer volume and surface integrals as midpoint sums, C14) · face normals ±e₁ (C04's faces)
  · orders of smallness (primer, C05) · e₁·∂Q/∂x₁ = ∂Q₁/∂x₁ (gloss in step 9) · (2.23) (C10).
- **Assumptions.** Q twice differentiable (Taylor remainder O(Δx²); steps 2–3) · the box shrinks in all three
  directions (step 8).
- **Steps.**
  1. *did:* Name the two faces perpendicular to x₁ · *tex:* $\text{EADH}:\ \mathbf n = +\mathbf e_1\ \text{at}\ x_1 + \tfrac
     {\Delta x_1}2;\qquad \text{FBCG}:\ \mathbf n = -\mathbf e_1\ \text{at}\ x_1 - \tfrac{\Delta x_1}2$ · *why:* Outward
     normals of a box are ±e_i (C04's faces); each of these two faces has area $\Delta x_2\Delta x_3$. We start with
     one pair because the other two are copies. · *plain:* A "front" face facing +x₁ and a "back" face facing −x₁.
  2. *did:* Taylor-expand Q to the two face centres · *tex:* $[\mathbf Q]_{\rm EADH} = \mathbf Q(\mathbf x) + \tfrac{\Delta
     x_1}2\,\dfrac{\partial\mathbf Q}{\partial x_1} + \dots,\qquad [\mathbf Q]_{\rm FBCG} = \mathbf Q(\mathbf x) - \tfrac{
     \Delta x_1}2\,\dfrac{\partial\mathbf Q}{\partial x_1} + \dots$ · *why:* First-order Taylor expansion from the box
     centre along ±x₁ (P26); the dots are O(Δx₁²). We expand from the centre so the two faces get opposite first-order
     terms. · *plain:* The front face sees Q plus half a step of slope; the back face sees Q minus it.
  3. *did:* Integrate n·Q over each face by the midpoint rule · *tex:* $\displaystyle\iint_{\rm EADH}\mathbf n\cdot\mathbf Q
     \,dA = \mathbf e_1\cdot[\mathbf Q]_{\rm EADH}\,\Delta x_2\Delta x_3 + O(\Delta x^4)$ · *why:* Over a face, Q varies
     linearly to first order, and the integral of a linear function over a rectangle equals its centre value times the
     area (midpoint rule, exact for linear variation); the neglected variation is O(Δx²) relative, i.e. O(Δx⁴)
     absolute here. · *plain:* Face integral ≈ centre value × face area.
  4. *did:* Add the front and back faces · *tex:* $\big(\mathbf e_1\cdot[\mathbf Q]_{\rm EADH} - \mathbf e_1\cdot[\mathbf Q]_
     {\rm FBCG}\big)\Delta x_2\Delta x_3$ · *why:* The back face has n = −e₁, so its n·Q is minus e₁·Q. Pairing opposite
     faces is the whole trick: their zeroth-order terms are about to cancel. · *plain:* Front minus back, times the
     face area.
  5. *did:* Insert the Taylor values: the constants cancel, the slopes add · *tex:* $\Big(\mathbf e_1\cdot\dfrac{\partial
     \mathbf Q}{\partial x_1}\Big)\Delta x_1\Delta x_2\Delta x_3 + O(\Delta x^5)$ · *why:* $\mathbf Q(\mathbf x) - \mathbf Q
     (\mathbf x) = 0$; $\tfrac{\Delta x_1}2 + \tfrac{\Delta x_1}2 = \Delta x_1$. Only the first-order change of Q across
     the box survives, multiplied by the volume $V = \Delta x_1\Delta x_2\Delta x_3$. · *plain:* What leaks through this
     pair is the slope of Q₁ along x₁ times the volume.
  6. *did:* Repeat for the pairs perpendicular to x₂ and x₃ · *tex:* $\Big(\mathbf e_2\cdot\dfrac{\partial\mathbf Q}{\partial
     x_2}\Big)V + O(\Delta x^5),\qquad \Big(\mathbf e_3\cdot\dfrac{\partial\mathbf Q}{\partial x_3}\Big)V + O(\Delta x^5)$ ·
     *why:* Steps 1–5 with the axis relabelled: faces ABCD/EFGH (±e₂) and ABFE/DCGH (±e₃). Nothing in the argument
     used which axis was "1". · *plain:* Each pair of faces contributes its own slope times the volume.
  7. *did:* Add the six faces · *tex:* $\displaystyle\iint_A \mathbf n\cdot\mathbf Q\,dA = \Big(\mathbf e_1\cdot\dfrac{\partial
     \mathbf Q}{\partial x_1} + \mathbf e_2\cdot\dfrac{\partial\mathbf Q}{\partial x_2} + \mathbf e_3\cdot\dfrac{\partial
     \mathbf Q}{\partial x_3}\Big)V + O(\Delta x^5)$ · *why:* The closed surface of the box is exactly the six faces;
     the surface integral is the sum of the face integrals (steps 5–6). · *plain:* The total outflux is the three slopes
     times the volume, plus something much smaller.
  8. *did:* Divide by V and take the limit · *tex:* $\displaystyle\lim_{V\to0}\frac1V\iint_A \mathbf n\cdot\mathbf Q\,dA =
     \mathbf e_i\cdot\dfrac{\partial\mathbf Q}{\partial x_i}$ · *why:* $O(\Delta x^5)/V = O(\Delta x^2) \to 0$ as all three
     sides shrink (orders of smallness), while the slope terms do not depend on the box size. Summation convention on
     i. · *plain:* Outflux per unit volume tends to the sum of the three directional slopes.
  9. *did:* Recognise the Cartesian divergence · *tex:* $\displaystyle\nabla\cdot\mathbf Q = \dfrac{\partial Q_1}{\partial x_1}
     + \dfrac{\partial Q_2}{\partial x_2} + \dfrac{\partial Q_3}{\partial x_3} = \dfrac{\partial Q_i}{\partial x_i}$ · *why:*
     $\mathbf e_i\cdot\partial\mathbf Q/\partial x_i = \partial(\mathbf e_i\cdot\mathbf Q)/\partial x_i = \partial Q_i/\partial
     x_i$ because the unit vectors are constant (Cartesian). This is (2.23): the recipe returns from the definition. ·
     *plain:* Divergence = the three stretching rates added — now *because* that is what leaks out of a tiny box.
- **Result.** $\nabla\cdot\mathbf Q = \partial Q_i/\partial x_i$ (2.23) from (2.32) — *in words:* the Cartesian formula is
  the outflux per unit volume of a shrinking box; opposite faces cancel the value of Q and keep its slope.
- **Check.** Units: [Q]/m ✓. Limit: linear Q — the O(Δx⁵) terms vanish identically and the box formula is exact at any
  size ✓. Number: Q = (x², 0, 0), x₀ = 1, h = 0.2: front 1.21·0.04, back −0.81·0.04, net 0.016, /0.008 = 2.000 = 2x₀ ✓;
  Q = (x³, 0, 0): 3 + h²/4 → order 2 in h (`integral_divergence`, the C15 animation).
- **What it means.** Ch. 4 §4.2 runs this argument with Q = ρu on a fluid box to get the continuity equation; the same
  bookkeeping with n × Q gives the curl (Ex. 2.6, D26). In curvilinear coordinates the face areas themselves vary and
  extra terms appear — Appendix B.
- **Traps.** Cancelling zeroth-order terms between *non*-opposite faces. Treating a face value as exact rather than a
  face average (the midpoint rule is what makes step 3 honest). Which face carries +e₁. Writing e₁·∂Q/∂x₁ = |∂Q/∂x₁|.

### D25 · Gauss' theorem, from the fundamental theorem of calculus, Eq. (2.30) — ★★, 9 steps, in C14 (notebook · `gauss_flux_box`)
- **Goal.** Prove the theorem the book only states: the volume integral of a derivative equals the surface integral of
  the field with the outward normal — first for a box, then for any volume by tiling.
- **Start.** $\displaystyle\int_a^b f'(x)\,dx = f(b) - f(a)$ — *in words:* the fundamental theorem of calculus (primer,
  C14): integrating a derivative gives the values at the two ends.
- **Plan.** (1) On a box, write the volume integral of ∂Q/∂x₁ as an iterated integral and do the x₁ integral first. (2)
  Read the right side of (2.30) face by face: n₁ = ±1 on two faces, 0 on four. (3) Repeat for i = 2, 3. (4) Tile a
  general volume with boxes; interior faces cancel.
- **Tools.** Fundamental theorem of calculus (primer, C14) · iterated integrals (primer volume and surface integrals as
  midpoint sums, C14) · outward normals of a box (C04's faces) · additivity of integrals over disjoint regions (step 7)
  · vector area / opposite normals (C05 primer; step 8).
- **Assumptions.** Q continuously differentiable on V and its boundary (steps 2, 9) · the boundary A piecewise smooth,
  so that the tiling's outer faces approach it (step 9).
- **Steps.**
  1. *did:* Write the i = 1 volume integral on a box as iterated integrals · *tex:* $\displaystyle\iiint_V \frac{\partial Q}{
     \partial x_1}\,dV = \int_{a_3}^{b_3}\!\int_{a_2}^{b_2}\Big[\int_{a_1}^{b_1}\frac{\partial Q}{\partial x_1}\,dx_1\Big]
     dx_2\,dx_3$ · *why:* On a box $[a_1, b_1]\times[a_2, b_2]\times[a_3, b_3]$ a volume integral is three nested
     one-dimensional integrals in any order; we do x₁ innermost because the derivative is in x₁. · *plain:* Integrate
     along x₁ first, on every line of the box.
  2. *did:* Apply the fundamental theorem to the inner integral · *tex:* $\displaystyle\int_{a_1}^{b_1}\frac{\partial Q}{
     \partial x_1}\,dx_1 = Q(b_1, x_2, x_3) - Q(a_1, x_2, x_3)$ · *why:* Along a line with x₂, x₃ fixed, ∂Q/∂x₁ is an
     ordinary derivative of a function of x₁; the fundamental theorem gives the two end values (Q ∈ C¹). · *plain:* Each
     line contributes Q at its far end minus Q at its near end.
  3. *did:* Insert and split into two face integrals · *tex:* $\displaystyle\iiint_V \frac{\partial Q}{\partial x_1}\,dV =
     \iint Q(b_1, x_2, x_3)\,dx_2dx_3 - \iint Q(a_1, x_2, x_3)\,dx_2dx_3$ · *why:* Step 2 into step 1; the remaining
     double integrals run over the face $[a_2, b_2]\times[a_3, b_3]$ — the two faces perpendicular to x₁. · *plain:* The
     inside adds up to "Q on the far face minus Q on the near face".
  4. *did:* Read the right side of (2.30) face by face · *tex:* $n_1 = +1\ \text{on}\ x_1 = b_1,\quad n_1 = -1\ \text{on}\
     x_1 = a_1,\quad n_1 = 0\ \text{on the other four faces}$ · *why:* The outward normal of a box face is ±e_i; its
     first component is nonzero only on the two faces perpendicular to x₁. Four of the six faces contribute nothing to
     the i = 1 term. · *plain:* Only the two faces facing ±x₁ count for the x₁-derivative.
  5. *did:* Conclude (2.30) for i = 1 on the box · *tex:* $\displaystyle\iint_A n_1Q\,dA = \iint Q(b_1,\cdot)\,dA - \iint Q
     (a_1,\cdot)\,dA = \iiint_V\frac{\partial Q}{\partial x_1}\,dV$ · *why:* Step 4 makes the surface integral exactly
     the two face integrals of step 3 with the signs of n₁. Left and right sides of (2.30) coincide for i = 1. ·
     *plain:* Gauss' theorem holds for the x₁-derivative on a box.
  6. *did:* Repeat for i = 2 and i = 3 · *tex:* $\displaystyle\iiint_V\frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA
     \quad\text{on any box, } i = 1, 2, 3$ · *why:* Steps 1–5 with the roles of the axes relabelled (integrate along
     x_i first; the faces ⊥ x_i carry n_i = ±1). Nothing used which axis was first. · *plain:* The same for each of the
     three directions.
  7. *did:* Tile a general volume with small boxes and add · *tex:* $\displaystyle\sum_k\iiint_{V_k}\frac{\partial Q}{
     \partial x_i}\,dV = \sum_k\iint_{A_k} n_iQ\,dA,\qquad \sum_k\iiint_{V_k} = \iiint_{V_{\rm tiles}}$ · *why:* Step 6
     holds on each box V_k; adding the equations adds both sides. Integrals over disjoint regions add up to the
     integral over their union. · *plain:* Gauss on every tile, summed.
  8. *did:* Cancel the interior faces · *tex:* $\displaystyle\sum_k\iint_{A_k} n_iQ\,dA = \iint_{\text{outer faces}} n_iQ
     \,dA$ · *why:* A face shared by two tiles appears twice with opposite outward normals (n and −n) and the same
     values of Q (continuity), so the two contributions sum to zero. Only faces on the outside of the tiling survive. ·
     *plain:* Every inner wall is counted once from each side and cancels.
  9. *did:* Refine the tiling · *tex:* $\displaystyle\iiint_V\frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA$ · *why:*
     As the boxes shrink, the tiled region fills V and its staircase boundary tends to A with the outer normals tending
     to n (A piecewise smooth, Q continuous), so both sides converge to the integrals over V and A. Q may carry any
     other indices: nothing depended on its order. · *plain:* Gauss' theorem for any volume and any field.
- **Result.** $\displaystyle\iiint_V\frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA$ (2.30) for Q of any order —
  *in words:* a derivative integrated over the inside equals the field weighted by the outward normal on the boundary.
- **Check.** Units: [Q]·m² both sides ✓. Limit Q constant: left 0; right ∮ n_i dA = 0 (vector area of a closed surface)
  ✓. Number: Q = (x, y, z) on the unit cube: ∭ 3 dV = 3; faces +1 + 0 + 1 + 0 + 1 + 0 = 3 ✓ (`divergence_theorem_box`);
  sphere benchmark F = (2x, y², z²): 8π/3 both sides (`divergence_theorem_sphere`); `divergence_theorem_tiled` reports
  interior faces summing to 1e-15.
- **What it means.** Every conservation law of Ch. 4 passes through this line: the rate of change of what is inside a
  volume equals the flux through its boundary, and shrinking the volume (D21) turns the integral law into a partial
  differential equation. The theorem needs Q differentiable *inside* — a point source or a shock breaks it and adds a
  surface term.
- **Traps.** Thinking all six faces contribute to each i. Cancelling interior faces without opposite normals and equal
  Q. Forgetting that the tiles' boundary is only an approximation of A until the limit.

### D26 · Stokes' theorem for a planar surface: one rectangle, then tiling, Eq. (2.34) — ★★, 10 steps, in C16 (notebook · `stokes_circulation_loop`)
- **Goal.** Prove Stokes' theorem for a flat surface — the circulation round the boundary equals the curl integrated
  inside — by doing one tiny rectangle exactly (Ex. 2.6's bookkeeping run backwards) and tiling; then read off the
  corollary ∇×∇φ = 0.
- **Start.** A rectangle $\Delta x_1\times\Delta x_2$ centred at $(x_1, x_2)$ in the plane $x_3 = $ const, with n = e₃
  and t counterclockwise — *in words:* the smallest loop we can integrate around by hand.
- **Plan.** (1) Fix the orientation of the four sides. (2) Midpoint values on each side. (3) Add bottom + top and right
  + left: Taylor differences. (4) Recognise the curl. (5) Tile the surface: interior edges cancel; refine. (6) Corollary
  for gradient fields.
- **Tools.** Line integral of a vector field round a loop (primer, C16) · orientation t counterclockwise about n (N69,
  right-hand-rule primer, C08) · first-order Taylor expansion (P26 reminder) · midpoint values on a side (C14 primer,
  as in D22) · the curl's third component (2.25) (C11, D12) · Ex. 2.6's bookkeeping (N73) · the chain rule for dφ along
  a curve (P49 reminder; step 10).
- **Assumptions.** u continuously differentiable on A (steps 2–3, 8) · A a flat surface with piecewise smooth boundary
  C (steps 7–9; a curved A is done patch by patch with t = n_c × n — stated) · every tile oriented by the same n
  (step 8).
- **Steps.**
  1. *did:* Fix the four sides and their tangents · *tex:* $\text{bottom}\ (x_2 - \tfrac{\Delta x_2}2):\ \mathbf t = +\mathbf
     e_1;\quad \text{right}:\ +\mathbf e_2;\quad \text{top}:\ -\mathbf e_1;\quad \text{left}:\ -\mathbf e_2$ · *why:*
     Counterclockwise about n = e₃ (right-hand rule: thumb along n, fingers along t — N69). The orientation decides
     every sign below, so we fix it first. · *plain:* Walk the rectangle anticlockwise seen from above.
  2. *did:* Write each side's integral with its midpoint value · *tex:* $\displaystyle\int_{\rm side}\mathbf u\cdot\mathbf t
     \,ds = (\mathbf u\cdot\mathbf t)_{\rm mid}\times\text{length} + O(\Delta^3)$ · *why:* Along a short side u varies
     linearly to first order and the integral of a linear function equals its midpoint value times the length (as in
     D22 step 3); the error is O(Δ³) per side. · *plain:* Each side contributes "velocity along it at its middle times
     its length".
  3. *did:* Add bottom and top · *tex:* $\big[u_1(x_1, x_2 - \tfrac{\Delta x_2}2) - u_1(x_1, x_2 + \tfrac{\Delta x_2}2)\big]
     \Delta x_1$ · *why:* On the bottom $\mathbf u\cdot\mathbf t = +u_1$, on the top $-u_1$ (step 1); both sides have
     length Δx₁ and midpoints directly below/above the centre. · *plain:* The horizontal sides see u₁ below minus u₁
     above.
  4. *did:* Taylor-expand the difference · *tex:* $-\,\dfrac{\partial u_1}{\partial x_2}\,\Delta x_1\Delta x_2 + O(\Delta^4)$ ·
     *why:* $u_1(x_2 \mp \tfrac{\Delta x_2}2) = u_1 \mp \tfrac{\Delta x_2}2\,\partial_2u_1 + \dots$ (P26); the two u₁'s
     cancel and the halves add to Δx₂ with a minus sign. · *plain:* If u₁ grows upward, the top side walks against it:
     a clockwise contribution.
  5. *did:* Add right and left the same way · *tex:* $\big[u_2(x_1 + \tfrac{\Delta x_1}2, x_2) - u_2(x_1 - \tfrac{\Delta x_1}
     2, x_2)\big]\Delta x_2 = +\dfrac{\partial u_2}{\partial x_1}\,\Delta x_1\Delta x_2 + O(\Delta^4)$ · *why:* On the right
     $\mathbf u\cdot\mathbf t = +u_2$, on the left $-u_2$; Taylor along x₁; this time the plus side is at $+\Delta x_1/2$
     so the sign is positive. · *plain:* If u₂ grows to the right, the right side walks with it: an anticlockwise
     contribution.
  6. *did:* Add the four sides and recognise the curl · *tex:* $\displaystyle\oint_{\rm rect}\mathbf u\cdot\mathbf t\,ds =
     \Big(\dfrac{\partial u_2}{\partial x_1} - \dfrac{\partial u_1}{\partial x_2}\Big)\Delta A + O(\Delta^4) = (\nabla\times
     \mathbf u)\cdot\mathbf n\,\Delta A$ · *why:* Steps 4 and 5 added, with $\Delta A = \Delta x_1\Delta x_2$; the bracket
     is $(\nabla\times\mathbf u)_3$ by (2.25) (D12), and n = e₃ picks that component. Dividing by ΔA and letting Δ → 0
     gives (2.35) — Ex. 2.6 in the x₃-plane. · *plain:* One small loop's circulation is the normal curl times its area.
  7. *did:* Tile the surface with small rectangles and add · *tex:* $\displaystyle\sum_k\oint_{C_k}\mathbf u\cdot\mathbf t\,ds
     = \sum_k(\nabla\times\mathbf u)\cdot\mathbf n\,\Delta A_k \to \iint_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA$ · *why:*
     Step 6 holds for every tile, all oriented counterclockwise about the same n; the right side is a Riemann sum of a
     continuous function, which tends to the surface integral as the tiles shrink. · *plain:* The sum of all the little
     circulations is the total curl inside.
  8. *did:* Cancel the interior edges · *tex:* $\displaystyle\sum_k\oint_{C_k}\mathbf u\cdot\mathbf t\,ds = \oint_{\text{outer
     edges}}\mathbf u\cdot\mathbf t\,ds$ · *why:* An edge shared by two neighbouring tiles is traversed once in each
     direction (both tiles anticlockwise), so u·t flips sign while u is the same: the two contributions cancel. Only
     the edges on the outside survive. · *plain:* Every inner edge is walked twice, opposite ways, and drops out.
  9. *did:* Refine: the outer edges become C · *tex:* $\displaystyle\oint_C\mathbf u\cdot\mathbf t\,ds = \iint_A(\nabla\times
     \mathbf u)\cdot\mathbf n\,dA$ · *why:* As the tiles shrink the staircase of outer edges tends to the boundary curve
     C with tangent t (u continuous, C piecewise smooth), so the left side of step 8 tends to the circulation round C.
     With step 7 this is (2.34). For a curved A: the same on each flat patch, with t = n_c × n at the boundary. ·
     *plain:* Stokes' theorem: circulation round the edge equals curl flux through the surface.
  10. *did:* Corollary for a gradient field · *tex:* $\displaystyle\oint_C\nabla\phi\cdot\mathbf t\,ds = \oint_C d\phi = 0\
      \Rightarrow\ \nabla\times\nabla\phi = 0$ · *why:* Along the curve, $\nabla\phi\cdot\mathbf t\,ds = d\phi$ (chain
      rule, P49): the integral of an exact change round a closed loop is φ(end) − φ(start) = 0. Then (2.34) forces
      $\iint(\nabla\times\nabla\phi)\cdot\mathbf n\,dA = 0$ for every A, so the integrand vanishes (Exercise 2.20). ·
      *plain:* A field that is the slope of something has no circulation and no curl.
- **Result.** $\displaystyle\iint_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA = \oint_C\mathbf u\cdot\mathbf t\,ds$ (2.34), and
  $\nabla\times\nabla\phi = 0$ — *in words:* the circulation round a loop is the total curl inside it; gradient fields
  have none.
- **Check.** Units: (1/s)·m² = (m/s)·m = m²/s ✓. Limit: u uniform: both sides 0 ✓. Number: u = b × x, b = e₃, unit disc:
  circulation 2πR = 6.283, curl 2 × area π = 6.283 ✓ (`stokes_theorem_check`); shear u₁ = Γx₂, unit square: −Γ both
  sides ✓; u = ∇(x² − y²): circulation 1e-16 ✓; the irrotational vortex with the core inside: 2πK vs "curl flux 0" —
  the hypothesis (u differentiable on A) fails, and the code says so.
- **What it means.** Circulation is "total spin inside", not "curved streamlines": a straight shear flow has it, a
  whirlpool with a singular core has it only through the core. Kelvin's theorem (Ch. 5) follows the circulation of a
  material loop; the lift of a wing (Ch. 6, 14) is proportional to it; a vortex tube's strength is the same through
  every cross-section because ∇·(∇×u) = 0. The theorem fails when A contains a point where u is not differentiable
  (the vortex core) — then the loop's circulation is fixed by the core, whatever the loop.
- **Traps.** The orientation of the horizontal sides (the top runs in −e₁). Cancelling interior edges when tiles are
  oriented differently. Applying (2.34) across a singular core. Confusing Γ (circulation) with Γ (shear rate) — the
  code names differ.

---
*End of the design. Counts: CORE blocks 16 · IDs placed 94 · explainers 5 + backup · derivations 15 (steps: D01 6, D02
8, D03 4, D05 9, D06 8, D09 11, D12 5, D14 6, D15 7, D17 15 + sympy, D18 8, D21 7, D22 9, D25 9, D26 10 = 122) ·
primers 25 new + 20 ch01 reminders + 9 glosses · ledger rows 123.*
