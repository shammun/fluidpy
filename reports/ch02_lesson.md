# Chapter 2 — lesson review (round 1 → round 2 below)          2026-09-16

coverage_check: **OK** (0 errors, 28 warnings — all "ledger primer not found among the notebook's primers"; every one is
a *reminder* or *gloss* row of design Part E, checked by hand below: 26 are present as 🔁 reminders / glosses, 2 are
missing — f-strings and tuple unpacking — see Should fix) · sections covered **14/14** (§2.1, 2.2, 2.3, 2.4+2.6, 2.5,
2.7+2.8, 2.9, 2.10, 2.11, 2.12, 2.13, 2.14, plus S01/S02 pointer lines at cells 402–403) · CORE (A) blocks **16/16**
with code + visual (C01 fig 31 + plotly 35 · C02 plotly 64 + animation 68 · C03 fig 89 + explainer 93 · C04 plotly 111 ·
C05 plotly 135 + slider 139 + explainer 143 · C06 fig 157 + slider 159 · C07 fig 180 · C08 fig 206 · C09 fig 228 +
slider 230 + fig 233 · C10 fig 245 + animation 248 · C11 fig 265 + explainer 268 · C12 animation 290 + explainer 294 ·
C13 slider 315 + live 320 · C14 plotly 339 + slider 342 + explainer 346 · C15 animation 365 + slider 368 · C16 plotly
388 + animation 391 + fig 394) · derivations **15/15** present with the design's step counts (D02 8, D03 4, D01 6, D05 9,
D06 8, D18 8, D09 11, D12 5, D14 6, D15 7, D17 15 + sympy cell 308, D25 9, D21 7, D22 9, D26 10) · explainers **5/5**,
each embedded exactly once (cells 93, 143, 268, 294, 346; `stokes_circulation_loop` sits in C11 with a pointer back
from C16 at cell 396, as the design planned) · the executed notebook has **405 cells (96 code, 309 markdown)**, 0 error
outputs, 11 PNG figures, 12 plotly figures, 5 animations, 5 explainer frames. I read all 405 cells in order (dumped to
`outputs/ch02/lesson_scratch/dump.txt`), looked at all 11 PNGs and at the first and last frame of all 5 animations
(`outputs/ch02/lesson_scratch/anim_*.png`), and re-ran the questionable computations in a scratch shell.

```
$ .venv/Scripts/python.exe tools/coverage_check.py ch02 --nb outputs/ch02/executed.ipynb
COVERAGE OK for ch02 (0 errors, 28 warnings)
```

Spot-checked numbers (text vs printed output): Ex. 2.1 u_r = 1.866, u_θ = 1.232 (cells 82–83) ✓ · Ex. 2.2 f = (0.5,
0.866) Pa, |f| = 1, 60°/240°, τ' = [[0.866, 0.5], [0.5, −0.866]] (cells 128–132, 153–154) ✓ · Ex. 2.3 3a and 2b on the
grid (cells 242, 262) ✓ · Ex. 2.4 λ = ±Γ, 45°, S' = diag(1, −1) (cells 311–312) ✓ · invariants (9, 24, 18) and roots
1.268/3/4.732 (cells 176–177, 312) ✓ · sphere 8π/3 = 8.3776 both sides (cell 336) ✓ · tiled interior faces 0.0 (cell
336) ✓ · circulation 2π = curl flux 2π, shear −Γ, ∇×∇φ 5.9e-13 < 1e-11 (cell 385) ✓. Conventions stated as required:
C_ij = e_i·e'_j passive with x' = Cᵀx (cells 42, 51, 63, 79) ✓ · f_i = τ_ji n_j contracts the first index (cells 124,
127, 134) ✓ · the book's R = G − Gᵀ carries ω = ∇×u while A = ½R carries ½∇×u, never called ω (cells 273, 280–283,
286–288, 404) ✓ · Stokes n_c into A, t = n_c × n counterclockwise about n, (n_c, n, t) right-handed (cells 378, 379,
389) ✓ · book typos taught corrected: Ex. 2.4 Γ ≡ S₁₂ (cell 310), Ex. 2.3 a x₃e₃ (cell 239), Ex. 2.6 u_y (cell 381),
(2.27) lower limit i = 1 (cells 281, 282 step 5) ✓ · orthogonality cited as "Exercise 2.8 / D02", never as Eq. (2.6)
(cells 52, 54, 60) ✓.

## Must fix (numbered: cell number · what · the concrete change)

1. **Cells 242 and 262 · the "sympy twin" prints the wrong answer and the text claims the right one.** Cell 242's last
   line prints `(0, [0, 0, 0])` while its comment says "(3a, [0, 0, 0])" and cell 243 item 6 says "`exact_div_curl`
   computes the same symbolically: (3a, 0)"; cell 262 prints `(0, [0, 0, 0])` while its comment says "(0, [0, 0, 2])"
   and cell 263 item 4 says "The sympy twin gives (0, [0, 0, 2])". Cause (reproduced in a scratch shell):
   `ch02.coordinates(3)` creates `Symbol('x1', real=True)` etc., so `sp.symbols('x1 x2 x3')` (no assumptions) are
   *different* symbols and every derivative is 0. Fix in `notebooks/build_ch02.py`: build the field from the library's
   coordinates — `X3 = sp.Matrix(ch02.coordinates(3))` in cell 242 (cell 262 then inherits it) — and ideally make
   `core/fields.VectorField` accept a `coords=` argument or default to the expressions' free symbols so this cannot
   silently happen again (a one-line `assert div != 0` in the builder would have caught it). Re-execute and check the
   printed lines read `(3*a, [0, 0, 0])` and `(0, [0, 0, 2])`.
2. **Cell 189 · primer P73 states the `np.transpose` semantics wrongly.** It says `np.transpose(eps, (1, 2, 0))` gives
   `new[i, j, k] = eps[j, k, i]`. With numpy's convention (new axis 0 = old axis 1, …) the correct statement is
   `new[i, j, k] = eps[k, i, j]` (verified with `np.arange(27).reshape(3,3,3)`: the primer's formula fails, the
   corrected one holds). The demo in cell 190 passes only because ε is invariant under *both* cyclic shifts, so the
   error is invisible there and would bite the first time a reader moves an index of a non-cyclic array (N28's
   fourth-order tensors). Fix: correct the sentence to "`np.transpose(a, axes)` puts old axis `axes[p]` in position p,
   so `np.transpose(eps, (1, 2, 0))[i, j, k] = eps[k, i, j]`", and make the demo honest by also showing it on a
   non-symmetric array (e.g. `a = np.arange(27).reshape(3,3,3); np.transpose(a,(1,2,0))[0,1,2] == a[2,0,1]`).
3. **Cell 391 (animation, frames player) · the bookkeeping text is clipped at the right edge of the figure.** The line
   reads "side h = 2.00 m: Γ_circ = ∮u·t ds = −4.0000 m²/s, A = 4.0000 m², Γ_circ/A =" and stops — the ratio
   "−1.000 = (∇×u)₃ = −Γ", which is the whole point of the animation and which cell 393 tells the reader to watch
   ("the ratio stays at −1.000"), is cut off in every frame. Also the grey shear-profile arrows at x₁ ≈ −1.2 overlap the
   start of the text. Fix: split the readout into two lines (or use `fig.text` below the axes as cells 248/290 do), drop
   the "m²/s"/"m²" units into the axis title, and move the profile quiver to x₁ = −1.3 or below the text's y.

## Should fix

1. **Cell 128 · `arctan2` used one cell before its primer.** The tiny example's step 5 writes "θ = arctan2(0.866, 0.5) =
   60°" and P70 (cell 129) explains `arctan2` afterwards. Swap cells 129–130 before 128.
2. **Cells 159–161 · I₁ and "Mohr's circle" used before they are defined.** Cell 159's comment says "centre I₁/2 = 0"
   and cell 160 "centre I₁/2 (here 0)"; I₁ is defined in C07 (cell 169). Write "centre (τ₁₁ + τ₂₂)/2 — the trace/2, an
   invariant (C07)". Mohr's circle itself is only ever described as "the circle every such point lies on"; add one
   sentence at cell 160 saying what it is (the set of all (σ_n, τ_s) pairs over all planes is a circle of centre
   (τ₁₁ + τ₂₂)/2 and radius √(((τ₁₁ − τ₂₂)/2)² + τ₁₂²), Ch. 4) — the reader is told to "find the shear-free planes" with it.
3. **Cell 170 (D18) step 8 and cell 169 item 3 · "eigenvalues", "roots", "principal frame A' = diag(λ)" used two blocks
   before the eigenvalue primer P80 (cell 300).** The forward pointer is there ("D17 fact 3, forward pointer"), but a
   one-clause gloss would make D18 self-contained: "(an eigenvalue λ and eigenvector b satisfy A·b = λb — a direction A
   only stretches; for symmetric A three perpendicular such directions exist and, taken as axes, make A diagonal — C13)".
4. **Cell 52 (D02) step 7 and cell 170 (D18) step 5 · det(AB) = det A · det B is cited as "Ch. 1 primer P53", but P53
   (knowledge/primers.md) covered "determinant zero ⇔ dependent columns; minors; 3×3 by cofactors" only.** Add a
   one-line gloss the first time (D02 step 7): "the determinant is the volume-scale factor of a linear map, and scale
   factors multiply — so det(AB) = det A det B and det Cᵀ = det C".
5. **Cell 307 (D17) step 15 · the variance bound "a variance of numbers confined to an interval of length L is at most
   (L/2)²" is asserted, not justified.** One line does it: Var ≤ E[(X − c)²] for any constant c; take c the midpoint of
   [λ_min, λ_max], then every |X − c| ≤ L/2. (Popoviciu's inequality; the reader is stats-strong, so name it.)
6. **Cell 266 and cell 267 · "the irrotational vortex u_θ = K/r" is used before it is introduced** (it is defined only
   in cell 385's code comment). Add a one-sentence gloss at cell 266: "a flow circling the origin with speed K/r —
   curl-free everywhere except at r = 0, where the speed blows up".
7. **Cell 60 · "Rodrigues' formula" is named in a comment without a gloss** (ledger says "code comment"). Add "(the
   closed form for a rotation by θ about a unit axis: cos θ I + sin θ [k×] + (1 − cos θ) k kᵀ)" or point to the
   docstring.
8. **Cell 328 · comment vs printed number.** The comment says "(0.3329 with n = 20)"; the cell prints 0.333125.
   That value is right (the midpoint sum of x² over 20 cells is 1/3 − (b − a)h²·f''/24 = 1/3 − 1/4800 = 0.333125);
   change the comment to "0.3331 with n = 20 — the midpoint rule undershoots a convex integrand by ∝ 1/n²".
9. **Cell 336 · comment "(quadrature error ~1e-4)" for the sphere benchmark, but the printed sides agree with 8π/3 to
   ~1e-14** (8.377580409572767 / 8.377580409572744 / 8.377580409572781). Cell 337 item 3 repeats "to quadrature
   error". Say what is true: the spherical product rule integrates these low-degree polynomials essentially exactly.
10. **Cell 55 / cell 56 (D03) · D03 starts from (2.5) before D01 derives it (D01 is in the next block).** The notebook
    says so honestly ("here we only need that it holds"), but a first-time reader has not yet seen what (2.5) says.
    Either move D03 to C03 right after D01 (the curation places D03 in C02, so alternatively) add one sentence before
    cell 56 stating (2.5) in words: "each new component is the old components weighted by the cosines of one column of
    C — the projection you saw in Fig. 2.2".
11. **Derivation titles · "(Eq. Exercise 2.8)", "(Eq. Exercise 2.9)", "(Eq. §2.10)", "(Eq. §2.11 facts (1)–(4))",
    "(Eq. 2.32 → 2.23)" (cells 52, 170, 278, 307, 359).** `nb.derivation` prefixes "Eq." to every label; let the
    builder pass a label without the prefix (or make nbkit skip "Eq." when the label does not start with a digit).
12. **Missing ch01 reminders planned in design Part E:** f-strings (P04; first used in cell 3/28 — a one-line
    "(Ch. 1 P04)" comment) and tuple unpacking (P14; first used at cell 50 `Q, _ = np.linalg.qr(...)`, the reminder
    appears only at cell 301). Also `np.linspace` is first used at cell 68, its reminder comes at cell 138; `np.geomspace`
    (cell 391) and `ax.imshow` (cell 206, "heatmap" in 245) have no gloss — one clause each in the comment.
13. **Cell 365 (animation) · the last frames' text line runs past the axes frame** ("… (1/V)∮ → ∂Q₁/∂x₁ = 2x₀ = 2.000
    Δx = 0.130: (1/V)∮ = 2.000") — legible but untidy; wrap onto two lines.
14. **Cell 216 (P76) · "x on the last axis … built with `np.meshgrid(z, y, x, indexing='ij')`" and cell 217 prints
    `g3.h` as "(hx, hy, hz)"** — fine, but say explicitly in the primer that `Grid.h` is ordered (x, y, z) while arrays are
    [k, j, i] = (z, y, x); readers will mix the two orders (cells 227, 244 rely on it).
15. **Cell 354 (D21) "what it means" · "a delta function, not a derivative"** — jargon for a first-time reader; either
    gloss ("an infinitely concentrated source whose total is finite") or drop it.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / symbol / tool | First cell | Explained at |
|---|---|---|
| index letters i, j, k ∈ {1,2,3}, ≡, boldface vs components, λ^k as label | 6 | 6 (notation cell) ✓ |
| scalar / vector / second-order tensor (orders 0, 1, 2) | 9–12 | 9–12 (N01, N02, N06, N24) ✓ |
| summation convention, free vs dummy index | 13–15 | 14–15 (C01) ✓ |
| position vector (2.1), column form, transpose ᵀ, `x.T` | 16 | 16 (N03, N04) ✓ |
| `np.einsum` index strings | 18 | 17 (P62) ✓ |
| `@`, `.T`, `np.eye`, matrix product by hand | 20 | 19 (P63) ✓ |
| `np.arange(...).reshape(3, 3)` | 18 | comment in 18 ("a 3×3 matrix 0…8"); numpy arrays Ch. 1 P03 ✓ |
| Kronecker δ (2.16), substitution rule (2.17), δ_ii = 3 | 21 | 23–24 (N40, N41) ✓ |
| `assert np.allclose` | 30 | 30 comment "(Ch. 1 primer P15)" ✓ |
| f-strings | 3, 28 | **missing reminder** (Ch. 1 P04) — Should fix 12 |
| `go.Cone`, `go.Scatter3d` lines, `go.Mesh3d`, `plotly_arrow` | 35 | 33 (P64) ✓ |
| rotated frame, e'_j, x'_j (2.3) | 40–44 | 42–44 (C02 idea, N09, N10) ✓ |
| orthonormal basis, projection, completeness | 46, D02 | 45 (P65) ✓ |
| `rotation_matrix_3d` (Rodrigues) | 46 | 46 comment; 60–61 name it "Rodrigues' formula" — gloss missing (Should fix 7) |
| `np.deg2rad`, cos(π/2 − θ) = sin θ | 46 | 47 (P66) — one cell late, trivial ✓ |
| `np.linalg.norm`, `np.linalg.qr`, random rotation | 50 | 49 (P67) ✓ |
| `np.random.default_rng` | 50 | 50 comment "(Ch. 1 primer P10)" ✓ |
| tuple unpacking `Q, _ = …` | 50 | 301 "(Ch. 1 P14)" — late; Should fix 12 |
| det(AB) = det A det B, det Cᵀ = det C | 52 (D02 step 7) | cited as P53, which did not cover it — Should fix 4 |
| continuity argument for det C = +1 | 52 (D02 step 8) | gloss in step 8 ✓ |
| bilinearity of the dot product | 52 (D02 step 3), 76 | gloss in D02 step 3 and D01 step 2 ✓ |
| (2.5) forward rule used by D03 | 55–56 | stated at 55, derived at 76 (D01) — Should fix 10 |
| passive vs active rotation | 63 | 63 (⚠️ callout) ✓ |
| `animate`, `show_animation(player="frames")`, `fig.set_layout_engine("none")` | 68 | 67 reminder (P16); layout comment in 68 ✓ |
| `np.linspace` | 68 | 138 reminder (P06) — late but a Ch. 1 primer; Should fix 12 |
| projection (2.4), formal definition (2.8), residual test | 75, 80 | 75 (N11), 80 (N18) ✓ |
| `dict` of `lambda` candidates, `vector_params` | 87 | 86 reminder (P23, P29); `vector_params` explained in 88 ✓ |
| log axis / `ax.set_yscale("log")`, round-off floor | 89 | comments in 89; log axes Ch. 1 P13 ✓ |
| `show_viz` | 93 | 91 reminder (P18) ✓ |
| matrix forms x' = Cᵀx, x = Cx' (§2.3) | 79, 96 | 79 (N16), 96 ✓ |
| stress τ_ij, face/force indices, sign convention | 100–105 | 101–105 (C04) ✓; stress P05 reminder at 101 ✓ |
| Newton III at a point (opposite faces) | 102 | 103 (R01 recap) ✓ |
| pressure as −pδ_ij | 105 | 105 item 4 (forward to Ch. 4) ✓ |
| `cube_face_tractions`, `stress_component_meaning` | 108 | comments + 109 ✓ |
| free-body diagram / Newton II | 117, 124 | 117 reminder (P09) ✓ |
| orders of smallness h² vs h³ | 118 | 118 (P68) ✓ |
| vector area of a closed surface, dA_i = n_i dA | 120 | 120 (P69), 123 (N35) ✓ |
| `np.arctan2` | 128 (tiny example) | 129 (P70) — one cell late; Should fix 1 |
| σ_n, τ_s split; double-angle forms | 128 | 128 step 6 (stated, Mohr preview) ✓ |
| `slider_figure`, `np.linspace` | 139 | 138 reminder (P17, P06) ✓ |
| `fig.update_yaxes(scaleanchor=…)` | 139 | comment "equal axes" ✓ |
| tensor transformation (2.12), tensor ≠ matrix, fourth order (2.13), outer product | 145–152 | 147 (idea), 148 (D06), 150–152 (N27–N29) ✓ |
| "true for every n' ⇒ coefficients agree" | 148 (D06 step 7) | gloss in step 7 ✓ |
| Mohr's circle, I₁ | 159–161 | I₁ defined at 169 (C07); Mohr's circle only described — Should fix 2 |
| contraction, trace, invariants I₁ I₂ I₃, closed chains | 164–169 | 166, 169 (C07) ✓ |
| Vieta's formulas, `sp.expand` | 168 | 167 (P71); sympy Ch. 1 P40 ✓ |
| characteristic polynomial, cofactor expansion, principal minors | 170 (D18 steps 6–7) | steps 6–7 (P53 for cofactors) ✓ |
| eigenvalues, roots, principal frame diag(λ) | 169 item 3, 170 step 8 | 300 (P80) — forward pointer only; Should fix 3 |
| (2.14) contractions as matrix products; A·u vs Aᵀ·u; double dot conventions | 173–175 | 173–175 (N31–N33) ✓ |
| `np.multiply.outer` | 172 | 172 (N30) ✓ |
| ε_ijk (2.18), index moves, cyclic order, parity, `itertools.permutations` | 186–188 | 186 (idea), 187 (P72) ✓ |
| three-axis arrays, `np.transpose(a, axes)` | 190 | 189 (P73) — **stated wrongly**; Must fix 2 |
| right-hand rule, orientation | 191–194 | 191 (P74) ✓ |
| cross product (2.20), determinant form, index form (2.21) | 193–196 | 194–196 (N46–N50) ✓ |
| epsilon–delta (2.19), 2δ, 6, triple product | 193, 199 | 200 (D09) ✓ |
| isotropic tensor, pseudotensor | 198 | 198 (N42) ✓ |
| `np.count_nonzero`, `np.array_equal`, `itertools.product` | 190, 203, 205 | comments ✓ |
| `ListedColormap`, `ax.imshow` | 206 | comments ("a colormap from a list of colours"; imshow only "heatmap" at 245) — gloss missing (Should fix 12) |
| ∇ (2.22), partial derivative, finite differences | 210 | 210 (R02 recap of P21, P22, P25) ✓ |
| level sets, directional derivative, chain rule | 213–214 | 214 (P75, P49 reminder) ✓ |
| `np.meshgrid`, project grid layout [k, j, i], `np.stack` on axis 0, `grid`, `grid2d` | 216–217 | 216 (P76) ✓ (ordering of `Grid.h` — Should fix 14) |
| numpy broadcasting, neighbour slices | 218–219 | 218 (P77) ✓ |
| `plt.contour`, `quiver`, `streamplot`, `ax.clabel`, `ax.annotate` | 221, 228 | 220 (P78); `clabel`/`annotate` commented ✓ |
| one-sided second-order stencil | 222 item 4 | 222 item 4 + from-scratch 227 ✓ |
| observed order of convergence, `tools.convergence.observed_order` | 225 | 224 gloss (P13 reminder) ✓ |
| `operator_convergence` "sin/cos field" | 225 | 226 item 5 ✓ |
| divergence (2.23), velocity gradient G_ij, divergence of a tensor (second index) | 236–240 | 239–240 (C10, N52 ⚠️) ✓ |
| ∂x_i/∂x_j = δ_ij, product rule | 239 | gloss in 239 ✓ |
| `radial_field`, `solid_body_rotation_field`, `VectorField`, `exact_div_curl`, `coordinates` | 242, 262 | comments + 243, 263 ✓ (but see Must fix 1) |
| `np.einsum('ii...', G)` (ellipsis) | 242 | comment "trace of G" only — minor |
| `scipy.linalg.expm`, `linear_flow_map` | 248 | 247 forward pointer; 274 (P79) — planned in the ledger ✓ |
| curl (2.24), components (2.25), operator ordering, paddle wheel | 252–256 | 254–256 (C11, D12) ✓ |
| solenoidal / irrotational | 242 (`is_solenoidal`), 259 | 259 (N54) — used in code one block earlier with a comment "solenoidal" ✓ |
| irrotational vortex u_θ = K/r, K | 266–267 | 385 code comment — Should fix 6 |
| `plt.Circle`, `connectionstyle="arc3"` | 265 | comments ✓ |
| symmetric / antisymmetric, S and A, the book's R vs A | 271–277 | 273, 276–277 (C12, N56) ✓ |
| (2.26)–(2.27), R·x = ω × x | 276, 280–282 | 280–282 (N57, N58, D15) ✓ |
| (2.28)–(2.29), τ:A = 0, even × odd analogy | 276 item 4, 284–285 | 284–285 (N59–N61) ✓ |
| `velocity_gradient_preset`, `rotation_tensor`, `deform_square` | 287, 290 | comments ✓ |
| eigenvalues/eigenvectors, `np.linalg.eigh`, `np.roots` | 299–301 | 300 (P80) ✓ |
| complex conjugate, |z|² | 302 | 302 (P81; P45 reminder) ✓ |
| quadratic form, Rayleigh quotient, Gram–Schmidt | 304 | 304 (P82, gloss) ✓ |
| variance bound (Popoviciu) | 307 step 15 | asserted only — Should fix 5 |
| principal angle ½ arctan2(2S₁₂, S₁₁ − S₂₂), principal values | 306 item 4 | 306 (stated) ✓ |
| quadratic formula, `np.c_`, `try/except ValueError` | 312–314 | comments ✓ |
| `live(...)` | 320 | 318 reminder (P47) + 319 note ✓ |
| Gauss (2.30), divergence theorem, outflux reading | 324–334 | 326, 331, 334 (C14, N64) ✓ |
| volume/surface integrals as midpoint sums, iterated integrals | 327 | 327 (P83) ✓ |
| fundamental theorem of calculus, `np.trapezoid` | 329–330 | 329 (P84; P27/P37 reminders) ✓ |
| tiling, cancellation of interior faces | 326, 332 | 332 (D25 steps 7–8) + 336 `divergence_theorem_tiled` ✓ |
| `flux_through_faces`, `gauss_gradient_box`, `smooth_test_field`, `volume_integral_box` | 336, 342 | comments ✓ |
| mean-value theorem for integrals | 351 | 351 (P85) ✓ |
| (2.31)–(2.33) generalised derivative, integral curl | 353–357 | 354 (D21), 356–357 (N66, N67, N72) ✓ |
| Ex. 2.5 face letters EADH/FBCG… | 359 | 359 start line ✓ |
| Riemann sum | 382 step 7 | not glossed; P83 covers the idea — minor |
| line integral round a loop, parametrising a circle/rectangle | 376–377 | 376 (P86; P35 reminder) ✓ |
| circulation, Stokes (2.34), (2.35), orientation n_c, t | 373–381 | 375, 378–381 (C16, N69–N73) ✓ |
| `planar_loop`, `planar_disc`, `stokes_theorem_check`, `hypothesis_ok`, `rectangle_loop` | 385, 391 | comments + 386 ✓ |
| `np.geomspace` | 391 | not glossed (P06 covered linspace/logspace) — minor, Should fix 12 |
| Green's theorem | 382 what-it-means | named as the planar case ✓ |
| comma notation (2.36), `comma_to_partial` | 399–400 | 399 (N74) ✓ |

Missing or late (the rows that are findings): f-strings (3) · tuple unpacking (50) · det(AB) rule (52) · (2.5) used by
D03 before D01 (56) · `np.linspace` (68) · `np.arctan2` (128) · I₁ / Mohr's circle (159) · eigenvalues / principal
frame in D18 (169–170) · `np.transpose` semantics wrong (189) · `ax.imshow` (206) · irrotational vortex (266) ·
variance bound (307) · `np.geomspace` (391).

## Derivation audit (D id · steps · every step follows? · why/in-words ok · check ok · verdict)

I worked every line by hand (and re-ran the numbers in the notebook's own outputs).

| D | Cell | Steps | Every step follows by its stated move? | *why* / *in words* | Check + meaning | Verdict |
|---|---|---|---|---|---|---|
| D02 orthogonality, det C = +1 | 52 | 8 | yes — completeness → dot → δ (rows), swapped bases (columns), det, continuity | yes; step 8's continuity gloss is good; step 7 cites P53 for det(AB) which P53 did not cover (Should fix 4) | 30° numbers, C = I limit, 200 random rotations (cell 60) ✓ | OK |
| D03 inverse (2.7) | 56 | 4 | yes — multiply by C_kj, D02, δ, rename | yes | round trip 1.866/1.232 → (1, 2) ✓ | OK (order vs D01: Should fix 10) |
| D01 transformation (2.5) | 76 | 6 | yes — both spellings, dot with e'_k, δ_jk, substitute, name C, rename | yes; bilinearity glossed | units, θ = 0, x = (1, 2) → (1.866, 1.232), lengths √5 ✓ | OK |
| D05 Cauchy (2.15) | 124 | 9 | yes — the −τ_ji on the −e_j faces, dA_j = n_j dA from P69, ρ(g − a)dV/dA ∝ h → 0 | yes; every assumption named where used | Pa, n = e₃, Ex. 2.2 (0.5, 0.866), τ = −pδ ✓ | OK |
| D06 tensor rule (2.12) | 148 | 8 | yes — f' = f_i C_in, insert (2.15), n_j = n'_m C_jm from (2.7), equate, "true for all n'", rename; τ' = CᵀτC | yes; passive/active trap called out | τ'₁₁ = 0.866, τ'₁₂ = 0.5, τ'₂₂ = −0.866 (cell 154) ✓ | OK |
| D18 invariants | 170 | 8 | yes — C-pair → δ for trace and chain; det; cofactor expansion; M = ½(I₁² − A_ijA_ji) verified by expanding; Vieta | yes; step 8 uses "principal frame" ahead of C13 (Should fix 3) | (9, 24, 18), roots 1.268/3/4.732, A = δ limit ✓ | OK |
| D09 ε–δ (2.19) + contractions + triple product | 200 | 11 | yes — i = j, l = m cases, unique k*, {l, m} = {i, j}, the two signs, 3δ − δ = 2δ, 2δ_rr = 6, two-place moves in the triple product | yes; δ_ii = 3 vs 1 trap called out | 81 cases residual 0, 2I and 6 printed, a × (b × c) number ✓ | OK |
| D12 curl components (2.25) | 256 | 5 | yes — nine terms, seven die, ε_132 = −1, cyclic relabelling justified by N43 | yes | 1/s, b × x → 2b₃, shear −Γ (cell 262) ✓ | OK |
| D14 split S + A, uniqueness, tensor | 278 | 6 | yes — add/subtract ½Bᵀ, symmetries, X = −X ⇒ 0, linearity of (2.12) with dummy swap | yes | B = [[0,2],[0,0]] → S, A; 6 + 3 = 9 ✓ | OK |
| D15 vector of R (2.26)–(2.27), R·x = ω × x | 282 | 7 | yes — antisymmetry, entries, contract with ε_ijl → −2ω_l, rename, one-place swap, (2.21) two-place form | yes; sign and misprint traps called out | ω = (0,0,−1), x = e₁ → (0,−1,0) both ways; round trip ✓ | OK |
| D17 ★★★ symmetric eigenproblem | 307 + 308 | 15 + sympy | yes — conjugate trick (symmetry used at 4 and 7, said), C from unit b's, CᵀτC = diag, n in eigenbasis, weighted-mean bounds, variance bound for shear | yes, except step 15's variance bound is stated without its one-line reason (Should fix 5) | sympy re-runs the construction for a symbolic 2×2 (discriminant ≥ 0, b¹·b² = 0, CᵀτC diagonal, σ_n = Σλc², Σc² = 1) ✓; Ex. 2.4 numbers ✓ | OK |
| D25 Gauss (2.30) | 332 | 9 | yes — iterated integral, FTC, two faces, n₁ = 0 on four faces, relabel, tile, interior cancel, refine | yes; hypotheses named | Q = (x, y, z) 3 = 3, sphere 8π/3, tiled interior 0.0 ✓ | OK |
| D21 (2.31)–(2.33) as limits | 354 | 7 | yes — MVT, ÷V, x* → x₀, contract, cross | yes | units, constant Q, (x², 0, 0) → 2.000, cubic order 2 ✓ | OK ("delta function" jargon — Should fix 15) |
| D22 Ex. 2.5 | 359 | 9 | yes — normals, Taylor, midpoint, front − back, constants cancel, other pairs, ÷V, e_i·∂Q/∂x_i = ∂Q_i/∂x_i | yes; the "face value = face average" trap is exactly right | (1.21 − 0.81)·0.04/0.008 = 2.000; x³ → 3 + h²/4 ✓ | OK |
| D26 Stokes planar + ∇×∇φ = 0 | 382 | 10 | yes — orientation of four sides, midpoint, bottom − top → −∂₂u₁ΔA, right − left → +∂₁u₂ΔA, tile, interior edges cancel, refine, dφ round a loop | yes; the top-side orientation trap called out | 2π = 2π, shear −Γ, ∇φ circulation 0, vortex core flags hypothesis_ok=False, ∇×∇φ 5.9e-13 ✓ | OK |

All 15 derivations: every line follows from the previous one by the stated move; assumptions are named at the step
where they enter; each has a units/limit/number check and a "what it means"; D17 has a sympy cell that rebuilds the
construction rather than only testing the result. Equation numbers and signs match the rendered pages as recorded in
the curation (I compared (2.5), (2.7), (2.12), (2.15), (2.19), (2.21), (2.25), (2.26)–(2.27), (2.30)–(2.35)).

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **The "one arrow, two shadows" thread** runs from Fig. 2.1 (cell 35) through the two-frame plotly (64), the frame
  animation (68, whose last frame literally shows x'₁ = 2, x'₂ = −1) to the residual bar chart (89): passive rotation
  becomes a picture before it is a formula. The passive/active ⚠️ callout (63) with Wikipedia/scipy named is exactly
  the "Which p_o?" pattern from ch01.
- **Deriving (2.15) before (2.12)** (cells 99, 147–148) and saying so turns the book's "Sommerfeld's tetrahedron" hand-wave
  into an eight-move proof; the −e_j face sign trap is called out three times (124, 125, 127) and the non-symmetric
  counter-example (134: n·τ ≠ τ·n) makes the index order a number, not a footnote.
- **Vector test as a residual table** (87–90): "a non-vector fails by order one, not by a rounding error" — a reusable
  device for every later "is this a tensor?" question (ch03 S, R; ch04 τ).
- **From-scratch loops next to every einsum** (30, 62, 85, 110, 134, 156, 179, 205, 227, 244, 264, 289, 314, 338, 364,
  387) with an assertion — the summation convention *is* the triple loop.
- **Book typos taught, not hidden**: Γ ≡ S₁₂ with the factor-2 argument (310), a x₃e₃ (239), u_y (381), i = 1 (281) —
  each in a ⚠️ or 📝 with the reason, never silently "fixed".
- **R vs A kept apart on purpose** (273, 283, 288, 404) with the numerical pin "vector(A) = ½∇×u on the grid" — this
  is the convention ch03 §3.4 needs; worth a line in knowledge/notation.md.
- **The integral definitions as experiments**: the face-bookkeeping frames player (365), the secant→tangent slider
  (368) and the two convergence plots (233, 394) make "(2.32) is a limit, approached at second order" measurable — the
  ch01 lesson "make an approximation measurable" applied again.
- **Stokes' hypothesis failure as a status, not an assertion** (`hypothesis_ok=False` with a note, 385–386): the code
  reports when a theorem does not apply instead of asserting through it. Candidate rule for `verify-implementation`.
- Teaching-style lesson candidates: (i) when a code cell has a "symbolic twin", assert it against the numeric result in
  the builder so a silent zero cannot pass (Must fix 1 would have been caught); (ii) view every animation's *last* frame,
  not only the first — clipped readouts (391) appear only once the numbers grow long; (iii) a Python idiom primer must be
  demoed on an object where a wrong statement would *fail* (189/190).

## Verdict (round 1): FAIL

Three Must-fix items, all small and local: the two `exact_div_curl` cells print `(0, [0, 0, 0])` while the text claims
`(3a, 0)` / `(0, [0, 0, 2])` (a real symbol-mismatch bug, not a display issue); primer P73 states `np.transpose` wrongly;
the Stokes shrinking-loop animation clips the very number its caption tells the reader to watch. Everything else —
coverage, all 16 CORE blocks, all 15 derivations, conventions, typo corrections, the spot-checked numbers, the figures —
passes. After the builder fixes cells 242/262, 189 and 391 (and ideally Should-fix 1–3), re-execute, re-run
`tools/coverage_check.py`, and this review can be flipped to PASS on a re-read of those cells.

---

## Round 2 (2026-09-16, after the notebook-builder's fixes)

Re-read from the fresh `outputs/ch02/executed.ipynb` (rebuilt and re-executed: **405 cells, 96 code / 309 markdown, 0
error outputs**; 11 PNG, 12 plotly, 5 animations, 5 explainer frames — the gauss_flux_box and stokes_circulation_loop
iframes are now the built explainers, 278 kB and 270 kB). `tools/coverage_check.py ch02 --nb outputs/ch02/executed.ipynb`
→ **COVERAGE OK (0 errors, 28 warnings)** — the same 28 name-matching warnings as round 1 (ledger reminders the tool cannot
match by name; all present). I diffed the round-1 and round-2 cell dumps (`outputs/ch02/lesson_scratch/round2.diff`,
170 changed lines): every change is one of the 18 requested fixes; nothing else changed, so no regression. I also
re-extracted first/last frames of all five animations myself and compared them with the builder's dumps
(`outputs/ch02/anim_stokes_loop_*_s.png`, `anim_ex25_box_*_s.png`) — identical content.

**Must fix — all three verified fixed**

1. Cells 242 / 262 · the sympy twins now use `X3 = sp.Matrix(ch02.coordinates(3))` and print **`3*a [0, 0, 0]`** and
   **`0 [0, 0, 2]`**, each followed by assertions that the symbolic value is non-zero and equals the grid value
   (`div.mean()` = 3; `w[:, 5, 5, 5]` = (0, 0, 2)). Cells 243 item 6 and 263 item 4 describe exactly what is printed. ✓
2. Cell 189 (P73) · now reads "`np.transpose(a, axes)` puts old axis `axes[p]` in position p, so
   `np.transpose(eps, (1, 2, 0))[i, j, k] = eps[k, i, j]`", and cell 190 demonstrates it on `np.arange(27).reshape(3,3,3)`:
   printed `19 19 15` = new[0,1,2], a3[2,0,1], a3[1,2,0] — the correct reading matches and the wrong one visibly does not. ✓
3. Cell 391 (Stokes shrinking-loop animation) · the readout is a two-line `fig.text` below the axes; units moved into
   the title; the profile quiver moved to x₁ = −1.3. First frame: "side h = 2.00 m: Γ_circ = ∮ u·t ds = −4.0000, A = h² =
   4.0000 / Γ_circ / A = −1.000 = (∇×u)₃ = −Γ"; last frame (h = 0.10): "−0.0100, 0.0100 / −1.000 = (∇×u)₃ = −Γ" — the
   ratio cell 393 tells the reader to watch is visible in every frame, nothing clipped. The `constrained_layout` switch is
   restored at the end of the cell; cell 394's figure is unchanged. ✓

**Should fix — 15/15 verified at the builder's cell map**

1 arctan2 primer P70 is now cell 128–129, before the tiny example (130) ✓ · 2 Mohr's circle: cell 159's comment gives
centre (τ₁₁ + τ₂₂)/2 and the radius formula; cell 160 item 2 defines Mohr's circle in a sentence and calls the centre
"half the trace, an invariant (C07)" ✓ · 3 eigen gloss at cell 169 item 3 (A·b = λb, principal frame) and D18 step 8's
*why* now names it ✓ · 4 det(AB) gloss in D02 step 7 ("volume-scale factor … scale factors multiply", P53 credited only
for the cofactor recipe) and D18 step 5 points back to it ✓ · 5 D17 step 15 now states Popoviciu's inequality with the
one-line reason (Var X ≤ E[(X − c)²], c the midpoint) ✓ · 6 irrotational-vortex gloss at cell 266 (u_θ = K/r, K in m²/s,
speed blowing up at r = 0) ✓ · 7 Rodrigues' formula glossed in cell 60's comment (cos θ I + sin θ [k×] + (1 − cos θ) k kᵀ)
✓ · 8 cell 328 comment "0.3331 … undershoots a convex integrand by ∝ 1/n²" ✓ · 9 cell 336 comment and cell 337 item 3
now say "agree to ~1e-14 — the spherical product rule integrates these low-degree polynomials essentially exactly" ✓ ·
10 cell 55 states (2.5) in words before D03 uses it ✓ · 11 derivation titles read "(Exercise 2.8)", "(Exercise 2.9)",
"(§2.10)", "(§2.11 facts (1)–(4))", "(2.32 → 2.23)" — no "Eq." prefix on non-numbered labels ✓ · 12 ch01 reminders:
f-string (cell 4), tuple unpacking (cell 50), np.linspace (cell 68), imshow gloss (cells 206, 245), np.geomspace (cell
391) ✓ · 13 cell 365's readout is now two lines top-left inside the axes, fully visible in the last frame ✓ · 14 P76 (cell
216) states the two orders: arrays [k, j, i] = (z, y, x) vs `grid.h` and component index (x, y, z) ✓ · 15 D21's "delta
function" glossed as "an infinitely concentrated source whose total is finite" ✓. Bonus: D26 step 7 glosses "Riemann sum";
cell 242 glosses the `'ii...'` einsum ellipsis.

**Regression skim.** All spot-checked numbers of round 1 unchanged (1.866; (0.5, 0.866); ±Γ at 45°; (9, 24, 18); 8π/3;
tiled interior 0.0; ∇×∇φ 5.9e-13). Cell 134's output is now split into two stream chunks ("[0. 1.7321]" / " [1. 0.]") —
same content, cosmetic. Cell indices of the CORE blocks, derivations and explainer embeds are unchanged.

**Remaining items:** none that block. (Optional, not required: the 28 coverage warnings are the tool matching ledger
reminder rows by exact name — a `reminder:` marker in `nb.md`/the ledger would silence them.)

## Final verdict: PASS
