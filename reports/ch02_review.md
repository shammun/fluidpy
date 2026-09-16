# ch02 — Cartesian Tensors · derivation review (phase 6)

Independent, read-only, line-by-line review of the coded chapter against the rendered book pages
(chapters/pages/ch02/p067–p089, plus ch03 p104–p106 for the rotation-tensor cross-check). Reviewer re-ran
`pytest tests/test_ch02.py -q` → 91 passed, 1 xfailed (strict, O1) before any fix.

Verdict at review time: **2 must-fix** (both real physics/labeling errors that the tests shared with the code),
7 should-fix, everything else verified. Resolution status is tracked at the bottom.

## Must fix

**M1 — `rotation_tensor` is half the book's rotation tensor** (`fluidpy/core/tensors.py` ≈ l.948).
The code returns the antisymmetric part A = ½(G − Gᵀ) and calls it R. The book (§2.10 p082; ch03 (3.15), (3.17) p106)
defines R_ij = ∂u_i/∂x_j − ∂u_j/∂x_i = −ε_ijk ω_k with ω = ∇×u the *vorticity* (no ½): ∂u_i/∂x_j = S_ij + ½R_ij,
and "ω and R represent twice the fluid element rotation rate". Proof: for u = b × x,
`vector_from_antisymmetric(rotation_tensor(G))` gave b = ½∇×u instead of 2b.
Fix: `rotation_tensor` returns G − Gᵀ; docstring restated; discriminating test (vector of R for b × x = 2b).
Upstream wording to correct: analysis/ch02.md F18, design Part C 2.13; the notebook must say
"A = ½R_book, vector(A) = ½ω_book = ½∇×u", never "R = A".

**M2 — `velocity_gradient_preset` mislabels its flows** (`fluidpy/ch02_cartesian_tensors.py` ≈ l.424–448).
Code: pure_strain = [[0,Γ],[Γ,0]], uniaxial_extension = diag(Γ,−Γ) (traceless, i.e. biaxial planar strain),
irrotational_strain = alias of pure_strain. Design Part C 6.6 and the E3 storyboard: pure_strain = diag(Γ,−Γ),
uniaxial_extension = diag(Γ,0) (∇·u = Γ ≠ 0), irrotational_strain = [[0,Γ],[Γ,0]]. A learner pressing
"uniaxial extension" in E3 would get a divergence-free flow. Fix: matrices renamed to the design's; tests updated.

## Should fix

- O2 (from verify, confirmed + extended): `boundary_tangent` and `planar_loop` docstrings and analysis §9 item 7 say
  n_c is "outward"; Fig. 2.10 draws n_c *into* A (with outward n_c, t = n_c × n would be clockwise). Code is right.
  Record "n_c points into A" in knowledge/notation.md.
- `scripts/ch02_fig2_2_rotated_axes.py` l.87, l.151 label the orthogonality C_ij C_ik = δ_jk as Eq. (2.6); (2.6) is
  x'_i = x_k C_ki; orthogonality is unnumbered (Exercise 2.8 / D02).
- `material_line_angle` docstring self-contradictory ("turns at ω₃/2 (= the vector of A)"); one consistent statement.
- `operators.divergence` raises for a 3-component z-independent field on a 2-D grid while `curl` and
  `vector_gradient` accept it; sum over the first nd components.
- `integral_gradient` returns [i,j] = ∂Q_j/∂x_i while `vector_gradient` returns G[i,j] = ∂u_i/∂x_j (both documented,
  the book itself is split); cross-reference "= transposed" in both docstrings.
- `traction_2d` docstring restricts σ_n = τ'₁₁, τ_s = τ'₁₂ to symmetric τ; it holds for any τ.
- O1 (from verify): `gauss_gradient_box` docstring promises a vector-Q (3,3) result but raises; strict xfail exists.

## Verified (page by page)

- (2.5), (2.7), (2.8) p069–p070: transform_vector = Cᵀx, inverse = Cx', direction_cosines row old / column new,
  rotation_matrix_2d = Ex. 2.1's C applied passively, rotation_matrix_3d Rodrigues, polar_components = Ex. 2.1.
- (2.12) p073: transform_tensor = C_im C_jn τ_ij = CᵀτC; example_2_2 τ'₁₁, τ'₁₂ match p077 line by line.
- Stress sign convention / Fig. 2.4 p073: STRESS_CUBE_FACES letters and ± face rows match; cube_face_tractions = traction(τ, ±e_i).
- (2.15) p075: traction = τ_ji n_j = n·τ, discrimination test with non-symmetric τ; double_dot "book" = A_ij B_ji vs "frobenius".
- (2.18)–(2.21) p078–p079: permutation_sign on all six permutations, epsilon_delta_residual pairing, cross and cross_einsum.
- (2.24)–(2.25) p081: curl einsum on G[k,j], curl_components, 2-component branch; Ex. 2.3 (a x₃e₃) symbolic.
- (2.26)–(2.29) p082: antisymmetric_from_vector = −ε_ijk ω_k, vector_from_antisymmetric = −½ε_ijk R_ij (i = 1),
  R·x = ω × x; symmetric_double_contraction as (2.28) writes.
- §2.11 / Ex. 2.4 p083–p084: principal_axes refuses non-symmetric input, det B = +1; example_2_4 (Γ, −Γ), 45°, S' = diag(Γ,−Γ),
  Γ ≡ S₁₂ with the book's factor-2 inconsistency documented.
- (2.30)–(2.33) p085–p086: gauss_gradient_box pairs n_i with the free index; _cube_face_sum weights; faces at x0 ± h/2.
- (2.34)–(2.35), Fig. 2.10, Ex. 2.6 p088–p089: _plane_basis e₁ × e₂ = n; loops counterclockwise about n
  (= n_c × n with n_c into A); integral_curl_component = circulation/h²; u_y correction discriminated (u = (0,z,0) → −1).
- Stencils / grid: partial axis = ndim − 1 − direction, interior central, one-sided second-order edges, second_partial
  edge stencil (2,−5,4,−1)/h², periodic grid drops the duplicate node; observed orders 1.97–2.03, Simpson 4.00.
- Labels: no unverified/qualitative labels in ch02 modules; V5 entries cited; discrimination proofs match the tests read.

## Resolution

(filled in by the orchestrator after the implementer / verifier round)
