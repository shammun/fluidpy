# Chapter 2 — Cartesian Tensors: curation
(from `analysis/ch02.md` — 94 inventory rows, 36 numbered equations (2.1)–(2.36), 26 derivations D01–D26, 26
implementation rows F1–F26; `book.yaml` ch02, `policy.tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`.
Curated 2026-09-16, concept-curator. Read: `knowledge/CUMULATIVE.md`, `concept_map.md`, `primers.md` (61 ch01 primers),
`notation.md`, `viz_patterns.md`, `analysis/ch01_curation.md` (format), the chapter text and rendered pages.)

Counts: A 16 · B 60 · C 18 · RECAP 2 · SKIP 2 · derivations written out 15 (★ 4 ★★ 10 ★★★ 1) · demoted to statements 11

Tier words (parsed by `tools/nbkit.py`): **CORE 16 · NOTE 74 · RECAP 2 · SKIP 2** = 94 rows. Depth: A 16 (all CORE) ·
B 60 (59 NOTE + R02) · C 18 (15 NOTE + R01 + S01 + S02). Every inventory row number `[#n]` appears exactly once in §2.

**Decisions that shape this chapter.**
1. **The chapter is a maths review, so "SEEN (pre-book)" is not "RECAP".** The analyst marked 20 rows SEEN; only the
   ones a *previous chapter of this project* taught become RECAPs (R01 Newton's third law on opposite faces, ch01 D05
   free-body reasoning; R02 the ∇ operator, ch01 primer P25). Pre-book knowledge (dot product, matrix product, identity
   matrix, cross product, eigenvalues, the divergence and Stokes theorems as statements) is tiered like any other row
   and its tools are primed (§4), because the promise is that nothing is used unexplained.
2. **Cauchy's traction (2.15) is taught before the tensor transformation rule (2.12)** (analysis §9): the book states
   (2.12) from a tetrahedron balance it never shows; our D06 derives (2.12) *from* (2.15) + (2.8) + (2.7). The notebook
   says so in one sentence. Teaching order §2.4 (stress tensor) → §2.6 (Cauchy) → §2.4 (2.12) → §2.5.
3. **Sixteen A items** (analyst's ★LB list, merged where two rows are one idea): the epsilon–delta identity (#53) is B
   inside the ε block but keeps its derivation D09 written out because the book never proves it; orthogonality of C
   (#18) is B inside the direction-cosine block for the same reason (D02, D03). The vector ↔ antisymmetric-tensor map
   (#71, #72) is B inside the decomposition block with D15 written out (the sign R·x = ω × x matters for ch03/ch05).
4. **Book typos are handled in our own words** (analysis §9): Ex. 2.4's Γ is taken as the off-diagonal element S₁₂
   (the notebook flags that for u₁(x₂) alone S₁₂ = ½ du₁/dx₂); Ex. 2.6's second bracket uses u_y; Ex. 2.3's third term is
   a x₃ e₃; (2.27)'s lower limit is i = 1; the "[?]" in Ex. 2.2 is ignored. Book-quoted numbers stay in the git-ignored
   `tests/book_values_ch02.json`.
5. **Conventions stated where first used** (⚠️ callouts, ch01 "Which p_o?" pattern): passive C (columns = new axes,
   x' = Cᵀx) vs Wikipedia/scipy active R; traction contracts the *first* index (n·τ) and equals τ·n only for symmetric
   τ (Ch. 4); tensor divergence contracts the *second* index; vector gradient G[i, j] = ∂u_i/∂x_j; double dot A:B =
   A_ij B_ji (book) vs A_ij B_ij (Frobenius); grid layout [k, j, i] = (z, y, x) for the whole project.
6. **Placement rule for sections without an A item** (§2.3, §2.8, §2.14): each still gets its own short notebook
   section whose B/C paragraphs are tagged with their A parent and say which A block develops them (§2.3 → C01/C03,
   §2.8 → C08, §2.14 → C10). Where the parent comes earlier, the paragraph sits inside the parent's block.

## 1. Teaching order (A IDs grouped by book section, B/C IDs under each; one sentence each: "once you see X, Y follows")
A items in **bold**; B and C items listed where they are taught, with depth in brackets.

**§2.1 Scalars, Vectors, Tensors, Notation**
- **C01 The summation convention** [#7]: once a repeated index means "add over 1, 2, 3", every formula of the book
  shrinks to one line and a machine (`expand_indices`) can write the hidden sums back out. Inside: N01 [B] scalar,
  N02 [B] vector, N03 [B] the basis form (2.1), N04 [C] column/transpose notation, N05 [B] our 3-D drawing of the
  position vector (Fig. 2.1), N06 [B] second-order tensor = one number per pair of directions, N07 [B] dot product
  (2.2), N08 [B] boldface vs index notation and "order = number of free indices", N12 [B] free vs dummy index (§2.2),
  N21 [B] matrix product as an index sum (2.9) (§2.3), N22 [C] the single-dot notation (2.10), N23 [B] the explicit
  3 × 3 product (2.11) with its nine sums printed, N24 [C] the 1 → 3 → 9 hierarchy (§2.4), N40 [B] the Kronecker delta
  (2.16) and N41 [B] its substitution rule (2.17) (§2.7, taught here because D02 needs δ), N45 [C] u·v = uv cos θ (§2.8).

**§2.2 Rotation of Axes: Formal Definition of a Vector**
- **C02 The direction-cosine matrix C_ij = e_i·e'_j** [#13]: once you see that the columns of C are the new axes
  written in old components, orthogonality C·Cᵀ = δ (N15, D02) and the inverse transformation (2.7) (N14, D03) follow.
  Inside: N09 [B] the rotated system with the same origin, N10 [B] the same vector in the primed basis (2.3), N17 [C]
  Fig. 2.2 (our drawing is N09's), N14 [B] (2.7), N15 [B] orthogonality and det C = +1.
- **C03 How components transform: x'_j = x_i C_ij (2.5) — and the formal definition of a vector (2.8)** [#14]: once the
  projection trick (2.4) is done once, every vector obeys the same rule, and any triple that does not is not a vector.
  Inside: N11 [B] the projection (2.4), N13 [C] (2.6) as (2.5) relabelled, N16 [B] the matrix forms x' = Cᵀx and x = Cx'
  (§2.3), N18 [B] the formal definition (2.8) with pass/fail examples, N19 [B] Ex. 2.1 polar components (worked
  number), N20 [C] Fig. 2.3.

**§2.4 Second-Order Tensors (part 1)**
- **C04 The stress tensor τ_ij and its sign convention** [#28]: once "first index = face normal, second = force
  direction" and "on the +e_i face positive components point along +e_j" are fixed, the cube of Fig. 2.4 reads itself
  and Ch. 4 can balance forces on it. Inside: N25 [B] our 3-D cube with nine arrows (Fig. 2.4), R01 [C] opposite faces
  carry equal and opposite stresses (Newton III, ch01), N26 [C] the 3 × 3 matrix of τ.

**§2.6 Force on a Surface** (taught here, before (2.12) — decision 2)
- **C05 Cauchy's traction formula f_i = τ_ji n_j (2.15)** [#43]: once the tetrahedron's volume terms shrink faster than
  its face terms, the force per area on *any* plane is τ contracted with n on the first index. Inside: N34 [B] the
  surface element as a vector dA = n dA, N35 [B] the tetrahedron balance and dA_i = n_i dA, N36 [B] our tetrahedron
  figure (Fig. 2.5), N37 [C] the analogy with u_n = u·n, N38 [B] Ex. 2.2 shear-flow traction (worked number), N39 [C]
  Fig. 2.6 (channel element; the E2 preset).

**§2.4 Second-Order Tensors (part 2)**
- **C06 The tensor transformation rule τ'_mn = C_im C_jn τ_ij (2.12)** [#32]: once f and n are both vectors, writing
  (2.15) in two frames forces (2.12) — one C per index — and that rule *is* the definition of a tensor. Inside: N27 [B]
  tensor ≠ matrix (a fixed array fails the test), N28 [B] fourth order (2.13), N29 [B] examples: u_i v_j and ∂u_i/∂x_j.

**§2.5 Contraction and Multiplication**
- **C07 Contraction: trace, invariants and the double dot** [#36]: once a contracted index pair is seen as "one C times
  its transpose = δ" (D18), every fully contracted quantity is frame-independent. Inside: N30 [C] multiplication raises
  the order, N31 [B] the four single contractions (2.14) as einsum patterns, N32 [B] A·u vs Aᵀ·u, N33 [B] the double dot
  in the book's convention (⚠️ callout).

**§2.7 Kronecker Delta and Alternating Tensor** (+ §2.8 Vector, Dot, and Cross Products)
- **C08 The alternating tensor ε_ijk** [#51]: once ε encodes "cyclic +1, anticyclic −1, repeat 0", the cross product,
  the curl, the rotation tensor and every vector identity become index sums, and the epsilon–delta identity (N44, D09)
  is the one tool that closes them. Inside: N42 [B] isotropic tensors (δ, ε up to a factor; ε only for proper
  rotations), N43 [B] index moves on ε, N44 [B] epsilon–delta (2.19) and its contractions, N46 [B] cross-product
  definition, N47 [B] components (2.20), N48 [C] determinant form, N49 [B] index form (2.21), N50 [C] the k = 1 check.

**§2.9 Gradient, Divergence, and Curl**
- **C09 The gradient ∇φ** [#61] (opens with R02, the ∇ recap): once ∇φ is the vector perpendicular to the level
  surfaces whose length is the steepest slope, ∂φ/∂n = ∇φ·n follows — and the grid + stencil tools used by every later
  field computation are introduced here. Inside: R02 [B] the del operator (2.22) in index form, N51 [B] Fig. 2.7
  (contours ⊥ ∇φ, direction n).
- **C10 The divergence ∇·u (2.23)** [#63]: once ∂u_i/∂x_i is read as the sum of the stretching rates along the axes,
  "solenoidal" has meaning and Ch. 4's continuity equation is one line away. Inside: N52 [B] the velocity gradient
  ∂u_i/∂x_j and the divergence of a tensor (⚠️ which index), N74 [B] comma notation (2.36) (§2.14).
- **C11 The curl ∇×u (2.24)–(2.25)** [#65]: once the curl is ε contracted with the velocity gradient, its three
  components (D12) and the solid-body-rotation result ∇×(b × x) = 2b follow. Inside: N53 [B] the components (2.25),
  N54 [B] solenoidal and irrotational fields, N55 [B] Ex. 2.3 (worked number for C10 and C11).

**§2.10 Symmetric and Antisymmetric Tensors**
- **C12 Every tensor is a symmetric plus an antisymmetric part** [#70]: once B = ½(B + Bᵀ) + ½(B − Bᵀ) is seen to be
  unique and frame-independent (D14), the antisymmetric part is a hidden vector (N57, N58, D15) and a symmetric tensor
  cannot see it (N59–N61). Inside: N56 [B] definitions and component counts, N57 [B] R from ω (2.26), N58 [B] the
  two-way map (2.27), N59 [B] (2.28), N60 [B] (2.29), N61 [C] τ_ij A_ij = 0 → ch04 dissipation.

**§2.11 Eigenvalues and Eigenvectors of a Symmetric Tensor**
- **C13 Principal axes of a symmetric tensor** [#76]: once τ·b = λb has real λ and orthogonal b (D17), rotating to the
  eigenvectors makes τ diagonal and the eigenvalues bound every normal stress. Inside: N62 [B] Ex. 2.4 principal axes
  of a plane shear strain rate (worked number, traced by hand), N63 [B] Fig. 2.8 (deforming square with the 45° axes).

**§2.12 Gauss' Theorem**
- **C14 Gauss' theorem (2.30)** [#79]: once the fundamental theorem of calculus is applied along one axis of a box
  (D25), "derivative inside = normal outside" holds for a field of any order, and the divergence theorem (N64) is the
  vector case. Inside: N64 [B] divergence theorem with the sphere benchmark 8π/3, N65 [B] Fig. 2.9 (our box with
  coloured faces).
- **C15 Divergence as outflux per unit volume (2.32)** [#83]: once (2.30) is applied to a shrinking box (D21), the
  Cartesian formula (2.23) comes back face by face (D22, Ex. 2.5) and the operators no longer depend on coordinates.
  Inside: N66 [B] the generalised derivative (2.31), N67 [B] integral curl (2.33), N68 [B] Ex. 2.5 (statement; its
  derivation is the block's D22), N72 [C] the coordinate-free remark (§2.13).

**§2.13 Stokes' Theorem**
- **C16 Stokes' theorem and circulation (2.34)** [#88]: once one small rectangle's circulation equals its curl times its
  area (D26), tiling any surface gives Stokes, and the normal curl is circulation per unit area (N71). Inside: N69 [B]
  the orientation rule t = n_c × n, N70 [B] Fig. 2.10, N71 [B] (2.35), N73 [B] Ex. 2.6 (with the u_y correction).

**§2.14 Comma Notation** — N74 [B → C10] in its own short notebook section (one paragraph, (2.36), pointer to §5.6).

**Exercises and literature** — S01, S02 (one pointer line each at the end of the notebook).

## 2. Chapter map (depth) — every inventory row exactly once
One row per inventory row (94, `[#n]` = analysis §2 row). `A parent` names the A block a B or C item is written in.

| ID | Item | § | Depth | Tier | A parent | Reason (A) / treatment (B) / pointer (C, SKIP) / source chapter (RECAP) |
|---|---|---|---|---|---|---|
| N01 | Scalar = zero-order tensor: one number, the same in every frame (p, T, ρ) [#1] | 2.1 | B | NOTE | C01 | stated in a paragraph: "zero free indices"; the ch01 fields as examples |
| N02 | Vector = first-order tensor: magnitude and direction, components change with the frame [#2] | 2.1 | B | NOTE | C01 | stated with the boldface / index pair u, u_i; the sharp definition is C03 |
| N03 | Position vector in the basis, x = e_i x_i (2.1) [#3] | 2.1 | B | NOTE | C01 | stated with Eq. (2.1) and one code line (`unit_vectors`, `vector_from_components`); x = (1, 2, 3) as the running example |
| N04 | Single free index i ∈ {1,2,3}; column form; transpose; triplet notation [#4] | 2.1 | C | NOTE | C01 | named in one sentence next to N03 (`x.T`); used in every code cell of the chapter |
| N05 | Fig. 2.1: position vector OP with its components and the unit vectors [#5] | 2.1 | B | NOTE | C01 | our own plotly 3-D figure of N03's vector: arrow, dashed projections, e₁ e₂ e₃ (rotatable) |
| N06 | Second-order tensor = one component per pair of directions, 3 × 3 = 9 [#6] | 2.1 | B | NOTE | C01 | stated in a paragraph: two free indices; the first physical example is C04, the definition is C06 |
| C01 | **The summation (Einstein) convention**: a repeated index is summed 1..3; free vs dummy indices [#7] | 2.1 | A | CORE | – | load-bearing: every index expression of the book (and of the ch03–ch05 notebooks) is read with it; `expand_indices` makes the hidden sums visible |
| N07 | Dot product as an implied sum, a·b = a_i b_i (2.2) [#8] | 2.1 | B | NOTE | C01 | stated with Eq. (2.2), (1,2,3)·(4,5,6) = 32 by hand, `dot` = `np.einsum('i,i')` |
| N08 | Boldface (meaning) vs indicial (manipulation); order = number of free indices; ab ≠ ba for tensors [#9] | 2.1 | B | NOTE | C01 | stated in a paragraph; `tensor_order("A_ij B_kl") == 4` as the one code line |
| N09 | Rotated system O1'2'3' with the same origin; unit vectors e'_j; same x, new components x'_j [#10] | 2.2 | B | NOTE | C02 | stated with our 3-D figure of two frames sharing an origin and one arrow with both component sets |
| N10 | The same vector in the primed basis, x = e'_j x'_j (2.3) [#11] | 2.2 | B | NOTE | C02 | stated with Eq. (2.3); `vector_from_components(xp, Ep)` reconstructs x — asserted equal to N03 |
| N11 | Projection of x on e'_1 gives x'_1 (2.4) — the first move of D01 [#12] | 2.2 | B | NOTE | C03 | stated with Eq. (2.4) and the number for θ = 30°, then generalised inside D01 |
| C02 | **The direction-cosine matrix** C_ij = e_i·e'_j: row = old axis, column = new axis; columns are the new axes [#13] | 2.2 | A | CORE | – | load-bearing: the matrix behind every rotation, principal-axis frame (ch03), stress in rotated frames (ch04) and rotating frame (ch13); orthogonality D02 and the inverse D03 live here |
| C03 | **Transformation of components** x'_j = x_i C_ij (2.5) — and the formal definition of a vector (2.8) [#14] | 2.2 | A | CORE | – | load-bearing: (2.5)/(2.8) is the definition of a vector; velocity, vorticity and ∇p pass this test in later chapters; D01 derived here |
| N12 | Free vs dummy index: the summed letter is arbitrary, the free letter must match on both sides [#15] | 2.2 | B | NOTE | C01 | stated with x_i C_ij = x_k C_kj; `classify_indices`, `rename_dummy` print the same expansion |
| N13 | (2.5) relabelled: x'_i = x_k C_ki (2.6) [#16] | 2.2 | C | NOTE | C03 | named; the demo expands (2.5) and (2.6) and asserts they coincide |
| N14 | Inverse transformation x_j = x'_i C_ji (2.7); the summed index of C changes slot [#17] | 2.2 | B | NOTE | C02 | stated with Eq. (2.7) and the round trip x → x' → x in code; derived in D03 (book never writes it) |
| N15 | Orthogonality of C: C_ij C_kj = C_ji C_jk = δ_ik, det C = +1 for a rotation (Exercise 2.8) [#18] | 2.2 | B | NOTE | C02 | stated with the identity and `is_orthogonal`, `is_proper_rotation` on `rotation_matrix_2d(θ)`; derived in D02 (book never writes it) |
| N16 | Matrix forms: x' = Cᵀ·x and x = C·x' [#19] | 2.3 | B | NOTE | C03 | stated with the matrix equation; `C.T @ x` asserted equal to `transform_vector` — ⚠️ passive C vs Wikipedia/scipy active R |
| N17 | Fig. 2.2: original and rotated axes, one vector [#20] | 2.2 | C | NOTE | C02 | named; drawn as N09's figure (our own) and animated in §6 |
| N18 | Formal definition of a Cartesian vector: any triple with u'_j = u_i C_ij (2.8) [#21] | 2.2 | B | NOTE | C03 | stated with Eq. (2.8) and `transforms_as_vector` on x, a x, b × x (pass) and (x₁², x₂², x₃²), (\|x\|, 0, 0) (fail) — a residual table |
| N19 | Ex. 2.1: Cartesian → polar components u_r, u_θ; C for the polar frame [#22] | 2.2 | B | NOTE | C03 | the block's worked number: u = (1, 2), θ = 30° traced by hand, then `polar_components`; C = [[cos, −sin],[sin, cos]] is (2.5) with j ∈ {r, θ} |
| N20 | Fig. 2.3: u resolved in (x₁, x₂) and (r, θ) [#23] | 2.2 | C | NOTE | C03 | named; our own 2-D drawing appears as the polar mode of the `rotation_of_axes` explainer |
| N21 | Matrix inner product P_ij = A_ik B_kj summed over the adjacent index (2.9) [#24] | 2.3 | B | NOTE | C01 | stated with Eq. (2.9); `inner` = `np.einsum('ik,kj->ij')` asserted equal to `A @ B` |
| N22 | Boldface form P = A·B; the single dot = one index summed (2.10) [#25] | 2.3 | C | NOTE | C01 | named next to N21; the dot-as-contraction reading recurs in C05 (n·τ) and C07 |
| N23 | Explicit 3 × 3 product with row/column boxes; P₁₂ = A₁₁B₁₂ + A₁₂B₂₂ + A₁₃B₃₂ (2.11) [#26] | 2.3 | B | NOTE | C01 | stated with `expand_indices("A_ik B_kj")` printing the nine sums and a row-i / column-j highlight figure (E5 pattern) |
| N24 | Hierarchy: scalar (1) → vector (3, one free index) → second-order tensor (9); orders 0–2 suffice for Newtonian fluids [#27] | 2.4 | C | NOTE | C01 | named; the fourth-order exception is N28 (ch04 §4.5 viscosity tensor) |
| C04 | **The stress tensor τ_ij**: first index = normal of the face, second = direction of the force; sign convention (tensile positive; on the +e_i face positive components point along +e_j) [#28] | 2.4 | A | CORE | – | load-bearing: the object Ch. 4 builds Cauchy's equation from; the sign convention decides every stress sign in ch04–ch09; ch01 C03 gave only the normal/shear split |
| N25 | Fig. 2.4: cube with the nine stress components on three faces [#29] | 2.4 | B | NOTE | C04 | the block's visualization: our plotly 3-D cube with nine arrows from `cube_face_tractions`, toggle for the hidden faces |
| R01 | Opposite faces carry equal and opposite stresses (Newton's third law at a point) [#30] | 2.4 | C | RECAP | C04 | reminded in one sentence: ch01 §1.7 C20 (D05 box faces, primer P09); `cube_face_tractions` sums to zero; the *difference* across a finite cube is Ch. 4's net force |
| N26 | The 3 × 3 matrix of τ (unnumbered) [#31] | 2.4 | C | NOTE | C04 | named; the array the code carries everywhere |
| C06 | **Transformation rule of a second-order tensor** τ'_mn = C_im C_jn τ_ij, τ' = Cᵀ·τ·C (2.12); a quantity obeying it *is* a tensor [#32] | 2.4 | A | CORE | – | load-bearing: defines "second-order tensor"; used for S and R (ch03), τ and σ (ch04), Reynolds stress (ch12); D06 derives it from (2.15) (book states it only) |
| N27 | Tensor ≠ matrix: an array's entries are tensor components only if they obey (2.12) [#33] | 2.4 | B | NOTE | C06 | stated with `transforms_as_tensor`: u_i v_j passes, a fixed array "the same in every frame" fails (residual printed) |
| N28 | Fourth-order transformation A'_mnpq = C_im C_jn C_kp C_lq A_ijkl (2.13) [#34] | 2.4 | B | NOTE | C06 | stated with Eq. (2.13) and one line `transform_tensor(A4, C)` (81 components); pointer: the viscosity tensor of ch04 §4.5 |
| N29 | Examples of second-order tensors: τ_ij, ∂u_i/∂x_j, the outer product u_i v_j (Exercise 2.10) [#35] | 2.4 | B | NOTE | C06 | stated; `outer(u, v)` transforms per (2.12) in one assertion; the velocity gradient is N52 |
| C07 | **Contraction**: A_jj = trace, a scalar invariant; the three invariants I₁, I₂, I₃ (Exercise 2.9) [#36] | 2.5 | A | CORE | – | load-bearing: ∇·u = S_ii (ch03), p = −τ_ii/3 (ch04), Reynolds-stress invariants (ch12), the characteristic polynomial of C13; the einsum patterns (N31–N33) are used in every later notebook; D18 (book gives hints only) |
| N30 | Multiplication raises order: P_ijkl = A_ij B_kl is fourth order [#37] | 2.5 | C | NOTE | C07 | named with `tensor_product` = `np.multiply.outer`; used with N28 in ch04 |
| N31 | The four single contractions of A_ij B_kl and their matrix forms (2.14) [#38] | 2.5 | B | NOTE | C07 | stated as a 4-row table "index string → matrix form" (`contract(A, B, "ij,ki->kj")` = B·A, …), each asserted |
| N32 | Tensor·vector: A_ij u_j = (A·u)_i vs A_ij u_i = (Aᵀ·u)_j [#39] | 2.5 | B | NOTE | C07 | stated with both forms and `dot_tensor_vector(A, u, index)`; the distinction is the one C05 needs (n·τ vs τ·n) |
| N33 | Double contraction: A:B = A_ij B_ji (book) vs A_ij B_ij (Frobenius) [#40] | 2.5 | B | NOTE | C07 | stated with both sums and a ⚠️ convention callout with numbers (they agree for symmetric operands); `double_dot(A, B, convention=)`; pointer: dissipation τ_ij ∂u_i/∂x_j in ch04 §4.5 |
| N34 | Surface element as a vector dA = n dA [#41] | 2.6 | B | NOTE | C05 | stated in a paragraph with a sketch; reused in C14 (n dA on a closed surface) |
| N35 | Tetrahedron force balance f₁ dA = τ₁₁ dA₁ + τ₂₁ dA₂ + τ₃₁ dA₃; geometry dA_i = n_i dA [#42] | 2.6 | B | NOTE | C05 | stated with the two relations and `tetrahedron_face_areas(n, dA)`; they are steps 3–6 of D05 |
| C05 | **Cauchy's traction formula** f_i = τ_ji n_j, f = n·τ (2.15) [#43] | 2.6 | A | CORE | – | load-bearing: every wall force, drag and lift (ch04 §4.3, ch08, ch09, ch14) is ∮ n·τ dA; the route to (2.12) (D06); D05 derived here |
| N36 | Fig. 2.5: tetrahedron with faces dA₁, dA₂, dA₃, slanted face dA, normal n, traction f [#44] | 2.6 | B | NOTE | C05 | the block's visualization: our plotly 3-D tetrahedron with stress arrows on the coordinate faces and f(n) on the slanted face |
| N37 | Analogy: f = n·τ is to τ what u_n = u·n is to u, but f is a vector [#45] | 2.6 | C | NOTE | C05 | named in one sentence after (2.15) |
| N38 | Ex. 2.2: shear flow τ = [[0, a],[a, 0]], n at 30°: f = (a/2, a√3/2), \|f\| = \|a\|, direction 60° (a > 0) or 240° (a < 0); τ'₁₁ = √3a/2, τ'₁₂ = a/2 via (2.12) [#46] | 2.6 | B | NOTE | C05 | the block's worked number, traced by hand with a = 1 Pa, then `example_2_2`; the general σ_n = a sin 2φ, τ_s = a cos 2φ stated (Mohr preview for C13); the (2.12) route is repeated inside C06 |
| N39 | Fig. 2.6: channel profile, element with normal at 30°, f and the rotated axes [#47] | 2.6 | C | NOTE | C05 | named; the channel case is a preset of the `cauchy_traction_principal_axes` explainer and the §6 slider figure |
| N40 | Kronecker delta δ_ij (2.16) = the identity matrix [#48] | 2.7 | B | NOTE | C01 | stated with Eq. (2.16) early in C01 (D02 needs it): `kronecker_delta()` = `np.eye(3)` |
| N41 | δ substitutes its summed index: δ_ij u_j = u_i (2.17) [#49] | 2.7 | B | NOTE | C01 | stated with Eq. (2.17) and the general rule "δ replaces its summed index" (δ_ij A_jk = A_ik, δ_ii = 3); `expand_indices("delta_ij u_j")` |
| N42 | Isotropic tensors: δ'_ij = δ_ij under any rotation (Exercise 2.11); δ the only one of order 2, ε of order 3 (up to a factor, proper rotations) [#50] | 2.7 | B | NOTE | C08 | stated with `is_isotropic` on δ (pass), ε (pass for det +1, sign flip under a reflection), a random tensor (fail); pointer: isotropic Newtonian stress model, ch04 §4.5 |
| C08 | **The alternating tensor** ε_ijk (2.18): +1 cyclic, −1 anticyclic, 0 repeated [#51] | 2.7 | A | CORE | – | load-bearing: cross product, curl, Coriolis 2Ω × u (ch04), vorticity (ch03, ch05) all run on ε; with (2.19) (N44, D09) it closes every vector identity |
| N43 | Index moves on ε: two places keep the sign, one place flips it [#52] | 2.7 | B | NOTE | C08 | stated with the two identities and `np.transpose(eps, …)` assertions |
| N44 | Epsilon–delta relation ε_ijk ε_klm = δ_il δ_jm − δ_im δ_jl (2.19); contractions 2δ_ij and 6; the vector triple product (Exercises 2.5, 2.7) [#53] | 2.7 | B | NOTE | C08 | stated with Eq. (2.19), `epsilon_delta_residual() == 0` (81 cases) and `triple_product`; **derived in D09** because the book only says "verify by choosing values"; used by D15, N55, ch05 |
| N45 | Dot product restated; u·v = uv cos θ; u·v = trace of u_i v_j (Exercises 2.12–2.13) [#54] | 2.8 | C | NOTE | C01 | named next to N07 with `angle_between`; the trace link points to C07 |
| N46 | Cross product w = u × v: magnitude uv sin θ, perpendicular, right-handed; u × v = −v × u; e₁ × e₂ = e₃ [#55] | 2.8 | B | NOTE | C08 | stated in a paragraph (pre-book reminder, right-hand rule primer) before the index form N49 |
| N47 | Cross product in components (2.20) (Exercise 2.14) [#56] | 2.8 | B | NOTE | C08 | stated with Eq. (2.20); `cross` (explicit) asserted equal to `np.cross`; D10 demoted (§4c) |
| N48 | Determinant form of u × v [#57] | 2.8 | C | NOTE | C08 | named; sympy `Matrix.det` reproduces (2.20) in one demo line |
| N49 | Index form (u × v)_k = ε_ijk u_i v_j = ε_kij u_i v_j (2.21) [#58] | 2.8 | B | NOTE | C08 | stated with Eq. (2.21); `cross_einsum` = `np.einsum('ijk,i,j->k')` asserted equal to N47; `expand_indices("eps_ijk u_i v_j")` prints (2.20); D11 demoted (§4c) |
| N50 | Check k = 1 of (2.21): only (2,3) and (3,2) survive [#59] | 2.8 | C | NOTE | C08 | named; the printed expansion of N49 shows it |
| R02 | The del operator ∇ = e_i ∂/∂x_i (2.22) [#60] | 2.9 | B | RECAP | C09 | reminded from ch01 primer P25 (∂, ∇ "points uphill") and P21/P22 (finite differences); new here: the index form and the grid + stencil tools `core.grids.grid`, `core.operators.partial` (📎 primer, §4) |
| C09 | **The gradient** (∇φ)_i = ∂φ/∂x_i: perpendicular to level surfaces, direction and size of the fastest change; ∂φ/∂n = ∇φ·n [#61] | 2.9 | A | CORE | – | load-bearing: ∇p (ch04), ∇φ potential flow (ch06), ∇T; the first field computed on the project grid — the stencil, layout and convergence test every later chapter reuses |
| N51 | Fig. 2.7: level curves of φ, ∇φ perpendicular to them, an arbitrary direction n [#62] | 2.9 | B | NOTE | C09 | the block's visualization: our contour + quiver figure with a slider for n (§6) |
| C10 | **The divergence** ∇·u = ∂u_i/∂x_i (2.23) [#63] | 2.9 | A | CORE | – | load-bearing: the most used operator in the book — continuity (ch04 §4.2), volume strain rate (ch03), incompressibility everywhere |
| N52 | Divergence of a tensor (∇·τ)_i = ∂τ_ij/∂x_j (contracts the second index); gradient of a vector ∂u_i/∂x_j raises the order [#64] | 2.9 | B | NOTE | C10 | stated with both formulas and a ⚠️ index-order callout (G[i, j] = ∂u_i/∂x_j fixed for ch03; (∇·τ)_i contracts the *second* index, (2.15) the first — consistent only for symmetric τ); `vector_gradient`, `tensor_divergence`; pointer ch03 §3.4, ch04 §4.4 |
| C11 | **The curl** (∇×u)_i = ε_ijk ∂u_k/∂x_j (2.24), components (2.25) [#65] | 2.9 | A | CORE | – | load-bearing: vorticity ω = ∇×u is the subject of ch05 and half of ch13; D12 derived here |
| N53 | The three curl components (2.25) [#66] | 2.9 | B | NOTE | C11 | stated as the result of D12; `curl_components` asserted equal to `curl` |
| N54 | Solenoidal (∇·u = 0) and irrotational (∇×u = 0) fields [#67] | 2.9 | B | NOTE | C11 | stated with both definitions, `is_solenoidal`, `is_irrotational`; pointer: incompressible flow ch04, potential flow ch06 (both at once), why "irrotational" in ch03 |
| N55 | Ex. 2.3: u = a x → ∇·u = 3a, ∇×u = 0; u = b × x → ∇·u = 0, ∇×u = 2b [#68] | 2.9 | B | NOTE | C11 | the worked number for C10 and C11: a = 1, b = e₃ by hand (index route in one line), then `radial_field`, `solid_body_rotation_field` on the grid vs sympy; the book's "a x₂ e₂" typo written as a x₃ e₃; pointer: vorticity = 2Ω (ch03); D13 demoted (§4c) |
| N56 | Symmetric (B_ij = B_ji, 6 independent) and antisymmetric (B_ij = −B_ji, zero diagonal, 3 independent) tensors [#69] | 2.10 | B | NOTE | C12 | stated with the definitions and `independent_components` (6 / 3 / 9) |
| C12 | **Decomposition** B_ij = ½(B_ij + B_ji) + ½(B_ij − B_ji) = S_ij + A_ij, unique [#70] | 2.10 | A | CORE | – | load-bearing: the spine of ch03 §3.4 (strain rate S + rotation R from ∂u_i/∂x_j) and of ch04 §4.5 (dissipation sees only S); D14 and D15 derived here |
| N57 | The antisymmetric tensor R associated with a vector ω (2.26) [#71] | 2.10 | B | NOTE | C12 | stated with the matrix (2.26) and `antisymmetric_from_vector(ω)`; pointer: R is the rotation tensor of the vorticity in ch03 |
| N58 | The two-way map R_ij = −ε_ijk ω_k, ω_k = −½ ε_ijk R_ij (2.27); R·x = ω × x [#72] | 2.10 | B | NOTE | C12 | stated with Eq. (2.27) (lower limit i = 1, the book's "i−1" is a misprint) and the round trip in code; **derived in D15** (inverse and the sign of R·x = ω × x are never derived by the book) |
| N59 | Double contraction of a symmetric τ with any B splits: τ_kl B_kl = τ_ij S_ij + τ_ij A_ij (2.28) [#73] | 2.10 | B | NOTE | C12 | stated with Eq. (2.28) and `symmetric_double_contraction(τ, B)` → (P, P_S, P_A) with a number |
| N60 | The same P with the A term negated (2.29) ⇒ τ_ij A_ij = 0 [#74] | 2.10 | B | NOTE | C12 | stated with Eq. (2.29) and the one-line lemma in words (X = −X ⇒ X = 0); P_A printed as 0 to machine precision; D16 demoted (§4c) |
| N61 | Result τ_ij B_ij = τ_ij S_ij = ½ τ_ij (B_ij + B_ji); the even × odd analogy [#75] | 2.10 | C | NOTE | C12 | named with the analogy in one sentence; pointer: viscous dissipation τ_ij ∂u_i/∂x_j = τ_ij S_ij in ch04 §4.5 |
| C13 | **Eigenvalues and eigenvectors of a real symmetric tensor**: real λ^k, orthogonal principal axes b^k, τ' = diag(λ) in the eigen-frame, λ bound the normal stresses [#76] | 2.11 | A | CORE | – | load-bearing: principal strain rates (ch03), principal stresses and Mohr's circle (ch04), normal modes (ch11) rest on the symmetric eigenproblem; D17 (★★★, book: "can be proved") derived here |
| N62 | Ex. 2.4: S = [[0, Γ],[Γ, 0]] for u = (u₁(x₂), 0); λ = ±Γ; b¹ = (1,1)/√2, b² = (−1,1)/√2; C = 45° rotation; S' = diag(Γ, −Γ) [#77] | 2.11 | B | NOTE | C13 | the block's worked number traced by hand (2 × 2 characteristic polynomial, normalisation, det C = +1 choice), then `example_2_4`; ⚠️ Γ is taken as S₁₂ (for u₁(x₂) alone S₁₂ = ½ du₁/dx₂ — the book's line "2S₁₂ = Γ" is inconsistent with its matrix); interpretation (stretch along b¹, squeeze along b²) previews ch03 §3.4; D19 demoted (§4c) |
| N63 | Fig. 2.8: original axes and the 45° eigen-axes [#78] | 2.11 | B | NOTE | C13 | the block's visualization: a square in the shear flow deforming, eigen-axes and Mohr's circle, slider Γ (§6) |
| C14 | **Gauss' theorem** ∭_V ∂Q/∂x_i dV = ∯_A n_i Q dA for a field of any order (2.30) [#79] | 2.12 | A | CORE | – | load-bearing: every conservation law of ch04 goes from integral to differential form through (2.30) (mass, momentum flux tensor, energy); D25 (book never proves it) derived here |
| N64 | Divergence theorem ∭ ∇·Q dV = ∯ n·Q dA — "volume integral of the divergence = net outflux" [#80] | 2.12 | B | NOTE | C14 | stated as the vector case of (2.30) with `divergence_theorem_box` (Q = (x, y, z): 3 = 3) and the cited sphere benchmark F = (2x, y², z²) → 8π/3; the outflux reading; D20 demoted (§4c) |
| N65 | Fig. 2.9: volume V, surface A, element n dA, element dV [#81] | 2.12 | B | NOTE | C14 | the block's visualization: our plotly 3-D box with faces coloured by n·Q and bars for both sides (mesh slider, §6) |
| N66 | Generalised field derivative 𝒟Q = lim_{V→0} (1/V) ∯ n_i Q dA (2.31) — coordinate-free gradient of any order [#82] | 2.12 | B | NOTE | C15 | stated with Eq. (2.31) and `integral_gradient` converging to `gradient` (order 2 in h, printed); **derived in D21** with (2.32), (2.33) (book: "limiting form") |
| C15 | **Integral definition of the divergence** ∇·Q = lim_{V→0} (1/V) ∯ n·Q dA (2.32): outflux per unit volume [#83] | 2.12 | A | CORE | – | load-bearing: the coordinate-free picture behind continuity and the control-volume forms of ch04; D21 and D22 (Ex. 2.5) derived here |
| N67 | Integral definition of the curl ∇×Q = lim (1/V) ∯ n × Q dA (2.33) [#84] | 2.12 | B | NOTE | C15 | stated with Eq. (2.33) and `integral_curl` → `curl` (one line); pointer ch05 |
| N68 | Ex. 2.5: Cartesian divergence from (2.32) — Taylor to the six face centres, opposite faces cancel the zeroth order, divide by the volume, take the limit [#85] | 2.12 | B | NOTE | C15 | the example is stated here; its derivation is the block's D22 (written out, animated face by face in §6 and in the `gauss_flux_box` explainer); pointer: the same box argument gives continuity in ch04 §4.2 |
| N69 | Stokes orientation: choose the outside → n; t = n_c × n counterclockwise seen from outside; (n_c, n, t) right-handed [#86] | 2.13 | B | NOTE | C16 | stated with the rule, a right-hand-rule primer and `boundary_tangent`, `planar_loop` (∮ (x − c) × t ds = 2A n as the check) |
| N70 | Fig. 2.10: open surface A, boundary C, n dA, n_c, t [#87] | 2.13 | B | NOTE | C16 | the block's visualization: our plotly 3-D cap with the three unit vectors at a boundary point |
| C16 | **Stokes' theorem** ∬_A (∇×u)·n dA = ∮_C u·t ds (2.34); the right side is the circulation [#88] | 2.13 | A | CORE | – | load-bearing: circulation is the quantity of Kelvin's theorem (ch05), of Kutta–Joukowski lift (ch06, ch14) and of vortex tubes; D26 (book never proves it) derived here |
| N71 | Normal curl as circulation per unit area: n·(∇×u) = lim_{A→0} (1/A) ∮ u·t ds (2.35) [#89] | 2.13 | B | NOTE | C16 | stated with Eq. (2.35) and `integral_curl_component` converging to n·curl (order 2 printed); pointer: vorticity = twice the angular velocity (ch03), ch05 |
| N72 | Remark: the integral definitions of ∇, ∇·, ∇× do not depend on the coordinate system (Exercises 2.16–2.18) [#90] | 2.13 | C | NOTE | C15 | named; pointer: cylindrical and spherical forms in ch03 / Appendix B |
| N73 | Ex. 2.6: Cartesian curl from (2.35) — rectangles in the three coordinate planes, midpoint values on each side [#91] | 2.13 | B | NOTE | C16 | stated with the x-component recipe written correctly (the second bracket contains u_y, not the printed u_z), the orientation of the four sides, and `integral_curl_component` vs (2.25); its bookkeeping is reused by D26; D23 demoted (§4c) |
| N74 | Comma notation A_,i = ∂A/∂x_i; ∇·u = u_i,i; (∇×u)_i = ε_ijk u_k,j (2.36) [#92] | 2.14 | B | NOTE | C10 | stated in its own short §2.14 section with Eq. (2.36), `expand_indices("u_i,i")`; a comma index is a tensor index (Cartesian only); pointer §5.6; D24 demoted (§4c) |
| S01 | Exercises 2.1–2.20 [#93] | Ex. | C | SKIP | – | pointer line: "the identities we need from them (2.5, 2.7–2.11, 2.19, 2.20) are tested in `tests/test_ch02.py`; exercise text is private" |
| S02 | Literature (Sommerfeld 1964) and supplemental reading (Aris 1962, Prager 1961) [#94] | Lit. | C | SKIP | – | pointer line: bibliography; Sommerfeld's tetrahedron argument is our D05/D06 |

### 2a. What each A block contains (for the lesson-designer)
- **C01** picture: three boxes labelled 1, 2, 3 being added · question: "why does the book never write Σ?" · maths: (2.2),
  (2.9), (2.17) expanded by hand · number: (1,2,3)·(4,5,6) = 32; δ_ij u_j for i = 2 · code: `expand_indices`,
  `dot`, `inner`, `kronecker_delta` + from-scratch triple loop · figure: N23 highlight · B/C: N01–N08, N12, N21–N24, N40,
  N41, N45.
- **C02** picture: two rulers through one origin · question: "what does one number C_ij measure?" · D02, D03 · number:
  θ = 30°, C = [[0.866, −0.5],[0.5, 0.866]], CᵀC = I by hand · code: `rotation_matrix_2d/3d`, `direction_cosines`,
  `is_orthogonal`, `random_rotation` · figure: N09 (3-D, two frames) + §6 animation · B/C: N09, N10, N14, N15, N17.
- **C03** picture: the same arrow, two shadows · D01 · number: x = (1, 2) at 30° → x' = (1.866, 1.232) · code:
  `transform_vector`, `inverse_transform_vector`, `transforms_as_vector` + from-scratch loop · figure: residual table of
  N18, polar figure of N19 · explainer E1 · B/C: N11, N13, N16, N18–N20.
- **C04** picture: the cube with arrows · question: "which index is the face?" · maths: sign convention as a rule
  table · number: τ = [[−p, a, 0],[a, −p, 0],[0, 0, −p]] read face by face · code: `cube_face_tractions`,
  `stress_component_meaning` · figure: N25 3-D cube · B/C: N25, R01, N26.
- **C05** picture: a slanted cut through the cube · D05 · number: Ex. 2.2 at 30° (N38) · code: `traction`,
  `normal_shear_stress`, `example_2_2` + from-scratch loop · figure: N36 tetrahedron, §6 slider · explainer E2 · B/C:
  N34–N39.
- **C06** picture: the stress element seen by two observers · D06 · number: τ' for Ex. 2.2's C, τ'₁₁ = √3a/2 · code:
  `transform_tensor`, `transforms_as_tensor` + from-scratch double loop vs `C.T @ tau @ C` · figure: §6 slider (τ'(θ),
  Mohr) · explainer E1 tensor mode · B/C: N27–N29.
- **C07** picture: a chain of indices closing on itself · D18 · number: I₁, I₂, I₃ of a 3 × 3 before and after a random
  rotation · code: `trace`, `invariants`, `contract`, `double_dot` + from-scratch · figure: invariants flat over 50
  rotations · B/C: N30–N33.
- **C08** picture: the 3 × 3 × 3 cube of ε · D09 · number: ε₁₂₃ ε₃₁₂ … the 81-case table summarised · code:
  `levi_civita`, `permutation_sign`, `epsilon_delta_residual`, `cross_einsum`, `is_isotropic` + from-scratch parity loop
  · figure: ε cube · B/C: N42–N44, N46–N50.
- **C09** picture: a hill's contour map · question: "which way is steepest?" · maths: ∇φ ⊥ level sets via the chain
  rule along a contour · number: φ = x² + y² at (1, 2): ∇φ = (2, 4), ∂φ/∂n along n = (1, 0)/1 → 2 · code: `grid`,
  `partial`, `gradient`, `directional_derivative` + from-scratch stencil · figure: N51 + §6 slider · B/C: R02, N51.
- **C10** picture: a box with more leaving than entering · maths: (2.23) · number: Ex. 2.3 a x → 3a · code:
  `divergence`, `vector_gradient`, `tensor_divergence`, `is_solenoidal` + from-scratch · figure: heatmaps of N55,
  order plot · B/C: N52, N74.
- **C11** picture: a paddle wheel in a shear flow · D12 · number: b × x → 2b · code: `curl`, `curl_components`,
  `is_irrotational` + from-scratch · figure: N55 quiver + tracer animation · explainer E5 · B/C: N53–N55.
- **C12** picture: a square that stretches and spins at once · D14, D15 · number: G = [[0, Γ],[0, 0]] → S, A, ω = −Γ e₃/2
  (hmm sign: ω from A = ½(G − Gᵀ) gives ω₃ = −Γ/2 with G[i,j] = ∂u_i/∂x_j — pinned by test F19) · code: `symmetric_part`,
  `antisymmetric_part`, `antisymmetric_from_vector`, `vector_from_antisymmetric`, `symmetric_double_contraction` +
  from-scratch · figure: §6 three-panel animation · explainer E3 · B/C: N56–N61.
- **C13** picture: the ellipse n·τ·n over all directions · D17 · number: Ex. 2.4 (N62) by hand · code:
  `principal_axes`, `diagonalize`, `characteristic_polynomial`, `normal_stress_bounds`, `example_2_4` + from-scratch
  2 × 2 quadratic formula · figure: N63 + Mohr slider · explainer E2 · B/C: N62, N63.
- **C14** picture: a box cut into slabs · D25 · number: Q = (x, y, z) on the unit cube: 3 = 3; sphere 8π/3 · code:
  `gauss_gradient_box`, `divergence_theorem_box`, `divergence_theorem_sphere` + from-scratch six-face sum · figure:
  N65 3-D box, mesh slider · explainer E4 · B/C: N64, N65.
- **C15** picture: the box shrinking to a point · D21, D22 · number: (1/V)∮ for Q = (x², 0, 0) at x₀ = 1 with h = 0.2 →
  2.000 (exact 2) · code: `integral_divergence`, `integral_gradient`, `integral_curl` · figure: §6 face-bookkeeping
  animation, h slider · explainer E4 · B/C: N66–N68, N72.
- **C16** picture: a loop in a whirlpool · D26 · number: b × x, disc R = 1, b = e₃: circulation 2π = curl flux 2π ·
  code: `circulation`, `curl_flux`, `stokes_theorem_check`, `integral_curl_component` + from-scratch midpoint sum ·
  figure: N70 3-D cap, shrinking-loop animation · explainer E5 · B/C: N69–N71, N73.

## 3. Section coverage
| § | Title | A | B | C | RECAP | SKIP |
|---|---|---|---|---|---|---|
| 2.1 | Scalars, Vectors, Tensors, Notation | C01 | N01, N02, N03, N05, N06, N07, N08 | N04 | — | — |
| 2.2 | Rotation of Axes: Formal Definition of a Vector | C02, C03 | N09, N10, N11, N12, N14, N15, N18, N19 | N13, N17, N20 | — | — |
| 2.3 | Multiplication of Matrices | — (inside C01, C03) | N16, N21, N23 | N22 | — | — |
| 2.4 | Second-Order Tensors | C04, C06 | N25, N27, N28, N29 | N24, N26 | R01 | — |
| 2.5 | Contraction and Multiplication | C07 | N31, N32, N33 | N30 | — | — |
| 2.6 | Force on a Surface | C05 | N34, N35, N36, N38 | N37, N39 | — | — |
| 2.7 | Kronecker Delta and Alternating Tensor | C08 | N40, N41, N42, N43, N44 | — | — | — |
| 2.8 | Vector, Dot, and Cross Products | — (inside C08, C01) | N46, N47, N49 | N45, N48, N50 | — | — |
| 2.9 | Gradient, Divergence, and Curl | C09, C10, C11 | N51, N52, N53, N54, N55 | — | R02 | — |
| 2.10 | Symmetric and Antisymmetric Tensors | C12 | N56, N57, N58, N59, N60 | N61 | — | — |
| 2.11 | Eigenvalues and Eigenvectors of a Symmetric Tensor | C13 | N62, N63 | — | — | — |
| 2.12 | Gauss' Theorem | C14, C15 | N64, N65, N66, N67, N68 | — | — | — |
| 2.13 | Stokes' Theorem | C16 | N69, N70, N71, N73 | N72 | — | — |
| 2.14 | Comma Notation | — (inside C10) | N74 | — | — | S01, S02 (end of chapter) |

Totals: A 16 · B 60 · C 18 = 94 (R01 counted under C, R02 under B, S01/S02 under C). No section is empty; §2.3, §2.8 and
§2.14 have no A item and follow decision 6 (own short notebook section, paragraphs tagged with the A parent).

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not A/B/RECAP)
Rows marked **gloss** are one plain sentence where the item is used; all others are full 📎 primers (one or two
sentences, what it means here, a 2–4-line runnable demo with easy numbers). Items already in `knowledge/primers.md` get
a one-line reminder naming the ch01 primer (listed at the end).

| Concept or tool | Needed by | Why it is not A/B/RECAP |
|---|---|---|
| `np.einsum` index strings (`'i,i'`, `'ik,kj->ij'`, `'ijk,i,j->k'`): the code form of the summation convention | C01, C07, C08, C11 | Python tool; not in primers.md |
| matrix multiplication by hand and `@`; transpose `.T`; identity `np.eye` | C01 (N21), C03 (N16), C06 | P53 primed determinants and minors only |
| subscript indices i, j, k ∈ {1, 2, 3} and the "≡ defined as" symbol | C01 | notation **gloss** |
| orthonormal basis, projection onto a unit vector, completeness Σ_j e'_j e'_jᵀ = I ("a vector is the sum of its projections") | C02 (D02), C03 (D01) | pre-book linear algebra |
| bilinearity of the dot product (distributes over sums, scalars pull out) | D01 | **gloss** inside D01's *why* |
| cosine of the angle between unit vectors; cos(π/2 − θ) = sin θ, cos(θ + π/2) = −sin θ; radians vs degrees (`np.deg2rad`) | C02, N19 | pre-book trigonometry |
| "true for every n ⇒ the coefficients agree" (compare coefficients of an identity) | D06 | logic move, **gloss** |
| limits and orders of smallness: face terms ∝ h², volume terms ∝ h³, so volume terms vanish first | D05, D22 | ch01 only glossed it (§1.7); D05 hinges on it → full primer |
| vector area of a closed surface Σ n dA = 0, hence dA_i = n_i dA for the tetrahedron | D05, N35 | pre-book geometry |
| `np.arctan2(y, x)`: the four-quadrant angle (handles a < 0 in Ex. 2.2) | N38, `angle_between` | Python tool |
| permutations, cyclic order and parity (even / odd) | C08 | pre-book combinatorics |
| numpy arrays with three axes; `np.transpose(a, axes)`; `arr[i, j, k]` | C08 (N43), N28 | Python tool |
| `np.meshgrid`; the project grid layout `[k, j, i]` = (z, y, x) with x on the last axis; vector fields stacked on axis 0 (`u[c, k, j, i]`) | C09, C10, C11, C14, C15, C16 | Python tool; the project-wide convention is fixed here (analysis §9.8) |
| numpy broadcasting (an array op applied to every grid point at once) | C09–C11 | Python tool (P03 primed arrays, not broadcasting) |
| level sets (contour lines / surfaces) and the directional derivative ∂φ/∂n | C09 | pre-book calculus |
| `plt.contour`, `plt.quiver`, `plt.streamplot` | C09, C10, C11 | matplotlib tools (P01 primed `plot` only) |
| product rule for partial derivatives; ∂x_i/∂x_j = δ_ij | N55, D22 | P38 primed differentials; the δ_ij form is new — **gloss** |
| eigenvalues and eigenvectors: A·b = λb, det(A − λI) = 0, `np.linalg.eigh`, `np.roots` | C13 | pre-book (the book says "the reader is assumed familiar") |
| complex conjugate z̄ and \|z\|² = z z̄ | D17 fact (1) | P45 primed i² = −1 and Euler only |
| quadratic form n·τ·n (Rayleigh quotient) and its extreme values; Gram–Schmidt for repeated eigenvalues (**gloss**) | D17 facts (2), (4) | pre-book linear algebra |
| Vieta's formulas: coefficients of a cubic ↔ sum, pair-sums and product of its roots | D18 | pre-book algebra |
| `np.linalg.norm`; `np.linalg.qr` (random rotations — **gloss** in a code comment) | C02, C07, N18 | Python tools |
| volume and surface integrals as sums: midpoint rule in 2-D and 3-D, iterated integrals | C14, C15 | P27/P37 primed 1-D integrals only |
| fundamental theorem of calculus ∫_a^b f′ dx = f(b) − f(a) | D25 | P27 primed the definite integral, not the theorem |
| mean-value theorem for integrals (∫_V f dV = f(x*) V for some x* in V) | D21 | pre-book analysis |
| line integral of a vector field ∮ u·t ds along a parametrised curve; parametrising a circle and a rectangle | C16, D26 | P35 primed ∫p dv along a path; the vector form and parametrisation are new |
| right-hand rule and the orientation of a boundary | N46, N69, C16 | pre-book geometry |
| `scipy.linalg.expm` (deforming square x(t) = e^{Gt} x₀) | C12 figure, E3 | Python tool |
| plotly 3-D arrows and lines (`go.Cone`, `go.Scatter3d` lines, `go.Mesh3d`) | C04, C05, C14, C16 figures | P41 primed `go.Surface`/`Scatter3d` points only |
| observed order of convergence from a log–log slope | C09 (F16 test), C15, C16 | P13 + ch01 `tools/convergence.observed_order` — reminder sentence |

Reminders only (already primed in ch01): partial derivative P25 · finite differences P21/P22 · first-order Taylor P26 ·
definite integral P27 · trapezoid rule P37 · chain rule P49 · matrices/determinants/minors P53 · sympy P40 · sympy Matrix
P61 · functions as arguments / lambda P29 · `assert np.allclose` P15 · `slider_figure` P17 · `animate` P16 · `show_viz`
P18 · f-strings P04 · dictionaries P23 · Newton's second law / free-body diagram P09 · stress P05 · np.random P10.

## 4b. Derivations written out (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Written out because (a) the result belongs to an A item, or (b) the book never writes it out and the lesson needs it
(D02, D03, D09, D15 belong to B items and are kept under rule (b)). One small move per step; the book's skipped moves
(analysis §2b) are filled in and marked. 15 rows.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | vector transformation x'_j = x_i C_ij (2.5) from (2.1) and (2.3) | C03 | ★ | 6 | dot product (N07), bilinearity of the dot product (gloss), orthonormality e'_i·e'_j = δ_ij (N40, primer orthonormal basis), Kronecker substitution (N41), naming C_ij (C02) | forgetting that e'_i·e'_j = δ_ij kills every primed term but one; writing C_ji for C_ij (row = old axis); doing j = 1 and saying "similarly" — all three projections are one indexed line | notebook · `rotation_of_axes` |
| D02 | orthogonality C_ij C_kj = δ_ik and C_ji C_jk = δ_ik; det C = ±1, +1 for a rotation (Exercise 2.8; book never writes it) | C02 | ★★ | 8 | completeness Σ_j e'_j e'_jᵀ = I (primer), matrix product `@` (primer), determinant of a product (P53 reminder), continuity to the identity for the sign | proving CᵀC = I and assuming CCᵀ = I follows without the second completeness argument; det = −1 is a reflection, excluded by "rotation" | notebook · `rotation_of_axes` |
| D03 | inverse transformation x_j = x'_i C_ji (2.7) (Exercise 2.2; book: "it can be shown") | C02 | ★ | 4 | multiply (2.5) by C_kj and sum on j (N12 dummy renaming), D02, Kronecker substitution (N41) | the summed index of C moves from the first to the second slot; C⁻¹ = Cᵀ holds only because C is orthogonal | notebook |
| D05 | Cauchy's traction formula f_i = τ_ji n_j (2.15) from the tetrahedron | C05 | ★★ | 9 | free-body balance (P09 reminder), stress sign convention on the −e_j faces (C04), vector area of a closed surface dA_j = n_j dA (primer), orders of smallness h² vs h³ (primer), limit h → 0 | sign of the tractions on the coordinate faces (outward normal −e_j gives −τ_j·); dropping body-force and inertia terms without saying why (they are ∝ h³); writing τ_ij n_j, which needs symmetry (Ch. 4) | notebook · `cauchy_traction_principal_axes` |
| D06 | tensor transformation rule τ'_mn = C_im C_jn τ_ij, τ' = Cᵀ·τ·C (2.12) from (2.15), (2.8), (2.7) (book states it via Sommerfeld's tetrahedron) | C06 | ★★ | 8 | f is a vector (2.8) (C03), (2.15) (C05), (2.7) (N14), dummy renaming (N12), "true for every n' ⇒ coefficients agree" (gloss) | the book's order (2.12 before 2.15) is reversed here and the notebook says so; the final rename i ↔ j; matrix form CᵀτC, not CτCᵀ (passive C) | notebook · `rotation_of_axes` |
| D09 | epsilon–delta relation ε_ijk ε_klm = δ_il δ_jm − δ_im δ_jl (2.19); contractions ε_pqi ε_pqj = 2δ_ij, ε_pqr ε_pqr = 6; triple product a × (b × c) = (a·c)b − (a·b)c (Exercises 2.5, 2.7; book: "verify by choosing values") | C08 | ★★ | 10 | antisymmetry of ε in any index pair (N43), case enumeration (both sides vanish unless {i, j} = {l, m}, i ≠ j), Kronecker substitution (N41), δ_ii = 3, (2.21) (N49) | δ_ii = 3, not 1 (so 3δ − δ = 2δ); the sign when moving the summed k next to its partner; the 81-case table is a check, not the proof | notebook |
| D12 | curl components (2.25) from the index form (2.24) | C11 | ★ | 5 | ε enumeration (C08), (2.21) with u_i → ∂/∂x_j (N49), operator ordering (∂ acts on u_k) | writing u_k ∂_j instead of ∂_j u_k; the minus from ε_132 = −1; doing i = 1 and saying "similarly" — i = 2, 3 are cyclic relabellings and that is said | notebook · `stokes_circulation_loop` |
| D14 | decomposition B_ij = S_ij + A_ij with S = ½(B + Bᵀ), A = ½(B − Bᵀ), unique; S and A are tensors | C12 | ★ | 6 | transpose (primer), linearity of (2.12) (C06), the "both symmetric and antisymmetric ⇒ zero" argument | thinking the split depends on the frame (it does not: both parts transform by (2.12)); forgetting the uniqueness half; component count 6 + 3 = 9 | notebook · `strain_vs_rotation_split` |
| D15 | R_ij = −ε_ijk ω_k and ω_k = −½ ε_ijk R_ij (2.26)–(2.27); R·x = ω × x | C12 | ★★ | 7 | (2.19) contraction 2δ (D09), index moves on ε (N43), (2.21) (N49), antisymmetry of ε ⇒ R antisymmetric | the sign: R_ij x_j = −ε_ijk ω_k x_j = +ε_ikj ω_k x_j = (ω × x)_i needs a one-place move; the book's lower limit "i−1" is a misprint for i = 1 | notebook · `strain_vs_rotation_split` |
| D17 | §2.11 facts for a real symmetric τ: (1) real eigenvalues; (2) orthogonal eigenvectors for distinct λ; (3) τ' = Cᵀ·τ·C = diag(λ¹, λ², λ³) with C = [b¹ b² b³]; (4) n·τ·n ∈ [λ_min, λ_max] for every unit n and max shear (λ_max − λ_min)/2 (book: "can be proved") | C13 | ★★★ | 14 | eigen-equation (primer), complex conjugate (primer), symmetry τ_ij = τ_ji (N56), orthonormality and D02, (2.12) (C06), expansion of n in the eigenbasis (primer orthonormal basis), quadratic form / Rayleigh quotient (primer), Gram–Schmidt for repeated λ (gloss), sympy check cell | symmetry is used twice (realness and orthogonality) — say where; repeated eigenvalues: orthogonal eigenvectors exist but are not unique (`eigh` picks one); flip one b if det C = −1; the book's "λ bound the elements τ_ij" is loose — the sharp statement is about normal stresses n·τ·n and the shear bound | notebook · `cauchy_traction_principal_axes` |
| D18 | invariants I₁ = A_ii, I₂ = ½(I₁² − A_ij A_ji), I₃ = det A are unchanged under (2.12); characteristic polynomial λ³ − I₁λ² + I₂λ − I₃ = 0; in the principal frame I₁ = Σλ, I₂ = Σ_{k<l} λ^k λ^l, I₃ = Πλ (Exercise 2.9; book gives hints) | C07 | ★★ | 8 | (2.12) (C06), D02 orthogonality as "a contracted C-pair is δ", determinant of a product (P53 reminder), Vieta's formulas (primer), closed index chains | I₂ uses A_ij A_ji (a closed chain), not A_ij A_ij (= tr(A Aᵀ): also frame-independent, but not the λ-coefficient of the cubic; equal to A_ij A_ji only for symmetric A); the sign pattern of the polynomial; λ^k is a label, not a power | notebook |
| D21 | integral definitions (2.31)–(2.33) as small-volume limits of Gauss' theorem (2.30) (book: "the limiting form") | C15 | ★★ | 7 | (2.30) (C14), mean-value theorem for integrals (primer), limits (primer), contraction with e_i for (2.32), ε-cross for (2.33) (N49) | forgetting the 1/V; assuming the limit is independent of the shape of V without continuity of ∂Q; (2.32) and (2.33) are (2.31) contracted / crossed, not new theorems | notebook · `gauss_flux_box` |
| D22 | Ex. 2.5: Cartesian divergence (2.23) recovered from (2.32) on a box Δx₁Δx₂Δx₃ | C15 | ★★ | 9 | first-order Taylor (P26 reminder), midpoint rule — a face integral of a linear term equals the centre value (primer), face normals ±e₁ (C04's faces), orders of smallness (primer), (2.23) (C10) | zeroth-order terms cancel only between *opposite* faces; a "face value" is really a face average (exact for linear variation, error ∝ Δx⁵ against the kept Δx³); which face carries +e₁; e₁·∂Q/∂x₁ = ∂Q₁/∂x₁ | notebook · `gauss_flux_box` |
| D25 | Gauss' theorem (2.30): for a box from the fundamental theorem of calculus, then for any V by tiling (book never proves it) | C14 | ★★ | 9 | fundamental theorem of calculus (primer), iterated integrals (primer), outward normal n_i = ±1 on two faces and 0 on four, cancellation of interior faces (opposite normals), hypotheses Q ∈ C¹ | four of the six faces have n₁ = 0 and contribute nothing to the i = 1 term; interior faces cancel only because the tiles share faces with opposite normals; the boundary of the tiling approximates A | notebook · `gauss_flux_box` |
| D26 | Stokes' theorem (2.34) for a planar surface: one small rectangle, then tiling; corollary ∮ ∇φ·t ds = 0 ⇒ ∇×∇φ = 0 (Exercise 2.20) (book never proves it) | C16 | ★★ | 10 | line integral of a vector field on a parametrised curve (primer), orientation t = n_c × n (N69), midpoint values / Taylor (P26 reminder), Ex. 2.6 bookkeeping (N73), tiling with shared edges traversed twice in opposite directions | orientation of each side of the rectangle (the u_y sides at z ± Δz/2 run opposite ways); interior edges cancel only when every tile is oriented by the same n; a field singular inside the loop (irrotational vortex core) breaks the hypothesis — the explainer's ⚠️ status | notebook · `stokes_circulation_loop` |

Counts: ★ 4 (D01, D03, D12, D14) · ★★ 10 (D02, D05, D06, D09, D15, D18, D21, D22, D25, D26) · ★★★ 1 (D17, shown in
`cauchy_traction_principal_axes` with a sympy check cell that re-runs the construction: builds C from `eigh`, checks
CᵀτC diagonal symbolically for a symbolic 2 × 2, and n·τ·n = Σ λ^k (n·b^k)²).

## 4c. Derivations demoted to statements (first column is the A parent in bold, e.g. **C20** — never a bare ID or a D id: the parser reads a bare first-cell ID as an item and blanks its tier)
Results **given, not derived**, inside the named A block (a paragraph, the equation, a number). 11 rows.

| A parent | Analysis §2b item | Result stated (Eq.) | Stated in (B item) | Why not written out |
|---|---|---|---|---|
| **C03** | D04 Ex. 2.1 polar components | u_r = u₁cos θ + u₂sin θ, u_θ = −u₁sin θ + u₂cos θ; C = [[cos θ, −sin θ],[sin θ, cos θ]] | N19 | the book writes it out; it is (2.5) with j ∈ {r, θ} and becomes C03's worked number (traced with u = (1, 2), θ = 30°), not a separate derivation |
| **C05** | D07 Ex. 2.2 traction and rotated stress in a shear flow | f = (a sin φ, a cos φ); \|f\| = \|a\|; θ = 60° / 240°; τ'₁₁ = a sin 2φ, τ'₁₂ = a cos 2φ | N38 | the book writes it out; it is C05's worked number (by hand at 30°, then `example_2_2`); the general double-angle forms are stated and reappear in C13 (Mohr) |
| **C01** | D08 (2.17) δ_ij u_j = u_i | δ replaces its summed index; δ_ij A_jk = A_ik; δ_ii = 3 | N41 | one-line enumeration the book writes; `expand_indices` prints it |
| **C08** | D10 (2.20) cross-product components from its defining properties | (2.20) | N47 | pre-book result the lesson does not need to derive (bilinearity would have to be assumed anyway); stated with the one-sentence reason e_i × e_j = ε_ijk e_k; `cross` = `np.cross` asserted |
| **C08** | D11 (2.21) index form of the cross product | (u × v)_k = ε_ijk u_i v_j = ε_kij u_i v_j | N49 | the book does the k = 1 check; the printed expansion of `expand_indices("eps_ijk u_i v_j")` shows all three components — a check, not a derivation |
| **C11** | D13 Ex. 2.3 div and curl of a x and b × x | 3a, 0; 0, 2b | N55 | the book writes every partial derivative; stated as the worked number of C10/C11 with the index route in one line (uses (2.19) from D09) and the a x₃ e₃ correction |
| **C12** | D16 (2.28) → (2.29) ⇒ τ_ij A_ij = 0 | τ_ij B_ij = τ_ij S_ij | N59, N60, N61 | the book writes the rename-and-swap chain in full; stated with the one-line lemma (X = −X ⇒ X = 0) and P_A = 0 printed |
| **C13** | D19 Ex. 2.4 principal axes of a plane shear strain rate | λ = ±Γ; b¹, b²; C = 45° rotation; S' = diag(Γ, −Γ) | N62 | the book writes it out; it is C13's worked number traced by hand (2 × 2 determinant, quadratic formula, normalisation, det C = +1 choice); the general angle ½ atan2(2S₁₂, S₁₁ − S₂₂) is stated |
| **C14** | D20 divergence theorem from (2.30) | ∭ ∇·Q dV = ∯ n·Q dA | N64 | one move (set Q → Q_i and sum on i, allowed because (2.30) is linear in Q) — said in a sentence |
| **C16** | D23 Ex. 2.6 Cartesian curl from (2.35) | (∇×u)_x = ∂u_z/∂y − ∂u_y/∂z (second bracket with u_y), the others by cyclic permutation | N73 | the book writes it out (with the u_z misprint corrected in our statement); D26 runs the same rectangle bookkeeping in the other direction, so the moves are shown there |
| **C10** | D24 (2.36) comma-notation forms of ∇·u and ∇×u | u_i,i and ε_ijk u_k,j; a comma index transforms as a vector index | N74 | notation; the one non-trivial fact (∂/∂x'_j = C_ij ∂/∂x_i by the chain rule) is stated in a sentence with a pointer to P49 |

## 5. Interactive explainers (4–5 + backup)
Five explainers on the A items where manipulation teaches most; one backup. Each has the required live **Explain** tab
("Explanation & interpretation", numbered sections computing every number on screen with the reader's settings, then
"Reading the current setting"), a synced **Code** tab, a **Derivation** tab for its D ids, and ≥ 2 more depth features.
Physics mirrors `fluidpy` functions that accept plain floats/lists (analysis §6 parity list). Colours: old frame /
Cartesian teal, new frame / rotated orange, normal stress blue, shear rose, symmetric part teal, antisymmetric part
orange, flux in blue / out orange, circulation purple.

### E1 · rotation_of_axes
- **A:** C02, C03, C06 (also shows N15 orthogonality residual, N16 x' = Cᵀx, N18 the vector test, N19 Ex. 2.1 as the
  "polar" preset, N27 tensor ≠ matrix as the "fixed array" preset, N38 τ' of Ex. 2.2 as the "pure shear" preset)
- **Confusion removed:** "rotating the axes" is confused with "rotating the vector" (passive C vs active R = Cᵀ; the
  book's C equals Wikipedia's R(θ) but is applied as Cᵀ); and (2.12) looks like a new rule when it is (2.8) applied once
  per index.
- **Why interactive:** the whole idea is that one object has different numbers in different frames. Dragging θ while
  the arrow stays put, and watching the C matrix, the components x'_j and the tensor components τ'_mn move together, is
  the experience a static pair of frames cannot give; the "what rotates?" toggle shows the same C doing opposite jobs.
- **Stage:** (1) the plane with the fixed arrow x (or the stress element in tensor mode), the old axes (teal) and the
  rotated axes (orange), both component sets drawn as projections; (2) the direction-cosine matrix view (E5's
  matrix-view pattern): C with row i / column j lit as x'_j = Σ_i x_i C_ij is assembled term by term, plus CᵀC printed
  live; (3) components vs θ: x'_1(θ), x'_2(θ) (vector mode) or τ'₁₁, τ'₁₂, τ'₂₂ vs θ (tensor mode) with a moving dot
  and the special angles marked (θ = 45° where τ'₁₂ = 0 for pure shear — the bridge to E2).
- **Controls:** rotation angle θ (−180° … 180°) · drag the arrow tip (x₁, x₂) · mode: vector / tensor / polar (Ex. 2.1
  with the point moving on a circle) · tensor preset chips (pure shear a, uniaxial, hydrostatic −p, fixed non-tensor
  array) · toggle "what rotates: axes (passive, book) / vector (active)" · optional: 3-D axis of rotation (Rodrigues C)
  for the vector mode.
- **Equations shown (live):** (2.4), (2.5)/(2.6), (2.7), (2.8), (2.12), C_ij C_kj = δ_ik; the polar C of Ex. 2.1.
- **Mirrors:** `core.tensors.rotation_matrix_2d`, `rotation_matrix_3d`, `transform_vector`,
  `inverse_transform_vector`, `transform_tensor`, `orthogonality_residual`, `ch02.polar_components`,
  `ch02.example_2_2` (τ_rot).
- **Derivations:** D01, D02, D06 (Derivation tab: D01's step "dot with e'_1" lights row 1 of C; D06's last step sets
  the tensor mode to pure shear so τ'₁₁ = a sin 2θ appears on the curve).
- **Depth features:** Explain tab (C from θ with numbers → x' term by term → CᵀC = I check → τ' = CᵀτC entry by entry →
  "the arrow's length is unchanged: \|x'\| = \|x\|" → reading the current setting), synced Code tab, + linked views (3),
  presets (θ = 0, 45°, 90°, the polar frame at θ, the pure-shear tensor), status verdict ("✅ orthogonal: CᵀC − I = 1e-16,
  det C = +1" / "⚠️ det C = −1: a reflection, not a rotation" for a mirrored preset / "❌ not a tensor: (2.12) residual
  0.83" for the fixed array), inspector (click a matrix cell: C_ij = e_i·e'_j = cos(angle) with the two unit vectors
  drawn), modes (vector / tensor / polar).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (windows linked by one state, numbered live derivation
  in the explanation) with E5 `buckingham_pi_machine`'s matrix view for the C matrix.
- **Aha:** the arrow never moved — only the ruler did; x' = Cᵀx and τ' = CᵀτC are the same ruler applied once per index.

### E2 · cauchy_traction_principal_axes
- **A:** C05, C13, C04 (also shows N38 Ex. 2.2 as the "channel shear" preset, N39 the channel element, N62 Ex. 2.4 as
  the "strain rate" preset, N32 n·τ vs τ·n toggle for a non-symmetric τ, N33-free; the Mohr circle previews ch04)
- **Confusion removed:** "stress is a vector" (it is nine numbers: the force per area depends on which plane you cut);
  "shear and normal stress are properties of the material point" (they are properties of the plane; two planes carry no
  shear at all and they are the eigenvectors); and the loose "eigenvalues bound the elements" (sharp: they bound the
  normal stress on every plane, the shear is bounded by half their spread).
- **Why interactive:** the phenomenon is a function of the cut direction. Dragging the normal n around the point and
  watching f, σ_n, τ_s and the Mohr point move at once is the definition of (2.15) made visible; finding the shear-free
  directions by hand and then seeing them coincide with `principal_axes` is the eigenvector idea *felt* before it is
  proved in the Derivation tab.
- **Stage:** (1) the 2-D stress element (a square with its τ₁₁, τ₂₂, τ₁₂ arrows, C04's sign convention) and a cut plane
  whose normal n (drag, angle φ) carries the traction arrow f = n·τ decomposed into σ_n (blue) and τ_s (rose); the
  eigen-axes drawn faintly as a ghost; (2) σ_n(φ) and τ_s(φ) over 0–180° with a moving dot, the principal angles and
  the max-shear angle marked, λ_min / λ_max as horizontal dashed bounds; (3) Mohr's circle (σ_n, τ_s) with the current
  point, centre (I₁/2), radius (λ_max − λ_min)/2 (hidePortrait; its numbers repeated in view 2's title).
- **Controls:** τ₁₁, τ₂₂, τ₁₂ (Pa or 1/s in strain mode) · plane angle φ (drag or slider) · presets (pure shear Ex. 2.2 —
  a > 0 and a < 0 halves, uniaxial, hydrostatic −p, Couette τ₁₂ = μU/h with water numbers, Ex. 2.4 strain rate) · toggle
  "contract first index n·τ (book) / second index τ·n" with a non-symmetric demo tensor (shows they differ) · optional:
  3-D tetrahedron view of the same n (Viz.three; hidePortrait).
- **Equations shown (live):** (2.15), dA_i = n_i dA, σ_n = n·f, τ_s = \|f − σ_n n\|, det(τ − λδ) = 0, (τ − λδ)·b = 0, τ' =
  Cᵀ·τ·C, σ_n = a sin 2φ and τ_s = a cos 2φ for the shear preset, the 2-D principal formula (σ₁₁ + σ₂₂)/2 ± √(…).
- **Mirrors:** `core.tensors.traction`, `normal_shear_stress`, `principal_axes`, `diagonalize`, `invariants`,
  `normal_stress_bounds`, `ch02.example_2_2`, `ch02.example_2_4`, `ch02.principal_angle_2d`, `ch02.traction_2d` (§8).
- **Derivations:** D05 (★★, each step draws the tetrahedron's face it balances; the "volume terms vanish" step shrinks
  the element), D17 (★★★, 14 steps; fact (2) sets φ to a principal angle so τ_s = 0 on screen; fact (4) sweeps φ and
  shows the σ_n curve inside the dashed bounds).
- **Depth features:** Explain tab (n from φ → f_i = τ_ji n_j term by term with numbers → σ_n, τ_s → the invariants →
  λ from the quadratic → principal angle → the bounds → reading the current setting: "this plane carries 87 % of the
  maximum shear"), synced Code tab, + linked views (3), presets, status verdict ("🎯 principal plane: τ_s = 0" / "🔪
  maximum shear plane" / "⚠️ non-symmetric τ: n·τ ≠ τ·n — Ch. 4 shows τ is symmetric" / "🌊 hydrostatic: every plane
  is principal"), term bars (f decomposed: τ₁₁n₁, τ₂₁n₂ | τ₁₂n₁, τ₂₂n₂ summing to f₁, f₂), inspector (click the Mohr
  point or a curve point: the arithmetic at that φ).
- **Follows reference:** `forced_damped_vibrations.html` (system + graphs on one parameter, the explanation panel
  with boxed numbers and a regime-dependent reading).
- **Aha:** the same nine numbers push differently on every plane through the point; the two planes with no shear are
  the eigenvectors, and their normal stresses bound everything else.

### E3 · strain_vs_rotation_split
- **A:** C12 (also shows N52 the velocity gradient G[i,j] = ∂u_i/∂x_j, N56 symmetric/antisymmetric definitions, N57–N58
  the ω ↔ R map, N59–N61 "a symmetric tensor cannot see A" as a term-bar check, N62–N63 Ex. 2.4's 45° axes, N55 Ex. 2.3
  b × x as the "solid-body rotation" preset; previews ch03 §3.4)
- **Confusion removed:** "simple shear is a rotation" (it is half pure stretching along 45° and half solid-body spin);
  "the antisymmetric part is just leftover algebra" (it is a vector: the local angular velocity ω/2… precisely, the
  rotation rate is ½∇×u = the vector of A); and the sign of R·x = ω × x.
- **Why interactive:** the decomposition is a statement about *motion*. A square carried by u = G·x deforms and spins at
  once; the same square under S alone only stretches, under A alone only spins, and playing all three on one clock —
  while sliders change G — makes B = S + A a thing you watch, not an identity you accept.
- **Stage:** (1) three panels (or one panel with toggles): a square of tracer points advected by the linear field
  u = G·x, u = S·x, u = A·x from the same start, with the velocity field's streamlines faint behind and the principal
  axes of S drawn on the S panel (stretch arrows blue, squeeze arrows rose), a small paddle wheel spinning at ω₃/2 on
  the A panel; (2) the matrix view G = S + A with the four (2-D) or nine entries as bars, S_ij = ½(G_ij + G_ji) assembled
  live; (3) the angle of a material line vs time: under S it tends to the stretching axis, under A it rotates uniformly,
  under G the sum — with the analytic ghost e^{Gt}.
- **Controls:** G₁₁, G₁₂, G₂₁, G₂₂ (1/s) · presets (simple shear u₁ = Γx₂ — Ex. 2.4, solid-body rotation b × x — Ex.
  2.3, pure strain (stretch + squeeze), uniaxial extension, irrotational strain with rotation removed) · show: G / S / A /
  all three · transport (t, play/step/scrub, loop) · optional: a 3-D preset (G 3 × 3 with the ω vector drawn; matrix
  view only).
- **Equations shown (live):** B_ij = ½(B_ij + B_ji) + ½(B_ij − B_ji), (2.26), (2.27), R·x = ω × x, (2.28)–(2.29) as the
  "τ:A = 0" bar, S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i) (Ex. 2.4), x(t) = e^{Gt} x₀.
- **Mirrors:** `core.tensors.symmetric_part`, `antisymmetric_part`, `antisymmetric_from_vector`,
  `vector_from_antisymmetric`, `symmetric_double_contraction`, `principal_axes` (2 × 2), `ch02.example_2_4`,
  `ch02.deform_square` / `ch02.linear_flow_map` (§8).
- **Derivations:** D14 (uniqueness step sets the preset to pure shear and shows both parts), D15 (the step R_ij x_j =
  (ω × x)_i sets the A-only view and draws ω × x on the wheel).
- **Depth features:** Explain tab (G from the preset → S and A entry by entry → ω from A with (2.27) → principal axes
  and rates of S → "after t = … s the square has stretched by e^{λt} = … along … and turned by … rad" → reading the
  current setting), synced Code tab, + linked views (3), transport (one clock for the three panels, end-of-run card with
  the stretch factor and the turned angle), presets, term bars (S:S, A:A and the cross term S:A = 0 for a symmetric
  probe tensor — the (2.29) check live), modes (show G / S / A), status verdict ("🔄 pure rotation: S = 0, shape kept" /
  "↔ pure strain: A = 0, no spin" / "🌀 simple shear: equal parts stretch and spin (Ex. 2.4)").
- **Follows reference:** `angular_frequency_explorer_1.html` (linked views on one clock, modes, presets, "Right now"
  notes, end-of-run summary).
- **Aha:** simple shear u₁ = Γx₂ is half pure stretching along the 45° axes and half solid-body spin at rate Γ/2 — the
  same square, two motions added.

### E4 · gauss_flux_box
- **A:** C14, C15 (also shows N64 the divergence theorem, N65 Fig. 2.9's coloured faces, N66 (2.31) as the limit view,
  N68 Ex. 2.5 face bookkeeping, N54 solenoidal fields as the "zero net flux" preset, N55 Ex. 2.3 fields as presets, C10
  the divergence heatmap)
- **Confusion removed:** "the divergence is a formula with partial derivatives" (it is a fact about a tiny box: what
  leaks out per unit volume; the formula is what that becomes in Cartesian coordinates) and "Gauss' theorem is a trick
  for computing integrals" (it is the statement that interior faces cancel when you add up boxes).
- **Why interactive:** the theorem has a *left side that lives inside* and a *right side that lives on the boundary*;
  moving and resizing the box over a field and watching the face fluxes (bars) re-sum to the same total as the
  integrated divergence, then shrinking the box until (1/V)∮n·Q dA settles on ∇·Q(x₀), is the (2.30) → (2.32) chain as
  an experiment. Subdividing the box into tiles and watching interior faces cancel is D25 on screen.
- **Stage:** (1) a 2-D field Q (x–y plane; the 3-D box has unit depth so the theorem's structure is identical — said in
  the subtitle) as a heatmap of ∇·Q (blue in / orange out) with streamlines and a draggable, resizable box whose four
  faces carry arrows n·Q; optional "tile it" split into 2 × 2 or 4 × 4 sub-boxes with interior arrows drawn and
  cancelling; (2) term bars: flux through each face (signed) summing to the total, next to the bar ∫∇·Q dV — equal to
  the pixel; (3) the limit view: (1/V)∮n·Q dA vs box side h on log axes approaching the dashed line ∇·Q(x₀), with the
  order-2 slope as a ghost.
- **Controls:** field preset (a x — Ex. 2.3, b × x — solenoidal, Q = (x², 0), a point source, (2x, y²) sin/cos smooth,
  a "sink in the corner") · box centre (drag) · box side h (log slider) · tiles (1 / 2 × 2 / 4 × 4) · toggle "show the
  Taylor face values of Ex. 2.5" (face value = centre ± ½h ∂Q/∂x).
- **Equations shown (live):** (2.30), the divergence theorem form, (2.31), (2.32), (2.23), Ex. 2.5's face expansion
  [Q]_{face} = Q(x) ± (Δx₁/2) ∂Q/∂x₁ + …, the 8π/3 sphere benchmark in the Equations tab.
- **Mirrors:** `core.integral_theorems.divergence_theorem_box` (with a z-independent field and unit depth),
  `divergence_theorem_rect2d` (§8), `integral_divergence`, `flux_through_faces` (§8), `core.operators.divergence`,
  `ch02.radial_field`, `ch02.solid_body_rotation_field`.
- **Derivations:** D25 (step "tile V and cancel interior faces" sets tiles = 4 × 4 and highlights the interior
  arrows), D22 (each face pair's step lights the two faces and shows their Taylor values with numbers), D21 (the
  1/V → limit step scrubs h down and points at the limit view).
- **Depth features:** Explain tab (the field at the box centre with numbers → ∇·Q by the formula → each face's flux as
  n·Q × length → the sum → the volume integral → their difference and why it shrinks with h → reading the current
  setting: "this box is a net source of 0.42 m³/s per m of depth"), synced Code tab, + linked views (3), term bars
  (face fluxes → total vs ∫∇·Q), presets, status verdict ("🟠 net outflow: source inside" / "🔵 net inflow" / "⚖️
  solenoidal: every box has zero net flux" / "📐 box small enough: (1/V)∮ within 1 % of ∇·Q"), inspector (click a
  face: its midpoint values, normal, n·Q, length, flux), transport-like scrub on h.
- **Follows reference:** `fid_formula_lab.html` (term-by-term bars summing to the total, clickable terms) with the
  `forced_damped_vibrations.html` explanation panel.
- **Aha:** divergence is what leaks out of a tiny box per unit volume; Gauss' theorem just adds the boxes up, and every
  interior face cancels.

### E5 · stokes_circulation_loop
- **A:** C16, C11 (also shows N69 orientation, N71 (2.35) as the limit view, N73 Ex. 2.6 rectangle bookkeeping, N54
  irrotational fields, N55 Ex. 2.3's b × x, N53 the curl components; previews ch03 "vorticity = 2 × angular velocity"
  and ch05 Kelvin)
- **Confusion removed:** "curl means the flow goes round in circles" (a straight shear flow has curl; an irrotational
  vortex goes round with zero curl everywhere except its core) and "circulation is about the loop" (it is the total
  curl inside the loop — as long as the field is smooth there).
- **Why interactive:** the theorem equates a loop quantity with a surface quantity; dragging and resizing the loop over
  different fields, flipping its orientation, and shrinking it until Γ/A → (∇×u)·n is the experiment that turns
  (2.34)–(2.35) into a picture. The irrotational-vortex preset with the loop enclosing vs excluding the core shows the
  hypothesis failing in a way no static figure does.
- **Stage:** (1) a 2-D velocity field as a heatmap of (∇×u)₃ (purple positive / green negative) with streamlines and
  tracer particles, a draggable/resizable loop (circle or rectangle) with u·t arrows along it and a paddle wheel at its
  centre turning at ½(∇×u)₃; (2) u·t vs arc length s with the signed area shaded = circulation; (3) bars: ∮u·t ds vs
  ∬(∇×u)·n dA, and the limit view Γ/A vs loop size → (∇×u)₃(x₀) (hidePortrait; the number repeated in view 2's title).
- **Controls:** field preset (solid-body rotation b × x — Ex. 2.3, simple shear u₁ = Γx₂, irrotational vortex
  u_θ = K/r, potential flow ∇φ (∮ = 0, Exercise 2.20), source, a smooth sin/cos field) · loop centre (drag) · loop radius
  or rectangle sides · loop shape (circle / rectangle — Ex. 2.6) · orientation toggle (flip n: both sides change sign) ·
  optional transport for the tracer particles and the wheel.
- **Equations shown (live):** (2.34), (2.35), t = n_c × n, (2.24)/(2.25) for the z-component, Ex. 2.6's four-side sum,
  Γ_circ = 2\|b\|πR² for the solid-body preset, ∮∇φ·t ds = 0.
- **Mirrors:** `core.integral_theorems.circulation`, `curl_flux`, `stokes_theorem_check`, `integral_curl_component`,
  `planar_loop`, `rectangle_loop`, `boundary_tangent`, `core.operators.curl`, `ch02.solid_body_rotation_field`,
  `ch02.shear_field`, `ch02.irrotational_vortex_field` (§8).
- **Derivations:** D26 (the rectangle step sets loop shape = rectangle and lights each side with its u·t value and
  sign; the tiling step splits the loop's interior into tiles with cancelling shared edges; the corollary step sets the
  potential-flow preset), D12 (the z-component step shows ∂u₂/∂x₁ − ∂u₁/∂x₂ evaluated at the centre).
- **Depth features:** Explain tab (the field at the loop centre → (∇×u)₃ by (2.25) with numbers → the four side sums
  (rectangle) or the midpoint sum (circle) → the circulation → the curl flux → their difference → Γ/A vs the point
  value → reading the current setting), synced Code tab, + linked views (3), presets, status verdict ("🌀 positive
  circulation: counterclockwise net spin" / "↔ shear flow: straight streamlines, yet Γ ≠ 0" / "0️⃣ irrotational: Γ → 0
  as the loop shrinks" / "⚠️ singular core inside the loop: Stokes does not apply — Γ = 2πK regardless of size"), term
  bars (the four sides), inspector (click the loop: u, t, u·t, ds at that point), transport (tracers + wheel).
- **Follows reference:** `templates/viz_example_field.html` field stage (heatmap + streamlines + tracers + probe) with
  the `forced_damped_vibrations.html` explanation panel; wheel idea from `angular_frequency_explorer_1.html`.
- **Aha:** circulation is the total spin inside the loop, not the shape of the streamlines — a straight shear flow has
  it, and a whirlpool with a singular core has it only because of the core.

### B1 · index_machine (backup)
- **A:** C01, C08 (also shows N12 free vs dummy, N21/N23 the matrix product (2.9)/(2.11), N41 δ substitution, N44 ε–δ,
  N49 the cross product (2.21), N74 comma notation)
- **Confusion removed:** index expressions read as gibberish; which index is summed; why renaming a dummy is free and
  renaming a free index is not; what ε_ijk "does".
- **Why interactive:** choosing an expression from chips (a_i b_i, A_ik B_kj, δ_ij u_j, ε_ijk u_i v_j, ε_ijk ε_klm,
  C_im C_jn τ_ij, u_i,i) and watching the summation convention unfold it into its 3, 9 or 27 explicit terms — each term
  lighting the matrix/ε-cube cells it multiplies — turns notation into mechanics.
- **Stage:** (1) the expression with its indices coloured (free blue, dummy orange) and the expanded sum below; (2) the
  ε cube (3 × 3 × 3, Viz.three or an unfolded 3 × 9 grid) and the operand matrices with the active cells lit as the sum
  is stepped; (3) the numeric result for the reader's operand values.
- **Controls:** expression chips · step through the terms (transport) · operand values (small sliders or presets) · dummy
  letter (rename i → k live) · dim = 2 / 3 toggle.
- **Equations shown:** (2.2), (2.9), (2.17), (2.18), (2.19), (2.21), (2.12), (2.36).
- **Mirrors:** `core.index_notation.expand_indices`, `classify_indices`, `tensor_order`, `core.tensors.levi_civita`,
  `epsilon_delta_residual`, `cross_einsum`, `double_dot`.
- **Derivations:** D09 (the case enumeration, stepping the ε cube).
- **Depth features:** Explain tab, synced Code tab, + transport (term by term), presets, inspector (click a term), modes
  (dim 2 / 3), status ("order 2 tensor: 2 free indices" / "❌ an index appears three times").
- **Follows reference:** `stride_padding_playground.html` (a formula with numbers plugged in, badges, click a cell) and E5.
- **Aha:** a repeated index is a loop; the free indices are the shape of the answer.

## 6. Python animations and interactive figures (A ID → what, why, player/figure kind)
A items get their own visuals; B items reuse the parent's figure as an overlay or panel. FAST budget: grids 24³ (full
48³), sphere quadrature 24 × 48 (full 48 × 96), convergence sequences n = 8, 16, 32 (full + 64); slider figures
precomputed (≤ 40 steps × ≤ 4 traces × ≤ 400 points); animations ≤ 60 frames at dpi 80.

| Kind | CORE | What moves / what the slider controls | Why | Player / figure |
|---|---|---|---|---|
| Animation | C02, C03 (N09, N17, N16) | the frame O1'2' rotating through 0 → 90° while the arrow x stays fixed; component bars x_i (teal) and x'_j (orange) update; the live C matrix in a corner with CᵀC = I | "the vector does not change, its components do" is motion; the E1 idea for the page | `player="frames"` (step 10°) |
| Animation | C10, C11 (N55) | tracer particles in a x (spreading from the origin) and in b × x (turning as a rigid disc), side by side, with div and curl printed | divergence = spreading, curl = spinning, seen rather than computed | `player="video"` |
| Animation | C12 (N57, N58, N62, N63) | a square of tracers under G, S, A in three panels on one clock for the simple-shear preset (stretch along 45° + spin) and for solid-body rotation | the decomposition is a statement about motion; the E3 idea for the page | `player="video"` |
| Animation | C15 (N68) | Ex. 2.5's six-face bookkeeping: face values appear, the zeroth-order terms of opposite faces cancel (fade), the first-order terms remain, the box shrinks and (1/V)∮ converges | the limit (2.32) is a process with a cancellation inside it | `player="frames"` |
| Animation | C16 (N71) | a square loop shrinking around a point in the shear field u₁ = Γx₂; the running Γ/A next to (∇×u)₃; a paddle wheel spinning at half the curl | curl as circulation per unit area is a limit | `player="frames"` |
| Plotly slider | C05 (N38, N39) | φ (0–180°, 37 steps): channel element with the normal, the traction arrow f = n·τ, its σ_n and τ_s parts; traces for a > 0 and a < 0 | traction depends on the cut direction; the double-angle structure appears | `slider_figure` |
| Plotly slider | C06 (N27, N38) | θ (0–180°): τ'₁₁, τ'₁₂, τ'₂₂ vs θ for the pure-shear and uniaxial presets, with the Mohr circle drawn as the second panel | (2.12) as a curve: τ'₁₂ = 0 at 45° links to C13 | `slider_figure` |
| Plotly slider | C09 (N51) | direction angle of n (36 steps): contour plot of φ = x² + y²/4, the ∇φ quiver, the n arrow at a probe point and ∂φ/∂n = ∇φ·n printed in the trace name | maximum when n ∥ ∇φ, zero along the contour | `slider_figure` |
| Plotly slider | C13 (N62, N63) | Γ (or the pure-shear amplitude, 30 steps): the deforming square with the eigen-axes, the eigenvalues ±Γ and Mohr's circle for S | principal axes stay at 45° while the rates scale with Γ | `slider_figure` |
| Plotly slider | C14 (N64, N65) | mesh count n = 4 … 32 (8 steps): bars ∫∇·Q dV vs ∮n·Q dA for a smooth Q on the unit cube and their gap on a log axis; the 3-D box with faces coloured by n·Q as the second panel | both sides converge to the same number as the quadrature refines | `slider_figure` |
| Plotly slider | C15 (N66) | box side h (20 log steps): (1/V)∮n·Q dA and (1/V)∮n Q dA (gradient) vs the point values, error vs h on log axes | the (2.31)–(2.32) limit with its order-2 error | `slider_figure` |
| Plotly 3-D | C04 (N25), C05 (N36), C14 (N65), C16 (N70) | the stress cube with nine `go.Cone` arrows and a toggle for the hidden faces; the tetrahedron with the coordinate-face stresses and f(n) on the slanted face; the box with faces coloured by n·Q; the cap surface with n, n_c, t at a boundary point | 3-D geometry that a 2-D sketch flattens | `go.Cone`, `go.Scatter3d`, `go.Mesh3d`, height 520 |
| Static figures | C01 (N23), C07, C08, C09–C11, C15, C16 | (2.11) row/column highlight of P₁₂; invariants I₁, I₂, I₃ flat across 50 random rotations while the entries scatter; ε as a 3 × 3 × 3 coloured cube (unfolded); heatmaps of ∇·u and (∇×u)₃ for the Ex. 2.3 fields; log–log convergence of `gradient`/`divergence`/`curl` and of `integral_divergence`, `integral_curl_component` (slope 2) | one relationship each | matplotlib |
| Live widget (kernel only; paired with the C13 slider figure) | C13 | S₁₁, S₂₂, S₁₂ sliders → eigen-axes on the element, Mohr's circle, λ and the principal angle | free exploration beyond the precomputed figure | `live(...)` |

Optional cheap overlays for B items: the N18 residual table as a bar chart; N42 isotropy residuals (δ, ε, random) as
three bars; N44's 81-case check as a 9 × 9 heatmap of the residual (all zero).

## 7. From-scratch moments
A items only; every section with a computable A item has at least one; each is followed by
`assert np.allclose(mine, library)`.

| § | CORE | Hand-written version | Tested function it must agree with |
|---|---|---|---|
| 2.1 | C01 | explicit `for i in range(3)` loop for a_i b_i and a double loop with an inner `for k` for A_ik B_kj | `core.tensors.dot`, `inner` (and `np.einsum`) |
| 2.2 | C02 | build C entry by entry as `E_old[i] @ E_new[j]` (nine dot products) for a frame rotated by θ | `core.tensors.direction_cosines`, `rotation_matrix_2d` |
| 2.2 | C03 | `xp[j] = sum(x[i] * C[i, j] for i in range(3))` | `core.tensors.transform_vector` |
| 2.4 | C04 | list the six faces ±e_i and read the traction ±τ[i, :] by the sign rule | `core.tensors.cube_face_tractions` |
| 2.6 | C05 | `f[i] = sum(tau[j, i] * n[j] for j in range(3))` — and the *wrong* index order shown to differ for a non-symmetric τ | `core.tensors.traction` |
| 2.4 | C06 | double loop `tp[m, n] = sum(C[i, m] * C[j, n] * tau[i, j] …)` next to `C.T @ tau @ C` | `core.tensors.transform_tensor` |
| 2.5 | C07 | I₁ as a loop over the diagonal, I₂ = ½(I₁² − Σ_ij A_ij A_ji), I₃ by cofactor expansion (P53) | `core.tensors.trace`, `invariants` |
| 2.7 | C08 | fill a 3 × 3 × 3 array by the parity of (i, j, k) with `itertools.permutations`; the 81-case ε–δ check as a quadruple loop | `core.tensors.levi_civita`, `epsilon_delta_residual` |
| 2.9 | C09 | central difference (φ[j, i+1] − φ[j, i−1])/(2h) on the grid interior, one-sided second order at the edges | `core.operators.partial`, `gradient` |
| 2.9 | C10 | sum of the three hand-written partials of u₁, u₂, u₃ | `core.operators.divergence` |
| 2.9 | C11 | the three components of (2.25) from hand-written partials | `core.operators.curl`, `curl_components` |
| 2.10 | C12 | `0.5 * (B + B.T)`, `0.5 * (B − B.T)` and reading ω from A's three entries | `core.tensors.symmetric_part`, `antisymmetric_part`, `vector_from_antisymmetric` |
| 2.11 | C13 | 2 × 2 characteristic polynomial, quadratic formula, eigenvector from one row of (S − λI), normalise, order and sign-fix for det C = +1 | `core.tensors.principal_axes`, `diagonalize`, `ch02.example_2_4` |
| 2.12 | C14 | six-face midpoint sums of n·Q on the unit cube vs a 3-D midpoint sum of ∇·Q | `core.integral_theorems.divergence_theorem_box` |
| 2.12 | C15 | the six-face sum on a small cube of side h around x₀ divided by h³ for three h values, with the observed order | `core.integral_theorems.integral_divergence` (+ `tools/convergence.observed_order`) |
| 2.13 | C16 | midpoint sum of u·t Δs around a circle parametrised by hand vs a polar-grid sum of (∇×u)₃ over the disc | `core.integral_theorems.circulation`, `curl_flux` |

(§2.3, §2.8 and §2.14 have no A item; their B items are exercised inside C01, C08 and C10.)

## 8. Notes for the implementer
Analysis §4 (F1–F26) covers the A items' physics. Additions the figures, from-scratch cells and explainers need:

1. **2-D helpers for E2 and the C05/C06 slider figures** (scalar-callable for parity rows): `ch02.traction_2d(tau, phi)`
   → dict(n, f, sigma_n, tau_s) for any symmetric 2 × 2 τ (Ex. 2.2 is the pure-shear case); `ch02.stress_vs_angle(tau,
   phi_array)`; `core.tensors.mohr_circle_2d(tau) -> (center, radius)`; `ch02.principal_angle_2d` (already F22).
   `traction` must accept 2 × 2 as well as 3 × 3 (the same einsum).
2. **Linear-flow kinematics for E3 and the C12 animation**: `ch02.linear_flow_map(G, t)` = `scipy.linalg.expm(G t)` and
   `ch02.deform_square(G, t, n_side=10)` → tracer positions; `ch02.velocity_gradient_preset(name)` for simple shear,
   solid-body rotation, pure strain, uniaxial extension (2-D and 3-D forms); `ch02.material_line_angle(G, t, theta0)`.
   Pin the sign convention with a test: for u = b × x, `vector_from_antisymmetric(antisymmetric_part(vector_gradient(u)))`
   equals b (so the vector of A is ½∇×u), with G[i, j] = ∂u_i/∂x_j.
3. **2-D versions of the integral theorems for E4/E5 parity and the 2-D notebook figures** (the 3-D functions stay the
   reference): `core.integral_theorems.divergence_theorem_rect2d(Q_fn, bounds, n) -> (lhs, rhs, face_fluxes)`,
   `flux_through_faces(Q_fn, bounds, n) -> dict[face, flux]`, `integral_divergence_2d(Q_fn, x0, h, n_side)`,
   `circulation` and `curl_flux` accepting 2-D fields (z-component); `rectangle_loop` and `planar_disc` in 2-D.
4. **Field presets as callables with sympy twins** (`ch02.exact_div_curl` already F17): `ch02.shear_field(Gamma)`,
   `ch02.irrotational_vortex_field(K)` (u_θ = K/r; flag `singular_at=(0,0)` so the E5 status and a test can detect a
   loop enclosing the core: circulation = 2πK independent of radius, `stokes_theorem_check` reports the hypothesis
   failure instead of a false mismatch), `ch02.potential_field(phi_expr)`, `ch02.point_source_field(m)`,
   `ch02.smooth_test_field()` (sin/cos products for the convergence tests).
5. **Grid and stencil**: fix the layout `[k, j, i]` = (z, y, x) in `core.grids` once, with `grid2d` using
   `indexing='xy'`; `core.operators.partial` second-order one-sided at the edges (never `np.gradient`, ch01 lesson);
   vector fields on axis 0. Expose `laplacian` now (ch04+). FAST: 24³ grids; the operator convergence test uses n = 8,
   16, 32 in FAST.
6. **`principal_axes` raises for non-symmetric input** (analysis §9); returns λ ascending, B with det +1, and a stable
   result for repeated eigenvalues (test with projectors, not eigenvectors). `example_2_4(Gamma)` takes Γ as S₁₂ and its
   docstring states the ½ du₁/dx₂ factor for a u₁(x₂) flow.
7. **Convention flags in signatures**: `double_dot(A, B, convention="book"|"frobenius")` with no silent default in the
   notebook text; `tensor_divergence(T, h, index=1)`; `rotation_matrix_2d/3d` docstrings say "passive: columns = new
   axes; x' = Cᵀx; equals the active R(θ) of Wikipedia/scipy applied as Rᵀ"; `is_isotropic(..., proper=True)`.
8. **Promotion**: move `ch01.traction_components` into `core/tensors.py` (`normal_shear_stress` calls it) with a
   re-export from `ch01_introduction` so ch01 tests keep passing.
9. **Book values** (private, git-ignored `tests/book_values_ch02.json`): keys per analysis §7 for Examples 2.1–2.4 and
   Exercise 2.1a; `reference/ch02/benchmarks.json` holds the public constants (sphere flux 8π/3, ε contractions 2, 6,
   3) and `SOURCES.md` the Wikipedia citations from analysis §8.
10. **Explainer parity functions must take plain floats/lists**: `rotation_matrix_2d`, `transform_vector`,
    `transform_tensor`, `orthogonality_residual`, `polar_components`, `traction`, `normal_shear_stress`,
    `traction_2d`, `principal_axes`, `example_2_2`, `example_2_4`, `symmetric_part`, `antisymmetric_part`,
    `vector_from_antisymmetric`, `linear_flow_map`, `divergence_theorem_rect2d`, `integral_divergence_2d`,
    `circulation`, `curl_flux`, `integral_curl_component`, `levi_civita`, `epsilon_delta_residual`, `expand_indices`
    (string result for B1).
11. **Notebook budget**: every 3-D plotly figure ≤ 300 kB (few cones, coarse meshes); animations ≤ 60 frames; the two
    convergence studies cached in one cell and reused; total target < 4 min on Colab CPU with FAST.
