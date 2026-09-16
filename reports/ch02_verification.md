# Chapter 2 verification — Cartesian Tensors                     2026-09-16, commit ebc80e3 + review round (M1/M2 code fixes re-verified)

Verifier: math-verifier. Suite: `tests/test_ch02.py` (82 test functions, 92 collected items), figures/metrics:
`tests/ch02_verify_figures.py` → `outputs/ch02/verify/` (git-ignored), cited data: `reference/ch02/` (`make_refs.py`,
`benchmarks.json`, `SOURCES.md`). Nothing in `fluidpy/` or `scripts/` was edited by the verifier.

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Windows, `.venv`).

Runs (after the review round): `pytest tests/test_ch02.py -q -p no:cacheprovider` → **92 passed, 0 failed, 0 xfailed in 37 s**
(with the private book file; without it 4 V6 tests skip). Full suite `pytest -q` → **206 passed, 1 failed in 30 s** — the
one failure is `tests/test_machinery.py::test_viz_library_inlined_and_template_lints` ("STALE viz/ch02/rotation_of_axes.html":
the viz-builders are mid-build and `viz/` is outside the verifier's scope; every ch01/ch02 physics test passes; ch01's 107
tests still pass with `traction_components` re-exported from `core.tensors`). First round (before the fixes): 91 passed,
1 strict xfail; full suite 206 passed, 1 xfailed. All 14 runnable
`scripts/ch02_*.py` exit 0 (`ch02_drawings.py` is an imported helper module without `__main__`, Open item O5).
`tools/check_public.py` → OK (165 files).

Collected items per evidence tag: **V1 41 · V2 13 · V3 15 · V4 1 · V5 6 · V6 4 · V7 12** (= 92).

**Discrimination proofs** (skill lesson ch01): eight wrong variants were monkey-patched into the live modules (scratchpad
script) and the convention tests re-run — every variant fails at least one test:
τ_ij n_j instead of τ_ji n_j (2 tests fail) · active C x instead of Cᵀx (4) · R_ij = +ε_ijk ω_k (1) · book/Frobenius
double dot swapped (1) · G[i, j] transposed (2) · tensor divergence on the first index (1) · clockwise loops (2) ·
`np.gradient` (first-order edges) in place of the one-sided stencil (4 of 5 operator-order tests).

## Validation table — A items (CORE, ≥ 2 independent levels, one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
|---|---|---|---|---|---|---|
| C01 summation convention (2.2), free/dummy indices, (2.9)–(2.11), (2.16)–(2.17), (2.36) | CORE | `expand_indices`, `expand_indices_str`, `classify_indices`, `tensor_order`, `rename_dummy`, `comma_to_partial`, `dot`, `inner`, `kronecker_delta` | V1 a_i b_i expansion (symbolic + numeric 32), exact string, free/dummy/rename, three-repeat error, δ substitution (δ_ii = 3), P_12 of (2.11); V2 `"C_im C_jn tau_ij"` = (CᵀτC)_mn all 9, `"eps_ijk u_i v_j"` = cross, `"u_i,i"` = div, `"eps_ijk u_k,j"` = (2.25), `"-1/2 eps_ijk R_ij"` = (2.27) | all exact (sympy `expand == 0`) | analytic, symbolic | |
| C02 direction cosines C_ij = e_i·e'_j; N14 (2.7), N15 orthogonality; D02, D03 | CORE | `direction_cosines`, `rotation_matrix_2d/3d`, `rotation_angle`, `random_rotation`, `is_orthogonal`, `is_proper_rotation`, `orthogonality_residual`, `inverse_transform_vector` | V1 columns = new axes, R_z embedding, composition, angle round trip, axis eigenvector; V5 `rotation_matrix_3d` == scipy `Rotation.from_rotvec` (20 random axes) and the passive-vs-active discrimination; V7 200 random rotations orthogonal/proper, reflections/scaled/random rejected; V2 D02 all 8 steps for a general Euler rotation | entries 1e-15; scipy 1e-13; residual < 1e-14; det − 1 < 1e-12 | analytic, symbolic, benchmark (library) | the book's C *equals* R(θ) and is applied as Rᵀ |
| C03 transformation of components (2.5), (2.8); N16, N18, N19 (Ex. 2.1); D01 | CORE | `transform_vector`, `transform_tensor` (order 1), `transforms_as_vector`, `polar_components`, `cartesian_from_polar`, `polar_components_deg`, `vector_from_components` | V1 = Cᵀx and the hand loop, |x'| = |x|, passive sign (old e₁ has a negative 2'-component after +30°), polar == rotation for 13 angles, round trips; V2 D01 x·e'_j = x_i C_ij and completeness; V7 pass/fail table (x, 2.5x, b×x, invariant ∇φ pass; x², (|x|,0,0), a frame-specific ∇ formula fail); V6 Ex. 2.1 | 1e-14; residuals < 1e-13 vs > 1e-2 | analytic, symbolic, book-value | |
| C04 stress tensor τ_ij and sign convention; N25, N26, R01 | CORE | `cube_face_tractions`, `stress_component_meaning`, `STRESS_CUBE_FACES` | V1 face +2 = τ[1,:], −2 = −τ[1,:], every face == `traction(τ, n)`, Ex. 2.5 lettering, meaning text; V1 six tractions sum to 0 and opposite faces cancel (R01); V7 face tractions transform as vectors under 10 random rotations | 1e-14; 1e-13 | analytic | |
| C05 Cauchy's traction (2.15); N34–N38 (Ex. 2.2); D05 | CORE | `traction`, `traction_components` (promoted), `normal_shear_stress`, `tetrahedron_face_areas`, `shear_flow_stress`, `traction_2d`, `example_2_2`, `stress_vs_angle`, `mohr_circle_2d` | V1 first-index discrimination (n@τ == τᵀn ≠ τn, gap > 0.1), pressure f = −pn, 2×2/lists, n = e₃ returns row 3; V2 D05 tetrahedron balance solved symbolically and h → 0 limit, τ_ji n_j ≠ τ_ij n_j symbolically, vector-area closure; V5 Wikipedia T_j = σ_ij n_i index placement; V7 covariance f' = Cᵀf for non-symmetric τ; V1 Ex. 2.2 general forms a sin 2φ / a cos 2φ over 50 φ, two routes agree, extremes = eigenvalues; V5 Mohr formula; V6 Ex. 2.2 | 1e-13/1e-14; V6 1e-12 | analytic, symbolic, benchmark, book-value | |
| C06 tensor transformation rule (2.12), (2.13); N27–N29; D06 | CORE | `transform_tensor`, `transforms_as_tensor`, `tensor_product`, `outer` | V1 order 2 = CᵀτC (≠ CτCᵀ), order 4 vs four-fold loop, product rule (2.13); V2 D06 steps 1–8 symbolically (vector rule + (2.15) + (2.7) → τ' = CᵀτC, coefficient matrix vanishes, dummy rename, C = I limit); V7 pass/fail (u_i v_j, x_i x_j, δ pass; a fixed array fails) | 1e-12–1e-14; fail residual > 1e-2 | analytic, symbolic | |
| C07 contraction, invariants, double dot; N30–N33; D18 | CORE | `contract`, `dot_tensor_vector`, `double_dot`, `trace`, `invariants`, `characteristic_polynomial` | V1 four (2.14) patterns, book == tr(AB) vs Frobenius == tr(ABᵀ) (differ unless symmetric), Vieta on EXAMPLE_TENSOR (9, 24, 18), roots == eigvalsh, 2×2 case; V7 invariants/trace/double dots unchanged under 100 rotations; V2 D18 I₁, chain, det invariant (symbolic z-rotation), cubic −λ³ + I₁λ² − I₂λ + I₃ exact for a general 3×3, M = Σ principal minors = I₂, Vieta | 1e-12; polynomial identities exact | analytic, symbolic | see Open item O3 (A_ij A_ij *is* invariant) |
| C08 alternating tensor (2.18), ε–δ (2.19), cross product (2.20)–(2.21); N42–N50; D09 | CORE | `levi_civita`, `permutation_sign`, `epsilon_delta_residual`, `triple_product`, `cross`, `cross_einsum`, `angle_between`, `is_isotropic` | V1 27 values (0/1-based), index moves, residual = 0, 2δ and 6, triple product; V2 D09 81 cases by sympy with the step 2–7 case split, contractions via the expander, BAC−CAB symbolically; V5 Wikipedia identities (a)–(e) + pseudotensor sign flip; V1 cross on 1000 pairs (Lagrange identity, right-handedness, angle stability); V7 isotropy of δ, λδ, ε (proper) and non-isotropy under reflection (residual 2) | exact / 1e-13 | analytic, symbolic, benchmark | |
| C09 gradient (2.22), ∂φ/∂n = ∇φ·n; N51; R02 | CORE | `partial`, `gradient`, `directional_derivative`, `grid`, `grid2d`, `evaluate_field` | V1 exact on quadratics including one-sided edges (1e-12), ∇φ ⊥ level curves of x² + y²/4, |∂φ/∂n| ≤ |∇φ|, extremes along/across ∇φ; V3 one-sided order 1.970 (pairwise 1.945, 1.976, 1.989), periodic 1.997; V7 rotational invariance (rotated field's gradient = Cᵀ∇φ, order 2.0) | see convergence table | analytic, converged | |
| C10 divergence (2.23); N52 (∂u_i/∂x_j, ∇·τ), N55 (Ex. 2.3), N74 | CORE | `divergence`, `vector_gradient`, `tensor_divergence`, `is_solenoidal`, `is_irrotational`, `radial_field`, `solid_body_rotation_field`, `example_2_3`, `exact_div_curl` | V1 ∇·(a x) = 3a and ∇·(b×x) = 0 exactly on every node, plane 2a; V2 symbolic 3a / 0 / 2b with symbolic a, b (and the a x₃ e₃ correction); V1 G[i, j] = ∂u_i/∂x_j pinned (u = (x₂, 0, 0)), div = tr G, tensor divergence index 1 vs 0 on a linear tensor field vs A_ijj / A_jij, quadratic outer-product field vs sympy; V3 one-sided 1.997, periodic 1.997; V6 factors 3 and 2 | 1e-12; orders in table | analytic, symbolic, converged, book-value | |
| C11 curl (2.24)–(2.25); N53, N54; D12 | CORE | `curl`, `curl_components`, `potential_field` | V1 curl == (2.25) components (1e-14), 3-component field on a 2-D grid, plane scalar curl, ∇×(b×x) = 2b exactly; V2 D12 via the expander; V1 ∇×∇φ = 0 and ∇·(∇×u) = 0 to round-off (discrete partials commute) on one-sided and periodic grids; V3 one-sided 2.029 (pairwise 2.048, 2.026, 2.013), periodic 1.997; V7 rotational invariance (curl u' = Cᵀ curl u, order 2.0) and mirror antisymmetry | 1e-11 identities; orders in table | analytic, symbolic, converged | |
| C12 symmetric + antisymmetric split; N56–N61 (2.26)–(2.29); D14, D15 | CORE | `symmetric_part`, `antisymmetric_part`, `strain_rate_tensor`, `rotation_tensor`, `independent_components`, `antisymmetric_from_vector`, `vector_from_antisymmetric`, `symmetric_double_contraction`, presets, `linear_flow_map`, `deform_square`, `material_line_angle` | V1 S + A = B, counts 6/3/9, uniqueness (sympy solve), stacked fields; V7 parts and the hidden vector transform under rotation; V1 (2.26) entries, R·x = ω×x (negated variant ≠), round trip, ½∇×u for b×x (pins G order + sign), plane scalar case; V2 D15 steps 1–7 symbolically (antisymmetry, entries, −2ω_l contraction, inversion with i = 1, cross product; opposite sign gives −ω×x); V1 P_A = 0 for symmetric τ, ≠ 0 otherwise; V1 kinematics (rotation map = C(ωt), pure strain diag(e^{±Γt}), uniaxial diag(e^{Γt}, 1), d/dt = G, area preserved for traceless G and grows by e^{Γt} for uniaxial, S₁₂ = Γ/2 for simple shear, traces Γ / 0 and spin content of the five presets); **V1 (review M1) `rotation_tensor(G)` = G − Gᵀ = 2A, ≠ A; for u = b × x: vector(R) = 2b = ∇×u, vector(A) = b; R = −ε·(∇×u); G = S + ½R** | 1e-13–1e-15 | analytic, symbolic | review M1/M2 re-verified from p082 |
| C13 principal axes of a symmetric tensor; N62 (Ex. 2.4), N63; D17 ★★★ | CORE | `principal_axes`, `diagonalize`, `normal_stress_bounds`, `example_2_4`, `principal_angle_2d` | V1 τB = Bdiag(λ), BᵀB = I, det +1, τ' diagonal (20 random), EXAMPLE_TENSOR λ = 3 ∓ √3, 3, hydrostatic, repeated λ via projectors, non-symmetric refused; V5 2-D Mohr formula (20 random) + τ_max = ½|σ₁ − σ₂|, (σ_n, τ_s) on the circle; V7 σ_n ∈ [λ_min, λ_max], τ_s ≤ Δλ/2 (4000 normals), σ_n = Σλ_k c_k²; V2 D17 steps 1–5 (im(b̄τb) = 0 for a general real symmetric 3×3 — and ≠ 0 without symmetry), 7, 5/8/9/11/13/14/15 on the symbolic 2×2 (discriminant, Vieta, orthogonality, CᵀC = I, CᵀτC = diag(λ), σ_n = Σλ_k c_k², τ_s² = Δλ² sin²2φ/4); V1 Ex. 2.4 with Γ = 1, 2.5, −0.4 (λ = (Γ, −Γ), 45°, S' diag, det C = +1); V6 Ex. 2.4 | 1e-12–1e-14 | analytic, symbolic, benchmark, book-value | Open item O4 (design's sympy cell uses `.norm()`) |
| C14 Gauss' theorem (2.30); N64, N65; D25 | CORE | `divergence_theorem_box`, `gauss_gradient_box` (scalar and vector Q), `flux_through_box_faces`, `volume_integral_box`, `divergence_theorem_sphere`, `flux_through_sphere`, `divergence_theorem_rect2d`, `flux_through_faces`, `divergence_theorem_tiled`, quadrature nodes | V1 (x, y, z) → 3 = 3, (x², 0, 0) → 1 = 1, face keys/sum, scalar xyz gradient form ¼, **vector Q = (xy, yz, zx): both (3, 3) sides equal the exact ∫∂Q_j/∂x_i dV to 1e-12 (FD) / 1e-13 (exact dQ), index order = transposed `vector_gradient`**, 2-D 2 = 2 and 1 = 1, Simpson/Gauss node exactness; V3 midpoint 2.002, Simpson 4.00, 2-D 2.0, FD-divergence fallback 1e-10; V5 F = (2x, y², z²) through the unit sphere = 8π/3 (both sides rel 4e-16 / 8e-16; FD fallback 1e-9; shifted centre 8π/3(1 + c_y + c_z)); V4 solenoidal b×x: zero net outflux through boxes, spheres, rectangles and tiles (< 1e-13); point source m through enclosing boxes = m, non-enclosing 0; V2 D25 FTC face pair for a symbolic polynomial Q (all three i), shared-face cancellation symbolically, numeric tiling interior −5.6e-17 | see tables | analytic, symbolic, converged, conserved, benchmark | O1 closed in the review round |
| C15 integral definitions (2.31)–(2.33); N66–N68 (Ex. 2.5), N72; D21, D22 | CORE | `integral_gradient`, `integral_divergence`, `integral_curl`, `integral_divergence_2d`, `integral_definition_convergence` | V1 exact for linear fields at h = 1, 0.3, 0.05 (3a, 2b, linear ∇, plane 2a), Q = (x², 0, 0) → 2.000 exactly, Q = (x³, 0, 0) → 3 + h²/4, index order [i, j] = ∂Q_j/∂x_i pinned; V3 orders 1.998 / 1.997 / 1.998 / 1.998; V1 agree with the grid operators at a node (3e-3 at h = 0.05); V2 D21/D22 symbolic box outflux for a symbolic cubic Q: zeroth-order terms cancel per face pair, ratio → ∂Q_i/∂x_i, remainder exactly c₅Δ₁²/4 + c₈Δ₃²/4, limit 0; n·Q and (n×Q)_k index forms; mean-value theorem model | exact / orders | analytic, symbolic, converged | |
| C16 Stokes' theorem (2.34)–(2.35); N69–N71, N73 (Ex. 2.6); D26 | CORE | `boundary_tangent`, `planar_loop`, `rectangle_loop`, `planar_disc`, `planar_rectangle`, `circulation`, `curl_flux`, `stokes_theorem_check`, `integral_curl_component`, `irrotational_vortex_field`, `shear_field` | V1 t = n_c × n with n_c into A (Fig. 2.10) is right-handed, ∮(x − c)×t ds = 2A n for circles and rectangles (10 random planes), t = n × r̂, tangents = derivative of points, disc/rectangle areas, 2-D convenience; V1 b×x: circulation = curl flux = 2|b|πR² (rel 7e-15) for R = 0.5, 1, 2, FD-curl fallback 1e-8, flipped n flips both sides, plane 2πb₃, shear −Γ; V1 ∮∇φ·t ds = 0 (2-D and 3-D, < 1e-10); V3 rectangle loop 2.000 (2-D), 3-D tilted rectangle 2.0 ± 0.25; V1 (2.35) exact for b×x in three orientations, Ex. 2.6 u_y discrimination (u = (0, z, 0) → −1, u = (0, 0, y) → +1), 2-D form, order 2 for a smooth field; V7 singular vortex: hypothesis flag False, Γ = 2πK for every R, curl flux 0; non-enclosing loop OK; V2 D26 one rectangle exactly (steps 1, 3–6 signs, exact Stokes for polynomial u, (2.35) limit, shared-edge cancellation, corollary ∮∇φ·t = 0) | 1e-12–1e-15; orders | analytic, symbolic, converged | Open item O2 (docstring wording "outward") |

Every computable A item has ≥ 2 independent levels with at least one of V1/V2/V3/V5.

## Validation table — coded B/C items (NOTE, ≥ 1 level)
| Concept | fluidpy target | Evidence (test) | Label |
|---|---|---|---|
| N03/N10 basis form (2.1), (2.3) | `unit_vectors`, `vector_from_components` | V1 `test_inverse_transform_V1_round_trip`, `test_direction_cosines_V1_columns_are_new_axes` | analytic |
| N07 dot product (2.2) | `dot` | V1 `test_summation_convention_V1_dot_product_expansion` (32; float out) | analytic |
| N08 order = free indices | `tensor_order` | V1 `test_summation_convention_V1_free_dummy_and_errors` | analytic |
| N12/N13 free vs dummy, (2.6) | `classify_indices`, `rename_dummy` | V1 same test | analytic |
| N14 inverse (2.7) | `inverse_transform_vector` | V1 round trip (20 rotations), V2 D03 in `test_orthogonality_V2_derivation` | analytic, symbolic |
| N15 orthogonality, det +1 | `is_orthogonal`, `is_proper_rotation`, `orthogonality_residual`, `random_rotation` | V7 200 draws + counterexamples; V2 D02 | analytic, symbolic |
| N16 matrix forms | `transform_vector` | V1 = Cᵀx, V5 passive vs scipy active | analytic, benchmark |
| N18 formal definition (2.8) | `transforms_as_vector` | V7 pass/fail | analytic |
| N19 Ex. 2.1 | `polar_components*`, `cartesian_from_polar` | V1 (13 angles), V6 | analytic, book-value |
| N21–N23 matrix product (2.9)–(2.11) | `inner`, `expand_indices("A_ik B_kj")` | V1 | analytic |
| N25/N26/R01 stress cube | `cube_face_tractions`, `STRESS_CUBE_FACES`, `stress_component_meaning` | V1 ×2, V7 | analytic |
| N27–N29 tensor ≠ matrix, (2.13), u_i v_j | `transforms_as_tensor`, `tensor_product`, `outer` | V1, V7 | analytic |
| N30–N33 contractions, A·u vs Aᵀ·u, double dot | `contract`, `dot_tensor_vector`, `double_dot`, `trace` | V1, V7 | analytic |
| N35 dA_i = n_i dA | `tetrahedron_face_areas` | V1 (in `test_cauchy_traction_V2_derivation`) | analytic |
| N38 Ex. 2.2 | `shear_flow_stress`, `example_2_2`, `traction_2d`, `stress_vs_angle`, `mohr_circle_2d` | V1, V5, V6 | analytic, benchmark, book-value |
| N40/N41 δ (2.16)–(2.17) | `kronecker_delta`, expander | V1 | analytic |
| N42 isotropic tensors | `is_isotropic` | V7 | analytic |
| N43/N44 ε moves, ε–δ (2.19) | `levi_civita`, `permutation_sign`, `epsilon_delta_residual`, `triple_product` | V1, V2 D09, V5 | analytic, symbolic, benchmark |
| N45 angle, trace of u_i v_j | `angle_between`, `trace(outer)` | V1 | analytic |
| N46–N50 cross product (2.20)–(2.21), k = 1 check | `cross`, `cross_einsum`, expander | V1 (1000 pairs), V2 | analytic, symbolic |
| N51 Fig. 2.7 | `gradient`, `directional_derivative` | V1 | analytic |
| N52 ∂u_i/∂x_j, (∇·τ)_i | `vector_gradient`, `tensor_divergence` | V1 (index pinned, discrimination), V3 (vector_gradient 1.970) | analytic, converged |
| C10 `divergence` of a 3-component z-independent field on a 2-D grid (review should-fix) | `divergence` | V1 linear field → 2 − 3 exactly, = tr(`vector_gradient`), = 2-component result; 2 components on a 3-D grid still refused | analytic |
| N53 (2.25) | `curl_components` | V1 (1e-14) | analytic |
| N54 solenoidal / irrotational | `is_solenoidal`, `is_irrotational` | V1 | analytic |
| N55 Ex. 2.3 | `radial_field`, `solid_body_rotation_field`, `example_2_3`, `exact_div_curl` | V1, V2, V6 | analytic, symbolic, book-value |
| N56–N58 S, A, (2.26)–(2.27) | `symmetric_part`, `antisymmetric_part`, `independent_components`, `antisymmetric_from_vector`, `vector_from_antisymmetric` | V1, V2 D15, V7 | analytic, symbolic |
| N59–N61 (2.28)–(2.29) | `symmetric_double_contraction` | V1 (+ discrimination) | analytic |
| N62 Ex. 2.4 | `example_2_4`, `principal_angle_2d` | V1, V6 | analytic, book-value |
| N64/N65 divergence theorem, Fig. 2.9 | `divergence_theorem_*`, `flux_through_*`, `divergence_theorem_tiled` | V1, V3, V4, V5 | analytic, converged, conserved, benchmark |
| N66–N68 (2.31)–(2.33), Ex. 2.5 | `integral_gradient/divergence/curl`, `integral_divergence_2d` | V1, V2 D22, V3 | analytic, symbolic, converged |
| N69/N70 orientation, Fig. 2.10 | `boundary_tangent`, `planar_loop`, `rectangle_loop`, `planar_disc`, `planar_rectangle` | V1 | analytic |
| N71/N73 (2.35), Ex. 2.6 | `integral_curl_component` | V1 (u_y discrimination), V3 1.998 | analytic, converged |
| N74 comma notation (2.36) | `comma_to_partial`, `expand_indices("u_i,i")` | V1, V2 | analytic, symbolic |
| grids (convention 8) | `grid`, `grid2d`, `axis_of_direction`, `spacing`, `evaluate_field` | V1 layout test | analytic |
| test fields (Part C 6.5) | `smooth_test_field`, `periodic_test_field`, `periodic_scalar_field`, `smooth_scalar_field`, `potential_field`, `point_source_field`, `irrotational_vortex_field`, `shear_field`, `VectorField`, `ScalarField` | V1/V2 smoke + periodicity + error paths (`test_part_c_contract_V1_…`), used throughout | analytic, symbolic |
| kinematics presets (Part C 6.6) | `velocity_gradient_preset`, `linear_flow_map`, `deform_square`, `material_line_angle` | V1 | analytic |
| convergence helpers (Part C 6.7) | `operator_convergence`, `integral_definition_convergence` | V3 (the studies themselves) | converged |

## Functions used by the notebook and explainers (design Part C) — test name each
`test_part_c_contract_V1_every_callable_exists_and_parity_rows_return_floats` asserts that all **143 Part C names**
exist in `ch02.__all__` and that the parity-row expressions of the design (`traction_2d(...)["sigma_n"]`,
`principal_axes(...)[0][1]`, `polar_components`, `double_dot`, `circulation`, `mohr_circle_2d`, `rotation_angle`,
`epsilon_delta_residual`, `integral_divergence`, `vector_from_antisymmetric`, `material_line_angle`) return Python
floats. Every name is additionally exercised with a physical assertion in the tables above (the table's "test"
columns); `symbol_tensor`, `field_functions`, `coordinates`, `periodic_*` fields and the `VectorField` error paths are
covered in the contract test itself. Machinery (`style`, `anim`, `interact`, `embed`) is covered by `tests/test_machinery.py`.

## Convergence studies (scheme | grids | observed order | design order)
| Scheme | Grids | Observed order (pairwise) | Design |
|---|---|---|---|
| `gradient`, one-sided edges, smooth φ on [−1, 1]³ | n = 16, 32, 64, 128 | **1.970** (1.945, 1.976, 1.989) | 2 |
| `divergence`, one-sided | same | **1.997** (1.993, 1.998, 2.000) | 2 |
| `curl`, one-sided | same | **2.029** (2.048, 2.026, 2.013) | 2 |
| `laplacian`, one-sided | same | **2.019** (2.030, 2.018, 2.010) | 2 |
| `vector_gradient`, one-sided | same | **1.970** (1.945, 1.976, 1.989) | 2 |
| `gradient` / `divergence` / `curl` / `laplacian`, `bc="periodic"`, `periodic_test_field` on [−π, π]³ (`grid(..., periodic=True)`) | n = 16 … 128 | **1.997 / 1.997 / 1.997 / 1.998** | 2 |
| `integral_divergence` (2.32) | h = 0.4, 0.2, 0.1, 0.05 | **1.998** (err 1.3e-2 → 2.1e-4) | 2 |
| `integral_curl` (2.33) | same | **1.997** | 2 |
| `integral_gradient` (2.31) | same | **1.998** | 2 |
| `integral_curl_component` (2.35), n = e₃ | same | **1.998** (7.7e-4 → 1.2e-5) | 2 |
| `integral_curl_component`, n = e₁ | h = 0.4, 0.2, 0.1 | 2.0 ± 0.15 | 2 |
| `divergence_theorem_box`, midpoint | n = 8 … 64 | **2.002** | 2 |
| `divergence_theorem_box`, Simpson | n = 9 … 65 (h = 2/(n−1)) | **4.00** (errors 1.2e-4 → 3.0e-8) | 4 |
| `divergence_theorem_rect2d`, midpoint | n = 8 … 64 | 2.0 ± 0.15 | 2 |
| `stokes_theorem_check`, rectangle loop + rectangle surface (2-D) | n = 8 … 64 per side | **2.000** (4.9e-6 → 7.7e-8) | 2 |
| same, tilted 3-D rectangle | n = 8 … 64 | 2.0 ± 0.25 | 2 |
| gradient / curl rotational invariance (V7) | n = 33, 65, 129 / 17, 33, 65 | 2.0 ± 0.25 | 2 |

The implementer's reported orders (gradient 1.994, divergence 1.988, curl 2.052, laplacian 2.026, integrals ≈ 1.998)
are reproduced within the grid-range dependence (ours use n = 16 … 128); the non-periodic field under `bc="periodic"`
was not used (it would falsely show order ≈ −1, as the implementer noted).

## Conservation / invariant residuals
| Quantity | Residual |
|---|---|
| ε–δ relation over 81 cases | 0.0 (exact) |
| Six stress-cube tractions summed (Newton III at a point) | < 1e-14 |
| Net outflux of b × x through boxes / sphere / rectangle / 3 × 3 tiles | < 1e-13 |
| Interior tile faces in `divergence_theorem_tiled` (4 × 4 tiles) | −5.6e-17 |
| ∇×∇φ and ∇·(∇×u) on the grid (one-sided and periodic) | < 1e-11 (round-off: discrete partials commute) |
| ∮∇φ·t ds (2-D polynomial φ; 3-D sin x₁ x₂ + x₃²) | < 1e-13 / < 1e-10 |
| Stokes on discs, b × x, R = 0.2 … 2 | max \|Γ − curl flux\| = 7.1e-15 |
| Invariants I₁, I₂, I₃, trace, both double dots under 100 random rotations | < 1e-12 |
| Length under `transform_vector`; det B of `principal_axes` | < 1e-13; \|det − 1\| < 1e-12 |

## Benchmarks used (value | our value | source + URL | date verified)
| Benchmark | Value | Ours | Source | Verified |
|---|---|---|---|---|
| Divergence theorem example F = (2x, y², z²), unit sphere | 8π/3 = 8.37758… | rhs rel 8.5e-16, lhs rel 4.2e-16 (Gauss–Legendre/midpoint, n = 24); FD-divergence fallback 1e-9 | Wikipedia "Divergence theorem" https://en.wikipedia.org/wiki/Divergence_theorem | 2026-09-16 |
| Levi-Civita identities ε_ijk ε_imn = δ_jm δ_kn − δ_jn δ_km; ε_jmn ε_imn = 2δ_ij; ε_ijk ε_ijk = 6; (a×b)_i = ε_ijk a_j b_k; det A = ε_ijk a_1i a_2j a_3k; pseudotensor | exact | exact / 1e-13 | Wikipedia "Levi-Civita symbol" https://en.wikipedia.org/wiki/Levi-Civita_symbol | 2026-09-16 |
| Cauchy T_j = σ_ij n_i (first index contracted); σ₁,₂ = (σ_x+σ_y)/2 ± √(((σ_x−σ_y)/2)² + τ_xy²); τ_max = ½\|σ₁ − σ₃\|; invariants in principal values | exact | same index placement; Mohr formula 1e-13 on 20 random 2×2 | Wikipedia "Cauchy stress tensor" https://en.wikipedia.org/wiki/Cauchy_stress_tensor | 2026-09-16 |
| R(θ) = [[cos, −sin],[sin, cos]] active; passive = Rᵀ; Rodrigues formula | exact | `rotation_matrix_2d/3d` entries 1e-15; == scipy `Rotation.from_rotvec` 1e-13 (library cross-check) | Wikipedia "Rotation matrix" https://en.wikipedia.org/wiki/Rotation_matrix ; scipy docs | 2026-09-16 |
| Stokes' theorem statement + right-hand rule (no numeric example on the page) | — | our V1 case 2\|b\|πR² | Wikipedia "Stokes' theorem" https://en.wikipedia.org/wiki/Stokes%27_theorem | 2026-09-16 |

## Numbers from the text (book vs ours) — private values redacted to percentages
| Item | Result |
|---|---|
| Example 2.1 (u_r, u_θ formulas, C matrix) at 5 sampled angles | reproduced to 1e-12 relative (exact formulas) |
| Example 2.2 (n, f, \|f\|, angle for a > 0 and a < 0, τ'₁₁, τ'₁₂ at φ = 30°) | every value reproduced to 1e-12; the book's stray "[?]" ignored |
| Example 2.3 (factors 3 and 2; zeros) | exact (linear fields; central stencils exact) |
| Example 2.4 (S, λ, b¹, b², C, 45°, S') with Γ ≡ S₁₂ | exact to 1e-12 |
| Exercise 2.1a dot product | exact |
Tests: `test_polar_components_V6_book_example_2_1`, `test_example_2_2_V6_book_values`, `test_example_2_4_V6_book_values`,
`test_example_2_3_V6_book_factors_and_exercise_2_1a` (skip without `tests/book_values_ch02.json`).

## Figures reproduced with our code (`outputs/ch02/verify/`, local, git-ignored)
| Figure | File | Visual verdict |
|---|---|---|
| Fig. 2.3 / Ex. 2.1 (our numbers u = (1, 2), θ = 30°) | `fig2_3_polar_components.png` | u decomposes onto the rotated (r, θ) axes with u_r = 1.866, u_θ = 1.232 and the two dotted projections closing a rectangle; \|u\| = 2.236 preserved. |
| Ex. 2.2 stresses vs plane angle + Mohr | `example2_2_traction_mohr.png` | σ_n = a sin 2φ and τ_s = a cos 2φ sit exactly on the dotted closed forms, meet at 30° at 0.866/0.5; all (σ_n, τ_s) points lie on the unit Mohr circle centred at 0 with the eigenvalues ±a at its ends. |
| Fig. 2.7 gradient | `fig2_7_gradient.png` | Elliptical level curves of x₁² + x₂²/4, quiver arrows everywhere perpendicular to them and pointing outward (uphill), longer along x₁ where the slope is steeper; ∂φ/∂n = 1.792 at (1, 0.5) for the drawn n. |
| Fig. 2.8 / Ex. 2.4 principal axes | `fig2_8_principal_axes.png` | The unit square carried by S = [[0,Γ],[Γ,0]] stretches along b¹ = (1,1)/√2 and shrinks along b² = (−1,1)/√2 into a rhombus at 45°, area preserved (traceless S). |
| Convergence orders | `convergence_orders.png` | All operator curves are straight lines parallel to the ∝ h² guide (one-sided and periodic), the Laplacian with a larger constant; the four integral-definition curves are parallel h² lines (curl_component with the smallest constant). |
| Fig. 2.9 / Gauss on a tiled box | `fig2_9_gauss_tiled.png` | Streamlines of the 2-D test field cross the box; the four face fluxes sum to ∮n·Q = 5.66482 = Σ tiles; the interior faces cancel to −5.6e-17; ∬∇·Q = 5.66505 differs by the O(h²) midpoint error. |
| Fig. 2.10 / Stokes | `fig2_10_stokes_circulation.png` | For b × x the circulation dots sit on the 2bπR² parabola; the singular vortex's circulation is flat at 2πK for every radius — the hypothesis failure the code flags. |

## Deviations & justifications
No `# DEVIATION` markers in the ch02 code. Book typos handled as decided by the curation and tested: Ex. 2.4 Γ ≡ S₁₂
(`test_example_2_4_V1…`, `strain_rate_tensor(G)` with G = [[0, 2Γ],[0, 0]]), Ex. 2.6 u_y (`test_stokes_V1_curl_component_2_35_uses_u_y…`,
u = (0, z, 0) → −1 whereas the printed u_z would give 0), (2.27) lower limit i = 1 (`test_antisymmetric_vector_V2_derivation`,
round trip), Ex. 2.3 a x₃ e₃ (`test_divergence_V2_example_2_3_symbolic`). The book's order (2.12) → (2.15) is reversed in
D06 as designed; the symbolic test follows the design's route.

## Review round (derivation review → implementer fixes → re-verification, 2026-09-16)
`reports/ch02_review.md` found two physics errors that the first-round tests had **shared with the code** (they asserted
the code's convention instead of the book's). Re-verified from the rendered page p082 (§2.10) and the review's ch03
citations ((3.15) R_ij = −ε_ijk ω_k with ω = ∇×u; (3.17) R_ij = ∂u_i/∂x_j − ∂u_j/∂x_i):

| Item | Was | Now (code, fluidpy only) | Test change | Result |
|---|---|---|---|---|
| **M1** `rotation_tensor(G)` | returned A = ½(G − Gᵀ) and called it R; the first-round test asserted `rotation_tensor(B) == antisymmetric_part(B)` | returns the book's R = G − Gᵀ = 2A (p082: "R is the rotation tensor corresponding to the vorticity vector ω"; for u = b × x, ω = ∇×u = 2b and R_ij = −ε_ijk ω_k = 2[b]ₓ = G − Gᵀ) | `test_decomposition_V1_unique_parts_and_component_counts`: asserts R = 2A = B − Bᵀ, R ≠ A, and the discriminating check vector(`rotation_tensor`(G)) = 2b = ∇×u vs vector(`antisymmetric_part`(G)) = b for u = b × x; R = `antisymmetric_from_vector`(∇×u); G = S + ½R | pass (1e-14) |
| **M2** `velocity_gradient_preset` | pure_strain = [[0,Γ],[Γ,0]], uniaxial_extension = diag(Γ, −Γ), irrotational_strain alias of pure_strain; the first-round test enshrined those matrices (diag(e^{1.4}, e^{−1.4}) for "uniaxial") | design Part C 6.6: pure_strain = diag(Γ, −Γ), uniaxial_extension = diag(Γ, 0) (∇·u = Γ), irrotational_strain = [[0,Γ],[Γ,0]] ≠ pure_strain; simple_shear, solid_body_rotation unchanged | `test_linear_flow_V1_kinematics_presets_and_material_lines`: matrices per design, e^{Gt} = diag(e^{1.4}, e^{−1.4}) for pure strain and diag(e^{1.4}, 1) for uniaxial, trace Γ for uniaxial and 0 for the other four, spin content (R ≠ 0 only for simple shear and solid-body rotation), area preserved (pure strain) vs × e^{Γt} (uniaxial) | pass |
| **O1** `gauss_gradient_box` vector Q | raised in `_eval` (strict xfail) | works; (3, 3) with [i, j] = ∂Q_j/∂x_i | `test_gauss_theorem_V1_vector_Q_gradient_form` (xfail removed): both sides equal the exact matrix; index order = transposed `vector_gradient` | pass (1e-12 / 1e-13) |
| should-fix `divergence` 3-on-2-D | raised | sums the first nd components (∂u₃/∂x₃ = 0), like `curl` / `vector_gradient` | new V1 row in `test_divergence_V1_example_2_3_linear_fields_exact` | pass (exact) |

Lesson recorded for the skill: a test that asserts `f(B) == g(B)` between two *code* functions is not evidence for the
book's convention — the factor 2 between A and R passed unnoticed because both sides were ours. Convention tests must
pin the value to the book's equation on a field with a known answer (here ∇×(b × x) = 2b).

## Open items
- **O1 — closed (review round):** `gauss_gradient_box` now handles a vector Q (`_eval` broadcasts any number of leading
  component axes); `test_gauss_theorem_V1_vector_Q_gradient_form` (formerly the strict xfail) passes with both sides equal
  to the exact ∫∂Q_j/∂x_i dV (1e-12 FD, 1e-13 exact dQ).
- **O2 (implementer docstrings; knowledge-keeper):** `boundary_tangent` (integral_theorems.py l. 527) says "n_c = outward
  in-surface normal to C" and `planar_loop` (l. 566) says "n_c the outward radial direction". With n_c outward, t = n_c × n is
  the *clockwise* tangent about n. The book's Fig. 2.10 (rendered page, printed p. 61) draws n_c up the cap, i.e. **into
  A**, and the code's loops are counterclockwise about n (right-hand rule, verified ∮(x − c)×t ds = 2A n; Stokes lhs = rhs
  positive). The code and the design (N69 ⚠️, E5 parity row with n_c = −r̂) are right; the two docstrings and analysis §9
  item 7 / inventory row #86 ("outward") carry the wrong word. `test_stokes_V1_orientation_right_hand_rule_and_geometry`
  pins the inward convention and shows the outward choice reverses t. `knowledge/notation.md` should record "n_c points
  into A".
- **O3 (designer / derivation-reviewer / notebook-builder):** design Part F D18 ("What it means" and "Traps"), curation
  §4b D18 traps and the analysis claim that A_ij A_ij "is not invariant for a non-symmetric A". It **is** invariant:
  A_ij A_ij = tr(A Aᵀ) and tr(CᵀACCᵀAᵀC) = tr(A Aᵀ) (checked numerically to 1e-12 in
  `test_invariants_V7_unchanged_under_rotation_and_the_chain_remark`). The correct statement is that A_ij A_ij is not the
  λ-coefficient of the characteristic polynomial — I₂ needs the closed chain A_ij A_ji (proved symbolically in
  `test_invariants_V2_derivation`: with A_ij A_ij the cubic's roots are wrong). The derivation's steps and end result are
  right; only the remark should be reworded before it reaches the notebook.
- **O4 (notebook-builder / designer):** the D17 `check_src` cell in design Part F normalises with `b1.norm()`; in sympy
  1.14 `.norm()` introduces `Abs(...)` and `sp.simplify(D[0, 0] - lam1)` does not reduce to 0 (the cell's assert would
  fail). Using `b / sp.sqrt(b.dot(b))` makes every step simplify to 0 (`test_principal_axes_V2_derivation`, which also
  proves CᵀC = I symbolically that way).
- **O5 (orchestrator, informational):** `scripts/ch02_drawings.py` is an imported helper (no `__main__`); run standalone
  it exits 1 with `ModuleNotFoundError: fluidpy` because it does not set `sys.path`. The 14 demo scripts that import it all
  exit 0. Not a defect unless the pipeline expects every `scripts/ch02_*.py` to be runnable.
- **O6 (informational, for the notebook text):** the discrete identities ∇×∇φ = 0 and ∇·(∇×u) = 0 hold to round-off
  (< 1e-11) on the project grid because the difference stencils along different axes commute exactly; the analysis F16
  says "to truncation order". The stronger statement is the one the notebook can make.
- Labels: the scipy comparison of `rotation_matrix_3d` is a library cross-check counted as V5 (`SOURCES.md` says so);
  no `qualitative` or `unverified` label remains in the chapter.

## Verdict: PASS
(review round M1, M2, O1 and the `divergence` should-fix re-verified and passing; remaining Open items O2 docstring wording
and O3/O4 design text/check-cell corrections are routed to the knowledge-keeper, designer and notebook-builder. The single
full-suite failure is the viz-staleness machinery check while the viz-builders are mid-build, outside this phase.)
