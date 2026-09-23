# ch03 — Kinematics: independent derivation review

Reviewer: `derivation-reviewer` (fresh context, read-only), 2026-09-23. Scope: `fluidpy/ch03_kinematics.py`,
`fluidpy/core/{kinematics,coords,vortices,transport}.py`, the ch02 re-exports, `scripts/ch03_*.py`, `tests/test_ch03.py`,
`reference/ch03/`, against the rendered book pages. Re-ran `tests/test_ch03.py`: 120 passed. Seven extra wrong-convention
variants of the reviewer's own (cylinder direction in the fluid frame, sign of the u_r term in (3.23), warm/cold advection
labels, (3.6) without |u|, halved Gaussian ω_z, solid-body ω₀ vs 2ω₀, sector-leg sign) were each caught by the suite.

## Must fix
None.

## Should fix
1. `fluidpy/core/kinematics.py` `pathline` — labelled "converged" but the test asserts no order → relabel `analytic` (V1 vs expm and the Ex. 3.1 circle) or measure an order.
2. `fluidpy/core/transport.py` `rtt_ellipse_2d` — "converged" from a one-Δt finite-difference agreement → `analytic`. `fluidpy/core/vortices.py` Wikipedia forms labelled "benchmark" (V5) are form cross-checks, not published numbers → relabel or cite a textbook.
3. `fluidpy/core/kinematics.py` `velocity_gradient_preset` — "simple_shear (Example 2.4's flow; S₁₂ = Γ/2)" conflates Γ = du₁/dx₂ (= ch03's γ = 2S₁₂) with Ex. 2.4's Γ ≡ S₁₂; spell out the split. `knowledge/notation.md` Γ rows need the same split and must note ch03 already uses Γ for circulation.
4. `fluidpy/ch03_kinematics.py` `parallel_shear_kinematics` — for γ < 0 returns diag(|γ|/2, −|γ|/2) and −45°, not the documented diag(γ/2, −γ/2) → document "stretching axis first".
5. `fluidpy/core/kinematics.py` `material_derivative_terms` — default `ht = h` reuses a length step [m] as a time step [s] → own documented default.
6. `fluidpy/core/kinematics.py` `volumetric_strain_rate` — `float(np.trace(G))` fails for (d, d, N) gradients → trace over axes 0, 1 and scalar-if-0d.
7. `scripts/ch03_fig3_18_rtt.py` — boundary loop omits the last segment (gap in the swept band) → wrap the index.
8. `tests/test_ch03.py` — compares with the book's rounded 1.12091 → use the Lambert-W value (rule 9: book numbers only in the private JSON).
9. `reports/ch03_verification.md` open item O4 is out of date (docstrings now say V3).

## Verified
- (3.9) and the Fig. 3.2 frames: x = x′ + Ut + x′₀, u′ = u − U; the cylinder moves at −U e_x in the fluid frame; worked number local_y = −0.592593 recomputed.
- R = G − Gᵀ (no ½), ω_k = −½ε_ijk R_ij, ω₃ = R₂₁ in 2-D, spin = ½ω, material-line rate −γ sin²θ, closing rate 2n₁·S·n₂, γ = 2S₁₂, ω′ = ω − 2Ω, (3.19) sign as printed.
- (3.23) polar vorticity (non-axisymmetric field: 1.0508701 both ways); (3.28)/(3.29) symbol by symbol; the Gaussian maximum from eˣ = 1 + 2x and its Lambert-W branch.
- (3.30)–(3.35) Leibniz and RTT signs; Ex. 3.2 cone geometry (side velocity, normal, area element, Jacobian, b·n = 0 on the base); box, cylinder, sphere and ellipse normals/velocities.
- Ex. 3.1 circle centres and parametric forms; coordinate conventions (θ from +z, e_r × e_θ = e_φ).
- ch02 re-exports byte-identical, no import cycle; the simply-connected caveat for u = ∇φ is flagged and the path dependence reported.

## Resolution
Should-fix 1–7 → `concept-implementer`; 8–9 → `math-verifier`; the notation.md split → `knowledge-keeper` (phase 9).
