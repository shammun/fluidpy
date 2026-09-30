# Chapter 9 — derivation review (independent, read-only)

Verdict: no physics bug that changes a headline result. One wrong transcription of (9.8) the tests cannot see, one mislabelled "slip", several misleading defaults and label overclaims. Findings below are from the `derivation-reviewer` reply; status column filled in by the orchestrator after fixes.

## Must fix
1. `fluidpy/core/boundary_layer.py` (~:241, docstring ~:253) — (9.8) coded with the wrong overall sign: page reads (1/Re)(u*v*_x* + v*v*_y*) = −dp*/dy* + (1/Re²) v*_x*x* + (1/Re) v*_y*y*, so dp*/dy* = −(1/Re)(u*v*_x* + v*v*_y*) + (1/Re²)v*_x*x* + (1/Re)v*_y*y*. Code has the exact negative. `delta_p` sign wrong (ratio/order use |·|). Test (`tests/test_ch09.py` ~:182) asserts abs(delta_p) and did not kill the sign mutant. Fix code + docstring, add a sign test against a manufactured field.
2. Slip R5 is not a slip: the book explicitly says "one side of the plate". `book_slips()` R5, module docstrings and analysis R5 misrepresent it. Relabel as a ⚠️ trap/note (the `sides=` argument stays).

## Should fix
- Unrecorded printed slip (R17): for n < 0 the book says (d²u/dy²)_{y=0} < 0 in adverse gradient; (9.9) at the wall gives u_yy = (dp/dx)/μ > 0. Code right; record so notebook text does not copy it.
- R11 stated too strongly: only "no attached bounded (0 ≤ f′ ≤ 1) solution below the fold m = −0.09043" is established.
- R8 dimension: Ψ is m⁵/s³ (code right); fix analysis text (m⁴/s³).
- `thwaites` default closure separates at λ = −0.0681 vs book Example 9.2 criterion −0.090 (x_sep/L 0.1256 vs 0.1583); `LAMBDA_SEP_FS` typed rounded (−0.0681) vs computed −0.068148 — compute from `_fs_fold_state()`; notebook Example 9.2 should show the book criterion first.
- `falkner_skan` default eta_max=10 not converged near the fold (m=−0.09: 0.018909 vs 0.018872 at 16/30); use 16 for m ≤ −0.05 (also `falkner_skan_fields`).
- `separated_pressure_drag` default cp_base=None gives C_D,p ≈ 2.59 at 82°; default to −1.2 or make required.
- `drag_crisis_pair` docstring "factor 3–4" vs model 0.65 — say "weaker".
- `magnus_sign`: docstring wrong about which side moves with the flow (book: larger relative velocity on the side moving against the flow); returns "+" for both-subcritical (unsupported). Sphere Re_cr 3e5 vs book 5e5 (G9) still open.
- Validation labels overclaim: `karman_pohlhausen` ("converged" is 6 % agreement → approximate), `plate_drag_coefficient` ("benchmark" but secondary-sourced 0.074), `blasius_fields`/`falkner_skan_table` ("converged" with no asserted order), `thwaites_l` (ANSYS handout is not a benchmark).
- `displacement_thickness`/`momentum_thickness` tail=True assumes exponential decay; Blasius deficit is Gaussian — document (effect small).

## Verified against rendered page images
Blasius (9.27), Falkner–Skan (9.34)–(9.36) incl. Töpfer scaling, constants (f″(0) = 0.3320573362, η99 = 4.90999, δ* = 1.72079, θ = 2f″(0), H = 2.5911), Hiemenz 1.2325877, m_sep = −0.0904286; erfc far field; Thwaites (9.43)–(9.50) incl. L = 2l − 2(2+H)λ, λ(0) = 0.075, cylinder closed form (103.11° at −0.09), Example 9.2 closed form; free jet (3f‴ + ff″ + f′² = 0, sech², C = 4√6/3, ṁ = (36Jρ²νx)^{1/3}, entrainment, Re_x/Re_h99); wall jet (R3 genuine slip: correct 4f‴ + ff″ + 2f′² = 0; (9.83); f″(0) = f∞³/72; Ψ = C²ν f∞⁴/40; ṁ ∝ x^{1/4}); R6 genuine slip (5.6152 is the 4 % point; 1 % coefficient √6·arccosh10 = 7.3319); C_D,p = sin φ_s(1 − (4/3)sin²φ_s − C_b); street speed (Γ/2a)tanh(πb/a), b/a = arccosh(√2)/π; von Mises marching, σ² mapping, wall-shear extrapolation; sign conventions (dp/dx = −ρU_eU_e′, adverse dp/dx > 0, u = ψ_y, v = −ψ_x).

Note: the reviewer's own mutation run covered 6 of 25 planted mutants (5 killed; the p_y sign mutant survived — Must-fix 1). The verification report's "25 mutants killed" claim is therefore unconfirmed by this review.
