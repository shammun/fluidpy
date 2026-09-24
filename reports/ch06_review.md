# ch06 — derivation review (Phase 6)

Reviewer: `derivation-reviewer` (fresh context, read-only), 2026-09-23. Tests at review time: `tests/test_ch06.py` 133 passed;
all 12 `scripts/ch06_*.py` exit 0. No sign, factor, branch or boundary bug found in the physics code.

## Must fix
- **M1 · (6.82) is not a book slip.** The book prints
  $\frac{1}{r}\frac{\partial}{\partial r}(r^2u_r)+\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}(u_\theta\sin\theta)=0$ (6.82),
  which is exactly $r\times$ the Appendix-B divergence
  $\frac{1}{r^2}\partial_r(r^2u_r)+\frac{1}{r\sin\theta}\partial_\theta(u_\theta\sin\theta)$ — correct as printed. The test built a
  hybrid (book's first factor, App. B's second) that the book never prints. Fix: remove (6.82) from every slip list
  (analysis, curation, design, verification report, code docstring); `spherical_continuity_residual` keeps the App. B
  normalisation and says the book's form is $r$ times it; the test asserts the printed form vanishes for the sphere flow
  and equals $r\times$ the residual.

## Should fix
- **S1 · `blasius_state("tilted_ellipse")`**: D, L are x/y components; for α = 15° they are −6.21 / 23.18 N/m while
  |F| = ρUΓ = 24 N/m ⟂ stream (from (6.60) with $c_0 = Ue^{-i\alpha}$: $D - iL = -i\rho U\Gamma e^{-i\alpha}$). Document, return
  complex c₀, `F_perp`, stream angle; `alpha=0` silently becoming 15° must be explicit. E5 text must be mode-aware.
- **S2 · Example 6.2 V4 flux check is a tautology** (Σ diff ψ telescopes to the BCs). Relabel V1 or use a real invariant
  (zero discrete circulation round interior loops / ∮∂ψ/∂n = 0 round a sub-box).
- **S3 · Public tests do not pin Example 6.2 geometry/BCs** (a planted non-uniform outlet passed all public tests). Add a
  public test: 24 unknowns, solid x > 5, y < 2, inlet ψ = Qy/5, outlet ψ = Q(y−2)/3, walls 0 and Q.
- **S4 · Inconsistent ρ defaults** (1.0 / 1.2 / 1000) across force functions — `lift_per_span(U=10, Gamma_cw=2)` = 20 N/m vs
  `surface_pressure_force` 24 N/m. One default per medium or keyword-required ρ.
- **S5 · `superposition_state` status for a net sink** says "no body"; it is an open body facing upstream.
- **S6 · `sphere_potential_vector`** reads a scalar `d_vec` as the radius a — add `a=` and reject scalar `d_vec`.
- **S7 · V5 citations are Wikipedia only** — add primary/textbook sources (Lamb, *Hydrodynamics* §92; Rayleigh 1917, Phil.
  Mag. 34, 94).
- **S8 · Wording**: `harmonic_polynomials` iz² matches (6.25) only up to sign; Exercise 6.42 V6 check should assert the one
  correct cubic form (z/a)²(1 − z/a).

## Verified (against rendered pages)
Γ convention (clockwise in (6.36)–(6.40), (6.52); L = +ρUΓ_cw, Γ_ccw flips); contour forces (6.54)–(6.56) with ccw contour;
Blasius (6.60) residue → −iρUΓ; (6.61) slip is real (−(Ud/π + Γ²/4π²), extra outer square); held source D = −ρmU; 2-D
doublet (6.49) = dipole −d e_x, cylinder d = −2πUa² e_x; 3-D doublet (6.86)–(6.89) and Stokes ψ sign (6.75), flux (6.78);
line sink (6.93)–(6.94) and axial influence matrix; moving sphere (6.96)–(6.105) re-derived, (6.104) middle-bracket slip
real, added mass 2.0944 kg for a = 0.1 m by pressure, energy and volume, bubble a₀ = 2g; log/Joukowski/corner branches;
circle theorem; Example 6.1 (image sign, drift Γ/4πh, wall pressure ρΓ²(η² − h²)/4π²s²); Example 6.2 matrix, sweep order,
geometry, BCs, Fig. 6.25 reproduced to 2 decimals; half-body C_p zero at 113.218°, Rankine oval, Rankine vortex, spheroid K,
Rayleigh collapse, cylinder added mass ρπa², panel self term ½. Scripted validation labels honest except S2.

## Resolution
All items resolved (2026-09-23).
- M1: (6.82) removed from every slip list (analysis §9, curation, design Part G6, code docstrings, test header, verification
  report); new sympy test: printed form vanishes for (6.90) and equals r × App. B.
- S1: `blasius_state` returns F_perp, F_par, stream_angle_deg, c0_im; `alpha=None` convention; test checks
  D − iL = −iρUΓe^{−iα} at four angles. Design G7 makes E5 text mode-aware.
- S2: flux assertion relabelled V1; new V4 test on discrete cell circulation (3.3e-10 converged, falls with sweeps).
- S3: public geometry/BC test for Example 6.2 (Q = 1); planted non-uniform outlet now caught.
- S4: 2-D force helpers default ρ = 1.2; test pins 24 N/m across seven helpers.
- S5: net-sink status "open body (net sink): extends upstream" + test.
- S6: `sphere_potential_vector(..., a=)`; positional scalar kept as deprecated radius (warning pinned).
- S7: Lamb (1932) §92 and Rayleigh (1917) added as primary sources.
- S8: docstring "up to sign"; Exercise 6.42 check asserts the single correct cubic.
Tests after fixes: tests/test_ch06.py 136 passed; full suite 749 passed.
