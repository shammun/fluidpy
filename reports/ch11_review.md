# Chapter 11 review — Instability
2026-10-01 · reviewed state `9f78058` (verify PASS) · reviewer: `derivation-reviewer` (fresh context, read-only) · written by
the orchestrator from its reply.

## Verdict: 2 Must fix (wrong verdict labels in two classifier functions) — no equation is transcribed wrongly

Every numbered equation compared against its page image matches the code symbol for symbol; the twelve claimed printed
slips are real (S1–S4, S6–S8, S11, S12) or correctly described convention changes (S9, S10); S5 cites a ch07 page that was
not reopened. `tests/test_ch11.py -m "not slow"`: 125 passed at review time.

## Must fix

### M1 · `salt_finger_regime` labels viscously stable top-heavy layers "overturning" (`fluidpy/ch11_instability.py:1092`)
- Book: the set (11.45) reduces to the Bénard problem with Ra → Rs − Ra, so the layer is unstable only when
  $\frac{gd^4}{\nu}\left[\frac{\beta}{\kappa_s}\frac{dS}{dz}-\frac{\alpha}{\kappa}\frac{d\bar T}{dz}\right] > 657$ (11.46),
  i.e. Rs − Ra > 27π⁴/4 = 657.5. A top-heavy layer below that threshold is stable.
- Code: `if not sf["density_stable"]: regime = "overturning"` — the computed growth rate is never consulted.
- Measured (d = 5 cm, defaults):

  | dT/dz [K/m] | dS/dz | Rs − Ra | σ_max | label returned |
  |---|---|---|---|---|
  | −0.005 | 0 | 437.8 | −0.159 | "overturning" (wrong) |
  | −0.0074 | 0 | 647.9 | −0.159 | "overturning" (wrong) |
  | −0.0076 | 0 | 665.5 | +0.156 | "overturning" (right) |

- Shows in: the regime map of `scripts/ch11_double_diffusion.py` (strip −0.0075 < dT/dz < 0), explainer E4 status text.
- Missed by the suite because `test_salt_finger_regime_V7_four_regimes_and_eos` only uses a supercritical point.
- Fix: density-unstable → "overturning" only if the margin is positive (Re σ_max > 0), otherwise "stable".

### M2 · `gradient_richardson` / `miles_howard_stable` call a statically unstable, shear-free layer "stable" (`:1490`, `:1502`)
- Book: stability is guaranteed only if $N^2 > \tfrac14 (dU/dz)^2$ everywhere, i.e. $\mathrm{Ri} > \tfrac14$ (11.67) with
  $\mathrm{Ri} \equiv N^2/(dU/dz)^2$ (11.66).
- Code: `Ri = np.where(Up == 0, np.inf, N2v / Up**2)` — +inf wherever U′ = 0, whatever the sign of N².
- Measured: `miles_howard_stable(z, dUdz=0, N2=−1)` → `guaranteed_stable: True, Ri_min: inf`.
- Matters at every shear extremum (jet centre) with N² < 0 — the convective case.
- Fix: sign of N² at U′ = 0 (N² = 0 = U′ does not satisfy the inequality), or test `N2 > 0.25 * Up**2` pointwise.

## Should fix
1. `poiseuille_neutral_curve` docstring: unstable band at Re = 10⁴ is k = 0.797–1.095, not "0.80 to 1.07".
2. `sin_profile_max_growth`: docstring "b = 1.7 → ≈ 0.0076 (k ≈ 0.1)" — the function returns 0.01139 (N = 80), 0.01187
   (N = 120), maximum near k ≈ 0.2; 4 % N-dependence near the 2b → π margin, so "converged" is not earned there.
3. `lorenz_largest_lyapunov` is labelled "converged" with only a 0.7 < λ < 1.2 assertion — relabel `qualitative`;
   `potential_well_demo` similarly over-labelled.
4. `salt_finger_unstable` docstring quotes 61 254 (the g = 9.81 value); the default G0 gives 61 233.19.
5. `benard_growth_rate` with Ra < 0 returns the real part of a complex mode unless `return_complex=True`;
   `benard_free_free_sigma` says "two real growth rates" but returns a complex pair for Ra < 0 — document.
6. KH degenerate inputs: `kh_critical_k(1, 1, 1, 1)` → 0.0 and `kh_unstable_band(0.0, 1, 1)` → (0, inf) for a uniform
   fluid at rest, which is neutral at every k.
7. Taylor numbers for μ ≠ 1 have no cited benchmark row (independent evidence is the μ = 1 Bénard limit + Galerkin);
   ours at μ = −1: 18 662.6 at k = 3.9991 — a cited Chandrasekhar row would make it V5.
8. `analysis/ch11_curation.md:185` writes `rayleigh_number(g, alpha, dT, d, kappa, nu)`; the code is
   `(alpha, dT, d, kappa, nu, g=G0)` — the notebook must call it with keywords.
9. From the verification report (still open): `benard_marginal_Ra_det` error hint (odd mode needs K ≥ 4 at the default
   `Ra_max`; "≥ 1e5" is not enough for K ≲ 1.2), stale "[0.05, 0.5]" comment in `tg_growth`, stale sentence about
   `decay_map_scale` in the `core/stability.py` module docstring, `scripts/ch11_stratified_shear.py` semicircle panel
   not passing `decay_map_scale`.

## Verified (read on page images and compared with the code)
- Kelvin–Helmholtz: (11.13)–(11.18), the criterion, (11.19), (11.20), Ex. 11.1 (a, b), Ex. 11.2, mixing energies.
- Bénard: (11.21), (11.34)–(11.38) and the solver's block matrices, (11.39)–(11.42), the 3 × 3 determinant, (11.43)–(11.44).
- Double diffusion: (11.45), Ra, Rs′, Rs, (11.46); the cubic of `double_diffusive_sigma` re-derived by hand.
- Taylor: ΔE, (11.50)–(11.54), matrix signs, narrow-gap limit.
- Stratified shear: (11.55)–(11.61), companion linearisation, (11.63)–(11.67), (11.68)–(11.72), Ex. 11.10 budget.
- Viscous: (11.77)–(11.80), (11.83)–(11.86), Table 11.1 scales (half-width for Poiseuille, δ* for Blasius), (11.88).
- Chaos: Ex. 11.11, (11.89)–(11.91), r_H = 24.7368.
- Numbers re-run: Ta_c(μ = 0) = 3389.900 at 3.12660; tanh layer max kc_i = 0.1897021 at k = 0.444918; `tg_growth` at
  R ≠ 1 agrees N = 100 vs 200 to 3e-7; docstring examples Ra = 3677.49, KH k_c = 327.0 m⁻¹, ΔU_min = 6.7026 m/s.
- Deviations challenged and accepted: BC elimination, tan and linear maps, `decay_box`, `decay_map_scale`, the Tollmien
  breakpoint, `os_3d_eigs` taking U′.

## Resolution (2026-10-01) — all Must-fix items closed; verification PASS (`reports/ch11_verification.md`)
- **M1 closed.** `salt_finger_regime` labels a top-heavy layer "overturning" only if the (11.46) margin is positive or
  a root of the cubic grows; otherwise "stable". At d = 5 cm, dS/dz = 0: dT/dz = −0.005 → stable (σ_max = −0.1586),
  −0.0076 → overturning (+0.1562); the label flips at dT/dz = −0.007509. The verifier then found the verdict *text*
  wrong where the growing root is real ((−0.03, −0.0016): roots +16.890, +7.322, −142.8); a "monotonic overturning"
  branch was added. Pinned by tests and mutants M26, M27.
- **M2 closed.** Ri at U′ = 0 is sign(N²)·∞ (nan for 0/0) and `miles_howard_stable` tests N² > ¼(dU/dz)² pointwise;
  (U′, N²) = (0, −1) → −∞, not stable. Independent check: for U ≡ 0, N² = −1 between walls the Taylor–Goldstein
  solver returns the exact growing modes. Mutants M28–M30 killed.
- **Should-fix 1–6, 8, 9 done**; 7 partly (co-rotating Taylor branch has a cited asymptotic check, agreement ≤ 1.1e-4;
  no citable value for μ < 0 — still open).
- **Found later in the same phase (by explainer builders and the lesson reviewer) and fixed:**
  1. Rayleigh table false zeros near neutral wavenumbers (Bickley jet k = 1.8, 1.9: 0 → c_i = 0.02432, 0.01163) — new
     complex-path solver `rayleigh_eigs_contour`, agreeing with shooting to 1e-10; `sin_profile_max_growth` now
     converged (b = 1.7 → 0.0118854).
  2. Orr–Sommerfeld far-field box: the Bickley lower neutral branch was a fixed-box artefact (true picture: one unstable
     band to Re ≈ 17.5, then a stable gap between a long-wave band and the main band); Blasius lower branch +9 % at
     Re = 6000; three grids regenerated with the box max(profile box, 12/k).
  3. Blasius critical point now box-converged: Re_c(δ*) = 519.0601 (was 519.0765 in a box of 20 δ*).
  4. `benard_growth_rate` could return the second eigenvalue (rigid–free K = 1.84, Ra = 1.6 Ra(K): −39.61 instead of
     +7.1195); it now returns the leading one or raises.
  5. `rayleigh_eigs_contour(bc="semi_infinite")` returned no mode for growing modes; the option is now refused.
  6. Printed slip S13: p. 503 prints "(7.128)" beside N²; Chapter 7 defines it as (7.127).
- Final state: `tests/test_ch11.py` 160 tests (154 + 6 slow) pass; full suite 1620 passed; 61 mutants planted over all
  loops, 1 survivor (equivalent on every reachable case).
