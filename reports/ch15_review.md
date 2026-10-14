# Chapter 15 — Compressible Flow: derivation review

**Reviewer:** `derivation-reviewer` (fresh context, read-only), 2026-10-13. Written up by the orchestrator from the reviewer's reply.
**Scope:** `fluidpy/core/gasdynamics.py`, `fluidpy/core/acoustics.py`, `fluidpy/ch15_compressible_flow.py`, `scripts/ch15_*.py`,
`tests/test_ch15.py`, `reference/ch15/`, `reports/ch15_verification.md`.
**Pages rendered and read:** p761, 762, 765, 767, 772, 777–782, 788–791, 795, 798, 800–802.
**Tests:** `tests/test_ch15.py` without the script-runner test — 106 passed, 1 deselected, 50.5 s.

**Verdict: no wrong book equation; 1 Must fix (gas-constant mix in one function), 9 Should fix.**

## Must fix

1. **`GD.mass_flow_vs_back_pressure` (gasdynamics.py ≈1014–1041) mixes two gases.** With `dimensionless=False` the critical
   ratio, exit Mach number and flux function use the `gamma` argument while R comes from `gas`. A caller passing
   `gas=PerfectGas(1.25, 320)` and the default `gamma` gets γ = 1.4 with the rocket gas's R: 20.966 kg/s instead of
   20.149 kg/s (4.1 % high) at p_B/p₀ = 0.3, p₀ = 3 MPa, T₀ = 3000 K, A = 0.01 m². The correct form uses one γ:
   $\dot m = \dfrac{p_0 A}{\sqrt{R T_0}}\sqrt{\gamma}\,M\left[1+\tfrac{\gamma-1}{2}M^2\right]^{-(\gamma+1)/(2(\gamma-1))}$.
   No test reaches it (the only dimensional test uses air). Fix: take γ from `gas` when given; refuse an explicit
   disagreeing `gamma`. Test: a γ = 1.25 gas on the plateau equals `choked_mass_flow(p0, T0, A, gas)`.

## Should fix

1. `ch15.stagnation_pressure` re-exports ch04's incompressible p + ½ρ|u|² under a bare name that, in this chapter, means
   $\dfrac{p_0}{p}=\left[1+\dfrac{\gamma-1}{2}M^2\right]^{\gamma/(\gamma-1)}$ (15.29); 5 % low at M = 0.8. Rename or wrap with a warning docstring.
2. Three air-constant sets share one namespace (`CP_BOOK` = 1004.5, `PerfectGas()` R = 287.058 → C_p = 1004.70,
   `AIR_BOOK` R = 287.0): `stagnation_temperature(300, 400)` = 379.6416 K vs 300·`stagnation_temperature_ratio` = 379.6255 K.
   State it at call sites / pass `cp` explicitly in the notebook (trap T4).
3. `ch15.rocket_nozzle(…, gas=None)` silently defaults to air for a rocket (A* 114.16 → 116.82 cm²).
4. Argument order differs between siblings: `ackeret_pressure(delta, M_inf)`, `ackeret_cp(delta, M_inf)` vs
   `GD.weak_wave_pressure(M, delta)` for the same $\dfrac{p-p_\infty}{p_\infty}=\dfrac{\gamma M_\infty^2\delta}{\sqrt{M_\infty^2-1}}$ (15.57).
   A swap returns NaN (not silent), but it is a learner trap.
5. `GD.normal_shock` (admissible = M₁ ≥ 1) and `GD.shock_admissible` (M₁ > 1) disagree at M₁ = 1. Use one rule (≥ 1).
6. `GD.influence_coefficients` docstring: `printed_15_47` / `printed_15_48` change more keys than the first bullet says;
   add that the book's increments are inlet-referred (df_book = (p/p₁)·df_local, dq_book = (h/h₁)·dq_local), which is
   why $\dfrac{du}{u}=\dfrac{p_1}{p}\dfrac{df}{1-M^2}$ (15.47) carries p₁/p.
7. `ch15.mach_regime` docstring misstates the book's hypersonic class as "M ≫ 1"; the book sets a definite lower
   threshold than the default 5 (value private, T3).
8. Label over-claims (confirming the verifier, plus one): `AC.breaking_time` → V1 · analytic; `GD.isentropic_state` →
   V1, V5; the four table functions → V7 (structure); level text V7 → V1, V7 for `duct_role`,
   `rayleigh_temperature_sign`, `laval_area`, `strong_shock_limits`, `flat_plate_wave_pattern`; `ch15.ackeret_pressure`
   to cite `test_flat_plate_V1_linear_theory_formulas`; cite the six second-pass tests in their functions' docstrings;
   `GD.prandtl_meyer` drop "conserved".
9. `GD.deflection_angle` docstring: "negative for σ < μ₁" holds only on 0 < σ < μ₁; the sign of δ is not a validity
   test — point to `oblique_shock` / `oblique_shock_components`.

## Verified

- **The three corrected slips**, re-derived from the conservation laws and compared with the page images:
  (15.32) third member is −c²dρ/ρ, the only form that gives $d\rho/\rho=-M^2du/u$ (15.33) and
  $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34); (15.47) numerator is 1 + (γ−1)M² (the printed form misses by exactly
  2(γ−1)M²/(1+(γ−1)M²)); (15.48) second relation is (γ−1)M²/(1−M²), the only sign consistent with the book's own text.
  Each `printed_15_*` option reproduces the page exactly.
- All other `influence_coefficients` entries re-derived by hand, including the three that are finite at M = 1.
- (15.45)–(15.46): the quadratic, discriminant 1 − 2A(γ+1), conjugate form, root selection ("+" always subsonic),
  f_max / q_max; 30 independent brentq cases agree to 4e-13.
- Normal shock (15.39)–(15.43): standard values at M₁ = 2; the (15.43) series; station 1 upstream; (15.40) involution guarded.
- Oblique shock (15.50)–(15.53), (15.55); δ_max and sonic locus; weak/strong roots; M₂ cross-checked to 1e-10.
- Prandtl–Meyer (15.56): ν(2) = 26.380°, ν(3) = 49.757°; `turn = −δ_book` is the book's sign; fan geometry and staircase step.
- Ackeret (15.57)–(15.59): surfaces, chord b, lift and drag split; flat-plate wave-pattern angles.
- Becker shock structure: ODE, implicit relation, thickness formula; book slips #7 and #8 confirmed.
- Acoustics (15.12), (15.16), (15.17), linearisation error, breaking time.
- Fanno/Rayleigh "ours" forms against standard values; never mixed with A/A*.
- Exact Riemann solver vs Toro (p* = 0.303130, u* = 0.927453); MacCormack + artificial viscosity conservative.
- Tests are independent where it matters (finite differences, sympy on the finite duct laws, flux tests written from the jump laws).

## Disposition

- Must fix 1 and Should fix 1–9 → `concept-implementer` (fluidpy only, public signatures stable); the γ ≠ 1.4 test → `math-verifier`.
