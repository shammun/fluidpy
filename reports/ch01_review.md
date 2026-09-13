# Chapter 1 — Introduction: derivation review

Reviewer: `derivation-reviewer` (fresh context, read-only), 2026-09-13. Written up by the orchestrator from its reply.
Scope: `fluidpy/ch01_introduction.py`, `fluidpy/core/{thermo,dimensional,stratification,statics,diffusion,units,_util}.py`,
`scripts/ch01_*.py`, `tests/test_ch01.py`, `reference/ch01/`, `reports/ch01_verification.md`.
Suite re-run by the reviewer: 101 passed in 7.2 s.

**Verdict: 2 Must fix, no sign or factor error in any transcribed equation.**

## Must fix
1. **Boundary-layer default is unstable, not neutral** (verification F1, upgraded). `synthetic_boundary_layer_profile` /
   `synthetic_boundary_layer_column` default `mixed_dT_dz = -0.0098` K/m is steeper than Γ_a = −g/C_p = −9.7607 K/km, so
   by (1.32) the "near-neutral" mixed layer has N² = −1.354e-6 s⁻² and is labelled `unstable` (book p. 19, Fig. 1.9:
   the mixed layer is drawn at the adiabatic slope and is neutral). Used by the lapse-rate script, the C55 slider and
   the E4 inversion preset. The test hard-codes −0.0098 and checks only above 810 m.
   Fix: default `mixed_dT_dz = adiabatic_lapse_rate()`; test that the stability below 790 m is `neutral`.
2. **Book wording in committed docstrings** (CLAUDE.md rule 9; `tools/check_public.py` does not catch it):
   `fluidpy/core/dimensional.py:225-226` (definition of rank, p. 23), `fluidpy/core/stratification.py:232` (parcel
   forces phrase, p. 18), `fluidpy/core/dimensional.py:17` (kmole phrase, Example 1.2, p. 27). Fix: own words.

## Should fix
- `adiabatic_lapse_rate`: `T` given without `alpha` is silently ignored (perfect-gas −g/C_p returned even for water).
- `blast_energy`: K = 0.856 is for a free spherical blast; the book's Fig. 1.11 is a ground hemisphere (energy 2E
  equivalent) — state the geometry.
- `wall_shear_history`: returns +μ ∂u/∂y at both walls; state the direction (stress on the moving plate is −μ ∂u/∂y).
- `standard_atmosphere(z)` takes geopotential altitude but the argument is called `z`; the notebook plan passes
  geometric heights (80 km).
- `ftcs_diffusion_1d`: unknown `bc` strings silently become Dirichlet; the periodic node convention is undocumented.
- `couette_startup_profile` at t = 0: absolute `np.isclose` tolerance fails for tiny h.
- `scripts/ch01_parcel_oscillation.py:66`: ocean excess gradient 0.1 kg/m⁴ is ~10× a strong thermocline.
- Argument-order traps between siblings: `irreversible_process(kind, T1, v1…)` vs `path_heat_work_totals(kind, v1, T1…)`
  / `perfect_gas_entropy_change(T1, v1, …)`; `molecular_pressure(n, m, velocities)` vs
  `wall_impact_pressure(velocities, m, n)`.
- Test gaps: the nonlinear parcel's division by parcel density is not pinned (a patched environment-density version
  still passes); the EOS-80 UNESCO check values are reproduced (1027.67547 at S = 35, 5 °C; 1023.34306 at S = 35,
  25 °C; 999.96675 at S = 0, 5 °C) and should become a cited V5 test, closing O2/F4.
- Labels: every docstring still says "Validation (planned) … Label: pending" (F2). Wedge (C18) is analytic, not
  converged; the C50 ζ₀ → 0 limit is V7, not V3; "conserved" on C51/C54/C69 is a consistency sweep. Claims not supported:
  `water_viscosity` ±1 % (−2 % at 0 °C), `water_density` NIST V5 (no stored reference), `molecular_pressure` 0.5 %
  (test allows 1 %), `van_der_waals_pressure` sympy check (test is numerical).
- Public repo: `tests/test_ch01.py:789` asserts the book's own worked scale-height example at 250 K — move to the
  private `book_values_ch01.json`.

## Verified (symbol by symbol against the rendered pages)
- (1.5) concave-side pressure and signed radii; Fig. 1.5 wedge balances (p₂ − p₁ = ½ρg dz, p₃ = p₁).
- (1.7)–(1.9), z up; `layered_pressure`, `net_pressure_force_on_box`, `linear_lapse_pressure`, isothermal law with H = RT/g.
- Example 1.1: α = π/2 − θ_c; rise consistent with the Laplace jump through a meniscus of radius R/sin α.
- (1.10)–(1.11) work done on the system (−p dv); (1.17)–(1.18); (1.19); sign of (1.20); (1.21)–(1.28) kmol constants;
  path corner states, q = Δe − w, stirring and free-expansion entropy.
- (1.29), (1.30) signs; exact equal-pressure buoyancy g(T_p − T_e)/T_e; parcel cooling −(g/C_p)T_p/T_e.
- (1.31)–(1.34); (1.35) as "same sign" (dρ_θ/dz = (ρ_θ/ρ)(dρ/dz − dρ_a/dz)).
- (1.39) and its minors; pipe exponents (−2, 0, −1), Π₁ = 32 Π₂Π₄; Examples 1.2–1.5.
- Couette series at t = 0; zero-flux ghost node conserves the trapezoid integral; all end stencils second order.
- Kell, IAPWS, USSA composition, Sutherland coefficients; CO₂ van der Waals per-mass conversion; Tait c = √(K₀/ρ₀);
  EOS-80 coefficients.
- A flipped sign or wrong factor in Γ_a, N², (1.35), θ, ρ_θ, the Laplace jump, capillary rise, the wedge ½, first-law
  signs, the Neumann factor 2, the top-wall stress, ν = μ/ρ or kmol/mol would fail the suite.

## Actions
- Must fix 1–2 and the code/docstring Should-fix items → `concept-implementer` (signatures stable).
- Test changes (boundary-layer neutral check, parcel first integral, EOS-80 V5, private scale-height value, labels)
  → `math-verifier` after the implementer.
