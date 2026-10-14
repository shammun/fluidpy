# Chapter 15 verification — Compressible Flow                       2026-10-13 · phase 5, third pass (after the derivation-review fixes), working tree above `7472cbb` (review fixes uncommitted)

## Verdict: PASS — `tests/test_ch15.py`: 111 passed, 0 failed, 0 skipped (127.9 s, slow script test and the five private V6 tests included); review Must-fix 1 and Should-fix 1, 3, 5 pinned by 4 new tests and 1 changed expectation, each new pin shown to fail on the pre-fix code where there was a pre-fix behaviour; label audit: no over-claim left (citation-only residuals listed); full suite not re-run in this pass (merge gate)

### Third pass (derivation-review fixes, `reports/ch15_review.md`) — summary
Code changed since `7472cbb` (AST without docstrings): `GD._one_gamma` (new), `GD.mass_flow_vs_back_pressure`,
`GD.shock_admissible`, `ch15.rocket_nozzle`, `ch15.stagnation_pressure` / `stagnation_temperature` (now thin wrappers of
ch04's functions) and the alias `ch15.stagnation_pressure_incompressible`; the rest is docstrings.

| Review item | Test | What it pins | Discriminates? |
|---|---|---|---|
| Must-fix 1 (one γ per call) | **new** `test_convergent_nozzle_mass_flow_V1_one_gamma_per_call` (V1) | for `PerfectGas(1.25, 320)`, p₀ = 3 MPa, T₀ = 3000 K, A = 0.01 m², p_B/p₀ = 0.3: choked, ṁ = `choked_mass_flow` (< 1e-14) = the closed form p₀A√(γ/RT₀)[2/(γ+1)]^{(γ+1)/(2(γ−1))} (< 1e-13) = 20.149037 kg/s; 60 subsonic back pressures against a hand-written single-γ ρuA (< 1e-12); explicit `gamma=1.25` identical; `gamma=1.4` with that gas → `ValueError`; default/`gamma=1.4`/`gas=AIR` calls bit-identical, parity value 0.56066 kept | on the pre-fix module the plateau is 20.966 kg/s, a 4.05 % error → fails. Pre-fix vs post-fix on every air or single-γ call (4 γ × 2001 back pressures, dimensionless and dimensional): max difference **0.0** |
| Should-fix 5 (M₁ = 1 admissible) | **changed** `test_shock_entropy_V1_equals_the_entropy_state_function_of_the_jump_state` | expectation `[False, False, True]` → `[False, True, True]` for M₁ = 0.5, 1, 1.5; `shock_admissible`, `normal_shock(...)["admissible"]` and sign(ΔS ≥ 0) agree on 29 Mach numbers × 4 γ including M₁ = 1 exactly and its two neighbouring doubles | fails on the pre-fix `>` rule |
| Should-fix 1 (incompressible p₀) | **new** `test_stagnation_pressure_V1_incompressible_wrapper_is_not_15_29` (V1) | alias is the wrapper (and not in `__all__`); both equal ch04's function and p + ½ρU² on M = 0…0.95; deficit 1 − (p + ½ρU²)/p₀ against (15.29) written out (< 1e-12); **0.050081 at M = 0.8** (the trap pinned: the incompressible value is 5 % low), 0.13 % at M = 0.3; ratio of excesses = `pitot_compressibility_factor` | the alias did not exist before |
| Should-fix 3 (rocket gas) | **new** `test_rocket_nozzle_V7_warns_without_a_gas` (V7) | `gas=None` (positional default or keyword) → `UserWarning`; silent with a gas; fallback = `PerfectGas()`; fallback throat area differs from the exhaust gas's | no warning before |
| Should-fix 8 (labels) | **new** `test_small_helpers_V1_against_independent_closed_forms` (V1) | `strong_shock_limits` against the five closed forms for 4 γ (1e-15); `laval_area` against the two parabolas on 201 points (1e-15); `duct_role` against the sign of d(A/A*)/dM (39 M × 4 γ × 2 signs); `rayleigh_temperature_sign` against the Rayleigh-line T/T* and T₀/T₀* (not (15.48)) | gives the V1 the docstrings now claim a test of its own |

Physics ruling on Should-fix 5: (15.43) gives S₂ − S₁ = 0 at M₁ = 1 exactly (the jump is the identity, a sound wave);
the second law demands ΔS ≥ 0, not > 0, so M₁ = 1 is admissible. Agreed; the expectation was updated.

### Second pass (fix loop 1) — summary
The implementer changed code in six functions of `core/gasdynamics.py` (and docstrings in all three modules). Both
first-pass failures are fixed: `test_influence_coefficients_V7_sonic_point_is_infinite_as_documented` and
`test_oblique_shock_V7_shock_angle_below_the_mach_angle_is_refused` pass unchanged. Six tests were added (V2 1, V1 1,
V7 4) and four existing ones tightened (no tolerance loosened): the 3 ulp/|M₁² − 1| allowance of the entropy test is
removed; the M₁ = 1 return trip of `test_prandtl_meyer_turn_V4_isentropic_and_composes` is no longer skipped; the
Ackeret quadrature order is now asserted two-sided (|p − 2| < 0.15 instead of p > 1.85); `test_area_ratio_V2_derivation`
now ties the code to the sympy expression. Rulings 6–10 below; label audit and second-pass mutation check further down.
The first-pass text that follows is kept as the record of that pass (its numbers describe the code before the fix).

### First pass
`tests/test_ch15.py`: **99 passed, 2 failed, 0 skipped in 61.4 s** (101 tests; the five private V6 tests ran — the book
file is present on this machine). Full suite `.venv/Scripts/python.exe -m pytest -q tests`: **2078 passed, 2 failed,
0 skipped, 1 warning in 2777.5 s (46 min 17 s)** — the two failures are the two ch15 tests below; chapters 1–14 and the
machinery tests are unchanged (1979 passed). Mutation check: **103 of 103** planted variants killed, and the 24 that were
first caught by a contract test were killed again by a physics test with the five contract tests deselected.

**No wrong equation was found**: every formula of `core/gasdynamics.py`, `core/acoustics.py` and
`ch15_compressible_flow.py` agrees with an independent check (conservation laws written in the test, sympy re-derivation,
NACA Report 1135 read first-hand, quadrature, a second numerical route). The verdict is FAIL because two tests fail, and
both are behaviour at the edge of a function's range — one contradicts the function's own docstring, the other the
module's stated convention and the guard that `normal_shock` already has for the same trap:

| # | Function | Call | Expected | Got | Diagnosis | Suggested fix (implementer) |
|---|---|---|---|---|---|---|
| F1 | `GD.influence_coefficients` | `influence_coefficients(1.0, 1.4, dA_over_A=0.1)["du_over_u"]` (also `df=0.1` alone, `dq=-0.1` alone, `dA_over_A=0.1, df=0.03`) | ±inf — docstring: "±inf at M = 1 unless the driving combination −dA/A + df + dq vanishes (then NaN)"; `GD.area_velocity(1.0, 0.1)` returns −inf for the same equation $\dfrac{du}{u}=-\dfrac{1}{1-M^2}\dfrac{dA}{A}$ (15.34) | `nan` for every key (only a call with all three increments non-zero gives inf) | the three parts are multiplied by `inv = 1/(1 − M²) = inf` separately, so each zero increment gives `0·inf = nan`, which poisons the sum | form the driver first, `du = (−a + f_ + q_) * inv` (the printed variants handled before the multiplication), or mask the zero increments |
| F2 | `GD.oblique_shock` | `oblique_shock(2.0, sigma=0.75·μ₁)` (σ = 0.3927 rad < μ₁ = 0.5236 rad) | `ValueError` or NaN with `attached=False` — module docstring: "outside its range a function returns NaN or raises `ValueError`"; argument documented as "sigma … in [μ₁, π/2]"; `normal_shock` refuses the same thing (trap T6) | a silent "expansion shock": `p_ratio` = 0.5168, `dS_over_cv` = −0.01162, `attached=True`, no flag | with `sigma=` the normal Mach number M₁ sin σ < 1 goes straight into the bare normal-shock formulas | raise for M₁ sin σ < 1 (to round-off), or accept `allow_expansion=` as `normal_shock` does and return `attached=False` |

Tests: `test_influence_coefficients_V7_sonic_point_is_infinite_as_documented`,
`test_oblique_shock_V7_shock_angle_below_the_mach_angle_is_refused`. Neither touches a number the notebook or an
explainer needs away from M = 1 or below the Mach angle; both are one-line fixes. Nothing was loosened to pass.

## Environment
python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · mpmath 1.3.0 · matplotlib 3.11.2 · pint 0.25.3 · pytest 9.1.1
(Windows 11; `.venv`). Run: `.venv/Scripts/python.exe -m pytest tests/test_ch15.py -q -p no:cacheprovider`.

## What the orchestrator asked about — rulings
1. **The three corrected slips** — confirmed three ways each. The corrected defaults agree with (a) finite differences of
   $\dfrac{A}{A^*}=\dfrac1M\left[\dfrac{2}{\gamma+1}\left(1+\dfrac{\gamma-1}{2}M^2\right)\right]^{(1/2)(\gamma+1)/(\gamma-1)}$ (15.31) (area part, to 1e-6) and of
   $M_2^2=\dfrac{-(1-2A\gamma)\pm\left[1-2A(\gamma+1)\right]^{1/2}}{(\gamma-1)-2A\gamma^2}$ (15.46) (friction and heat parts, central differences with Richardson, to 1e-6, eight
   quantities each, six Mach numbers, two γ); (b) a sympy solution of the four linearised laws (all nine returned keys, to
   1e-12); (c) sympy differentials of the finite laws (15.45) themselves. Each printed form **fails the same checks**:
   `printed_15_32=True` misses du/u by 18 % and more (and would make a widening duct slow a supersonic stream);
   `printed_15_47=True` misses dp/p by exactly 2(γ − 1)M²/(1 + (γ − 1)M²) — 7 % at M = 0.3, 18 % at 0.5, 49 % at 0.9, and
   the wrong sign at M = 2; `printed_15_48=True` has the opposite sign of u du/h at every Mach number.
2. **`duct_outlet_mach` exactly at f_max / q_max** — honest and documented, not a defect. At the limit the discriminant
   D = 1 − 2A(γ + 1) is zero up to round-off (measured |D| ≤ 2 ulp) and the outlet Mach number is a square-root branch point:
   |M₂ − 1| = ½(γ + 1)√|D|. Measured worst case over 11 inlet Mach numbers × 3 γ × both limits × both roots:
   **|M₂ − 1| = 1.26e-8**, against the bound ½(γ + 1)√(4 ulp) = 3.6e-8 that the test asserts (a conditioning bound, not a
   relaxed tolerance); the backward error is at round-off (the returned M₂ satisfies A(M₂) = A to 1e-14). Just beyond
   (f_max(1 + 1e-9)) both roots are NaN; just inside there are two roots, one on each side of 1. Two remarks for the
   implementer (not failures): (i) the clamp is one-sided — D ∈ (−1e-13, 0) is set to 0 (so M₂ = 1 exactly), but a
   positive round-off of the same size is not, which is where the 1e-8 comes from; a symmetric clamp would return 1 on
   both sides; a side effect of the clamp is that f up to about f_max(1 + 1e-13) still returns ≈ 1 instead of NaN (the
   docstring says so); (ii) the docstring of `duct_choking_limits` says the outlet "returns 1 there from both sides" — it
   should say "1 to √(round-off) ≈ 1e-8".
3. **Removable 0/0 of (15.46)** — handled. Where 2Aγ² = γ − 1 the numerator and denominator of the "+" root vanish
   together (|num|, |den| < 1e-14 at the two points tested, f = 0 and f = 0.1); the code returns a finite root that
   satisfies the definition of A to 1e-12 and joins its neighbours smoothly (second difference < 1e-9). With f = q = 0 the
   point is M₁ = √((γ − 1)/2γ), where the shock partner of (15.40) runs off to infinity. The planted variant "printed form
   only" is killed by this test.
4. **Flat plate, linear against exact** — the linear theory is the small-angle limit of the shock–expansion result, but
   the miss of C_L and C_D is **second order in α** (measured 2.02, 2.00, 1.99 for C_L at M∞ = 1.6, 2, 3), not first order
   as the docstring of `flat_plate_supersonic` and analysis T29 say: the second-order (Busemann) pressure terms are equal
   on the two surfaces and cancel in p₃ − p₂. A single surface is first order (measured 1.03, 1.02, 1.02) with the
   coefficient of NACA Report 1135's expansion series.
5. **The stray `outputs/ch14/_resume_scratch/mut/tests/conftest.py`** — confirmed: bare `pytest -q` from the repo root
   collects 2080 tests and then stops with `ImportPathMismatchError: ('tests.conftest', …/tests/conftest.py,
   …/outputs/ch14/_resume_scratch/mut/tests/conftest.py)`. The whole folder `outputs/ch14/_resume_scratch/mut/` (9 MB) is
   a copy of the package made by the ch14 mutation check; it is git-ignored (`outputs/`), nothing refers to it, and it is
   **safe to remove** (it also holds a copy of the private ch14 book file, another reason to delete it). I did not remove
   it. The ch15 mutation trees were created under `outputs/ch15/_verify_scratch/` and deleted by the script when it
   finished; what remains there is the mutation script, its two result files, a measuring script and the logs.

## Second pass — review of fix loop 1 (rulings 6–10)
Code changed in six functions of `core/gasdynamics.py` only (`influence_coefficients`, `normal_shock_entropy`,
`oblique_shock`, `oblique_shock_components`, `oblique_shock_curve`, `prandtl_meyer_turn`, plus the two constants
`_SIGMA_RTOL`, `_NU_TOL`) — established by comparing the syntax trees of the three modules with their docstrings removed
against `HEAD`; `core/acoustics.py` and `ch15_compressible_flow.py` changed in docstrings only.

6. **`influence_coefficients` at M = 1 — correct, and as documented.** The six keys that carry 1/(1 − M²) (`drho_over_rho`,
   `du_over_u`, `dp_over_p`, `dT_over_T`, `dM_over_M`, `u_du_over_h`) are ±inf with the sign of the subsonic side whenever
   the driver −dA/A + df + dq is non-zero — one increment alone included — and NaN when it is exactly zero. The three keys
   now returned finite are right: re-derived with sympy from the finite laws of a short duct (mass $\rho_1u_1A_1=\rho_2u_2A_2$ (15.18),
   momentum with the wall's pressure force and a friction force (15.24), energy $h_2+\tfrac12u_2^2=h_1+\tfrac12u_1^2+Q$ (15.22))
   and the *definitions* of entropy, T₀ and p₀ — not from the code's linear system:
   $ds/C_p = dq + \frac{\gamma-1}{\gamma}df$, $dT_0/T_0 = dq/(1+\frac{\gamma-1}{2}M^2)$,
   $dp_0/p_0 = -df - \frac{\gamma M^2}{2}\,dq/(1+\frac{\gamma-1}{2}M^2)$; none contains 1/(1 − M²), the area change drops out of all
   three, and they satisfy $ds/C_p = dT_0/T_0 - \frac{\gamma-1}{\gamma}dp_0/p_0$. Away from M = 1 nothing changed: 138 240 values
   (10 Mach numbers × 3 γ × 64 increment triples × the 8 combinations of the printed options × 9 keys) agree with the
   first-pass code to 1.1e-15 of max(1, |value|); the three printed options fail the same checks as before (the three
   first-pass tests are unchanged and pass) and keep their own behaviour at M = 1 (pinned). One remark, not a failure:
   "the driver vanishes" means *exactly* in doubles — `dA_over_A=0.1, df=0.04, dq=0.06` sums to −6.9e-18 and returns
   infinities whose sign is that of the round-off; a caller who wants the 0/0 must pass increments that cancel exactly.
7. **`oblique_shock` with a given σ — correct.** Scalars raise `ValueError("not a shock: …")` for M₁ sin σ < 1 − 1e-12 or
   σ > π/2 + 1e-12 (σ ≤ 0 and M₁ ≤ 1 with σ < π/2 included); arrays return NaN in all ten numerical keys and
   `attached=False` exactly at those entries; σ = μ₁ within the allowance is the Mach wave (M_n1 = 1, δ = 0, p₂/p₁ = 1,
   ΔS = 0 exactly), which matters because `asin(1/M₁)` in doubles gives M₁ sin σ one ulp below 1 for many M₁ (more than 20 of
   the 800 tested). `oblique_shock_components` returns NaN in every key outside the range; `oblique_shock_curve` raises for
   M₁ < 1 (and NaN). In range the σ route equals the first-pass code to 4e-15. NaN inputs propagate quietly
   (`attached=False`), which is consistent with the module convention.
   **`deflection_angle` left unguarded — acceptable.** It is the bare formula $\tan\delta = 2\cot\sigma\,\dfrac{M_1^2\sin^2\sigma-1}{M_1^2(\gamma+\cos2\sigma)+2}$ (15.53),
   the root finder brackets it at σ = μ₁ and σ_max, a first-pass test already relies on its negative value below the Mach
   angle, it has the same status as `normal_shock_mach` beside the guarded `normal_shock`, and the docstring now names the
   two guarded callers. One imprecision in that docstring (Open item 3): the negative sign marks the forbidden root only for
   0 < σ < μ₁ (and π/2 < σ < π − μ₁); the formula is odd in σ and antisymmetric about π/2, so for π − μ₁ < σ < π (and
   −μ₁ < σ < 0) it returns a *positive* angle although M₁ sin σ < 1. The sign of δ is therefore not a validity test; the
   guarded callers are (pinned by `test_deflection_angle_V7_bare_formula_outside_the_range_and_its_guarded_callers`). Every
   call in the chapter module and the scripts is inside [μ₁, π/2].
8. **`normal_shock_entropy` — the round-off edge is gone.** With m = (M₁ − 1)(M₁ + 1) the relative error against 50-digit
   arithmetic is ≤ 3e-15 for |M₁ − 1| from 1e-12 to 9e-3 (first-pass code: 1.5e-9 at |M₁ − 1| = 1e-9, 1.2e-10 at 1e-7) and
   at most 6e-12 on the logarithm branch just above the switch |m| = 0.02. ΔS = 0 exactly at M₁ = 1, with the right sign
   one ulp to either side.
9. **`prandtl_meyer_turn` — the round-off edge is gone.** ν(M₁) + turn in (−1e-12, 0) rad is set to 0 (M₂ = 1 exactly, sonic
   ratios); an overshoot of 2e-12 rad and more is still NaN; a small positive ν is not touched; the ν_max end is
   unchanged. Out from M = 1 and back by the same angle: 0 of 1600 return trips NaN (first-pass code: 683).
10. **Docstring labels — audited, see "Label audit".** No "planned"/"pending" remains; the two untrue statements of the
    first pass (`flat_plate_supersonic`, `duct_choking_limits`) and the blast-constant sentence are corrected.

## Validation table
Levels: V1 analytic · V2 symbolic · V3 convergence · V4 conservation/invariant · V5 cited benchmark · V6 book value
(private) · V7 limits/guards. "Label" is the label each function has earned.

| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers | Label | Notes |
|---|---|---|---|---|---|---|
| C01 wave equation (15.12) | CORE | `AC.convected_wave_residual`, `wave_equation_residual`, `ch15.linearised_acoustics_sympy` | V2 (D01, D02 in three dimensions, independent of the engine; plane wave in any direction) · V3 (stencil residual of a solution) · V7 (a reversed stream is not a solution; U → 0) | residual order 1.96; exact (round-off) for U = 0 | symbolic, converged | |
| C02 d'Alembert (15.15)–(15.17) | CORE | `AC.dalembert`, `plane_wave_velocity`, `velocity_from_pressure`, `dalembert_fg`, `dalembert_initial_value`, `dalembert_in_stream`, `wave_speeds`, `xt_characteristics` | V2 (D03) · V1 ((15.16) = quadrature of (15.14); initial data reproduced; Galilean shift) · V4 (acoustic energy) | quadrature 1e-7 of the scale; energy ptp/E = 4.0e-16 (planted (f + g)/ρ₀c: > 5 %) | symbolic, analytic, conserved | |
| C03 duct laws (15.18)–(15.26) | CORE | `GD.mass_flow`, `total_enthalpy`, `energy_residual`, `wall_force`, `nozzle_isentropic`, `ch15.duct_equations_sympy` | V2 (D06, D07) · V4 (ρuA, h₀, entropy along four isentropic solutions; wall force = ∫p dA by quadrature) | ρuA ptp < 1e-9; h₀ < 1e-12; force 1e-8 | symbolic, conserved | |
| C04 stagnation ratios (15.27)–(15.30) | CORE | `GD.stagnation_*_ratio`, `mach_from_pressure_ratio`, `isentropic_state`, `isentropic_table`, `pitot_compressibility_factor`, `entropy_change_from_stagnation` | V1 (from the ch04 energy law, the isentrope and the entropy state function) · V5 (NACA 1135 Tables I, II) · V2 (series) · V7 (γ → 1, M → 0) | 1e-13; table: 62 of 63 entries to half a unit, all to one unit | analytic, benchmark | |
| C05 area ratio (15.31) | CORE | `GD.area_ratio`, `mach_from_area_ratio`, `mass_flux_function`, `choked_mass_flow`, `sonic_ratios` | V1 (ρ*c*/(ρu) from the three ratios; two inverses) · V2 (D10, D11) · V5 (NACA A/A*; cited 0.528 and gas table) · V7 (√ conditioning at the throat) | 1e-13; round trips 1e-9 | analytic, symbolic, benchmark | |
| C06 area–velocity (15.32)–(15.34) | CORE | `GD.influence_coefficients`, `area_velocity`, `duct_role` | V1 (finite differences of (15.31)) · V2 (sympy linear system; second pass: ds/C_p, dT₀/T₀, dp₀/p₀ re-derived from the finite laws of a short duct with area change, friction and heat — no 1/(1 − M²) in them) · V7 (sign pattern, M → 0; **the sonic point, key by key — F1 fixed**: six keys ±inf with the sign of the subsonic side when the driver ≠ 0, also for one increment alone; NaN when the driver is exactly 0; three keys finite; arrays through M = 1; the printed forms' own behaviour; no numpy warning leaks) | 1e-6 against differences; 1e-12 against sympy; the three regular keys 1e-14 at M = 1 | analytic, symbolic | printed (15.32) fails; values away from M = 1 unchanged by the rewrite (138 240 values, all option combinations, ≤ 1.1e-15) |
| C07 normal shock (15.35)–(15.42) | CORE | `GD.normal_shock*`, `normal_shock_residuals`, `normal_shock_fluxes`, `normal_shock_table`, `hugoniot_pressure_ratio`, `strong_shock_limits`, `ch15.normal_shock_sympy` | V4 (three jump laws written in the test) · V2 (D13 line by line, D14, D15) · V5 (NACA 1135 rows and forms (105), (106), (111), (127); secondary table) · V7 (limits, involution, guard) | flux residual 3.3e-16; table: 36 of 36 to half a unit | symbolic, conserved, benchmark | |
| C08 entropy jump (15.43)–(15.44) | CORE | `GD.normal_shock_entropy`, `weak_shock_entropy`, `shock_admissible`, `normal_shock_stagnation_pressure_ratio`, `ch15.weak_shock_series` | V1 (entropy state function; 40- and 50-digit arithmetic for |M₁ − 1| from 0.3 down to 1e-12; cubic law) · V2 (D17 line by line, D18; second pass: the code against the cubic + quartic series sympy derives from (15.43), remainder of fifth order) · V5 (NACA (146), (147) and its series in ξ) | slope 2.997; relative error ≤ 1e-13 for |M₁² − 1| < 0.02 (measured ≤ 3e-15), ≤ 1e-10 elsewhere (measured 6e-12 just above the switch) — the allowance 3 ulp/|M₁² − 1| of the first pass is removed | analytic, symbolic, benchmark | m = (M₁ − 1)(M₁ + 1) |
| C09 convergent nozzle | CORE | `GD.mass_flow_vs_back_pressure`, `convergent_nozzle`, `critical_pressure_ratio` | V1 (piecewise; value and slope continuous) · V4 (ṁ = ρuA at the exit) · V2 (D19) · V5 (cited ratio) | leaves the plateau at order 2.00 | analytic, conserved, symbolic | |
| C10 Laval nozzle | CORE | `GD.laval_critical_pressures`, `laval_regime`, `laval_shock_location`, `laval_exit_state`, `laval_nozzle`, `ch15.nozzle_with_shock` | V1 (thresholds, 12 regime probes × 15 cases) · V2 + V4 (D20: shock found independently by mass conservation) · V4 (mass flow at every station, across the shock) · V7 (ends of the shock range, continuity) | M₁ 1e-8; mass-flow ptp 2.1e-15 | analytic, conserved | |
| C11 duct, friction and heat (15.45)–(15.46) | CORE | `GD.duct_parameter_A`, `duct_outlet_mach`, `duct_state_ratios`, `duct_entropy_change`, `duct_residuals`, `duct_curve` | V4 (three laws on both roots, 91 roots of 50 cases) · V2 (D21 line by line) · V1 (f = q = 0: inlet or shock) · V7 (removable 0/0) | residuals < 1e-11 | symbolic, conserved | planted "one sign for all M₁" fails |
| C12 friction and heating theorems, choking | CORE | `GD.influence_coefficients` (df, dq), `duct_choking_limits`, `duct_choking_locus`, `rayleigh_temperature_sign`, `fanno_ratios`, `rayleigh_ratios`, `fanno_length_to_f`, `fanno_line_Ts`, `rayleigh_line_Ts` | V1 (finite differences of (15.46); monotone towards M = 1; limits) · V2 (D22, D23, D24) · V5 (cited forms, 0.845) | 1e-6; |M₂ − 1| ≤ 1.26e-8 at the limit | analytic, symbolic, benchmark (forms) | printed (15.47), (15.48) fail |
| C13 Mach angle (15.49) | CORE | `GD.mach_angle`, `mach_from_angle`, `AC.source_wavefronts`, `in_zone_of_action`, `mach_cone_delay`, `doppler_factor` | V1 (tangent from the source to 200 fronts; zone of action against 3000 fronts) · V5 (NACA μ) | 1e-12; μ to 0.005° | analytic, benchmark | |
| C14 oblique shock (15.50)–(15.55) | CORE | `GD.oblique_shock`, `oblique_shock_components`, `deflection_angle`, `max_deflection`, `sonic_deflection`, `detachment_mach`, `oblique_shock_curve`, `weak_wave_pressure`, `weak_wave_error`, `ch15.wedge_flow` | V4 (four laws with velocity vectors, 135 shocks) · V2 (D26, D27, D28) · V5 (NACA (124), (132), (137), (138), (139b), (142), (145), (167), (168), (169), (170), 45.6°) · V1 (two roots, δ_max against a minimiser) · V7 (**the guard — F2 fixed**: a given σ with M₁ sin σ < 1 − 1e-12 or σ > π/2 + 1e-12 raises "not a shock" for scalars and gives NaN in all ten numerical keys with `attached=False` for arrays; σ = μ₁ to round-off is the Mach wave with no jump at all; `oblique_shock_components` NaN in every key outside; `oblique_shock_curve(M₁ < 1)` raises; over 4101 angles × 3 Mach numbers × 2 γ no returned entry has p₂/p₁ < 1 or S₂ < S₁; the bare `deflection_angle` pinned outside its range) | residuals < 2e-11; δ_max 1e-10; in range the σ route equals the first-pass code to 4e-15 | symbolic, conserved, benchmark (forms) | ruling 7 |
| C15 Prandtl–Meyer (15.56) | CORE | `GD.prandtl_meyer`, `prandtl_meyer_max`, `mach_from_prandtl_meyer`, `prandtl_meyer_turn`, `ch15.expansion_fan`, `fan_staircase`, `curved_wall_mach_lines` | V1 (quadrature in two variables) · V2 (D29, D30 line by line) · V5 (NACA ν) · V4 (isentropic, composes; mass flux between rays) · V3 (staircase) · V7 (second pass: the return to M = 1 — 1200 out-and-back turns, none NaN (683 of 1600 were NaN in the first-pass code); ν + turn in (−1e-12, 0) is the sonic state exactly, beyond it still NaN, a small positive ν untouched) | 1e-11; 5 of 6 rows to half a unit, all to one unit | analytic, symbolic, benchmark | |
| C16 flat plate (15.57)–(15.59) | CORE | `ch15.ackeret_pressure`, `ackeret_cp`, `ackeret_flat_plate`, `flat_plate_supersonic`, `ackeret_airfoil`, `flat_plate_wave_pattern`, `lift_slope_vs_mach` | V1 (formulas; double wedge 4τ²/β, half-sine π²τ²/2β) · V7 (limit of the exact result, order measured) · V5 (form: first term of NACA's expansion series) | limit order 2.00 | analytic | see ruling 4 |
| N: the gas, regimes | NOTE | `GD.PerfectGas`, `AIR_BOOK`, `mach_from_speed`, `ch15.mach_regime` | V1, V2 (units) | 1e-13 | analytic | |
| N: small-amplitude family | NOTE | `AC.isentropic_pressure_perturbation`, `acoustic_density`, `acoustic_temperature`, `condensation`, `linearisation_error`, `newton_sound_speed`, `acoustic_impedance`, `sound_pressure_level` | V1 | error order 2.00 with the coefficient ½γ(γ − 1)p₀/ρ₀² | analytic | |
| N: steepening (ours) | NOTE | `AC.simple_wave_speed`, `simple_wave_evolve`, `breaking_time` | V4 (Riemann invariant) · V1 (first crossing; γ = 2 equals ch07) | 1e-10; t_b 2e-4 | analytic, conserved | |
| N: moving shock, blast, thickness | NOTE | `GD.moving_shock`, `ch15.blast_constant_crude`, `shock_thickness_estimate` | V4 (laboratory-frame balances) · V1 | 1e-12 | conserved; blast: **qualitative** | Open item 4 |
| N: shock structure | NOTE | `ch15.shock_structure_sympy`, `shock_structure_constants`, `becker_shock_*` | V2 · V1 (implicit relation, ODE) · V3 (`solve_ivp`) · V4 | 1e-8 | symbolic, analytic | |
| N: rocket nozzle | NOTE | `ch15.rocket_nozzle`, `nozzle_thrust`, `GD.max_speed` | V1 (thrust by two routes; ṁ = choked) | 1e-12 | analytic | |
| N: Riemann solver, Euler run (ours) | NOTE | `GD.riemann_exact`, `ch15.euler1d_run`, `AC.wave_equation_step` | V4 (jump laws, global balance, round-off conservation) · V3 | see below | conserved, converged | DEVIATION |
| N: tables, sketches | NOTE | `ch15.illustrative_inputs`, `conventions_table`, `book_slips`, `traps`; `scripts/ch15_drawings.py` (12) | V7 | — | analytic | |

Counts per validation label (primary level in the test name): **V1 37 · V2 29 · V3 3 · V4 11 · V5 7 · V6 5 · V7 19** = 111
tests (second pass: 34 · 29 · 3 · 11 · 7 · 5 · 18 = 107; first pass: 33 · 28 · 3 · 11 · 7 · 5 · 14 = 101); counting every
level named in a test's header line: V1 39 · V2 31 · V3 6 · V4 16 · V5 8 · V6 5 · V7 27. Every computable CORE row has at least two independent levels, one of them V1, V2, V3
or V5. "Label" is now also what the docstrings say (label audit below).

## Functions used by the notebook and explainers (design Part C)
`test_contract_V7_every_public_name_is_exercised_here` enforces that every one of the 163 public names (82 of
`core.gasdynamics`, 30 of `core.acoustics`, 35 own callables of `ch15`, 16 reused) is called by at least one test in the
file; `…_reexports_are_the_same_objects`, `…_signatures_carry_their_conventions` (keyword-only `branch=`, `method=`,
`printed_*`), `…_scalars_in_floats_out`; `test_contract_V7_design_parity_expressions_evaluate` extracts **every `py:`
expression of design Part B** (70 distinct), evaluates it and checks 31 numbers the design prints beside them (6e-5).
`test_drawings_V7_every_sketch_returns_a_figure` calls the 12 sketches; `test_scripts_V7_every_chapter_script_runs_headless`
runs the 8 scripts (`--no-show --fast`) and `reference/ch15/make_refs.py` (idempotent).

## Derivations (curation §4b, design Part F) — sympy, independent of the chapter's engines
23 `test_*_V2_derivation` tests cover every two- and three-star row: D01, D02, D03 (+ D04), D06, D07, D10, D11, **D13**,
D14, D15, **D17**, D18, D19, D20, **D21**, D22, D23, D24, D26, D27, D28, D29, **D30**. The four three-star rows are
checked line by line (D13: the factorised quadratic (r − 1)[(γ + 1)r − (2γa + γ − 1)], the second root, 1 − r, the return
to momentum; D17: the two brackets, the vanishing first and second orders, the cubic and the quartic coefficients, the
conversion to the pressure jump; D21: the three Mach-number forms, the two sides of A, the quadratic and its
discriminant, f = q = 0; D30: the substitution, the bracket, the partial fractions, the prefactor 1, the two
arctangents). All lines follow. One number in a "Check" is wrong (below, "Found wrong in design / analysis").

## Convergence studies
| Scheme / limit | Grids | Observed order | Design order |
|---|---|---|---|
| central-difference residual of (15.12) on a solution of (15.17) (`convected_wave_residual`) | h = 0.08 … 0.01 | 1.96 | 2 |
| leapfrog (15.13), `wave_equation_step`, C = 0.5 / 0.8 | n = 100 … 800 | 2.002 / 2.001 (pairs 2.006, 2.001, 1.999); exact at C = 1 (1e-11 after 150 steps); raises for C > 1 | 2 |
| `euler1d_run`, smooth simple wave to 0.5 t_b, no smoothing, against `simple_wave_evolve` | n = 100 … 800 | 1.98 (pairs 1.85, 1.97, 1.995) | 2 |
| `euler1d_run`, Sod's problem at t = 0.2, L¹ against `riemann_exact` | n = 100 … 800 | p 1.01 · u 0.99 · ρ 0.85 (the smeared contact) | ≈ 1 (± 0.25) |
| `fan_staircase` (explicit Euler in the turning angle) | N = 20 … 320 | M 1.00 · p 1.00 · weak-wave product 1.01 | 1 |
| `ackeret_airfoil` thickness integral (second-order slopes, Simpson) | n = 51 … 401 | 1.97 | 2 |
| Becker profile: implicit relation against `solve_ivp` (rtol 1e-11) | — | agreement 1e-8 | — |
| limits measured as orders (closed forms, not discretisations): linearisation error 2.00 · entropy jump against M₁² − 1: 2.997 · plateau departure of ṁ(p_B): 2.00 · Hugoniot − isentrope: 3.0 · plate C_L, C_D against shock–expansion: 2.00 · one surface: 1.02 | | | |

## Conservation / invariant residuals
| Check | Residual |
|---|---|
| normal shock, three fluxes (4 γ × 8 M₁, up to M₁ = 40) | 3.3e-16 |
| isentropic nozzle: ρuA / h₀ / entropy along x | < 1e-9 / < 1e-12 / < 1e-11 |
| wall force (15.24) against ∫p dA (quadrature) | < 1e-8 of p₁A₁ |
| Laval nozzle, mass flow at 801 stations, across the shock | 2.1e-15 |
| duct with friction and heat, (15.45) on 91 roots | < 1e-11 |
| oblique shock, four laws in two dimensions (135 shocks) | < 2e-11 |
| moving shock, laboratory frame (mass / momentum / energy) | < 1e-12 / 1e-9 / 1e-9 (relative to the jump) |
| Becker profile: momentum flux with viscous stress | < 1e-10 |
| expansion fan: r ρ c between rays; streamline direction | 1e-12; 1e-6 rad |
| acoustic energy while two pulses cross | ptp/E = 4.0e-16 |
| simple wave: u − 2c/(γ − 1) | < 1e-10 c₀ |
| `riemann_exact`: jump laws at each shock; global balance of mass, momentum, energy (Gauss–Legendre between the wave edges) | < 1e-10; < 1e-11 |
| `euler1d_run`, periodic: mass / energy / momentum | 4.4e-16 / 2.2e-16 / 4.2e-18 (also with a captured shock: < 1e-13) |

## Benchmarks used (`reference/ch15/benchmarks.json`, `SOURCES.md`; all read on 2026-10-13)
| Value | Our value | Source | Verified |
|---|---|---|---|
| Table I rows M = 0.50, 0.55, 0.60 (seven columns) | 20 of 21 entries to half a unit of the last digit; A/A* at M = 0.55 printed 1.2550, ours 1.254948 (one unit) | Ames Research Staff, NACA Report 1135 (1953), Table I, p. 633 — public domain; https://www.nasa.gov/wp-content/uploads/2023/03/equations-tables-charts-compressibleflow-report-1135.pdf | first-hand, page image |
| Table II rows M = 1.60, 1.65, 1.70, 2.00, 2.05, 2.10 (15 columns: isentropic, ν, μ, normal shock) | 89 of 90 entries to half a unit; ν at M = 1.60 printed 14.861°, ours 14.86035° (one unit) | same report, Table II, p. 634 | first-hand, page image |
| Equations (105), (106), (111), (124), (127), (132), (137), (138), (139b), (142), (145), (146), (147), (167), (168), (169), (170), (171c); weak-shock series in ξ; expansion series; wedge limit 45.6° | agree at our own inputs to 1e-10 … 1e-13 (the two near-sonic approximations to their stated order) | same report, pp. 620–626 | first-hand, page images |
| normal-shock rows M₁ = 1.5, 2.0, 2.5, 3.0 | 5e-5 absolute | Wikipedia, "Normal shock tables" (secondary) | fetched |
| p*/p₀ = 0.528; minimum pressure ratios of eight gases; ṁ form | 5e-4; 6e-4 relative; 1e-12 | Wikipedia, "Choked flow" (secondary) | fetched |
| Fanno forms (dM²/M², 4fL*/D, ratios, Δs/c_p) | 1e-12 (forms); friction length = quadrature of the corrected (15.47) to 1e-9 | Wikipedia, "Fanno flow" (secondary) | fetched |
| Rayleigh forms; maximum of T at M = 0.845 | 1e-13 (forms); 1/√1.4 = 0.84515 | Wikipedia, "Rayleigh flow" (secondary) | fetched |
| Taylor's similarity constant (reused) | crude sphere constant 9.4 % low | `reference/ch01/SOURCES.md` | (ch01) |

The analyst's row M = 1.90 was not re-read and is not stored. The report's charts were not digitised (its equations give
the same information exactly).

## Numbers from the text (private values redacted to relative differences)
- Table 15.1 (251 rows): every entry within half a unit of the last printed digit except the two misprints the analyst
  found (slip #9) — confirmed.
- Table 15.2 (101 rows): one misprint (slip #10, 2.3 units) — confirmed; **six further entries are 0.51–0.56 of a unit
  off** (the digit is truncated rather than rounded). The analyst's "100 of 101 rows" holds at 0.6 of a unit, not at
  strict rounding.
- Example 15.1: the seven printed results within 1 % (the book reads a table to three digits); the analyst's exact values
  to 1e-9. Example 15.2: six printed results within 1 %; exact values to 1e-9; the same case is recovered from its back
  pressure by `laval_shock_location` (M₁ to 1e-7).
- §15.1–§15.2: gas constants, transonic band, the incompressible limit, ρ₀c² (2 %), c at 15 °C (0.5 %) — reproduced; the
  book's hypersonic threshold differs from the customary default of `mach_regime` and is reached with the keyword (T3).
- §15.6: asymptote of Fig. 15.8 (5e-4); the value read off Fig. 15.9 (2 %); the plateaus of Fig. 15.11 are those of
  M₁ = 2 (1e-6, 1e-5) and not of the captioned Mach number (off by > 0.02) — slip #8 confirmed; the thickness estimate.
- §15.8: the analyst's choking limits (1e-10) and differentials; Fig. 15.16's supersonic heating labels are the decimal
  values (slip #15 confirmed: with them the curves turn at M₁ = 1.4 … 1.9, with the printed whole numbers at 2.6 and
  beyond, the last outside the axis).
- §15.9: δ_max at M₁ = 2 (0.5°); the analyst's δ_max, σ, ν and exact/linear ratios (1e-8 … 1e-10).
- Exercises 15.3, 15.7, 15.8, 15.9, 15.12, 15.13, 15.14, 15.15, 15.17: the analyst's values to 1e-9 … 1e-12; Exercise
  15.7's quoted loss is 3.5 % above ours (slip #17, unresolved, pinned between 3 % and 4 %); Exercise 15.9's crude
  constant is 4.0 % above the constant the exercise quotes.
- Rule 9: none of the strings or patterns the private file lists as forbidden occurs in `tests/test_ch15.py`,
  `tests/ch15_verify_figures.py` or `reference/ch15/`; `tools/check_public.py` is clean on the tracked tree (769 files) and on the six new public files of this phase.

## Printed slips — confirmed or refuted
| # | Slip | Verdict |
|---|---|---|
| 4 | (15.32) third member lacks its minus sign | **confirmed** (finite differences, sympy; the printed option fails) |
| 11 | (15.47) numerator 1 − (γ − 1)M² | **confirmed**: 1 + (γ − 1)M² (sympy from (15.45); finite differences of (15.46)) |
| 12 | (15.48) second relation, sign | **confirmed**: (γ − 1)M²/(1 − M²) |
| 8 | Fig. 15.11 caption | confirmed (private) |
| 9, 10 | table misprints | confirmed (private); six truncated entries in Table 15.2 in addition |
| 15 | Fig. 15.16 labels | confirmed (private) |
| 17 | Exercise 15.7 | discrepancy reproduced, unresolved |
| 5, 7 | constant-density remark; structure equation before scaling | consistent with the code (M → 0 gives du/u = −dA/A; the quadratic is V²(U₁ − U)(U − U₂)) |

## Mutation check (do the tests discriminate?)
The package was copied to `outputs/ch15/_verify_scratch/mut*/` (three workers; the private book file not copied, so only
public tests can kill) and one wrong variant planted at a time — **103** in all: 59 in `core/gasdynamics.py` (specific
heats, sound speed, the exponents of (15.29)–(15.31), signs of (15.33), (15.34), each of the three printed slips made the
default, (15.39)–(15.44), the upstream guard removed, exp(−ΔS/R), the Hugoniot, the moving shock, the choking threshold,
p_d, the forgotten growth of A* behind the shock (three places), A, the discriminant, one sign of (15.46) for all M₁, the
removable 0/0 unhandled, signs of f and q, f_max, q_max, the heating band, Fanno and Rayleigh forms, tan for sin in
(15.49), (15.52), (15.53), the maximum-deflection and sonic loci, detachment not raised, (15.55), (15.56), the sign of
the turn, ν_max, the Riemann solver), 19 in `core/acoustics.py` (Newton's c, the cross term of (15.12), the leapfrog
step, the direction of g, the sign of g in (15.16), of (15.14), the initial-value formulas, (15.17), wavelet speed,
breaking time, front radius, cone delay, zone of action, Doppler, decibels) and 25 in `ch15_compressible_flow.py`
(regime thresholds, thrust, the rocket and shock chains, the blast constant, Becker's roots, ODE, profile and thickness,
fan geometry, staircase sign, curved wall, Ackeret's forms, the plate's two surfaces and its force resolution, the lift
slope, the Euler energy flux and smoothing, wedge detachment). The unmutated copy passes (93 passed, 5 skipped, 3
deselected). **103 of 103 were killed in the first pass.** 24 of them were stopped first by a contract test (23 by the
check of the design's printed numbers, one by the float-return check); those 24 were re-run with the five contract tests
deselected and **24 of 24 were killed by a physics test**. No assertion was changed between or after the passes (afterwards one test was renamed from V3 to V7 and one literal argument of an error-path call was changed to satisfy `tools/check_public.py`).

**Second pass (after fix loop 1): 45 of 46 killed; the survivor is equivalent.** Same harness
(`outputs/ch15/_verify_scratch/mutate2.py`, three workers, public tests only, slow script test deselected, nothing else
deselected; the unmutated copy passes: 101 passed, 5 skipped, 1 deselected), all 46 mutants in the six changed functions
of `core/gasdynamics.py`: `influence_coefficients` 18 (13 physics — the sign of (15.34) in the driver, the heat sign in the
energy driver, each printed slip made the default, the signs of the mass and momentum laws, γ for γ − 1 in the energy law,
dM/M, the three regular keys' factors, u du/h — and 5 of the new sonic-point contract: 0/0 returned as 0, F1 put back for
du/u, the infinity taking the supersonic side's sign, the printed (15.32) pure area change NaN at M = 1, the three regular
keys NaN at M = 1); `oblique_shock` 8 guard mutants (lower limit not checked, upper limit not checked, scalar does not
raise, no round-off allowance, allowance 10⁶ times too wide, `attached` judged from σ alone, refused entries keep their
σ, raise only for 0 < σ < μ₁); `oblique_shock_components` 7 (no lower limit, σ > π/2 accepted, the Mach wave keeping
M_n1 = 1 − ulp, keeping its negative round-off δ, geometry keys not NaN, (15.52) divided by sin σ, u₂ = u₁ρ₂/ρ₁);
`oblique_shock_curve` 2 (M₁ < 1 not refused, first point not the Mach wave); `normal_shock_entropy` 5 (the old M₁² − 1,
m = 2(M₁ − 1), the (15.43) exponent, the series without its cubic term, the series cut after the quartic);
`prandtl_meyer_turn` 6 (no clamp, allowance 10⁶ times too wide, two-sided clamp, every over-compression clamped, the sign
of the turn, the pressure exponent). 19 were first stopped by one of the six new tests, 5 by the contract check of the
design's printed numbers — those 5 re-run with the contract tests deselected were **5 of 5 killed by a physics test**
(`test_oblique_shock_V4_jump_laws_in_two_dimensions`, `test_normal_shock_V5_naca_1135_and_cited_table_rows` (2),
`test_weak_wave_V1_exact_shock_and_expansion_approach_the_linear_rule_from_both_sides` (2)). The survivor, "refused
entries keep their σ" (`sig = np.where(bad, np.nan, sig)` removed from `oblique_shock`), is **equivalent**: the guard of
`oblique_shock_components` uses the complementary condition (`ok` = not `bad` for every non-NaN input), so M_n1 is NaN
exactly there and the next line `sig = sig + 0.0 * Mn1` makes σ NaN anyway — the removed line is redundant. Note that
the first-pass F2 test alone would have let three of the guard mutants through (it accepts NaN *or* `attached=False`);
the new test pins the documented behaviour (raise for scalars). The mutation trees were deleted when the runs finished;
the scripts, `mut2_results.json`, `mut2_results_strict.json` and the log stay in the scratch folder (git-ignored).

## Label audit (second pass) — do the 145 rewritten docstring labels hold?
Every public function and class of the three modules (81 + 29 + 35 = 145) now has a "Validation: V… — tests/test_ch15.py:
… Label: …" block; none says "planned" or "pending". Audit, two layers:

1. **Mechanical, all 145** (`outputs/ch15/_verify_scratch/label_audit.py`): every cited test exists (145/145); for 144 of
   145 at least one cited test calls the function by name (the 145th, `ch15.ackeret_flat_plate`, is called by its other
   cited test — one citation is wrong, below); every claimed level appears in the name or header line of a cited test
   that calls the function, except two whose V3 lived in a test whose header did not say so (`ch15.fan_staircase`,
   `ch15.ackeret_airfoil` — the body measures the order; I added "V3" to the two headers and made the Ackeret order
   assertion two-sided); every label word matches a claimed level except nine "analytic" labels on V7-only claims.
2. **By reading the assertions, 88 functions**: all 52 that the nine explainers mirror (every `ch15.<name>` in the
   design's Part B `py:` rows) and 36 more (all six changed functions, the numerics `riemann_exact`, `euler1d_run`,
   `wave_equation_step`, `convected_wave_residual`, the Becker profile, the blast constant, the nozzles, the duct, the
   oblique-shock family, the tables). For each I checked that a cited test asserts the function's output against
   something independent at the claimed level (a closed form written in the test, sympy, a conservation law written in
   the test, a cited value, an order measured).

**Over-claims (claimed → supported):**

| Function | Claimed | Supported | Why |
|---|---|---|---|
| `AC.breaking_time` | V1, V4 · analytic, conserved | V1 · analytic | in the cited V4 test it is only compared with the closed form 2/((γ + 1)ak); no invariant or conservation law tests it |
| `GD.isentropic_state` | V1, V2, V5 · analytic, benchmark, symbolic | V1, V5 · analytic, benchmark | in the cited V2 test (D11) it only supplies the static state for a ρuA comparison; no sympy expression is compared with its output |
| `ch15.illustrative_inputs`, `conventions_table`, `book_slips`, `traps` | V7 · analytic | V7 (structure) | tables of inputs and strings: the test checks shapes, columns, IDs and that the inputs lie in their regimes — there is no analytic truth to check; the vocabulary has no word for V7, so the honest label is the level itself |

**Level text that understates (label earned, V-level mis-stated — not over-claims):** `GD.duct_role`,
`GD.rayleigh_temperature_sign`, `GD.laval_area`, `GD.strong_shock_limits`, `ch15.flat_plate_wave_pattern` say "V7"
but are labelled "analytic"; the tests assert exact values (the classification from the theory, the area at the inlet,
throat and exit and its zero slope at the throat, the limits against `normal_shock(1e6)`, the leading-edge waves against
the exact oblique shock and Mach angle), so "analytic" is earned and the level should read V1 (+ V7).
**Citation slips:** `ch15.ackeret_flat_plate` cites `test_flat_plate_V7_…` which does not call it;
`ch15.ackeret_pressure`'s V1 rests on `test_flat_plate_V1_linear_theory_formulas` (not cited) — the cited weak-wave test
compares it with another chapter function only; `GD.area_ratio`'s V2 rested on the D18 route only — the D10 test now
ties the code to the sympy expression (added this pass). The new tests of this pass are not yet cited in the docstrings
of the functions they pin (`influence_coefficients`, `normal_shock_entropy`, `oblique_shock`,
`oblique_shock_components`, `oblique_shock_curve`, `deflection_angle`, `prandtl_meyer_turn`).

### Label audit (third pass, after review Should-fix 8)
Re-run of `outputs/ch15/_verify_scratch/label_audit.py` on the working tree: 147 of 147 functions carry a Validation line;
cited tests that do not exist **0**; label words without a claimed level **0**; functions no cited test calls **0**.
All second-pass over-claims are gone (`AC.breaking_time` → V1 · analytic; `GD.isentropic_state` → V1, V5; `prandtl_meyer`
drops V4/conserved; the four tables say "V7 (structure)"). Residuals, none an over-claim:
- **Vocabulary**: the four tables' label "V7 (structure)" is not one of the eight vocabulary words; it is honest (there is
  no analytic truth in a table) and the vocabulary has no word for V7. Cosmetic.
- **Citation only**: `duct_role`, `rayleigh_temperature_sign`, `laval_area`, `strong_shock_limits` now claim "V1, V7"; the
  V1 is carried by the new `test_small_helpers_V1_against_independent_closed_forms`, which their docstrings do not cite
  yet. `flat_plate_wave_pattern`'s V1 is carried, in substance, by exact comparisons (to 1e-14 with
  `flat_plate_supersonic`, `oblique_shock` σ and `mach_angle`) inside a test named V7. `ch15.ackeret_flat_plate` still cites
  `test_flat_plate_V7_…`, which does not call it (its other cited test does).

## Figures reproduced with our code (`tests/ch15_verify_figures.py` → `outputs/ch15/verify/`, local)
| Figure(s) | File | Visual verdict |
|---|---|---|
| Fig. 15.1; steepening; Table 15.1 as curves; Figs. 15.5–15.6; 15.13; 15.12 | `A_sound_and_nozzles.png` | A 20 Pa pulse released from rest splits into two 10 Pa halves with antisymmetric velocity; the finite wave reaches a vertical tangent at t_b at mid-wavelength and folds after; the three ratios fall monotonically from 1 and A/A* has its minimum 1 at M = 1; four solutions cross or avoid M = 1 exactly at the throat; the pressure curves leave the supersonic branch through a shock that moves downstream as p_B falls and coincide below p_d; ṁ(p_B) meets its plateau tangentially at p*. |
| Figs. 15.8, 15.9, 15.11, 15.15, 15.16; T–s lines | `B_shocks_and_ducts.png` | M₂(M₁) passes (1, 1), tends to 0.378 and diverges at M₁ = 0.378; the entropy jump is negative below 1, touches zero flat at 1 and reaches about 1.1 at M₁ = 5, the cubic hugging it near 1; the Becker profiles fall between their two plateaus and steepen with M₁; friction curves turn round on M₂ = 1 and leave a choked gap that widens with f; heating curves do the same, cooling curves avoid (1, 1); both T–s lines have their largest entropy at the sonic point and the heating line a temperature maximum before it. |
| Figs. 15.19, 15.17, 15.24; ν(M); linear against exact | `C_planar_waves_and_plate.png` | δ(σ) starts at the Mach angle, ends at 90°, peaks at 23° for M₁ = 2 and 45.6° for M₁ → ∞, with the sonic locus just left of the maximum locus; ν rises from 0 to 130.5° (γ = 1.4), 180° (1.25), 90° (5/3); the exact/linear ratios leave 1 linearly, upward for shocks and downward for expansions; C_L and C_D follow the dashed linear laws up to a few degrees; the fronts of the supersonic source are tangent to the drawn cone. |

## Deviations & justifications
| `# DEVIATION` | Where | Verified by |
|---|---|---|
| leapfrog scheme for (15.13) is ours | `AC.wave_equation_step` | order 2.00, exact at C = 1, stability limit raised |
| MacCormack's scheme for the Euler equations is ours (ch10's step) | `ch15.euler1d_run` | order 1.98 on a smooth wave, conservation to round-off |
| pressure-switched artificial viscosity is ours | `ch15.euler1d_run` | L¹ order 1.0 at a captured shock, shock position within 4 cells at n = 800, still conservative |

Formulas labelled "ours" (thresholds of the Laval nozzle, choking limits, Fanno/Rayleigh forms, δ_max and the sonic
locus, the fan streamline, the blast constant) are each derived in a test or equal a first-hand form of NACA Report 1135.

## Open items
Resolved in this pass: first-pass items 1 (F1, F2 — fixed and pinned), 2 (labels written; audit above), 3 (both untrue
docstring statements corrected), 5a (`normal_shock_entropy` precision) and 5b (`prandtl_meyer_turn` return to M = 1).

Third pass: second-pass items 1 (over-claims), 2 (`deflection_angle` sentence) and 3 ("exactly zero") are resolved by
the review fixes; what is left of item 1 is citation housekeeping.

Remaining (none blocks the verdict; none is a failing test):
1. **Docstring citations (implementer, docstrings only).** Cite `test_small_helpers_V1_against_independent_closed_forms`
   in `duct_role`, `rayleigh_temperature_sign`, `laval_area`, `strong_shock_limits`;
   `test_convergent_nozzle_mass_flow_V1_one_gamma_per_call` in `mass_flow_vs_back_pressure`;
   `test_stagnation_pressure_V1_incompressible_wrapper_is_not_15_29` in `ch15.stagnation_pressure`;
   `test_rocket_nozzle_V7_warns_without_a_gas` in `rocket_nozzle`; drop the V7 citation of `ch15.ackeret_flat_plate`.
2. (resolved) 3. (resolved)
4. **`blast_constant_crude` is `qualitative`.** Its sphere constant is 9.4 % below the cited similarity value (the
   docstring now says about 9 %); the closed form itself is verified against an integration of its own ODE. Note for the
   notebook: `ch01.blast_radius` writes D = (E t²/(Kρ))^{1/5}, so its K is this chapter's K to the power −5.
5. **`wedge_flow` detached case and `flat_plate_wave_pattern` trailing-edge closure** are sketches, as their docstrings
   say: only the attached case, the centre-line normal shock and the leading-edge waves are verified (p₄ returns to p∞
   to O(α³)).
6. **`duct_outlet_mach` at f_max / q_max** — one-sided clamp of the discriminant, |M₂ − 1| ≤ 1.3e-8 (ruling 2; now stated
   correctly in the docstring). Could-fix only.
7. **Housekeeping (orchestrator):** `outputs/ch14/_resume_scratch/mut/` is still on disk (ruling 5); it breaks a bare
   `pytest -q` from the repo root and holds a copy of a private book file.
8. **Full suite** not re-run in this pass (the orchestrator runs it at the merge gate); the first-pass run had chapters
   1–14 and the machinery unchanged, and this pass touched only `tests/test_ch15.py` and `reports/ch15_verification.md`.

## Found wrong in design / analysis
- **Design Part F, D17, "Check"**: "at M₁ = 1.1 … 6.0 × 10⁻⁴ against the exact 4.2 × 10⁻⁴" — the exact value is
  4.290 × 10⁻⁴ (it should read 4.3 × 10⁻⁴). The other numbers of that line (7.0 × 10⁻⁵ against 5.9 × 10⁻⁵, the quartic
  coefficient) are right.
- **Analysis §6 T29 and the docstring**: the plate's linear theory errs at second order in α, not first (ruling 4).
- **Analysis §6 T19**: "crude blast constant within a few per cent (5 %)" — 9.4 % against the cited sphere constant.
- **Analysis §9a #10 / T17**: "100 of 101 rows of Table 15.2" needs the qualifier "to 0.6 of a unit"; six entries are
  truncated.
- **Design Part C counts**: on disk `GD.__all__` has 82 names, `AC.__all__` 30 and `ch15` 37 own (35 callables + `AC`,
  `GD`); the design says 81, 29, 35. All are covered.

## Rest of the suite
`.venv/Scripts/python.exe -m pytest -q tests`: 2080 collected; **2078 passed, 2 failed (both in `tests/test_ch15.py`),
0 skipped, 1 warning (a scipy `IntegrationWarning` in `tests/test_ch13.py`) in 2777.5 s**. `tools/check_public.py`: OK (769 tracked files; the six new files checked separately).
Scripts: the eight `scripts/ch15_*.py` and `reference/ch15/make_refs.py` exit 0. The full-suite run preceded the two cosmetic edits of `tests/test_ch15.py` named under "Mutation check"; the file was re-run alone afterwards (99 passed, 2 failed, 61.4 s).

Second pass: `.venv/Scripts/python.exe -m pytest -q tests/test_ch15.py -p no:cacheprovider` → **107 passed in 85.8 s**
(the slow script test included, so the eight scripts and `make_refs.py` ran again; the five private V6 tests ran).
`tests/test_machinery.py` → **9 passed in 23.4 s**. `tools/check_public.py` → OK (769 tracked files) and OK on
`tests/test_ch15.py` and this report. The full suite was not re-run (Open item 8).

Third pass: `.venv/Scripts/python.exe -m pytest -q tests/test_ch15.py -p no:cacheprovider` → **111 passed in 127.9 s**
(0 skipped: the slow script test and the five private V6 tests ran). `tools/check_public.py` → OK (775 tracked files; a slow full scan of about 10–20 min while other jobs were running)
and OK on `tests/test_ch15.py` + this report.
Full suite not re-run (merge gate).

## Verdict: PASS
