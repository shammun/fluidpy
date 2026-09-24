# ch07 — Gravity Waves: derivation & code review (phase 6)

Reviewer: `derivation-reviewer` (fresh context, read-only), 2026-09-24. Report written by the orchestrator from the reviewer's reply.
Inputs: `fluidpy/core/waves.py`, `fluidpy/ch07_gravity_waves.py`, `fluidpy/core/similarity.py`, `scripts/ch07_*.py`,
`tests/test_ch07.py`, `analysis/ch07_design.md` Parts C/F, `analysis/ch07.md` §9, `reports/ch07_verification.md`.
Pages compared symbol by symbol: (7.62), (7.65)–(7.67), (7.82)–(7.87), (7.93)–(7.96), (7.103)–(7.119), (7.143)–(7.158), (4.105).
Re-run: `pytest tests/test_ch07.py -m "not slow"` → 129 passed, 1 deselected (33.6 s).

**Verdict: 2 must-fix, then PASS** (status of each fix recorded at the end).

## Must fix

**M1 — negative k in the complex-step group velocity** (`core/waves.py` `_abs_k`, reached by `group_velocity_numeric`,
`_disp_fn`, `gaussian_packet`). For complex input `_abs_k` returned k unchanged, so for k < 0: deep water
`group_velocity_numeric(k=-1)` = 3.1e20 m/s, `gaussian_packet(k0=-0.1)` c_g = 9.9e19 (no warning); finite depth the sign
disagreed with `group_velocity`. From $\omega^2 = gk\tanh kH$ (7.28) with $|k|$,
$d\omega/dk = \operatorname{sgn}(k)\,\tfrac{c}{2}\left(1+\tfrac{2|k|H}{\sinh 2|k|H}\right)$, i.e. (7.69) with a sign.
Fix: `np.where(np.real(k) < 0, -k, k)` for complex input; one signed convention for `group_velocity`; tests at k < 0.

**M2 — γ = 1 attributed to Exercise 7.2** (`ch07_gravity_waves.py` `stokes_wave_speed` docstring,
`scripts/ch07_nonlinear.py` label, `analysis/ch07.md` D40/V37). Set up as printed (φ = a(ω/k)e^{kz} sin θ,
η = a cos θ + αka² cos 2θ, kinematic condition truncated at (ka)¹) the exercise gives γ = 3/8 in
$c^2 = \tfrac{g}{k}(1+\gamma k^2a^2)$. The book's (7.83), γ = 1, is right: it needs the potential's first-harmonic
correction −5/8·(ka)² relative to aω/k (−1/8·(ka)² relative to a√(g/k)), fixed only by the O((ka)³) kinematic balance.
Teaching decision: state (7.82)–(7.83) (Stokes 1847), show the consistent third-order sympy expansion, and add one trap
line explaining why the truncated Exercise 7.2 ansatz misses γ. Analysis text fixed by the orchestrator.

## Should fix
1. Overclaimed labels: `wavenumber_from_omega` "benchmark" (evidence is V1; FM/Guo are approximations checked against it);
   FM/Guo error figures → 1.7 % (Fenton & McKee 1990) and 0.79 % (Guo, exponent 5/2); `hydraulic_jump` "benchmark" → V1;
   "converged" labels on ka-scaling checks → analytic + V7 (asymptotic order).
2. Stokes limiting steepness 0.1410633: citation. On re-reading, HandWiki footnotes the digits to Dyachenko, Lushnikov & Korotkevich (2016) and "about one-seventh" to Schwartz & Fenton (1982); the docstring now states both (primary source of the digits unchecked).
3. `internal_pe_interface_limit`: "≈ 1.2 kε, 1.2 %" → 1.37 %.
4. Book slips added to `analysis/ch07.md` §9 as T10 (p. 288 interfacial E_p middle form /(2λ) gives ⅛; code ¼ is right)
   and T11 (text cites Ex. 7.16 for interfacial E_k; it is 7.18); module header to list them.
5. `ch07.reduced_gravity` re-exports ch04's ρ₁ form — not (7.117) $g' = g(\rho_2-\rho_1)/\rho_2$; point learners at
   `reduced_gravity_book`; fix the "(ρ₁+ρ₂)/2ρ₂" factor note for the mean form.
6. Sibling signature trap: `omega_capillary_gravity(k, H, sigma, rho, g)` vs `phase_speed/group_velocity(k, H, g, sigma, rho)`
   and different σ defaults. Public signatures stay stable this chapter; documented as a warning (candidate for a
   keyword-only change in a later machinery pass).
7. Docstring precision: `dyed_line` drift at the mean depth; `trace_velocities` returns a tuple.
8. Verification O1 (81 "Validation (planned)" docstrings), O2c–f, O9 (`np.errstate(over="ignore")` on complex tanh).

## Verified
- (7.145) for k < 0: sign-safe $\mathbf c_g = \nabla_{\mathbf K}(N|k|/K)$; printed form contradicts Fig. 7.31; all four quadrants tested.
- g′ with ρ₂ (7.117) forced by (7.113) as kH → 0; ch04's ρ₁ definition reachable via `ref=`, backward compatible.
- Interfacial E_p = ¼Δρga², E_k = E_p; two-layer constants (7.106)–(7.109), (7.114); interface (7.93)–(7.95), vortex sheet 2ωa cos θ.
- Polarization and energy (7.153)–(7.159), F = c_g E (also k < 0); standing wave (7.62); seiche (7.65); beats with ½Δω t.
- Orbits clockwise; Stokes drift (7.86) overflow-safe form identical; jump energy line; KdV (7.87)/(7.88); cnoidal speed.
- c_g for capillary–gravity waves; c_g,min from 3x² + 6x − 1 = 0; inverse-dispersion bracket; packet chirp; FFT sign.
- Tests catch flipped signs, halved coefficients, swapped densities in every CORE formula (35/35 planted variants).

## Fix status
Implementer fix round (2026-09-24): all fixed. M1 — `_abs_k` complex-safe, `group_velocity` returns signed dω/dk
(unchanged for k > 0; `phase_speed` stays |c|), 2 new V1 tests (odd symmetry, complex step = analytic to 1e-12,
closed forms −½√g, −(c/2)(1 + 4/sinh 4), −√(gH); packet mirror image and FFT agreement to 1e-12). M2 — docstrings and
script output corrected. Should-fix 1–3, 5–7, T10/T11, O1 (108 docstrings now name their tests), O2c–f, O9 — done;
`core/similarity.py` ratio note fixed by the orchestrator. `internal_pe_interface_limit` leading error now given
analytically as −2 ln 2·kε ≈ −1.39 kε. Should-fix 6 (keyword-only siblings) documented only; signatures unchanged.
Result: `tests/test_ch07.py` 132 passed; other chapters 750 passed (882 total); 12 scripts exit 0. **Review: PASS.**
