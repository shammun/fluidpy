# Chapter 16 — derivation review                                        2026-10-14

Independent, fresh-context, read-only review of `fluidpy/ch16_biofluid_mechanics.py`, `fluidpy/core/womersley.py`,
`fluidpy/core/rheology.py`, `scripts/ch16_*.py` and `tests/test_ch16.py`, checked line by line against the rendered pages
p817–p876. Re-run: `tests/test_ch16.py -m "not slow"` 90 passed, 1 deselected (118 s).

**Verdict: 1 Must fix (regime misuse), no formula, sign, factor or unit error found.**

## Must fix

1. **Collapsible-tube law used outside its fitted regime (α > 1), unguarded.** (16.188), $-P\approx\alpha^{-n}-1$, is
   Shapiro's fit to the collapse branch (P < 0); Fig. 16.19 rises almost vertically for α ≥ 1, and design D28 says the fit
   fails there. `scripts/ch16_collapsible.py:58` uses `al_t = 0.6`, so the frictionless subcritical branch upstream of the
   throat is distended (α = 1.501 at x = −0.12 m, 1.428 at −0.06 m, 1.208 at −0.03 m): P is capped at 1 and
   $c^2=nK_p\alpha^{-n}/\rho$ falls as the tube inflates — non-physical. The tube-law plot (`:38`) runs to α = 1.6.
   `collapsible_tube_law` (`ch16…py:2401–2415`), `collapsible_frictionless_alpha` (`:2705`) and `collapsible_terms`
   (`:2572`) have no range guard. **Fix:** keep every demo state at α ≤ 1 — upstream α ≤ 1 iff
   $H=K_pP(\alpha_t)+\tfrac12\rho c_t^2\le0$ iff $(1-n/2)\alpha_t^{-n}\ge1$ iff $\alpha_t\le0.25^{2/3}=0.397$ for n = 3/2
   (α_t = 0.35 gives H = −41.5 Pa, upstream α = 0.876); warn (arrays: NaN + warning, T16 pattern) for α > 1; docstring
   "valid for α ≤ 1 (collapse branch, Fig. 16.19)"; a test pins the warning. E9's throat mode must use the same inputs.

## Should fix

1. `womersley_flow_rate(t, …, 100.0)` treats a scalar as the n = 1 amplitude (4.58e-7 m³/s) while
   `womersley_flow_rate_harmonics(t, …, 100.0)` treats it as the mean G₀ (1.12e-4 m³/s) — same scalar, 245× apart
   (`core/womersley.py:384–389` vs `:451–461`). Make the meaning agree (or raise) and test it.
2. Positional order differs between siblings: `vascular_resistance(mu, L, a)` vs `poiseuille_flow_rate(a, dp, L, mu)`
   (a ↔ μ swap gives 1.7e7 instead of 8.9e4); `tube_radius_hookean(p, pe, a0, E, h)` vs `tube_law_hookean(A, A0, E, h, a0)`,
   `tube_law_linear`, `tube1d_pressure_flux`. Keyword-only or aligned order (deferred while the notebook and explainers
   are being built against the current signatures; see below).
3. `cardiac_output` docstring (`:200`): "multiply by 6e7 for L/min" — ×6e7 gives mL/min; L/min is ×6e4.
4. S10 justification (`core/womersley.py:298`): "J₁/J₀ → −i[1 − i/(2z) + …]" — numerically it is **+i**[1 − i/(2z)]
   (J₁/J₀ = −0.00712 + 0.99293i at z = i^{3/2}·50). The coded asymptote is right (slope 3); only the stated reason is wrong.
5. Labels that claim "converged" without an asserted refinement order: `elastic_tube_flow_rate` `:850`,
   `elastic_tube_pressure` `:873`, `windkessel_ode` `:990`, `elastic_tube_determinant_numeric` `:1636`, `bifurcation_point`
   `:1805`, `network_pulse_run` `:1934`, `collapsible_march` `:2666`, `starling_state` `:2816`, `flow_limitation_curve`
   `:2869` — these are independent-route agreements (V1 cross-checks). Tests named V3 without an order:
   `test_distensible_flux_V3_…` (673), `test_windkessel_V3_…` (748), `test_murray_V3_…` (1410), `test_junction_V3_…`
   (1482), `test_elastic_tube_V3_independent_4x4_determinant` (1273, asserts only large-α slopes), Dean c₄ "V3" (1581,
   single resolution; its strength is the 9e-6 agreement with the exact rational, V2). Relabel or add a resolution study.
6. Docstring regime notes: `tube_radius_hookean` accepts p < p_e (collapse) without a note; `pulse_state` (`:1294`)
   treats the pulse pressure as the full transmural pressure (a₀ = radius at mean pressure); `elastic_tube_wave` defaults
   c0 = 1, omega = 1 make k₁, k₂ meaningless unless both are set.

## Verified

- **Womersley (16.73)–(16.100):** S6 (right side −G/μ); U = Ga²/(iμα²)(1 − F₁); Q via ∫J₀(λr)r dr = aJ₁/λ → (1 − F₂);
  τ_w = (Ga/2)F₂ (S8); φ = arctan2(X, Y); small-α series re-derived; large-α (16.88); corrected S10 slope 3; Crank–Nicolson
  axis row, wall and forcing.
- **Elastic-tube waves (16.116)–(16.136):** β̂ = e^{iπ/4}α under e^{−iωt} (T1); numeric determinant reproduces (16.131);
  X = k²c₀²/ω² matches (16.135)/(16.136); (16.137) via E_ratio.
- **Tubes and pulse waves:** (16.41)/(16.55), compliance, (16.46) with a₀⁴ (S4); Moens–Korteweg; 1-D tube flux, c(A), and
  inlet Riemann invariant sign.
- **Two-layer model:** (16.27)–(16.31) with corrected (16.30); thin-layer (16.34).
- **Casson pipe:** (16.215), u_c, flux bracket, a_c ≥ a guard (T21); Herschel–Bulkley and power-law fluxes.
- **Windkessel:** systolic particular solution, p(T_s), p₀; p_D sampling error 6e-3 Pa.
- **Murray, junctions, network:** a⁶ = 16μQ²/(π²K) (S15); pumping = upkeep/2; weighted-Fermat angles, θ = arccos 2^{−1/3};
  R, T and the energy identity; junction node formula.
- **Dean:** K = 2κ², D = 4√2κ, λ inversion −0.01106, w₀ = 2(1 − r²); S34 confirmed independently (ζ = −∇²ψ gives + on the left).
- **Shapiro:** (16.194), (16.202), M = 2 − n, laminar R, α*.
- **Plants:** (16.221) SI and book units (S30); xylem root pressure; (16.222).
- **Units:** MMHG, PRU, cP (S27), dyn/cm², (16.5) ×0.1, Whitmore τ_y = 4μ_p.

## Resolution

**PASS (2026-10-14).**

- **Must fix 1 — fixed.** `collapsible_tube_law` and its callers warn for α > 1 (`ALPHA_MAX_COLLAPSE = 1.0`; scalar →
  RuntimeWarning + extrapolated value, array → NaN + RuntimeWarning); docstrings say "valid for α ≤ 1 (collapse branch,
  Fig. 16.19)". Demo throat α_t = 0.35 (limit $(1-n/2)^{1/n}=0.3969$): H = −41.47 Pa, upstream α 0.8758 at x = −12 cm,
  Q* = 47.646 cm³/s; tube-law panel stops at α = 1. Tests pin the guard and prove H ≤ 0 ⇔ α_t ≤ 0.3969 (sympy).
- **Should fix 1, 3, 4, 5, 6 — fixed.** `womersley_flow_rate_harmonics` raises for a bare scalar (scalar = n = 1
  amplitude elsewhere); ×6e4 for L/min; S10 reason +i[1 − i/(2z)]; "converged" labels → V1 cross-checks, tests renamed,
  and four genuine order studies added (quadrature 2 / 4, Windkessel CN 2.000, elastic-tube ODE determinant 6.0, Dean
  spectral resolution); regime notes added.
- **Should fix 2 — deferred.** Positional order unchanged while the notebook and explainers are built; each docstring
  states its order and names the sibling whose order differs.
- **Found during the re-test (L3-1) — fixed.** `collapsible_frictionless_alpha` raised for wide upstream sections
  (A₀/A₀,t ≳ 24, α* < 1e-6 outside the brentq bracket); now uses the closed form
  α* = (ρQ²/(A₀²K_p n))^{1/(2−n)} (both law forms share ∂P/∂α), raises clearly for n ≥ 2.
- `tests/test_ch16.py -m "not slow"`: 98 passed, 1 deselected; the slow script test passed in loop 3; `-W
  error::RuntimeWarning` on the collapsible and elastic-tube tests: 16 passed.
