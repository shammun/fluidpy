# Chapter 12 review — Turbulence                       2026-10-07 · code at `d158c86`

Independent, fresh-context, read-only review by the `derivation-reviewer`; written up by the orchestrator from its reply.

## Verdict: no wrong equation and no wrong physics — 4 Must-fix items (a default, table constants, docstring labels, one citation)

Every numbered equation compared against the rendered pages matches symbol by symbol, and all six `# DEVIATION` sites are
justified. None of the Must-fix items changes a formula.

**Runs.** `tests/test_ch12.py -m "not slow"` under a line tracer: 114 passed, 3 deselected (179 s). Lines never executed by
those tests: `ch12_turbulence.py` 82 of 1487, `core/turbstats.py` 38 of 368, `core/wall_turbulence.py` 8 of 228.
`tools/check_public.py` OK (646 files) — but its forbidden list holds no table constants, so it cannot see M2.
Pages read as images: p573–583, 585–592, 594–597, 599–607, 610–622, 624–634 (not read: p584, 593, 598, 608, 609, 623 —
figures and prose without a numbered equation — and p635–647, exercises and references).

## Must fix

### M1 — `gradient_richardson_thermal` defaults to $\Gamma_a=0$
`fluidpy/ch12_turbulence.py` — `gradient_richardson_thermal(dTdz, dUdz, alpha, g=G0, Gamma_a=0.0, …)`. The docstring promises
Ri from the in-situ temperature gradient, but the default treats the input as a potential-temperature gradient. With
$\mathrm{Ri}=N^2/(dU/dz)^2=\alpha g\,(d\bar T/dz)/(dU/dz)^2$ (12.108), $\bar T$ potential, the in-situ form is
$N^2=g\alpha\,(dT/dz-\Gamma_a)$ with $\Gamma_a=-g/C_p\approx-9.76$ K/km (meteorological: $\Gamma_d=+9.76$ K/km, stable when
$\Gamma_{met}<\Gamma_d$). The standard atmosphere, `gradient_richardson_thermal(-6.5e-3, 0.01, 1/288)`, returns Ri = −2.213,
"unstable", and prints "Γa = 0.0 K/km"; with `Gamma_a=-9.76e-3` it gives Ri = +1.110, "stable". The test suite already
names this call a mutant. **Fix:** make `Gamma_a` a required keyword; callers holding a θ-gradient pass `Gamma_a=0.0`
explicitly; parity rows pass `adiabatic_lapse_rate()` rather than a typed −0.0098.

### M2 — Book Table 12.1 constants in committed files
(Values deliberately not repeated here — this report is public.)

| Where | Which "illustrative" constant | Coincides with |
|---|---|---|
| `scripts/ch12_free_shear_table.py` | the `C_Y`, `C_U` pair | the round-jet pair |
| `scripts/ch12_drawings.py` | `dxi80_coeff` | the shear-layer coefficient |
| `scripts/ch12_common.py` | `ILLUSTRATIVE_JET` `C_Y`, `xi_half_Y` | plane-jet scalar constant; planar-plume scalar half-width |
| `analysis/ch12_design.md` (A.15a) | `xi_half_Y` of the mass-fraction row | plane-jet scalar half-width |
| `analysis/ch12_design.md` (A.15a) | the B typed into `nagib_chauhan_kappa(…)` | the ZPG (κ, B) pair |

**Fix:** numbers visibly not the table's (C_U = 7, C_Y = 4, dxi80 = 0.1, xi_half_Y = 0.15), a round trip
`nagib_chauhan_kappa(nagib_chauhan_B(0.384))` for $\kappa B=1.6[\exp(0.1663B)-1]$ (12.92), and the patterns added to the
public checker. Coinciding but cited (a decision, not a violation): `WT.LOG_LAW_CONSTANTS["classical"]` = (0.41, 5.0) and
κ = 0.384 (Lee & Moser). No table constant was found in the tests, the three modules or `reference/ch12/`.

### M3 — Validation statements in docstrings that are false or stale
1. `mixing_length_intercept`: "within 5 %" is false — B(0.41, A⁺ = 26) = 5.277, 5.5 % above 5.0. Say "approximate (model)".
2. `skin_friction_zpg`: the three laws agree within 10 % to Re_x = 10⁸ and 12.0 % at 10⁹.
3. `stress_partition`: equal shares at y⁺ ≈ 9.9, not 11; label it approximate.
4. `spalding_uplus`: no "approximate" label (≤ 0.85 wall units against DNS).
5. `k_epsilon_channel`: qualitative — converges, no grid study, no DNS comparison.
6. `model_spectrum`: c_L = 6.78 is verified; drop "unverified".
7. `explainer_tables`: "to the 6 digits stored" holds only at un-rounded heights (≈ 3e-5 at the stored ones).
8. `rans_eddy_viscosity_residual`: the printed (12.97) with a lone j on $\partial P/\partial x_j$ is ill-formed; the sum is one arbitrary realisation, not "the summation convention".
9. `channel_mixing_length`: production peak at y⁺ ≈ 10.2–10.5 with value 0.222–0.249, not "≈ ¼ near y⁺ ≈ 12".
10. `dispersion_local_slope`: 1.370 at x = 3.2 and 1.301 at 4.0, not "1.3 at 3.2".
11. `dispersion_rate_from_particles`: cites "(12.115) dX/dt = u"; the page's (12.115) is $\frac{d}{dt}\overline{X_\alpha^2}=2\,\overline{X_\alpha\,dX_\alpha/dt}$ — $u_\alpha=dX_\alpha/dt$ is the unnumbered definition after it.

### M4 — The Businger–Dyer "benchmark" rests on a citation nobody has read
`core/wall_turbulence.py` attributes $\phi_m=(1-16\zeta)^{-1/4}$ with stable coefficient 4.7 to the AMS Glossary;
`reference/ch12/SOURCES.md` says the same entry gives $1+5\,z/L$ and that the page returned HTTP 403. The pair (16, 4.7)
belongs to neither published set as far as the reviewer recalls (Businger et al. 1971: 15 with 4.7; Dyer 1974: 16 with 5) —
unverified, from memory. **Fix:** a human reads the source. Until then the test is relabelled V1 (it proves ψ_m is the
integral of the coded φ_m — analytic, not a benchmark) and the coefficient is marked "literature value, not verified
first-hand". The coded physics is self-consistent.

## Should fix
1. The key `uv_plus` holds **minus** $\overline{uv}$ (in `mixing_length_wall_profile`, `channel_mixing_length`, `channel_energy_budget_at`, `explainer_tables.json`), while `mean_energy_budget` and `tke_budget` take the true $\overline{uv}$. Add the alias `minus_uv_plus`.
2. `inertial_spectrum_1d`: `C1` is the one-sided constant and the default `two_sided=True` halves it; that combination is untested. Document or rename.
3. `WT.surface_layer_wind(unstable="log_linear")` returns numbers (even negative winds) where `dimensionless_shear` returns NaN (ζ ≤ −1/β). Return NaN there.
4. `TS.integral_scale(upto="first_zero")` is not $\Lambda=\int_0^\infty r\,d\tau$ (12.18) when r has a negative lobe — exactly the transverse g(r). `scripts/ch12_isotropic.py` prints Λ_g/Λ_f = 0.62 beside the taught ½. Use `upto="all"` for Λ_g and warn in the docstring.
5. The meteorological verdict string writes "Γa" for +9.8 K/km where ch12 calls it Γ_d; explain once beside the badge.
6. `K_EPSILON_CONSTANTS` comment cites a page that returned 403; `SOURCES.md` cites others. Align.
7. `SOURCES.md` lacks rows for A⁺ = 26 (van Driest 1956), the core cap l/δ = 0.09 and the "classical" preset.
8. The Lee & Moser subset is redistributed without a licence statement (cited) — the user's decision.
9. `test_k_epsilon_V5_constants_…` and `test_inertial_spectrum_V5_…` check constants, not code: call them "constant consistency". The PASS is unaffected.
10. The mutant campaign ran on the default tests only; `isotropy_report`, `TS.longitudinal_transverse_correlation` and the fresh `explainer_tables()` call live in `@slow` tests.
11. Never executed by a default test (all behave correctly by hand): `periodogram(window="hann")`, `make_ensemble` with non-uniform t, `check_averaging_rules` with 3-D input, the numeric branch of `_derivative_of`, the r = 0 branch of `isotropic_correlation_tensor`, the linspace branch of `_wall_grid`, two fallbacks in `integral_scale` / `correlation_time`, and about 45 `raise` lines.
12. `WT.log_law_crossing` surfaces brentq's raw message when no crossing exists (e.g. undamped B = −1.23).

## Verified (against the rendered pages)
- **§12.3–12.4**: (12.1)–(12.23), including $\Lambda_t=\int_0^\infty r_{11}\,d\tau$ (12.18), $\lambda_t^2=-2/[d^2r_{11}/d\tau^2]_0$ (12.19), $S_e=\frac1{2\pi}\int R_{11}e^{-i\omega\tau}d\tau$ (12.20), $\overline{u_1^2}=\int S_e\,d\omega$ (12.22).
- **§12.5**: (12.24)–(12.35); the mean momentum equation (12.30) with Reynolds stress $-\rho_0\overline{u_iu_j}$; $Q_j=-k\,\partial\bar T/\partial x_j+\rho_0C_p\overline{u_jT'}$ (12.32).
- **§12.6**: (12.36)–(12.45); $\bar\varepsilon=-15\nu\overline{u^2}f''(0)=30\nu\overline{u^2}/\lambda_f^2=15\nu\overline{u^2}/\lambda_g^2$ (12.43); on an independent 48³ field ε̄ = 24.03 against 24.06.
- **§12.7**: (12.46)–(12.55); the budget closes by hand in wall units; $(18/55)(1.5)/2=0.2455$.
- **§12.8**: (12.56)–(12.75) and the Table 12.1 *forms* for all seven flows.
- **§12.9**: (12.76)–(12.93); h is the full height in $dP/dx=-2\tau_0/h$ (12.90); $\kappa B=1.6[\exp(0.1663B)-1]$ (12.92).
- **§12.10 (incl. p621)**: (12.94)–(12.105); undamped $B=(\ln4\kappa-1)/\kappa=-1.2324$; $\nu_T=C_\mu\bar e^2/\bar\varepsilon$ (12.104).
- **§12.11**: (12.106)–(12.114); $L_M=-u_*^3/\kappa\alpha g\overline{wT'}$ (12.110) gives +83.01 m for H = −30 W/m² and −16.60 m for +150 W/m².
- **§12.12**: (12.115)–(12.129); 2×10⁴ particles within one standard error of Taylor's closed form.
- **DEVIATION sites**: slip #6 ($\partial P/\partial x_i$ in (12.97)), slip #1 ($k_1^{-5/3}$ in (12.54)), $C_3G=-\tfrac12F\int_0^\xi F\,d\xi$ from (12.63) (independent outward integration agrees to 5.4e-7), signed mixing-length stress, slip #15 (the ½ in (12.112)), slip #12 ($t\gg\Lambda_t$ in (12.129)) — all correct. Slips #2–#5, #8–#11, #13 confirmed on the pages; #7 and #14 concern exercises not read.
- **Conventions**: Reynolds-stress sign, buoyant production sign, Rf and L_M signs, both lapse-rate verdicts (with Γ_a supplied), Ri = Pr_T·Rf, spectrum normalisation, one-component $\overline{u^2}$, λ_g = λ_f/√2, Λ_g = Λ_f/2, σ² = 2νt.

## Status
Fixes: see the "review fixes" entry in `progress.json` and the commit following this report.
