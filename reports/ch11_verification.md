# Chapter 11 verification — Instability
2026-10-01 · repository state `8ae4f15` (design + implement) + uncommitted loop-1 and loop-2 implementer changes
(`fluidpy/core/stability.py`, `fluidpy/ch11_instability.py`, `scripts/ch11_double_diffusion.py`,
`scripts/ch11_inviscid_criteria.py`, `scripts/ch11_stratified_shear.py`, `reference/ch11/tg_growth_map.csv`,
`reference/ch11/explainer_tables.json`, `reference/ch11/rayleigh_spectra.json`) + uncommitted verifier files
(`tests/test_ch11.py`, `tests/ch11_verify_figures.py`, `reference/ch11/SOURCES.md` — the two "Verifier re-check" sections)
· verifier: `math-verifier` (loop 2 of max 3)

## Verdict: PASS — both Taylor–Goldstein defects (loop 0: decay box at small k; loop 1: false zeros under the neutral curve at k ≥ 0.65) are fixed and independently confirmed; 0 failing tests

- **`tests/test_ch11.py`: 128 test functions; full file (slow included): 128 passed, 0 failed.** Default run
  (`-m "not slow"`): **125 passed, 3 deselected** (384 s). Slow run (`-m slow`): **3 passed, 125 deselected** (197 s).
  Labels in the test names: V1 55 · V2 22 · V3 17 · V4 7 · V5 8 · V6 2 · V7 31 (a test that carries two labels counts twice).
- The file as the implementer saw it in loop 2 (119 test functions) was re-run unchanged first: **118 passed, 1 failed**
  (637 s, slow included). The loop-1 failing test `test_tg_growth_V3_no_false_zero_under_the_neutral_curve_at_large_k`
  **passes unchanged**. The one failure was the verifier's own pinned count (`== 299` positives in the table), which had
  encoded the loop-1 defect; it is corrected on evidence (see "Test changes in loop 2"), not loosened.
- **Full suite (`pytest tests`): 1588 passed, 0 failed in 1225 s — 1460 tests of the other chapters and the machinery + 128 of ch11 (slow included).****
- **Mutants: 25 of 25 planted wrong variants killed** (15 from loop 1 re-run + 10 new for the loop-2 changes; details at the end).
- `tools/check_public.py`: OK (591 files). Book values are used only through the git-ignored `tests/book_values_ch11.json`.
- All 10 runnable scripts `scripts/ch11_*.py --no-show --fast` run clean (slow test; `ch11_common.py` is a helper and
  `ch11_tables.py` is run only by `make_refs.py` because it rewrites `reference/ch11`); `reference/ch11/make_refs.py
  --no-tables` rewrites a byte-identical `benchmarks.json`. `scripts/ch11_inviscid_criteria.py` also runs clean with
  `RuntimeWarning` promoted to an error (the `_sech2` overflow of loop 1 is gone).

### Loop history
| Loop | Defect | Evidence | State |
|---|---|---|---|
| 0 → 1 | `tg_growth` imposed ψ̂ = 0 at a fixed y_max = 30 (1.5 e-folds of e^{−k|z|} at k = 0.05) | `test_tg_growth_V3_converged_in_the_decay_box_at_small_k`, `test_tg_growth_map_V1_cached_csv_matches_live_and_the_exact_tongue` failed | **fixed in loop 1** (`decay_box(k)` = max(30, 12/k)); both tests pass unchanged |
| 1 → 2 | the table reported 0 in a strip 0.0375–0.06 wide under the neutral curve for k ≥ 0.65 (critical layer not resolved) | `test_tg_growth_V3_no_false_zero_under_the_neutral_curve_at_large_k` failed | **fixed in loop 2** (`decay_map_scale(k)` = min(0.5, max(0.035, 0.025/k)) — the tan-map scale follows the wave; N = 100 unchanged); the test passes unchanged |
| 2 | — | 128 passed | **closed** |

The verifier's loop-1 hypothesis ("resolution in N") named the symptom; the implementer's fix shows the cheaper lever is
where the nodes sit (map scale), not how many there are: with the rule, N = 100 reproduces N = 200 / 240 and an
independent shooting solution to ≤ 1.7e-6 absolute.

**Loop 0 → 1, decay box** (kc_i; "converged" = box 240 … 960, N = 100 … 240, shooting on the unbounded layer — all equal
to the digits shown):

| k | J | loop 0 (y_max = 30) | loop 1 (`decay_box`, s = 0.5) | loop 2 default | converged / shooting | loop-0 error | loop-2 error |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.00 | 0.0455051 | 0.0457335 | 0.0457335 | 0.0457335 | −0.50 % | < 1e-7 |
| 0.05 | 0.02 | 0.0351371 | 0.0340758 | 0.0340758 | 0.0340758 | +3.1 % | < 1e-7 |
| 0.05 | 0.04 | 0.0205913 | 0.0163047 | 0.0163047 | 0.0163046 | +26.3 % | 1e-7 |
| 0.05 | 0.045 | 0.0152825 | 0.0082157 | 0.0082157 | 0.0082156 | +86 % | 1e-7 |
| 0.05 | 0.05 | 0.0077967 (spurious) | 0 | 0 | 0 (J > k(1 − k) = 0.0475) | spurious growth | exact |
| 0.10 | 0.08 | 0.0217285 | 0.0211428 | 0.0211430 | 0.0211428 | +2.8 % | 2e-7 |
| 0.10 | 0.085 | 0.0135108 | 0.0126792 | 0.0126794 | 0.0126792 | +6.6 % | 2e-7 |
| 0.15 | 0.10 | 0.0423453 | 0.0423204 | 0.0423205 | 0.0423204 | +0.06 % | 1e-7 |
| 0.15 | 0.125 | 0.0058937 | 0.0058231 | 0.0058248 | 0.0058244 | +1.2 % | 4e-7 (loop 1: 1.3e-6) |
| 0.30 | 0.10 | 0.1089743 | 0.1089743 | 0.1089743 | 0.1089743 | 2e-9 | < 1e-7 |

**Loop 1 → 2, false zeros at large k** (kc_i; N = 200 and N = 240 with the rule, half and twice the map scale at N = 200,
twice the box at N = 160 and the shooting solution all agree to the digits of the "converged" column; N = 240 and shooting
differ by ≤ 1.2e-10 at all 35 points probed):

| k | J | distance below J = k(1 − k) | loop 1 (s = 0.5, N = 100) | loop 2 default (N = 100) | converged / shooting | loop-2 error |
|---|---|---|---|---|---|---|
| 0.90 | 0.04 | 0.050 | **0** | 0.0332469 | 0.0332455 | 1.4e-6 |
| 0.90 | 0.05 | 0.040 | **0** | 0.0267117 | 0.0267103 | 1.4e-6 |
| 0.90 | 0.06 | 0.030 | **0** | 0.0201206 | 0.0201192 | 1.4e-6 |
| 0.90 | 0.08 | 0.010 | **0** | 0.0067666 | 0.0067652 | 1.4e-6 |
| 0.85 | 0.08 | 0.0475 | **0** | 0.0326887 | 0.0326878 | 9e-7 |
| 0.85 | 0.12 | 0.0075 | **0** | 0.0052635 | 0.0052627 | 8e-7 |
| 0.80 | 0.12 | 0.040 | **0** | 0.0286322 | 0.0286323 | 1e-7 |
| 0.80 | 0.13 | 0.030 | **0** | 0.0215916 | 0.0215917 | 1e-7 |
| 0.80 | 0.15 | 0.010 | **0** | 0.0072775 | 0.0072776 | 1e-7 |
| 0.95 | 0.01 | 0.0375 | **0** | 0.0242820 | 0.0242804 | 1.6e-6 |
| 0.95 | 0.04 | 0.0075 | **0** | 0.0049144 | 0.0049133 | 1.1e-6 (2.2e-4 relative) |
| 0.75 | 0.18 | 0.0075 | **0** | 0.0056917 | 0.0056929 | 1.2e-6 |
| 0.70 | 0.20 | 0.010 | **0** | 0.0079147 | 0.0079164 | 1.7e-6 (largest of the 35 points) |
| 0.65 | 0.22 | 0.0075 | **0** | 0.0062373 | 0.0062384 | 1.1e-6 |
| 0.60 | 0.23 | 0.010 | **0** | 0.0087377 | 0.0087383 | 6e-7 |
| 0.50 | 0.24 | 0.010 | **0** | 0.0098401 | 0.0098401 | < 1e-7 |
| 0.93 | 0.055 | 0.0101 (off grid) | **0** | 0.0066980 | 0.0066964 | 1.6e-6 |
| 0.97 | 0.02 | 0.0091 (off grid) | **0** | 0.0058840 | 0.0058825 | 1.5e-6 |
| 0.4449 | 0.00 | — | 0.1897021 | 0.1897017 | 0.1897021 | 4e-7 (Michalke point unchanged) |

Zero strip under the exact neutral curve in the published table, per column (the J step is 0.01):

| k | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 | 0.85 | 0.90 | 0.95 |
|---|---|---|---|---|---|---|---|---|
| loop 1 | 0.0300 | 0.0375 | 0.0400 | 0.0475 | 0.0500 | 0.0575 | 0.0600 | 0.0475 |
| loop 2 | 0.0100 | 0.0075 | 0.0100 | 0.0075 | 0.0100 | 0.0075 | 0.0100 | 0.0075 |

(0.0075 or 0.01 = the first grid point at or above the curve; in the columns k = 0.1, 0.2, …, 0.9 the curve passes
through a grid point, where kc_i = 0 exactly.) Positive entries: loop 0 300 (some of them spurious), loop 1 299, **loop 2
335 = exactly the 335 grid points with J < k(1 − k)**; none on or above the curve.

### Independent checks of the implementer's loop-2 claims (all confirmed unless stated)
1. **Large-k fix.** The verifier's failing test passes unchanged. Own N-doubling / map-scale / box variation / shooting at
   35 (k, J) points (tables above): default vs N = 200: ≤ 1.70e-6 absolute; vs N = 240: ≤ 1.70e-6; vs half the scale
   (N = 200): ≤ 1.67e-6; vs twice the scale: ≤ 1.70e-6; vs twice the box (N = 160): ≤ 1.71e-6; vs shooting: ≤ 1.70e-6
   (largest at (0.7, 0.2), 2.1e-4 relative).
2. **Whole map, N = 100 vs N = 200** (620 points each, 162 s and 965 s): identical zero pattern (335 positives); largest
   difference 2.05e-6 absolute (at (0.7, 0)), largest relative 2.2e-4 (at (0.95, 0.04)) — the implementer's numbers. Because
   the table is rounded to 4 s.f., the N = 200 map rounds to the published entry in 604 of 620 entries; the other 16 sit on a
   rounding boundary (e.g. (0.7, 0): 0.1460481 → 0.146, 0.1460502 → 0.1461), i.e. the table is exact to ±1 in its last digit.
3. **No spurious growth at or above the neutral curve.** 504 points with J = k(1 − k) + {0, 1e-4, 1e-3, 3e-3, 0.01, 0.03}
   and J = 0.25, 0.2501, 0.26, 0.3, 0.5, 1.0 at 41 wavenumbers 0.03 … 1.3 (32 of them off the table's grid; k ≥ 1 also at J = 0,
   0.01, 0.1): **0 non-zero**. Miles–Howard (11.67): a further 56 points with J = 0.25 … 2.0, k = 0.03 … 1.2: 0 non-zero.
4. **Michalke point**: `tg_growth(0.44492, 0)` = 0.1897017, `tanh_max_growth()` = (0.4449181, 0.1897021) — unchanged.
5. **N² = 0 parity with Rayleigh on the rule grid**: |c_TG − c_Rayleigh| ≤ 5.1e-11 at k = 0.05, 0.2, 0.4449, 0.7, 0.9, 0.95,
   0.99 (both solvers with `decay_map_scale(k)`, `decay_box(k)`).
6. **Howard semicircle**: every computed mode at (0.9, 0.04), (0.95, 0.04), (0.65, 0.22), (0.05, 0.04), (0.8, 0.15),
   (0.5, 0.24), (0.4, 0.1) and all 37 modes of the regenerated semicircle figure (5 J × 9 k: exactly one mode wherever
   J < k(1 − k), none elsewhere) lie inside.
7. **The stated near-neutral limitation is true as written.** Bisection in 13 columns: `tg_growth` returns 0 below the
   curve only within 0.00019 (k = 0.05), 0.00035 (0.1), 0.00070 (0.2), 0.00110 (0.3), 0.00123 (0.33), 0.00156 (0.4),
   0.00207 (0.5), 0.00264 (0.6), 0.00324 (0.7), 0.00423 (0.8), 0.00502 (0.87), 0.00536 (0.9), 0.00583 (0.95) of it
   (docstring: 0.0002 … 0.0059); the growth at the edge of that sliver is 0.00109 … 0.00378 (docstring: "at most 0.0011 …
   0.0038", "kc_i < 0.004"); the true growth in the middle of the sliver (shooting) is 0.0006 … 0.0019 — positive, so the
   0 there is a resolution limit, as the docstrings and the CSV header say. 25 log-spaced probes per column: no holes
   (0 only inside the sliver, growth everywhere outside).
8. **Tables.** `tg_growth_map.csv` = `explainer_tables.json → tg_map` entry for entry; a full live recomputation rounds to
   the CSV in 620 of 620 entries (largest relative rounding difference 4.1e-4). `rayleigh_spectra.json` =
   `explainer_tables.json → rayleigh` (8 profiles); a full live recomputation rounds to the file in every c_r and c_i entry.
   Against `HEAD`: only the keys `tg_map` and `rayleigh` of `explainer_tables.json` changed — in `rayleigh` exactly four
   physical entries (jet k = 0.1: 0.09273 + 0.2149i → 0.0928 + 0.215i; shear layer k = 0.1: c_i 0.836 → 0.8364) plus three
   round-off-sized c_r of the shear layer (|c_r| < 1e-11); `benard`, `taylor`, `lorenz_sweep`, `note`, `benchmarks.json`,
   `critical_points.json`, `os_modes.json` are identical.
9. **Rayleigh table box.** With max(profile box, `decay_box(k)`) the k = 0.1 rows equal the value in twice and four times
   the box to ≤ 3.5e-10 relative (shear layer c_i 0.8364428; jet 0.0927971 + 0.2149559i); the fixed boxes gave 0.8360159 and
   0.0927330 + 0.2148749i (5.1e-4, 4.4e-4). `scripts/ch11_inviscid_criteria.py`: c_i(k = 0.02) = 0.9649425 in the box 600,
   0.9649424 in 1200, 0.9469873 in the old box 30 (−1.9 %).
10. **`ValueError` paths.** `benard_growth_rate(…, N = 6, 8, 10)` with the filter raises `ValueError` (message names the
    filter; `filter=False` and N ≥ 12 work); `benard_marginal_Ra_det` raises `ValueError` instead of returning NaN.
    **Correction to the claim "K = 0.5 even and odd K ≤ 1.5 now raise":** measured on K = 0.3 … 9.0 (step 0.1) the default
    bracket (Ra_max = 2×10⁴) finds the even root for K ≥ 0.6 and the odd root only for **K ≥ 4.0** — the odd mode raises
    for every K < 4.0 (e.g. K = 2: root 47 006; K = 3: 26 147), and "pass Ra_max ≥ 1e5" is not enough for K ≲ 1.2 (K = 1:
    163 128; K = 0.5: 629 153). With `Ra_max=1e6` the determinant equals the Chebyshev value to ≤ 2e-9 in all those cases.
    Behaviour is correct (a clear error, never NaN); the docstring and the message are imprecise — Should-fix 1.
11. **Double-diffusion script**: the "diffusive" row is now classified diffusive, σ_max = 9.65 + 70.39i (= the tested
    (Ra, Rs) = (−2×10⁴, −1.9×10⁶) root 9.650 ± 70.389i).

### Test changes in loop 2 (verifier)
- **Corrected, not loosened:** `test_tg_tables_V1_csv_equals_explainer_json_and_small_k_columns_are_live` asserted
  `positives == 299`. That number was read off the loop-1 table and so pinned its defect (36 unstable grid points at
  k ≥ 0.65 reported as 0); it contradicted the verifier's own large-k test (a strip ≤ 0.03 needs ≥ 313 positives). It is
  now derived from the grid: the set of positive entries must **equal** the set {J < k(1 − k)} (335 points, counted
  independently by a live recomputation) — a stricter statement than any pinned count.
- **Tightened:** `test_rayleigh_spectra_V1_cached_json_matches_live` additionally recomputes the k = 0.1 and 0.2 rows of the
  two unbounded profiles (file = 4-s.f. rounding of live; live = twice-the-box value to 1e-8; fixed box 3e-4 … 7e-4 off; file
  = explainer key); the slow box test additionally requires box 120 vs 480 at k = 0.1 to agree to 1e-8 (the old 1e-3
  statement about the profiles' fixed boxes is kept as it is).
- **New (9):** `test_decay_map_scale_V1_formula_guards_and_limits`,
  `test_decay_map_scale_V3_tg_growth_unchanged_by_N_doubling_and_scale_halving`,
  `test_tg_growth_V1_shooting_on_the_unbounded_layer_agrees_incl_near_neutral` (independent route: no box, no Chebyshev,
  no filter; 10 points, 3e-6 absolute),
  `test_tg_growth_V7_no_growth_at_or_above_the_neutral_curve_off_grid` (70 points),
  `test_tg_growth_V3_near_neutral_sliver_is_as_documented` (5 columns; also that docstrings and CSV header state it),
  `test_tg_growth_V7_default_is_the_map_scale_rule_and_explicit_scale_overrides` (incl. the cache key),
  `test_taylor_goldstein_V1_N2_zero_equals_rayleigh_on_the_rule_grid`,
  `test_parallel_profile_V7_sech2_far_field_is_finite_and_the_clip_is_invisible` (incl. the k = 0.02 sweep box of
  `scripts/ch11_inviscid_criteria.py`), `test_benard_guards_V7_unresolved_or_unbracketed_cases_raise_value_error`.
- No existing test needed a change because of the Bénard behaviour change (the determinant tests use even K ≥ 2 and odd
  K ≥ 4.5). `tests/ch11_verify_figures.py`: the semicircle figure now uses `decay_map_scale(k)` like `tg_growth`.
- Unchanged and passing: the two loop-0 tests, the loop-1 large-k test, all 7 other loop-1 tests.

## Environment
python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · matplotlib 3.11.2 · pint 0.25.3 (Windows, `.venv`).

## Validation table (CORE rows; every computable row has ≥ 2 independent levels incl. one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers | Label | Tests |
|---|---|---|---|---|---|---|
| C01 normal mode (11.1), σ = −i∣K∣c, stability vocabulary | CORE | `ST.normal_mode`, `sigma_from_c`, `c_from_sigma`, `stability_class`, `stability_verdict`, `marginal_type`, `ch11.normal_mode_growth`, `potential_well_demo` | V1 + V7 + V4 | two forms agree to 1e-12 on a (x, y, t) field for 5 random (k, m, c); σ(2, 1 + 0.5i) = 1 − 2i; ch07 `real_field` parity 1e-12; "kh" σ = 1.3229 − 1.5i exact; "interface" |σ_i| = ch07 ω to 1e-12; potential-well energy non-increasing (1e-9) | analytic | `test_normal_mode_V1_…`, `test_stability_vocabulary_V7_…`, `test_normal_mode_growth_V1_…`, `test_potential_well_V4_…` |
| C02 KH (11.18), criterion, Ex. 11.1/11.2 | CORE | `kh_phase_speed`, `kh_discriminant_terms`, `kh_growth_rate`, `kh_critical_k`, `kh_stability_boundary`, `kh_min_shear`, `kh_unstable_band`, `vortex_sheet_c`, `rayleigh_taylor_cutoff`, `kh_amplitudes`, `kh_fields`, `kh_residuals`, `kh_mixing_energy` | V1 + V2 + V7 + V4 | both roots solve ρ₁(U₁ − c)² + ρ₂(U₂ − c)² = (g/k)(ρ₂ − ρ₁) to 1e-12 on 400 random k × 4 density pairs; (11.19) = ch07 (7.95) to 1e-13; (11.20); Galilean shift 1e-12; c_r sign flips with the moving stream, c_i does not; c_i = 0 at k_c and > 0 at k_c(1 + 1e-6); boundary ΔU_min(k) has zero discriminant to 1e-11 (h = ∞, 1 cm, 10 cm); min shear 6.703 m/s at λ = 1.727 cm (= numerical minimum to 1e-10); band [149.15, 887.44] m⁻¹ at 8 m/s; RT cut-off 1.727 cm (neutral there); residuals of (11.3), (11.9), (11.13) ≤ 1.1e-15; mixing E_f/E_i = 2/3 (1e-14), momentum conserved, 10 random erf smoothings lower E | analytic, symbolic | `test_kh_*` (10 tests incl. D02, D03, D04) |
| C03 Bénard amplitude problem (11.36)–(11.37), Ra (11.21), σ real | CORE | `rayleigh_number`, `gamma_conventions`, `benard_base_state`, `benard_scales`, `benard_growth_rate`, `benard_free_free_sigma` | V1 + V2 + V3 + V7 | Chebyshev σ = D12 quadratic to 1e-8 (4 (K, Ra, Pr)); 11.0155, −40.6243 at K = π/√2, Ra = 2000, Pr = 1; all |σ_i| < 1e-10·|σ| over 3 cases; σ₁ = 0 (≤ 7e-9) at the marginal Ra, sign change at ±1 %; Pr leaves the margin unchanged (Pr = 0.01, 1, 100); convergence 1.2e-6 (N = 12) → 1e-12 floor (N ≥ 16); Ra pint-dimensionless and invariant under (mm, min, °F-sized) unit changes to 1e-13; heated from above → Ra < 0; (11.23) residual 1e-7 rel | converged, symbolic | `test_benard_growth_rate_*`, `test_rayleigh_number_V2_…`, `test_benard_base_state_V1_…`, D06, D08 |
| C04 rigid–rigid neutral curve, Ra_c | CORE | `benard_marginal_Ra`, `benard_critical`, `benard_neutral_curve`, `benard_char_roots`, `benard_determinant`, `benard_marginal_Ra_det`, `benard_eigenfunction`, `planform`, `benard_neutral_table` | V5 + V1 (two routes) + V3 + V2 | Ra_c = 1707.7618 at K_c = 3.11632 (Chandrasekhar 1707.762/3.117: 1.0e-7, 2.2e-4 — the source rounds 3.1163); rigid–free 1100.6496/2.68234 (1100.65/2.682: 3e-7, 1.3e-4); odd 17610.394/5.36465 (17610.39/5.365); Chebyshev = determinant at K = 2, 3.1163, 5, 7.5 (≤ 1e-9) and odd K = 4.5, 5.3647, 7 (≤ 1e-8); convergence 3e-6 (N = 10) → 2e-10 (N = 14) → 1e-12; det purely imaginary (|Re| < 1e-9|det|), slip #1 breaks it and moves the root 1.6 % (K = 2), 2.4 % (K = 5), 0.44 % (K_c); eigenfunction BCs < 1e-10, parity to 5e-9, w = ∂ψ/∂x, u = −∂ψ/∂z; ∇_H²f = −K²f to 1e-4 (FD) | benchmark, converged | `test_benard_critical_V5_…`, `test_benard_marginal_V1_…`, `…_V3_…`, D10, `test_benard_eigenfunction_V1_…`, `test_planform_V1_…`, `test_benard_tables_V1_…` |
| C05 free–free (11.44), 27π⁴/4 | CORE | `benard_free_free_Ra`, `benard_free_free_critical`, `benard_free_free_mode`, `benard_free_free_sympy` | V1 + V2 | Ra(1, 2, 3, 4) = 1284.23, 667.01, 746.53, 1082.06; Chebyshev free walls = (11.44) to ≤ 1.6e-9 (4 K); sympy: sin nπ(z + ½) kills every even derivative at ±½ (n = 1, 2, 3), solves (11.40), root K² = π²/2, Ra_c = 27π⁴/4 exactly; slip #2 (W(±½) = ∓1) and slip #3 (no root) fail | analytic, symbolic | `test_free_free_*` (D11, D12) |
| C06 salt fingers (11.46) | CORE | `linear_eos`, `thermal_rayleigh_signed`, `salinity_rayleigh`, `salinity_rayleigh_prime`, `double_diffusive_margin`, `salt_finger_unstable`, `double_diffusive_sigma`, `salt_finger_regime` | V1 + V2 + V7 | lhs ≡ Rs − Ra (1e-13); thinnest finger-unstable layer 1.6095 cm (flip at ±0.1 %); κ_s = κ, no salt → (11.21) Ra = 27π⁴/4 (1e-13); τ = 1, Rs = 0 → Bénard quadratic + σ = −a² (1e-9); every root solves the cubic (1e-9); fingers +0.0330081 (real), diffusive 9.650 ± 70.389i; sympy: the 3 × 3 dimensional determinant is −(Pr κ³/d⁸) × our cubic; four regimes classified | analytic, symbolic | `test_salt_finger_*`, `test_double_diffusive_sigma_V1_…`, D13 |
| C07 Taylor (11.51)–(11.54) | CORE | `taylor_marginal_Ta`, `taylor_critical`, `taylor_critical_approx`, `taylor_critical_table`, `taylor_growth_rate`, `taylor_eigenfunction`, `taylor_galerkin_Ta`, `taylor_number`, `taylor_number_narrow_inner`, `ring_interchange_energy`, `rayleigh_circulation_criterion`, `couette_rayleigh_line`, `taylor_stability_boundary` | V1 (two routes, Bénard limit) + V3 + V2 + V7 | μ = 1: Ta_c = Ra_c to 1.5e-9 and Ta(k, 1) = Bénard determinant to 2.5e-8; Galerkin (Ex. 11.9) vs Chebyshev ≤ 1.1e-4 (M = 4), 1.4e-6 (M = 8); μ = 0: 3389.900 at 3.12660; μ = 0.5: 2275.095 at 3.11749; μ = −0.5: 6413.722 at 3.19845; (11.54) +0.77 % (μ = 0), +0.10 % (μ = 0.5); convergence 9.5e-6 (N = 10) → 1.5e-9 (N = 14), round-off 6.5e-9 at N = 60; Ta = 6097.56 (narrow 6250, |ratio − 1| = d/2R); Ta = −4AΩ₁d⁴/ν² with ch08 `circular_couette` A (1e-12); Ta invariant under unit changes; Ta < 0 beyond Rayleigh's line; ΔE = E_f − E_i | converged, symbolic | `test_taylor_*`, `test_ring_interchange_…`, D14 |
| C08 Taylor–Goldstein (11.61) | CORE | `ST.taylor_goldstein_eigs`, `ST.decay_box`, `ST.decay_map_scale`, `richardson_profiles`, `tg_growth`, `tg_growth_map`, `tg_tanh_neutral_J`, `tg_tanh_neutral_mode`, `stratified_shear_sympy` | V1 (three routes) + V5 + V2 + V7 + V3 | N² = 0 equals `rayleigh_eigs` to 3e-12 observed, 1e-8 asserted (k = 0.2, 0.4449, 0.7 in `decay_box(k)`; and on the rule grid at k = 0.05, 0.2, 0.7, 0.9, 0.95: ≤ 5.1e-11); exact neutral mode residual 1e-6 (FD) and 40-digit exact; **shooting on the unbounded layer (no box, no Chebyshev, no filter) = `tg_growth` to ≤ 1.7e-6 absolute at 35 points (3e-6 asserted at 10), incl. the 20 points under the neutral curve at k ≥ 0.5 that the loop-1 scale reported as 0**; growth inside / zero outside the tongue; conjugate pairs to 1e-6; Michalke k = 0.44492 (7e-4), kc_i = 0.18970 (1e-5); box: kc_i unchanged to 1e-6 absolute between 8, 12, 16, 24, 48 e-folds; map scale: default (N = 100) = N = 200 = half the scale to ≤ 1.7e-6 (3e-6 asserted); whole map N = 100 vs 200: same zero pattern, ≤ 2.05e-6 absolute, ≤ 2.2e-4 relative; table: positives = {J < k(1 − k)} exactly (335 of 620), CSV = JSON `tg_map` = live (620 / 620 after 4-s.f. rounding); no growth at 504 + 56 points on or above the curve / above Ri = ¼; `decay_box`, `decay_map_scale` formulas and guards; stated limitation measured: 0 returned only within 0.0002 (k = 0.05) … 0.0058 (k = 0.95) of the curve, where kc_i < 0.004 | **converged** (a 0 within 0.006 of the neutral curve means "kc_i < 0.004", not "stable" — documented in the docstrings and the CSV header) | `test_taylor_goldstein_*`, `test_tg_growth_*` (8), `test_decay_box_*`, `test_decay_map_scale_*`, `test_tg_tables_V1_…`, `test_richardson_profiles_V7_…`, `test_decay_callers_V7_…`, D17 |
| C09 Miles–Howard (11.67) | CORE | `gradient_richardson`, `miles_howard_stable`, `richardson_identity_check`, `stratified_energy_budget` | V4 + V2 + V7 + V1 | (11.65) residual < 1e-8 for 3 computed modes, Im parts equal (1e-7) and negative; J = 0.26, 0.3, 0.5: no growth (5 k each); 6 random Ri > ¼ walled profiles: no converged unstable mode; Ri = J cosh²z (1e-12 with U′, 1e-6 central difference); Ex. 11.10 budget residual < 1e-9 | conserved, symbolic | `test_taylor_goldstein_V4_…`, `…_V7_…`, `test_richardson_V1_…`, D18 |
| C10 Howard's semicircle (11.72) | CORE | `ST.howard_semicircle`, `in_howard_semicircle`, `howard_identity_check` | V1 + V7 + V4 + V2 | boundary points on the circle; (11.69), (11.70) residuals < 1e-8, c_r = ∫UQ/∫Q to 1e-8, ∫(U² − c_r² − c_i²)Q ≥ 0; every unstable c of tanh, Bickley, sin (Rayleigh) and tanh/J (TG) inside, kc_i ≤ kΔU/2 | analytic, conserved | `test_howard_semicircle_*`, D19 |
| C11 Orr–Sommerfeld (11.79), Squire (11.78) | CORE | `ST.orr_sommerfeld_eigs`, `os_mode`, `squire_transform`, `ch11.os_3d_eigs`, `os_derivation_sympy`, `ST.cheb*`, `apply_bc_rows`, `generalized_eigs`, `constrained_eig`, `converged_eigs` | V5 + V3 + V2 + V7 + V1 (two BC routes) | Orszag c = 0.23752649 + 0.00373967i: |Δc| < 1e-8 (N = 100); errors 9e-6 (N = 40), 2e-7 (50), 1e-8 (60), 1e-10 (80) — exponential; row-replacement route = elimination route to 1e-9; Re → ∞: |c − c_Rayleigh| 3.7e-5 at Re = 1e5, monotone; viscous spectrum not closed under conjugation; 3-D primitive spectrum contains the OS mode at (k̄, Re̅) to 1e-8 (2 cases); Rebar = Re cos 30° exactly | benchmark, converged | `test_orr_sommerfeld_*`, `test_os_mode_V1_…`, `test_squire_V1_…`, D21 |
| C12 Rayleigh's inflection-point theorem, Fjørtoft | CORE | `ST.rayleigh_eigs`, `inflection_points`, `rayleigh_criterion`, `fjortoft_criterion`, `rayleigh_identity_check`, `fig_11_21_verdicts`, `sin_profile_max_growth`, `critical_layer`, `cats_eye_*`, `piecewise_*`, `rayleigh_spectrum_table`, `parallel_profile`/`inviscid_profile` | V1 + V4 + V2 + V7 + V5 | analytic neutral modes (tanh k = 1, c = 0; Bickley k = 2, 1, c = 2/3) exact at 40 digits, computed branches end there (c_r within 0.007 of 2/3); (11.83), (11.84) residuals < 1e-9; no unstable mode for Poiseuille, Couette, Blasius; Fig. 11.21 verdicts (a–c none, d Rayleigh only, e–f both); sin y stable at 2b = 3, growth 0.0114, 0.0596, 0.157 at b = 1.7, 2, 3; cat's-eye expansion error order 3.00; piecewise layer kh = 1.278465, max 0.2012 at 0.797, Ex. 11.11(d) | analytic, conserved | `test_rayleigh_*`, `test_inflection_and_fjortoft_V1_…`, `test_sin_profile_V7_…`, `test_critical_layer_…`, `test_piecewise_…`, D22 |
| C13 plane Poiseuille, Table 11.1 | CORE | `poiseuille_critical`, `poiseuille_neutral_curve`, `poiseuille_spectrum`, `couette_max_growth`, `blasius_*`, `bickley_critical`, `tanh_*`, `falkner_skan_neutral_curve`, `table_11_1`, `tollmien_profile`, `ST.critical_point`, `ST.neutral_curve` | V5 + V3 + V1 (two routes) + V7 | **live** (cache=False) Re_c = 5772.2223 (N = 60), 5772.2218 (N = 80), k_c = 1.020547, 1.020545, c_r = 0.264000 (Orszag 5772.22, 1.02056 ± 1e-5, 0.26400: 4e-8, 1.5e-5, 2e-7); Re = 10⁴ band 0.797–1.0947 (bracketed by the row-replacement route); Couette c_i < 0 on 4 k × 3 Re; Blasius live N = 60 519.005, cached 519.077, k 0.3038, ω 0.1205 (Thomas 519.2/0.303/0.120: 0.02 %, 0.26 %, 0.41 %); Bickley 4.0170 at 0.1728 (Tatsumi & Kakutani 4.0 at 0.2: 0.42 %, k 14 % — a half-domain independent solve confirms Re_n(0.1728) = 4.01705 < Re_n(0.2) = 4.0366); tanh k_u(2, 10, 200) = 0.2732, 0.6544, 0.9702 → 1 | benchmark, converged | `test_poiseuille_*`, `test_os_tables_V1_…`, `test_couette_V7_…`, `test_blasius_critical_…`, `test_bickley_*`, `test_tollmien_profile_V1_…` |
| C14 disturbance-energy equation (11.88) | CORE | `ST.disturbance_energy_budget`, `ts_wave_fields`, `ts_mode`, `energy_equation_sympy` | V4 + V2 + V7 | dE/dt − (P − Λ) < 1e-8·max(P, Λ) for 5 OS modes (Poiseuille × 3, tanh, Blasius; observed 4.5e-12 at Re = 10⁴); sign(P − Λ) = sign(c_i) for those and every |c_i| > 1e-4 entry of the cached 600-point grid; P/Λ = 1.6163 (Re 10⁴), 0.7654 (Re 5000); ⟨uv⟩ = ½Re(ûv̂*) (sympy) and = numerical x-average (1e-14) | conserved, symbolic | `test_energy_budget_V4_…`, D23 |
| C15 Lorenz (11.91), chaos | CORE | `lorenz_*`, `pendulum_*`, `phase_portrait`, `hopf_normal_form`, `limit_cycle_amplitude`, `logistic_*`, `cobweb`, `bifurcation_diagram`, `superstable_points`, `period_doubling_points`, `feigenbaum_estimate` | V1 + V5 + V4 + V3 + V2 | fixed points exact; eig(C) = 0.09396 ± 10.19451i, −13.85458; origin 11.8277, −2.6667, −22.8277; r_H = 24.73684 = root of Re λ(r) (1e-8) (Wikipedia 24.74); tr J = −(Pr + 1 + b) everywhere; RK4 order 4.00 ± 0.15; RK4 = DOP853 at t = 5 (1e-6); separation slope 0.897, predictability 29.14; Lyapunov (t = 60) 0.7–1.2; r = 0.5 → origin, r = 15 → |X| = √(b(r − 1)); r(Ra_c(k), k) = 1; pendulum energy ptp < 1e-9; Hopf radius √μ (1e-6); logistic A₂ = 1 + √6, A₃, A₄, A₅ to 2e-7, δ₈ = 4.669191 | analytic, benchmark | `test_lorenz_*` (D24, D25), `test_pendulum_hopf_…`, `test_logistic_…` |

NOTE items coded and tested (≥ 1 level each): N09/N12/N14 (KH residuals), N15–N19, N22 (mixing), N25/N28/N31, N40, N44, N46,
N48–N50, N51, N53–N56, N59, N62–N68, N73, N75–N82, N88–N97, N99, N100, N102, N104, N105, N107, N108, N113, N114, N116, N117,
N120–N127, N129 — see the test names above. RECAP parities: R09 = ch07 (7.95), R16 = ch08 `circular_couette`, R18–R20
(Ri), R23 (Bickley), R24 (Falkner–Skan: slow test, adverse m = −0.05 unstable at Re 300 where m = +0.05 is not).

## Derivations (curation §4b) — every ★★/★★★ D row re-derived with sympy, independently of the chapter's own engines
| D | Result | ★ | Test | What the sympy re-derivation checks (Part F lines) | Result |
|---|---|---|---|---|---|
| D02 | (11.9), (11.13) | ★★ | `test_kh_V2_derivation_D02_…` | step 5 (the √(1 + ζ_x²) cancels from all three members, exact), step 7 (O(ε) Taylor transfer), step 9 (O(1) balance (11.12)), step 12 (U∂φ/∂x survives); second-order remainder non-zero | all lines ✓ |
| D03 | (11.18), criterion | ★★ | `…_D03_…` | dsolve of A″ = k²A, steps 5–6 amplitudes, step 8 line, step 9, step 12 identity, both roots, step 13 inequality ×k(ρ₁ + ρ₂)²; coded root = derived root at 20 random points (1e-12) | ✓ |
| D04 | Ex. 11.1, ΔU_min | ★★ | `…_D04_…` | cosh k(z + h) meets the bottom, step 3 A₊, step 6 line (with σ_s k²), coded root solves step 6 (15 random cases), step 8–9 minimum, coth → 1 | ✓ |
| D06 | (11.29) | ★★ | `test_benard_V2_derivation_D06_…` | built on a divergence-free vector-potential field: div of (11.26) is the Poisson part only, ∇²R_z − ∂_z(∇·R) ≡ (11.29); planted ∇² ≠ | ✓ |
| D08 | σ real | ★★ | `…_D08_…` | steps 2, 6, 7–8 by parts on complex polynomial trial functions; steps 10–12; **and step 10 on three computed Chebyshev eigenmodes (1e-6)** | ✓ |
| D10 | roots, det, Ra_c | ★★★ | `…_D10_…` | step 1–4 (each root of (11.42)), step 9 operator identities, step 12 (det imaginary, 4 (Ra, K); slip #1 not), step 13 sign change at the Chebyshev root, slip #1 shifts | ✓ |
| D11 | (11.43)–(11.44), π²/2 | ★★ | `test_free_free_V2_derivation_D11_…` | even derivatives to order 6 vanish (n = 1–3), (11.40) residual, step 9 quotient rule, roots, 27π⁴/4, slip #3 root set empty | ✓ |
| D12 | free–free growth | ★★ | `…_D12_…` | steps 3, 4, 6, 7 (discriminant = a⁴(1 − 1/Pr)² + 4RaK²/(Pr a²)), 8; 10 random numeric root comparisons (1e-10) | ✓ |
| D13 | (11.45)–(11.46) and our cubic | ★★ | `test_double_diffusion_V2_derivation_D13_…` | dimensional 3 × 3 determinant (heat, salt, w) with steps 4–5 scales = −(Pr κ³/d⁸) × our cubic; σ = 0 ⇒ Rs − Ra = a⁶/K²; 27π⁴/4; root 0.0330081 (and 0.0307930 if τ = 0.01) | ✓ (see doc correction 1) |
| D14 | (11.50)–(11.52) | ★★★ | `test_taylor_V2_derivation_D14_…` | steps 2–4 linearisation (factor 2), step 6 operators, step 7 û_z, step 8 p̂ satisfies the axial equation exactly, step 10 (k² × radial = the printed line), step 9 2A, steps 13–15 (narrow gap, rescale, Ta = −4AΩ₁d⁴/ν²), (11.52) identity, slip #12 by L-homogeneity | ✓ |
| D17 | (11.61) | ★★ | `test_tg_V2_derivation_D17_…` | steps 2–4 ((11.58)–(11.60) from (11.57)), steps 5–9 elimination; slip #4 (∂p/∂x) does not reach (11.61) | ✓ |
| D18 | Miles–Howard | ★★★ | `test_miles_howard_V2_derivation_D18` | steps 3, 4 (both derivative formulas), step 6 (the "rearrangement" line), step 8 (11.64), step 12 Im parts; step 13 sign on a computed mode | ✓ |
| D19 | Howard | ★★★ | `test_howard_V2_derivation_D19` | step 2, steps 3–4 (U″ cancels), step 5 divergence form; slip #11 form ≠; step 8; steps 9–11 moment algebra; step 14 completing the square | ✓ |
| D21 | (11.79) | ★★ | `test_os_V2_derivation_D21_…` | steps 1–8 from (11.77) (continuity automatic, p̂′ cancels, step 6 line); planted v̂ = +ikφ ≠ | ✓ |
| D22 | (11.83)–(11.84) | ★★ | `test_rayleigh_V2_derivation_D22_…` | step 4 by parts (complex polynomial), step 6 conjugate trick; steps 7–9 on a computed tanh mode (positive and negative U″ parts cancel to 1e-9) | ✓ |
| D23 | (11.88) | ★★ | `test_energy_V2_derivation_D23_…` | steps 1–7 pointwise identity in index form (2-D, ∂_ju_j = 0 by ψ), step 12, ⟨uv⟩ = ½Re(ûv̂*), step 10 periodic integral | ✓ |
| D24 | (11.91) | ★★★ | `test_lorenz_V2_derivation_D24_…` | step 4 (J(ψ, ∇²ψ) = 0), step 8 Jacobian line, steps 7, 9, 10 Galerkin projections, steps 11–15 scalings (α_L, β_L, γ_L, r) → (11.91), b(π/√2) = 8/3 | ✓ |
| D25 | fixed points, r_H | ★★ | `…_D25_…` | steps 1–4 solve, step 5 Jacobian, step 6 factorised origin polynomial, step 9 cubic, steps 10–11 r_H, the λ = iω condition a₂a₁ = a₀ | ✓ |

No intermediate line of Part F was found wrong. The chapter's own engines (`kh_sympy`, `kh_depth_tension_sympy`,
`benard_perturbation_sympy`, `exchange_of_stabilities_sympy`, `benard_free_free_sympy`, `taylor_perturbation_sympy`,
`stratified_shear_sympy`, `os_derivation_sympy`, `energy_equation_sympy`, `lorenz_sympy`, `derive_all`) are tested separately:
all residuals 0, all planted variants non-zero.

## Functions used by the notebook and explainers (design Part C) — each is called by at least one test
`test_part_c_V7_every_contract_function_exists_and_is_exercised` parses design Part C (C.1 + C.2: **149 names**) and asserts
each exists in `ch11` and is called as `ch11.<name>` or `ST.<name>` somewhere in this file (plus 9 helpers the scripts use:
`parallel_profile`, `fig_11_21_verdicts`, `ts_mode`, `tg_tanh_neutral_J`, `tg_tanh_neutral_mode`, `superstable_points`,
`logistic_fixed_point`, `tollmien_coefficients`, `lorenz_divergence`), and that 12 explainer-mirrored functions are
scalar-callable (selftest parity rows). Writers that touch `reference/ch11` are exercised with `write=False` / a temporary
folder (`neutral_curve_tables(tmp_path, fast=True, cache=False)` in the slow test). New since Part C was written:
`ST.decay_box` / `ch11.decay_box` (loop 1) — V1 formula, guards and e-fold count, V3 box and N doubling
(`test_decay_box_V1_…`, `test_decay_box_V3_…`), label `converged`; `ST.decay_map_scale` / `ch11.decay_map_scale` (loop 2)
— V1 formula, limits, guards, effect on the grid (`test_decay_map_scale_V1_…`), V3 N doubling and scale halving
(`test_decay_map_scale_V3_…`), V1 shooting route through `tg_growth`, label `converged`; the new `ValueError` paths of
`benard_growth_rate` and `benard_marginal_Ra_det` (`test_benard_guards_V7_…`). The heavy cached functions are each
compared with a live recomputation at sample points (Bénard table, Taylor table, TG map, Rayleigh spectra, OS grid,
critical_points.json, Poiseuille/Blasius/Bickley critical points, Poiseuille neutral curve, Lorenz sweep).

## Convergence studies
| Quantity | Grids | Observed | Design |
|---|---|---|---|
| Chebyshev −u″ = λu (first 3 λ) | N = 8, 16, 32 | error 2e-6 → < 1e-9 → < 1e-10 | spectral |
| Bénard σ (K = 3, Ra = 2500, Pr = 0.7) | N = 12, 16, 20, 24 vs 30 | 1.2e-6, then 1e-12 floor | spectral |
| Bénard marginal Ra (K = 3.1163) vs determinant | N = 10, 14, 20, 30 | 3.0e-6, 2.1e-10, 3e-12, 5e-12 | spectral |
| Taylor marginal Ta (μ = 0, k = 3.1266) | N = 10, 14, 20 vs 24 | 9.5e-6, 1.5e-9, 3e-12; round-off grows to 6.5e-9 at N = 60 | spectral |
| Orr–Sommerfeld TS mode (Orszag case) | N = 30, 40, 50, 60, 80 vs 140 | 9e-6 (40), 2e-7 (50), 1e-8 (60), 1e-10 (80); log-error slope < −0.2 per N | spectral |
| Poiseuille Re_c | N = 60, 80 (live) | 5772.2223, 5772.2218 (8e-8) | converged |
| Blasius Re_c | N = 60 / 100 (y_max 20), 100 (y_max 40) | 519.005 / 519.077 / 519.056 | converged (1.4e-4, 4e-5) |
| Galerkin (Ex. 11.9) vs Chebyshev | M = 4, 8 | 1.1e-4, 1.4e-6 | converging |
| Lorenz RK4 | dt = 0.02, 0.01, 0.005 (t = 0.5) | order 4.00 ± 0.15 | 4 |
| cat's-eye expansion error | Δy = 0.04, 0.02, 0.01 | order 3.00 ± 0.15 | 3 |
| TG growth vs box size (loop 0, fixed y_max = 30) | y_max = 30, 120, 240 | 26 % at (0.05, 0.04); spurious 0.0078 at (0.05, 0.05) | FAILED in loop 0 |
| TG growth vs box size (loop 1, `decay_box(k)` = max(30, 12/k)) | box, 2 × box, 4 × box; 8, 16, 24, 48 e-folds | (0.05, 0.04): 0.0163047, 0.0163048, 0.0163046; (0.05, 0.045): 0.0082157, 0.0082158, 0.0082156; (0.1, 0.08), (0.1, 0.085), (0.15, 0.1), (0.2, 0.15), (0.3, 0.1): identical to 7 decimals; (0.05, 0.05): 0 in every box ≥ 160 | converged (≤ 2e-7 absolute observed; 1e-6 asserted) |
| TG growth vs N at small k (default box; loop-2 default scale) | N = 100, 200, 240 | ≤ 4e-7 absolute at (0.02, 0.01), (0.03, 0.02), (0.05, 0 … 0.045), (0.1, 0.08), (0.1, 0.085), (0.15, 0.1), (0.2, 0.15), (0.3, 0.1), (0.3, 0.2); (0.15, 0.125), 0.0025 below the neutral curve: 0.0058248 vs 0.0058244 (loop 1, s = 0.5: 0.0058231) | converged |
| TG growth vs N at k ≥ 0.65, near neutral (loop 1, fixed map scale 0.5) | N = 100, 160, 240 (filtered `tg_growth`) | (0.9, 0.04): 0 / 0.03324 / 0.03325; (0.9, 0.06): 0 / 0 / 0.02010; (0.8, 0.13): 0 / 0.02142 / 0.02159; (0.95, 0.01): 0 / 0.02418 / 0.02428; zero strip per column at N = 100: 0.0375 … 0.06 | FAILED in loop 1 |
| TG growth vs N, map scale and box (loop 2, `decay_map_scale(k)` = min(0.5, max(0.035, 0.025/k))) | N = 100 (default), 200, 240; scale × ½, × 2 (N = 200); box × 2 (N = 160); shooting (unbounded) | 35 points: (0.9, 0.04) 0.0332469 / 0.0332455 / 0.0332455, shooting 0.0332455; (0.95, 0.04) 0.0049144 / 0.0049133 / 0.0049133, shooting 0.0049133; (0.7, 0.2) 0.0079147 / 0.0079164 / 0.0079164 (largest difference, 1.7e-6); N = 240 vs shooting ≤ 1.2e-10; whole 20 × 31 map N = 100 vs 200: same zero pattern, ≤ 2.05e-6 absolute, ≤ 2.2e-4 relative; zero strip per column 0.0075 / 0.01 (one J step) | converged (≤ 1.7e-6 absolute observed, 3e-6 asserted) |
| TG near-neutral sliver (where 0 is still returned under the curve) | bisection in J, 13 columns | width 0.00019 (k = 0.05), 0.00110 (0.3), 0.00207 (0.5), 0.00324 (0.7), 0.00536 (0.9), 0.00583 (0.95); growth at its edge 0.0011 … 0.0038; shooting growth at mid-sliver 0.0006 … 0.0019 | as documented (resolution limit of the N-filter, kc_i < 0.004) |
| Rayleigh c vs box size (profile boxes 30 tanh / 40 Bickley against 2 × `decay_box`) | k = 0.02, 0.05, 0.1, 0.2, 0.3 | fixed boxes: tanh 1.9e-2, 5.0e-3, 5.1e-4, 3.3e-6, 1.6e-8; Bickley 2.0e-1, 2.1e-2, 4.4e-4, 1.9e-7, 1e-10 (relative); loop 2 — table and script use max(profile box, `decay_box(k)`): k = 0.1 rule vs 2 × and 4 × the box ≤ 3.5e-10; k = 0.02 (script sweep) box 600 vs 1200: 1.7e-7 | converged (the E7 table's k = 0.1 row is now right to its 4 s.f.) |

## Conservation / invariant residuals
KH interface conditions ≤ 1.1e-15 · (11.65) < 1e-8 (observed 1e-11) · (11.69) < 1e-8, (11.70) < 1e-8 · Ex. 11.10 budget
< 1e-9 relative (observed 2e-14) · (11.83), (11.84) < 1e-9 · (11.88) dE/dt − (P − Λ) < 1e-8 relative (observed −4.5e-12 at
Re = 10⁴) · Lorenz tr J = −(Pr + 1 + b) to 1e-14 · pendulum energy ptp < 1e-9 · potential wells: energy non-increasing (1e-9)
· KH momentum conserved by mixing (exact) · D08 step 10 on computed modes < 1e-6.

## Benchmarks used (`reference/ch11/SOURCES.md`, re-read 2026-09-30)
| Benchmark | Value | Ours | Source | Status |
|---|---|---|---|---|
| plane Poiseuille | Re_c 5772.22, α_c 1.02056 ± 1e-5, c_r 0.26400 | 5772.2218, 1.020545, 0.264000 | Orszag, JFM 50, 689 (1971), PDF read | ✓ (α_c 5e-6 outside Orszag's ±1e-5 band edge: 1.5e-5 relative, within 1 %) |
| Orszag eigenvalue | 0.23752649 + 0.00373967i | 0.2375264888 + 0.0037396706i | same | ✓ 1e-9 |
| rigid–rigid / rigid–free / free–free Bénard | 1707.76/3.117, 1100.65/2.682, 657.511/2.2214 | 1707.7618/3.11632, 1100.6496/2.68234, 657.5114/2.22144 | Chandrasekhar (1961) via Nek5000 examples PDF table, arXiv:nlin/0302057 | ✓ |
| odd Bénard mode | 17610.39 at 5.365 | 17610.394 at 5.36465 | Chandrasekhar p. 39 (not reachable today) | two independent routes + V6 |
| Blasius (parallel) | 519.2, 0.303, 0.120 | 519.077, 0.30378, 0.12049 | Thomas via Gallagher, Griffiths & Stephen, Phys. Fluids 28, 074107 (2016), Table II | ✓ ≤ 0.41 % |
| Bickley jet | Re_c ≈ 4.0 at k ≈ 0.2 | 4.0170 at 0.1728 | Tatsumi & Kakutani (1958) via a 2017 abstract | Re ✓ 0.42 %; k 14 % (approximate 1958 value; independent half-domain solve agrees with ours) |
| tanh shear layer | k 0.4446, kc_i 0.1897 | 0.44492, 0.18970 | Michalke (1964) (digits secondary) | approximate ✓ 7e-4 |
| Lorenz Hopf | 24.74 | 24.736842 | Wikipedia "Lorenz system" | ✓ |
| Feigenbaum | δ 4.669201609; a_n 3.4494897, 3.5440903, 3.5644073, 3.5687594 | δ₈ 4.669191; a_n to 2e-7 | Wikipedia "Feigenbaum constants" | ✓ |

## Numbers from the text (book vs ours) — private values redacted to relative differences
Rigid–rigid Ra_c 0.014 %, K_c 0.12 %; cell width λ/d 0.81 % (the book gives one significant figure); free–free 0.002 %;
rigid–free Ra_c 0.032 %, K_c 0.088 %; odd mode (Ex. 11.7) Ra 0.002 %, K 0.007 %; Taylor numerator of (11.54) 0.014 %, k_c 0.12 %;
Ex. 11.9 k_c (μ = 0) 0.21 %; double-diffusive threshold 0.08 %; mixing energies exact; Table 11.1: jet 0.42 %, Blasius
0.18 %, plane Poiseuille 0.13 % (the book rounds 5772.22 up), shear layer exact (0), pipe/Couette "stable" ✓; Tollmien profile
with the book's coefficients: corrected middle branch continuous to 0.49 % of the inner value (book rounding), the printed
branch jumps by 99.8 %; tanh most-amplified wavelength / thickness 0.87 % (book "≈", one figure); Lorenz parameters exact;
Feigenbaum δ 2e-4 %. All within 0.5 % except the two one-significant-figure values (explained).

## Figures reproduced with our code (`outputs/ch11/verify/`, local; `tests/ch11_verify_figures.py`)
| Figure | PNG | Visual verdict |
|---|---|---|
| Fig. 11.10 | `fig_11_10_benard_neutral_curves.png` | Four U-shaped curves ordered rigid > rigid–free > free–free with minima at 1707.8/3.12, 1100.6/2.68, 657.5/2.22 and the odd branch an order higher (17610 at 5.36); all rise ∝ K⁴ at large K. |
| KH (implied by (11.18), Ex. 11.1) | `kh_growth_and_boundary.png` | Growth switches on at k_c ≈ 128 m⁻¹ without surface tension and becomes a closed band [149, 887] m⁻¹ with it; the boundary ΔU_min(k) is a U with minimum 6.70 m/s at λ = 1.73 cm; a 5 mm lower layer changes nothing visible (ρ₂/ρ₁ ≫ tanh kh). |
| Taylor vs (11.54) | `taylor_critical_vs_11_54.png` | Exact Ta_c falls from 6414 (μ = −0.5) to 1707.8 (μ = 1); (11.54) lies above it, coinciding for co-rotation and diverging for counter-rotation. |
| TG growth map (E6/F5) | `tg_growth_map.png` | (regenerated in loop 2 and looked at) A smooth symmetric tongue that fills the exact neutral parabola J = k(1 − k) on both sides — the lowest contour hugs the cyan curve from k = 0.05 to 0.95, the ragged right edge of loop 1 is gone; nothing above Ri = ¼; maximum 0.19 at k ≈ 0.45, J = 0; growth → 0 at k = 1. |
| Fig. 11.20 | `fig_11_20_howard_semicircle.png` | (regenerated in loop 2 with `decay_box(k)` and `decay_map_scale(k)`) All 37 TG eigenvalues (c_r = 0, 0.009 < c_i < 0.84; exactly one wherever J < k(1 − k)) inside the unit semicircle and all Bickley Rayleigh eigenvalues inside the [0, 1] semicircle. |
| OS spectrum + thumb | `os_poiseuille_spectrum_and_neutral_curve.png` | Y-shaped spectrum with one unstable mode 0.237526 + 0.003740i; thumb-shaped neutral curve with its nose at (5772.22, 1.0205), upper branch peaking ≈ 1.10 near Re 9000, both branches falling as Re → 10⁶. |
| Fig. 11.23, 11.26 | `fig_11_23_tanh_and_fig_11_26_blasius.png` | tanh upper neutral k rises from 0.27 at Re = 2 toward 1 (Re_c = 0); Blasius neutral curve in F closes at Re_δ* ≈ 519 with F ≈ 2.3e-4 and narrows as Re grows. |
| Fig. 11.25 / p. 520 | `fig_11_25_energy_budget_profiles.png` | Production −⟨uv⟩U′ peaks at the critical layers y ≈ ±0.85; dissipation is concentrated in the wall layers; P/Λ = 1.6163. |
| Figs. 11.29–11.31 | `fig_11_29_31_chaos.png` | Butterfly attractor around C± = (±8.49, ±8.49, 27); separation grows exponentially (slope 0.897) from t ≈ 12 and saturates near t ≈ 32; logistic tree doubles at 3, 3.449, 3.544, 3.564 and becomes chaotic near 3.57 with the period-3 window. |
| Fig. 11.21 | `fig_11_21_profiles_and_criteria.png` | Six stand-ins; inflection markers only on (d)–(f); verdicts (a–c) none, (d) Rayleigh only, (e), (f) both. |

## Deviations & justifications (every `# DEVIATION` / documented choice in the code)
| Deviation | Where | Evidence it is harmless |
|---|---|---|
| Chebyshev collocation (method choice) | `ST.cheb` | V1 exactness, spectral convergence (tables above) |
| BC elimination instead of row replacement | `ST.constrained_eig` and every solver | the row-replacement route (`apply_bc_rows` + `generalized_eigs`) gives Orszag's c to 2e-8 and equals the elimination route to 1e-9 |
| ∞ truncated at y_max | `cheb_grid`, Rayleigh/OS/TG decay and semi-infinite | passes for OS tanh/Bickley/Blasius (slow box test); `tg_growth` / `tg_growth_map` take `decay_box(k)` = max(30, 12/k) and are box-converged to 2e-7 (loop-0 defect fixed); since loop 2 `rayleigh_spectrum_table` and the k-sweep of `scripts/ch11_inviscid_criteria.py` take max(profile box, `decay_box(k)`) (k = 0.1 rows converged to 3.5e-10); `taylor_goldstein_eigs` and `rayleigh_eigs` themselves still default to 30 — every script call either passes `y_max` or uses a literal k ≥ 0.2 (3.3e-6 at k = 0.2) |
| cosh argument clipped at ±350 | `richardson_profiles` (N², U″); since loop 2 `_sech2` (`parallel_profile` shear layer / jet) | closed forms reproduced to 1e-14 (1e-12 for the jet's U″) on |z| ≤ 20; finite and < 1e-300 at |z| = 356 … 1000 with overflow raised as an error (`test_richardson_profiles_V7_…`, `test_parallel_profile_V7_…`); `Up` and `Ri` of `richardson_profiles` are not clipped (they warn beyond |z| ≈ 355, where no solver evaluates them) |
| linear map for semi-infinite OS | `orr_sommerfeld_eigs(bc="semi_infinite")` | Blasius Thomas benchmark 0.02 %; needs N ≈ 100 at y_max = 40 (N = 60 gives 522.0) |
| tan map (s = 1) for Rayleigh decay | `rayleigh_eigs` | Michalke, analytic neutral modes, TG at N² = 0 equal to 1e-8 |
| TG map scale: `decay_map_scale(k)` = min(0.5, max(0.035, 0.025/k)) in `tg_growth` / `tg_growth_map` (loop 2; loops 0–1: fixed 0.5; Part C: 3.0); `taylor_goldstein_eigs` itself still defaults to 0.5 | `decay_map_scale`, `tg_growth` | N = 100 with the rule = N = 200 / 240 = half / twice the scale = shooting to ≤ 1.7e-6; N² = 0 equals Rayleigh on the same grid to 5e-11; no spurious mode at 560 points on or above the neutral curve; the fixed 0.5 stays reachable (`map_scale=0.5`) and reproduces the loop-1 zeros. Side effect to know: the rule grid is coarser in the far field, so the integral identities (11.65), (11.69) of a mode computed on it close only to ≈ 1e-7 … 2e-4 (against < 1e-8 on the s = 0.5 grid the V4 tests use) — the eigenvalue is unaffected (Should-fix 3) |
| N-convergence filter, ci_min | solvers | TS mode kept, spurious ones dropped (`test_converged_eigs_V3_…`); near-neutral TG modes that fail the filter are reported as 0: with the loop-2 map scale only within 0.0002 (k = 0.05) … 0.0058 (k = 0.95) of the neutral curve, kc_i < 0.004 (measured; documented in `tg_growth`, `tg_growth_map`, `taylor_goldstein_eigs` and the CSV header; loop 1: 0.0375–0.06 at k ≥ 0.65) |
| Tollmien breakpoint ours (η₁ = 0.2), coefficients from continuity | `tollmien_profile` | V1 continuity of value and slope; book coefficients (V6) continuous to 0.49 % |
| `os_3d_eigs` takes U′ (Part C wrote Upp) | `os_3d_eigs` | passing U″ changes c from 0.2536 − 0.0011i to 0.794 + 0.411i (wrong physics) — the code is right, Part C is the slip |

## Open items (verbatim for the orchestrator)
Nothing blocks. No row is `qualitative` or `unverified`.

Closed:
- **Closed in loop 1 — `ch11.tg_growth` decay box.** Loop 0: with `y_max = 30`, kc_i at k = 0.05 was 26 % too high at
  J = 0.04 and spuriously positive (0.0078) at J = 0.05 > k(1 − k). Loop 1: `decay_box(k)` = max(30, 12/k).
- **Closed in loop 2 — false zeros under the neutral curve at k ≥ 0.65** (loop-1 Must-fix). `decay_map_scale(k)` =
  min(0.5, max(0.035, 0.025/k)); table regenerated; zero strip 0.0375–0.06 → 0.0075–0.01 (one J step); 299 → 335
  positives = every grid point under J = k(1 − k); confirmed by N doubling, scale and box variation and an independent
  shooting solution (≤ 1.7e-6 absolute).
- **Closed in loop 2 — Rayleigh table box** (loop-1 item 2): `rayleigh_spectrum_table` and the k-sweep of
  `scripts/ch11_inviscid_criteria.py` use max(profile box, `decay_box(k)`); k = 0.1 rows now right to 4 s.f.; `_sech2` clipped.
- **Closed in loop 2 — Bénard robustness** (loop-1 item 3): `ValueError` instead of `IndexError` / NaN.
- **Closed in loop 2 — "diffusive" demo row** (loop-1 item 4): now σ_max = 9.65 + 70.39i, classified diffusive.
- **Closed in loop 2 — stale "Validation (planned …)" module docstring** (loop-1 item 5): rewritten (one sentence of the
  new text is already out of date, see Should-fix 2).

Should-fix (implementer / notebook-builder; documentation and cosmetics, none changes a number):
1. `benard_marginal_Ra_det`: the docstring and the error message say "the odd mode at K ≲ 1.5 has its root above the
   default: pass `Ra_max` ≥ 1e5". Measured (K = 0.3 … 9.0, step 0.1): with the default `Ra_max = 2e4` the even root is found
   for K ≥ 0.6 and the odd root only for **K ≥ 4.0**; the odd mode raises for every K < 4.0 (K = 3: root 26 147; K = 2:
   47 006) and needs `Ra_max` > 1.6×10⁵ at K = 1 (163 128) and > 6.3×10⁵ at K = 0.5 (629 153). The even-mode failure
   (K ≤ 0.5) prints the hint about the odd mode. Either widen the scan automatically (e.g. to 10 × the Chebyshev value) or
   state "even: K ≥ 0.6; odd: K ≥ 4.0 with the default; otherwise pass `Ra_max` above the root". The notebook must not call
   `benard_marginal_Ra_det(K, mode="odd")` with K < 4 without `Ra_max` (tested: `test_benard_guards_V7_…`).
2. Stale text after loop 2: (a) `fluidpy/ch11_instability.py`, `tg_growth`, the comment "s = 0.025/k within [0.05, 0.5]" —
   the floor is 0.035; (b) `fluidpy/core/stability.py` module docstring, "`decay_map_scale` — … (measured by the
   implementer in verify loop 2; the verifier's large-k test exercises it through `tg_growth`)" — it now has its own tests
   (`test_decay_map_scale_V1_…`, `test_decay_map_scale_V3_…`, the shooting route `test_tg_growth_V1_shooting_…`).
3. `taylor_goldstein_eigs` itself keeps the defaults s = 0.5, box 30 (documented). Direct callers that want near-neutral
   modes must pass `map_scale=decay_map_scale(k)`, `y_max=decay_box(k)` (as `tg_growth` does):
   `scripts/ch11_stratified_shear.py` passes the box but not the scale, so its semicircle panel (N = 80, k up to 0.8) can
   miss a near-neutral point (cosmetic; nothing wrong is plotted). Conversely, for integral identities / energy budgets of a
   mode ((11.65), (11.69)–(11.70), Ex. 11.10) keep the s = 0.5 grid away from the neutral curve: on the rule grid the
   quadrature closes only to 1e-7 … 2e-4 (e.g. (11.69) at (0.95, 0.04): 2.2e-4; (11.65) at (0.4, 0.1): 1.2e-6 against
   < 1e-8 on the s = 0.5 grid).
4. The published table is rounded to 4 s.f. while the N = 100 values are converged to ≤ 2.05e-6 absolute: 16 of the 620
   entries sit on a rounding boundary and would print one unit different in the last digit at N = 200 (e.g. (0.7, 0): 0.146
   vs 0.1461). Explainer E6 parity rows that compare against `tg_growth` must use a tolerance ≥ 6e-4 relative (the tests use
   6e-4) — not exact equality.
5. Explainer E6 / the notebook text must carry the stated limitation in words: an entry 0 (or a `tg_growth` value 0) within
   0.006 of the neutral curve means "kc_i < 0.004", not "stable"; on the table's grid this affects no point (every grid
   point under the curve is positive).
6. Michalke (1964) digits and Chandrasekhar's odd-mode value could not be read in a primary source (labelled
   "approximate" / two-route evidence in SOURCES.md) — unchanged from loop 1.

## Documentation corrections required (analysis / design / docstrings; the verifier does not edit them)
1. Design Part C 4.5 (`double_diffusive_sigma`) expect "+0.0308" at τ = 0.0107: the root of the cubic (re-derived from the
   dimensional equations, D13) is **+0.0330081**; 0.0308 is the root for τ = 0.01 exactly. The implementer's value is right.
2. Design Part C 7.12 expect "Re = 10⁴: unstable k from ≈ 0.80 to ≈ 1.07": computed **0.797 to 1.0947** (c_i(1.07, 10⁴) =
   +0.00128; sign change between 1.09 and 1.095, confirmed by the independent row-replacement solver).
3. Design Part C 7.2: `os_3d_eigs(k, m, Re, U, Up, …)` — the fifth argument is U′ (Part C says Upp).
4. Design Part C 4.2/4.4 expect values (Ra = 875.9, Rs = 62 130, lhs = 61 254) use g = 9.81; the code default is G0 = 9.80665
   (875.6, 62 109, 61 233). The notebook must state which g it uses.
5. Design Part C 7.11 / curation C13: Tatsumi & Kakutani's k ≈ 0.2 is approximate; our k_c = 0.1728 (Re_c 4.017) is confirmed
   by an independent half-domain solve — say "Re_c ≈ 4 at k ≈ 0.17 (0.2 in the 1958 paper)".
6. Design Part C 3.13 says the determinant scan starts at "1.0001 K⁴ + 1"; the code starts at max(1.02K⁴, 50) — docstring and
   design should agree (behaviour verified: Chebyshev = determinant to 1e-8).
7. Chandrasekhar's K_c "3.117" is a rounding of 3.1163 (2.2e-4); quote 3.1163 (ours) beside 3.117 (source).

## Discrimination (planted wrong variants on a scratch copy of `fluidpy/`)
25 planted variants, each patched into a scratch copy of the repository (`fluidpy/`, `tools/`, `scripts/`,
`tests/test_ch11.py`, `reference/ch11`, no cache) and run with the targeted `-k` subset of the default tests. **25 of 25
killed** (all re-run in loop 2 against the 128-test file); baseline on the unpatched copy: 68 passed (233 s).

| # | Planted wrong variant | Where | Tests that fail |
|---|---|---|---|
| M01 | `decay_box` returns the fixed box (30) | `core/stability.py` | 10: `test_tg_growth_V3_converged_in_the_decay_box_at_small_k`, `test_decay_box_V1_…`, `test_decay_box_V3_…`, `test_tg_growth_V7_default_is_the_decay_box_…`, `test_tg_tables_V1_…`, `test_decay_map_scale_V1_…` (y_max/s = 480), `test_tg_growth_V1_shooting_…`, `test_tg_growth_V7_no_growth_at_or_above_…`, `test_tg_growth_V3_near_neutral_sliver_…`, `test_parallel_profile_V7_…` |
| M02 | `tg_growth` default box back to 30 (the loop-0 defect) | `ch11_instability.py` | 8 (the box tests, the tables test, shooting, no-growth, sliver, N² = 0 on the rule grid) |
| M03 | cache key of `tg_growth_map` without the box rule | `ch11_instability.py` | 1: `test_tg_growth_V7_default_is_the_decay_box_and_explicit_y_max_overrides` |
| M04 | cosh clip at ±3.5 instead of ±350 | `richardson_profiles` | 7 (incl. `test_richardson_profiles_V7_…`, shooting, map-scale V3, sliver) |
| M05 | sign of N² in (11.61) | `_tg_solve` | 14 (N² = 0 parity on the rule grid passes by construction; everything else on the TG solver fails, incl. shooting) |
| M06 | sign of U″ in the linear-in-c matrix of (11.61) | `_tg_solve` | 12 (incl. both N² = 0 parity tests and shooting) |
| M07 | sign of the viscous term of (11.79) | `_os_solve` | 4: Orszag eigenvalue, row-replacement route, spectral convergence, energy budget |
| M08 | Howard radius = full range instead of half | `howard_semicircle`, `in_howard_semicircle` | 1: `test_howard_semicircle_V1_boundary_and_membership` (the V7 theorem check cannot see a larger circle, as expected) |
| M09 | sign of the shear term in (11.18) | `kh_phase_speed` | 9 KH tests |
| M10 | K⁴ instead of K² in (11.44) | `benard_free_free_Ra` | 1: `test_free_free_V1_formula_examples_and_chebyshev_route` |
| M11 | κ and κ_s swapped in (11.46) | `double_diffusive_margin` | 1: `test_salt_finger_V1_criterion_is_Rs_minus_Ra_…` |
| M12 | factor ½ dropped in (11.54) | `taylor_critical_approx` | 1: `test_taylor_critical_V1_11_54_accuracy_and_cached_table_is_live` |
| M13 | sign in the denominator of r_H | `lorenz_hopf_r` | 2: `test_lorenz_V1_…_V5_hopf_point`, `test_reference_files_V7_…` |
| M14 | d³ instead of d⁴ in (11.21) | `rayleigh_number` | 2: `test_rayleigh_number_V2_…`, `test_salt_finger_V7_single_component_limit` |
| M15 | dissipation doubled in (11.88) | `disturbance_energy_budget` | 1: `test_energy_budget_V4_closes_for_computed_modes_…` |
| M16 | `decay_map_scale` returns the fixed 0.5 (the loop-1 defect) | `core/stability.py` | 7: `test_tg_growth_V3_no_false_zero_under_the_neutral_curve_at_large_k`, `test_decay_map_scale_V1_…`, `test_decay_map_scale_V3_…`, `test_tg_growth_V1_shooting_…`, `test_tg_growth_V7_no_growth_at_or_above_…`, `test_tg_growth_V3_near_neutral_sliver_…`, `test_tg_growth_V7_default_is_the_map_scale_rule_…` |
| M17 | `tg_growth` ignores the rule (default scale 0.5) | `ch11_instability.py` | 7: the same without the formula test, plus `test_taylor_goldstein_V1_N2_zero_equals_rayleigh_on_the_rule_grid` |
| M18 | floor s_min = 0.05 instead of 0.035 | `decay_map_scale` | 5: formula, shooting, map-scale V3, sliver, wiring |
| M19 | s·k = 0.25 instead of 0.025 | `decay_map_scale` | 6: formula, shooting, map-scale V3, no-growth, sliver, wiring |
| M20 | cache key of `tg_growth_map` without the map-scale rule | `ch11_instability.py` | 1: `test_tg_growth_V7_default_is_the_map_scale_rule_and_explicit_scale_overrides` |
| M21 | no cap: s = max(0.035, 0.025/k) (long waves get s > 0.5) | `decay_map_scale` | 1: `test_decay_map_scale_V1_…` only — the growth-rate tests cannot see it (the cap acts for k < 0.05, below the table; (0.03, 0.02) and (0.02, 0.01) agree with shooting to 4e-7 with the cap) |
| M22 | `benard_marginal_Ra_det` returns NaN instead of raising | `ch11_instability.py` | 1: `test_benard_guards_V7_…` |
| M23 | guard removed in `benard_growth_rate` (bare `IndexError` again) | `ch11_instability.py` | 1: `test_benard_guards_V7_…` |
| M24 | Rayleigh table back to the fixed profile boxes | `rayleigh_spectrum_table` | 1: `test_rayleigh_spectra_V1_cached_json_matches_live` |
| M25 | `_sech2` without the clip | `ch11_instability.py` | 1: `test_parallel_profile_V7_sech2_far_field_is_finite_and_the_clip_is_invisible` |
