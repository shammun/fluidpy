# Chapter 10 — derivation review (independent, read-only)

Verdict: four Must-fix items. One (M1) is a sign error the suite could not see, because the tests compared the live code with
tables written by the same buggy code. No transcription error in the schemes themselves. Findings are from the
`derivation-reviewer` reply (pages p451–p455, p467–p468, p471, p473, p478–p479, p481–p482 re-read); status filled in by the
orchestrator.

## Must fix
1. `fluidpy/core/mac.py` `primary_vortex_centre` (~:765): the fitted ψ_min adds −¼(f₂ − f₀)s* where the parabola vertex value
   is f₁ + ¼(f₂ − f₀)s* — the result is pushed above the discrete minimum (paraboloid with ψ_min = −0.1 → −0.09505; 64² cavity
   −0.102920 vs discrete −0.102953). Contaminates the C15 ψ_min grid study (p 1.946 → ≈ 1.998, GCI_fine 0.19 % → ≈ 0.16 %),
   the `psi_min` headers of `reference/ch10/cavity_mac_*.csv`, `explainer_tables.json`, the verification report C15 row and
   the E7 parity row. Fix sign, regenerate tables, add a V1 paraboloid test. — *sent to concept-implementer; tests to
   math-verifier.*
2. `fluidpy/core/fd.py` `rod_heating_exact` (Eq. (10.199)) (~:1103): the series stops when the actual term is small, which
   happens where sin((2m − 1)πx/L) = 0 (x = 0.2, t = 0.001 → −0.110; x = 1/3, t = 0.01 → 0.00097 vs 0.0184). Stop on the
   envelope 4/((2m − 1)π)·exp(−D(2m − 1)²π²t/L²). — *sent to concept-implementer.*
3. `fluidpy/core/maccormack.py` `cavity_maccormack` (~:509): default `dt_rule="inviscid"` drops the viscous factor
   1/(1 + 2/Re_Δ) of (10.110), undocumented and unguarded (Re = 5 on 64², Re = 2 on 32² → NaN). Default to the full rule,
   raise on violation, mark DEVIATION if it differs from the book. — *sent to concept-implementer.*
4. Rule 9 leak: book run parameters (cavity Mach/grid, block geometry and Mach, Re = 1000 smoke-line case, FE-cylinder channel,
   σ, the book's unbounded St) in `analysis/ch10.md`, `analysis/ch10_curation.md`; a value revealing the book's printed period in
   `reports/ch10_verification.md`; `tools/check_public.py` did not catch it. — **fixed by the orchestrator**: values redacted to
   "(book value, private)" with our own values; `tools/check_public.py` now scans tracked text files for strings listed under
   `_forbidden_public` in the local, git-ignored `tests/book_values_*.json`; the unpushed ch10 commits that contained the values
   are squashed before the push so they never reach the public history.

## Should fix (all sent to concept-implementer)
- `fd.stability_verdict`: reason text wrong for upwind / Lax–Wendroff with β > 0 (true limits |C| + 2β ≤ 1, C² + 2β ≤ 1).
- `fd.numerical_diffusivity`: `scheme="upwind"` with default C = 0 silently returns the steady value — make C required.
- `fd.convergence_study` docstring promises asymptotic orders on the default grids (measured 0.76, 1.80).
- Validation labels "converged" without an asserted order (`lax_demo`, `body_forces`, `cylinder_forces`, `march_unsteady`,
  `cylinder_steady`, `primary_vortex_centre`); "benchmark" on the MacCormack cavity (1.43 % at 64²).
- `mac.step`/`run` default `lid=1.0` turns any walled grid's top into a lid — default `None`.
- `fd.stretched_grid` refines only at x = L (layer at x = 0 for R < 0) — document; `make_cavity_bc(wall_density="momentum")`
  corner treatment — document.
- Stale test comments `tests/test_ch10.py:1918–1919` (F3 fixed), `:1953` (0.47 → 0.43). — *math-verifier.*

## Verified against rendered page images
Taylor stencils; FTCS (10.10)–(10.11), BTCS (10.13) with the Neumann ghost; Lax–Wendroff 2α² = C²/2; truncation error (10.17);
G (10.24), |G|² (10.26), Noye (10.27); BTCS, CN, upwind (both signs of u, R11), Lax–Wendroff limits; Courant with |u|;
layer thickness; uneven-grid weights; discrete roots (1 + P/2)/(1 − P/2), 1 + P, 1/(1 − R); FE element matrices (h/6)[[2,1],[1,2]],
(u/2)[[−1,1],[−1,1]] + (D/h)[[1,−1],[−1,1]], f¹ = −g k¹_{a1}, assembly A = e − 1, e (R1, R2), θ-scheme and MOL Jacobian;
MacCormack predictor/corrector directions, a₁–a₁₁ (a₁₀ = 2(a₅ + a₆), c₅ = μΔt/(12ΔxΔy)), cavity wall densities
(10.139)–(10.146), block-face densities (10.151)–(10.154) (8/9, 1/18), drag/lift traction, (10.155) and R12; the five DEVIATIONs
justified; MAC staggered index map, face gradient (10.126), divergence (10.123), Poisson matrix link by link, pin validity,
FFT eigenvalues, SOR, stability limits (10.127)–(10.128), streamfunction sign, Taylor–Green and Poiseuille; P2 shapes (10.185),
7-point Radon rule, Newton blocks (10.189)–(10.197) incl. R6 (B_vpᵀ) and R10, drag reaction sign, Kovasznay; Θ-scheme
substeps, Marchuk–Yanenko, artificial-compressibility sign, both GCI formulas; Ghia and Hou data public and cited;
sign-sensitive tests present (drag from a linear pressure, cylinder trend, projection and cross-viscous mutants).
