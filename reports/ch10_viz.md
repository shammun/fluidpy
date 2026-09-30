# Chapter 10 — explainer review                                  2026-09-30

Re-run for the Phase 7 gate after the fluidpy changes (stability_verdict reason strings, ψ_min fix, MacCormack
`dt_rule="additive"`, `numerical_diffusivity` needing C, `mac.step/run` lid default `None`).

Evidence gathered this round:
- `tools/viz_lint.py --chapter ch10` → **ok** for all 8 files.
- `tools/shot.py --chapter ch10` (full: 8 sizes, every tab, step and pager page) → **7 PASS, 1 FAIL** (cell_peclet_wiggles, tablet).
- `audit.json` selftest: **317/317 rows ok** (272 fluidpy parity rows evaluated by shot.py). Tolerances are honest: 1e-9…1e-14
  for closed forms and bit-level ports; the loose ones only where the reference itself is rounded (6-digit Ghia/INFSUP tables
  2e-6/1e-5, 4-digit Kovasznay orders 5e-4, 6-digit centreline curves 2e-3, min-norm lstsq of a singular P1–P1 system 1e-6).
- A browser pass (Playwright/Edge) that changes the state *while the Explain tab is open* and compares it with a fresh page
  set to the same state (25 cases, transport paused): **no stale Explain text anywhere**; status lines checked in each case.
- A scan of every config string and live-function output, on every tab, for control characters (lost `\f`, `\t`, `\b`, `\v`)
  and `.katex-error` nodes: **0 in all 8**. KaTeX rendered at every size (no fallback).
- `tools/eq_refs.py --chapter ch10`: 12 hits, all `ref:` metadata of an equation/derivation/code card whose TeX is displayed
  with it → accepted. A wider scan of my own found one real bare number (mixed_fe_lbb) and one weak pointer (fd_stencil_order).
- `tools/check_public.py` → OK (541 files). Because `viz/ch10/` is still untracked, I also searched the 8 files for every
  string under `_forbidden_public` in `tests/book_values_ch10.json`: **0 hits**. The embedded Ghia/Hou values are cited public
  benchmark data (`reference/ch10/`), not book numbers.
- Derivation tabs against the notebook (`notebooks/build_ch10.py --dump`, i.e. PF after all `pf_sub`/`pf_why_add` edits):
  **every D id of the storyboards is present, every step count matches** (D01 9, D03 11, D04 11, D05 8, D06 13, D07 7, D08 9,
  D18 14, D15 8, D16 12, D17 8, D09 9, D11 13, D12 10, D14 11, D19 9, D20 12, D21 7, D22 11, D23 8) and every `did` title
  matches; the only differences are titles in which the explainer replaces a bare "(10.xx)" by words.
  Walkthrough `derive: {id, step}` links all point at the step they quote.
- Equations checked against the rendered pages p453 (10.13–10.17), p455 (10.26–10.31), p466 (10.93–10.94), p471 (10.121–10.125),
  p473 (10.127–10.128), p474 (10.134–10.136): signs and factors match in the JS, the Equations tab and the derivations.
- The fluidpy changes: the lid-driven cavity's `vortexCentre` uses the corrected ψ_min vertex formula and its steady/unsteady
  ψ_min rows pass at 1e-11; the cached MacCormack table rows (max dev, mass drift) re-run against the new `dt_rule="additive"`
  default and pass at 1e-12; every Code tab that calls `numerical_diffusivity` passes C or uses `"upwind_steady"`; the cavity
  Code tab passes `lid=1.0` explicitly. **Only the stability_verdict strings are stale** (von_neumann: a real, visible bug;
  upwind_cfl: dead code).

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| fd_stencil_order | ok | PASS 8 / 280 | 25/25 (22 py) | 7 §, live, regime-dependent (stencil + FTCS modes) | D01 9/9, D03 11/11 ✓ | strong; 1 quiz item not answerable in the UI | good | **PASS** |
| von_neumann_amplification | ok | PASS 8 / 472 | 41/41 (36 py) | 8 §, table, live | D04 11, D05 8, D06 13, D07 7 ✓ | strong | status wrong for upwind β > 0; "1 − −1" in §2 | **FAIL** |
| upwind_cfl_advection | ok | PASS 8 / 312 | 32/32 (29 py) | 9 §, table, live, end card | D08 9, D18 14 ✓ | strong | good | **PASS** |
| cell_peclet_wiggles | ok | **FAIL** 8 / 360 (tablet tour-step2) | 34/34 (30 py) | 8 §, table, inspector | D15 8, D16 12, D17 8 ✓ | strong | one clipped derivation mini-card | **FAIL** |
| fem_hat_assembly | ok | PASS 8 / 504 | 32/32 (27 py) | 8 §, scatter-add table, inspector | D09 9, D11 13, D12 10, D14 11 ✓ | strong | good | **PASS** |
| mac_projection_staggered | ok | PASS 8 / 368 | 38/38 (35 py) | 8 § (+ checkerboard mode 6 §) | D19 9, D20 12, D21 7 ✓ | strong | 16² illegible on phones | **PASS** |
| lid_driven_cavity | ok | PASS 8 / 176 | 41/41 (36 py) | 7 §, benchmark table (row lit), live | D23 8/8 ✓ | strong | small legend gap | **PASS** |
| mixed_fe_lbb | ok | PASS 8 / 200 | 57/57 (47 py) | 8 §, β table, live | D22 11/11 ✓ | strong | bare "(10.184)" in Explain | **FAIL** |

## fd_stencil_order
**Must fix** — none.

**Should fix**
1. Check yourself Q2 ("Why is the central stencil exact for f = x² but not for x³?") cannot be answered by experimenting: the
   test functions are sin x, eˣ, e^(−x²). Add a `poly` chip (x² / x³, exact derivatives are trivial in `fk`) or reword the
   question to one of the offered functions (e.g. "why is the Taylor-terms bar c₃ the first rose one for sin x").
2. D01 step 6 *why* ends "…the third line of (10.6)" without the equation; write "…the third line of
   $\big[\frac{\partial T}{\partial x}\big]_i=\frac{T_{i+1}-T_{i-1}}{2\Delta x}+O(\Delta x^2)$ (10.6)". (The step titles
   "Repeat with the left neighbour (10.5)", "Subtract (10.5) from (10.4)", "Add (10.4) and (10.5)" are the notebook's and the
   start page shows both series — acceptable.)
3. FTCS mode colours: the stencil nodes are orange while the "time" term is amber; on small screens the two read as the same
   colour. Use a clearly different hue for one of them.

**What works** — log–log |error|(h) with the rose round-off band, the ◆ best-h marker and the faint "so far" curve on one
transport (h shrinks tenfold every 0.8 s); Taylor-term bars with exact fractions (c₃ = 1/6, "cancels" labels) next to the
measured error; the ×1/h zoom inset that keeps the nodes visible below h = 10⁻²; the FTCS mode whose three truncation terms
add up to the measured one-step residual (0.0317 vs 0.0317 1/s) and the β = 1/6 cancellation preset.

## von_neumann_amplification
**Must fix**
1. **Stale stability_verdict reason for upwind with diffusion** (JS `verdict()`, line 2413):
   `reason = Math.abs(2 * a) > 1 + 1e-12 ? 'unstable: C > 1, …' : 'unstable: 2b > 1, the zigzag grows';`
   fluidpy now returns `"unstable: |C| + 2b > 1, the zigzag grows"` here. Visible effect (browser-verified): scheme = upwind,
   α = 0.3, β = 0.3 → the status reads **"💥 unstable: 2β > 1, the zigzag grows"** although 2β = 0.6 (the real failure is
   C + 2β = 1.2 > 1); the Code tab's `# {{reason}}` comment prints the same wrong string, so it no longer equals what the shown
   Python prints. Change the fallback to `'unstable: |C| + 2b > 1, the zigzag grows'`, extend `pretty()` with
   `.replace('|C| + 2b > 1', '|C| + 2β > 1')` (and the phone shortening), and add the parity row
   `{ name: 'verdict text: upwind |C|+2b', js: eqs(verdict('upwind', 0.3, 0.3).reason, 'unstable: |C| + 2b > 1, the zigzag grows'), py: 'ch10.stability_verdict("upwind", 0.3, 0.3)["reason"] == "unstable: |C| + 2b > 1, the zigzag grows"', rtol: 0 }`
   so the next string change is caught.

**Should fix**
1. Explain §2 (upwind) prints `=1-0.6(1--1+0i)-1.2` at θ = π (screenshot `scratchpad/ch10_stale/von_neumann_amplification_case0_A.png`):
   the template `1-${t3(n.C)}(1-${t3(cs)}+${t3(sn)}\mathrm i)` does not handle a negative cos θ or a round-off sin θ. Write the
   bracket with a sign-aware helper (e.g. `(1-(${t3(cs)})+${t3(sn)}\,\mathrm i)` or `ctex`), as the FTCS branch does.
2. Explain §3 for upwind could add the Noye-style two-row table (C ≤ 1, C + 2β ≤ 1) so the failing edge is lit, as FTCS has.

**What works** — three linked views on one clock (G(θ) curve with the outside part rose and the Re/Im split of the probe
wave, Noye's lens drawn from the closed form *and* checked by a brute-force scan of rose dots, the live 10⁻¹⁰ kick marched to a
×10⁶ cap with the predicted slope max|G|ⁿ dashed); the draggable (α, β) dot; the classic-settings table; the exact
fluidpy reason strings pinned by selftest rows for FTCS (a pattern to extend to every scheme).

## upwind_cfl_advection
**Must fix** — none.

**Should fix**
1. Mirror of stability_verdict is stale (JS `verdict()`, line 2480: final `else reason = 'unstable: 2b > 1, the zigzag grows'`).
   With β = 0 the branch is unreachable, so nothing wrong is displayed, but the function claims to be `ch10.stability_verdict`;
   change it to fluidpy's wording (`'unstable: |C| + 2b > 1, the zigzag grows'` for upwind, `'unstable: C^2 + 2b > 1, the zigzag grows'`
   for Lax–Wendroff/MacCormack) so the next copy is not inherited wrong.
2. At C > 1 the Code tab comment shows `D_num = -2.632e-4 m^2/s`; add "(negative: anti-diffusion, C > 1)" to the comment or
   Explain §2 so a negative diffusivity is not read as a bug.

**What works** — the x–t stencil view with the characteristic's foot ◆ and the interpolation weights (×1.05 / ×−0.0526 when
CFL fails) — the CFL idea made visible; upwind vs MacCormack vs FTCS on one ring with an end-of-run card; |G| and phase-speed
panels with the pulse's main wave marked; the "u = −1, printed" mode that demonstrates the book's sign slip (R11); D18 ★★★ in
14 short steps with live local-error coefficient that vanishes at C = 1.

## cell_peclet_wiggles
**Must fix**
1. **shot.py FAIL**: `[tablet] tour-step2: viz-step-card page 1/1: viz-der-mini 'Derivation · step 7: Let R be large…' 17px too wide`
   (`reports/viz/ch10/cell_peclet_wiggles/CLIP__tablet__tour-step2__pg1.png`: the equation number "(10.88)" is cut off at the
   right edge of the mini-card). D15 step 7 tex is `T\approx\frac{e^{Rx/L}}{e^R}=e^{-R(1-x/L)}\quad(10.88)`; split it into
   `\begin{aligned}T&\approx\frac{e^{Rx/L}}{e^R}\\&=e^{-R(1-x/L)}\ \ (10.88)\end{aligned}` (or drop the `\quad` and put the
   number on the second line). Re-run `tools/shot.py viz/ch10/cell_peclet_wiggles.html`.

**Should fix**
1. The derivation header "EQS. (10.90)–(10.93)" on the step pages is only a number range; acceptable as metadata, but the goal
   page could name the two schemes' equations as it already does for (10.91)/(10.93) in the steps.

**What works** — the profile + discrete-root + diffusivity/error views: the root plot with the rose r < 0 band and the ± sign
strip of rʲ makes "wiggles iff R_cell > 2" a one-glance fact; the sweep transport that stops at the first wiggle (n = 19,
pinned by an invariant row); the inspector that traces T_j = (rʲ − 1)/(rⁿ − 1) with the tiny example (−0.05, 0.1, −0.35);
the honest stretched-grid reading ("shrink, not removed") instead of the storyboard's over-claim; the dashed modified-equation
ghost with its measured gap.

## fem_hat_assembly
**Must fix** — none.

**Should fix**
1. On phones the matrix view hides the numbers for n ≥ 4, so walkthrough step 4 ("Click any matrix cell…") depends on the
   inspector; consider printing the three entries of the highlighted row under the matrix on phones (the title already has
   "row 2: (−1.5, 2, −0.5)" — good; keep it on every step).

**What works** — the element loop as a transport (each element's 2 × 2 block outlined amber as it lands), the scatter-add
table with the current element lit, the matrix-cell inspector ("K₂₂ = k⁽²⁾₂₂ + k⁽³⁾₁₁ = 1.5 + 0.5 = 2"), FE = FD to 10⁻¹⁷ as a
status and an invariant row, the natural-end mode whose learned slope (0.667 on 4 elements) teaches "natural = approached";
four derivations (D09, D11, D12, D14) that match the notebook step for step.

## mac_projection_staggered
**Must fix** — none.

**Should fix**
1. Phones: the storyboard caps n at 12 on phones; at 16² (preset "16² squeeze", tour step 7, `phone-tall__explain.png`)
   the face arrows are 1–2 px and the divergence borders a speckle — the picture carries no information. Cap n at 12 in
   portrait (or draw every other arrow) and say so in the hint.
2. Tour step 1 says "Press ▶" but the step leaves the transport paused (fine), while step 1 of the storyboard wanted the 8²
   preset; the default shear field is used — consistent with the text, just note it in the design.

**What works** — one projection in slow motion (u* → Poisson → push → uⁿ⁺¹) with the cell-border divergence colours turning
white; the cell inspector with the wall face marked "wall" (never corrected); the 3-cell pipe tiny example (p − p₁ = 0, 1, 2);
SOR history vs direct solve and "divergence left = Δt × residual" pinned by an invariant row; the checkerboard mode with the
row view and the null-space bars (collocated 4, staggered 1).

## lid_driven_cavity
**Must fix** — none.

**Should fix**
1. Grid-convergence panel: the purple diamonds (MacCormack deviations) are not in the view title's legend nor in Explain §0;
   add "◆ MacCormack" to both.
2. While a live grid is spinning up on a phone, three clocks disagree in one frame (`phone-tall__tour-step3.png`: transport
   0.6 L/U, view title t 0.26, readout t 0.51). Either show "computing … t = 0.26" in the status whenever the picture lags the
   transport by more than one snapshot, or hold the transport until the run catches up.

**What works** — a real MAC solver running live (16²–32²) with a banded LU factorised once, the centreline against Ghia's
17 points, cached 64²/128² and MacCormack runs read from a parity-checked table, the three-grid Richardson study (p = 2.00,
ψ_min ≈ −0.10353, GCI) with the "128² is further from Ghia than 64² — the comparison measures Ghia" interpretation, and the
benchmark table with the current row lit. All ψ_min rows use the corrected vertex formula.

## mixed_fe_lbb
**Must fix**
1. **Bare equation number in the Explain tab**: with "one element" on, Explain §0 reads "the parent triangle of **(10.184)**"
   (line 2978) — the equation is never shown there. Write it: "the parent triangle and its isoparametric map
   $x(\xi,\eta)=\sum_{a=1}^6x^e_a\phi_a(\xi,\eta),\ y(\xi,\eta)=\sum_{a=1}^6y^e_a\phi_a(\xi,\eta)$ (10.184)". (`tools/eq_refs.py`
   did not catch this one; it sits in a template string with other `$…$` further along.)

**Should fix**
1. Tour step 5 claims "the noise is gone; only the two lid corners stand out"; at n = 4 the P2–P1 picture (colour saturated at
   ±3.77) still shows a broad violet band along the lid. Either say "a smooth pattern plus the two singular lid corners" or
   show the lower-half spread number (already computed in Explain §4) in the step text.

**What works** — counting that becomes an argument (n_u vs n_p − 1, "at least 6 invisible patterns" at n = 2 matching the 6
spurious modes found), the "invisible pattern" picture with max|Bq| = 4×10⁻¹⁶, β_h vs n with P1–P1's zeros on the floor and its
falling next eigenvalue, live dense FE up to 6 × 6 with a Cholesky + Jacobi generalised eigen-solver (57 parity rows), the
one-element mode with the six P2 shapes and the 7-point rule, and the Kovasznay order check (3.00 / 2.23).

### Round-1 verdict: FAIL
Three explainers need one small edit each before the gate passes:
- **von_neumann_amplification** — Must 1: update the upwind fallback reason to `"unstable: |C| + 2b > 1, the zigzag grows"` (the
  status is wrong today for upwind with β > 0) and add a parity row for it.
- **cell_peclet_wiggles** — Must 1: split the D15 step 7 formula so the tablet mini-card is not clipped (shot.py FAIL).
- **mixed_fe_lbb** — Must 1: write the isoparametric map next to "(10.184)" in the one-element Explain text.

fd_stencil_order, upwind_cfl_advection, fem_hat_assembly, mac_projection_staggered and lid_driven_cavity **PASS** (Should-fix
items above are optional polish). After the three edits, re-run `tools/viz_lint.py --chapter ch10` and
`tools/shot.py viz/ch10/{von_neumann_amplification,cell_peclet_wiggles,mixed_fe_lbb}.html`.

---

## Round 2 — 2026-09-30 (after the builders' fixes)

Re-checked: `tools/viz_lint.py --chapter ch10` → ok ×8; `tools/eq_refs.py` → the same 12 `ref:` metadata hits (accepted);
`tools/check_public.py` → OK. The four untouched files (fd_stencil_order 14:18, fem_hat_assembly 14:50,
mac_projection_staggered 15:05, lid_driven_cavity 15:20) have the same mtimes as in round 1 → round-1 PASS stands.
`tools/shot.py` re-run on the four edited files (von_neumann twice, upwind twice), plus a browser pass on the states named
in round 1 (`scratchpad/ch10_r2/*.png`).

| slug | round-1 item | round-2 check | shot (this run) | parity rows | verdict |
|---|---|---|---|---|---|
| von_neumann_amplification | Must: stale upwind reason | **fixed**: upwind α = β = 0.3 → status "💥 unstable: \|C\| + 2β > 1, the zigzag grows · max\|G\| = 1.400 at θ = π"; Explain §3 "\|C\| + 2β = 0.6 + 0.6 = 1.2 > 1"; 3 exact-text parity rows + 2 `pretty()` rows. Should (θ = π §2): **fixed**, now "= 1 − 0.6(2 + 0 i) − 1.2 → G = −1.4" with the (1 − cos θ) + i sin θ note; FTCS line also sign-clean | **FAIL ×2 (reproducible)**: `[phone] equations: pager page 3/4: viz-eq 'Pure diffusion Eq. (10.28) …' clipped 26px (cut through a line)` | 47/47 (40 py) | **FAIL** |
| cell_peclet_wiggles | Must: tablet tour-step2 clip | **fixed**: D15 step 7 is two aligned lines | PASS | 34/34 | **PASS** |
| mixed_fe_lbb | Must: bare (10.184) | **fixed**: element-mode Explain §0 writes x(ξ,η) = Σ xᵃφₐ, y(ξ,η) = Σ yᵃφₐ (10.184); tour step 5 now quotes the lower-half spread | PASS | 57/57 | **PASS** |
| upwind_cfl_advection | Should: stale verdict mirror | **fixed**: `verdict()` mirrors stability_verdict incl. β > 0 (\|C\| + 2b / C^2 + 2b), 4 new parity rows | run 1 **FAIL** `[phone-land] explore: pager page 5/6: viz-notes 'Stable, but not free…' clipped 1px`; run 2 PASS → intermittent | 36/36 (33 py) | **PASS** (flaky 1 px — see Should) |
| fd_stencil_order, fem_hat_assembly, mac_projection_staggered, lid_driven_cavity | — | unchanged | PASS (round 1) | 25/25, 32/32, 38/38, 41/41 | **PASS** |

### Must fix (round 2)
**von_neumann_amplification**
1. `tools/shot.py` FAIL, reproduced in two runs: on the 360×640 phone the Equations tab, page 3/4, clips the "Pure diffusion"
   card (10.28) by 26 px through its live line (`reports/viz/ch10/von_neumann_amplification/CLIP__phone__equations__pg3.png`).
   The audit reaches Equations with the huge-step state (β = 100, from the BTCS preset/quiz), so the (10.26) and (10.27)
   live lines are long ("1.59 × 10⁵ + 0 = …", "0 ≤ 0.04 ≤ 200 > 1") and the page packs one card too many. Fix: shorten the
   (10.28) card — put the two relations on one line `0\leeta\le	frac12\iff\Delta t\le	frac12rac{\Delta x^2}{D}` and the
   live text on one line (`\Delta t\le0.005	ext{ s}` with Δx, D in the note) — or mark the long live lines `optional`; then
   re-run `tools/shot.py viz/ch10/von_neumann_amplification.html` until PASS.

### Remaining Should-fix (compact)
- upwind_cfl_advection: the phone-land Explore "Right now" note clips by 1 px intermittently (FAIL in one of two runs) — trim
  one clause or allow one more line of slack so the audit is stable.
- von_neumann: cosmetic "−0 i" / "+0 i" at θ = π in §2 (round-off sin θ); print the real number alone when |Im| < 1e-12.
- fd_stencil_order: quiz Q2 (x² vs x³) not answerable in the UI; D01 step 6 *why* "the third line of (10.6)" without the formula;
  amber vs orange too close in FTCS mode.
- mac_projection_staggered: cap n at 12 on phones (16² illegible).
- lid_driven_cavity: MacCormack ◆ missing from the convergence legend / Explain §0; spin-up clocks disagree on phones.
- fem_hat_assembly: keep the row-of-node-A numbers visible on phones for tour step 4.

## Round-2 verdict: FAIL
7 of 8 explainers PASS; **von_neumann_amplification FAILS** on one reproducible phone clip (Equations page 3, the (10.28)
card, 26 px). All round-1 Must-fixes are resolved and every parity row passes (fd 25, vn 47, upwind 36, cell 34, fem 32,
mac 38, lid 41, mixed 57 = 310/310).

## Round-2 resolution (orchestrator)
- von_neumann_amplification (last fix round): (10.28) card on one line with a one-line live value; (10.26)/(10.27) live
  lines shortened when β > 1; round-off imaginary parts dropped ("−1.04", no "± 0i"). Full shot.py PASS twice (8 sizes,
  472 views, 47/47 rows); own 360×640 screenshots of every Equations page with β = 100 under BTCS and FTCS: no overflow.
  Re-run by the site-publisher: PASS.
- upwind_cfl_advection: phone-landscape "Right now" note shortened; full shot.py PASS three runs in a row (36 rows).
- cell_peclet_wiggles: an intermittent 11 px phone clip of the (10.88) Large-R card under machine load (merge-gate run)
  was fixed with margin (card and live line shortened, neighbouring cards trimmed); full shot.py PASS three runs in a row
  and PASS under load in the merge gate.
- Merge gate: tools/viz_lint.py --chapter ch10 ok (8/8); tools/shot.py --chapter ch10 --quick 8/8 PASS; published page
  audit (shot.py --page) PASS.
- Remaining Should-fix items above are carried as open items to knowledge/viz_patterns.md.

## Verdict: PASS
