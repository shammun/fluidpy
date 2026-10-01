# Chapter 11 — lesson review                                   2026-10-01

coverage_check: **OK** (0 errors, 0 warnings) · sections covered **14/14** (11.1–11.14) · CORE blocks **15/15** with
code + visual · notebook read in order: 613 cells (462 markdown, 151 code), 0 execution errors · derivations audited
**25/25 (258 steps)** · new primers checked 25 (P255–P279) · explainer cells 9/9 with "why interactive" + "What to try" ·
figures looked at 22 of 33 · book equations compared with page images: **49 numbered + 2 unnumbered, all match**.

Reviewed: `outputs/ch11/executed.ipynb` (source `notebooks/ch11_instability.ipynb`). Scratch files of this review:
`outputs/ch11/review_dump.txt`, `outputs/ch11/review_img/`, `outputs/ch11/_scan2.py`, `_scan3.py`, `_chk.py`.

## Final verdict: PASS (round 3, at the end of this report)

## Round 1 verdict: FAIL — 7 small, local Must-fix items; no equation, sign or physics statement is wrong

The physics, the numbers and every derivation step are right. What blocks PASS is one broken formula, one figure whose
text describes a curve that is not on the plot, two primers placed after their first use, two terms never explained, and
one ★★★ check cell that does not test the result it claims to test. One builder pass should close all seven.

## Must fix

1. **Cell 34 · broken formula.** Item 2 of "What does the code above do?" reads `$sigma_i` + line break + `eq0$` (a `\n`
   was eaten as a newline and the backslash of `\sigma` is missing). Write `$\sigma_i\neq0$` (raw string in the builder).
2. **Cells 239 / 241 · C05 figure, panel (a).** The legend lists "n = 2" (dashed) and *What you see* says "one half-wave
   (n = 1, purple) and two (dashed)", but the n = 2 curve never enters the panel: its minimum is 108π⁴ ≈ 10 520 and
   `ylim=(0, 6000)`. Raise the limit (or use a log axis) so the curve is visible, or drop the curve, its legend entry and
   the words.
3. **Cell 183 (D10) uses P263 before it is explained.** "Tools we use (each explained before this point)" lists the
   neutral-curve recipe P263, and step 14 cites it, but the primer is at cells 188–189. Move 188–189 in front of D10 (after
   the cube-root primer, cells 181–182).
4. **Cell 340 (D17) uses "critical layer" / P269 before it is explained.** Step 7 ("not at a critical layer, P269") and
   *What it means* use the term; the singular-point primer is at cells 346–347. Move 346–347 in front of D17.
5. **"self-adjoint" is never explained** (C09 idea sketch cell 362; D18 plan, step 2 *why*, step 8 *in words*, cell 367;
   note N76 cell 369). Add one sentence where it first appears: the form (p φ′)′ − q φ = 0, whose derivative term
   integrates by parts into −∫p|φ′|² with no boundary term — which is the only property D18 uses.
6. **"pitchfork" is never explained** (D25 title and step 7 "Read the pitchfork", cell 574; N116, cell 575). Add one
   sentence in step 7 or in P274: one fixed point loses stability and two symmetric ones (C₊, C₋) appear.
7. **D24 (★★★), cell 570 · the check does not test the result derived.** The sympy cell verifies steps 9–10 and b = 8/3
   only; step 7 (the vorticity projection) and steps 11–15 (the rescaling that yields
   $\dot X=\Pr(Y-X),\ \dot Y=-XZ+rX-Y,\ \dot Z=XY-bZ$ (11.91)) are not checked, and the *Check it* line cites
   `ch11.lorenz_sympy()`, which no cell runs. Add the projection on cos πz sin kx and the substitution X = α_L A,
   Y = β_L B, Z = γ_L C, or print `ch11.lorenz_sympy()["scaled_residuals"]`. (Done independently here in
   `outputs/ch11/_chk.py`: step 7 residual 0, the three rescaled equations reproduce (11.91) with residuals 0, 0, 0 —
   the derivation is right; only the notebook's check is incomplete.)

## Should fix

1. **Cell 135** — "Air is a hundred times harder to set convecting than water" sits next to 16 K (1 cm of air) and
   0.93 K (5 mm of water), a factor 17. At equal depth the factor is about 140. Say which comparison is meant, with a
   computed number.
2. **N² equation number (cells 328, 336 step 6).** The notebook writes the equation with "(Ch. 7 §7.8)" and no number.
   The number is **(7.127)** (see below); add it. `analysis/ch11.md` (4 places), `analysis/ch11_curation.md` row R18 and
   `analysis/ch11_design.md` (5 places) say (7.128) — wrong, correct them upstream.
3. **K_c = 3.117 vs 3.1163.** The C04 heading, D10 result and the summary quote Chandrasekhar's 3.117; the code prints
   3.1163 and the figure labels 3.116. One sentence: our converged value is 3.1163; the 1961 table rounds to 3.117.
4. **Check lines that no cell backs.** D22 cites `rayleigh_identity_check` (never run); D25 cites the r = 20 eigenvalues
   −13.36, −0.155 ± 8.71i (correct — re-run here — but not computed in the notebook). Add one-line cells or drop the claims.
5. **Compressed steps.** D13 steps 1, 3, 4, 6 ("as D05 steps 4–8", "the D06 moves", "as D07 step 2") and D04 step 7
   ("D03 steps 10–12 with ρ₂ replaced by ρ₂ coth kh") point back instead of showing the line; write the salt-equation
   line and the quadratic once.
6. **D14 step 14** — "$\hat u_\varphi\to\frac{\nu}{2k^2d^2\Omega_1}\hat u_\varphi$" does not say which side is old and
   which is new. Write û_φ(new) = (2k²d²Ω₁/ν) û_φ(old).
7. **2.5 % vs 2.4 %** — C07 tiny example step 2 says the shortcut is 2.5 % higher, D15 *Check it* says 2.4 %. Pick one
   (6250/6097.6 = 1.025).
8. **C06 regime map is for a 25 cm layer, the tiny example and slider F3 for 5 cm.** Cell 266 hints at it in a
   parenthesis; say outright why (at 25 cm the threshold offset is invisible and the thin diffusive band is resolvable;
   F3 uses 5 cm to make the offset visible).
9. **Figures.**
   - Cell 586: the 3-D panel's "Z" label collides with the "X(t)" label of panel (b). The text says the two runs stay
     together "until about t = 25"; the plot and the printed 29.14 say about 30. The separation is flat for the first
     ~12 time units and the reading notes do not say why.
   - Cell 485 (b): the lower branch of the Bickley curve has an unexplained kink near Re ≈ 20.
   - Cell 520 (b): one stray phase point at y = 0 (amplitudes vanish there) — mask it or mention it.
   - Cell 355: Howard's semicircle shows only as two fragments at the panel edges and is not mentioned in *What you see*.
   - Cell 453: the six profiles have no axis labels.
   - Cell 604 (a): the A = 2.8 "single point" is not drawn (the late cobweb is a point) — add a marker.
   - Cell 553 (b): "a single decaying line" is in fact a broad low bump near 1 cycle per time unit.
   - Cell 308 (b): the curve is a 9-point polygon and "hugs the Rayleigh line" overstates it (the gap is ~1500 units).
10. **Generic code comments.** 116 × "show the numbers computed above", 68 × "draw the curve(s)", 52 × "the panel's
    title…", 26 × "repeat the indented lines…", 10 × "the values for this call"; 150 of 1391 code lines carry no comment.
    On the lines that carry the idea (the ◆ of the minimum wind, the critical-layer lines, the separatrix) say what is drawn.
11. **Notation.** N59 and recap R14 use (r, U_θ); D14 uses (R, U_φ). One line that they are the same.
12. **Volume contraction −(Pr + 1 + b)** (cell 577 step 4, D25 *What it means*) is stated without its one-line origin:
    ∂Ẋ/∂X + ∂Ẏ/∂Y + ∂Ż/∂Z, the trace of the Jacobian.
13. **Terms without a defining clause:** "eigenfunction" (first at cell 164), "capillary length" (cell 94).
14. Cell 0 says equations "are cited by their numbers" — they are shown with their numbers; reword. Recap R23's heading
    reads "(Ch. 9 (free jets))".
15. Cell 498 (Table 11.1, recomputed) has the same rows as the book's table with our values and cited benchmarks — have
    `tools/check_public.py` confirm it is acceptable.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| term | first used | explained |
|---|---|---|
| neutral-curve recipe (P263) | 183 (D10 tools, step 14) | 188 — **after use** |
| critical layer / singular point (P269) | 340 (D17 step 7) | 346 — **after use** |
| self-adjoint | 362 | **missing** |
| pitchfork | 574 | **missing** |
| eigenfunction | 164 | only implied by P257 (cell 146) |
| capillary length | 94 | **missing** (named only) |
| volume contraction rate | 577 | **missing** (number given, origin not) |
| Reynolds stress | 1, 108 | 504 / 514 (forward pointer given — acceptable) |
| Tollmien–Schlichting wave | 0 | 434 (acceptable) |
| Galerkin | 318 | inline at 318, primer P275 at 567 (acceptable) |
| P259 | 4 (code comment "primer P259 below") | 151 (forward pointer — acceptable) |

All other primers P255–P279 sit at or before their first use; Python functions new to this chapter (`np.roots`,
`scimath.sqrt`, `linalg.eig(A, B)`, `ST.cheb`, Clenshaw–Curtis, `np.genfromtxt`, `symlog`, `twinx`, `ListedColormap`,
`Chebyshev.fit`, `np.polyfit`, 3-D axes, pandas) are explained where they appear.

## Derivation audit

| D | ★ | steps | every step follows? | why / in words | check | verdict |
|---|---|---|---|---|---|---|
| D01 | ★ | 5 | yes | ok | units, number | ok |
| D02 | ★★ | 12 | yes | ok | units, Ch. 7 limit, residuals | ok |
| D03 | ★★ | 14 | yes | ok | sympy cell 71, np.roots, limits | ok |
| D04 | ★★ | 9 | yes (re-derived: signs of the tension jump right; ΔU_min = 6.70 m/s) | ok | units, limits, number; no sympy (not required) | ok — step 7 compressed (S5) |
| D05 | ★ | 8 | yes | ok | units, sign | ok |
| D06 | ★★ | 8 | yes | ok | sympy + planted mistake | ok |
| D07 | ★ | 9 | yes | ok | dimensionless | ok |
| D08 | ★★ | 12 | yes | ok | sympy cell 161 | ok |
| D09 | ★ | 6 | yes | ok | no Pr, free–free | ok |
| D10 | ★★★ | 14 | yes | ok | sympy cell 184, two routes, slip #1 | ok — **P263 after use (M3)** |
| D11 | ★★ | 11 | yes | ok | sympy cell 228, slips #2, #3 | ok |
| D12 | ★★ | 9 | yes — the added step 9 is right (Pr = 1 gives σ₊ = −a² + √(RaK²/a²) = 11.0155) | ok | number vs Chebyshev | ok |
| D13 | ★★ | 12 | yes (scalings re-derived) | ok | limit, number | ok — compressed (S5) |
| D14 | ★★★ | 15 | yes (steps 8, 10, 13, 15 re-derived) | ok | sympy cells 293, 296 | ok — step 14 notation (S6) |
| D15 | ★ | 6 | yes | ok | number, μ → 1 | ok |
| D16 | ★ | 7 | yes | ok | units, U = 0 | ok |
| D17 | ★★ | 9 | yes | ok | sympy cell 353 | ok — **P269 after use (M4)** |
| D18 | ★★★ | 14 | yes | ok — **"self-adjoint" unexplained (M5)** | sympy cell 368, identity to 1e-11 | fix M5 |
| D19 | ★★★ | 14 | yes | ok | sympy cell 393, identities | ok |
| D20 | ★ | 7 | yes | ok | limit | ok |
| D21 | ★★ | 10 | yes | ok | sympy + planted sign, Orszag | ok |
| D22 | ★★ | 9 | yes | ok | cited identity check not run (S4) | ok |
| D23 | ★★ | 12 | yes | ok | sympy cell 512, budget to 1e-12 | ok |
| D24 | ★★★ | 15 | yes (independently re-derived) | ok | **incomplete (M7)** | fix M7 |
| D25 | ★★ | 11 | yes | ok — **"pitchfork" unexplained (M6)** | r = 28 run; r = 20 not (S4) | fix M6 |

The 26 "why" sentences the builder added (seen in D03/11, D06/7, D07/4, 9, D08/4, D09/5, D12/3, 4, 9, D14/7, D15/1,
D16/2, 4, D17/1, 8, D19/2, D20/6, D21/8, D22/2, 8, D23/2, 6, 7, 11, D25/1, 4) are all correct.

## Equations compared with page images (all match)

(11.1) p476 · (11.13)–(11.18) p479 · instability inequality, (11.19), (11.20) p480 · (11.21) and the Boussinesq set as
printed (4.10, 4.86, 4.89) p484 · (11.39)–(11.42) p488 · (11.43), (11.44) p490 · (11.45), Ra, Rs′, Rs, (11.46) p494 ·
(11.51)–(11.53) p499 · (11.54) p501 · (11.58)–(11.63) p504 · (11.64)–(11.67) p505 · (11.71), (11.72), the unnumbered
semicircle p507 · (11.74)–(11.78) p509 · (11.79), (11.80) p510 · (11.83)–(11.86) p512 · (11.88) p519 · (11.89), (11.90)
p526 · (11.91) with r and b p528. Ch. 7: (7.127) p294, (7.96) p289.

- Slips confirmed on the pages: #2 (sin nπz), #3 (the extra 3), #5 ("see (7.96)", which is the energy), #8 (no weight Q in
  the first inequality). (11.72) is printed with ">"; the notebook derives "≥" and says so. Only (11.72) is numbered; the
  notebook says the semicircle itself carries no number.
- **(7.127) vs (7.128):** $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ is **(7.127)** (printed p. 294). (7.128) is the
  x-momentum equation that follows.
- Other cross-chapter citations agree with the earlier builders: (1.29), (3.5), (4.10, 4.86, 4.89), (5.8), (7.18),
  (7.21), (7.95), (7.96), (8.9), (9.51), (9.52), (9.71).

## Known risks from this run — findings

- **Recomputed numbers.** The double-diffusive roots +0.033 and 9.80 ± 70.25i appear only as printed output of cell 251,
  next to their parameters (K² = π²/2, Pr = 7, τ = 0.0107, (Ra, Rs) = (1000, 2000) and (−2×10⁴, −1.9×10⁶)); re-run here:
  9.796 ± 70.245i. The verifier's 9.650 ± 70.389i belongs to a different parameter set and is not quoted. Poiseuille band
  at Re = 10⁴ is printed as "0.80 to 1.09" (cell 481). Fastest tanh mode: k = 0.4449, kc_i = 0.1897, c = 0.4264i.
- **Necessary vs sufficient** — all stated correctly: Ri > ¼ everywhere guarantees stability and Ri < ¼ only allows
  instability (D18 result, N78, summary); Rayleigh and Fjørtoft necessary, with sin y as the counter-example; top-heavy
  but below 27π⁴/4 is stable (cell 261, fourth row); (11.46) is flagged as not valid for the oscillatory diffusive
  regime; Squire flagged as not valid with rotation or stratification; linear onset vs transition (N100, N101, N111);
  the Taylor–Goldstein map caveat (cell 378); the Lyapunov slope called qualitative (cell 579).
- **The 12 slips** are each taught in corrected form with a ⚠️ box.
- **Sign conventions.** Three-row Γ table in the conventions cell and again in C03; the two Ra signs in D13 step 5 and
  the slip #10 box; N² with its sign. Good.
- **Departures from the storyboard.** C06 25 cm vs 5 cm — only partly explained (S8). C12's six library profiles — fine,
  each named and described. D04 without sympy — fine for ★★. P270 describes the tan map, which is what the solvers use.
- Review-report items M1 and M2 (`salt_finger_regime`, `miles_howard_stable`) are fixed in what the notebook runs.

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- Slip boxes with a planted-mistake check: the printed form is built into the sympy or numeric check and shown to fail.
- "Trust what does not move with N": raw spectra at two resolutions overlaid (cells 164, 355, 432).
- The conventions cell: one-letter-several-meanings table, the four scalings, the three Γ conventions.
- Two independent routes to every headline number (determinant vs Chebyshev; Galerkin vs Chebyshev; 3-D primitive
  variables vs Squire + Orr–Sommerfeld).
- The RK4-vs-DOP853 disagreement used as the lesson on chaos (cell 581).
- Every live widget is paired with a slider figure that works on the page.
- Necessary/sufficient logic primed once (P255) and reused by name in every theorem.

---

# Round 2 — re-review after fix round 1                          2026-10-01

coverage_check: **OK** (0 errors, 0 warnings) · executed notebook now 615 cells (463 markdown, 152 code), 0 execution
errors · sections 14/14 · CORE 15/15 with code + visual · primers P255–P279 all at or before first use (re-scanned) ·
equations compared with page images this round: **13 more numbered + 2 unnumbered, all match**.

## Verdict (round 2): FAIL — one Must-fix left, a single sentence (cell 489)

All seven round-1 Must-fix items are closed. The one remaining item is the explanation of the kink in the Bickley
neutral curve, which states a cause that the solver does not support. Rewording that sentence is enough; a grep of cell
489 can confirm it, no third full review is needed.

## Round-1 Must-fix items — all verified closed

| # | claim | checked in the new notebook |
|---|---|---|
| 1 | `$\sigma_i\neq0$` at cell 34 | yes |
| 2 | C05 panel (a) on a log Ra axis, n = 2 visible | yes — figure looked at (cell 239); text at 241 matches it |
| 3 | P263 before D10 | yes — primer 183, D10 185 |
| 4 | P269 before D17 | yes — primer 340, D17 342; no earlier use of "critical layer" |
| 5 | "self-adjoint" explained | yes — cell 362, the form (pφ′)′ − qφ = 0 and the one property D18 uses |
| 6 | "pitchfork" explained | yes — D25 step 7 (cell 576) |
| 7 | D24 check complete | yes — cell 572 prints 0, 0, 8/3, 0 (step 7), 0 0 0 (rescaling to (11.91)) and `[0, 0, 0]` from `ch11.lorenz_sympy()["scaled_residuals"]` |

## Must fix (round 2)

1. **Cell 489 · the kink in the lower Bickley branch near Re ≈ 20.** The sentence says it "comes from our stored table —
   few points, and near-neutral modes that converge slowly — not from the physics". Tested here
   (`outputs/ch11/_bick.py`: leading sinuous eigenvalue at Re = 14.58, 19.32, 25.62 for k = 0.04–0.12, with
   (y_max, N) = (60, 80), (150, 140), (300, 200)):
   - The modes on the right of the kink are **converged**: e.g. Re = 25.62, k = 0.065 gives c = −0.1068 + 0.0092i in all
     three boxes. Slow convergence is not the cause.
   - Left of the kink the stored branch **depends on the box**: at Re = 14.58, k = 0.04 the table's box (y_max = 60) gives
     c_i = −0.0026 (decaying), the two larger boxes give +0.00125 (growing). The true lower neutral wavenumber there is
     below 0.04, not the stored 0.065.
   - Across the kink the neutral point **changes mode**: the stored c_r on the lower branch jumps from −0.015 (Re = 14.6) to
     −0.063 (Re = 19.3), and the box test shows two different leading modes on either side.

   Reword to what is known, for example: "The lower branch of the jet is indicative only. For Re ≲ 17 its position
   depends on the size of the computational box (our table uses y_max = 60; a larger box moves it to smaller k), and
   near Re ≈ 20 the neutral point passes from one mode to another, which makes the kink." Nothing else in the notebook
   depends on this branch (the critical point 4.017 at k = 0.173 is box-independent), so publishing does not need the
   table recomputed — see Should-fix 1.

## Should fix (round 2)

1. **Recompute the lower Bickley branch** in `reference/ch11/os_neutral_bickley.csv` with a box that grows as k falls
   (the rule `ST.decay_box(k)` already used for the Taylor–Goldstein solver), or drop the lower branch for Re < 17. For
   the verifier / implementer, not the notebook builder.
2. **Code comments (the partial item).** The rule is met on the lines that carry an idea (the ◆ of the minimum wind, the
   separatrix, the critical layers now say what is drawn), and uncommented lines fell to 67 of 1420. Two things remain:
   - 114 machine-made `# print: …` comments that echo the format string, 95 of them with "…" gaps and some garbled
     (`# print: )])`, `# print: np.round(mR[ … ][0], 4), … … … imaginary part … …`). These are noise; delete them or say
     what the number means. The printed lines are already explained in "What does the code above do?".
   - Plumbing comments ("the figure and its panels" × 38, "display the figure" × 33, "repeat the indented lines for each
     item listed" × 26, "the values for this call" × 10) restate syntax. Acceptable; not worth another round.
3. **Animation A1 (cell 12), last frames:** the "cap" panel is empty because both balls have left the plotted range. Hold
   the markers at the panel edge, or say in the reading notes that an empty panel means "escaped".
4. `notebooks/build_ch11.py` still contains "7.128" three times (in the substitution patterns that replace it). Harmless
   once the design files are corrected — they are (0 occurrences now) — so the dead patterns can go.

## Round-1 Should-fix items — status

Verified done: equal-depth Ra ratio (cell 135, "about 140"); (7.127) on the N² equation (328, 336) and in the three
analysis files; K_c sentence (192); D22 and D25 check cells run (447: residuals 1.6e-13 and 3.1e-13; 580: r = 20
eigenvalues and the real part −4.4e-16 at r_H); D13 steps 1, 3, 4, 6 and D04 step 7 now show their lines (re-read: all
correct); D14 step 14 written as new = factor × old; 2.5 % in both places; the 25 cm vs 5 cm paragraph (265); notation
line (284); volume contraction as the trace of the Jacobian (579); eigenfunction (146) and capillary length (94)
defined; cell 0 wording; R23 heading. Figures looked at again: 239, 455 (axis labels added), 522 (stray point gone),
588 (label collision gone, "about t = 28", flat start explained in 590), 606 (fixed-point markers), 355 with its new
sentence in 357, 310 ("rises with the Rayleigh line and stays above it"), 557 ("one low bump").

## Re-checks asked for

- **Nothing new is used before it is explained.** The primer scan shows every P255–P279 defined at or before its first
  mention (the only earlier mention is the forward pointer to P259 in a code comment of cell 4). The added text uses only
  terms already explained.
- **Animations at dpi 60, 18–22 frames.** Last frames extracted and looked at for A1 (cell 12), A2 (102) and A5 (592);
  A3 and A4 extracted. Labels are readable, and each still shows its point: A1 the four verdicts, A2 the linear wave
  leaving the panel while the sheet rolls up (printed: 0.067 vs 0.056 m at t = 0.82 s, 2.68 vs 0.197 m at the end), A5
  the two runs on different lobes with the separation saturated.
- **More equations against page images (not sampled in round 1), all match:** (11.28), (11.29), (11.30),
  (11.31)–(11.33) p486 · (11.34)–(11.38) p487 · (11.50) p498 · (11.68), (11.69), (11.70) and the divergence form p506 ·
  (11.87) and its expansion near the critical layer p514. Slips #11 (−k²(U − c)F as printed) and #12 (continuity without
  1/R as printed) confirmed on the pages.

---

# Round 3 — targeted check after the Bickley recomputation        2026-10-01

coverage_check: **OK** (0 errors, 0 warnings) · executed notebook 618 cells (465 markdown, 153 code), 0 execution errors.
Targeted check only (cells 12–14, 202, 392, 459, 487–492, 545–546, the Table 11.1 cell, comment scans); not a full re-read.

## Final verdict: PASS — no Must-fix left; four Should-fix items, none blocking

## What was checked

1. **Bickley jet (cells 487–492; both figures looked at).** Everything drawn and said agrees with the recomputed tables
   (`reference/ch11/os_neutral_bickley.csv`, `os_neutral_bickley_longwave.csv`):
   - one band up to the tabulated Re = 16.78, the gap first tabulated at Re = 19.32 spanning k = 0.0447 … 0.0759 (printed);
     the notes say "up to Re ≈ 17" and "between the two tabulated Reynolds numbers printed below" — correct;
   - ticks on k = 0.02 at Re = 8.3–16.8 and Re ≥ 91, exactly where `k_lower` is NaN, labelled "unstable down to k = 0.02
     (edge not resolved)" — correct;
   - **shading k = 0.02 … k_lower as stable gap for Re ≈ 29–79 is right**: there the long-wave band's upper edge has dropped
     below 0.02, so everything between 0.02 and `k_lower` is stable;
   - **"the gap … leaves the computed range near Re ≈ 80" is right**: `k_lower` = 0.02025 at Re = 79.11 and NaN from 91.08;
   - upper branch 1.976 at Re = 1000 (printed), critical point 4.017 at k = 0.173 unchanged;
   - "from Re ≈ 19 the two bands are carried by two different modes" matches the table's `same_mode` = 0 column.
   The round-2 Must-fix (the false cause of the kink) is gone with the kink.
2. **Blasius.** 519.06 everywhere: cell 545 (note N108, with the 519.08 fixed-box value explained), cell 546 (prints
   Re_c = 519.06, k_c = 0.30377, c_r = 0.39664, ω_c = 0.12049), Table 11.1 cell (519.060118 / 0.303771). "0.03 % below"
   519.2 is right (0.027 %). No other "519" in the notebook.
3. **Round-2 Should-fix items.** `# print:` comments: 0 left. A1 last frame (extracted and looked at): the escaped balls
   are now held at the panel edge of "cap". "7.128" in `build_ch11.py`: 0. Sin-profile caveats removed and b = 1.7 →
   0.0119 printed (cell 459). D19 step 3 is titled "Insert into the Taylor–Goldstein equation" with (11.61) written out.
   The live widget (cell 202) catches `ValueError`.
4. **(7.128) on p. 503 — confirmed on the page image.** The book prints "(7.128)" beside
   $N^2\equiv-\frac g{\rho_0}\frac{d\bar\rho}{dz}$ on p. 503; Ch. 7 p. 294 defines it as (7.127). (The same page also
   confirms slip #4.)

## Should fix (none blocks publishing)

1. **Add the (7.128) back-reference as a printed slip.** A reader following along in the book sees (7.128) on p. 503 and
   (7.127) in the notebook. One line in recap R18 (cell 328) — "the book's p. 503 labels this (7.128); the definition in
   Ch. 7 is (7.127)" — and a thirteenth row in the slips table of the conventions cell.
2. **Soften the "flag" sentence (cell 492).** The picture is defensible as a heuristic: the long sinuous wave moves the jet
   sideways almost as a whole, so the strain it adds — and with it the viscous damping, of order k²/Re — is small. But it
   is our reading, not a result shown in the notebook, and for waves this long (k = 0.02 is a wavelength of about 300 jet
   widths) the parallel-flow assumption itself is doubtful, because a real jet spreads within one wavelength. Say both:
   "one way to see it … (our reading); for such long waves the parallel-flow model is itself questionable (N109)".
3. **Two small gaps in the Bickley notes.**
   - Cell 492 does not say that for Re ≥ 29.5 the long-wave band is still there, below k = 0.02; the orange region closing
     to a point near Re ≈ 27 can be read as the band ending.
   - Cell 489 does not mention the two short blue pieces at the bottom of panel (b) (the jet's lower edges at Re 4–7 and
     19–79); only a printed line points to the next figure.
4. **Echo comments.** "draw the curve(s)" was replaced by 33 × `# draw <y> against <x>` and 54 × `# the curve labelled: …`,
   which repeat the code. Same remedy as before: say what the curve is, or leave the line to the numbered explanation.
