# Chapter 9 — lesson review                                   2026-09-30

coverage_check: OK (0 errors, 0 warnings) · sections covered 11/11 (§9.1–§9.11, exercises as S01) · CORE blocks 14/14 with code+visual · derivations 22/22 present (★★★: D14, D16, D19, D20, D21 all have a sympy cell) · explainers 9/9 embedded, each with "Why interactive" + "What to try"

Notebook read: `outputs/ch09/executed.ipynb` (454 cells, 0 errors; all 22 static figures looked at). Review-driven corrections from `reports/ch09_review.md`: (9.8) sign ✓ (cell 33–35), R5 relabelled as a trap ✓ (cells 6, 127, 133), R17 ✓ (cell 161), R11 wording ✓ in cell 162 **but undone in cell 169** (Must 4), Example 9.2 book criterion −0.090 first, then the exact FS fold −0.068148 ✓ (cells 227–229), Ψ in m⁵/s³ ✓ in cell 419 **but still m⁴/s³ in D19** (Must 3).

## Must fix

1. **Cell 49 (reading notes of the C01 figure, cell 48 panel b).** "the cross-stream pressure bar is smaller still" is false: in every group the orange cross-stream-pressure bar is *larger* than the grey dropped-diffusion bar (≈1.7e−3 vs 1.3e−3 at Re_L = 10³; same ratio at 10⁵, 10⁷). Change to: "the cross-stream pressure term is of the same small order as the dropped x-diffusion (both ∝ 1/Re, ~10⁻³ of U²/L at Re_L = 10³), which is why ∂p/∂y = 0 (9.10) is safe."
2. **Cell 123, D06 step 6 (θ/δ = ∫f′(1−f′)dη = 2f″(0)).** (a) Integration by parts is used here for the first time, but its primer (P218) only comes at cell 399 and no earlier chapter primes it (not in `knowledge/primers.md`); "integration by parts" is listed under *Tools* without a source. (b) The step packs four moves into one line and names an undefined quantity ("E′ = 0"). Add a 📎 primer on integration by parts before D06 (2-line demo), and split step 6 into: 6a write f f″ = f (f′−1)′ (because (f′−1)′ = f″); 6b integrate (9.27) from 0 to η: f″(η) − f″(0) = −½∫₀^η f (f′−1)′ dη; 6c by parts: ∫f (f′−1)′ = f(f′−1) − ∫f′(f′−1), so f″(0) = ½I_θ(η) − ½f(1−f′) + f″(η); 6d let η → ∞ (both boundary terms decay like a Gaussian) ⇒ θ/δ = 2f″(0). Either delete "checked at η = 1, 2, 4 to 1e-11" or add the cell that does it.
3. **Cell 402, D19 "Check it".** "Units: u³L² = m⁴/s³" is wrong: u³L² = (m/s)³·m² = **m⁵/s³** (the derivation-review item R8; cell 419 already says m⁵/s³, so the notebook contradicts itself). Correct the line.
4. **Cell 169 (What would change if… n = −0.1).** "only if the profile may have reverse flow (amber) and then f′ → 1 no longer holds…" suggests reversed-flow solutions exist at n = −0.1 — the overstatement R11 removed in cell 162 (the amber branch lives only on −0.0904 < n < 0 and joins the attached one at the fold). Replace with: "…nothing: below the fold there is no attached bounded solution with f′ → 1; the amber reversed branch also ends at the fold (it exists only for −0.0904 < n < 0). A marched layer at n = −0.1 separates (cell 260)."
5. **Cell 232 (primer P211 code).** The comment says `# 0.0958 twice`; the printed output is 0.1104 twice. Fix the comment to 0.1104.
6. **Cells 271 and 287 (C09 idea box and reading notes).** "drag ∝ area between the curves" / "The shaded area between the real and the ideal C_p is the missing pressure recovery = the form drag" is not true: C_D,p = ½∮(C_p − C_p,ideal) cos φ dφ, a **cos φ-weighted** area (between 82° and ≈110° the shaded region lies where the real pressure is *above* the ideal one and cos φ changes sign at 90°). Say "the shaded gap, weighted by cos φ, is the form drag" (or shade (C_p − C_p,ideal)cos φ instead).
7. **Cell 185, D08 "Check it".** "`BL.momentum_integral_residual` ≈ 1e-8" contradicts the printed residuals in cell 190 (2.9e-7, 1.9e-7, 3.6e-8 of max τ₀) and cell 191 ("about 10⁻⁷"). Quote the printed numbers.

## Should fix

- **Cell 33 (N10).** The sentence "It fails at separation and where the surface curves as sharply as the layer is thick: $\frac{\partial p}{\partial y}=0$." is garbled (the colon attaches (9.10) to the place where it *fails*). Rewrite as "…as the layer is thick; elsewhere ∂p/∂y = 0 (9.10)."
- **Unexplained terms / idioms** (see the table below): von Mises variables (cell 46), lambda default-argument capture `lambda y_, U0_=U0_, d0_=d0_:` (cell 260), "Lamb's asymptote" (cell 327); add one sentence each. Add broadcasting P77 (`xs[None, :]`, first in cell 47) to the C01 recap list (cell 13).
- **Numbers in prose that no cell computes:** cell 18 "stretched by a factor 15" (the panel aspect is ≈ 10; compute or drop); cell 136 "η_max = 5 … off by 4×10⁻³, at 6 by 5×10⁻⁴" (add to the truncation loop of cell 130); D10 check "reproduces them to 2e-8" (cell 215). Cell 132 prints τ₀ as `[0.0049 0.0022 0.0015]` while cell 133 quotes 1.543×10⁻³ Pa — print with `.4g`.
- **"What does the code above do?" missing after 11 from-scratch cells:** 79, 134, 166, 192, 230, 258, 285, 314, 385, 424, 444. The line comments are good, but each needs the numbered block (what is built by hand, what it is compared with, the assert).
- **Slider figures without reading notes:** cells 145, 172, 264, 288 (one "What you see / How to read it" line each; 84, 217 and 388 have "What to try" and are fine).
- **Symbol clash σ (C10).** D14 uses σ = sinh(πb/a)/cosh²(πb/a) and λ, μ for eigenvalues; cells 316, 318 and the axes of the cell-315 figure (panel c "Re σ", "Im σ") use σ for the eigenvalue/growth rate ("σ = π/4", "σ = 0.48"). Rename to λ (or "growth rate") there.
- **Cell 308, D14 step 10** mentions "Part F went straight to…" — an internal design-document reference the reader cannot see; delete the parenthesis.
- **Cell 281 (D13 trap)** "(ch06 measured from downstream, D09 there)" is ambiguous now that this chapter also has a D09 (Thwaites); write "Ch. 6's derivation D09".
- **Cell 69, D03 check** "(equal to (9.39), …'s v_∞ of C06)" reads badly; write "the same v_∞ that appears in (9.39) … (C06)".
- **Cell 215, D10 check**: "so a claim of '< 5 %' holds for H, not for l" — whose claim? Name it (the book's accuracy statement) or drop it; the paragraph is also too long for a Check.
- **Cell 453 (What should have clicked, C11)**: "A turbulent layer separates later (82° → 125°) and the drag falls" re-attributes the drop to later separation, contradicting C09's measured lesson (at fixed C_b later separation *raises* the model drag). Add "…and the wake pressure rises, so the drag falls".
- **D19 ★★★ check (cell 403)** tests the moves of steps 6 and 11 on a test ψ, not the result (9.80). Cell 425 already computes the invariant at four stations (`w_inv`) but only plots it; print and assert its constancy (and the falling ρ∫u²dy) and point D19's check to that cell.
- **D22 (cell 439)** step 5 writes ρu_e²/R as "1000×0.04/0.04" (0.04 is both u_e² and R) — write 1000×0.2²/0.04; step 6 cites the Cartesian (6.2) for an axisymmetric flow — say "mass conservation (the axisymmetric form of (6.2))".
- **Figures:** cell 48 (b) legend sits on the bars; cell 207 (a) "Hiemenz" arrow and "−0.0681 (exact)" label cover data points; cell 425 (c) legend covers the bars; cell 445 (a) the axis is not drawn and the tea leaf sits away from the centre (the whole point is that leaves collect at the axis); cells 168 and 262 say "dots: inflections" but the n = 0 inflection at the wall is not dotted (`if ie:` skips 0.0 — use `if ie is not None`); cell 17 left panel: the "U = 1 m/s" label overlaps the streamlines.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / tool | First used | Explained |
|---|---|---|
| integration by parts | 123 (D06 step 6) | **missing at use** (P218 only at cell 399; not in earlier primers) |
| "E" in "E′ = 0" | 123 (D06 step 6) | **missing** (never defined) |
| von Mises variables (x, ψ) | 46 | **missing** |
| lambda default-argument capture (`U0_=U0_`) | 260 | **missing** |
| Lamb's asymptote | 327 | **missing** (named only) |
| numpy broadcasting `a[None, :]` | 47 | Ch. 2 P77 exists; not recalled in cell 7/13 lists (add) |
| "Part F" | 308 | internal reference — remove |
| δ₉₉ | 18 | forward pointer to C02 (cell 68) ✓ |
| `observed_order` | 34 | 35 ✓ |
| backward Euler | 46 | Ch. 8 P192 ✓ |
| `CubicSpline` | 134 | comment gloss ✓ (acceptable) |
| parabolic/elliptic, marching | 44 | primer P201, cell 42 ✓ |
| Töpfer scaling, shooting vs BVP | 117 | P207, P208 (cells 118, 120) ✓ |
| fold / saddle-node | 162 | P209, cell 157 ✓ |
| integrating factor | 211 | P210, cell 209 ✓ |
| ∫sin⁵ by c = cos φ | 233 | P211, cell 231 ✓ |
| inflection point | 247 | P212, cell 248 ✓ |
| cotangent lattice sums, linear stability | 308 | P213, P214 (cells 304, 306) ✓ |
| sech, arccosh | 356 | pointer, P215 cell 370 ✓ |
| total-derivative move | 374 | P216, cell 372 ✓ |
| momentum / mass flux | 357 | P217, cell 359 ✓ |
| partial fractions, `apart`, implicit inversion | 414 | P219, cell 412 ✓ |
| radial force balance in a swirl | 435 | P220, cell 437 ✓ |
| six Reynolds numbers | 19 | P200 ✓ |

## Derivation audit (D id · steps · every step follows? · why/in-words ok · check ok · verdict)

| D | ★ | steps | every step follows? | why / in words | check | verdict |
|---|---|---|---|---|---|---|
| D01 | ★★ | 14 | yes (verified the Re powers of steps 6–12 by hand; (9.8) sign correct) | ok | units, numbers, `bl_nondim_sympy` ✓ | OK |
| D02 | ★ | 5 | yes | ok | ✓ | OK |
| D03 | ★ | 6 | yes (step 6 continuity with fixed h) | ok | ✓ (wording of the (9.39) clause: Should) | OK |
| D04 | ★★ | 9 | yes | ok | 3.086e−3 N/m both sides ✓ | OK |
| D05 | ★★ | 13 | yes | ok | sympy residual 0 ✓ | OK |
| D06 | ★★ | 10 | **step 6 no** — needs f f″ = f(f′−1)′ and a by-parts move not yet taught; "E" undefined | step 6 overloaded | θ = 2f″(0) ✓ numerically | **FAIL** (Must 2) |
| D07 | ★★ | 12 | yes (checked steps 3–12 by hand) | ok | sympy ✓, δ = 1.22 mm ✓ | OK |
| D08 | ★★ | 11 | yes (steps 4–7 checked) | ok | residual quoted ≈1e−8 vs printed 1e−7 | OK after Must 7 |
| D09 | ★★ | 10 | yes (steps 2–10 checked; diffuser λ closed form re-derived) | ok | ✓ | OK |
| D10 | ★★ | 9 | yes (L values 0.441, ≈0, 0.818 re-computed) | ok | "2e-8" not shown; "< 5 %" unclear | OK (Should) |
| D11 | ★★ | 8 | yes (λ(30°) = 0.0720 re-computed) | ok | ✓ | OK |
| D12 | ★ | 6 | yes | ok | f‴(0) = −n ✓ | OK |
| D13 | ★ | 7 | yes (normal n_x = −cos φ, ½∮ factor checked) | ok | 2.667, 0.884, 0.578, 2.59 re-computed ✓ | OK |
| D14 | ★★★ | 16 | yes (steps 5, 7–11 checked: π²/3a², −π²/6a², P real, Q imaginary) | ok (σ symbol clash outside D14: Should) | sympy rebuilds steps 8/10 and tests 13, 15, 16 ✓; numeric Jacobian ✓ | OK |
| D15 | ★★ | 7 | yes | ok | J = 1.00000000 at 4 stations ✓ | OK |
| D16 | ★★★ | 14 | yes (steps 8–14 checked by hand) | ok | sympy residual = −ODE ✓ | OK |
| D17 | ★★ | 11 | yes ((2√6)³/C = 36 exactly) | ok | Bickley coefficients ✓ | OK |
| D18 | ★ | 5 | yes | ok | 1.515 vs 1.160 mm ✓ | OK |
| D19 | ★★★ | 13 | yes (T₁, T₂ bookkeeping and the factor 2 checked) | ok | **units line wrong (m⁴/s³)**; sympy tests moves 6 and 11 only | **FAIL** (Must 3; Should: assert the invariant) |
| D20 | ★★★ | 12 | yes (steps 4–11 checked) | ok | sympy: correct 0, printed −3f‴ ✓ | OK |
| D21 | ★★★ | 15 | yes (4.2896 tail constant, f″(0) = f∞³/72, C → C/λ² checked) | ok | sympy steps 4, 6, 10, 12, 14 ✓; IVP ✓ | OK |
| D22 | ★★ | 7 | yes | step 5 arithmetic ambiguous, step 6 continuity form | ✓ | OK (Should) |

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **Printed-vs-correct boxes** for every book slip (R1, R2, R3, R4, R6, R7, R17, Magnus, chapter cross-reference), each with a code option that reproduces the printed form and fails (`printed=True`, `coeff=1.0`, `printed_9_7=True`). The distinction slip vs trap (R5 "one side") is handled cleanly.
- **Example 9.2 order**: the book's criterion (−0.090, x/L = 0.1583) first, then the exact Falkner–Skan fold (−0.068148, computed, not typed), both closed form and marched, plus a from-scratch sign-change search. Good template for "two criteria, both reported".
- **Measured accuracy instead of quoted percentages** (Thwaites θ error on the exact family, cell 220; Pohlhausen within 3 %, cell 193) — the "make an approximation measurable" lesson applied again.
- **★★★ sympy cells that rebuild the construction** (D14 builds the 4×4 matrix from the conjugated lattice sums and tests the factorisation; D16/D20 show the residual *equals minus the ODE*; D20 shows the printed ODE leaves −3f‴).
- **"A bracket, not a number"** for marching near the fold (cells 259–261) — an honest way to present an inlet-dependent result.
- **Every qualitative input labelled** (rounded thresholds, illustrative C_b, schematic C_D curve "no dataset", illustrative teacup profile).
- The (9.8)-solved-for-∂p*/∂y* tiny example with easy numbers (cell 34) is a model for checking a transcribed sign.
- Candidate lessons for the skill: (i) when a derivation uses a tool that a later primer introduces (P218 integration by parts), move the primer to the first use; (ii) check "Check it" unit lines against the notebook's own later statements (Ψ m⁵/s³); (iii) `if x:` on a float that may be 0.0 silently drops the Blasius case — use `is not None`; (iv) a shaded area between two C_p curves is not the drag unless weighted by cos φ.

## Round-1 verdict: FAIL

Seven Must-fix items, all local text/comment edits except Must 2 (one new primer and D06 step 6 split into four steps). No physics or code errors were found in the executed results; the structure, coverage, derivations and explainer integration are otherwise at PASS level.

(Round-1 verdict above; superseded by the round-2 verdict below.)

---

# Round 2 — re-review of the rebuilt notebook                   2026-09-30

coverage_check: OK (0 errors, 0 warnings) · executed notebook `outputs/ch09/executed.ipynb`, 471 cells, **0 error outputs** · sections 11/11 · CORE 14/14 with code+visual · 22 derivations (D06 now 14 steps) · 9 explainers. Cell numbers below are round-2 numbers.

## Round-1 Must-fix — status

| # | Round-1 item | Round-2 evidence | Status |
|---|---|---|---|
| 1 | C01 bars: "cross-stream pressure bar smaller still" | cell 49: "of the same small order as the grey dropped x-diffusion bar (both ∝ 1/Re, ≈10⁻³ of U²/L at Re_L = 10³)"; figure legend moved off the bars | fixed |
| 2 | D06 step 6 compressed, integration by parts unexplained, "E′ = 0" | new primer **P218a** (cell 124) with a runnable `quad` demo (0.38177 twice) placed right before D06 (cell 126); D06 now steps 6–10: f f″ = f(f′−1)′ → integrate (9.27) 0…η → by parts → solve for f″(0) → η → ∞. I re-derived each line: step 7 f″(η) − f″(0) = −½∫f(f′−1)′ ✓; step 8 boundary term −f(1−f′), integral +I_θ ✓ (P218a sign convention matches); step 9 f″(0) = ½I_θ − ½f(1−f′) + f″(η) ✓; step 10 limits ✓. "E" and the uncomputed "η = 1, 2, 4" claim are gone; Tools line cites P218a | fixed |
| 3 | D19 units m⁴/s³ | D19 check: "(m/s)³·m² = m⁵/s³"; no m⁴/s³ left anywhere | fixed |
| 4 | cell 169 n = −0.1 reverse-flow overstatement | now: "no attached bounded solution … the amber reversed branch also ends at the fold (exists only for −0.0904 < n < 0); a layer marched into n = −0.1 separates" | fixed |
| 5 | P211 comment 0.0958 | `# 0.1104 twice`, output 0.1104 twice | fixed |
| 6 | "drag ∝ area between the curves" | idea box: "drag = ½∮(C_p − C_p,ideal) cos φ dφ"; reading notes say "weighted by cos φ … not its plain area"; panel title "shaded: gap to ideal flow" | fixed — **but the new sentence contains a wrong number (Must R2-1)** |
| 7 | D08 check "≈ 1e-8" | "about 10⁻⁷ … prints 2.9e-7, 1.9e-7, 3.6e-8" — matches the output | fixed |

## Round-1 unexplained-first-use rows — status

integration by parts → P218a before D06 ✓ · "E" → removed ✓ · von Mises variables → explained where used (cell 46: ψ as the cross-stream coordinate, removes v, grid follows streamlines) ✓ · lambda default-argument capture → explained in the note after the bracket cell ✓ · Lamb's asymptote → named with its formula C_D ≈ 8π/[Re(2.002 − ln Re)] (checked: 0.5 − γ_E + ln 8 = 2.002) ✓ · broadcasting `xs[None, :]` → added to the tools list (cell 7) ✓ · "Part F" → removed ✓. P218 (cell ~399 of round 1) now opens with a recap line pointing to P218a ✓.

## Round-1 Should-fix — status

Fixed: N10 garbled sentence; cell 18 "about ten times"; truncation study now prints η_max = 5, 6 (4.1e-3, 5.1e-4) matching the prose; τ₀ printed to 4 significant figures (0.001543); every from-scratch cell except one now has "What does the code above do?"; the four slider figures have reading notes; D13 trap says "Ch. 6's derivation D09"; D03 check wording; D10 check rewritten (no unattributed "< 5 %", numbers are the printed rows); C11 summary now credits the higher wake pressure; D19 invariant printed and asserted at four stations (constant to 1e-6, ρ∫u²dy ∝ x^(−1/4) asserted); D22 step 5 uses 0.2², step 6 says "axisymmetric form of (6.2)"; figures: C01 legend, Thwaites closure labels, wall-jet bar legend, teacup schematic (axis drawn, leaf at the centre, mirrored loop), `is not None` for the n = 0 inflection.

## Must fix (round 2)

1. **Cell 299 (reading notes of the C09 C_p figure, new text).** "between 82° and about 110° the real pressure sits *above* the ideal one" — wrong crossing angle (my round-1 wording carried the same slip). With C_b = −1.2, 1 − 4 sin²φ = −1.2 gives sin²φ = 0.55, φ ≈ **132°**; the wake plateau lies above the ideal curve from 82° to ≈ 132° and below it from 132° to 180° (visible in the cell-298 figure: the dashed ideal curve crosses the rose plateau near 130°). Change "about 110°" to "about 132° (where 1 − 4 sin²φ = C_b = −1.2)".

## Should fix (round 2, new or remaining)

- **C10 symbol clash moved, not removed (cells 330, 332).** The reading notes now write the eigenvalue quartet as "s = ±γ ± iω (growth rate γ …)" and quote "growth rate γ = 0.48" (b/a = 0.15) and "γ = π/4". In D14, γ = ½ − sech²(πb/a) is the *dimensionless coefficient* (γ = −0.307 at b/a = 0.15; the growth rate is (πΓ/2a²)|γ| = 0.48). Use D14's own names: eigenvalue λ = (πΓ/2a²)(±γ ± iσ), growth rate Re λ; write "growth rate Re λ = 0.48", "Re λ = π/4 (non-staggered)", and label the panel-(c) axes "Re λ", "Im λ".
- **Cell 47** (the C01 from-scratch terms-of-(9.1) cell) is still followed directly by the figure cell with no "What does the code above do?" (it was also missing in round 1; I missed it then).

## Quick sweep of new/changed cells

No other problems found: P218a's demo and sign convention are consistent with D06 step 8; the new D06 steps each carry one move with *why* and *in words*; Lamb's formula is correct; the von Mises and default-argument glosses are accurate; the D19 assert tests the derived result (9.80) and the x^(−1/4) decay; the teacup schematic now shows the right loop. The 0.2805, 103.11°, Example 9.2 (0.15829 / 0.12569) and all Blasius/jet numbers are unchanged.

## Round-2 verdict: FAIL

One Must-fix remains (R2-1, a single wrong angle in new prose: "about 110°" → "about 132°"). All seven round-1 Must-fix items and all round-1 unexplained-first-use rows are resolved. Once R2-1 is corrected, the notebook meets PASS; the γ/λ naming in C10 and the note after cell 47 are recommended but not blocking.

---

# Round 3 — final check                                          2026-09-30

coverage_check: OK (0 errors, 0 warnings) · `outputs/ch09/executed.ipynb` 472 cells, 0 error outputs.

| Round-2 item | Evidence in the executed notebook | Status |
|---|---|---|
| Must R2-1: crossing angle in the C09 C_p reading notes | cell ~299 "How to read it": 1 − 4 sin²φ = C_b = −1.2 ⇒ sin²φ = 0.55, φ ≈ 132°; real above ideal from 82° to 132°, below beyond; cos φ changes sign at 90° so parts of the gap push forward and parts backward (checked: arcsin √0.55 = 47.9°, 180 − 47.9 = 132.1°) | fixed |
| Should: cell 47 lacked "What does the code above do?" | cell 48: four numbered points (broadcast grid, `np.gradient` derivatives, residual assert of (9.9), dropped/kept ∝ 1/Re_x assert) — all match the code | fixed |
| Should: C10 γ/σ clash | cells 331, 333 and the spectrum axes now use D14's eigenvalue λ (growth Re λ, frequency Im λ), with a note that it is unrelated to Thwaites' λ; numbers consistent (Re λ = 0.48 at b/a = 0.15 ⇒ e-folding ≈ 2.1; 0.54 at 0.5; π/4 non-staggered) | fixed |

Nothing new introduced by these edits; the rest of the notebook is unchanged from round 2 (same headline numbers).

## Verdict: PASS
