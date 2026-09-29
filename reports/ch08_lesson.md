# Chapter 8 — lesson review                                   2026-09-28

Reviewed: `outputs/ch08/executed.ipynb` (496 cells, 0 execution errors, 16 PNG figures, 5 animations, 8 plotly figures,
1 live widget, 9 explainers), built by `notebooks/build_ch08.py`; against `analysis/ch08_curation.md` and
`analysis/ch08_design.md` (Parts A, E, F). Cell numbers are 0-based indices of the executed notebook.

**coverage_check: OK (0 errors, 0 warnings)** · **embed_check: FAILED (6 errors, viz side — see S11)** · sections covered
7/7 (8.1–8.7 + S01 pointer) · CORE blocks 15/15 with code + visualization · derivations 33/33 written out (296 + 2
inserted steps) · ★★★ sympy checks 7/7 run and pass (D10 c165, D13 c195, D14 c218, D22 c328, D28 c420, D30 c449,
D31 c453; plus D27 c414) · primers P185–P199 15/15 present, each before its first use.

Verified corrections required by the caller — present and right: backflow for U > 0 iff dp/dx > 2μU/h² (D05, c80);
net flow reversed beyond 6μU/h² (N10 c83, figure c85); engine film ρ = 870, μ = 0.05, ε²Re_L = 4.35 × 10⁻³ (c172–c175);
far-field ratio ½Re_a(r/a), crossover ≈ 2a/Re_a (D32 step 7, c483–c487); ∫ω dy = +U (N53 c304, c295); 2 erf⁻¹(0.95) =
2.772 (N59 c330, c331); rear minimum −3μU/2a (D30, N88); Oseen C_D = (24/Re)(1 + 3Re/16), Re on the diameter (N99);
ice μ ~ 10¹³ Pa s (N63). **Missing or incomplete:** the Couette–Poiseuille extremum y* (M6) and the slider
recirculation criterion (M5).

## Must fix

1. **Garbled inline maths `$e^(…)$` in derivation steps (20 places).** Cells 286 (D19 steps 4, 9: `$e^(ln G)$`,
   `$e^(-ζ^2)$`), 327 (D22 step 13), 361 (D24 tools, steps 1–3, 4, 7, 8, 9, in-words of the result — 13×), 362 (D24
   traps), 366 (D25 step 5 `$e^(-2π)$`), 479 (D33 tools, step 3). KaTeX renders only "(" as the superscript, so the
   reader sees e⁽iωt) — the key exponentials of Stokes' second problem are unreadable. **Cause:** `tidy_raw_tex()`
   (`build_ch08.py` l. 697–698) first turns `e^{…}` into `$e^{…}$`, then the second `re.sub(r"\^\{([^{}$]*)\}", r"^(\1)", seg)`
   runs on the same segment and rewrites the brace *inside* the maths it just created. **Fix:** apply the `^{…}→^(…)`
   rewrite only outside `$…$` (re-split `seg` with `_PROSE_SPLIT`/`\$[^$]+\$` after the first substitution), and add to
   `self_check_prose()` an assertion that no `$…^(…$` survives.

2. **The "printed" form is displayed as the corrected equation (self-contradictory text).** `show_eqs()` appends the
   corrected equation `EQ[n]` after every first mention of a number, including where the text names the *book's* slip:
   - c7 (slips table), row (8.13b): the "The book prints" column shows the correct ∂p/∂**y** equation, then "with −(1/ρ)∂p/∂x";
   - c7, row (8.19): the "book prints" column shows the consistent $…+U_0(1-\frac yh)$, then "ending in +U₀";
   - c186 (D12 step 8 *why*): "The printed (8.19), $…+U_0\big(1-\frac yh\big)$ adds U₀ instead of U₀(1 − y/h)";
   - c189 (N34): "The printed (8.19), $…+U_0(1-\frac yh)$ adds U₀ to the Couette part";
   - c168 (D11 goal): "restoring the ν the printed (8.17a), $0\cong…+\nu\frac{\partial^2u}{\partial y^2}$ omits";
   - c164 (D10 check) and c166 (D10 traps): "the printed ∂p/∂x in (8.13b), $…\partial p/\partial y…$".
   A first-time reader cannot tell which form is the slip. **Fix:** in `show_eqs()` skip the insertion when the 30
   characters before the number contain "printed" or "book prints"/"The book prints"/"The text"; where the printed form
   matters, write it explicitly (N34/D12 step 8: printed $u\cong-\frac{h^2}{2\mu}\frac{\partial p}{\partial x}\frac yh(1-\frac yh)+U_h\frac yh+U_0$;
   D11 goal: printed $0\cong-\frac1\rho\frac{\partial p}{\partial x}+\frac{\partial^2u}{\partial y^2}$; slips-table rows:
   put only the printed form in the left column and the corrected equation in the "We use" column).

3. **Wrong radius for Re = 2000 in the pipe figure text.** c117 ("How to read it": "the amber line, a ≈ 2.5 mm here") and
   c119 ("past about 2.5 mm Re would exceed 2000"). With dp/dz = −1000 Pa/m, V = a²|dp/dz|/(8μ), Re = 2aV/ν = 2000 ⇒
   a³ = 8 × 10⁻⁹ m³, **a = 2.0 mm** (check: a = 2 mm → V = 0.5 m/s, Re = 2000); the figure itself labels 2.1 mm (grid).
   **Fix:** `build_ch08.py` l. 1563 and l. 1587 → "a ≈ 2.0 mm" (better: print `a_cross` from the figure cell and quote it).

4. **Wrong sand fall speed.** c464 "How to read it": "a 10 µm sand grain about 0.2 mm/s". `terminal_velocity(10e-6, 2650,
   1000, 1e-3)` = 3.6 × 10⁻⁴ m/s = **0.36 mm/s** (2 × 1650 × 9.81 × 10⁻¹⁰ / 9 × 10⁻³), which is what panel (c) plots.
   **Fix:** `build_ch08.py` l. 3729 → "about 0.4 mm/s".

5. **Slider recirculation criterion incomplete (only α > 1 is taught).** c224 item 2 says "no inlet recirculation (α ≤ 1)";
   c234 (live widget) prints "inlet backflow"; c235 "past α = 1 a small recirculation appears at the wide inlet" and the
   What-to-try bullet only pushes α past 1; c230/c231 sweep α down to −0.8 without saying that recirculation reappears.
   The verified rule (reports/ch08_verification.md, `slider_bearing_state` docstring, O9): **recirculation at the wide end
   iff the wide/narrow gap ratio exceeds 2, i.e. α > 1 or α < −½.** **Fix:** state it once in C07 (N35 or after D14:
   pad-frame profile u/U = −s + 3(1 − h_m/h)s(1 − s), reversed where h > 1.5h_m, h_m = 2(1 + α)h₀/(2 + α) ⇔ gap ratio > 2),
   rename "inlet backflow"/"inlet recirculation" to "recirculation at the wide end" in c224, c234 (`build_ch08.py`
   l. 2216, 2305, 2312) and add a bullet "set α = −0.7: recirculation at the wide end, now at x = 0".

6. **The Couette–Poiseuille velocity extremum is never taught.** The trap c81 contrasts backflow with "when the velocity
   maximum moves", and the E1 explainer shows u_max and its height (reviewed fix M1, `couette_poiseuille_state`), but the
   notebook never derives where the extremum is. **Fix:** add one step to D05 (or a sentence to N09): du/dy = U/h −
   (1/2μ)(dp/dx)(h − 2y) = 0 ⇒ **y* = h/2 − μU/(h·dp/dx)** (an interior extremum only when 0 < y* < h; U = 0 gives h/2);
   number: h = 1 cm, U = 1 cm/s, dp/dx = −0.5 Pa/m ⇒ y* = 7 mm; and print `u_max`, `y_umax`, `u_min`, `y_umin` from
   `couette_poiseuille_state` in c72.

7. **D14 step 15 gives a wrong reason, and two step pointers are off by one.** c217 step 15 *why*: "For α ≪ 1, (2 + α)(1 +
   αx/L)² ≈ 2; (8.19) … is only valid for small slopes anyway." The gap slope is dh/dx = αh₀/L (≈ 10⁻⁴ here), small for
   *any* α of order one, so lubrication does not require α ≪ 1; the linearisation is a separate approximation that the
   notebook itself shows costs 15.6 % at α = 0.1 and ≈ 90 % at α = 0.5 (c224, c227). Also the Assumptions line (c217) and
   the trap (c219) say α ≪ 1 enters "at step 14" — after the inserted step 12 it is step 15. **Fix:** `pf_sub("D14",
   "step15.why", …)` → "An extra approximation, not required by lubrication (the wall slope αh₀/L is tiny for any α): for
   α ≪ 1, (2 + α)(1 + αx/L)² ≈ 2 — the error is 15.6 % at α = 0.1 (N35)"; `pf_sub` the assumptions and traps to "step 15".

8. **Two derivation headings lost their equation.** c286: "#### 🧮 Derivation — Solving the ODE: and $\frac uU=…$ `D19`";
   c393: "#### 🧮 Derivation — From to the Stokes equations …`D26`". `_title_no_numbers()` (l. ~524) deletes "(8.29)" and
   "(8.39)" instead of showing them — the heading is ungrammatical and the equation it names is not shown (house rule).
   **Fix:** in `_title_no_numbers`, replace a bare number that is followed by "and"/"to" with `$EQ[n]$ (n)`, or `pf_sub`
   the two titles: "Solving the ODE: $F(\eta)=A\int_0^\eta e^{-\xi^2/4}d\xi+B$ (8.29) and …" and "From
   $\rho\mathbf u\cdot\nabla\mathbf u+\nabla p=\mu\nabla^2\mathbf u$ (8.39) to the Stokes equations …".

## Should fix

1. **C01 figure (c31, left panel) does not show "dye fills the tube".** The seeded random walk drifts to −0.28, is clipped
   and then runs flat along the lower wall from x ≈ 0.33 m — it looks like dye sticking to the wall. Use a reflecting walk
   or 5–8 streaks from different heights (or a shaded band) so the turbulent tube is visibly filled.
2. **c116 panel (b):** the "Re = 2000 / a ≈ 2.1 mm" label overlaps the legend; move the legend to upper left or the label
   above the curve.
3. **c176 panel (a):** pressure (orange) and cross-gap friction (rose) both sit at 1, so the orange line is hidden, while
   c177 says "two flat lines". Dash one of them (or offset by ×1.05 as c399 does); the "ε²Re_L = 1" text overlaps the
   legend.
4. **c434 "What you see" (b)** says the fluid-frame streamlines "arrive along the axis ahead of the moving sphere … and leave
   along the axis behind". The plotted fluid-frame Stokes streamlines are open U-shaped curves (y² ∝ r far away),
   mirror-symmetric front to back; reword to what is drawn. Same cell, (c): "closed loops of a dipole" is right.
5. **Local reminders placed after the first use** (cell 8 lists them by id, so not a Must): `special.erf` first used c54
   (P195 primer at c282; P123 covers erfc only) · `np.trapezoid` c54 (P37 reminder at c114) · `np.meshgrid` + `ax.contour`
   c207 (P76/P78 reminder at c432) · `observed_order` c210 (reminder at c299). Move each one-line reminder into the
   "Tools from earlier chapters used here" cell before c54 (c53), before c207 (c206) and before c210 (c209).
6. **"What does the code above do?" missing after non-trivial from-scratch cells:** c254 (explicit finite-volume thin
   film), c300 (hand-made Crank–Nicolson, FTCS comparison, convergence table), c459 (midpoint + Gauss–Legendre drag),
   c485 (hand finite differences of u·∇u and ν∇²u). The comments are good, but a numbered 3–4-line block says what the
   printed numbers prove (e.g. c300: "observed order 2.00 = second order in space and time; FTCS needs 5528 steps").
7. **D27 step 5 *in words*** (c413): "vorticity neither diffuses nor is advected here" is misleading — in steady creeping
   flow advection is absent and diffusion is *in balance*; say "no advection and no change in time: the diffusion of ω
   balances itself, so each Cartesian component is harmonic".
8. **D29 step 4** (c423) packs two moves (write E²g = 0 with g = f″ − 2f/r²; expand (f/r²)″) — split into two steps.
9. **★★★ checks D10, D13, D14 call cached library engines** (`lubrication_nondim_sympy`, `reynolds_equation_sympy`,
   `slider_bearing_sympy`) and assert their outputs. They do test the derived results, but the teaching-style lesson
   (ch01) asks the check to re-run the construction; add 3–5 lines per cell that build one key step in the cell (D10:
   substitute the scalings into one term and read the coefficient; D14: `sp.dsolve`/integrate step 5's ODE and compare).
10. **c438 What-to-try** says the wake appears "once r reaches a/Re"; use a/Re_a (≈ 2a/Re_a with the prefactor), as C15
    does, so the two statements agree.
11. **Explainer derivation tabs (embed_check errors):** `stokes_drag_settling.html` lacks the D30, D31 steps and
    `stokes_sphere_flow.html` lacks D28, D29, D32, D33 (listed in `viz:derivations` but not in the script). Owner: the
    viz-builder (phase 7); the notebook's pointer "Open the Derivation tab at D29 step 8" (c438) depends on it.

## Unexplained-first-use list

| term / tool | first used | explained at | status |
|---|---|---|---|
| ν as momentum diffusivity, L²/ν | c14 | C01, P185 c19 | ok |
| erf (scipy.special.erf) | c54 (entrance profiles) | P195 c282 (erfc via Ch. 4 P123 listed in c8) | late (S5) |
| np.trapezoid | c54 | P37 reminder c114 (listed c8) | late (S5) |
| brentq | c54 | P108 reminder c53 | ok |
| partial derivative, continuity | c57 | P25 reminder c56, recap c34 | ok |
| np.diag, np.r_ | c76 | gloss c75 | ok |
| cylindrical Laplacian | c100 | P186 c101 | ok |
| Euler–Cauchy ODE | c129 | P187 c130 | ok |
| np.gradient | c145 | P22 reminder c144 | ok |
| quiver | c148 | P78 (listed c8) | ok (listed) |
| anisotropic scaling, ε bookkeeping | c158 | P188 c160, reminders c162 | ok |
| Leibniz rule (moving limit) | c194 | P189 c191 | ok |
| cumulative_trapezoid | c201 | P190 c198 | ok |
| meshgrid, contour | c207 | P76/P78 reminder c432 (listed c8) | late (S5) |
| observed_order | c210 | reminder c299 (P13, listed c8) | late (S5) |
| substitution / Taylor in α | c217 | reminders c216 | ok |
| minimize_scalar, log1p | c223 | reminders c222 | ok |
| nonlinear diffusion, implicit/Picard | c243–c252 | P191 c245, P192 c248 | ok |
| Π groups (pi_groups) | c274 | Ch. 1 §1.11 named in D17 tools | ok |
| Gaussian integral, erf/erfinv | c286 | P194 c280, P195 c282 | ok |
| solve_bvp | c295 | P196 c288 | ok |
| Crank–Nicolson, solve_banded | c300 | P193 c297 | ok |
| exponent matching | c323 | P197 c319 | ok |
| complex amplitude, √i | c361 | reminders c358 | ok |
| dominant balance | c393 | P198 c390 | ok |
| ε_ijk, Schwarz, curl commuting | c413 | reminders c411 | ok |
| Stokes operator E² | c419 | P199 c416 | ok |
| traction σ·n | c448 | gloss c447 | ok |
| Gauss–Legendre | c459 | reminder c456 | ok |
| solve_ivp (many tracers) | c436 | reminder c435 | ok |

Missing entirely: none. Late (explained after first use, but listed by id in the chapter's tool list c8): erf, np.trapezoid,
meshgrid/contour, observed_order.

## Derivation audit

| D | steps | every step follows? | why / in words | check | verdict |
|---|---|---|---|---|---|
| D01 ν = μ/ρ, t ~ L²/ν | 5 | yes | ok | units + computed ratio 15.01 (c26) | ok |
| D02 v ≡ 0, (8.4a,b) | 7 | yes | ok | units; sympy (c72) | ok |
| D03 dp/dx const | 5 | yes | ok | p = x² counter-example | ok |
| D04 (8.5) | 9 | yes | ok | walls, units, number −0.0125 m/s | ok |
| D05 limits, backflow | 7 | yes | ok | units, 5 mm reversal | ok — add y* (M6) |
| D06 (8.6) | 10 | yes | ok | parity Ch. 3/Ch. 4 (c112) | ok |
| D07 τ, τ₀, Q, f | 10 | yes | ok | numbers, force balance (c115) | ok |
| D08 (8.10) | 10 | yes | ok | walls, rigid limit, 3.9 mm/s | ok |
| D09 (8.11), (8.12) | 6 | yes | ok | R₂ = 10⁶R₁ parity | ok |
| D10 (8.15), (8.16) ★★★ | 14 | yes (verified coefficients) | ok | sympy c165 passes | ok (S9) |
| D11 (8.17) | 7 | yes | ok; "printed … omits" contradictory (M2) | units | fix M2 |
| D12 (8.18), (8.19) | 8 | yes | step 8 shows the wrong form as "printed" (M2) | walls | fix M2 |
| D13 Reynolds eq. ★★★ | 12 | yes | ok | sympy c195 passes | ok (S9) |
| D14 slider ★★★ | 16 | yes (numerator identity and C₁, C₂ re-checked) | step 15 *why* wrong; step pointers off by one (M7) | sympy c218, brentq c225 | fix M7 |
| D15 thin film | 10 | yes | ok | units, volume conserved | ok |
| D16 (8.20) | 6 | yes | ok | units | ok |
| D17 (8.24)→(8.25) | 8 | yes | ok | pi_groups c274 | ok |
| D18 (8.26)–(8.28) | 9 | yes | ok | solve_bvp c295 | ok |
| D19 (8.29), (8.30) | 12 | yes (inserted step 10 fixes the skip) | ok; steps 4, 9 garbled (M1); heading (M8) | F(0), F(∞), erfc(0.5) | fix M1, M8 |
| D20 (8.31) | 5 | yes | ok | 95 % gives 2.772 | ok |
| D21 Example 8.4 | 8 | yes | ok | units | ok |
| D22 vortex sheet ★★★ | 14 | yes | ok; step 13 garbled (M1) | sympy c328 passes | fix M1 |
| D23 bead n = m = 1/5 | 9 | yes | ok | slope 0.198 (c343) | ok |
| D24 (8.35)–(8.38) | 10 | yes | ok; 13 exponentials garbled (M1) | sympy c363, CN c371 | fix M1 |
| D25 reading (8.38) | 6 | yes | ok; step 5 garbled (M1) | units | fix M1 |
| D26 (8.39)→(8.43) | 9 | yes | ok; heading (M8) | sympy coefficients c396 | fix M8 |
| D27 ∇²ω = 0 | 6 | yes | step 5 *in words* misleading (S7) | sympy curl-curl c414 | ok (S7) |
| D28 (8.44) ★★★ | 12 | yes (identity checked for generic A) | ok | sympy c420 passes | ok |
| D29 (8.48) | 10 | yes (cubic factor re-checked) | step 4 two moves (S8) | f(a), f′(a), E⁴ψ | ok (S8) |
| D30 (8.50) ★★★ | 11 | yes (all components re-derived) | ok | sympy c449 passes | ok |
| D31 D = 6πμaU ★★★ | 13 | yes (step 7 derivative re-checked) | ok | sympy c453, midpoint + Gauss c459 | ok |
| D32 far field | 7 | yes (inertia (3U²a/16r²)(8cos²θ − 4sin²θ) re-checked with sympy) | ok | ratio ÷ (Re_a r/a) → 0.5 (c483) | ok |
| D33 Oseen → Stokes | 7 | yes | ok; tools/step 3 garbled (M1) | sympy c481 | fix M1 |

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **Slips table up front (c7) plus the "as printed" ghost curves** (c226, c230) and a `form="book"` switch that a test must
  fail — the reader sees the slip, the correction and the proof side by side. (Once M2 is fixed.)
- **Three-way parity cells** (c112: this chapter's formula, Ch. 3's developing profile far downstream, Ch. 4's exact
  solution) — reuse of earlier chapters as evidence, not just as recaps.
- **Measuring an order-of-magnitude claim**: c483 prints ratio ÷ (Re_a r/a) in three directions and turns "r ~ a/Re" into
  "prefactor ½, crossover 2a/Re_a, dip near 55°". Good pattern for every "~" in later chapters.
- **Inserted derivation steps flagged in the text** ("Part F went straight to …; this step is added") — honest and
  teaches that the book skips moves.
- **Collapse measured, not just shown** (c303 spread 1e-16; c341 right vs wrong exponents 5.6e-15 vs 1.02).
- **Tractions figure** (c463 a): pressure, friction and their constant-x sum as arrows — makes "⅓ + ⅔ = uniform 3μU/2a"
  visible at a glance.
- Candidate lesson for `teaching-style`/`nbkit`: *automatic equation insertion must skip "printed/book prints" mentions*,
  and *post-processing regexes must never run inside maths they just created* (M1, M2 are both builder-level).

## Verdict (round 1): FAIL

8 Must-fix items (all small text/builder fixes; the physics, derivations and code are otherwise correct and complete).
Re-run after fixing: `notebooks/build_ch08.py`, `tools/run_notebook.py ch08 --save`, `tools/coverage_check.py ch08`, and
grep the executed notebook for `\$[^$]*\^\(` (must be empty).


---

# Round 2 — re-review                                          2026-09-29

Fresh `outputs/ch08/executed.ipynb`: 501 cells, 0 execution errors, 17 figures. **coverage_check: OK (0 errors, 0
warnings)** · **embed_check: OK** (the D28–D33 errors of round 1 are gone; no error on the notebook side) · 9 explainers
embedded once each. Cell numbers below are round-2 indices.

## Must-fix items, checked in the new cells

| # | round-1 item | round-2 evidence | status |
|---|---|---|---|
| M1 | garbled `$e^(…)$` (20×) | regex scan of all markdown cells for `$…^(…$`: 0 hits; D19, D22, D24, D25, D33 now render $e^{i\omega t}$ etc. | fixed |
| M2 | printed form shown as corrected | c7 slips table: the left column shows printed (8.13b) with ∂p/∂x, printed (8.17a) without ν, printed (8.19) ending in +U₀; the right column the corrected forms. D11 goal, D12 step 8, N34, D10 check and traps now show or describe only the printed form | fixed |
| M3 | pipe radius 2.5 mm | c116 computes `a_cross` with brentq and prints 2.00 mm (a³ = 8 × 10⁻⁹ m³ ✓); c117/c119 quote 2.00 mm; the label is clear of the legend | fixed |
| M4 | sand 0.2 mm/s | c467 prints U_t = 0.36 mm/s; "How to read it" quotes 0.36 mm/s | fixed |
| M5 | recirculation rule | new c231 (pad-frame profile) + c232 (table α = 0.1 … −0.7). Re-derived: $C_1=-Uh_m/2$, $h_m=2(1+\alpha)h_o/(2+\alpha)$, $(u-U)/U=(1-\eta)[-1+3(1-h_m/h)\eta]$, reversal next to the pad iff $h>\tfrac32h_m$ ⇔ α > 1 (at x = L) or α < −½ (at x = 0) ⇔ gap ratio > 2 — correct; the printed verdicts flip at 0.99/1.01 and −0.45/−0.55 as they must. "inlet backflow" renamed "recirculation at the wide end" in c225, the live widget and the E3 What-to-try (new α = −0.7 bullet) | fixed |
| M6 | extremum y* | D05 now 9 steps: step 8 $\frac Uh-\frac1{2\mu}\frac{dp}{dx}(h-2y^*)=0$, step 9 $y^*=\frac h2-\frac{\mu U}{h\,dp/dx}$ (re-derived ✓; "interior only for 0 < y* < h" stated); c72 prints y* = 7.00 mm by hand = library y_umax 7.00 mm, u_max = 1.225 cm/s (0.007 + 250 × 0.007 × 0.003 = 0.01225 m/s ✓) | fixed |
| M7 | D14 step 15 why + pointers | step 15 *why*: "an extra approximation, not required by lubrication (the wall slope αh₀/L is tiny for any α of order one) … the linear load is 15.6 % high at α = 0.1"; Assumptions and trap say step 15 | fixed |
| M8 | lost heading equations | D19 "Solving the ODE: $F(\eta)=A\int_0^\eta e^{-\xi^2/4}d\xi+B$ and …"; D26 "From $\rho\mathbf u\cdot\nabla\mathbf u+\nabla p=\mu\nabla^2\mathbf u$ to …"; D07, D10, D12, D25, D33 headings also read as sentences with their equations | fixed |

## Should-fix items from round 1

S1 dye sketch now fills the tube (c31, several streaks) ✓ · S2 label/legend overlap (c116) ✓ · S3 pressure line visible
(c176) ✓ · S4 fluid-frame streamlines described as open U-shaped curves, y² ∝ r (c438) ✓ · S5 reminders for erf and
np.trapezoid (c53), meshgrid/contour and observed_order (c206) now precede first use ✓ · S6 "What does the code above do?"
after c257, c303, c462, c489 ✓ · S7 D27 step 5 in words ✓ · S8 D29 split into 11 steps (step 4 applies E² to g, step 5
expands in f — re-checked, correct) ✓ · S10 "a/Re" → "order a/Re_a (about 2a/Re_a, C15)" ✓ · S9 (★★★ checks via library
engines) left as is — acceptable · S11 explainer tabs: embed_check now OK.

## New findings from the round-2 spot checks (Should fix; none blocking)

1. **Stale step pointer after the D29 split:** c430 code comment "D29 steps 6 and 8–10" — now steps 7 (roots) and 9–11
   (constants).
2. **D29 notebook vs explainer:** the notebook has 11 steps; `viz/ch08/stokes_sphere_flow.html` still has 10 (its comment
   "D29 10"). The c441 pointer now reads "the far-field step of D29 (A = 0)", so nothing breaks, but the viz-builder should
   mirror the split (same steps in both places).
3. **D05 step 9 *in words*:** "the wall drag shifts the Poiseuille peak toward the moving wall" is true for a favourable
   gradient (a maximum); for dp/dx > 0 the extremum is a minimum that moves toward the fixed wall. Suggest: "the wall drag
   moves the extremum off mid-gap — a maximum toward the moving wall (dp/dx < 0), a minimum toward the fixed wall
   (dp/dx > 0)".

No fix introduced a new error: the numbers printed by the new cells match the text (y* 7.00 mm, a 2.00 mm, U_t 0.36 mm/s,
recirculation table), the slips-table columns are consistent, and all eight sympy check cells still run and pass.

## Verdict: PASS
