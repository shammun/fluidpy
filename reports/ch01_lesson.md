# Chapter 1 — lesson review (round 2)                          2026-09-13

coverage_check: **OK** (0 errors, 0 warnings) · sections covered **11/11** (§1.1–§1.11, plus S01/S02 pointer lines) ·
CORE (A) blocks **15/15** with code + visual · derivations **12/12** present (D05, D10, D11, D14, D18, D19 ★★★, D20,
D28 ★★★, D34–D37) · explainers **5/5**, each embedded exactly once (E1 cell 94 in C06, E2 cell 133 in C12, E3 cell 314 at
the end of C45, E4 cell 415 at the end of C55, E5 cell 490 in C69) · the executed notebook has 495 cells (134 code) and 0
error outputs · 1049 of 1065 code lines are commented (the 16 bare lines are setup boilerplate, argument-continuation
lines, and two loop-body lines in cell 463).

```
$ .venv/Scripts/python.exe tools/coverage_check.py ch01 --nb outputs/ch01/executed.ipynb
COVERAGE OK for ch01 (0 errors, 0 warnings)
```

I used the same tiered-depth standard as round 1: A blocks get the full checklist, B notes are checked for being clear,
correct and self-contained, and C notes for one sentence plus a pointer. Cell numbers below refer to the new executed
notebook (495 cells), not the 485-cell notebook of round 1.

### Round-1 Must-fix status

| # | Round-1 item | Status | Evidence (new cell numbers) |
|---|---|---|---|
| 1 | ∂, ∂², ∇ used before explained | **fixed** | P25 (cell 101) covers ∂f/∂y, ∂²f/∂y² as curvature, and ∇ with "down the gradient", plus a 3-line numpy demo (102). It comes before the C12 idea table (103), C08 (104) and FTCS (119). C20 recaps it in one line (152). |
| 2 | P45/P44 after D18; Euler, cosh; N² = 0 double root | **fixed** | P45 (320) gives i, Euler's formula, cosh and cos(iσt) = cosh σt, checked in numpy (321). P44 (322) now covers the double root (ζ = A + Bt). Both come before D18 (324). D18 is now 14 steps: step 13 handles two distinct roots, and the new step 14 handles N² = 0 with ζ′(0) = B = 0. |
| 3 | ρ_θ used before defined | **fixed** | C60 (351) glosses potential density with a forward pointer to C58. |
| 4 | D20 "Which p_o?" contradicted the book | **fixed** | Cell 390 says the book takes p_o = p(0), the sea-level pressure (≈ 100 kPa), which I checked against chapters/ch01.txt at printed p. 20. Meteorology and `P_REF` use exactly 1000 hPa, and cell 390 explains why θ = 287.07 K at the ground. The D20 assumptions (391) and the trap callout (392) are reworded to "one fixed constant; not p, not a varying surface pressure". |
| 5 | D28: groups for non-repeating variables never shown to exist | **fixed** | D28 (462) is now 13 steps. The new step 9 solves A_rep a_j = −A_{·j}, which has a unique solution because det ≠ 0 (step 8). Each resulting group contains its own q_j, so the groups are independent, and n − r of them in an (n − r)-dimensional null space form a basis. The sympy cell (463) builds the four groups this way and checks A·k = 0 and rank 4. |
| 6 | D19 sympy check: implicit-function rule and `.subs().doit()` | **fixed** | Cell 373 writes out ds = (∂s/∂T)_p dT + (∂s/∂p)_T dp = 0 ⇒ dT/dp = −(∂s/∂p)_T/(∂s/∂T)_p (P49), and explains that `.doit()` evaluates the derivatives once the concrete G has been substituted. The output is 0, 0, −g/c_p. |
| 7 | Prose/output mismatches | **fixed** (all four) | P21 comment (120): "≈ −1.5e-05 = −cos(0.5)·Δ²/6", and the cell prints −1.46e-05. The 80 km density ratio is now computed (89 prints 78 007), and the prose says "tens of thousands". C40 now has an 84 kPa isobar (277), and the dot sits on it (image checked). C64 (425) draws zero exponents as stubs, a dashed outline and the label "L: 0 ≠ −1", and cell 426 describes exactly that (image checked). |
| 8 | C23 copied wording | **fixed** | Cell 193 now reads "a lump of matter we choose to follow whose molecules stay the same: energy may cross its boundary as heat or as work, but matter may not". A new 7-word overlap scan of all markdown cells against chapters/ch01.txt finds only one run, a string of maths symbols ("e = e(T) and h = h(T)") in D14. |

## Must fix (numbered: cell number · what · the concrete change)

Round 1 counted a prose number that disagreed with fluidpy as a Must fix. For consistency, these three are Must fixes
too. Each is a one-line edit.

1. **Cell 95 · the altitude is wrong.** "the gas were ten times thinner (about 40 km up)". USSA-1976 in fluidpy gives
   ρ(0)/ρ(z) = 10.2 at 18 km but 318 at 40 km. **Change:** "about 18 km up". The rest of the sentence still holds: the mean
   free path grows ≈ 9× there.
2. **Cell 350 · the explanation contradicts the code and cell 328.** Item 2 says the linear and unlinearised solutions
   "separate only for N² < 0, where the runaway leaves the range of the small-displacement steps". But the live cell (349)
   passes an incompressible parcel (`drho_a_dz = 0.0`) and a linear environment. As cell 328 itself states, the two
   equations are then *identical*. The only difference is that the unlinearised run stops at 2 km.
   **Change:** reword item 2 to "with an incompressible parcel the two equations coincide; the dashed curve only ends
   where the run is stopped at 2 km". Better still, give the live cell a parcel gradient (for example N21's −ρg/c², or
   −0.005 kg m⁻⁴ as in cell 327), so that the live comparison shows the degradation cell 328 describes.
3. **Cells 480 and 462 · step references left stale by the D28 renumbering.** Cell 480 says "That invariance is what step
   11 of D28 uses". But step 11 is now "choose the scales"; the unit-invariance of a Π value is used in **step 12**. The D28
   assumption line (462) says "The variable list is complete (step 12: nothing else can enter f)", but completeness is now
   invoked in **step 13**. **Change:** "step 12" in cell 480. In 462, write "complete (step 13)" and "dimensionally
   homogeneous (steps 12–13)".

I also checked every other "Dxx step N" reference by script (cells 53, 54, 274, 276, 308, 327, 349, 381, 397) and every
primer/block pairing: all are correct. The D18 references to "steps 13–14" were updated properly.

## Should fix

1. **Cell 111, C12 maths step 6** says the h²/D clock is confirmed because the animation "runs identically for any gap".
   The animation only runs h = 1 mm. Say "its clock is t/(h²/ν); a wider gap gives the same frames on that clock (cell 129)".
2. **Cell 322 (P44)** names sinh σt, but P45 defines only cosh. Add sinh x = ½(eˣ − e⁻ˣ) to P45. In cell 321, comment that
   `1j` is Python's imaginary unit i.
3. **D28 steps 6 and 9 (cell 462)** cite C71 and C70, which only appear after the derivation (cells 465, 467). Write
   "(worked with numbers in C70 below)" so the forward pointer is explicit. The step itself is self-contained.
4. **Cell 195** cites "working model P33" one cell before P33 (196). Swap the two primers or drop the ID from the comment.
5. **Figures (images checked):**
   - Book map (cell 8): the full titles are back, but arrows still cross the labels "Vorticity Dynamics", "Turbulence" and
     "Boundary Layers and Related Topics".
   - C20 profile (cell 183): the inset has a stray rose "α" overlapping its title "capillary rise (C22)", and point labels
     "E" and "F" that are never explained. Remove them or say what they mark.
   - C36 (cell 257): the teal ρ = 1.2 line runs through the rose "c → ∞" label. Move the label down or left.
   - C69 collapse (cell 481): the staggered Reynolds numbers reach ≈ 3500 with the laminar law, while the comment says
     "≈ 1 … 2000". Cap Re at 2000 or call the data "idealised laminar".
6. **Cell 463:** two loop-body lines have no comment (`k[col] = a_j[i]`, `groups_k.append(k)`).
7. **Duplicated reading lines.** Cells 84, 289, 385 and 407 still end *What does the code above do?* with a "Reading it:"
   line, and full *What you see / How to read it / What would change if* notes now follow. Drop the inline line or keep
   it to one clause.
8. **Cell 283:** "a van der Waals gas … sags at small v, where molecular attraction lowers the pressure" holds only where
   a/v² dominates. Very close to v = b, the b term makes p rise steeply. Say "at moderate v".
9. Primer IDs are still out of order (P43 in C06, P25 before P20, P45 before P44, P61/P60/P57). Round 2 marks renumbering
   as optional, so I accept this as it is.

### Round-1 Should-fix items reported done: verification

| Round-1 Should | Status | Where |
|---|---|---|
| 1 reading notes on six plotly/3-D figures | confirmed | cells 85, 187, 283, 290, 386, 408 (duplicate inline lines remain, Should 7) |
| 2 D05 step 3 "≈" | confirmed | cell 159 step 3 |
| 3 D37 step 7 split | confirmed | cell 166: step 7 projection, step 8 one column, step 9 integral (9 steps) |
| 4 D28 step 1 reason | confirmed | "We *look for* … steps 9–13 show that nothing more general is ever needed" |
| 5 derivation headings without "Eq." | confirmed | cells 51, 62, 166, 402 |
| 6 exponent-rules primer into C06 | confirmed | P43 at cell 59, before D35 (62); recap line at 301 |
| 7 glosses (Pr, boundary layer, thermocline, inversion, singular values) | confirmed | 132, 75/98, 318, 318, P56 reworded at 444 |
| 8 C43 γ for CO₂/H₂O | confirmed | cell 296: linear vs non-linear rigid values, real CO₂ ≈ 1.3 from vibrations |
| 9 C12 h²/D reason | confirmed (wording, Should 1) | cell 111 step 6 |
| 10 C50 linear-vs-unlinearised | confirmed in 327–328 (gap/amplitude 1.8e-5 → 3.7e-3 as ζ₀ grows); **live cell 350 contradicts it** (Must 2) | |
| 11 C36 strengthened | confirmed | cell 251: dimensional argument + piston mass/momentum plausibility (I checked the algebra: Δρ ≈ ρu/c, Δp = ρcu ⇒ c² = ∂p/∂ρ); figure 257 shows the 1 km crossing times (2.94 s, 0.67 s) |
| 12 symbol n callout | confirmed | cell 268; D11 start (270) points to it |
| 13 cosh 3 ≈ ½e³ | confirmed | cell 338 step 4 |
| 14 moist rate harmonised | confirmed | "about 4–7 K/km" in D19 (372) and cell 416 |
| 15 figure polish | mostly confirmed | C69 legend ✓, C45 plain ticks ✓, C64 stubs + text ✓, book-map truncation ✓ (overlaps remain, Should 5) |
| 16 primer renumbering | skipped (optional) | — |

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / tool | First used | Explained at | Status |
|---|---|---|---|
| partial derivative ∂/∂t, ∂/∂y | 101 (P25) | P25, cell 101 | ok (round-1 Must 1 fixed) |
| second partial derivative ∂²/∂y² | 101 (P25), 104 (C08) | P25 + demo, cells 101–102 | ok |
| gradient ∇ | 101 (P25), 103 | P25, cell 101 | ok |
| κ = k/ρC_p, C_p | 103 | glossed in the same cell | ok |
| i, Euler's formula | 320 (P45) | P45 + demo, 320–321 | ok (round-1 Must 2 fixed) |
| cosh | 320 | P45, cell 320 | ok |
| sinh | 322 (P44) | — | missing (Should 2; named only, never used) |
| Python `1j` | 321 | — | uncommented idiom (Should 2) |
| double root of ζ″ = 0 | 322, 324 step 14 | P44, cell 322 | ok |
| potential density ρ_θ | 351 (C60), 353 (1.35) | gloss in C60, cell 351; full in C58, 409 | ok (round-1 Must 3 fixed) |
| implicit-function rule | 373 | comments in cell 373 (via P49) | ok (round-1 Must 6 fixed) |
| sympy `.subs(...).doit()` | 373 | comment in cell 373 | ok |
| sympy `LUsolve`, `sp.zeros` | 463 | inline comments | ok |
| exponent rule (ab)^k | 62 (D35) | P43, cell 59 | ok |
| boundary layer | 75 | gloss in cell 75 | ok |
| Prandtl number | 132 | gloss in cell 132 | ok |
| thermocline, inversion | 318 | gloss in cell 318 | ok |
| singular values | — | removed from P56 (444) | ok |
| C70 / C71 cited inside D28 | 462 steps 6, 9 | 465, 467 | forward pointer (Should 3; step is self-contained) |
| primer P33 cited | 195 | 196 | one cell late (Should 4) |

Everything else checked was explained at or before first use: all 61 primers, recaps R01–R03, P34 before (1.10), P26
before D05, P50/P51 before D19, P58/P59/P61 before D28, P52 before D36, and P40 before the first sympy cell. No
unexplained first use rises to Must-fix level.

## Derivation audit (D id · steps · every step follows? · why/in-words ok · check ok · verdict)

I worked every step again by hand, including the new D18 step 14, the D28 step 9, the D37 split and the D19/D28 sympy
cells (outputs 0, 0, −g/c_p; 3 4 4, True ×4, "True 4"). The equations are unchanged from round 1, where they were compared
with the rendered pages. The book's p_o definition was re-checked in chapters/ch01.txt (printed p. 20).

| D | Steps | Every step follows? | Why / in words | Check | Verdict |
|---|---|---|---|---|---|
| D34 kinetic pressure | 8 | yes | ok | units ✓, 1.013e5 Pa ✓, sampled 1.002 / 1.012 ✓ | PASS |
| D35 density noise | 5 | yes; P43 is now before it | ok | dimensionless ✓, limits ✓, 6.3e-6 ✓, from scratch ✓ | PASS |
| D05 hydrostatics | 9 | yes; step 3 now "≈" with the dz → 0 reason | ok | units ✓, ρ → 0 ✓, 10 m ≈ 1 atm ✓ | PASS |
| D37 buoyancy | 9 | yes; the projection has its own step (7) | ok | units ✓, special cases ✓, six-face integral 9.807 N ✓ | PASS |
| D10 Gibbs | 7 | yes | ok | units ✓, 696.4 both forms ✓, sympy residuals 0 ✓ | PASS |
| D11 p = ρRT | 7 | yes | ok; the n-meaning callout is linked | units ✓, 1.225 ✓, constants chain ✓ | PASS |
| D14 isentropic law | 11 | yes | ok | units ✓, γ → 1 ✓, 351.3 K ✓, Euler 199 999.3 Pa ✓ | PASS |
| D18 parcel equation, N² | 14 | yes: P45/P44 come first; step 13 (distinct roots, a = b from ζ′(0) = 0) and step 14 (double root, B = 0) both follow | ok | units ✓, neutral limits ✓, 9.57e-5 s⁻², 642 s ✓; unlinearised gap ∝ ζ₀ relative ✓ | **PASS** (was FAIL) |
| D19 lapse rate ★★★ | 16 | yes | ok | sympy tests step 12 and (1.30) for arbitrary G, with the rule explained; perfect gas → −g/c_p ✓ | **PASS** (check cell fixed) |
| D20 potential temperature | 5 | yes | ok; p_o now consistent with the book and with the code | p = p_o ✓, 304.8 K ✓, from scratch ✓ | **PASS** (was FAIL) |
| D36 N² from θ | 8 | yes | ok | units ✓, neutral ✓, 3.83e-4 s⁻² ✓, USSA θ vs lapse to 1e-3 ✓ | PASS |
| D28 Π theorem ★★★ | 13 | yes: step 9 fills the missing move (unique solution by det ≠ 0; independence from the identity pattern; basis by dimension count), and steps 11–13 follow | ok; step 1's reason fixed | 7 − 3 = 4 ✓, invariance in three unit systems ✓, sympy nullspace + step-9 construction rank 4 ✓ | **PASS on the maths**; the assumption line's step pointers are stale (Must 3) |

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **Round-1 fixes went beyond the minimum.** P25 carries a curvature reading of ∂² with a numpy demo, and P45 turns
  "imaginary frequency means growth" into three printed numbers. D18's double-root step is a model of "one move, one
  reason".
- **"Which p_o?" is a reusable pattern** for any place where the book's convention and the field's (or the code's)
  convention differ. State both, give the size of the difference in numbers (287.07 vs 288.15 K), and say which one the
  code uses and why. Promote it next to the lapse-rate convention pattern in `teaching-style`.
- **The D28 sympy check now tests what is derived:** it builds the groups exactly as step 9 says and checks that they form
  a basis, not just that some null space exists. Use this as the template for ★★★ checks: "re-run the construction of the
  derivation symbolically".
- **C50's linear vs unlinearised table** (cells 327–328) is honest. It explains why the incompressible case would be a
  trivial comparison and shows the error growing with ζ₀. That makes the linearisation assumption measurable.
- **C36 now teaches rather than asserts.** A dimensional argument, then a piston mass/momentum plausibility, then 1 km
  crossing times on the figure.
- **The lapse-rate conventions remain exemplary** (C53 table, P48, `lapse_rate_stability` with the assert, both forms in
  every slider title and in the summary).
- **The tiny examples are all correct on paper.** New or re-checked numbers include the Euler error ≈ 0.02 Pa (182), Newton's
  isothermal sound speed 287.6 m/s (15 % slow), H per 10 K ≈ 0.29 km, halved-N² periods 444/889 s, Taylor's blast 4.59 K kt,
  and the Rayleigh rank 2.

**3 strongest A blocks:** C54 adiabatic lapse rate, C06 continuum hypothesis, C50 displaced parcel (now with the
double-root step and a real linearisation test).
**3 weakest A blocks:**
- C12 viscosity: the "any gap" claim overstates what its animation shows (Should 1).
- C51 N²: its maths still largely restates D18, although it is correct and now primed.
- C67 dimensional matrix: solid but plain; its heatmap is the least expressive figure of the chapter.

## Verdict: FAIL

All eight round-1 Must-fix items are fixed, and every derivation now passes on its mathematics. Three new, very small
Must-fix items remain. One prose number is wrong (cell 95: 40 km → 18 km). One explanation contradicts its own code and an
earlier cell (cell 350). Two step pointers went stale when D28 was renumbered (cells 480, 462). Each is a one-line edit in
`notebooks/build_ch01.py`. Once they are made, the notebook meets the gate; no re-review of anything else is needed.

## Orchestrator close-out (2026-09-13, after the notebook builder's final round)
The three round-2 Must-fix items were fixed by `notebook-builder` and checked directly against
`outputs/ch01/executed.ipynb` (496 cells, 0 error outputs, run 96 s) instead of a third full review (2-round limit):
1. Cell 95 now computes the ten-times-thinner height from `standard_atmosphere`: "about 17.9 km (geopotential height,
   USSA-1976)"; cell 96 refers to the printed height.
2. Live cell 350 gives the parcel its own density gradient (dρ_a/dz = −0.005 kg m⁻⁴, environment gradient from (1.29));
   cell 351 says the curves lie on top of each other on the 2 km axis and differ only where the runaway run is stopped.
3. D28 has 13 steps (cell 463, "Step 13 of 13"); cell 481 cites "step 12 of D28"; the D28 assumptions cite steps 12–13;
   D05 and D11 internal step references corrected; all 14 "Dxx step N" references checked by script.
Should-fix items from round 2 addressed except the book-map arrow overlaps (in `scripts/ch01_book_map.py`) and primer
renumbering (optional). `tools/coverage_check.py ch01 --nb outputs/ch01/executed.ipynb` → COVERAGE OK (0 errors, 0 warnings).

## Final verdict: PASS
