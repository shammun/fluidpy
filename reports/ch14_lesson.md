# Chapter 14 — lesson review (re-review, round 2)                    2026-10-12

**This is the re-review after the round-1 fixes.** Round 1 (14:27, 541-cell notebook) ended FAIL with 15 Must fix and
23 Should fix. This round checks every one of those items against the re-executed notebook and repeats the full
first-time-reader pass.

coverage_check: OK (0 errors, 0 warnings) · sections covered 9/9 (14.1–14.9) · CORE blocks 16/16 with code + visual
(plain words, idea, tiny example, from-scratch cell with an assert, at least one figure/slider/animation with reading
notes in every block) · derivations 22/22 present (D01–D22), 4/4 ★★★ with a sympy check · explainers 9/9 embedded once
each (cells 101, 159, 224, 226, 262, 353, 447, 505, 545; all nine are the real pages, 250–295 kB each, no "not in this
copy yet" card) · 0 execution errors, 150/150 code cells executed · 240 calls to `fluidpy` physics modules in the
notebook, 0 positional.

Reviewed: `outputs/ch14/executed.ipynb` (15:58, 549 cells: 399 markdown, 150 code). It is newer than
`notebooks/build_ch14.py` (14:51), `scripts/ch14_drawings.py` (14:53) and `fluidpy/ch14_aerodynamics.py` (11:34), so it
was not re-run. Read top to bottom from a text dump; 38 PNG outputs, 3 videos (4 frames each) and 2 frame players
extracted to `outputs/ch14/lesson_review4/`; every figure named in round 1 was looked at, plus 20 others. Pages p719,
p739 and p741 rendered again: (14.1), (14.2), (14.14), (14.18)–(14.24) and slips #1, #2, #7 match the notebook.

## Round-1 item status

### Must fix — 15 of 15 FIXED, 0 partly, 0 not fixed

| # | Round-1 item (old cell) | Status | Now at | What I saw |
|---|---|---|---|---|
| 1 | bird upstroke sign and sketches (509) | FIXED | 515–517 | `oscillating_foil(U=10, V=-4, L=0.2)` prints −0.074 N (backward), power in −0.74 W, mean 0.33 N. Image `cell516_1.png`: relative wind downward, lift up and back, rose "backward part = L sin\|θ\|". Both sketches print V/U = 0.40, 21.8°, the cell's numbers. Text says the upstroke is "≈ 0", sign set by feathering; callout explains why the fish differs. |
| 2 | ASCII drag arrow to windward (515) | FIXED | 522 | `D ↙` back-and-left (leeward is left), `L ↖`; agrees with D22 step 4 and the figure. |
| 3 | "lift" printed on "lateral" (525) | FIXED | 533 | Image: "lift (normal to the apparent wind)" right of the teal arrow head, "lateral, to leeward (the keel resists it)" right of the amber line, "resultant F" left; no overlap. Windward and leeward sides are labelled. |
| 4 | sailing vocabulary never explained | FIXED | 524 | Primer P353 (bow, windward, leeward, abeam / beam reach, run, tack, heel) sits before D22 (525), tied to the sketch's orientation. |
| 5 | D15 step 1 three moves in one | FIXED | 375–378 | New primer P352 (cot A + cot B, 2 sin A sin B, d ln\|sin u\|) with a numeric demo; D15 is now 14 steps (candidate, differentiate, add cotangents, product of sines). The sympy cell now differentiates the antiderivative on both sides of γ₁ and takes the ε → 0 limit to get I₀ = 0. |
| 6 | triangle OQP labels overlap (123) | FIXED | 123 | Image: "OQ = b tan β" and "90° − θ" on leader lines, "QP = a = b/cos β" clear of QP, "OP = R" beside OP. |
| 7 | cusp figure titles collide (150) | FIXED | 150 | Image: three two-line titles, none touching. |
| 8 | three force triangles unreadable (343) | FIXED | 345 | Image: two-line titles (y, w, ε), arrows named in the right panel, tilt drawn ×10, rose arrow visible in all three. |
| 9 | C04 animation title clipped (152) | FIXED | 152 | Frames 0, 13, 26, 39: "β = 8.0° (0.140 rad), e = 0.038" fully inside the frame. |
| 10 | roll-up title clipped (288) | FIXED | 290 | All 5 frames: "t = 6.1 s · centroids 9.42 m apart · 1.46 m down" fits; the "qualitative" caveat moved to cell 292. |
| 11 | "elliptic loading", `Gn` used before defined | FIXED | 276–278 | Primer P351 with formula, span angle and a five-station demo, before the first use (282). |
| 12 | span efficiency used before defined | FIXED | 6, 405 | One-line meaning in cell 6; definition $e_{\text{span}}\equiv\Gamma_1^2/\sum n\Gamma_n^2\le1$ in cell 405 where 0.954 is first interpreted. |
| 13 | broken "Tools from earlier chapters" lines | FIXED | 21, 107, 165, 268, 304, 364, 395, 487 | All twelve reminder cells read as whole sentences about their use here; no stray brackets, raw LaTeX or foreign text. |
| 14 | text describes the wrong view (10) | FIXED | 10 | "a wing seen from behind (the stream goes into the page), with y along the span, z up". |
| 15 | comment's numbers not the output (401) | FIXED | 406 | Comment "Γ1 = 0.2505, Γ3 = 0.0276" = printed output. |

### Should fix — 21 of 23 FIXED, 2 partly

| # | Round-1 item | Status | What I saw |
|---|---|---|---|
| 1 | stagnation dot moves less than its size | FIXED | Cell 65 has a zoom panel; its title counts "dot 17.0 mm … 0.0 mm from the edge". |
| 2 | dashed negative ψ levels | FIXED | `linestyles="solid"` in 65, 150, 183; said in cell 66. |
| 3 | "tangent" not visible | FIXED | Cell 150 left panel is thickness against distance on log axes with slopes 3/2 and 1. |
| 4 | saturated colour map, grid over maps/bars | FIXED | Cell 321 uses a symmetric-log scale with readable ticks; no grid in 321, 387, 439. |
| 5 | legends on curves (92, 95) | FIXED | Cell 92 names the curves with text; cell 95's legend sits in the empty band. |
| 6 | label collisions (8, 41, 46) | FIXED | All three images clean ("camber" in 41 touches the outline but reads). |
| 7 | duplicate "What would change if…" | FIXED | 187/188 and 469/470 now differ (reading note, then the hand-over to the next block). |
| 8 | derivation headers D10, D11, D22 | FIXED | Each now shows the derived line. (D21: see new Should 3.) |
| 9 | ASCII corner letters (75) | FIXED | D bottom-middle, C bottom-right, as in the figure of cell 82. |
| 10 | "after a chord or two" (78) | FIXED | Says order of magnitude; 65 % after one chord, 90 % after about six. |
| 11 | (4.107)/(4.108) pairing in cell 0 | FIXED | $C_D$ (4.107), $C_L$ (4.108). |
| 12 | geometric series not in D05 Tools | FIXED | Named in Tools with the series written out. |
| 13 | D03 step 12 dense | FIXED | Split; D03 has 16 steps and all "step N" pointers (assumptions, notes cell 127, sympy comments) follow the new numbers. |
| 14 | `chebyshevt`, `chebyshevu`, `cancel` | FIXED | Cell 380 explains the five sympy functions. |
| 15 | circulatory lift, quasi-steady undefined | FIXED | Both defined in cell 96. |
| 16 | repr under widgets; empty stdout chunks | PARTLY | 221 and 416 now end in `_ = live(...)`. Cells 90 and 408 still emit an empty stdout chunk. |
| 17 | three orientations in C16; V label | FIXED | Cells 532 and 534 say which orientation each uses; sketches label V as the relative wind "(the wing moves down)". |
| 18 | N70 "rises as α²" | FIXED | Cell 475: section drag flat until stall; the α² part is the finite wing's induced drag. |
| 19 | explainer order in cell 1 | FIXED | Table order equals order of appearance. |
| 20 | comment 0.0393 | FIXED | 0.0392. |
| 21 | thin links for a climate reader | PARTLY | Cell 348 adds the ship-wave-drag analogy and mountain wave drag. Wind-turbine blades still only in the closing list. |
| 22 | sheet leaves the frame (288) | FIXED | ylim −6.5 m; the last frame's lowest point is at −5.2 m. |
| 23 | C12 had no static loading figure | FIXED | Cell 413: three taper ratios, loading against its ellipse and downwash, with reading notes. |

## Must fix (cell · what · the concrete change)

None.

## Should fix

1. **Cells 90 and 408** still emit an empty stdout chunk (a blank output box on the page). Left over from round-1
   Should 16.
2. **Cells 519–521** use "bow", "trimmed close" and "keel" before primer P353 (cell 524). "Keel" is glossed in 521;
   "bow" and "trimmed close" are not. Move P353 up to the start of §14.9 or gloss the two words in place.
3. **Derivation headers that show the starting equation, not the result:** D21 (cell 488) shows
   $C_{D_i}=C_L^2/(\pi\Lambda)$ (14.26) above a best-glide derivation; D02 (5.8), D05 (14.3), D06 (6.68), D08 (6.62) and
   D14 (14.13) likewise. Show the derived line, as D10, D11 and D22 now do.
4. **D22 step 8 (cell 525)** names "the tangent of a difference of two angles (gloss)" but never writes it. Add the
   line $\dfrac{\tan\beta_a-\tan\varepsilon_A}{1+\tan\beta_a\tan\varepsilon_A}=\tan(\beta_a-\varepsilon_A)$ between the
   quotient and the result.
5. **`lambda`** is first used in cell 112 (then 308, 327) but its reminder is in the tools list of cell 395. Move the
   reminder to the tools list of C04 (cell 107).
6. **"collocation"** appears in D15's "Fails when" (cell 377) before its reminder (cell 395); **"analytic"** (D06 step 1,
   R08) and **"winglets"** (D18, D19 "Fails when") have no gloss. One clause each.
7. **Cell 304** says "primer P344 below"; P344 is two cells above (302).
8. **Cell 423** comment: "comparable to the weight acting at the wing tip"; the printed moment is 126 391 N m against
   54 000 N m, 2.3 times larger. Say so.
9. **Long paragraphs:** "How to read it" in cell 96 (eight sentences) and "What does this show?" in cell 517 (seven).
   Split each in two; in 96 put the two definitions (circulatory lift, quasi-steady) in their own short paragraph.
10. **§14.8 sits outside its CORE block.** N77, N79 and the tail-stroke tiny example (cells 508–518) belong to C16,
    whose 🧩 heading opens in §14.9 (cell 520). Add one line to the §14.8 opening: "these notes are the first half of
    C16, which opens in §14.9".
11. **`scripts/ch14_drawings.py` still calls `ch14` functions positionally in 28 places** (for example
    `ch14.zhukhovsky_circle(b, e, beta)`, `ch14.kutta_circulation(U, a, alpha, beta)`,
    `ch14.apparent_wind(I["U_wind"], I["V_boat"], I["theta"])`). The notebook itself has none, and
    `oscillating_foil` there is now called by keyword, but a reader who opens the helper the notebook imports sees
    the other style. Convert them.
12. **Wind-turbine and propeller blades** are still only in the closing list. One sentence in C16 would do: a
    rotating blade is the same velocity triangle with the blade's own speed in place of the tail's.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / tool | First used | Explained |
|---|---|---|
| bow, "trimmed close" | 519, 521 | P353 at 524 (bow); "trimmed close" — **missing** (Should 2) |
| tan(A − B) identity | 525 (D22 step 8) | named as a gloss, formula **not written** (Should 4) |
| `lambda` | 112 | reminder at 395 — **late** (Should 5) |
| collocation | 377 | 395 — **late** (Should 6) |
| analytic function | 166 (D06 step 1), 171 | Ch. 6 prerequisite, **no reminder** (Should 6) |
| winglets | 434, 452 | **missing** (Should 6) |

Round-1 rows now resolved: elliptic loading and `Gn` (P351, cell 276), span efficiency (cells 6, 405), sailing words
(P353), cotangent and sine-product identities (P352), geometric series (D05 Tools), Chebyshev functions (cell 380),
circulatory lift and quasi-steady (cell 96), dashed contours (drawn solid, said in cell 66).

Everything else checked at first use was explained by a primer (P336–P353), a recap (R01–R12), the conventions cell 7 or
an inline gloss: stagnation point, complex potential, conformal map, Mach number, dynamic pressure, planform words,
section words, vortex sheet, principal value, Kronecker delta, condition number, washout, displacement thickness,
adverse pressure gradient, feathering, drag angle, Rankine core, `np.argmin`, `np.interp`, `quad(weight="cauchy")`,
`np.linalg.solve`, `np.outer`, `np.hypot`, `@`, `semilogy`, `loglog`, `SymLogNorm`, `**SEC`.

## Derivation audit (D id · steps · every step follows? · why/in-words ok · check ok · verdict)

| D | Steps | Follows | Why / in words | Check | Verdict |
|---|---|---|---|---|---|
| D01 ★ | 5 | yes | ok | units, limits, number | ok |
| D02 ★ | 5 | yes | ok | units, number, code | ok |
| D03 ★★★ | 16 | yes (steps 12–15 re-derived; 3/4 example gives −0.583, 2.083) | ok | sympy rebuilds (14.4) → (14.5) → (14.6) → arc; limits; number | ok |
| D04 ★★★ | 11 | yes (re-derived to 3√3/4) | ok | sympy series, maximum, quarter chord; rule against measured | ok |
| D05 ★★ | 7 | yes | ok; series now in Tools | units, number, code | ok |
| D06 ★★ | 7 | yes | ok | far field, streamline, circulation, code | ok |
| D07 ★ | 4 | yes; matches (14.9) | ok | units, limits, number | ok |
| D08 ★★ | 5 | yes | ok | loop sums in both planes (−10.0000 six times) | ok |
| D09 ★★ | 10 | yes; matches (14.10)–(14.12) | ok | 15.96 m²/s, 1.058, 0.987 reproduced | ok |
| D10 ★★ | 7 | yes | ok | 1.054/cos 5° = 1.058 | ok |
| D11 ★ | 4 | yes | ok | staircase 2 + 3 + 5 = 10 | ok |
| D12 ★★ | 7 | yes; matches (14.13), slip #3 | ok | horseshoe limit, number | ok |
| D13 ★★ | 8 | yes; matches (14.14), (14.15), slip #2 (page p739) | ok | number | ok |
| D14 ★★ | 7 | yes; slip #1 confirmed (page p741) | ok | sign check +0.5 / −0.5 / +0.5 | ok |
| D15 ★★★ | 14 | yes (steps 2–4 and 6 re-derived: F(0) = F(π) = 0, the two ln sin(ε/2) cancel) | ok | sympy: identity, recursion, both starting values, dF/dγ on both sides, I₀ limit, n = 1…5 | ok |
| D16 ★★ | 8 | yes; matches (14.16), (14.17), (14.20) | ok | limits, 0.217 Uc | ok |
| D17 ★★ | 6 | yes; matches (14.21) | ok | number 6786 N | ok |
| D18 ★★★ | 9 | yes; matches (14.22) | ok | sympy integrates the 16 products | ok |
| D19 ★★ | 9 | yes; matches (14.23)–(14.26) | ok | 318 N two ways | ok |
| D20 ★★ | 8 | yes | ok | limits, solver returns one coefficient | ok |
| D21 ★★ | 8 | yes | ok | 16.5 at 34.5 m/s, brute force | ok — header (Should 3) |
| D22 ★★ | 8 | yes (vectors and the quotient re-derived) | ok | limits; 65 N / 80 N; 375 N / 981 N | ok — write the identity (Should 4) |

## The two upstream corrections

- **Wagner's function is a lift ratio:** holds. Cell 95 plots `ch14.wagner_function(s=2*s_ch)` only on the right-hand
  lift panel, against `run["lift_circulatory"]/(1.0*sec["U"]*run["Gamma_steady"][0])`, axis "lift / steady lift [–]",
  label "circulatory lift / steady lift (Jones' fit)". The left panel carries circulations only. Cell 96 says it is a
  lift ratio that starts at one half while the exact circulation starts from zero, and "The model reaches 90 % of the
  steady value after about 6.1 chords (with this step of 0.05 chords)"; the printed value is 6.13. Nowhere is it called
  the circulation growth.
- **Keyword calls:** hold in the notebook. An AST scan of all 150 code cells found 240 calls to `AF`, `LL`, `ch14`,
  `ch06`, `CONF`, `POT`, `BS`, `SIM` and to the `draw_…` / `fig_…` helpers with an argument, none positional. The only
  positional calls on physics objects are four evaluations of the callables that `CONF.circle_flow_zeta` returns
  (`w_c(point)`, cells 63, 169, 183). The helper script is a different matter (Should 11).

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- A "first look" primer (P351) for a concept a later block owns: formula, five numbers, and a pointer to where it is
  taught properly. It removes a forward reference without moving the block.
- A zoom panel whose title counts the quantity the reader should watch ("dot 17.0 mm from the edge").
- Naming arrows once, in the panel where they are largest, and stating the exaggeration ("tilt drawn ×10").
- The bird callout: the fish gives thrust on both half-strokes because its lift changes sides; a wing that must keep
  lifting upward cannot. One contrast replaces a wrong symmetry.
- Sympy checks that rebuild the derivation's own steps (D15: the antiderivative on both sides of the singular point,
  then the limit).
- From round 1, still true: the conventions cell before any physics, "slip #k" boxes with a numerical demonstration,
  from-scratch cells ending in an assert in all 16 blocks, "approximation made measurable" tables, the staircase of
  horseshoes converging to the principal-value integral.

## Verdict: PASS

0 Must fix, 12 Should fix. All 15 round-1 Must-fix items are fixed, and the figure items were confirmed on the
regenerated images and frames. No wrong equation, no broken derivation step, no number in the prose that disagrees
with its cell. The Should items are cosmetic or small glosses; none blocks or misleads a first-time reader.
