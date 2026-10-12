# Chapter 14 — lesson review                                   2026-10-12

coverage_check: OK (0 errors, 3 warnings — D19, D21, D22: explainer not built yet, Derivation tab not checked) ·
sections covered 9/9 (14.1–14.9) · CORE blocks 16/16 with code + visual (plain words, idea, tiny example, from-scratch
cell and at least one figure/slider/animation in every block) · derivations 22/22 present (D01–D22), 4/4 ★★★ with a
sympy check · explainers 9/9 embedded once each (cells 101, 159, 224, 226, 262, 351, 440, 498, 537; three show the
"not in this copy yet" card, which is not counted against the notebook) · 0 execution errors · all `fluidpy` calls by keyword.

Reviewed: `outputs/ch14/executed.ipynb` (541 cells), read top to bottom; all 37 PNG outputs looked at; sample frames
pulled from the 3 videos and 2 frame players (`outputs/ch14/lesson_review/`); 22 numbered equations compared with the
rendered pages p731, p732, p734, p738, p740, p741, p742.

## Must fix (cell · what · the concrete change)

1. **Cell 509 (N79, §14.8 bird flight) · the upstroke has the wrong sign and the two sketches contradict the code.**
   `ch14.oscillating_foil(U=10, V=4, L=0.2)` is reused for the upstroke and prints "+0.074 N thrust on the upstroke; mean
   0.41 N". With the wing moving *up*, the relative wind comes from above-ahead, so a small positive lift leans
   *backward* (thrust −0.074 N, mean 0.33 N); it is forward only if the lift is negative. `draw_velocity_triangle("bird_up")`
   draws the same triangle as the downstroke (V arrow up, lift leaning forward) under the label "wing moving up". The
   sketches also print V/U = 0.50 (26.6°) and 0.80 (38.7°) while the cell and cell 510 use V/U = 0.4 (21.8°), and carry
   the fish caption "thrust power = sideways power". Change: flip the V arrow and the lean of the lift in `bird_up`, give
   the upstroke the sign it has (or state "≈ 0, sign depends on feathering" and drop the mean), and pass the cell's own
   U, V to both sketches.
2. **Cell 515 (C16, "The idea" ASCII) · the drag arrow points to windward.** The sketch has the wind from ahead-right
   (leeward = left) and draws `↘ D drag ∥ U_r`. D22 step 4 and the figure of cell 525 have drag "back and to leeward",
   i.e. `↙`, parallel to the apparent-wind arrow `↓↙`. Change the arrow (and tilt the lift arrow forward, e.g. `↖`).
3. **Cell 525 (C16, `scripts/ch14_drawings.py::draw_sail_polygon`) · confirmed: the label "lift" is printed on top of
   "lateral (the keel resists it)"** (reads "liftateral"). The lift is the vector the block is about. Move "lift" to the
   outer side of the teal arrow (mid-arrow, offset to the right) and put "lateral" left of the amber arrow.
4. **Cells 517, 520, 526, 529 (D22 and its figures) · sailing vocabulary is used and never explained:** "to leeward"
   (D22 step 1), "windward side" (step 2), "from exactly abeam" / "beam reach" (cell 520), "forward of the beam" (526),
   "a run" / "running" (529, D22 check), "tacks" before N82, "heels". Only the code comment of cell 523 glosses leeward.
   Add a 📎 primer before cell 517: bow, leeward/windward, abeam (beam reach), run, tack, heel — one line each, with the
   sketch's orientation.
5. **Cell 373 (D15 ★★★, step 1) · three moves in one "why".** "Differentiate the right side: ½cot + ½cot = sin γ₁/(cos γ −
   cos γ₁)" skips (a) d/dγ ln|sin u| = u′ cot u, (b) cot A + cot B = sin(A + B)/(sin A sin B) with A + B = γ₁, (c)
   2 sin A sin B = cos(A − B) − cos(A + B) with A − B = γ — (b) and (c) are in no primer (P347 has the cosine identity
   only). Split into steps or add the three lines; and add to the check cell 374 a `sp.diff` test of this antiderivative
   (the cell assumes I₀ = 0 in a comment and never tests step 1–3).
6. **Cell 123 (`draw_triangle_OQP`, used by D03 steps 6–7) · the labels "90° − θ" and "b tan β" overlap** and the
   length OQ = b tan β cannot be read; the radius label "a = b/cos β" is struck through by the line QP. Move the three
   labels off the lines.
7. **Cell 150 (C04, cusp figure) · the three panel titles run into one another** ("…become tangentedge of 17°: speed ∝
   r^0.05 → 0 at the edgecusp: speed…"). Shorten or wrap the titles / widen the figure.
8. **Cell 343 (C10, three force triangles) · the three titles overlap and the read-out (y, w, ε) is unreadable;** the
   arrows have no labels and the rose induced-drag arrow is a dot in the first two panels. Wrap the titles on two lines,
   label the arrows once, draw the tilt ×10.
9. **Cell 152 (C04 animation) · the left title — the β, e read-out that the animation is about — is cut off at the top
   edge in every frame** (frozen layout). Reserve the top margin before `fig.set_layout_engine("none")` or move the
   read-out into an `ax.text`.
10. **Cell 288 (C08 roll-up frame player) · the title is clipped at both ends** ("alitative past the first turn … 1.46").
    Shorten it ("t = 6.1 s · centroids 9.42 m apart · 1.46 m down") and put "qualitative past the first turn" in the notes.
11. **Cells 280–285 (C08) · "elliptic loading" is used (code comments, cell 281, cell 285 "the smooth elliptic loading")
    long before it is defined** (formula first in cell 365, block C14 at cell 442); the sine coefficients `Gn` are used
    in cells 276 and 280 with only "(C11)". Add at cell 280: "elliptic loading: Γ(y) = Γ₀√(1 − (2y/s)²), a half-ellipse
    over the span — C14 shows why it is special; `Gn=[Γ₁]` is its amplitude (C11)".
12. **Cells 399–407 (C12) · "span efficiency" (`e_span` = 0.954, "peaks at 0.99") is printed and interpreted before it is
    defined** (N54, cell 431); cell 400 says only "explained in C13". Give the one-line definition there:
    $e_{\text{span}}=\Gamma_1^2/\sum n\Gamma_n^2\le1$, 1 for the cheapest loading.
13. **Cells 362, 480, 302, 390 (and 21, 107, 268) · the "Tools from earlier chapters" lines are broken or foreign text.**
    362: unclosed brackets and raw LaTeX ("(L = ρUΓ. unknown = …", "overline{R_n^2}=nL^2", "(even W one row of cells,
    odd W two."); 480: "used for and to expose slip #3. the minimum…"; 302: "a segment grows into an infinite line in.";
    390: "(`np.vander` as the classic bad matrix."; 21: `np.interp` "reads the envelope's height at a crest"; 107:
    "Howard's last step with x = c_r", "an imaginary frequency means growth". Rewrite each reminder as one whole
    sentence about its use in this chapter.
14. **Cell 10 · the text describes the wrong view.** "Right: a wing from above with y along the span" — the figure is
    titled "wing seen from behind (stream into the page)" with z up. Correct the sentence.
15. **Cell 401 (C12 from-scratch) · the comment's numbers are not the output:** comment "Γ1 = 0.2474, Γ3 = 0.0297",
    printed 0.2505 and 0.0276. Correct the comment.

## Should fix

1. Cell 65 (C02 frame player): the amber stagnation dot travels 1.7 % of the chord (x = 0.483 → 0.500 m) — less than
   its own size — so "the dot slides back to the edge" is not visible; add a zoom on the last 10 % of the plate or a
   ζ-plane panel where the point moves 15°.
2. Cells 65, 150, 183: negative stream-function levels are drawn dashed (matplotlib default) and nothing says so;
   use `linestyles="solid"` (cell 79 already does).
3. Cell 150, left panel: two straight lines meeting at about 2° do not show "tangent"; plot thickness against
   distance from the tail on log axes (slope 3/2) or stretch the vertical scale.
4. Cell 319: the colour scale is cut at ±0.6 m/s and almost the whole map is saturated; widen the cut or use a
   symmetric-log norm; remove the white grid over the map (also over the bars of cells 382, 432).
5. Cells 92 and 95: legends lie on the curves (the "their sum" line strikes through the legend text in the animation;
   the legend covers the shed curve in the left panel of 95).
6. Label collisions: cell 8 ("lift normal to the stream" on the y-axis; "+s/2 starboard" on "y (span)"), cell 41
   ("thickness t_max" on "camber line"), cell 46 ("yaw axis" on the fuselage).
7. Duplicate cells: 187/188 and 462/463 repeat the same "What would change if…".
8. Derivation headers show an equation that is neither the start nor the result: D10 (14.12), D11 (14.13),
   D22 (14.15 — lifting-line integrals above a sailing derivation). Show the derived line instead.
9. Cell 75 ASCII puts D at the bottom-left corner; the figure of cell 82 and D02 need D at the bottom-middle (the
   cut BD). The top-right corner has no letter in either.
10. Cell 78 (N14): "Γ = Γ_Kutta for t > c/U — after a chord or two" against the model's 65 % after one chord and
    90 % after 6.1 (cells 88, 96). Say that the book's sketch gives the order of magnitude.
11. Cell 0: "$C_L$ … and $C_D$ … (4.107), (4.108)" reverses the pairing (C_D is 4.107, C_L is 4.108, as cell 30 has it).
12. D05 step 4 (cell 147): the second-order expansion is justified by "the geometric series" while Tools lists
    "Taylor expansion (P26)", which is first order; name the geometric series in Tools with one line.
13. D03 step 12 (cell 125): two substitutions and a collection in one step; split in two.
14. Cell 374: `sp.chebyshevt`, `sp.chebyshevu`, `sp.cancel` are explained only by half a comment; add one sentence.
15. Cell 96: "circulatory lift" and "quasi-steady" are used without a definition (quasi-steady is explained at cell 531).
16. Cells 221 and 409 print the function repr under the widget; cells 90 and 210 emit an empty stdout chunk.
17. C16 uses three orientations (ASCII: boat heading up; cell 523: wind toward −x; cell 525: course along +x) and
    cell 526 does not say which; in the fish and bird-downstroke sketches the V arrow is the relative-flow component
    but is labelled as the body's motion ("wing moving down" on an upward arrow).
18. Cell 468 (N70) says the drag "rises roughly as α²"; the model curve of cell 469 is flat until stall.
19. Cell 1 lists the explainers in an order (builder, Kutta, starting vortex) different from where they appear
    (starting vortex, builder, Kutta).
20. Cell 454: comment 0.0393, output 0.0392.
21. Links for a climate-dynamics reader are thin: the book's own analogy between induced drag and a ship's wave drag
    (Ch. 7) is not used; wind-turbine blades appear only in the closing list.
22. Cell 288, last frame: the middle of the sheet has left the frame (ylim −4.5 m).
23. C12 has only a plotly slider and a live cell; a static loading figure (as cell 458 gives in C14) would survive
    a page without JavaScript.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / tool | First used | Explained |
|---|---|---|
| elliptic loading | 280 (comment), 281, 285 | 365 (formula), 445 (C14) — **missing at first use** |
| sine coefficients `Gn`, span angle γ | 276, 280 | 353–363 (C11); 280 glosses γ, `Gn` only "(C11)" — **thin** |
| span efficiency `e_span` | 5 (dict), 399–407 | 431 (N54) — **missing at first use** |
| leeward, windward, abeam / beam reach, run, tack, heel | 517, 520, 526, 529 | comment in 523 (leeward only) — **missing** |
| cot A + cot B, 2 sin A sin B identities; d ln sin | 373 (D15 step 1) | **missing** (P347 has the cosine identity only) |
| geometric series to second order | 147 (D05 step 4) | named in the step, not in Tools — thin |
| `sp.chebyshevt`, `sp.chebyshevu`, `sp.cancel` | 374 | half a comment — thin |
| circulatory lift, quasi-steady | 95, 96 | 531 (quasi-steady); circulatory lift **missing** |
| `LL.solve_lifting_line` | 343 | 385–404 (C12); forward pointer given — acceptable |
| dashed contour lines (negative ψ) | 65, 150, 183 | **missing** |

Everything else checked at first use was explained (primers P336–P350, recaps R01–R12, the conventions cell 7, or an
inline gloss): stagnation point, complex potential, conformal map, Mach number, dynamic pressure, vortex sheet, principal
value, Kronecker delta, condition number, washout, displacement thickness, adverse pressure gradient, Rankine core,
`np.argmin`, `np.interp`, `quad(weight="cauchy")`, `np.linalg.solve`, `np.outer`, `np.hypot`, `@`, `semilogy`.

## Derivation audit (D id · steps · every step follows? · why/in-words ok · check ok · verdict)

| D | Steps | Follows | Why / in words | Check | Verdict |
|---|---|---|---|---|---|
| D01 ★ | 5 | yes | ok | units, number | ok |
| D02 ★ | 5 | yes | ok | units, number, code | ok |
| D03 ★★★ | 15 | yes (re-derived; step 12 dense) | ok | sympy re-does (14.4)→(14.5)→(14.6)→arc; limits; number | ok — split step 12 |
| D04 ★★★ | 11 | yes (re-derived to 3√3/4) | ok | sympy series + maximum; table rule vs exact | ok |
| D05 ★★ | 7 | yes | ok; geometric series not in Tools | units, number, code | ok — Should 12 |
| D06 ★★ | 7 | yes | ok | far field, streamline, circulation, code | ok |
| D07 ★ | 4 | yes; matches (14.9) | ok | units, limits, number | ok |
| D08 ★★ | 5 | yes | ok | loop sums in both planes | ok |
| D09 ★★ | 10 | yes; matches (14.10)–(14.12) | ok | numbers reproduced (15.96 m²/s, 1.058, 0.987) | ok |
| D10 ★★ | 7 | yes | ok | number 1.054/cos 5° = 1.058 | ok — header equation |
| D11 ★ | 4 | yes | ok | staircase number | ok — header equation |
| D12 ★★ | 7 | yes; matches (14.13) and slip #3 | ok | horseshoe limit, number | ok |
| D13 ★★ | 8 | yes; matches (14.14), (14.15) | ok | number | ok |
| D14 ★★ | 7 | yes; slip #1 confirmed on p741 | ok | sign check, code | ok |
| D15 ★★★ | 11 | **step 1 skips three moves** | rest ok | sympy tests identity, recursion, n = 1…5; not step 1, not I₀ | **Must 5** |
| D16 ★★ | 8 | yes; matches (14.16), (14.17), (14.20) | ok | limits, number 0.217 Uc | ok |
| D17 ★★ | 6 | yes; matches (14.21) | ok | number | ok |
| D18 ★★★ | 9 | yes; matches (14.22) | ok | sympy integrates the 16 products | ok |
| D19 ★★ | 9 | yes; matches (14.23)–(14.26) | ok | number 318 N | ok |
| D20 ★★ | 8 | yes; starts from (14.27) | ok | limits, solver | ok |
| D21 ★★ | 8 | yes | ok | number 16.5, brute force | ok |
| D22 ★★ | 8 | yes (vectors and tan(β_a − ε_A) re-derived) | ok; vocabulary unexplained | limits, numbers reproduced (375 N, 981 N) | ok — Must 4, header equation |

## Points the orchestrator asked about

- **C03 / Wagner:** correct. Wagner's function is on the right-hand *lift* panel only, labelled "circulatory lift /
  steady lift (Jones' fit)", argument in semi-chords; the text says the exact circulation starts from zero and the
  circulatory lift at one half; "90 % after about 6.1 chords (with this step of 0.05 chords)" matches the printed 6.13.
- **Keyword calls:** no positional `fluidpy` call found (script scan of all 147 code cells plus reading).
- **C16 leeward vector:** `n_hat = (−sin θ, cos θ)` is right for the cell's axes (wind toward −x): U·n̂ = +U sin θ, as
  D22 step 2 has it; drive 374.6 N and side force 980.9 N equal `ch14.sail_forces`; the figure of cell 525 agrees in
  direction. Not consistent: the ASCII sketch (Must 2) and the label overlap (Must 3, confirmed).
- **Animations at 40 frames** (cells 92, 152, 530; 3.2–3.6 s each): all three still teach — the starting vortex and the
  flat amber sum, the two knobs one after the other, thrust and drive forward on both half-strokes/tacks. The defect is
  the clipped title of cell 152 (Must 9), not the frame count.
- **Equations against pages:** (14.4)–(14.8), the arc equation, (14.9)–(14.13), (14.15)–(14.27) all match the page
  images, with slips #1, #3, #4, #6 as described. (14.1)–(14.3), (14.14) and the Ch. 4–6 numbers were not rendered.

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- The conventions cell (clockwise Γ, two sets of axes, overloaded letters) placed before any physics, and "slip #k —
  the book prints … ; the correct form is …" boxes with a numerical demonstration (slip #1: +0.5 / −0.5 / +0.5).
- From-scratch cells that end in an `assert` against the library in all 16 blocks.
- "Make the approximation measurable" tables: thickness rule against measured, N-term convergence, taper against
  span efficiency.
- C09's staircase of horseshoes converging to the principal-value integral (0.637 → 0.9998): the integral explained by
  its discrete limit.
- D15's route to Glauert's integral (two starting values, a recursion, a candidate, induction) and its Chebyshev check.
- C13's table "net lift / induced drag per sine" and the bar chart behind it.
- Tiny examples with 3–4–5 triangles and sines of 30° that the code cells then reproduce digit for digit.

## Verdict: FAIL

15 Must fix, 23 Should fix. No wrong equation and no broken derivation result; the failures are one sign error with
its sketches (bird upstroke), one wrong arrow, six unreadable figure labels/titles, three terms and one vocabulary set
used before they are explained, one compressed ★★★ step, broken reminder lines, and two text/number mismatches.

Addendum (outside the notebook cells): `scripts/ch14_drawings.py::draw_velocity_triangle` calls `ch14.oscillating_foil(U, V, L, D)` positionally and hard-codes V/U = 0.5 and 0.8 for the two bird sketches; `oscillating_foil` assumes "the same triangle on both half-strokes" (right for a symmetric fish tail whose lift reverses, not for a bird upstroke with positive lift) - the source of Must 1.
