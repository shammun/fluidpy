# Chapter 15 — lesson review (round 2)                                   2026-10-14

coverage_check: OK (0 errors, 0 warnings) · sections covered 10/10 (15.1–15.10, plus S01/S02) · CORE blocks 16/16 with
code + visual · derivations 31/31 present (D01–D31), 4/4 ★★★ with a sympy check that rebuilds the derivation's steps ·
9/9 explainers embedded (none shows the "not in this copy yet" card) · 0 execution errors.

Reviewed: `outputs/ch15/executed.ipynb` (00:57, 606 cells), which is newer than `notebooks/build_ch15.py` (00:54). I dumped
it to text and extracted its 25 PNG outputs. I diffed the dump against my round-1 dump: 501 changed lines once blank lines
are ignored. I read every changed passage in context, and re-read the derivations whose "why" fields were rewritten (D14,
D15, D16, D20, D21, D27–D30). I looked at the six figures named in round 1: c585 (C16 plate and force inset), c403 (N76 jet
sketch), c523, c272, c120 and c067 (label collisions). I also looked at c287 (frame sketch).

## Verdict: PASS

All 3 round-1 Must-fix items and all 21 round-1 Should-fix items are fixed. Nothing new rises to Must fix. There are
4 new Should-fix items, all wording or cosmetic, plus 1 note for the publisher about stale embeds.

## Round-1 items

| # | Item | Status | Evidence |
|---|---|---|---|
| M1 | C16 plate figure had no force arrows; p₃ label crossed by shock | **fixed** | c585 now has an inset (lower left): N ⟂ plate (purple), L vertical (teal), D horizontal downstream (rose), with "tilt drawn ×5: 10° for α = 2°". The p₃ label is clear of the lower shock. Cell 586 describes this figure. |
| M2 | N76 jet sketch titles collided | **fixed** | c403: each panel has a two-line title ("over-expanded jet (a sketch)" / "$p_e<p_B$: oblique shocks, jet contracts"), and they do not touch. |
| M3 | $4f_FL/D$ never defined | **fixed** | Cell 7(a) defines $f_F=\tau_w/(\tfrac12\rho u^2)$, L and D (hydraulic). N86 (cell 445) adds "Reading the right-hand column" with the slice bridge $df=\dfrac{4f_F\,dx}{D}\cdot\tfrac12\gamma M^2\dfrac{p}{p_1}$, which I checked from $dF=\tau_w(4A/D)dx$, $f=F/(p_1A)$ and $\tfrac12\rho u^2=\tfrac12\gamma pM^2$. It also explains $T_0^*$ and warns that D is not drag and L is not lift. The T5/T11 callout (cell 409) points to it. |
| S1 | stray "> ⚠️" in mid-paragraph | fixed | every callout now starts its own line (the one remaining, cell 447, is the start of a list item and renders as a quote inside the list) |
| S2 | cell 33 printed (15.7) twice | fixed | |
| S3 | two M values with no reason | fixed | comment at cell 13 (R = 287.058 against 287.05, trap T4) |
| S4 | "parabola" for an absolute-value plot | fixed | "a W: two V-shaped dips" |
| S5 | `lambda` had no reminder | fixed | comment at first use: "Ch. 1 primer P29" |
| S6 | slider grids missing the case read | fixed | γ lists contain 1.4 and open on it (4 figures); s list contains 0.1 and opens on it; Mach list contains 2.0 and opens on it |
| S7 | "within 6 %" | fixed | "within about 6 % … at 0.939" |
| S8 | "two-row" against a side-by-side figure | fixed | "left … right", which matches c287 |
| S9 | "agree to third order" | fixed | "share value, slope and curvature … agree through second order and part only at third order" |
| S10 | cell 342 printout before m, V, I defined | fixed | N66 (cell 342) defines m, V, I with units; the code (345) follows, then a "What does the code above do?" saying that 0 is the residual of the reduction. I checked V = u₁ + p₁/m = 802.1 m/s and U₁ = 5.6/6.6. |
| S11 | Kutta and d'Alembert used before the recap | fixed | recap is now cell 575, right after the C16 heading, before the first use (576) |
| S12 | four label collisions | fixed | c523 "M₂ = 1.08" is now under the title on a white box. In c272 the label sits beside the band, clear of the curve. In c120 "F" is under its arrow, clear of the tail. In c067 the labels sit to the right of the dashed lines, on white. |
| S13 | copy-paste `return out` comment | fixed | "all artists of this frame (circles, sources, Mach lines)" |
| S14 | D20 step 4: $p_d<p_b$ unjustified | fixed | The argument is now $p_d/p_0=(p_{02}/p_{01})(p/p_{02})$: both factors are smaller (D18, and a smaller $A_e/A_2^*$ means a faster subsonic exit). It is correct. |
| S15 | padded "why" fields | fixed | each now names the rule and the purpose; a few are still a little long (new S4) |
| S16 | D29 start line = step 1 | fixed | starts from $v_1=v_2$ (D26 step 3) |
| S17 | duplicate note IDs N26, N27, N102 | fixed | the second use of each now reads "…, continued" |
| S18 | "twice the 2° plate" | fixed | "about twice … 0.00281: the diamond's slopes are 2.86°" |
| S19 | bare equation numbers | fixed | C12 heading names the differentials; cell 415 writes out (15.45) and (15.46); D31 step and trap write out (15.47) |
| S20 | "(verification report)" in prose | fixed | replaced by the reason (0.05 % at 2°, about 1 % at 10°) |
| S21 | identical boilerplate in 9 🎮 cells | fixed | each is tailored (p′/ρ₀c²; A/A*; station 1/2; p₀ and p_B; f and q; μ; σ/δ against β/θ; ν with the sign of the turn; α and b) |

**Round 1: 24 fixed / 0 not fixed.**

## Must fix

None.

## Should fix

1. **Cell 586 (C16, How to read it) puts α at the wrong corner of the triangle.** The text says the three arrows "form a right
   triangle with the angle α at the top". In c585, L and N both start at the plate, so α sits between them at the **bottom**
   (where the arrows start). The right angle is at the top of the teal arrow. Change "at the top" to "at the bottom, between
   the teal and purple arrows".
2. **Cell 448 (N87): the trap T18 callout is a fragment.** "> ⚠️ Common confusion (trap T18): a statement about one band of
   Mach number, and the entropy still rises." does not say what the confusion is. Suggested text: "heating *cools* the gas
   only for $1/\gamma<M^2<1$. That is a statement about one band of Mach number, not about heating in general, and the
   entropy still rises."
3. **Slip callouts in cell 38 (slip #1) and cell 584 (slip #13) show literal asterisks.** The helper writes
   `…gases**; the correct form is**` and `…plate**; the correct form is**`. A closing `**` that comes right after a letter
   and right before ";" cannot open bold in CommonMark, so the page prints `**`. The other four slip callouts work because
   their text before the `**` ends in ")" or "$". Fix it in the helper: emit `**the correct form is**` after the
   semicolon, i.e. `…gases; **the correct form is** …`.
4. **A few "why" fields still say the same thing twice** (round-1 S15 follow-up; minor):
   - D29 step 5: "Logarithmic differential of V = Mc. Taking the logarithm of V = Mc and differentiating is allowed…"
   - D14 step 6: the numerator algebra comes before "Subtracting 1 from both sides…". Put the move first and the algebra
     second.
   - D21 step 4 is just "Rearrangement." Name what is solved for, e.g. "divide by $p_1M_1$…".

**Note for publish (not a lesson defect).** `viz/ch15/oblique_shock_wedge.html` (00:59) and `prandtl_meyer_fan.html`
(01:00) were saved **after** this execution (00:57). The copies embedded in cells 536 and 571 may therefore be older
than the files. Re-execute (or re-embed) once the viz-reviewer has passed them.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / tool | First used | Explained |
|---|---|---|
| $f_F$, L, D in "$4f_FL/D$" | 7(a) | **7(a)**, then the bridge in N86 (445) — fixed |
| `lambda` | 15 | comment at 15 (Ch. 1 P29) — fixed |
| m, V, I (shock-structure constants) | 342 | 342 (N66), before the code at 345 — fixed |
| Kutta condition, d'Alembert's paradox | 576 | 575 (recap) — fixed |
| mean free path | 258 (N48) | inline gloss in N48 (about 7 × 10⁻⁸ m, Ch. 1 C06) — fixed |
| `**kw` / `*bars` unpacking | 456, 458 | not primed (acceptable for this reader) |
| p_rms | 52 | not glossed (acceptable for a statistics reader) |

No term is missing. Everything else was explained at first use in round 1, and the edits introduced no new symbol
without a gloss. $\tau_w$, Fanning, hydraulic diameter and wetted perimeter are each defined in the sentence that first
uses them.

## Derivation audit (D id · steps · every step follows? · why/in-words ok · check ok · verdict)

Unchanged derivations keep their round-1 verdicts (all ok). The rewritten fields were re-checked:

| D | Change | Follows | Why / in words | Check | Verdict |
|---|---|---|---|---|---|
| D14 ★★ | steps 3, 5, 6 why | yes. Step 6 numerator $(\gamma+1)(1+\gamma M_1^2)-[2\gamma M_1^2-(\gamma-1)]=2\gamma+\gamma(\gamma-1)M_1^2$ re-derived | ok (order in step 6, S4) | 1/3 at M₁ = 2 | ok |
| D15 ★★ | steps 1, 3 why | yes | ok | unchanged | ok |
| D16 ★ | steps 3, 4 why | yes ($C_v+R=\gamma C_v$) | ok | unchanged | ok |
| D20 ★★ | step 4 argument for $p_d<p_b$ | **yes, now complete**: smaller $A_e/A_2^*$ → larger subsonic $M_e$ → smaller $p/p_{02}$; $p_{02}<p_{01}$ | ok (forward reference to step 6 is marked) | 0.0939 < 0.5134 < 0.9372 | ok |
| D21 ★★★ | step 5 why | yes (divide step 1 by $p_2(1+\gamma M_1^2-f)$) | ok; step 4 terse (S4) | sympy unchanged | ok |
| D27 ★★ | step 1 why | yes | ok | unchanged | ok |
| D28 ★★ | steps 2, 3, 6, 7 why | yes (Pythagoras for cot μ = √(M²−1)) | ok | unchanged | ok |
| D29 ★★ | start line | yes ($v_1=v_2$, then step 1 writes $V_1\cos\sigma=V_2\cos(\sigma-\delta)$) | ok; step 5 padded (S4) | unchanged | ok |
| D30 ★★★ | steps 3, 6, 8, 10 why | yes (prefactor $\sqrt{(\gamma+1)/(\gamma-1)}\cdot\ldots=1$) | ok | sympy unchanged | ok |
| D31 ★ | step 2 why names slip #13 with (15.47) written out | yes | ok | 0.0806 | ok |

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **The force inset with a stated exaggeration** (c585): "tilt drawn ×5: 10° for α = 2°", plus the sentence saying the true
  drag arrow would be 3.5 % of lift and invisible. This is the honest way to draw a small-angle decomposition; reuse it
  wherever a component is O(α).
- **Bridge between rival parameters, done slice by slice** (N86): a one-line differential that converts the book's lumped f
  into the usual $4f_FL/D$. The cell then prints both numbers for the same duct (0.100 against 0.532). This is a
  `teaching-style` pattern for any "book's parameter against the field's parameter" difference.
- **Sliders open on the case the text reads** (`active=values.index(1.4)` and explicit value lists). It should be a
  python-viz rule: the reading note and the default slider position must agree.
- **An inequality proved by factoring a ratio** (D20 step 4): write $p_d/p_0$ as a product and show each factor is smaller.
- All round-1 "works well" items still stand: the slip settled by finite differences of the exact formula (cell 440),
  three routes to one number (330), the measurable-approximation tables, the non-additive f (446), "divide out the known
  root", and the A3 back-pressure animation.

## Verdict: PASS

There are no Must-fix items. All 24 round-1 items are fixed and checked in the executed notebook, including the six
figures I looked at. coverage_check is OK with 0 errors and 0 warnings, and all 9 explainers are embedded. Four small
Should-fix items remain: the α-corner wording, the T18 fragment, literal `**` in two slip callouts, and three verbose
"why" fields. Re-execute before publishing so the two explainers saved after this run are embedded fresh.
