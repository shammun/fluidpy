# Chapter 13 — lesson review, round 3 (final)                   2026-10-08

## Verdict: PASS — 0 Must fix, 3 Should fix (wording only; none blocks publishing)

coverage_check: OK (0 errors, 2 warnings — "explainer not built yet" for D25 and D26 only) · sections covered 18/18 ·
CORE blocks 17/17 complete (plain words, idea, derivation, tiny example, from-scratch code, at least 2 visuals each,
see / read / change) · derivations 29 (253 steps) · primers P307–P335 present, none cited before its cell · 821 cells,
0 execution errors, no code cell without output, no stderr, no symbol-initial sentence, no equation number without its
equation · `check_public.py` clean on the builder, the notebook and this report.

Reviewed file: `outputs/ch13/executed.ipynb` (0-based cell indices of the 821-cell notebook).

---

## Must fix

None.

## Should fix

1. **Parcel wording, cells 127, 128, 129 (tools line), 150, 153 against the figure note of cell 152.** The figure and its
   note are right: a frictionless parcel released from rest in a uniform pressure gradient traces **cycloid arches with
   cusps** — it stops for an instant at each cusp and its path never closes into a loop (checked on the image; and
   analytically, the velocity is the balanced value times 1 − e^{−ift}, which returns to zero once per inertial period).
   The surrounding prose still says the air "overshoots and loops", "loops about the balanced drift", "removes the loops",
   "it circles". The physics meant is the same (the *velocity* circles about the balanced value) and nothing stated is
   used wrongly, but a first-time reader meets "loops" in the text and "arches, not loops" under the figure. One sentence
   in cell 152 would reconcile them — "these arches are the 'loops' of the text: the velocity vector circles the balanced
   velocity; the *path* closes into true loops only if the parcel starts with extra speed (as in the explainer's circular
   low)" — or replace "loops" by "swings in arches" in 127/128/150/153.
2. **Cell 1, explainer table, row 1:** "the Coriolis force turns it until the only motion left runs along the isobars" is
   the last remnant of the settling picture; "…turns it; on average the motion runs along the isobars" matches cells
   127–152.
3. **Figure of cell 45 (ocean profile):** the label "thermocline" sits inside the amber band at the top-left where the band
   is only a few pixels tall at this figure size; the note is correct, the label is hard to see. Move it beside the band.

## Rounds

| round | Must fix | Should fix | what changed |
|---|---|---|---|
| 1 | 15 | 19 | parcel "settles"; mirrored Kelvin sketch; three figure notes describing other figures; one ratio and two rounded multiples not ours; pendulum day; a wrong *why* in D26; mode boundary conditions after D11; unexplained NaN and WKB; slip #13 box missing; an unsupported attribution; x–t prediction reversed |
| 2 | 2 | 9 | all 15 verified fixed, no regression; two more figure-note mismatches found in figures first viewed in round 2 (Poincaré dispersion panel, inertia–gravity hodographs) |
| 3 | 0 | 3 | both verified fixed; all nine round-2 Should-fix items verified applied; eight rewritten figure notes checked against their images |

Round-3 verification of the round-2 items (new text read on the executed notebook, figure opened beside it):

| item | status | cells |
|---|---|---|
| Poincaré dispersion: N87, code comment, note | fixed — one log–log panel in units of abs(f) and Λ, three curves, "slope 1 means non-dispersive"; matches the image | 457–459 |
| Inertia–gravity figure: N122, comment, note | fixed — two hodographs, sense of rotation per hemisphere read from the phase marks, axis ratio as plotted; c_p ⟂ c_g referred to N121/D21 | 636–638 |
| slip #13 box wording | fixed — spelling only, attribution not asserted, no date | 717 |
| `sech` before its definition | fixed — defined inside D24's Check | 705 |
| D09 step 6 symbol clash | fixed — complex velocity renamed | 251 |
| "several times less" | fixed — "about half as much"; panel titles neutral, same scale | 578–579 |
| eddy re-entering the periodic box | fixed — run re-timed to one year, periodicity stated, clean background | 696–697 |
| over-long arrows in the four-centre figure | fixed | 148 |
| first "What to try" bullet | fixed — overshoot; settles only with drag | 153 |
| D11 tools list | fixed — boundary conditions and P315 listed | 353 |
| caption touching the drag curve | fixed — moved to a figure title | 151 |

Figure notes rewritten by the builder after opening the images — each compared with its image:

| figure | cell | note agrees with the image? |
|---|---|---|
| parcel paths | 151–152 | yes — arches with cusps along the isobar, mirror image at 35° S, straight drag path toward low pressure |
| ocean profile | 45–46 | yes — fastest change just under the mixed layer, N largest at the top of the thermocline (label placement: Should fix 3) |
| Taylor column sketch | 184–185 | yes — side view, obstacle, shaded column, equal arrows at three heights, Ω |
| potential-vorticity columns | 545–546 | yes — three columns, counter-clockwise on the stretched one, clockwise on the squashed one, stated as f > 0 |
| wavenumber-5 map | 654–655 | yes — dashed contours on the poleward (low) side, five waves, wind along the contours |
| Eady wedge | 724–725 | yes — sloping surfaces, shaded wedge, two arrows flatter than the surfaces |
| cascade sketch | 788–789 | yes — two branches, stirring scale, opposite arrows |
| drifting eddy | 696–697 | yes — westward drift, weakening, short waves to the east |

## Unexplained-first-use list

Empty: every term, symbol, tool and function sampled in three rounds is explained at or before its first use (the round-1
and round-2 rows — mode boundary conditions, WKB, NaN, Courant number, Cramer's rule, hPa, spin-down time, sech — are all
closed).

## Derivation audit

All 29 derivations were re-worked step by step in round 1; the edited steps were re-read in rounds 2 and 3 (D02 step 6,
D09 step 6, D11 step 13 and tools, D13 step 6, D16 Check, D18, D24 Check, D26 step 8, D27 steps 3 and 5). Every step
follows from the previous one by the stated move; *why* and *in words* are correct; each has a check and a meaning; all
nine ★★★ derivations (D06, D11, D12, D16, D17, D19, D22, D25, D26) carry a sympy cell that rebuilds the steps.
Verdict: ok for D01–D29.

## What was verified across the three rounds

- Equations against the page images: 70 numbered displays on 24 pages; signs and factors agree; (13.112) only in its
  printed form; the one flagged departure is the cut-off inequality written with the symbol for the constant.
- Lapse rate in both conventions wherever it appears; hemisphere statements for geostrophy, thermal wind, Ekman spiral,
  transport and pumping, bottom layer, Kelvin wave (also run against `GFD.kelvin_wave`), adjustment jet, potential
  vorticity and flow over a step.
- Slips (14 rows, two loose, each with a box) and traps (18) kept apart; orthogonality with weight 1 for both lids;
  rigid-lid index; bottom-layer overshoot; the cubic's caveat; turbulence statements qualitative; "ours" labels.
- Public-repo rule: no forbidden string or pattern; no worked example pairs a book-quoted input set with a book-quoted
  result; the three places that once echoed book-quoted values now use our computed ones.
- Figures: 38 of 43 viewed over the three rounds; five never opened (cells 15, 35, 228, 308, 809; their
  notes were read and are plausible from the code, but the images were not compared). Reading notes of all viewed figures
  agree with the images.
- Code-explanation cells of C01–C17 read against the outputs above them; quoted numbers match.

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- Two routes to one result with the reason the second is stronger (D07); the conserved quantity shown not moving before
  it is derived (535); "period mean versus snapshot" with measured numbers (532–534).
- "Honest caveat" boxes with a criterion and computed cases (435); slip boxes backed by an engine that fails with the
  printed form (56, 176); boundary conditions signposted just before the derivation that needs them (346).
- ASCII sections that name the viewing direction and give the mirror rule for the other hemisphere (481).
- **Lesson for the pipeline:** a reading note must be written from the rendered figure, not from the storyboard. Six
  Must-fix items in rounds 1–2 were notes describing a planned figure that the figure function does not draw; the
  builder's round-3 habit of opening each PNG before writing its note removed the whole class. Worth a line in the
  `teaching-style` and `python-viz` skills, and a guard in the builder (panel count and axis scale asserted against the
  words "inset", "upper/lower", "left/right", "slope").
