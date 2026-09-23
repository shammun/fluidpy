# Chapter 3 — lesson review (round 2)                           2026-09-23

History: round 1 (2026-09-23): FAIL, with 7 Must-fix and 15 Should-fix items. The Must-fix items were the bare (3.24)/(3.26) in D18, the broken LaTeX of (3.32) in cell 392, the incomplete ψ gloss, the false far-field claim in C05, the false dye/particle reading in C04, and the non-existent checks cited by D13 and D19.

Reviewed: `outputs/ch03/executed.ipynb` (413 cells, 0 execution errors), built by the revised `notebooks/build_ch03.py`. I diffed the round-1 and round-2 dumps cell by cell (`outputs/ch03/nb_diff.txt`), re-read every changed cell and re-ran the whole-notebook scans.
coverage_check: **OK** (0 errors, 63 warnings: ledger primers matched by name, and explainers not built yet) · sections 6/6 · CORE blocks 15/15 with code + visual · 24/24 derivations.
Scans: no control characters in any cell; no ",." / ",," / "— —" left over from the equation injector; no bare equation numbers outside headings whose equation follows in the same cell.

## Round-1 Must-fix items — all resolved

| # | item | status |
|---|---|---|
| 1 | D18 bare "(3.24)/(3.26)" | fixed. The Tools line and step 4 now write $\Gamma=2\pi ru_\theta=2\pi r^2\omega_0$ (3.24) and $\Gamma=2\pi ru_\theta=2\pi B$ (3.26). Both match the rendered book page (printed p. 84). |
| 2 | cell 392: (3.32) broken by a form-feed character | fixed. The text reads `\Delta t\frac{…}` and no control characters are left anywhere. |
| 3 | ψ gloss | fixed. The gloss now gives $u=\partial\psi/\partial y$, $v=-\partial\psi/\partial x$ and shows $\mathbf u\cdot\nabla\psi=0$ (so ψ is constant along a streamline) and that adding $c\,y$ adds $u=c$. A sympy cell prints `U c`. |
| 4 | C05 "far-field streamlines point backwards" | fixed. The text now says the far fluid still moves toward +x, faster than U, and the cylinder drifts toward +x. |
| 5 | C04 dye/particle reading | fixed. The tagged particle's circle is directly above the port; the dye circle's centre goes round the port and meets the particle's circle at T/2. |
| 6 | D13 check cited a numeric curl that did not exist | fixed. A new code cell computes the numerical curl: ω = 1.88, ω′ = 0.48 = ω − 2Ω, with an assert (checked by hand). |
| 7 | D19 check cited a sympy line that did not exist | fixed. A new sympy cell prints `0` and `Gamma`. |

All 15 round-1 Should-fix items were also done:
- An idiom gloss was added to cell 8.
- The formulas for the developing pipe profile and the tanh front are now written out.
- "E3"/"E7" are replaced by the explainers' descriptions.
- "(analysis §8)" is removed.
- `material_volume_rate` is now run (both values 1.2566).
- The tiny example now gives the true 1.4° turn and explains why the net rate is γ/4.
- The `Gamma=` keyword is explained in a comment.
- N29 now states $u_\theta=B/r$.
- R05 was moved before C08.
- D10 step 5 now states the "axes coincide" assumption.
- D08 is split into 8 steps and D24 into 8 steps.
- There are 7 "Why interactive" sentences.
- Every same-cell equation reference now has its equation written.

## Must fix

1. **Step pointers left stale by the renumbering** (the ch01 lesson says to check them by script after any renumbering):
   - **Cell 390 (D24), step 8:** the title "Equate steps 4 and 6, cancel" should be "Equate steps 5 and 7". Step 5 is now the shrunk left side and step 7 the inserted swelling rate; step 4 is only the right side.
   - **Cell 390 (D24), *Check it*:** "F = 1: steps 4 and 6 both give…" should be "steps 5 and 7". "the left side of step 4 is 0" should be "step 5".
   - **Cell 199 (C07 code explanation):** "the antisymmetric part drops out of any n·G·n, D08 step 6" should be "D08 step 8". D08 step 6 is now "Differentiate the squared length".
   - D11 step 2's pointer to "D08 step 5" is still correct.

## Should fix

1. **Cell 394 ("What does the code above do?")** — the new item 7 (D23 on a sphere) is listed before item 6. Renumber, or put it last.
2. **Cells 40 and 395** each show an extra empty `stream` output after the real one. It is harmless but looks odd on the page; check where the stray output comes from.

## Unexplained-first-use list

All rows marked missing or late in round 1 are now explained at or before first use:
- ψ: cell 162.
- B: cell 258.
- Explainer names replace "E3"/"E7".
- Front formula: new cell before the front figure.
- Pipe profile formula: new cell before the pipe figure.
- `Gamma=`: comment in the C08 code cell.
- (3.15): R05 now comes before D10.
- Python/sympy idioms: gloss in cell 8.

No new missing rows.

## Derivation audit (round 2: changed derivations re-checked)

| D | steps | every step follows? | why / in-words ok | check ok | verdict |
|---|---|---|---|---|---|
| D08 | 8 | yes: step 6 is $\ell\dot\ell=\delta\mathbf x\cdot\mathbf G\cdot\delta\mathbf x$, step 7 divides by ℓ², step 8 uses $\mathbf n\cdot\mathbf R\cdot\mathbf n=0$ | yes | yes | pass (cell 199's pointer, Must 1) |
| D10 | 5 | yes | yes; step 5 now states the "axes coincide" assumption | yes | pass |
| D13 | 6 | yes | yes | yes (new numeric-curl cell) | pass |
| D18 | 5 | yes | yes; (3.24)/(3.26) shown and match the page | yes | pass |
| D19 | 7 | yes | yes | yes (new sympy cell: 0, Γ) | pass |
| D20 | 7 | yes | yes | yes (private reference removed) | pass |
| D23 | 5 | yes | yes | yes (`material_volume_rate` now run) | pass |
| D24 | 8 | yes (step 4 = right side, step 5 = left side) | step 8 and *Check* cite the old step numbers | yes otherwise | **fix (Must 1)** |
| others | — | unchanged since round 1; only the injector punctuation was cleaned | yes | yes | pass |

## What works well

Everything from round 1 still holds: honest callouts for the book's typos, γ / Γ / circulation kept apart, two independent routes to every key number, orders of smallness shown as fitted slopes, consistent colours, and the observer-dependence thread with its climate hooks. New strengths in this round:
- The ψ gloss proves its two claims in two lines of sympy.
- The D13 check uses a nonlinear field, not just the formula's arithmetic.
- The 1.4° tiny example now explains *why* the ellipse turns at about γ/4 and not γ/2 — a good "measurable approximation" lesson for the teaching-style skill.

## Verdict: FAIL

The only remaining Must-fix is the set of step pointers left stale by renumbering D08 and D24 (cells 199 and 390): three text edits in `notebooks/build_ch03.py`. After them the lesson passes; no re-review is needed beyond a grep for "step [0-9]" pointers.

## Round 3 (orchestrator, 2026-09-23)
Round-2 Must-fix (stale D08/D24 step pointers) applied by the builder; verified by grep of notebooks/build_ch03.py (D24 step 8 "Equate steps 5 and 7", check "steps 5 and 7" / "step 5", C07 note "D08 step 8"; no other pointers into D08/D24). Coverage OK for ch03 (0 errors, 63 warnings), 413 cells, 129 s. Should-fix: list order and stray outputs fixed; one stray "
" stream output in the D22 sympy cell remains (kernel stream split, cosmetic).

**Verdict: PASS**
