# Chapter 16 — lesson review (round 2)                                   2026-10-14

**Resolution (orchestrator, 2026-10-14): PASS.** M1 fixed: a 🔁 recap cell before the steepening figure derives the
tube's Riemann invariant ($c=c_0(A_0/A)^{1/4}$ ⇒ $\int c\,dA/A=-4c$ ⇒ $u+4c$ on $dx/dt=u-c$), with the figure note and the
code comment made consistent (checked by the orchestrator). Should-fix 1–6 applied (headings in C01, C04, C07, C09–C11,
C13–C15, C17; P375 `jv` limit α ≈ 1000 and α ≈ 30 labelled our habit; D15 "about 0.1 % (0.14 %)"; p_e gloss in the
conventions table; c080 label; comment leftovers). Executed: 577 cells, 133 s, 0 errors; coverage 0/0; check_public OK.

coverage_check: OK (0 errors, 0 warnings) · sections covered 4/4 (16.1–16.4, plus S01/S02) · CORE blocks 17/17 with
code + visual · derivations 34/34 present, 8/8 ★★★ with an executed sympy cell · explainers 10/10 embedded (all ten viz
files are older than the run, so no embed is stale) · 0 execution errors.

Reviewed: `outputs/ch16/executed.ipynb` (565 cells, 15:04), built from `notebooks/build_ch16.py` (14:56).
- I dumped it to text and diffed it, line by line, against the round-1 dump (kept in `outputs/ch16/review_r1/`).
- I re-read every changed passage and every changed derivation field.
- I looked at the 23 regenerated PNGs.
- I recomputed the corrected numbers with the module and with `outputs/ch16/scratch_checks.py`.
- I rendered and compared printed pp. 804, 813, 816, 828 and 830 (`--eq 16.45`, `16.94`, `16.112`, `16.166`, `16.179`).

## Verdict: PASS (with one small Must fix)

All 11 round-1 Must-fix items are resolved and correct. Every new equation display matches the rendered page. Every
corrected number agrees with the module.

One new Must-fix comes from the rewritten steepening note (cell 252). It now rests on a Riemann invariant, $u+4c$, but
gives no reminder and no reason for the 4. It is a one-sentence text fix, so the verdict is PASS on the condition that
M1 is applied. If you want a strict gate, treat this as FAIL until that sentence is in. Nothing else blocks.

## Round-1 status

| R1 | Item | Status in round 2 | Evidence |
|---|---|---|---|
| M1 | (16.31) displayed as "$Q_p+Q_c\neq$" | **fixed** | All 8 places now show $Q=Q_p+Q_c=\frac{\pi a^4\Delta p}{8\mu_pL}[1-(1-\delta/a)^4(1-\mu_p/\mu_c)]$. (16.179) in the conventions table and the N92 table is now the book's $D=\sqrt{2\delta}\,\hat Ga^2a/(\mu\nu)$, which matches p. 830. |
| M2 | D06 printed − corrected sign | **fixed** | Cells 138 and 141 now read $+\pi\Delta p\,a^2\delta(2a-\delta)/(4\mu_pL)>0$ and say the printed form overstates the core flux by 75 %, which matches the +0.753 the code prints. I re-derived $a^2(a^2-b^2)=a^2\delta(2a-\delta)$. |
| M3 | D15 "5 % below at 0.9a" | **fixed** | The text now gives 0.1 % in amplitude, 1.3 % complex, and 2 % at 0.97a. My values: 0.14 %, 1.34 %, 2.1 % (see S3). |
| M4 | Steepening mechanism | **fixed** | The note now gives $c^2=\frac{Eh}{2\rho a_0}(A_0/A)^{1/2}$ falling with A and $u+c$ rising: crest 6.47 m/s, foot 4.86 m/s. It names the ringing and has a pure right-going start. Figure c251 shows no left-going bump. The new wording needs M1 below. |
| M5 | "half as tall" → a quarter | **fixed** | The note gives the radius change $a_0^2(p-p_e)/(Eh)\propto1/E$; (16.41) linearised on p. 804 confirms it. |
| M6 | "four orders" | **fixed** | The note gives Re 6.1 decades and α 3.4 decades, and a new code line prints both. |
| M7 | "10¹⁵" | **fixed** | The note gives 6 × 10¹³, printed as 5.72e+13. |
| M8 | Murray "30 % either side" | **fixed** | The note gives +71 % / +24 %, printed by `murray_cost` and by the ⅓x⁻⁴ + ⅔x² formula. |
| M9 | Casson 3.6 % | **fixed** | The note gives 4.7 %; the printed bracket is 0.0474. |
| M10 | Bare equation numbers | **fixed** | (16.45) is written as printed and corrected, and (16.46) is added; I checked both against p. 804 and the integration. (16.94) is written as printed and in single-mode form (p. 813). (16.112) and (16.113) are written in full (p. 816). (16.166)–(16.170) appear as a display block in N86, and the slip-S34 box names $-w^2\cos\alpha$ and $+w^2\sin\alpha$ (p. 828). (16.180) appears as printed in N92 with the magnitude reading; p. 830 says "−Ĝ is the pressure gradient", which agrees. |
| M11 | p_e, transmural, a₀ | **fixed** | Cell 161 defines all three before D08, and cell 162 adds "The idea" with an ASCII sketch. |
| S1 | sympy inline in Check | **fixed** | All 8 ★★★ Check fields now point to the cell below. |
| S2 | Pipeline references | **fixed** | None of the 11 phrases remains. The P-numbers of primers remain, as a reader-facing convention. |
| S3 | Comments against output | **fixed** | Cell 70: 0.00774 / 0.00466. Cell 267: 1.21e11 / 0.0632. |
| S4 | Placeholder comments | **fixed** | 0 × "intermediate quantity", 0 × "loop over the cases", 0 × "draw the curve". The RK4, Euler, `quad`, `lapf` and `live` lines now carry their meaning and unit. A trivial residue is in S6. |
| S5 | Casson "stay on K_c²" | **fixed** | The note gives 19 % / 8 % and 41 % / 18 %, printed. |
| S6 | Tracer shift ~1 cm | **fixed** | The note gives 1.9 cm, printed. |
| S7 | Gauge / absolute wording | **fixed** | D34 step 5 and cell 546 give gauge −0.1 MPa at the ground and absolute tension from about 13 cm up. I checked it: 1325 Pa / (998 × 9.81) = 0.135 m. |
| S8 | Figure cosmetics | **fixed** | c034, c080, c206, c413, c449 and c469 are clean. One new small overlap is in S5. |
| S9 | Glosses | **fixed** | Diastole (N04, cell 34 comment), Kelvin–Voigt (N72), embolism (D34), `np.searchsorted` (P395 primer, cells 211–212, before its first use in 213), `np.unravel_index` (comment and reading note), `try/except` (cell 59, right after its first use), `np.arctan2` (comment, Ch. 2 P70) and `lambda` (Ch. 1 P29). |
| S10 | "The idea" headings | **open** (as announced) | C06 gained one. C01, C04, C09, C10, C11, C13, C14, C15 and C17 still lack one. C07 and C15 still lack "The problem in plain words". |
| S11 | "influence coefficients" | **fixed** | Replaced by "differential relations with a 1 − S² factor". |
| S12 | Stale embeds | **resolved** | All 10 embedded; `xylem_water_potential` is now built. Every viz file predates the 15:04 run. |

## Must fix

1. **Cell 252 (steepening, C08) uses a Riemann invariant without a reminder, and the "4" cannot be checked.** The note
   says the start is set so that "$u+4c$ is the same everywhere ($u=4(c_0-c)$, the Riemann invariant of this tube law)".
   - Riemann invariants were primed in Ch. 15 (P358), but only for a gas, as $u\pm2c/(\gamma-1)$.
   - Nothing here says where the factor 4 comes from, so the reader cannot check it.
   - The physics is right. I checked it: with $c=c_0(A_0/A)^{1/4}$ one gets $\int c\,dA/A=-4c$, so $u+4c$ is the quantity
     carried along the left-going characteristics $dx/dt=u-c$. Making it uniform removes the left-going wave.

   Add one reminder sentence, e.g.:

   > 🔁 Riemann invariant (Ch. 15 P358): $u\pm\int c\,dA/A$ keeps its value along the characteristics $dx/dt=u\pm c$ (Ch. 7 P174); here $c\propto A^{-1/4}$ gives $\int c\,dA/A=-4c$, so $u+4c$ rides the left-going ones — make it uniform and no left-going wave exists (the gas analogue is $u+2c/(\gamma-1)$).

   Put it either in a `🔁` line before cell 251 or in the "What you see" paragraph.

## Should fix

1. **Heading consistency (round-1 S10, still open).** Add a one-line "#### The idea" (or an ASCII sketch) to C01, C04,
   C09, C10, C11, C13, C14, C15 and C17. Add "The problem in plain words" to C07 and C15.
2. **Primer P375 (cell 266) overstates when `jv` overflows.** It says "beyond α ≈ 30 use `jve`". In fact `jv(0, i^{3/2}α)`
   is finite and accurate at α = 100 (≈ 2 × 10²⁹) and α = 500 (≈ 6 × 10¹⁵¹), and it becomes `inf` only near α ≈ 1000,
   where $e^{\alpha/\sqrt2}$ passes 10³⁰⁸. Write "overflows near α ≈ 1000; we use `jve` and ratios as a safe habit".
   Alternatively, keep α ≈ 30 and say it is a conservative choice of ours.
3. **D15 Check (cell ~281): "within 0.1 %".** My recomputation gives 0.14 %, which is not within 0.1 %. Write "about
   0.1 %". The "2 % at 0.97a" figure is the complex gap (2.1 %); the amplitude gap there is 1.1 %. Say "complex" so the
   reader knows which gap grows.
4. **The conventions table (cell 7, rows viii and α) uses $p-p_e$ before C06 defines p_e.** It is a reference table, so
   this is low-risk. Add "(p_e: pressure outside the vessel, C06)" to row viii.
5. **New small overlap in c080 (left panel).** The label "intercept √τ_y" sits on the blue Casson line near
   √γ̇ ≈ 4–7. Move `xytext` to about (6, 0.1), below the line.
6. **Comment residue.** "# a label with an arrow" sits on `annotate` calls that draw no arrow: cell 34 (corner names) and
   cell 120 (vessel names). "# (continued)" sits on the first line of the `np.where` profile in the C16 slider cell. Both
   are trivial.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / tool | First used | Explained |
|---|---|---|
| Riemann invariant $u+4c$ for the tube law | 251 (code comment), 252 | **missing**: no reminder of Ch. 15 P358 and no derivation of the 4 (M1) |
| p_e in the conventions table | 7 | 161 (C06). Acceptable in a reference table; a gloss would help (S4) |
| characteristic speed u + c | 251–252 | explained in place in 252 ("what a point of the pulse travels at"); Ch. 7 P174 ✓ |
| MacCormack's scheme / numerical ringing | 250, 252 | named as Ch. 10 in N54 ✓ |
| p_e, transmural pressure, a₀ | 159 (heading), 161 | 161, before D08 ✓ (round-1 M11 closed) |
| diastole | 34 | N04 cell before it and the cell 34 comment ✓ |
| Kelvin–Voigt wall | N72 | in place ✓ |
| embolism | D34 | in place ✓ |
| `np.searchsorted` | 213 | P395 primer, 211–212 ✓ |
| `np.unravel_index` | 449 | comment and reading note ✓ |
| `try` / `except ValueError` | 58 | 59 ✓ |
| `np.arctan2`, `np.hypot` | 275 | comment, Ch. 2 P70 ✓ |
| `lambda` | 70 | comment, Ch. 1 P29 ✓ |
| `DataFrame.set_index(...).loc[...]` | 98 | Ch. 15 P361 (pandas); acceptable for this reader |
| Bessel functions, `jv`/`jve`, Young's modulus, Laplace's law, integrating factor, convolution, perturbation series, `minimize`, osmosis | as in round 1 | primers before use ✓ (unchanged) |

## Derivation audit (changed derivations only; D id · steps · every step follows? · why/in-words ok · check ok · verdict)

No derivation line (the maths) changed. Only the Tools, Start, *why*, Check and What-it-means texts changed. I re-read
each of those texts against its step.

| D | Steps | What changed | Follows? | Why / in words | Check | Verdict |
|---|---|---|---|---|---|---|
| D01 | 6 | Tools (Ch. 1 / Ch. 5 reminders); Check ("the `heart` inputs above") | yes | ok | 0.9999 J, 1.0666 J ✓ | ok |
| D02 | 5 | Tools: $(x+y)^2$ written out | yes | ok | ✓ | ok |
| D03 | 4 | Tools: Stokes-layer recap (cell 91 is above) | yes | ok | ✓ | ok |
| D06 | 8 | Header, step 6 *why*, Result and Check now show the true (16.31); the sign is now $+\delta(2a-\delta)$ | yes (re-derived $a^4-a^2b^2$) | ok | ✓, consistent with +0.753 | ok |
| D07 | 6 | Tools (binomial written out); Start shows the true (16.31) | yes (step 5 identity re-expanded) | ok | 1.969 / 1.818 ✓ | ok |
| D10 | 9 | Check names `windkessel_summary` instead of "C.12" | yes | ok | ✓ | ok |
| D13, D19, D20, D21, D24, D27, D30, D33 | 13, 10, 10, 13, 11, 14, 15, 10 | Check text no longer inlines the code; it points to the executed cell (outputs unchanged: `0`, `1.3j` × 3, 25.7057°, 1541/50400) | yes | ok | ✓ | ok |
| D15 | 8 | Check: 0.1 % / 1.3 % / 2 % | yes | ok | true to rounding; S3 asks for "about" and "complex" | ok |
| D27 | 14 | Step 11 *why* ("the fastest fluid moves to the outer side"); What it means names the ¼ threshold | yes | ok | consistent with `threshold=0.25` in cell 449's code and with K ≈ 460 | ok |
| D30 | 15 | Steps 3 and 6 *why* (no pipeline references) | yes | ok | ✓ | ok |
| D34 | 5 | Step 5 *why* (gauge against absolute, 13 cm); What it means (embolism) | yes | ok; I checked 0.135 m | ✓ | ok |

The other 18 derivations are unchanged from round 1, where every step followed.

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **Every corrected adjective now carries a printing cell next to it:** decades, the 6 × 10¹³ resistance ratio, +71 % /
  +24 %, 19 % / 8 %, 1.9 cm, 4.7 %, and the steepening speeds. The reader can see the claim and its source together.
  Recommend this for teaching-style: "a quantitative word in a reading note gets a `print` in the cell above it".
- **The steepening block is now honest about numerics.** It separates the physics (the $u+c$ ordering, the real decay
  after breaking) from the artefacts (overshoot, wiggles, the dip ahead of the front), and it starts the run from an
  exact simple wave, so no spurious left-going wave appears. Once M1 is added, this is a good pattern for any shock-capturing
  demo.
- **Each slip is shown in both forms.** (16.45) and (16.94) appear as printed and as corrected, (16.180) as printed with
  the magnitude reading, and (16.46) is confirmed correct. The ch01 "which p_o?" pattern is applied to printed equations.
- **The definitions cell before D08 ("Three symbols for the rest of the chapter") is a reusable fix for symbols that
  the book introduces silently.**

## Verdict: PASS (conditional on M1)

- Round 1: all 11 Must-fix items are fixed and verified against rendered pages and module output. Ten of the 12
  Should-fix items are fixed; S10 (headings) is still open; S12 (stale embeds) is resolved.
- Remaining: **one Must-fix (M1)**, a one-sentence Riemann-invariant reminder in cell 252, plus six minor Should-fix
  items.
- Coverage is 0 errors and 0 warnings, 17/17 CORE blocks have code and a visual, all 10 explainers are embedded, and
  the run has no execution errors.
