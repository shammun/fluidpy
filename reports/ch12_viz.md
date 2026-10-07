# Chapter 12 — explainer review, round 2                          2026-10-07

Reviewer: `viz-reviewer` (read-only on `viz/`). Round 2 re-reviewed the five explainers that failed round 1, as they
stood at 15:40–15:56 local time (file times 14:38–15:14); the other five are unchanged since round 1 (file times
11:54–13:41). Storyboards: `analysis/ch12_design.md` Part B (E1–E10), derivations Part F.

## Verdict: PASS — all ten explainers PASS

All six round-1 Must-fix items are fixed, each confirmed in the state where it failed. No new Must fix. The lists of
Should-fix items below are polish and do not block the gate.

| slug | lint | shot (sizes/views) | parity rows ok | KaTeX errors (my probe) | explain | derivations (steps ok) | verdict |
|---|---|---|---|---|---|---|---|
| reynolds_averaging_window (E1) | ok | PASS 8 / 224 | 26/26 (16 py) | 0 in 52 state loads, 706 tab/page visits | 8 sections, live | D01 8/8 · D03 5/5 | **PASS** (round 2) |
| correlation_and_spectrum (E2) | ok | PASS 8 / 312 | 33/33 (26 py) | 0 in 33 states, 245 fragments | 9 sections | D03 5/5 · D04 7/7 · D05 9/9 | **PASS** |
| reynolds_stress_parcels (E3) | ok | PASS 8 / 208 | 32/32 (20 py) | 0 in 33 states, 340 fragments | 8 sections | D06 13/13 | **PASS** |
| energy_cascade_spectrum (E4) | ok | PASS 8 / 272 | 40/40 (30 py) | 0 in 36 state loads, 510 visits | 9 sections | D11 9/9 · D12 9/9 | **PASS** (round 2) |
| turbulent_energy_budget (E5) | ok | PASS 8 / 320 | 56/56 (41 py) | 0 in 42 state loads, 600 visits | 8 sections | D09 9/9 · D10 15/15 | **PASS** (round 2) |
| turbulent_jet_similarity (E6) | ok | PASS 8 / 528 | 97/97 (67 py) | 0 in 42 state loads, 972 visits | 9 sections | D13 10/10 · D14 14/14 · D15 10/10 · D16 12/12 | **PASS** (round 2) |
| law_of_the_wall (E7) | ok | PASS 8 / 432 | 72/72 (50 py) | 0 in 49 states, 671 fragments | 9 sections | D17 9/9 · D18 7/7 · D19 12/12 · D20 5/5 | **PASS** |
| mixing_length_closure (E8) | ok | PASS 8 / 192 | 54/54 (38 py) | 0 in 31 states, 441 fragments | 8 sections | D21 11/11 | **PASS** |
| stratified_surface_layer (E9) | ok | PASS 8 / 256 | 151/151 (134 py) | 0 in 30 states, 775 fragments | 9 sections | D24 8/8 · D25 7/7 | **PASS** |
| taylor_dispersion (E10) | ok | PASS 8 / 352 | 36/36 (22 py) | 0 in 46 state loads, 836 visits | 8 sections | D26 11/11 · D27 8/8 · D28 7/7 | **PASS** (round 2) |

Totals: 597 selftest rows, 597 ok; 24 derivation instances (220 steps), all with the step count and move titles of
Part F.

## Mechanical results, round 2 (my own runs)

`tools/viz_lint.py` on the five: `ok` each. `tools/check_public.py viz/ch12/<slug>.html` on the five: `public-repo check
OK (1 files)` each.

`tools/eq_refs.py` on the five: 0 hits in `reynolds_averaging_window`, `turbulent_energy_budget` and
`turbulent_jet_similarity`; `energy_cascade_spectrum` 3 (L2947 and L2999 "Eq. (12.54), exponent corrected" — `ref:`
labels on a card and a derivation that show the TeX; L3056 a fragment of a `py:` expression); `taylor_dispersion` 1
(L2980 "Eq. (12.129), condition corrected", a `ref:` label on a card that shows the TeX). All acceptable, as in round 1.

`tools/shot.py viz/ch12/<slug>.html`, full, one at a time, 15:40:21–15:51:10:
```
PASS  viz/ch12/reynolds_averaging_window.html  (8 sizes, 224 views, 6 steps, 26 selftest rows)
PASS  viz/ch12/energy_cascade_spectrum.html  (8 sizes, 272 views, 7 steps, 40 selftest rows)
PASS  viz/ch12/turbulent_energy_budget.html  (8 sizes, 320 views, 7 steps, 56 selftest rows)
PASS  viz/ch12/turbulent_jet_similarity.html  (8 sizes, 528 views, 7 steps, 97 selftest rows)
PASS  viz/ch12/taylor_dispersion.html  (8 sizes, 352 views, 7 steps, 36 selftest rows)
```

`tools/shot.py --chapter ch12 --quick` (parity sweep over all ten, 4 sizes), run last:
```
PASS  viz/ch12/correlation_and_spectrum.html  (4 sizes, 156 views, 7 steps, 33 selftest rows)
PASS  viz/ch12/energy_cascade_spectrum.html  (4 sizes, 136 views, 7 steps, 40 selftest rows)
PASS  viz/ch12/law_of_the_wall.html  (4 sizes, 216 views, 8 steps, 72 selftest rows)
PASS  viz/ch12/mixing_length_closure.html  (4 sizes, 96 views, 6 steps, 54 selftest rows)
PASS  viz/ch12/reynolds_averaging_window.html  (4 sizes, 112 views, 6 steps, 26 selftest rows)
PASS  viz/ch12/reynolds_stress_parcels.html  (4 sizes, 104 views, 6 steps, 32 selftest rows)
PASS  viz/ch12/stratified_surface_layer.html  (4 sizes, 128 views, 8 steps, 151 selftest rows)
PASS  viz/ch12/taylor_dispersion.html  (4 sizes, 176 views, 7 steps, 36 selftest rows)
PASS  viz/ch12/turbulent_energy_budget.html  (4 sizes, 160 views, 7 steps, 56 selftest rows)
PASS  viz/ch12/turbulent_jet_similarity.html  (4 sizes, 264 views, 7 steps, 97 selftest rows)
```
10 of 10 PASS, 0 failures; 597/597 selftest rows ok, 444 of them JS ↔ fluidpy; the worst relative error of a parity row is 8.0e-07 (law_of_the_wall: indicator y+ dU+/dy+ (exact derivative against differences)).

`tools/embed_check.py ch12` after the notebook rebuild (690 cells): `embed check OK for ch12` (exit 0).

**My own KaTeX probes** (scripts in my session scratchpad, written by me, not the builders' probes):
- *DOM probe* on the five fixed files at 1366×768 and 360×640: for the default state, every preset, every quiz state,
  both ends of every slider, every option of every chip control and every mode, open each tab (Explore, Explain,
  Equations, Code, Check), every derivation page in the default state and four pages per derivation in the others, and
  every walkthrough step, and count `.katex-error` nodes in the DOM (pager pages that are not on screen are in the DOM
  and are counted). Result: **0 nodes** in 218 state loads and 3 624 tab/page visits.
- *String probe* on all ten: every live TeX fragment (Explain, notes, status, inspector, equation and derivation live
  lines, quiz answers, view titles) in every preset, step, quiz and derivation state and at both ends of every slider,
  rendered with KaTeX in throwing mode: 335 states, about 4 240 distinct fragments, **0 failures** (the two hits in
  `mixing_length_closure` are my tag-stripper cutting a `$350<y^+<…$`; the page renders it, checked in the DOM in round 1).
- The same selector found the two round-1 failures, so the probe can see what it is looking for.
Limits: sliders are moved one at a time, not in combination.

Other round-2 passes: all 11 Code blocks of the five files run as displayed in three or four states each; config strings
re-scanned for control characters and TeX that lost a backslash (0 hits); derivations re-compared with Part F by script
(no deviation left); longest derivation *why* now 36 words (E6 D16 step 1), all others ≤ 35.

## Round-1 Must fixes → status

| # | round-1 item | status | evidence |
|---|---|---|---|
| E1-1 | Explain 4: boxed rms-error line failed in KaTeX at the default state | **fixed** | the line is now `\sqrt{0.226^{2}+\big(3.27\times10^{-5}\big)^{2}+0.267^{2}}`; rendered in `reports/viz/ch12/reynolds_averaging_window/REVIEW2__explain_rms_default.png` (default) and `REVIEW2__explain_rms_one_period.png` (preset "one whole period": `(2.01×10⁻⁴)²`); presets "too short a window", "one whole period" and check question 2 are in the DOM probe: 0 error nodes |
| E1-2 | Code block `prod` raised `NameError` | **fixed** | the three blocks run in order in four states (default; decaying, N = 64, t = 10; stationary, N = 16, r = 0.3; wave, N = 2, r = 0): `prod` prints 2e-16, 6e-17, 1e-16, 0.0 and its `assert` holds |
| E4-1 | Explain 5: boxed −5/3 arithmetic failed in KaTeX when ε̄ was in scientific notation | **fixed** | preset ocean now reads `0.245×(1×10⁻⁴)^{2/3}×110^{−5/3} = 2.09×10⁻⁷` (`reports/viz/ch12/energy_cascade_spectrum/REVIEW2__explain_law_ocean.png`); also rendered at ΔU = 0.01 and 50 m/s and at L = 10⁴ m (`REVIEW2__explain_law_mindU/maxdU/maxL.png`) |
| E5-1 | D10 step 4 was a placeholder | **fixed** | the line is the eight-term equation in three aligned rows, identical term by term to Part F D10 step 4; steps 5–12 pick its terms up in order (time derivative, advection by the mean, mean-stress term, shear term, triple term, pressure, buoyancy, viscous); seen at 1366×768, at 1000×700 (`REVIEW2__d10_step4_notebook.png`) and at 360×640 (`REVIEW2__d10_step4_phone.png`, shown row by row) |
| E6-1 | stale typed C₆ = 2.240 | **fixed** | the symbol row now reads "1/(C₅∫HF dξ) = 1/(2.577 × 0.1771) = 2.19", built from the page's constants; selftest rows pin 2.5774, 2.1903, 0.1505, 0.2129, 0.1771; no "2.24" is left in the file |
| E10-1 | "the page prints" for the book, five places | **fixed** | no "page prints" is left; slips #11, #12 and #13 are in the house wording with the equations written out (Explain 4 and 6, D27 check, D28 step and check, card `dtl`, check question 4 "…, as the book prints?") |

## New Must fixes in round 2

None.

---

## reynolds_averaging_window (E1) — PASS (round 2)

**Must fix** — none (E1-1 and E1-2 fixed).

The new `prod` block is honest and clear: it builds its own u and v runs (two `make_ensemble` calls, v made from u with
the chosen correlation), says "N -> inf here: 0.766 + 0.054 = 0.820" for those wave-free runs, and labels the page's own
sample separately ("the page's own 8 runs: 0.638 = 0.616 + 0.022"); the block's note says the runs are wave-free and
that the two parts only estimate their N → ∞ values.

**Should fix**
- `prod`, last comment with N = 2 and r = 0 prints "0.241 = 0.280 + -0.040": write "− 0.040", and round so that the sum
  closes on the printed digits.
- In mean + wave mode the block's N → ∞ covariance (0.054, wave-free) differs from Explain 5's (0.0835, with the wave).
  The note covers it; a comment "(no wave in these runs)" on the `N -> inf` line would put it where the reader looks.
- The OU window-noise term still has no fluidpy counterpart (ruling b): `ch12.time_average_variance_ou` at the next
  implementer pass.

**What works**
- Two estimates of one mean on one clock, with the bold "so far" part and the true mean as a ghost.
- The error view decomposes the window error into noise, drift and leak, and ◆ (formula) sits beside ◇ (measured).
- "A 30-year climate normal is neither 3 years nor 300" in the regime reading.

## correlation_and_spectrum (E2) — PASS

**Must fix** — none.

**Should fix**
- `phone-tall__tour-step4`: the "All steps →" button of the quoted derivation step is cut by the page slice
  ("continued on the next page ›"). Shorten the step's quote card or drop the button from the quote on phones.
- Colours: the negative part of the product is rose and the spectrum blue (rose and blue mean dissipation and buoyancy
  elsewhere in the chapter). The storyboard asked for this; say so in Explain section 0 in one clause.

**What works**
- Dragging the copy sideways builds one point of r(τ); the equal-area rectangle makes Λ_t a picture.
- Phones show two views at a time and swap them by step (signal + correlation, then correlation + spectrum).
- "A small Λ_t does not always mean a short memory" in the hidden-wave case, with the number (0.0138 s against 2 s).

## reynolds_stress_parcels (E3) — PASS

**Must fix** — none.

**Should fix**
- The corrected flux equation carries "⚠️ The book prints the middle term with ū …" without its number. It is now
  slip #16 in `ch12.book_slips()`: use `⚠️ slip #16 — the book prints … ; the correct form is …` (ruling c).
- Blue marks the parcels with uv > 0; blue is buoyancy elsewhere. A muted grey would keep the code.

**What works**
- One clicked parcel traced to its quadrant, then 2000 of them; the estimate-against-count curve with its ± 5 s.e. band.
- Reynolds against viscous stress on a log bar: 2667 times, and 833 times more again in water.
- Every sampled number is labelled an estimate; the selftest holds the page's own sample to 5 standard errors.

## energy_cascade_spectrum (E4) — PASS (round 2)

**Must fix** — none (E4-1 fixed).

Also checked: the new sentence in Explain 5 is right — the model's cut-off factor already acts well below the
Kolmogorov wavenumber, and $S_{11}(k_1)$ is an integral of the three-dimensional spectrum over all $K\ge k_1$, so a
missing high-wavenumber tail lowers it early. `keepTiny` works: with the probe at k₁ = 1.5×10⁵ rad/m the title reads
"law 5.8×10⁻¹⁰ · model 7.22×10⁻²⁰ m³/s²" (`REVIEW2__explain_law_far_probe.png`). The code comment reads "the book's
printed +5/3 (slip #1)".

**Should fix** — none new.

**What works**
- The ladder with pulses that speed up, and the table of real cases with the current row highlighted.
- "Switch air ↔ water: this number does not move" answers the step-1 question by experiment.
- The inspector checks the units of the law at the clicked wavenumber and shows the printed power failing them.

## turbulent_energy_budget (E5) — PASS (round 2)

**Must fix** — none (E5-1 fixed).

Also checked: D09's Check now says the 1/Re estimate holds where the mean gradient is ~ ΔU/L, not in the viscous wall
layer, and gives the direct share as 56.4 %, 44.7 %, 37.3 % at Re_τ = 180, 1000, 5200 — the same numbers as the table
of Explain 5 and the whole-channel view (44.7 % at Re_τ = 1000), and as the notebook (cell output "{180: 0.56, 1000:
0.45, 5200: 0.37, 20000: 0.33}"). The bar is labelled "sink (remainder)" ("remainder" in the notebook frame). The two
titles are "Channel · U⁺ · viscous + Reynolds stress · shaded 4×P" and "Budgets at y⁺ = 12 · production 0.238".

**Should fix**
- D10 step 4 takes 9 card pages on a 360×640 phone (the line is shown one row per page). Acceptable for a ★★★ line;
  hiding the step dots on that page would save two pages.
- At 1000×700 the whole-channel title still ends "… + via turbul…"; drop "(Re_τ = 1000, model)" there.
- The note for the undamped preset still ends "simulations put the peak near y⁺ ≈ 12" although the peak is then at
  y⁺ = 3.45; add "with damping".

**What works**
- One orange bar drawn twice with opposite signs, joined by a dashed link: the chapter's central idea in one glance.
- The whole-channel flow diagram (work → direct 44.7 % + via turbulence 55.3 %) and its table against Re_τ.
- Crossing (10.46) and exact maximum (10.39) both given, with the reason they differ by 0.075 wall unit.

## turbulent_jet_similarity (E6) — PASS (round 2)

**Must fix** — none (E6-1 fixed).

Also checked: a scan of the file for decimal literals in reader strings finds no typed result left; the constants 2.5774, 2.1903, 0.1505, 0.2129, 0.1771 occur only as `expect` values of selftest rows. Displayed values against
fluidpy at the page's station: `profile_integrals(0.1)` gives ∫F = 0.21289, ∫F² = 0.15054; U_CL(1.92 m) = 1.860 m/s
(page, at x = 1.9247: 1.858), V̇ = 0.7603 m²/s (0.7612), v_e = 0.0990 m/s (0.0989), Y_CL = 0.1581 (0.158);
`wrong_exponent_fluxes(4, −0.4, 1)` = (1.3195, 2.2974) and (4, −0.7, 1) = (0.5743, 1.5157), as in the wrong-decay
status and check question 2; D13's check reads "0.001 of the stress gradient, and V/U is 0.019" (fluidpy 0.00100,
0.0190). The cards `far`, `vol` and `scal` now say "(plane jet)" in their live lines. Phones have "y [m]", "x: 0 to
10 m" and "U [m/s]" inside the frames (`phone__tour-step1`). The repeated closing sentences of the derivations are gone.

**Should fix**
- Slip #18 (V does not vanish at the jet edge) appears as a "trap" in D13; number it when the others are numbered
  (ruling c).
- D16 step 1 *why* is 36 words.

**What works**
- Raw ↔ rescaled toggle; a wrong exponent that leaves each profile plausible and makes the invariant line turn rose.
- Four flows from one recipe, with the invariant of each in the live line (disc momentum, drag, buoyancy flux).
- Slip #5 worked with numbers for the printed family, and our e^{−ax/2} family passing both tests.

## law_of_the_wall (E7) — PASS

**Must fix** — none.

**Should fix**
- Walkthrough step 2 writes $u_*=\sqrt{\tau_0/\rho}$ beside (12.81); the book prints $u_*^2\equiv\tau_0/\rho$ (the
  equation card and Explain have the printed form). Use the printed form in the step too.
- Explore is 10 pages on a 360×640 phone (9 controls, 8 readouts). Cut the readouts shown there to five (U⁺, log law,
  viscous share, log layer [decades], fitted κ) to bring it to about 7.
- "Equal shares at y⁺ = 9.9" (Spalding's curve) against 10.46 in E5 (mixing length): add "(with Spalding's curve; the
  mixing-length channel of the energy-budget explainer gives 10.5)".
- Rough-wall view title prints "B -17.55" with a hyphen; use "−".
- The model curve fitted from y⁺ = 30 at Re_τ ≈ 5200 returns 0.382 / 3.94 here and 0.383 / 3.99 in the notebook's
  table (ruling j): align the wake strength or the fit grid, or say which is used.

**What works**
- Five DNS profiles collapsing in inner units and fanning out in outer units, by one toggle.
- The draggable fit window: the model built with 0.41 returns 0.382 from y⁺ = 30; the simulation moves the other way.
  "A quoted κ always belongs to a window."
- The inset of the same profile on linear axes ("sublayer + buffer: 3.0 % of δ").

## mixing_length_closure (E8) — PASS

**Must fix** — none.

**Should fix**
- Explain section 5 gives the DNS intercept as "about 4.29" (mean of U⁺ − ln y⁺/κ at fixed κ = 0.384); E7 and the
  notebook give 4.28 (least squares, same window). Use one estimator, or write "≈ 4.3" here.
- Explore is 10 pages on a 360×640 phone; drop the readouts `Ub` and `Cfl` from Explore.

**What works**
- κ turns the line, A⁺ slides it: two constants, two visibly different jobs, and the amber bracket that is B.
- The model's failure shown first (B = −1.23, six units low) before the fix.
- "Model off" returns the laminar parabola at the same pressure gradient: 21.6 times the centreline speed.

## stratified_surface_layer (E9) — PASS

**Must fix** — none.

**Should fix**
- The Code tab prints fluidpy's strings, in which the meteorological line reads "stable ⇔ Γ < Γa: 6.5 < 9.8 K/km" two
  lines under "dT/dz > Γa: −6.5 > −9.8 K/km": one symbol, two values. The badge, legends and Explain write Γd for
  +9.8 K/km. This is `verdict_met` in fluidpy (ch01 wording), not the explainer: for the orchestrator to decide; if
  fluidpy keeps Γa, add a comment line in the block ("fluidpy writes Γa for |g/C_p| in both conventions").
- (12.106) is shown with $-\overline{uw}\,dU/dz$; the book prints $-\overline{uw}\,\partial U/\partial z$ (L2838 via
  `E106I`, and D24). Use ∂ beside the number.
- The correction $-\overline{uw}=u_*^2$ is in the *watch* line of D25 step 1 only and has no number. It is now slip #17:
  move it to the step's *why* in the house wording (ruling c).
- Two values of Ri are on screen (0.12 in the readout and the (12.109) card, 0.0752 in the (12.108) card). Explain
  section 4 says why; label the readout "Ri = Pr_T·Rf (neutral shear)".
- Step 1 takes 4 card pages on a 360×640 phone (step text + notes). Drop `notes: true` from step 1.

**What works**
- H as the transport: the sweep through zero shows L_M jump through infinity and the wind cross the grey logarithm.
- The two-line badge: verdict · L_M · Rf, then both conventions with numbers — at every size, in the chosen order.
- The thermometer trap run in code (+1.110 against −2.213), and the log-linear curve ending at ✕ with the reason.

## taylor_dispersion (E10) — PASS (round 2)

**Must fix** — none (E10-1 fixed).

Also checked against the rendered pages once more: printed p. 604 has $(X_\alpha)_{rms}=(u_\alpha)_{rms}t$ for
$t\ll\Lambda_t$ (12.121), "τ/t in (11.119) is negligible" (slip #11; the page writes the equation meant,
$\overline{X_\alpha^2}=2\overline{u_\alpha^2}\,t\int_0^t(1-\tau/t)\,r_\alpha\,d\tau$ (12.119)) and
$(X_\alpha)_{rms}=(u_\alpha)_{rms}\sqrt{2\Lambda_tt}$ for $t\gg\Lambda_t$ (12.123); printed p. 607 has (12.126)–(12.128)
as shown and $D_T\cong\overline{u_\alpha^2}\Lambda_t$ "for $t\ll\Lambda_t$" as (12.129) (slip #12, shown with
$t\gg\Lambda_t$); the caption on printed p. 606 gives the width ∝ x^{1/2} near the outlet and ∝ x far away (slip #13,
shown the right way round). Step 1 opens at t = 0.3 Λ_t with rays already drawn and says "↺ replays from the release".
Plume mode: the dashed amber parabola is a darker amber than the pale concentration shading and carries its own label
"∝ √x"; legible at 1366×768 and at 360×640 (`REVIEW2__plume_desktop.png`, `REVIEW2__plume_phone.png`). Phones have
"Z [m] against x [m]" and "⟨X²⟩ [m²] against t [s]" inside the frames.

**Should fix**
- On a 360×640 phone the Explore intro ends with the full stop alone on a line after $u_{rms}$; end the sentence with a
  word.

**What works**
- Particle paths, the rms envelope and both limiting laws in one picture; the same story on log axes below.
- "A constant-D cloud starts too wide" — the failure of eddy diffusivity near a source, as a number.
- Regime boundaries on t/Λ_t (0.29 ballistic, 0.31 transition, 2.9 transition, 3.1 diffusive: probed in round 1).

---

## Rulings on the orchestrator's points (updated where round 2 changed them)

**a. turbulent_jet_similarity.** *(updated)* The page computes C₆ = 2.190 (C₅ = 2.577, ∫HF dξ = 0.1771) and now shows
it in every place: the live line and the symbol row of the (12.71) card, pinned by selftest rows; the row "mass
fraction, C6 from the invariant" = 0.10140 agrees with fluidpy. Slip #5 wording is right in Explain 7, card `gen` and
Code `slip` ("the invariant rules out some families, not all exponentials"; nowhere "only power laws"). D13's check
numbers are computed and equal `ch12.thin_shear_layer_terms(1.0, 0.1, 1.0, 1.0, 1.5e-5, C5="from_invariant",
xi_half=0.10)`: 0.00100 and 0.0190.

**b. reynolds_averaging_window.** *(updated)* Both are physically right, and now clearly explained: the line that adds
the OU noise to drift and leak renders. With a random phase per run the wave is a fluctuation (true mean = the decaying
exponential; σ_u² = σ² + B²/2; scatter 0.164 at N = 8); the storyboard's 0.106 is the stationary case — the design
should record the change. σ²·2(r − 1 + e^{−r})/r² is the exact variance of the window average of an OU signal (number
and both limits checked); no fluidpy counterpart yet — acceptable for this gate, follow-up noted.

**c. Book slips #16, #17, #18.** *(updated)* Both new slips I ruled real in round 1 are now in `ch12.book_slips()`
(18 rows), with a third: #16 printed p. 557, $\rho_0U\bar u$ for $\rho_0U\bar v$; #17 printed p. 598,
"$\overline{uw}=u_*^2$" for $-\overline{uw}=u_*^2$; #18 printed p. 575, V does not vanish at the jet edge.
**Leaving them unnumbered in E3, E9 and E6 is acceptable for this gate — Should fix, not Must fix.** In each explainer
the corrected form is the one shown, the reader is told the book prints otherwise (E3: a ⚠️ note on the card; E9: the
*watch* line of D25 step 1; E6: a trap in D13 and step 9's *why*), and nothing on the page is false. Only the label
"slip #k" is missing, and those numbers did not exist when the files were built. Number them in one edit at the
knowledge phase or the next time those files are opened.

**d. energy_cascade_spectrum.** *(updated)* Honest about what is drawn, and the reason the model leaves the −5/3 line
early is now given and is right. Normalisations stated (C₁ = 0.245 two-sided = ½·(18/55)·C with C ≈ 1.5; 0.491
one-sided). Slip #1 shown corrected everywhere; the arithmetic under it renders in every state probed.

**e. turbulent_energy_budget.** *(updated)* Crossing 10.46 and maximum 10.39 both reported and labelled model values.
The sink is described as the remainder in Explain and now on the canvas too ("sink (remainder)"). `uv_plus` is labelled
as minus the correlation everywhere.

**f. law_of_the_wall.** The fit window is clear and correct: DNS at Re_τ = 5186 gives κ = 0.400, B = 4.92 from y⁺ = 30
and 0.384, 4.28 from 350. The (12.92) card solving $\kappa B=1.6[\exp(0.1663B)-1]$ for the slider's κ (4.172 at 0.384)
is acceptable: computed from a printed equation and not labelled as any flow's constants; it coincides with one row of
the book's table of log-law constants, so never attach a flow name to it (the notebook prints the same line). Explore is
10 pages at 360×640 (16 is the Explain count): acceptable. D17's labels of (12.76) are right.

**g. mixing_length_closure.** (12.81), (12.82), (12.88), (12.94), (12.97), (12.98), (12.100), (12.101) are as printed.
The damped intercept is called approximate wherever it appears. Slip #6 on the (12.97) card. Explore 10 pages:
acceptable.

**h. stratified_surface_layer — the lapse-rate rule is kept.** Kundu's convention computed (Γa ≈ −9.8 K/km), the
meteorological one beside it, in the badge at 360×640, 844×390, 1000×700 and 1366×768 for stable, unstable and neutral
states, and in the H help, the convention control, the readout, Explain 6, the (12.108) card, steps 6 and 7, check
question 2 and the neutral note. On phones the temperature strip shows the verdict word and the two-convention line
stays in the badge above it. Trap: +1.110 with Γa, −2.213 with Γa = 0. Log-linear curve ends at ✕ |L_M|/5; continuation
"a commonly used form (often called Businger–Dyer)", no source cited. H as the transport is found by the reader.
Open for the orchestrator: fluidpy's `verdict_met` writes Γa where the page writes Γd.

**i. taylor_dispersion.** *(updated)* Boundaries on t/Λ_t. Slips #11, #12, #13 in corrected form and in the house
wording. The first frame now shows a cloud. The remaining `eq_refs.py` hit is a `ref:` label on a card that shows the TeX.

**j. Consistency.** *(updated)* Same numbers in explainers and notebook: L_M = 83.0 m, U(10 m) = 4.81 m/s; production
peak 0.245; undamped B = −1.23, damped 5.28; D_T = 0.952, 6.32, 10.0 m²/s; jet integrals and C₅; direct share 0.56 /
0.45 / 0.37. **The log-law fit is in the rebuilt notebook**: cell 476 prints "DNS Re_τ ≈ 5186, window 350 < y+ <
0.15 δ+: κ = 0.3837, B = 4.28 (35 points)" and a table with 0.384 / 4.28 from 350 and 0.400 / 4.92 from 30 — the same
as E7. I missed nothing in round 1: that cell is new. **The "B = None" line** (cell 470) is the print-out of the
`LOG_LAW_CONSTANTS` entry `channel_dns_lee_moser_2015`: the cited source quotes κ = 0.384 and no additive constant, so
the dictionary holds `None`. It is correct but reads like a bug to a novice: reword the print to "B: not quoted by the
source (fitted from their profile below: 4.28)" — a notebook Should fix for the lesson review, not an explainer matter.
Remaining differences: DNS intercept "about 4.29" in E8 against 4.28 in E7 and the notebook; the *model* curve fitted
from y⁺ = 30 gives 0.382 / 3.94 in E7 and 0.383 / 3.99 in the notebook's table; equal stresses at 9.9 (E7, Spalding)
and 10.46 (E5, mixing length), each labelled a model; κ = 0.4 in E9, stated. Slip wording: every numbered slip shown now uses the
house form; E3, E9 and E6 carry #16, #17, #18 unnumbered (ruling c).

**k. `tools/embed_check.py ch12`:** `embed check OK for ch12`, exit 0, after the notebook rebuild.

**l. Library and tool findings for the knowledge phase** (collected; nothing was edited):
1. `.viz-der-chain` is `align-items: center`, so a line the pager cuts over two pages is re-centred in its slice
   (E8 and E9 carry `.viz-der-chain.viz-sliced { align-items: flex-start; }`).
2. A `hidePortrait` view keeps its flex share of the row, so its neighbour does not grow (E8 and E9 hide the whole row
   and set `flex-grow` by hand).
3. `<sub>` in titles and status falls below 12 px (E2, E7, E8 add `sub { font-size: max(12px, .8em) }`).
4. No log axes in `v.plot()`: local `logFrame`/`L10` helpers in eight files (E1, E4, E5, E6, E7, E8, E9, E10).
5. `Viz.num.erf` is good to about 1e-7 only (E6 and E10 use it; their rows are at 1e-6).
6. No log-paced transport (E10 drives a 0…1 clock and maps it to time by hand).
7. No two-line status (E9 injects `<br><span class="ssl-conv">` into the status text).
8. **`Viz.tnum(x) + '^2'` is a KaTeX double superscript when x prints as `a\times10^{b}`.** It was live in E1 and E4 in
   round 1; both now carry a local bracketing helper. Move one helper (`Viz.tpow(x, n)`) into the library.
9. **`tools/shot.py` does not fail on `.katex-error` nodes** and audits Explain only in the default and step states; both
   round-1 failures passed it. Add the check and walk the presets with the Explain tab open.
10. A quoted derivation step's "All steps →" button can be cut by a page slice on phones (E2).
11. View titles are truncated with "…" instead of wrapping (E5 at 1000×700; E9 at 1000×700 and 844×390).
12. At least seven files hide the preset strip on phones outside Explore with the same six CSS lines — a library switch
    would replace them.
13. `Viz.fmt`/`Viz.tnum` print values below 1e-12 as 0 unless `keepTiny` is passed; a spectrum far beyond the cut-off
    showed "= 0" in E4 until the builder passed it. Worth a line in `viz_patterns.md`.
14. Process: two builders' probe scripts collided by file name in the shared scratch folder (reported by the
    orchestrator); probes belong in per-slug scratch subfolders, as Part B already says.

## Derivations

Every D id that Part F assigns to an explainer is present: 24 instances, 220 steps, each step with *did*, line, *why*
and *in words*; step counts and move titles equal Part F everywhere; lines differ from Part F only by `\textcolor`
wrappers and `aligned` splits. TeX inside *why* is rendered by the library. I did not re-derive the steps on paper
(Part F was the derivation-reviewer's gate); in round 2 I compared D10 step 4 with Part F term by term.

## The two best screenshots

- `reports/viz/ch12/stratified_surface_layer/desktop__tour-step7.png` — the thermometer trap: wind and temperature
  profiles, the budget at the cursor, the two-convention badge and the code that produces +1.110 and −2.213.
- `reports/viz/ch12/turbulent_energy_budget/desktop__tour-step3.png` — the production peak where the two stresses
  cross, the two budgets joined by the shared bar, and the whole-channel flow of energy.

## Round-1 history (14:06–14:22)

Round 1: **FAIL — 5 of 10** (E1, E4, E5, E6, E10); E2, E3, E7, E8, E9 passed and have not changed since. All ten passed
lint, `shot.py` (8 sizes, 3 092 views) and parity (583/583 rows) in round 1 too; the six Must fixes were found by
reading and by my KaTeX probe:

| # | file | what was wrong |
|---|---|---|
| E1-1 | reynolds_averaging_window L2741 | `…+3.27\times10^{-5}^2+…`: double superscript, boxed Explain line shown as red source text at the default state |
| E1-2 | reynolds_averaging_window L2876–2881 | Code block `prod` used undefined `ens_u`, `ens_v`, `k` |
| E4-1 | energy_cascade_spectrum L2792 | `…\times1\times10^{-4}^{2/3}…`: same failure on the preset "ocean" and at slider ends |
| E5-1 | turbulent_energy_budget L2873 | D10 step 4 was `\overline{u_i\times(\text{step 3})}` instead of its equation |
| E6-1 | turbulent_jet_similarity L3190 | typed "C₆ = 2.240" where the page computes 2.190 |
| E10-1 | taylor_dispersion, five lines | "the page prints" used for the book; slips not in the house wording |

Round-1 evidence screenshots are kept beside the round-2 ones (`REVIEW__explain_rms_error.png`,
`REVIEW__explain_law_ocean.png`).

## Verdict: PASS (all explainers PASS)
