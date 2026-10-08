# Chapter 13 — explainer review (round 2)                                  2026-10-08

## Verdict: PASS (all ten explainers PASS; no Must-fix item is open)

| Round | Reviewed | lint | shot | 844×345 audit | selftest rows (ok / JS ↔ fluidpy) | Must fix open | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | all ten | 10/10 | 9 PASS, 1 FAIL | 10/10 PASS | 474 / 474 ok, 370 py | 4 (vertical_modes 2, shallow_water_dispersion 1, rossby_waves 1) | FAIL |
| 2 | the three fixed files (the other seven are unchanged since round 1) | 3/3 | 3 PASS (vertical_modes twice, both PASS) | 3/3 PASS | 482 / 482 ok, 374 py | 0 | **PASS** |

**Round 1 (all ten).** Every file was re-linted and re-audited one at a time (`tools/viz_lint.py`, `tools/shot.py`, 8
sizes, every tab, step, derivation page and pager page), audited again by the same machinery at the in-page landscape
frame 844×345, probed for `.katex-error` / NaN / raw TeX at three sizes in the default state, under every preset and
with the hemisphere set to S, and given a click-and-drag pass over every view at two sizes. Screenshots were read at
phone, laptop, notebook-frame, desktop and 844×345 sizes.

**Round 2 (vertical_modes, shallow_water_dispersion, rossby_waves).** The same battery again on the three changed
files: lint ok; shot PASS (216, 264 and 280 views; vertical_modes run twice, PASS both times); 844×345 audit PASS
(567 pager pages); 0 `.katex-error`, 0 NaN / undefined / raw TeX in the default state, under every preset (including
the no-real-roots band) and with hemi S; click pass 0 errors; `tools/eq_refs.py` 0 hits; derivation step counts and
titles still equal to design Part F (D11 13, D12 10, D13 7, D22 10, D23 9); selftest 48 + 60 + 53 rows, all ok.
Screenshots were captured and read for every item the fixes touched (listed in the three sections below).

**Chapter totals after round 2.** lint 10/10 · shot 10/10 PASS · 844×345 audit 10/10 PASS · KaTeX everywhere, no
`.katex-error` · selftest 482 rows, 482 ok, 374 JS ↔ fluidpy · 14 derivations, 137 steps equal to design Part F.

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| geostrophic_balance | ok | PASS 8 / 152 | 39/39 (30 py) | 8 sections + reading | D03 6/6 | very good | good | **PASS** |
| thermal_wind | ok | PASS 8 / 168 | 31/31 (24 py) | 8 sections, both lapse-rate conventions | D04 7/7 | very good | good | **PASS** |
| ekman_spiral | ok | PASS 8 / 304 | 34/34 (25 py) | 8 sections | D06 14/14 · D07 8/8 | very good | good (phone Explore long) | **PASS** |
| ekman_force_balance | ok | PASS 8 / 176 | 49/49 (42 py) | 8 sections | D08 9/9 | very good | good | **PASS** |
| vertical_modes | ok | PASS 8 / 216 (twice) | 48/48 (39 py) | 9 sections | D11 13/13 | very good | good | **PASS** (round 2) |
| shallow_water_dispersion | ok | PASS 8 / 264 | 60/60 (47 py) | 10 sections incl. energy speed (ours) | D12 10/10 · D13 7/7 | very good | good | **PASS** (round 2) |
| kelvin_wave | ok | PASS 8 / 168 | 51/51 (37 py) | 10 sections | D15 8/8 | very good | good (phone step 2 long) | **PASS** |
| geostrophic_adjustment | ok | PASS 8 / 200 | 67/67 (51 py) | 8 sections | D16 11/11 | excellent | good | **PASS** |
| rossby_waves | ok | PASS 8 / 280 | 53/53 (36 py) | 11 sections; 16 pages at 360×640 | D22 10/10 · D23 9/9 | excellent | good | **PASS** (round 2) |
| eady_instability | ok | PASS 8 / 328 | 50/50 (43 py) | 10 sections | D25 13/13 · D26 12/12 | excellent | good (phone Explain 28 pages) | **PASS** |

## What was checked for all ten, and came out clean

- **Equation fidelity against the rendered pages.** Eleven pages were rendered and read for this review (printed pages
  630, 632, 635, 639, 645, 650, 656, 675, 682, 683 and the transport page), covering at least three numbered equations
  per explainer. Every one matches, with its number beside it: $f=2\Omega\sin\theta$ (13.8);
  $-fv=-\frac1{\rho_0}\frac{\partial p}{\partial x}$, $fu=-\frac1{\rho_0}\frac{\partial p}{\partial y}$ (13.11), (13.12);
  $\mathrm{Ro}=\frac{U^2/L}{fU}=\frac U{fL}$ (13.13); $0=-\frac{\partial p}{\partial z}-g\rho$ (13.14);
  $\frac{\partial v}{\partial z}=-\frac g{\rho_0f}\frac{\partial\rho}{\partial x}$, $\frac{\partial u}{\partial z}=\frac g{\rho_0f}\frac{\partial\rho}{\partial y}$ (13.15);
  $E=\frac{\nu}{fL^2}$ (13.18); $-fv=\nu_v\frac{d^2u}{dz^2}$, $fu=\nu_v\frac{d^2v}{dz^2}$ (13.22), (13.23);
  $\frac{d^2V}{dz^2}=\frac{if}{\nu_v}V$ (13.27); $V=Ae^{(1+i)z/\delta}+Be^{-(1+i)z/\delta}$, $\delta=\sqrt{2\nu_v/f}$ (13.28), (13.29);
  $-fv=\nu_v\frac{d^2u}{dz^2}$, $fu=\nu_v\frac{d^2v}{dz^2}+fU$ (13.33), (13.34); $\frac{d^2V}{dz^2}=\frac{if}{\nu_v}(V-U)$ (13.37);
  $\frac{d\psi_n/dz}{N^2\int_{-H}^z\psi_n\,dz}=\frac{\rho_0}g\frac{w_n}{\partial\rho_n/\partial t}\equiv-\frac1{c_n^2}$ (13.55);
  $\frac d{dz}\big(\frac1{N^2}\frac{d\psi_n}{dz}\big)+\frac1{c_n^2}\psi_n=0$ (13.56); the modal set (13.57)–(13.59);
  $w_n=\frac1{c_n^2}\frac{\partial p_n}{\partial t}$ (13.61); $c_n^2\equiv gH_e$ (13.62); (13.72)–(13.74);
  $\frac{\partial^3v}{\partial t^3}-gH\frac\partial{\partial t}\nabla_H^2v+f_0^2\frac{\partial v}{\partial t}-gH\beta\frac{\partial v}{\partial x}=0$ (13.75);
  $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76); the Kelvin set (13.84), (13.85) and $c=\sqrt{gH}$ (13.86);
  $c_x=\frac\omega k=-\frac{\beta}{k^2+l^2+f_0^2/c^2}$ (13.119); $c_x=U-\frac{\beta}{k^2+l^2+f_0^2/c^2}$ (13.120);
  the two lid conditions and $c=\frac{U_0}2\pm\frac{U_0}{\alpha H}\sqrt{\big(\frac{\alpha H}2-\tanh\frac{\alpha H}2\big)\big(\frac{\alpha H}2-\coth\frac{\alpha H}2\big)}$ (13.141).
  The unstable-band inequality is shown as $\frac{HN}f<\frac{\alpha_cH}k$ (13.142) with the cut-off kept symbolic (our
  five-digit value beside it) — deliberate, design header convention 12.
- **Derivations, on paper.** All 137 steps were followed line by line: the square root of $i$ and the polar form of $A$
  (D06), $M=-i\tau/(\rho f)$ (D07), the overshoot at $z=\tfrac34\pi\delta$ (D08), the Sturm–Liouville cross-multiplication
  and the boundary bracket for either lid (D11), the cancellation of $f_0\,gH\,\partial_x(\nabla\cdot\mathbf u)$ between
  steps 4 and 7 (D12), the factor $i$ and the AM–GM bound (D13), the sign of $d\hat\eta/dy$ (D15), the energy integrals
  $\tfrac12\rho g\eta_0^2\Lambda$ and $\tfrac32\rho g\eta_0^2\Lambda$ (D16), the circle and both group-velocity
  components (D23), the $w'$ formula and the cancelling pair (D25), the 2 × 2 determinant and the factorisation through
  $\tanh X+\coth X=2\coth2X$ (D26). No skipped move, no wrong sign or factor.
- **Hemisphere handedness in the drawing (f < 0), against the Python.** Screenshots at 35° S (60° S for the Ekman pair):
  `geostrophic_balance` — quiver clockwise round the low, Coriolis arrow to the left of the wind;
  `thermal_wind` — pole and cold air on the left, still ⊙ westerly, the $f$ and gradient bars swap sides, shear bars stay;
  `ekman_spiral` — surface current 45° to the left, summed arrow 90° to the left (north for an eastward wind);
  `ekman_force_balance` — high pressure drawn on the north side, wind leaning toward low to the right of the geostrophic
  wind, low circled clockwise with inflow, rising in both hemispheres;
  `kelvin_wave` — trapped only toward −x with the coast on the left, $e^{+y/\Lambda}$ shown and refused for +x;
  `geostrophic_adjustment` — jet along −y, drawn ⊙ with the correct right-handed sense for an x–η section.
  All correct, and every status sentence agrees with the picture.
- **Honest labels and house forms.** "ours — not in the book" is present on Ekman pumping, the thermal wind in
  temperature form, the whole of geostrophic adjustment, the Rossby group velocity and the Eady growth rate in physical
  units (the one gap of round 1, the energy speed in `shallow_water_dispersion`, is closed). Slips and traps are numbered as `ch13.book_slips()` (14) and
  `ch13.traps()` (18): slip #4, slip #7, T4, T5, T6, T8, T10, T14, T15, T17 all point at the right entry. "The book
  prints" everywhere; no sentence of tour or quiz text starts with a symbol; no book prose, figure or table.
- **Known points (a)–(g), judged.** (a) `ekman_spiral` — Should fix, below. (b) `vertical_modes` — was a Must fix, closed in round 2; the
  1e-3 row is acceptable (two discretisations; the same quantity is pinned at 1e-9 against the shooting route and the
  matrix mirror at 1e-8). (c) `rossby_waves` — the phone strip is acceptable (crest ▲ and group ◆ stay legible), the
  37-page Explain was a Must fix (closed in round 2: 16 pages), the free-surface / rigid-lid wavelengths are explained clearly (quiz 4, D23 check, code
  note, a parity row for each), columns without advection was a Should fix (closed in round 2). (d) `kelvin_wave` — Should fix. (e)
  `geostrophic_balance` — wording correct: "overshoots … loops for ever about a steady drift", "without drag the
  oscillation never dies"; "settles" appears only in the drag-on quiz answer. (f) **Eady numbers:** the explainer's
  ocean preset (N = 5 × 10⁻³ s⁻¹ over the top 1000 m, U₀ = 0.1 m/s → Λ_E 59.8 km, 22.3 days, 234 km) is exactly the
  "ocean twin" of the executed notebook (cell 762 prints the same three numbers); the atmosphere agrees too (1183 km,
  1.64 days, 4630 km, cut-off 3099 km). `scripts/ch13_instability.py` uses **different inputs** (N = 2.7 × 10⁻³ s⁻¹,
  U₀ = 0.11 m/s over 1000 m → about 11 days and 126 km). Nothing is wrong numerically; the script should be brought in
  line with the notebook (orchestrator: an implementer item, not an explainer item). (g) see "Across the chapter".

---

## vertical_modes — PASS (round 2)
**Must fix** — none open. Round 1: M1 (shot FAIL, equation card of walkthrough step 3 cut through a line at 360×640)
and M2 (stack rows crowded, values through the ◇ markers) are **closed**: shot PASS on two consecutive runs; step 3 at
360×640 is now 4 pages, the derivation quote replaces the equation card (captured and read, pages 1–3); on a short
stack only the selected row carries its value and every ◇ is clear, for uniform N and for the thermocline preset.

**Should fix**
- Walkthrough step 3 at 360×640: page 2 of 4 shows only the head of the derivation quote and its "All steps" button.
- Explain at 360×640 is 19 pages (text body about 100 px).
- The first-mode Rossby radius reads 43.1 km here (free surface, shooting) and 43.2 km in the three explainers that use
  $c_1=NH/\pi$; say so once in Explain §4.

**Closed since round 1** — the loose parity rows now carry their reason in the label.

**What works** — the three linked views (shape, stack, travelling mode) on one mode number; ◇ quick estimate labelled
ours; orthogonality with weight 1 for either lid shown as a live overlap in the title; a very thorough parity set (two
independent eigen-solvers mirrored in JS and each pinned to fluidpy).

## shallow_water_dispersion — PASS (round 2)
**Must fix** — none open. Round 1 M1 (energy speed shown but neither worked out nor labelled) is **closed**: Explain
§7 "Crest speed and energy speed (ours — not in the book)" differentiates
$\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76) to
$\big(3\omega^2-c^2K^2-f_0^2\big)\frac{\partial\omega}{\partial k}-2c^2\omega k-c^2\beta=0$ and solves it. The line was
re-derived on paper (correct, and the $l$ version has no β term), and the live numbers were recomputed by hand for
the slow root of the default state (35° N, external mode, 3100 km): numerator −7.12 × 10⁻⁷, denominator −1.76 × 10⁻⁷,
ratio 4.04 m/s eastward — the page shows the same three numbers, and the plan-view arrow and the readout agree. "ours"
is on the readout, the plan-view legend and title, the notes and the Explore callout. Read at desktop (page 4 of 5),
360×640 (pages 15–17 of 24) and 844×345 (page 12 of 17). The no-real-roots band still shows only cards and dashes
(no number, no NaN) in all three views and in Explain. Six new parity rows: the fast root against the Poincaré group
velocity at 1e-10, the slow root against the quasi-geostrophic formula with the approximation named in the label, and
two centred differences of the roots at 1e-6. Walkthrough step 2 at 360×640: 4 pages (was 6).

**Should fix**
- Explain at 360×640 grew from 18 to 24 pages with the new section (no heading-only page any more).
- Explain §7 on a desktop: the "why" note of the differentiation line wraps beside the second row of the formula; put
  it under the formula.

**What works** — the best "which term balances which" picture so far: log term bars for the selected root with
"smallest / largest", the three roots on one log–log diagram with the Kelvin line dashed as "not a root", β on/off, the
slip #7 toggle with its own quiz answer, the discriminant with its bound; and now the group velocity of any root from
one implicit differentiation, with "crest speed × energy speed ≈ c²" for a fast wave.

## rossby_waves — PASS (round 2)
**Must fix** — none open. Round 1 M1 (Explain 37 pages with a heading-only first page at 360×640) is **closed**:
Explain is 16 pages there, the text body is about 250 px, the first page carries the whole of section 0, no page has
fewer than 120 characters; Explore is 7 pages. On portrait phones these two tabs keep the ω(k) curve only, and Explain
says so in a hint.

**Judged in round 2**
- Columns with U > 0: they now drift east with the current and re-enter at the western end; the row stays evenly
  filled, the selected column keeps its ring, and Explain §8 says what happens at the edge. Two parity rows pin the
  drift and the wrap. Sensible.
- 844×345, Derivation tab: the title-only first page is gone (step 1 of D22 opens with the "we had" line). The curve's
  view title is cut to "group c_gx = +0.435 c…" there. **Acceptable, Should fix**: the full value with its unit stands
  in the status badge directly above it at that size.

**Should fix**
- 844×345 Derivation tab: shorten the curve title so the unit is not cut ("c_gx +0.435 cm/s").
- Explore at 360×640 shows the curve alone and leaves about 170 px of the text panel empty: the height strip with the
  crest ▲ and group ◆ would fit, and Explore would then show the phenomenon on a phone as the walkthrough does.

**Closed since round 1** — columns advected; title-only derivation pages; the reason for the 35° default is given.

**What works** — the chapter's best explainer on a desktop: crest ▲ against group ◆ on one clock, chord = phase speed
and tangent = group velocity on ω(k), the wavenumber circle with the energy arrow pointing at its centre, and the
clicked column whose arithmetic ($\beta Y$, stretching share, resulting spin) is traced in the walkthrough card
(`reports/viz/ch13/rossby_waves/desktop__tour-step2.png`). Stationary preset with both wavelengths explained.

---

## geostrophic_balance — PASS
**Should fix** — walkthrough steps 2–4 need 6–8 pages at 360×640 (step 3 carries an equation card and a derivation
quote of the same line: drop one on portrait phones) · check the D03 result page's live sentence "the leftover arrow is
… % of the Coriolis arrow" at t = 0, where the Coriolis arrow is still nearly zero (show it only after the first hour).
**What works** — cycloid arches on the map with the same four coloured arrows in the force view, tip-to-tail closure,
term bars per axis that sum to zero, the "balanced start" contrast in step 3, Ro written as a ratio of lengths (loop
radius over the size of the low) with the latitudes at which Ro reaches 0.3 and 1.

## thermal_wind — PASS
**Should fix** — step 4 is 7 pages at 360×640 · step 6 says "Switch to S" but leaves the chip on N (fine as an
invitation; say "use the chips below").
**What works** — "two sign changes, same westerly" as bars that swap sides while the shear bars stay; pressure surfaces
fanning out as the physical reason; the two-row lapse-rate table (Kundu's sign computed, the meteorological one beside
it); ocean mode with a level of no motion.

## ekman_spiral — PASS
**Should fix** — (a) Explore is 10 pages at 360×640: mark "Interior flow", "Coast" and "Wind toward" `optional` · code
comments wrap inside the walkthrough card even on a desktop (`desktop__tour-step6.png`): keep comment lines under about
46 characters or move the live value to its own line · the live comment "current at the cursor: (0.0001, 0.0001) m/s"
at 5 δ needs one more significant figure · Explain at 360×640 opens with a heading-only page.
**What works** — the depth cursor shared by the 3-D spiral, the hodograph and the profiles; the purple running sum that
swings to 90°; "four times the viscosity: twice as deep, half as fast, same transport" as a walkthrough step with the
code line that contains no ν; coast preset with upwelling / downwelling.

## ekman_force_balance — PASS
**Should fix** — step 3 is 8 pages at 360×640 · parity row "height of the overshoot / delta" has rtol 1e-5 without a
reason (the value agrees to 1e-8: tighten it or say how it is found).
**What works** — the closed force triangle that leans toward low pressure in both hemispheres; ▲ for the largest u; the
low and the high seen from above with rising / sinking; pumping clearly marked ours.

## kelvin_wave — PASS
**Should fix** — walkthrough step 2 is 9 pages at 360×640 and opens with the step title alone: drop the `eq` card on
portrait phones (the derivation quote shows the same third equation) · (d) 844×345, Derivation step 1: page 1 of 6 is
the title alone (`reports/viz/ch13/kelvin_wave/frame/frame__derive-d1p1.png`): shorten the title · the Equations card
writes $\Lambda\equiv c/f$ while every number uses $c/|f|$: add "(for f > 0; in general c/|f|)".
**What works** — three equations as three pairs of equal bars at the cursor; the wrong direction drawn growing offshore
and refused, with the status telling the reader to turn the wave round; channel mode with both waves; the straight
line ω = ck under the Poincaré curve.

## geostrophic_adjustment — PASS
**Should fix** — at 844×345 two view titles end in "…" and lose their number (`<scratch>/vizreview13/frame/geostrophic_adjustment/frame__derive-d1p3.png`):
use the short titles there · thirteen selftest rows are looser than 1e-6; most labels give the reason, "period mean …
vs numpy" (2e-6) does not · the end-of-run card covers the top of the curve near the front.
**What works** — no rotation first, then the same step on a turning earth; snapshot (orange) against its mean over one
inertial period (teal) against the prediction (dashed); the conserved quantity drawn flat in amber; the energy bill
with one third kept; narrow and wide bumps as presets. The clearest "aha" of the chapter in the first three steps.

## eady_instability — PASS
**Should fix** — Explain is 28 pages at 360×640 (text body about 50 px): same remedy as `rossby_waves` M1 · parity row
"tilt of the fastest wave is a quarter wavelength" has rtol 1e-4 while the value agrees to 5e-9.
**What works** — tilt, growth curve and heat flux on one wavelength cursor (`reports/viz/ch13/eady_instability/desktop__tour-step5.png`);
the cut-off ◆ and fastest ▲ as presets on both sides of the threshold; "growth rate is k c_i, not c_i" as a quiz with
both curves drawn; the wedge of exchange paths; an independent Chebyshev eigenvalue in the parity rows; numbers equal
to the executed notebook.

---

## Across the chapter (Should fix, for the round-2 builders and the knowledge-keeper)
- **Walkthrough length at 360×640.** Steps that carry both an equation card and a derivation quote take 7–9 pages
  (kelvin 2, geostrophic_balance 3, ekman_force_balance 3, ekman_spiral 3, thermal_wind 4, geostrophic_adjustment 5).
  One heavy extra per step on portrait phones.
- **Colour of velocity.** Forces follow convention 10 everywhere (pressure orange, Coriolis teal, friction rose,
  acceleration purple). Velocities do not have one colour: grey in `geostrophic_balance` and `ekman_force_balance`, teal
  in `ekman_spiral` and `kelvin_wave` (where teal is also the Coriolis force), purple in `thermal_wind` and
  `geostrophic_adjustment`. Pick one for round 2 or record the exception in `knowledge/viz_patterns.md`.
- **Defaults against `ch13.illustrative_inputs()`.** Latitude 35° (60° for the two Ekman explainers), ocean 4200 m and
  N = 2.7 × 10⁻³ s⁻¹, atmosphere 9 km / 1.1 × 10⁻² s⁻¹ / 27 m/s, τ and ν_v, U_g, the two-layer case and the 3100 km
  wavelength all agree. Exceptions, both labelled: the Eady ocean twin (notebook inputs) and the Rossby default latitude.
- **Library note (orchestrator).** At 844×345 the first page of a derivation's step 1 can be the title alone, and at
  360×640 a numbered Explain heading can be left alone on a page: a "keep heading with next block" rule in the pager
  would remove both for every chapter.

## Candidates for `knowledge/viz_patterns.md`
Snapshot vs period mean vs prediction for a ringing adjustment (E8) · chord and tangent on ω(k) with a wavenumber
circle whose centre the energy arrow points at (E9) · three equations as three pairs of equal bars at a cursor, with
the refused solution drawn (E7) · sign bars that swap sides while the result bar stays (E2) · log term bars with
"smallest / largest" for the selected root of a dispersion relation (E6) · one cursor shared by a 3-D spiral, a
hodograph and a running sum (E3) · reviewer tooling: `shot.audit_explainer` re-used at a ninth size (844×345) and a
"characters per pager page" scan that finds heading-only pages.

## Verdict: PASS — all ten explainers PASS after round 2; the Should-fix lists above stay open for the knowledge-keeper and any later polish pass
