# Chapter 11 — explainer review                                  2026-10-01

Reviewer: `viz-reviewer` (fresh context, read-only on `viz/`). Reviewed the nine files in `viz/ch11/` as they stood at
10:10–10:25 local time (file times 08:37–10:02). Storyboards: `analysis/ch11_design.md` E1–E9.

**Current verdict (round 2, below): PASS — all nine explainers PASS.**

## Round 1 verdict: FAIL — 4 of 9 need a fix (E2, E3, E4, E9); 5 PASS (E1, E5, E6, E7, E8)

All six Must-fix items are small text/CSS edits; no physics function is wrong and no layout fails.

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| normal_mode_growth (E1) | ok | PASS 8 / 224 | 36/36 (29 py) | 8 sections, regime reading right | D01 5/5 · D12 8 (notebook 9) | good | title legend colours lost; thin Bénard strip on 360×640 | **PASS** |
| kelvin_helmholtz_boundary (E2) | ok | PASS 8 / 424 | 40/40 (33 py) | 9 sections, very good | D02 12/12 · D03 14/14 · D04 9/9; **D03 result text broken** | very good | 14 of 17 code lines wrap | **FIX** |
| benard_neutral_curve (E3) | ok | PASS 8 / 392 | 38/38 (32 py; one row at rtol 0.25) | 8 sections; **σ accuracy claim wrong** | D09 6/6 · D10 14/14 · D11 11/11 | very good | clean | **FIX** |
| salt_fingers (E4) | ok | PASS 8 / 208 | 45/45 (40 py) | 8 sections; **"both decay" wrong in 2 regimes** | D13 12/12 | very good | **title legend colours lost** | **FIX** |
| taylor_couette_onset (E5) | ok | PASS 8 / 304 | 42/42 (34 py) | 10 sections, very good | D14 15/15 · D15 6/6 | good | 18 of 20 code lines wrap; Ta slider vs Ta readout past Rayleigh's line | **PASS** |
| richardson_shear_instability (E6) | ok | PASS 8 / 312 | 52/52 (45 py) | 9 sections, excellent | D17 9/9 · D18 14/14 | excellent | clean (legend colours fixed locally) | **PASS** |
| inviscid_shear_criteria (E7) | ok | PASS 8 / 312 | 51/51 (43 py) | 8 sections, very good | D22 9/9 · D19 14/14 | very good | clean | **PASS** |
| orr_sommerfeld_neutral_curve (E8) | ok | PASS 8 / 304 | 27/27 (21 py) | 7 sections, good | D21 10/10 · D23 12/12 | good | title legend colours lost (mild) | **PASS** |
| lorenz_attractor (E9) | ok | PASS 8 / 336 | 43/43 (35 py) | 9 sections, very good | D24 15/15 · D25 11/11 | very good | **title legend colours lost**; 16 of 23 code lines wrap | **FIX** |

Totals: 2 816 views, 374 selftest rows (312 JS ↔ fluidpy), 19 derivations (18 match the notebook's step count).

## Mechanical results (my own runs)

`tools/viz_lint.py --chapter ch11`:
```
ok   viz/ch11/benard_neutral_curve.html
ok   viz/ch11/inviscid_shear_criteria.html
ok   viz/ch11/kelvin_helmholtz_boundary.html
ok   viz/ch11/lorenz_attractor.html
ok   viz/ch11/normal_mode_growth.html
ok   viz/ch11/orr_sommerfeld_neutral_curve.html
ok   viz/ch11/richardson_shear_instability.html
ok   viz/ch11/salt_fingers.html
ok   viz/ch11/taylor_couette_onset.html
```
`tools/eq_refs.py --chapter ch11`: `0 string(s) cite an equation number without showing the equation` (0 in each file).

`tools/shot.py --chapter ch11`, run 2 (complete):
```
PASS  viz/ch11/benard_neutral_curve.html  (8 sizes, 392 views, 7 steps, 38 selftest rows)
PASS  viz/ch11/inviscid_shear_criteria.html  (8 sizes, 312 views, 7 steps, 51 selftest rows)
PASS  viz/ch11/kelvin_helmholtz_boundary.html  (8 sizes, 424 views, 7 steps, 40 selftest rows)
PASS  viz/ch11/lorenz_attractor.html  (8 sizes, 336 views, 7 steps, 43 selftest rows)
PASS  viz/ch11/normal_mode_growth.html  (8 sizes, 224 views, 6 steps, 36 selftest rows)
PASS  viz/ch11/orr_sommerfeld_neutral_curve.html  (8 sizes, 304 views, 7 steps, 27 selftest rows)
PASS  viz/ch11/richardson_shear_instability.html  (8 sizes, 312 views, 7 steps, 52 selftest rows)
PASS  viz/ch11/salt_fingers.html  (8 sizes, 208 views, 7 steps, 45 selftest rows)
PASS  viz/ch11/taylor_couette_onset.html  (8 sizes, 304 views, 8 steps, 42 selftest rows)
```
**Run 1 crashed and is not a result.** It stopped in the parity stage of the first explainer with
`AttributeError: module 'fluidpy.core.stability' has no attribute 'rayleigh_eigs_contour'` while
`fluidpy/core/stability.py` was being edited by another agent (file time 10:08:33; the import worked seconds later).
`fluidpy/ch11_instability.py` changed again at 10:24:42, after E1–E7 of run 2 had been audited. The row counts equal
the builders' reports and every row passed, so nothing moved between their runs and mine; `reference/ch11/rayleigh_spectra.json`
was not rewritten during the review (file time 08:22). If the implementer changes these modules again, re-run `shot.py`.

Extra passes I ran (scripts in the session scratchpad, not in the repo):
- every config string scanned for control characters and lost backslashes → 1 hit (E2, below);
- state changed with the Explain tab open vs a fresh page at the same state, all 50 presets → 0 stale panels;
- click and drag over every view at 1366×768 and 390×844 → 0 page errors; drag-to-orbit works on the E9 scene at both sizes;
- E4 status badge vs `ch11.salt_finger_regime` on 1 120 reachable states → 0 regime mismatches;
- derivation step counts and titles vs `notebooks/build_ch11.py --dump` → see "Derivations" below;
- book pages rendered and read for (11.13)–(11.18), (11.39)–(11.42), (11.45)–(11.46), (11.51)–(11.54), (11.64)–(11.67),
  (11.71)–(11.72), (11.79)–(11.80), (11.83)–(11.86), (11.87), (11.88), (11.90), (11.91), (11.96): every equation shown in
  the explainers matches the page.

Loosest parity rows (all pass): E3 one-mode σ estimate rtol 0.25 (measured 19.6 % off — see E3 M1); E6 weak mode near
the tongue top 0.02 (measured 0.33 %); E4 parcel-run growth 0.02 (0.14 %); E7 fastest-growth scans 0.01 (≤ 0.07 %);
E5 off-grid growth 0.005 (0.22 %); E8 interpolated production 0.005 (0.16 %). Only the E3 row is not an honest tolerance.

---

## normal_mode_growth (E1) — PASS

**Must fix** — none.

**Should fix**
1. View titles: `<b class="c-accent">purple axis</b>`, `σ_r`, `σ_i` render in plain text colour (base rule
   `.viz-view-title b` outranks `.c-*`). Add the five local CSS lines E6 uses, or wait for the library fix (L1).
2. D12 has 8 steps; the notebook has 9 (last: "Solve the quadratic for σ"). The result page shows σ± without the move
   that produces it. Add step 9 with the notebook's title; walkthrough link `derive: {id: 'D12', step: 8}` keeps its number.
3. Explain §0: "the grey dashed curve is the same flow no surface tension" — write "the same flow without surface
   tension" / "the same layer exactly at Ra_c" (the `ghostName` string is used as a sentence fragment).
4. `phone__explore.png` (360×640, Bénard): the layer is a ~45 px strip with seven tiny cells. Draw fewer wavelengths
   when the view is short.
5. Code tab: 5 of 14 lines wrap at 1366×768; keep lines ≤ 44 characters as E3 and E8 do.

**What works**
- One state, three systems (interface, shear, Bénard) with per-system controls hidden by a data attribute.
- σ(k) with the growing band shaded, ▲ cut-off, ◆ fastest mode and a "no surface tension" ghost; the σ-plane with both roots.
- Explain §4's stable/neutral/unstable table and the "stationary vs oscillatory" line computed from the same σ.
- Numbers checked: RT σ_r = 37.0 s⁻¹ at λ = 3.01 cm, cut-off 1.73 cm, Bénard band 0.754–5.36 at Ra = 2000, Ra(6) = 2681.

## kelvin_helmholtz_boundary (E2) — FIX

**Must fix**
1. `viz/ch11/kelvin_helmholtz_boundary.html:2902` — D03 `result.plain` is a normal quoted string with single
   backslashes: `'Eq. (11.18), $c=\bar U\pm[\,\cdot\,]^{1/2}$ with the density-weighted mean $\bar U$, …'`. JavaScript
   turns `\b` into a backspace and drops the others, so the D03 result page shows two red KaTeX errors
   (`c=⌫ar Upm[,cdot,]^{1/2}` and `⌫ar U`; DOM check: 2 `.katex-error` nodes on `#tab=derive&d=2`, result page).
   Make it a raw template (`R\`…\``) like the `tex` line above it, or double the backslashes. `shot.py` did not catch
   this (see L3).

**Should fix**
1. Code tab: 14 of 17 lines wrap at 1366×768 (`desktop__code.png`); reflow to ≤ 44 characters.
2. D03 step 9 `live` reads "(switch surface tension off and the depth to deep to test this line)" when tension is on.
   The derivation's own `set` turns tension off, so it only shows after the reader changes it; print the line with the
   tension term added instead of an instruction.
3. Empty `.viz-stage-badge` (16×4 px) is still in the DOM; not noticeable on the white card (L2).

**What works**
- Three term bars under the root (gravity, surface tension, shear) adding to the discriminant, with the sign deciding.
- The c-plane: two roots collide on the real axis exactly on the boundary (✕) and leave as a conjugate pair.
- Stability boundary with its minimum ◆ 6.70 m/s at λ = 1.73 cm; thermocline preset (k_c = 3.93 m⁻¹, λ < 1.6 m).
- Adaptive status (phone: verdict + one fact). Moving-frame toggle for the vortex sheet.
- D02–D04 checked line by line on paper: signs and factors right; (11.13)–(11.18) match the page.

## benard_neutral_curve (E3) — FIX

**Must fix**
1. Explain §7 says the rigid-wall growth rate is "about 10–20 % off away from the curve". Measured against
   `ch11.benard_growth_rate` (Pr = 7.14): +19.6 % at K = 3.1163, Ra = 2500 (9.905 vs 8.281); +24.7 % at Ra = 1000;
   **+35.1 % at K = 2, Ra = 2500** (2.632 vs 1.949); +14.6 % at Ra = 5000; +4.6 % at K = 5, Ra = 10⁴; rigid–free
   +8.5 %. The stated band is wrong and the parity row carries rtol 0.25. Preferred fix: embed a small σ table from
   `ch11.benard_growth_rate` for rigid–rigid and rigid–free (as E5 does for `taylor_growth_rate`, 4 s.f.), interpolate,
   and tighten the row to ≤ 5e-3. Minimum fix: change the text to the measured range ("up to about 35 % off; the sign is
   exact") and keep the "≈". The sign, hence every verdict and the band edges, is already exact.

**Should fix**
1. D10 `interpret`: "which is $1$ times the bottom of this valley" at the minimum — print "the bottom of the valley".
2. Readout "Im det" shows −0.25 or 137 "on the curve" because K is typed to four digits; show "≈ 0" there as tour
   step 4 already does.

**What works**
- (K, Ra) plane with four neutral curves, ◆ minima, the growing band as a rose bar with its edges labelled (1.76, 5.08).
- The Im det(Ra) view: the determinant's first zero is the marginal Ra — the derivation's step 13 is a picture.
- Explain §1 gives the gradient in three conventions (this chapter, Kundu Ch. 1, meteorology) with the same number.
- Explain §3 prints the 3 × 3 matrix entries with the reader's numbers. Slips #1–#3 quoted in D10/D11 `watch` lines.
- Determinant root vs Chebyshev eigenvalue to 1e-7 at K = 2 and 5; live curves vs the 4 s.f. table to 6e-4.

## salt_fingers (E4) — FIX

**Must fix**
1. `viz/ch11/salt_fingers.html:2781` — Explain §5 always labels the two non-leading roots "both decay". False whenever
   the leading root is one of a complex pair or a second real root is positive. Arctic preset:
   "leading root 9.032 + 72 i … the other two 9.032 − 72 i, −136.7 κ/d² — both decay" (the conjugate grows equally).
   At (dT/dz, dS/dz, d) = (−0.02, −0.001, 0.3): "the other two 3.278, −2132 — both decay". Derive the note from the
   signs: "its conjugate grows equally; the third decays" / "one more grows, one decays" / "both decay".
2. View-title legends lost their colours: `heat`, `salt`, `depth`, `density line`, `finger line`, `purple axis` are
   `<b class="c-…">` inside `.viz-view-title` and print in plain text colour (`desktop__tour-step4.png`). The history
   graph's three curves have no other label on the stage, and tour step 4 speaks of "the blue line" and "the orange
   line". Add the local CSS E6 uses (or L1).

**Should fix**
1. Tour step 3 and Explain §3 print (11.46) as an inequality ("… > 657 (11.46)"); the book prints the marginal equality
   "= 657". Write it as the Equations tab already does: the marginal state "= 657 (11.46); unstable when the left side
   is larger".
2. Status mirror of `salt_finger_regime` text: stripping " of (11.46)" for the badge is fine (the badge states
   "Rs − Ra = … < 27π⁴/4 = 657.5", and the full equation is in tour step 3, Explain §3 and Equations; five checksum
   rows pin the fluidpy text). Do not split the string to hide it from `eq_refs.py`; ask for an allow-marker (L6).
3. The "monotonic overturning" branch is reachable (top-heavy with dS/dz < 0 and d ≳ 0.12 m, e.g. −0.02, −0.001,
   0.3 m) and its wording is right. A preset for it is optional; the dT/dz range need not grow.

**What works**
- The parcel with two gauges (heat, salt) and a history graph beside a regime map with the density line, the finger
  line and the κ_s = κ ghost: "stable by density, unstable by diffusion" is a wedge the reader drags a dot into.
- "Top-heavy is not enough" step (Rs − Ra = 438 < 657.5, flips at −0.0076 K/m) — the review's M1 taught as a step.
- Two sign conventions on the slider (dT/dz and Γ = −dT/dz) and the ⚠ line on the two signs of Ra.
- Regime word equals fluidpy on 1 120 states; D13 checked on paper; (11.45), (11.46) match the page.

## taylor_couette_onset (E5) — PASS

**Must fix** — none.

**Should fix**
1. Past Rayleigh's line the slider reads "Taylor number Ta 1708" while the readout "Taylor number Ta" reads −175
   (`REVIEW__mu1_slider_vs_readout.png`; also tour step 6: −288.9 with the dot at log Ta 3.79). Explain §1 explains it
   ("the Ta slider cannot be met … speed that would give your slider value with the outer one at rest"), but the stage
   does not. Label the hollow dot ("slider value — not reachable here") and put the same clause in the slider's help line.
2. Code tab: 18 of 20 lines wrap at 1366×768 (`desktop__code.png`), comments detach from their lines.
3. Presets, μ range (−0.5 … 1.2), no Galerkin overlay, no Python parity row for the Couette profile: all acceptable.
   The two wall-speed invariants and `Ta = −4AΩ₁d⁴/ν²` cover the base flow.

**What works**
- Ring swap with draggable rings, E before/after bars and ΔE with its sign; Γ² profile with the "falls outward" strip.
- (μ, Ta) plane with the exact curve, the book's fit, Rayleigh's line and the Bénard end point ◆.
- "Rayleigh-unstable but viscously stable" status; the note that (11.54) is +6.5 % at μ = −0.5 because σ = 0 is assumed.
- D14 (15 steps) checked on paper including the 1/k² pressure step and the rescaling that leaves Ta; D15's Bénard map.
- Explain interprets μ = 1 correctly ("do not read μ = 1 as unstable at 1708").

## richardson_shear_instability (E6) — PASS

**Must fix** — none.

**Should fix**
1. Code tab: 6 of 15 lines wrap at 1366×768.
2. Fixed 0–40 L/U₀ clock: acceptable (the title prints "amplitude ×43.7 (drawn capped)" and Explain §7 gives the time
   the cap is reached).

**What works** (best of the chapter; candidates for `knowledge/viz_patterns.md`)
- Ri(z) on a log axis with the ¼ line and a rose band over the levels where N² − ¼U′² < 0 — the band is the set where
  D18's left integral can turn negative (step 12 `live` and `watch` point at it).
- Growth map with the exact neutral curve J = k(1 − k) and a slice at the reader's J with the J = 0 ghost.
- Three-way status: "🛡️ stable for every k" / "🌀 allowed and this wave grows" / "⚖️ allowed, not happening".
- "How far to trust it" paragraph: what a computed 0 means near the neutral curve.
- R = 3, J = 0.4 check: Ri at the centre 0.4 > ¼ and the layer still grows — "everywhere" taught by a counter-example.
- Thermocline and billow-cloud presets in metres and seconds (71.6 s; 2.09 min, 1428 m — arithmetic checked).
- Review M2 covered by nine rows (±inf, nan, shear-free levels). (11.64)–(11.67) match the page.

## inviscid_shear_criteria (E7) — PASS

**Must fix** — none.

**Should fix**
1. Tour step 7 is about the cat's eye, which is hidden on portrait phones; the readout carries the number, the picture
   is gone. Add the one-line hint used elsewhere ("turn the phone sideways to see the eye").
2. Explain numbers its sections 1–8; the other eight start at 0 with "What the views show".

**Howard's semicircle numbering — confirmed on the page image (p. 507).** (11.72) is
$\int[U^2-c_r^2-c_i^2]Q\,dz>0$; the semicircle inequality is unnumbered. E7 labels both correctly ("p. 507
(unnumbered)", "after (11.72)") in Equations, Explain §5, D19 result and tour step 4. E6 does not cite the semicircle.
No other explainer numbers it. D19 step 12's remark that the book first prints the bound without Q is also right.

**What works**
- Profile with U″ and the product (U − U_I)U″ shaded by sign; "Rayleigh ✓ Fjørtoft ✗" preset (sinh profile).
- Semicircle with the eigenvalue's track over k and ◇ fastest growth; Howard's bound beside the actual rate (43 %).
- Fjørtoft's two integrals evaluated on the computed mode as bars (−1 and 0), with the two lobes of (11.84) cancelling.
- sin y with 2b < π: "both conditions met, yet nothing grows"; the "≈" status just past the margin is honest.
- Supplementary Rayleigh tables: 12 rows against live `rayleigh_eigs` solves, all within 6e-4 (table) or 3e-3 (interpolated).

## orr_sommerfeld_neutral_curve (E8) — PASS

**Must fix** — none.

**Should fix**
1. View titles: `rose`, `P`, `Λ` in plain text colour (words still carry the meaning). Same CSS as E6, or L1.
2. Step 1 is a statement; end it with the question ("How can friction make a wave grow?").
3. Between stored modes the streamlines are a phase-aligned blend (source comment, line 2467); nothing on screen says
   so. Add "shape interpolated" to the wave title when the point is not a stored mode, as the eigenvalue already says
   "(table, interpolated)".
4. Imaginary unit is italic `i` in the equations and upright in the numbers; notation convention 4 asks for `\mathrm i`.

**Honesty of the table — acceptable.** Anchors: Re_c = 5772.2218 at k_c = 1.020545 (row vs `poiseuille_critical`, 1e-9);
band at Re = 10⁴ 0.7974–1.0944 against 0.7972–1.0947 (3e-4); Blasius 519.077 on δ* and Bickley 4.0170 (1e-9); Orszag
c = 0.237526 + 0.00374 i; P/Λ = 1.616; budget closes to 2.6e-4 of Λ. The readout says "(table)", Explain says
"computed directly" or "interpolated", and the notes say when another mode decays more slowly than the one followed.
Explain §5 lists our Re_c beside cited literature values, not the book's column.

**What works**
- Thumb with growth contours, ◆ Re_c and a Squire ghost point with its arrow.
- Production/dissipation bars adding to dE/dt = 2kc_iE, and where across the layer the production sits.
- "Onset of small waves, not the Reynolds number of transition" in tour step 7, Explain §5 and the Couette status.
- D21 and D23 checked on paper; (11.79), (11.80), (11.88), (11.96) match the pages.

## lorenz_attractor (E9) — FIX

**Must fix**
1. View-title legends lost their colours: `run 1`, `run 2`, `● C±`, `warm`, `cold` (8 `<b class="c-…">` items) print
   in plain text colour (`desktop__tour-step5.png`: "X(t): run 1, run 2"). The orange and blue traces have no other
   label on the stage. Add the local CSS E6 uses (or L1).

**Should fix**
1. Leftover from the string collision: preset `r = 24` has the tooltip "just below the Hopf point: a very slow spiral"
   (line 3065) while its status and notes say the run "still switches lobes at t = 50 … a chaotic set coexists". Change
   the tooltip to match.
2. The fitted λ depends on δ₀ because the window is 100 δ₀ < |δ| < 0.1: 0.89 at 10⁻⁸, 0.90 at 10⁻¹², **0.59 at 10⁻⁴**
   (preset "δ₀ = 10⁻⁴", `REVIEW__delta0_1e-4_slope.png`). A reader sees the "growth rate" change when only the start
   gap changed. Fit over δ₀ … 0.1 after the flat start, or print "window too short" when it spans under two decades.
3. Explain §5: "|δ| reaches 1 at t = 27.7; compare ln(1/δ₀)/λ = 20.6" leaves a 7-unit gap unexplained. Say that the
   distance does not grow for the first ≈ 12 time units (the flat part of the teal curve).
4. Code tab: 16 of 23 lines wrap at 1366×768.

**Judgements asked for.** The canvas projector with drag-to-orbit is fine (works at 1366×768 and 390×844, no CDN
dependency). 27.7 (fixed-step RK4) vs fluidpy's 29.1 is presented as the lesson in Explain §5 and the Code note —
keep. "≈ 0.89" is labelled a finite-time fit beside the long-time ≈ 0.9 — keep.

**What works**
- Two runs on one clock: orbit, X(t) of both, log₁₀|δ| with the fitted line and the |δ| = 1 mark, and the roll itself.
- Explain §8 "Weather and climate": the time-mean of Z of the two runs (24.07, 23.99) beside their unrelated states.
- r = 24 status: stable C± coexisting with a chaotic set, cited, not computed.
- D24 (15 steps) checked on paper including the Jacobian and the three scale factors; D25's Hopf condition a₂a₁ = a₀.
- (11.90), (11.91) and the steady states match the pages.

---

## Derivations against the notebook (`build_ch11.py --dump`)

19 derivations, 206 steps. Step counts equal the notebook's in 18; order equal in all. **D12: explainer 8, notebook 9**
(the notebook's last step is the quadratic formula; E1 Should 2). Title differences are rewordings that replace a bare
number by words ("Use (11.36)" → "Use the temperature equation", "Subtract (11.11) from (11.12)" → "Subtract step 8
from step 9", "Insert into (11.61)" → "Insert into the Taylor–Goldstein equation"); none changes a move. No walkthrough
`derive` link points at a wrong step.

## Library and tool items (for the orchestrator; builders cannot edit these)

- **L1 · `assets/viz_base.css`**: `.viz-view-title b { color: var(--viz-text) }` outranks `.c-teal` etc. Change the
  selector to `.viz-view-title b:not([class*="c-"])`, then `tools/viz_inline.py --all`. Fixes E1, E4, E8, E9 at once
  (E6 carries a local workaround; E2, E3, E5, E7 colour their legends another way).
- **L2 · empty `.viz-stage-badge`**: add `.viz-stage-badge:empty { display: none }` to the base CSS. Five files carry
  it locally; E2, E5, E7, E9 leave a 16×4 px pill (not noticeable on a white card).
- **L3 · `tools/shot.py`** passed E2 with two `.katex-error` nodes. Fail on any `.katex-error` and on control
  characters in rendered text, on every page of every pager.
- **L4 · deep links**: `#tab=…` works; `#<param>=…` is overridden in all nine because boot applies the hash and then
  `goStep(0)` applies step 1's `set` (8 of 9 probes; the ninth asked for the default value). Apply hash parameters after `goStep`.
- **L5 · code panel width**: about 44 characters at 1366×768. Either lint `code[].src` lines > 44 characters or shrink
  the panel's code font one step. Wrapping: E5 18/20, E9 16/23, E2 14/17, E6 6/15, E1 5/14, E4 3/15; E3, E7, E8 clean.
- **L6 · `tools/eq_refs.py`**: add an allow-marker for strings that mirror a fluidpy text byte for byte (E4).

Helpers worth promoting to `assets/viz_lib.js`: complex-pair formatting (`ctex`, `ctxt`, `cpy`, `csqrtReal`: E1–E4,
E9); `bisect` and `goldenMax` into `Viz.num` (E1); table interpolation on a rectangular grid, bilinear and cubic
(E5–E8); a real-coefficient cubic solver (E4, E9); `arrowPx` (E1, E2); a text-checksum row helper for exact-string
parity (E4, E6); the orthographic 3-D projector with drag-to-orbit as the offline path of `Viz.three` (E9).

## Round 1 verdict: FAIL (E2, E3, E4, E9 need the fixes above; E1, E5, E6, E7, E8 PASS)

---

# Round 2 — re-review after fix round 1                           2026-10-01 (13:30–14:00)

## Verdict: PASS (all nine explainers PASS)

Every round-1 Must fix is closed, checked in the browser and against live fluidpy calls, not from the builders'
reports. The new content (E3 growth table, E7 near-neutral interpolation, E8 two-band jet) is correct at every point I
sampled, honestly labelled, and fits at phone, landscape, notebook and laptop sizes. Five small Should-fix items remain.

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| normal_mode_growth (E1) | ok | PASS 8 / 232 | 36/36 (29 py) | ok | D01 5/5 · D12 9/9 | good | legend colours back; code clean | **PASS** |
| kelvin_helmholtz_boundary (E2) | ok | PASS 8 / 424 | 40/40 (33 py) | ok | D02 12 · D03 14 · D04 9; 0 KaTeX errors | very good | code clean | **PASS** |
| benard_neutral_curve (E3) | ok | PASS 8 / 392 | 54/54 (46 py) | σ from table, exact to 8e-5 | D09 6 · D10 14 · D11 11 | very good | 1 code line wraps | **PASS** |
| salt_fingers (E4) | ok | PASS 8 / 208 | 48/48 (40 py) | root notes right in every regime | D13 12 | very good | legend colours back; 3 code lines wrap | **PASS** |
| taylor_couette_onset (E5) | ok | PASS 8 / 304 | 42/42 (34 py) | ok | D14 15 · D15 6 | good | slider, readout, status and plane agree | **PASS** |
| richardson_shear_instability (E6) | ok | PASS 8 / 312 | 52/52 (45 py) | excellent | D17 9 · D18 14 | excellent | clean | **PASS** |
| inviscid_shear_criteria (E7) | ok | PASS 8 / 312 | 70/70 (62 py) | numbered from 0 | D22 9 · D19 14 | very good | clean | **PASS** |
| orr_sommerfeld_neutral_curve (E8) | ok | PASS 8 / 304 | 32/32 (26 py) | good | D21 10 · D23 12 | good | y-title clipped by one character at 360×640 | **PASS** |
| lorenz_attractor (E9) | ok | PASS 8 / 336 | 43/43 (35 py) | very good | D24 15 · D25 11 | very good | legend colours back; code clean | **PASS** |

Totals: 2 824 views, 417 selftest rows (350 JS ↔ fluidpy), 19 derivations, 207 steps.

## Mechanical results (my own runs, 13:29–13:55; explainer file times 10:41–13:25)

`tools/viz_lint.py --chapter ch11`: `ok` for all nine files.
`tools/eq_refs.py --chapter ch11`: `0 string(s) cite an equation number without showing the equation`.
`tools/shot.py --chapter ch11` (one complete run, exit 0):
```
PASS  viz/ch11/benard_neutral_curve.html  (8 sizes, 392 views, 7 steps, 54 selftest rows)
PASS  viz/ch11/inviscid_shear_criteria.html  (8 sizes, 312 views, 7 steps, 70 selftest rows)
PASS  viz/ch11/kelvin_helmholtz_boundary.html  (8 sizes, 424 views, 7 steps, 40 selftest rows)
PASS  viz/ch11/lorenz_attractor.html  (8 sizes, 336 views, 7 steps, 43 selftest rows)
PASS  viz/ch11/normal_mode_growth.html  (8 sizes, 232 views, 6 steps, 36 selftest rows)
PASS  viz/ch11/orr_sommerfeld_neutral_curve.html  (8 sizes, 304 views, 7 steps, 32 selftest rows)
PASS  viz/ch11/richardson_shear_instability.html  (8 sizes, 312 views, 7 steps, 52 selftest rows)
PASS  viz/ch11/salt_fingers.html  (8 sizes, 208 views, 7 steps, 48 selftest rows)
PASS  viz/ch11/taylor_couette_onset.html  (8 sizes, 304 views, 8 steps, 42 selftest rows)
```
Row counts equal the builders' claims. Loosest rows: E6 weak mode 0.02 (measured 0.33 %), E4 parcel run 0.02 (0.14 %),
E7 fastest-growth scans 0.01 (worst 0.73 %, sin y b = 1.6), E5 off-grid growth 0.005 (0.22 %), E8 interpolated
production 0.005 (0.16 %). E3's rtol 0.25 row is gone; its 15 new rows sit at 2e-3 (worst measured 4.4e-4).

## (a) Round-1 Must fixes — all closed

| item | check | result |
|---|---|---|
| E2 D03 result text | `.katex-error` count over every tab, tour step, derivation page and preset × Explain/Equations/Check | 0 in E2 (70 states) and 0 in the other eight (46–66 states each); no control characters in any config string |
| E3 growth-rate accuracy | page σ vs live `ch11.benard_growth_rate` at 60 off-grid points (rigid–rigid and rigid–free; K = 1.37, 2.03, 3.1163, 4.21, 6.3, 8.7; Ra = 613 … 21 000; Pr = 7.14) | worst relative difference 7.7e-5; round 1's worst case (K = 2, Ra = 2500) is now 2.190 vs 2.190; σ = 0 exactly at 0.999, 1 and 1.001 × Ra(K) — the page snaps to marginal within 0.1 % of the curve |
| E4 "both decay" | Explain §5 at all six presets and two extra states | Arctic and top-heavy oscillatory: "the first is the conjugate of the leading root and grows equally (one oscillating pair); the third decays"; (−0.02, −0.001, 0.3 m): "one more grows, one decays"; the other five: "both decay" (true) |
| E4, E9 legend colours | computed colour of every `<b class="c-…">` in view titles at 1366×768 and 390×844 | E4 6/6, E9 8/8, E1 3/3, E8 3/3, E6 7/7 carry their class colour (library L1 applied) |

Also closed from the Should lists: D12 has 9 steps; code tabs reflowed (0 wrapped lines in E1, E2, E5, E6, E7, E8, E9);
E4 prints (11.46) as the book's "= 657" and states the inequality in words; E5 shows "−175 (< 0)" on the slider, "slider
value: not reachable here" on the plane and the actual Ta in the status; E7 numbers Explain from 0; E8 step 1 ends with a
question and the wave title says "shape interpolated"; E9 preset tooltips match the status and the δ₀ = 10⁻⁴ preset prints
"no fitted line: window too short (2.7 decades)" instead of λ ≈ 0.59; empty stage badges are `display: none` in all nine (L2).

## (b) New content

**E3 growth table** — see the table above: correct to 8e-5 at 60 points that are not table nodes; the "≈" and the
"10–20 %" sentence are gone; Explain fits at 360×640, 844×390, 1000×700 and 1280×720 with no overflow.

**E7 near-neutral interpolation** — page c_i vs live `ch11.rayleigh_eigs_contour` at ten off-grid wavenumbers:

| profile | k | page c_i | fluidpy c_i |
|---|---|---|---|
| tanh | 0.97 | 0.01925 | 0.01925 |
| tanh | 0.33 | 0.54364 | 0.54364 |
| jet, sinuous | 1.85 | 0.01784 | 0.01784 |
| jet, sinuous | 1.97 | 0.00338 | 0.00337 |
| jet, varicose | 0.93 | 0.01276 | 0.01276 |
| jet, varicose | 0.65 | 0.06515 | 0.06514 |
| sin y, b = 1.6 | 0.15 | 0.00699 | 0.00699 |
| sin y, b = 1.7 | 0.30 | 0.03077 | 0.03077 |
| sin y, b = 1.7 | 0.36 | 0.00903 | 0.00903 |
| sin y, b = 2.4 | 0.55 | 0.19420 | 0.19418 |

Largest difference 1.6e-5. The modes near the neutral points (jet k = 1.8–1.97, varicose 0.9–0.95, sin b = 1.6) now grow
on the page as they do in the solver. The landscape-phone layout (profile | cat's eye) is legible.

**E8 two-band jet** — status word vs the sign of c_i from live `ch11.os_leading_mode` at 23 points (Re 3.5–500,
k 0.021–1.9): 23 of 23 agree, including below Re_c, both sides of the gap opening (Re = 17 and 18 at k = 0.068), the
gap (Re 19.5, 22.25, 25, 40, 60), the long-wave band (Re = 22.25, k = 0.025: +2.33e-6 vs +2.45e-6; Re = 25, k = 0.021)
and above the upper edge. Gap edges at Re = 22.25 are 0.0309 and 0.0697, equal to `bickley_neutral_curve` and
`reference/ch11/os_neutral_bickley*.csv`. Magnitudes agree within a few per cent, except inside the gap where c_i is
tiny (−1.69e-4 on the page vs −1.02e-4 live at Re = 18, k = 0.068: right sign, 65 % off in a number of size 10⁻⁴).
Labelling is honest: "its lower edge below k = 0.02, not resolved", "unstable below the map too", "(table)", "shape
interpolated".

The "within the table's resolution of a neutral curve" status is honest and clear. It says that the table cannot tell
the sign of c_i there and that the root-found curve puts the point just inside or outside; it does not claim growth or
decay. Keep it.

Fit: the jet-gap state passes the audit at 360×640, 844×390, 1000×700 and 1280×720 (no overflow, no text below 12 px).

## (c) Derivations against the notebook

All 19 derivations have the notebook's step count and order (D12 9 = 9; 207 steps). D19 step 3 reads "Insert into the
Taylor–Goldstein equation" in the explainer; `build_ch11.py --dump` still prints "Insert into (11.61)" for the notebook.
That is a wording difference, not a different move; align the notebook title if the lesson reviewer wants identical titles.

## (d) Rule 3 and fit

`eq_refs` 0; every equation I had compared with a page image in round 1 is unchanged. E6 now cites the buoyancy
frequency as (7.127): confirmed on the ch07 page image (p. 294) — round 1 should have caught "(7.128)". No text below
12 px and no scroll overflow in the shot run or in my 28 extra state × size screenshots (`REVIEW2__*.png`).

## Remaining Should fix (none blocks)

1. E8 `REVIEW2__jetgap__phone.png`: the rotated y title "wavenumber kL (log)" loses its last character at 360×640;
   shorten to "k L (log)".
2. E8 quiz 4 answer ends "Its lower edge below k = 0.02, not resolved." — add "is".
3. E8 status for the jet at high Re reads "the wave grows already at Re = 500"; drop "already" above, say, 10 × Re_c.
4. E4 code tab: 3 of 15 lines still wrap at 1366×768; E3: 1 of 41.
5. E4 keeps the split-string workaround for the mirrored fluidpy text until `eq_refs.py` has an allow-marker (L6).

## Library and tool items still open (orchestrator)

- L3: `tools/shot.py` still does not fail on `.katex-error`; my sweep found none, but the tool would not have.
- L4: `#<param>=…` deep links are still overridden by tour step 1 in all nine (8 of 9 probes; the ninth asked for the default).
- L6: allow-marker in `tools/eq_refs.py`.
- L1 and L2 are applied in the nine ch11 files; other chapters and templates need the re-inline at the merge gate.

## Verdict: PASS (all explainers PASS)
