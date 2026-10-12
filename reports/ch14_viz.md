# Chapter 14 — explainer review (RE-REVIEW after fix round 1)                2026-10-12

## Verdict: FAIL — all five round-1 Must-fix items are FIXED and the nine files need no further edit; the gate is held by one NEW Must-fix in the shared library (L1: the Explain tab is packed wrongly after any change of state made while it is open)

Round 1 (same day, commit `9a1e234`): 4 PASS, 5 FAIL with one Must-fix each. This round reviews the six files changed
since then (`9cd970d`, `844ea51`): E1 `zhukhovsky_airfoil_builder`, E2 `kutta_condition`, E3 `starting_vortex`,
E4 `airfoil_pressure_lift`, E5 `lift_curve_stall`, E9 `thrust_from_lift`. The other three were linted only, plus one
cheap probe of the library defect found here (L1). Everything below was re-run in this session; the builders' reports
were treated as claims.

| Check (six changed files unless stated) | Result |
|---|---|
| `tools/viz_lint.py --chapter ch14` | 9 / 9 ok |
| `tools/shot.py` on the six (8 sizes, every tab, step, derivation page and pager page) | 6 / 6 PASS, 1168 views, 2741 pager pages, 0 failures, exit 0 |
| selftest | 193 rows, 193 ok, 146 JS ↔ fluidpy (chapter: 328 rows, 240 JS ↔ fluidpy) |
| tolerances against round 1 (every row, old file from `git show 9a1e234`) | **0 loosened**; 11 tightened (E1: 5, E4: 6); 1 row added (E3); none removed. Largest error / tolerance after tightening: 0.29 (E1 `x_tmax/c measured`), 0.11 (E4 `separation, α = 6°`) |
| `.katex-error` probe: every tab, every walkthrough step, every derivation page and **every page of every pager**, plus Explain and Explore under every preset, at 1366×768, 360×640 and 844×345 | 624 views, 3098 pager pages: KaTeX loaded everywhere, **0 `.katex-error`**, 0 text fallback, 0 NaN / undefined / Infinity, 0 page errors |
| keyword calls, parsed with `ast` from the live selftest rows and the Code cards with their live values | 172 fluidpy calls in 146 `py:` rows: **all by keyword** (method and closure calls included). 34 calls in 12 Code cards: 32 by keyword, **2 with one leading positional** (E5 card `model`, Should fix) |
| Code tabs of E1 and E2 executed in Python with the live values | 3 / 3 cards run; every value equals its comment (E1: a 1.0826, μ −0.0800+0.0755j, chord 4.0222, t/c 0.0964, camber/c 0.0346, rules 0.1039 / 0.0349, mark 0.955 / 0.177 · E2: Γ_K 2.1283 twice, C_L 0.987 / 1.058, Γ 1.0641, rest point −4.486°, q 1.55 U, edge 0.9123 U, L 293.3 N/m) |
| E9 at the in-page frames, shot.py's own checks, sail + fish + bird | 844×345: 53 views, 254 pager pages, **0 failures**; 667×330: 53 / 369 / 0; 844×390: 53 / 208 / 0 |
| E4 at 360×640, 390×844, 844×390, 667×375, 844×345 — Explore, walkthrough step 6, and α = −6° with the plate ghost | 0 audit failures; screenshots read (one cosmetic clip at 667×375, Should fix) |
| `tools/eq_refs.py --chapter ch14` | 14 hits, the same ref labels and selftest row names as in round 1; nothing new |
| `tools/check_public.py viz/ch14/*.html` | OK (9 files) |
| **NEW — state changed while the Explain tab is open** (each preset, then shot.py's pager audit; all nine files, six sizes) | **7 of 9 explainers show clipped text** at one size or more; re-entering the tab or calling the pager's `layout()` repairs every case → library, **L1** |

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | round-1 Must-fix | verdict |
|---|---|---|---|---|---|---|---|---|---|
| zhukhovsky_airfoil_builder | ok | PASS 8 / 336 | 35/35 (30 py) | 9 sections + reading | D03 15/15 · D04 11/11 | very good | very good | FIXED | **PASS** |
| kutta_condition | ok | PASS 8 / 240 | 38/38 (29 py) | 9 sections + reading | D07 4/4 · D09 10/10 | very good | very good | FIXED | **PASS** (L1 seen) |
| starting_vortex | ok | PASS 8 / 144 | 21/21 (16 py) | 6 sections + reading | D02 5/5 | very good | good | FIXED | **PASS** (L1 seen) |
| airfoil_pressure_lift | ok | PASS 8 / 168 | 37/37 (29 py) | 9 sections + reading | D10 7/7 | very good | good | FIXED | **PASS** (L1 seen) |
| lift_curve_stall | ok | PASS 8 / 96 | 26/26 (19 py) | 8 sections + reading | none (by design) | very good | very good | none | **PASS** (L1 seen) |
| downwash_induced_drag | ok | not re-run (round 1: PASS 8 / 304) | round 1: 36/36 (23 py) | 8 sections + reading | D11 4/4 · D12 7/7 · D13 8/8 | very good | good (phone clutter) | none | **PASS** (L1 seen) |
| lifting_line_solver | ok | not re-run (round 1: PASS 8 / 296) | round 1: 61/61 (47 py) | 10 sections + reading | D15 11/11 · D18 9/9 | very good | very good | none | **PASS** (L1 seen) |
| elliptic_loading_aspect_ratio | ok | not re-run (round 1: PASS 8 / 264) | round 1: 38/38 (24 py) | 10 sections + reading | D19 9/9 · D21 8/8 | very good | very good | none | **PASS** (L1 seen) |
| thrust_from_lift | ok | PASS 8 / 184 | 36/36 (23 py) | 9 sections + reading | D22 8/8 | very good | good | FIXED | **PASS** |

"PASS" in the last column is the verdict on the file's own content. The chapter verdict is FAIL only because of L1.

---

## L1 — NEW Must-fix, shared library (`assets/viz_lib.js`, not any explainer file)

**What the reader sees.** Open the Explain tab and change the state while it is open — click a preset chip (they sit
right above the picture on that tab), drag the transport, or let the animation run. The text is cut off inside the
panel and stays cut until the reader leaves the tab and comes back.

- `reports/viz/ch14/_rereview/L1__starting_vortex__desktop__explain_after_preset_stop.png` — 1366×768, real click on
  the tab, real click on the preset "stop", 3 s wait: "Reading the current setting" is cut 33 px on the last page (3/3),
  so the regime sentence is nowhere on screen.
- `…/L1__lift_curve_stall__phone-tall__explain_pg8.png` — 390×844: the line "lift lost against the line = 9 % …" is
  cut through its second row. At 390×844 this is **new since round 1** (the sentence added to §6 moved the page
  break); at 360×640 it was already there in round 1 (19 px) and I missed it.
- `…/L1__airfoil_pressure_lift__tablet__explain_after_preset.png` — 768×1024: a line of section 5 is 47 px below the
  page (every one of the five presets).
- `…/L1__downwash_induced_drag__desktop__explain_after_preset.png` — 1366×768: all six presets, 2–47 px.

**Scope measured** (Explain tab opened paused; each preset applied; shot.py's `JS_PAGER_AUDIT`; presets with a clipped
page / presets):

| slug | 360×640 | 390×844 | 844×390 | 768×1024 | 1000×700 | 1366×768 |
|---|---|---|---|---|---|---|
| starting_vortex | 5/5 | 4/5 | 5/5 | 3/5 | 3/5 | 4/5 |
| downwash_induced_drag | 6/6 | 6/6 | 5/6 | 6/6 | 6/6 | 6/6 |
| lift_curve_stall | 5/5 | 5/5 | 5/5 | 0 | 0 | 0 |
| airfoil_pressure_lift | 0 | 5/5 | 0 | 5/5 | 0 | 0 |
| lifting_line_solver | 0 | 4/6 | 4/6 | 0 | 0 | 0 |
| elliptic_loading_aspect_ratio | 0 | 0 | 0 | 1/5 | 1/5 | 0 |
| kutta_condition | 0 | 0 | 0 | 1/6 | 0 | 1/6 |
| zhukhovsky_airfoil_builder, thrust_from_lift | 0 | 0 | 0 | 0 | 0 | 0 |

Calling the pager's own `layout()` afterwards repaired every one of these cases (all of them, at every size).

**A second symptom, same function.** While the transport is playing, the Explain pager jumps back to page 1 after
every click on "›" (E5 and E3: six clicks, the label stays "1 / 21" and "1 / 14"; paused, the same clicks give 2, 3, 4,
5). Walkthrough step 1 of E3 and E5 starts the animation, and E5's sweep loops, so a reader who goes from the Guide
to "Why" cannot turn the page until they find the pause button.

**Where.** `assets/viz_lib.js`, `updateCalc()` (lines 1230–1249): when the Explain HTML changes it calls
`calcPager.setItems(items)` — which sets `page = 0` (line 1113) — then `calcPager.layout()` (line 1241), and only
afterwards writes the `[data-live]` values (line 1246). The pages are therefore measured with the "…" placeholders
in place. I did not prove that this ordering explains the largest clips (47 px); what is certain is that one more
`layout()` after the update repairs all of them.

**Concrete change** (orchestrator or knowledge-keeper; then `tools/viz_inline.py --all`):
1. In `updateCalc`, fill the live values first and lay the pager out after that (when `app.tab === 'explain'`), and
   lay it out again when a live value changes an item's height.
2. Keep the reader's page across `setItems` when the number of items is unchanged (clamp to the new page count), so
   playing does not throw the reader back to page 1.
3. `tools/shot.py`: on the Explain tab, apply each preset (and one transport step) and repeat the pager audit — the
   present audit never changes the state on that tab, which is why 2032 + 1168 views passed.
4. Re-run `tools/shot.py --chapter ch14`. The library is inlined in every chapter's explainers; this was not tested
   outside chapter 14.

**Why it holds the gate.** Clicking the presets to watch the interpretation change is the intended use of that tab,
and "Reading the current setting" is the paragraph that gets cut. It predates this chapter and no builder caused it;
if the orchestrator prefers to track it as a library ticket and pass chapter 14 on the files alone, the per-file
verdicts above are the ones to use.

---

## zhukhovsky_airfoil_builder — PASS
**Round-1 Must fix**
1. Positional fluidpy calls — **FIXED.** Code card `shape`: 5 calls, all by keyword; 30 `py:` rows, 32 calls, all by
   keyword (`b=, e=, beta=, n=, mu=, zeta=, zeta_re=, zeta_im=, exact=`). 35/35 rows ok; the card runs and its
   numbers equal the comments.

**New Must fix** — none.

**Round-1 Should fix — all five done**
- D03 "What it means right now": `\tfrac12\beta=6.98\,\%` renders ("7.03 % of the chord against the rule's ½β = 6.98 %").
- D04 step 5 *your numbers* now compares with the β = 0 circle and says so ("mapped exactly (β = 0): 0.1490" against
  0.1432 at e = 0.1, θ = 75°, with β = 8° set).
- `phone-tall__derive-d1p7.png`: the triangle O–Q–P is clear; "μ · drag me" and the long edge labels are gone there.
- `notebook__equations.png`: both titles fit ("Circle · μ (−0.100, 0.000) b", "Section · t/c 11.8 % (rule 13.0 %)").
- Tolerances tightened: four "measured" rows 2e-3 → 1e-5, `x_tmax/c` 5e-3 → 2e-3.

**Should fix (new, small)**
- At 1366×768 in the Derivation tab the two view titles are still cut ("centre μ = −0.100 + 0.155 i,…" loses the
  radius; "t/c 12.0 % (rule 13.…"). `medium()` tests the window (< 1100 px); test the view's width instead.

**What works** (keep; candidates for `knowledge/viz_patterns.md`)
- A draggable centre with two guides ("up = β", "back = e") that *is* the idea; presets are the book's four cases plus
  "too fat".
- Rule-against-measured bars with the gap in words, and a small curve of the rule leaving the truth as e grows.
- The constant 3√3/4 pinned by a parity row ("computed, not typed").
- The D03 check with a number the reader can do by hand (tan β = ¾ → 1.5 = −0.583 + 2.083).

## kutta_condition — PASS
**Round-1 Must fix**
1. Positional fluidpy calls — **FIXED.** Code cards `kutta` and `edge`: 9 calls by keyword, the two siblings now read
   `kutta_circulation(U=, a=, alpha=, beta=)` and `kutta_circulation_root(U=, b=, e=, beta=, alpha=)`. 29 `py:` rows,
   48 calls, all by keyword, including `….speed(z=)`, `….velocity(x=, y=)`, `….stream(x=, y=)` and the closures
   `…[1](zeta=)`. 38/38 rows ok; both cards run, numbers equal the comments. No tolerance changed.

**New Must fix** — none in the file (L1 seen once at 768×1024 and once at 1366×768).

**Round-1 Should fix — all four done**
- `_fromLoop`: the comment at line 2665 names what the flag does in the library. Still a private API (library item).
  Ramp re-checked: from 0 it reaches exactly 1.000 after 5 s and stays there while playing; 0 errors.
- Explore intro says "▶ from below 1 runs the start-up to the Kutta value (above 1 it only moves the particles)".
- Flat-plate preset: status and title agree, "0.996 U" (three figures everywhere).
- `desktop__derive-d2p3.png`: the amber "α + β" label is back in the narrow circle view.

**What works**
- The hypothesis as an experiment: the amber rest point crawls to the tail, the readout says "unbounded" (never a
  number) until Γ = Γ_K, and the log plot shows the notch at three distances.
- Γ_K two ways, formula and a bisection "that knows no formula"; C_L three ways.
- The inset "the tail, ×9" with the rose patch; the status that says which surface and how far (0.71 % of chord).
- The ramp that ends in a hold card while the particles keep flowing.

## starting_vortex — PASS
**Round-1 Must fix**
1. "90 % after 6.2 chords" — **FIXED.** `crossing()` interpolates between the two bracketing steps; new row
   "90 % crossing (interpolated)": JS 6.129401661376997, Python 6.129401661376972, rtol 1e-9 (21 rows, 21 ok).
   Read off the page for every case (Explain §3 box, both readings, regime line, notes):

   | α | rigid wake: crossing → printed | free wake: crossing → printed |
   |---|---|---|
   | 2° | 6.144 → 6.1 | 6.145 → 6.1 |
   | 5° | 6.129 → 6.1 | 6.136 → 6.1 |
   | 8° | 6.102 → 6.1 | 6.123 → 6.1 |
   | 12° | 6.046 → 6.0 | 6.097 → 6.1 |

   At 12° with the rigid wake every place prints 6.0, consistently; no sentence anywhere states "6.1" as a fixed
   number (all five places take it from `S.s90`; the quiz answer says "About six"). The text reads correctly.

**New Must fix** — none in the file (L1 seen at every size; see the desktop screenshot above).

**Round-1 Should fix — all four done**
- "That 0.55 is quoted from the exact solution, not computed here."
- Rest state has its own sentence ("The foil is at rest. Nothing has been shed … bound 0, shed 0, total 0.").
- Step 1 ends with the question: "Where does it get its circulation?"
- `phone__explore.png`: the circulation graph has ticks 1, 0, −1.

**What works**
- The conserved total flat in amber while teal and rose grow apart; the status repeats "total 0.000" at every instant.
- The cut B–D: one click splits the circuit and each part reads what it holds (the derivation's steps 3–4 on screen).
- The lift panel that says in the plot what Wagner's function is a ratio *of*, and the table "model against Wagner".
- The "from scratch" code card: one step of the march as a 2 × 2 solve with the live numbers.

## airfoil_pressure_lift — PASS
**Round-1 Must fix**
1. Lift numbers missing on portrait phones — **FIXED.** Looked at 360×640 and 390×844 in Explore and in walkthrough
   step 6 (`reports/viz/ch14/_rereview/E4__360x640__tour-step6.png`, `E4__390x844__tour-step6.png`): a teal line
   under the x-axis, "C_L: loop 1.058 · Γ 1.058 · 2π rule 0.987", whole, 12.5 px, clear of the curves; also whole
   with negative values ("loop −0.236 · Γ −0.236 · 2π rule −0.219" at 360 px). On 844×390 it is two rows at the
   plot's top right. The numbers are out of the title, and the Explain hint now says where they are and is true at
   both orientations. 0 audit failures at five phone frames.

**New Must fix** — none in the file (L1 seen at 390×844 and 768×1024).

**Round-1 Should fix — three of four done**
- Tolerances tightened: three rows 1e-5 → 1e-9, three separation rows 2e-2 → 1e-3 (all still ok).
- Equations tab units: "m (here: b)", "m²/s (here: Ub)".
- Step 3: neutral badge "🔺 tall suction peak −4.94 …" (`desktop__tour-step3.png`).
- Left: `desktop__tour-step5.png`, the section uses about 40 % of the view and the lower arrows are 3–4 px.

**Should fix (new, small)**
- 667×375 (a small landscape phone, not one of the eight audit sizes) with negative lift:
  `_rereview/E4__667x375__explore_alpha-6.png` — the first of the two rows is wider than the 164 px view and loses
  its first letter ("_L: loop −0.236 · Γ −0.236"). The numbers are whole. When two rows do not fit, use three or
  drop the "C_L: " prefix.

**What works**
- Three pictures of one flow: arrows on the skin, the loop, and the two term bars that add up to C_y.
- "The part of the lift the loop does not show" — C_L − C_y cos α = −C_x sin α with the numbers: leading-edge suction
  made concrete.
- The inspector line "circle speed × map factor = speed → C_p" for any clicked point.
- A flat-plate ghost whose edge value is "unbounded", checked by two invariant rows.
- New: a number promised to phone readers is drawn where the plot has room (under the axis when the plot is short,
  rows at the top when it is tall), never squeezed.

## lift_curve_stall — PASS
**Must fix** — none in the file (L1 seen on phones; the 390×844 case is new, caused by the added sentence moving a
page break, and goes away with the library fix).

**Round-1 Should fix — all three done**
- Explain §6, looked at on phones (360×640 page 16 of 21; 390×844 page 7 of 9): "This is the section model's maximum
  applied to the whole wing; a finite wing's maximum is lower, and the elliptic-loading explainer gives this same
  aircraft a wing value of C_L,max = 1.4, hence its 26.5 m/s." It sits above the 22.7 m/s line, is whole and
  readable (shot.py reaches the tab after the walkthrough, at α = 5°, and its pager audit is clean there; see L1 for
  what happens once the angle is changed on that tab).
- §4: "an illustrative speed and chord, not the example aircraft's".
- Six `py:` rows now by keyword (`alpha=`, `Re_c=`).

**Should fix (new, small)**
- Code card `model` still has two calls with a leading positional argument: `ch14.stall_regime(Re, Re_laminar=…)`
  and `ch14.generic_lift_curve(al, K=…)`. Write `Re_c=Re` and `alpha=al`, as the rows now do.
- 360×640: the Explain tab is 21 pages of one or two lines each (about 60 px of text under two views). Hiding the
  section view on that tab on small portrait phones would halve it.

**What works**
- "Model" said in the title, in the plot, in the status, on the sketch and in the notes — the pattern for any drawn,
  non-predictive curve.
- The line beyond stall drawn faint with "theory no longer applies"; the zero-lift triangle on the axis.
- The stalled fraction σ written out with the numbers, so the cartoon is at least a transparent one.

## downwash_induced_drag — PASS (unchanged since round 1; lint ok)
**Must fix** — none in the file (L1 seen at every size, all presets).

**Should fix**
- `phone__explore.png`, `phone-land__explore.png`: at N = 3 and 7 the thin singular downwash lines of the staircase
  fill a 94 px plot. Hide them on small views and keep the dots at the panel centres.
- Station exactly at y₁ = ±0.5 with a constant loading: the status reads "upwash 2.78 U outboard of the tip" (0.794 U at
  N = 7). Show "unbounded at the tip filament", as the hygiene rule asks.
- Quiz 4 says "try the loading slider", which changes Γ as well (the answer admits it). Point to the ground-height
  slider: one horseshoe at h/s = 0.25 halves w at the same Γ.
- Code card `down`: the station is shown to four decimals (`y1=0.1429` for 3/21), so the card as written prints
  0.03495 while its comment says 0.03492. Show six decimals, or the fraction.
- The storyboard names `Viz.three`; the view is a canvas projection with orbit. Better than the plan (no WebGL
  dependency) — record it in the knowledge files.

**What works**
- "Build the sum": one horseshoe → 3 → 7 → 21 → ∞, with signed term bars per filament that add up to w.
- The honest staircase: "11 % less than an elliptic loading of the same lift (a coarse staircase underestimates it)".
- Upwash outboard shaded orange with "where a following bird flies"; the factor 0 · 1 · 2 ahead of, at and behind the
  wing, with "sheet only: the bound vortex is left out".
- ε drawn five times too large and labelled so.

## lifting_line_solver — PASS (unchanged since round 1; lint ok)
**Must fix** — none in the file (L1 seen at 390×844 and 844×390).

**Should fix**
- Explain §3 and the inspector list n = 1, 5, 7 at station 3 and skip n = 3 without a word (sin 3γ = 0 at γ = 60°).
  Add "(n = 3: sin 180° = 0)".
- With "perfect twist" on, code card `solve` uses `T(y)` that the *next* card defines ("# T: next block"); pasted in
  order it stops with a NameError. Swap the cards in that mode or define `T` first.
- Code comments cite numbers alone ("# (14.21): Gamma_1 only"). Write the formula in the comment.
- D15 step 10 and D18 step 9: put "(14.19)" and "(14.22)" beside the line, as E1 and E6 do.

**What works**
- The black box opened: one row of the matrix with the reader's numbers, left = right at the clicked station, and the
  condition number.
- Bars on a log scale with the share on top, the teal "lift" bar against the rose "drag only" bars; click a bar to see
  its sine on the span.
- "Make it perfect" by twist, checked against the closed form by two parity rows.
- The N = 1 preset whose badge says the efficiency of 1 "is the truncation speaking".

## elliptic_loading_aspect_ratio — PASS (unchanged since round 1; lint ok)
**Must fix** — none in the file (L1 seen once at 768×1024 and once at 1000×700).

**Should fix**
- Preset "infinite span" sets Λ = 30, where induced drag is still 24 %. Rename it "Λ = 30 (towards infinite span)".
- "Λ doubled" at unchanged speed is badged "fast"; say "above the new best-glide speed (29.0 m/s)".
- "a published estimate … about 37 % induced": name the source (`reference/ch14/SOURCES.md`).

**What works**
- The U-shaped curve with a stacked bar under the dot, a stall wall, and the grey ghost of the wing you started with.
- Trap T9 resolved by the badge: "best glide: induced = profile" / "slow: induced 66 %" / "fast: only 13 %".
- Coefficient and newtons shown as the same drag (C_Di·qA = D_i/e_span), with "span, not area".
- The hollow diamond when the best-glide C_L lies above C_L,max.

## thrust_from_lift — PASS
**Round-1 Must fix**
1. Derivation tab at 844×345, `curve` 59 px — **FIXED.** Re-audited with the round-1 method (shot.py's `VIZ.audit()`,
   text audit and pager audit on every tab, step and derivation page, for sail, fish and bird): 844×345 — 53 views,
   254 pager pages, 0 failures; 667×330 — 53 views, 369 pages, 0 failures; 844×390 — 0 failures. In the Derivation
   tab only `plan` is shown: 331 × 147 px (sail), 103 px (fish), 121 px (bird) at 844×345; 113 / 87 / 87 px at
   667×330. Screenshot `reports/viz/ch14/_rereview/E9__844x345__derive-d1p4.png`: one view, ticks no longer
   stacked, other tabs unchanged (two views side by side, 214 × 165 px).

**New Must fix** — none.

**Round-1 Should fix**
- θ label: done (`phone__explore.png`: θ sits inside its arc, clear of "Uᵣ"). At 844×345 "Uᵣ", "θ" and "drive" are
  close together but each is readable.
- Landscape titles: done ("From above · drive/L +0.36" fits whole).
- Bird at the top of the stroke: done — "thrust −0.6 N · … the lift does not lean: only drag acts along the flight";
  slow downstroke "the lift leans forward too little: drag wins"; lifting upstroke "the lift leans backward". Each
  agrees with the sign printed.
- Two hull models: done, Explain §6 ("Two hull models are in use, and they are not tied together …").
- Efficiency without thrust: "none: no thrust" in the readout, "none, no thrust (V/U < D/L)" in the bars.
- **Withdrawn:** "bird slow-stroke status 'thrust — N · upward force — N'". Not reproducible: bird and fish, three
  upstroke settings, L/D 1, 2, 10, 20, stroke speed 0.05–1.5, phase 0–1 in steps of 0.05 (15 120 states) at
  1366×768, 360×640, 844×390 and 768×1024, plus all eight walkthrough steps — no dash, NaN or undefined in the
  status or the readouts. I cannot give a state for it; treat the round-1 note as my error.

**What works**
- A force polygon that closes: lift ⟂ felt wind, drag along it, drive and side force, the keel's answer.
- No-go (rose) and cannot-hold (amber) zones with three distinct badges: "in irons", "too close to hold", "drive
  positive".
- The same triangle with other owners: fish tail and bird wing, with the upstroke chips "is feathered / keeps lifting
  / lifts downward" as an experiment.
- Selftest rows that test the *drawing* ("drawn drag arrow points to leeward", "drawn lift arrow is square to the felt
  wind").
- A glossary of sailing words in the notes, so no term is used unexplained.
- New: a status phrase chosen from the signs it prints, with a case for "the lift does not lean".

---

## For the orchestrator
- **L1 is the only open Must-fix**, and it is in `assets/viz_lib.js` (above). After it: `tools/viz_inline.py --all`,
  `tools/shot.py --chapter ch14`, and the preset-on-Explain probe; no explainer file of chapter 14 needs a change
  for the gate.
- **Design corrections the builders made, which the notebook should share** (from round 1, still valid). E8 quiz 4:
  halving C_D0 *raises* the best-glide speed by 2^(1/4) (design line 2326 says "down"). E3 step 5: "a third arrives
  at once" (design: "half"). E7 D15 step 11: each term tends to ±n²Γₙ/(2s) at the tips (design: ±n/(2s) Γₙ).
- **`tools/shot.py` gaps:** no state change on the Explain tab (L1), no `.katex-error` probe, no click pass, no
  844×345 size, walkthrough steps are screenshotted only at 390×844 and 1366×768 (E4's step 6 at 360×640 had to be taken by
  hand), no check that numbers promised for phones are visible there, and Code tabs are never executed.
- **Library:** `transport.max` / `holdAt` would remove E2's use of `app._fromLoop`.
- **For `knowledge/viz_patterns.md`:** E3's per-setting number ("6.0 at 12°, 6.1 elsewhere") is right only because
  every sentence reads it from one computed value; a typed "6.1" anywhere would now be wrong. E4's rule for a number
  that must survive on phones: draw it where the plot has room, and say in the Explain hint where that is.

## Verdict: FAIL — five of five round-1 Must-fix items FIXED, no regression in any file, no tolerance loosened; one new Must-fix, L1, in the shared library. With L1 fixed and re-audited the chapter passes without touching the nine files.
