# Chapter 14 — explainer review (round 1)                                  2026-10-12

## Verdict: FAIL (4 of 9 explainers PASS; 5 carry one Must-fix item each — all five are small)

Full review from scratch. The `audit.json` files left by the interrupted run were not used: `tools/viz_lint.py` and
`tools/shot.py` were re-run here (audits rewritten 16:47–16:57), and everything below comes from that run and from the
reviewer's own probes.

| Check | Result |
|---|---|
| `tools/viz_lint.py --chapter ch14` | 9 / 9 ok |
| `tools/shot.py --chapter ch14` (8 sizes, every tab, step, derivation page and pager page) | 9 / 9 PASS, 2032 views, exit 0 |
| selftest | 327 rows, 327 ok, 239 JS ↔ fluidpy |
| KaTeX | loaded at every size; **0 `.katex-error`**, 0 text fallback, 0 control characters in any TeX string — every tab, every walkthrough step, every derivation page at 1366×768 and 360×640, Explain under every preset, plus 59 mode / extreme states × 2 sizes |
| NaN / undefined / Infinity in rendered text, JS errors | 0 in the default state, under every preset and in the 59 extreme states (118 state × size cases) |
| click-and-drag pass over every view, two sizes | 0 errors |
| in-page landscape frame 844×345 (shot.py's own checks, 254 views, 1242 pager pages) | 8 PASS, **1 FAIL** (`thrust_from_lift`, Derivation tab) |
| Code tabs executed in Python with the live values substituted | 19 / 19 cards run without error (later cards reuse the first card's variables); every printed value equals its comment except one fourth-figure difference in E6 (Should fix) |
| derivations against design Part F | 14 derivations, 116 steps; step counts and move titles equal (one title reworded, D15 step 10) |
| `tools/eq_refs.py` | 14 hits, all ref labels or selftest row names; a second scan of the rendered text found no bare "(N.M)" in tour, Explain, Derivation, notes, status or quiz |
| equations against the rendered book pages (printed 704, 705, 707, 711, 713, 714, 715) | (14.4)–(14.8), (14.9)–(14.12), (14.13), (14.15)–(14.17), (14.18)–(14.26) all match |
| `tools/check_public.py viz/ch14/*.html` | OK (9 files) |

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| zhukhovsky_airfoil_builder | ok | PASS 8 / 336 | 35/35 (30 py) | 9 sections + reading | D03 15/15 · D04 11/11 | very good | very good | **FAIL** (1) |
| kutta_condition | ok | PASS 8 / 240 | 38/38 (29 py) | 9 sections + reading | D07 4/4 · D09 10/10 | very good | very good | **FAIL** (1) |
| starting_vortex | ok | PASS 8 / 144 | 20/20 (15 py) | 6 sections + reading | D02 5/5 | very good | good | **FAIL** (1) |
| airfoil_pressure_lift | ok | PASS 8 / 168 | 37/37 (29 py) | 9 sections + reading | D10 7/7 | very good | good (phone) | **FAIL** (1) |
| lift_curve_stall | ok | PASS 8 / 96 | 26/26 (19 py) | 8 sections + reading | none (by design) | very good | very good | **PASS** |
| downwash_induced_drag | ok | PASS 8 / 304 | 36/36 (23 py) | 8 sections + reading | D11 4/4 · D12 7/7 · D13 8/8 | very good | good (phone clutter) | **PASS** |
| lifting_line_solver | ok | PASS 8 / 296 | 61/61 (47 py) | 10 sections + reading | D15 11/11 · D18 9/9 | very good | very good | **PASS** |
| elliptic_loading_aspect_ratio | ok | PASS 8 / 264 | 38/38 (24 py) | 10 sections + reading | D19 9/9 · D21 8/8 | very good | very good | **PASS** |
| thrust_from_lift | ok | PASS 8 / 184 | 36/36 (23 py) | 9 sections + reading | D22 8/8 | very good | good (labels) | **FAIL** (1) |

Depth against the references: all nine Explain tabs have the structure of `forced_damped_vibrations.html` (what the views
show with colour words → every quantity from the controls as formula = substituted = result → the values now → a
regime-dependent reading that names its thresholds) and more sections than its seven. None is shallower than the
reference its storyboard names.

### The eight points raised by the build session
1. **KaTeX errors** — none anywhere (see the table). One cosmetic loss: an unescaped `%` inside TeX in E1 (Should fix).
2. **E3 "90 %"** — the page *computes* it (line 2482) but as the first grid step at or above 0.9, 6.15 chords, which
   `toFixed(1)` prints as "6.2". The verified value at its own time step is 6.13 → "6.1". **Must fix** (E3).
3. **Wagner's function** — only E3 mentions it; it is drawn only in the "Lift / steady lift" view, titled "dotted:
   Wagner (a lift ratio)" and labelled "circulatory lift / steady lift (Jones' fit)"; Explain §4 says it starts at one
   half while the exact circulation starts at zero. Nothing on a circulation axis. Clean.
4. **E2 ramp through `app._fromLoop`** — behaves: from 0 it reaches exactly 1.000 after 5 s and holds there with the
   particles running (also at 3× speed and at 360×640); pause → ▶ stays at 1; Restart ramps again; scrubbing to 0.8 while
   playing resumes the ramp; scrubbing to 1.3 holds 1.3; step buttons 0.95 → 1 → 1.05 → 1.1 → 1.05; `app.reset()` clears
   the hold; walkthrough steps 1 and 4 ramp and hold, 2, 3, 5, 6, 7 are still; going back from step 5 to 4 ramps again;
   0 errors. Still a private API (Should fix / library).
5. **E5 model** — labelled everywhere: view title "rose: a drawn model", in-plot "rose: a drawn model, not data" and
   "stall: C_L,max = 1.90 at 15.1° (model)", status "… — model", section view "separated flow — a sketch, not a
   computation", notes "the book's measured values are not reproduced". Its constants (K = 5.7, stall at 11° / 15°,
   plateau) are ours; none is the book's figure value.
6. **E6** — no "first order" claim; notes and reading say "the error falls roughly as N^(−1.5)", pinned by a selftest
   row (1.512) and reproduced in Python here (1.51–1.54). The 3-D view never uses WebGL: it is an orbitable projection
   on the 2-D canvas, so it draws identically with WebGL disabled (tested: canvas inked, drag orbits, 0 errors).
7. **Keyword calls** — E3–E9: every fluidpy call in Code tabs and `py:` rows is by keyword (a single leading positional
   `alpha` / `Re` in three E5 rows). **E1 and E2 are positional throughout — Must fix.**
8. **E4–E9** — reviewed in full, as E1–E3.

---

## zhukhovsky_airfoil_builder
**Must fix**
1. **Positional fluidpy calls** (build rule: sibling functions have different positional orders). Code card `shape`,
   lines 2992–2998: `ch14.zhukhovsky_circle(b, e, beta)`, `zhukhovsky_airfoil(b, e, beta)`, `zhukhovsky_shape(b, e,
   beta)` (twice), `map_point(ξ, η, b)`; and 32 positional calls in the 30 `py:` rows (`zhukhovsky_circle(1.0, 0.08, …)`,
   `zhukhovsky_circle_from_centre(mu, 1.0)`, `joukowski(z, 1.0)`, `map_point(0.62, 0.95, 1.0)`, `zhukhovsky_airfoil(1.0,
   0.1, …, 400)`, `zhukhovsky_shape(…)`, `circular_arc(1.0, …)`). Rewrite all with keywords (`b=, e=, beta=, n=, mu=,
   zeta_re=, zeta_im=, exact=`); re-run `tools/shot.py` (35/35 must stay ok). No number changes: every row passes and
   the card runs as written.

**Should fix**
- D03 "What it means right now": `$\tfrac12\beta=6.98 %$` — a bare `%` is a KaTeX comment, so the unit is dropped.
  Write `6.98\,\%`.
- D04 step 5 *your numbers* compares the β = 0 formula with the current contour. If the reader has set β ≠ 0 the two
  disagree badly (0.0064 against 0.0561 at β = 8°). Evaluate both at β = 0 or label the line "(β = 0)".
- `phone-tall__derive-d1p7.png`: the triangle O–Q–P that the *Watch* line points to is hidden under "μ · drag me" and
  "+b: the tail". Drop those labels in D03 steps 6–8 on small views.
- `notebook__equations.png`: both view titles are cut ("centre μ = …", "camber 0.0…"). Shorten them at this width.
- Tolerances wider than needed: the "measured" rows declare 2e-3 (5e-3 for x_tmax) and agree to 1e-6 (6e-4). Use 1e-5
  and 2e-3.

**What works** (keep; candidates for `knowledge/viz_patterns.md`)
- A draggable centre with two guides ("up = β", "back = e") that *is* the idea; presets are the book's four cases plus
  "too fat".
- Rule-against-measured bars with the gap in words, and a small curve of the rule leaving the truth as e grows.
- The constant 3√3/4 pinned by a parity row ("computed, not typed").
- The D03 check with a number the reader can do by hand (tan β = ¾ → 1.5 = −0.583 + 2.083).

## kutta_condition
**Must fix**
1. **Positional fluidpy calls.** Code cards `kutta` and `edge`, lines 2985–3016: nine calls, among them
   `kutta_circulation(1.0, a, al, be)` next to `kutta_circulation_root(1.0, 1.0, 0.08, be, al)` — the angle order is
   reversed between the two siblings, which is exactly the trap the rule exists for — and `stagnation_points`,
   `edge_speed`, `trailing_edge_speed`, `lift_per_span`, `zhukhovsky_lift_coefficient`. Also 41 positional calls in the
   29 `py:` rows (lines 3070–3105, including the shared strings `PC`, `FL`, `CF`). Rewrite all with keywords; re-run
   `tools/shot.py` (38/38 must stay ok). No number changes.

**Should fix**
- `pin()` (line 2666) and `onChange` (line 3231) use the private `app._fromLoop`. It works (point 4 above); replace it
  when the library gets `transport.max` / `holdAt`, and until then say in a comment which library lines it relies on.
- ▶ pressed with Γ/Γ_K > 1 holds the value instead of ramping; say so in the Explore intro ("▶ from below 1 runs the
  start-up; above 1 it only moves the particles").
- Flat-plate preset: the status says "leaves the edge smoothly at 1 U", the view title 0.996 U. Use three figures in
  the status.
- `desktop__derive-d2p3.png`: in the narrower circle view of the Derivation tab the labels α, β, α + β are dropped, yet
  the *Watch* line points to "the amber arc … down to B". Keep the α + β label there.

**What works**
- The hypothesis as an experiment: the amber rest point crawls to the tail, the readout says "unbounded" (never a
  number) until Γ = Γ_K, and the log plot shows the notch at three distances.
- Γ_K two ways, formula and a bisection "that knows no formula"; C_L three ways.
- The inset "the tail, ×9" with the rose patch; the status that says which surface and how far (0.71 % of chord).
- The ramp that ends in a hold card while the particles keep flowing.

## starting_vortex
**Must fix**
1. **"90 % after 6.2 chords" is not the verified number.** Line 2482 takes the first grid step with ratio ≥ 0.9
   (k = 122, 6.15 chords at dt = 0.05 c/U) and `f1()` rounds it to "6.2". The crossing by linear interpolation is
   6.13 chords (`tests/test_ch14.py` lines 724–726; `reports/ch14_verification.md` "6.1 chords, not 6.2"). The page
   shows it in four places — Explain §3 box (line 2817), both readings (lines 2692–2693), the notes, and the regime
   line (line 2830) — and prints 6.2 at α = 2°, 5°, 8° but 6.1 at α = 12°, a rounding artefact. Change line 2482 to
   interpolate, `s90 = s[k−1] + (0.9 − r[k−1])/(r[k] − r[k−1])·DT`, so that it prints 6.1, and pin it with a row such as
   `{name: '90 % crossing', js: S.s90, py: '(lambda r, np=np: float(np.interp(0.9, r["ratio"], r["s_chords"])))(ch14.impulsive_start_lumped_vortex(U=1.0, c=1.0, alpha=np.deg2rad(5.0), t_end=12.0, dt=0.05))', rtol: 1e-9}` (6.129).

**Should fix**
- Explain §3: "0.65 against about 0.55 after one chord" — 0.55 (the exact thin-airfoil circulation) is typed. Say
  "quoted from the exact solution, not computed here".
- Notes at s = 0: "the foil already carries 0.00 of its final circulation". Give the rest state its own sentence.
- Step 1 is a command ("Press ▶ and watch …"). End it with the question: "Where does it get its circulation?"
- `phone__explore.png`: the circulation graph is 76 px high with only the 0 tick; add ±1.

**What works**
- The conserved total flat in amber while teal and rose grow apart; the status repeats "total 0.000" at every instant.
- The cut B–D: one click splits the circuit and each part reads what it holds (the derivation's steps 3–4 on screen).
- The lift panel that says in the plot what Wagner's function is a ratio *of*, and the table "model against Wagner".
- The "from scratch" code card: one step of the march as a 2 × 2 solve with the live numbers.

## airfoil_pressure_lift
**Must fix**
1. **On portrait phones the three lift numbers are not on screen, and the Explain tab says they are.** The bars view
   is `hidePortrait`. The in-plot line "C_L: loop … · Γ … · 2π rule …" is drawn only when the plot is not `low`
   (line 2732: `if (barsHidden && !low)`), so it never appears at 360×640 or 390×844; and the title appends the numbers
   *after* the legend words (line 2853), where the ellipsis cuts them. `phone-tall__tour-step6.png` ("Three roads to
   the lift") and `phone__explore.png` show neither. Yet line 2894 tells the phone reader "its three lift numbers are
   written inside the graph". Put the numbers first in the `cp` title on non-wide layouts ("C_L loop 1.058 · Γ 1.058 ·
   2π 0.987", legend words after), or draw the line when `low` as well; make the hint say where they really are.

**Should fix**
- `desktop__tour-step5.png`: the section uses about 40 % of a 923 px view and the lower-surface arrows are 3–4 px.
  Scale the section (or the arrows) to the view.
- Tolerances wider than needed: "CL pressure route", "centre of pressure", "thick symmetric: CL" declare 1e-5 and agree
  to 4e-11; the three separation rows declare 2e-2 and agree to ≤ 1.1e-4. Use 1e-9 and 1e-3.
- Equations tab: the unit column gives "m" for ζ, z, b while the page works in units of b.
- Step 3: the badge is a green ✅ at a peak of −4.94; a neutral tone would match "a real section would be close to
  stall".

**What works**
- Three pictures of one flow: arrows on the skin, the loop, and the two term bars that add up to C_y.
- "The part of the lift the loop does not show" — C_L − C_y cos α = −C_x sin α with the numbers: leading-edge suction
  made concrete.
- The inspector line "circle speed × map factor = speed → C_p" for any clicked point.
- A flat-plate ghost whose edge value is "unbounded", checked by two invariant rows.

## lift_curve_stall
**Must fix** — none.

**Should fix**
- The stall speed uses the section model's C_L,max = 1.90 for the whole aircraft (22.7 m/s), while
  `elliptic_loading_aspect_ratio` gives the same aircraft C_L,max = 1.4 and 26.5 m/s. Add one sentence in Explain §6
  ("the other explainer uses a wing value of 1.4") or use the same number.
- Explain §4: Re = 60 × 1.5 / 1.5×10⁻⁵ belongs to no aircraft on the page; call it illustrative.
- Three `py:` rows pass the first argument positionally (`generic_lift_curve(np.deg2rad(4.0), K=…)`, `stall_regime(2e5,
  …)`); write `alpha=`, `Re=`.

**What works**
- "Model" said in the title, in the plot, in the status, on the sketch and in the notes — the pattern for any drawn,
  non-predictive curve.
- The line beyond stall drawn faint with "theory no longer applies"; the zero-lift triangle on the axis.
- The stalled fraction σ written out with the numbers, so the cartoon is at least a transparent one.

## downwash_induced_drag
**Must fix** — none.

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

## lifting_line_solver
**Must fix** — none.

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

## elliptic_loading_aspect_ratio
**Must fix** — none.

**Should fix**
- Preset "infinite span" sets Λ = 30, where induced drag is still 24 %. Rename it "Λ = 30 (towards infinite span)".
- "Λ doubled" at unchanged speed is badged "fast"; say "above the new best-glide speed (29.0 m/s)".
- "a published estimate … about 37 % induced": name the source (`reference/ch14/SOURCES.md`).

**What works**
- The U-shaped curve with a stacked bar under the dot, a stall wall, and the grey ghost of the wing you started with.
- Trap T9 resolved by the badge: "best glide: induced = profile" / "slow: induced 66 %" / "fast: only 13 %".
- Coefficient and newtons shown as the same drag (C_Di·qA = D_i/e_span), with "span, not area".
- The hollow diamond when the best-glide C_L lies above C_L,max.

## thrust_from_lift
**Must fix**
1. **Derivation tab at the in-page landscape frame 844×345: the `curve` view is 59 px high** (engine floor 60) on all
   ten pages — `VIZ.audit()` reports `view-too-small:curve`, 48 failures over the pager pages; 45 px at 667×330. The
   long three-line status and the mode chips leave 93 px for `plan` and a strip for `curve`, whose y ticks 1, 0, −1
   print on top of each other. The other eight explainers pass at this frame, and so do this one's other tabs. In the
   Derivation tab on landscape show only `plan` (the derivation's `view`), or use the short status there as on
   portrait; then re-audit at 844×345.

**Should fix**
- `phone__explore.png`, `phone-tall__tour-step2.png`: the label "U_r" sits under the θ arc; at 844×390 the plan title
  is cut inside its TeX ("closest …"). Move the label; shorten the title.
- Bird, slow stroke: the status reads "thrust — N · upward force — N". Say why there is no value.
- Bird at the top of the stroke: "the lift leans forward" with thrust −0.6 N. Choose the phrase by the sign.
- Explain §6: the steady speed uses a fixed hull drag area (0.05 m²) and the closest course uses the hull's L/D slider.
  Say that these are two hull models.

**What works**
- A force polygon that closes: lift ⟂ felt wind, drag along it, drive and side force, the keel's answer.
- No-go (rose) and cannot-hold (amber) zones with three distinct badges: "in irons", "too close to hold", "drive
  positive".
- The same triangle with other owners: fish tail and bird wing, with the upstroke chips "is feathered / keeps lifting
  / lifts downward" as an experiment.
- Selftest rows that test the *drawing* ("drawn drag arrow points to leeward", "drawn lift arrow is square to the felt
  wind").
- A glossary of sailing words in the notes, so no term is used unexplained.

---

## For the orchestrator
- **Design corrections the builders made, which the notebook should share.** E8 quiz 4: halving C_D0 *raises* the
  best-glide speed by 2^(1/4) (design line 2326 says "down"). E3 step 5: "a third arrives at once" (design: "half").
  E7 D15 step 11: each term tends to ±n²Γₙ/(2s) at the tips (design: ±n/(2s) Γₙ).
- **`tools/shot.py` gaps met again:** no `.katex-error` probe, no click pass, no 844×345 size (it would have caught
  the E9 item), no check that numbers promised for phones are visible there, and Code tabs are never executed. The
  reviewer's probes for all five are throw-away scripts.
- **Library:** `transport.max` / `holdAt` would remove E2's use of `app._fromLoop`.

## Verdict: FAIL — fix the five Must-fix items (E1, E2: keyword calls · E3: 6.1 · E4: phone lift numbers · E9: 844×345 Derivation), then re-review those five files
