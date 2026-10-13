# Chapter 14 — explainer review, round 4 (gate check of two library fixes) — Verdict: PASS — L1 stays fixed (0 clipped of 288 preset states), L2 is fixed (0 of 1 478 step changes open off page 1; round 3: 45 of 115), L3 is fixed (36 of 36 real-click runs walk 1, 2, 3, … n and end on Explore). 9 / 9 explainers PASS; nothing regressed.            2026-10-13

Round 3 (HEAD `8798d2a`): nine files PASS on their own content; gate held by two defects of the shared library —
L2 (a new step opened on the pager page the reader was on) and L3 (one click on "Next →" advanced two steps). Since
then an uncommitted change to `assets/viz_lib.js` and a re-inline into 121 files. This round treats that as a claim
and re-runs everything. All numbers below were produced in this session.

**Intermittent failures.** `tools/shot.py --chapter ch14` had none (first run, 0 failures, nothing to re-run). The
only failure of the session was outside the chapter: `viz/ch13/rossby_waves.html`, run as an extra regression check
while three other browser jobs were running, failed once (`[phone-land] tour-step2: viz-step-card page 4/4:
viz-inspect 'Column at x=0.375 λ…' clipped 10px`) and then **passed twice in a row** (280 views each) — intermittent,
not counted as a defect; noted under "Library, Should fix (d)".

| Check | Result |
|---|---|
| `git diff -- assets/viz_lib.js` | 3 hunks, 5 + / 4 − lines: `Pager.setItems(list, keepPage)` (line 1115), the "Next →" button built without a handler (line 1214), `calcPager.setItems(items, true)` (line 1242). Nothing else. |
| `git diff -- viz/ch14` | 9 files, 45 insertions, 36 deletions: **the same 5 + / 4 − lines in each file, identical to the library diff** (compared line by line by script). No selftest row, tolerance, text or picture changed. Across the repo: 121 explainer / template files changed, all 121 carry exactly the library diff and nothing else. HEAD `8798d2a`, nothing committed. |
| `tools/viz_lint.py --chapter ch14` | 9 / 9 ok, exit 0 |
| `tools/shot.py --chapter ch14` (8 sizes, every tab, step, derivation page, pager page) | **9 / 9 PASS, 2032 views, 4924 pager pages, 0 failures, 0 CLIP screenshots**, exit 0 (the same counts as round 3) |
| selftest | **328 rows, 328 ok, 240 JS ↔ fluidpy** (same as rounds 2–3). Largest error / tolerance 0.29 (E1 `x_tmax/c measured`, rtol 2e-3), next 0.13 (E1 `t/c measured`, rtol 1e-5); tolerances cannot have changed (see the diff) |
| **L2 probe** (turn a pager page with a real click, then change step; ch14, 9 files) | **0 off page 1** of 534 derivation changes at 360×640, 0 of 534 at 844×390, 0 of 48 at 1366×768; walkthrough 0 of 208, 0 of 150, 0 of 4. Every change also landed on the step asked for (0 wrong) and showed the step's first item; page errors 0 |
| L2 probe, outside ch14 at 360×640: `ch13/rossby_waves` · `ch09/karman_street_stability` · `ch05/point_vortex_lab` | derivation 0 of 86 · 0 of 69 · 0 of 54; walkthrough 0 of 22 · 0 of 22 · 0 of 18 (round 3: 6/20 · 5/17 · 5/12 and 1/5 · 1/4 · 0/4) |
| **L3 probe** (real `mouse.click` on "Next →"; 9 ch14 files + `templates/viz_example.html`, `ch13/rossby_waves`, `ch01/parcel_stability`; 390×844, 1366×768, 844×390) | **36 / 36 runs: steps 1, 2, 3, … n, then "Explore →" lands on Explore**; exactly one `goStep` call per click, none on the last click; 0 defects, 0 page errors |
| L3, the other ways to move (same 36 runs) | "← Back" −1, −1; "Next →" +1 after coming back from Explore; chip 1 → step 1; ArrowRight +1, +1; ArrowLeft −1; two quick clicks +2; `#step=3` opens step 3 and "Next →" gives 4; a touch tap +1; every Derivation tab walks goal → steps → result one page per click and its last button goes to the next derivation or Explore |
| Deep links (`starting_vortex`, `lifting_line_solver`; 390×844 and 1366×768) | `#tab=derive&d=1&ds=3` opens step 3 on page 1, "Next step →" gives step 4 on page 1; `#tab=explain&step=4` opens Explain, the walkthrough behind it is on step 4 and its "Next →" (never shown before) gives step 5; `#tab=tour&step=6` on a 6-step tour shows "Explore →" and it works |
| **L1 probe repeated** (Explain open, each preset, pager audit; 9 files × 6 sizes) | **0 clipped of 288 preset states**; entry 0 / 54; transport drag 0 of 42; 3.5 s of animation 0 of 42 (during and after); `.katex-error` 0, NaN 0, overflow 0, page errors 0; "›" while playing reaches pages 2, 3, 4 in every case with a transport |
| L1, Explain page kept across a preset | kept in E2, E3, E4, E5, E8, E9 at all six sizes (e.g. `kutta_condition` 360×640: 16/16 → 16/16 → 16/16; clamped when the count shrinks: `thrust_from_lift` 25/25 → 24/24); falls to page 1 in E1, E6, E7, whose presets change the number of items — exactly as in round 3 (library Should fix (a)) |
| `.katex-error` probe: every tab, walkthrough step, derivation page, every pager page, Explain + Explore under every preset; 1366×768 and 360×640 | **700 views, 3274 pager pages, 0 `.katex-error`**, 0 raw TeX, 0 NaN / undefined / Infinity, KaTeX loaded in all, 0 page errors |
| Regression, outside ch14: `tools/shot.py templates/viz_example.html` · `viz/ch13/rossby_waves.html` | PASS 200 views, 3 selftest rows · PASS 280 views, 53 rows (2 of 3 runs; see "Intermittent failures") |
| `tools/eq_refs.py --chapter ch14` · `tools/check_public.py viz/ch14/*.html` | 14 hits (the same `ref:` labels and row names as rounds 1–3) · OK (9 files) |

| slug | lint | shot (sizes/views/pages) | parity rows ok | explain | derivations (steps ok) | teaching | polish | L1 (clipped presets, 6 sizes) | L2 (off page 1 / step changes, 3 sizes) | L3 (Next clicks) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| zhukhovsky_airfoil_builder | ok | PASS 8 / 336 / 763 | 35/35 (30 py) | 9 sections + reading | D03 15/15 · D04 11/11 | very good | very good | 0 / 30 | 0 / 272 | 1…7, Explore | **PASS** |
| kutta_condition | ok | PASS 8 / 240 / 562 | 38/38 (29 py) | 9 sections + reading | D07 4/4 · D09 10/10 | very good | very good | 0 / 36 | 0 / 176 | 1…7, Explore | **PASS** |
| starting_vortex | ok | PASS 8 / 144 / 312 | 21/21 (16 py) | 6 sections + reading | D02 5/5 | very good | good | 0 / 30 | 0 / 84 | 1…6, Explore | **PASS** |
| airfoil_pressure_lift | ok | PASS 8 / 168 / 414 | 37/37 (29 py) | 9 sections + reading | D10 7/7 | very good | good | 0 / 30 | 0 / 104 | 1…7, Explore | **PASS** |
| lift_curve_stall | ok | PASS 8 / 96 / 235 | 26/26 (19 py) | 8 sections + reading | none (by design) | very good | very good | 0 / 30 | 0 / 42 | 1…7, Explore | **PASS** |
| downwash_induced_drag | ok | PASS 8 / 304 / 737 | 36/36 (23 py) | 8 sections + reading | D11 4/4 · D12 7/7 · D13 8/8 | very good | good (phone clutter) | 0 / 36 | 0 / 232 | 1…8, Explore | **PASS** |
| lifting_line_solver | ok | PASS 8 / 296 / 773 | 61/61 (47 py) | 10 sections + reading | D15 11/11 · D18 9/9 | very good | very good | 0 / 36 | 0 / 243 | 1…8, Explore | **PASS** |
| elliptic_loading_aspect_ratio | ok | PASS 8 / 264 / 673 | 38/38 (24 py) | 10 sections + reading | D19 9/9 · D21 8/8 | very good | very good | 0 / 30 | 0 / 205 | 1…7, Explore | **PASS** |
| thrust_from_lift | ok | PASS 8 / 184 / 455 | 36/36 (23 py) | 9 sections + reading | D22 8/8 | very good | good | 0 / 30 | 0 / 120 | 1…8, Explore | **PASS** |

Content columns (explain, derivations, teaching, polish) are the round-1/2 judgements; the files are byte-identical
outside the library block, so they were not re-judged. L2 column = walkthrough + derivation, all four ways of changing
step, 360×640 + 844×390 + 1366×768 (sum 1 478).

---

## The diff, judged

- **L2 — `setItems(list, keepPage)`.** `page` is kept only when the caller says so *and* the item count is unchanged;
  `layout()` still clamps. There are three callers in the library and only `updateCalc` (Explain, line 1242) passes
  `true`; `goStep` (line 1707) and `goDerive` (line 1362) pass nothing, which is the pre-`8798d2a` behaviour (page 1).
  No explainer calls `setItems` itself. `goStep` / `goDerive` are called only on a step change, a deep link or at
  start-up — never on a parameter change or an animation tick — so a reader's page inside a step is not reset while
  they drag a slider (live notes and inspector re-pack through `layout()`, which keeps the page).
- **L3 — one handler.** The `h()` helper registers an `on…` attribute with `addEventListener`; `goStep` also assigned
  `tourEls.next.onclick`, so two handlers ran. The button is now built without one, as the Derivation tab's always
  was. *Is it clickable before the first `goStep`?* The start-up code calls `app.goStep(...)` whenever a walkthrough
  exists (line 1888), before the first tab is shown, also when the page opens on another tab — verified by the first
  click after a fresh load in all 36 runs, by the tap, and by `#tab=explain&step=4` → Walkthrough → "Next →".
- **Last step.** Before, the last step's button ran `goStep(n)` (clamped: a needless re-render) and then
  `setTab('explore')`; now only the latter (0 `goStep` calls on the last click, 36 / 36).
- **Not touched and still right:** chips, "← Back", ArrowLeft / ArrowRight (walkthrough and Derivation), "All steps →"
  from a quoted derivation step, the derivation picker, hash links.

## L2 — FIXED (verified)

**The probe** (`reports/viz/ch14/_round4/scripts/r4_keep.py`, the round-3 `r3_behave.py keep` widened): for every
walkthrough step and every derivation page that has a second pager page — click "›" (real click), check the pager
is on page 2, then change step in one of four ways: the panel's "Next →" / "Next step →" (real click), "← Back" (real
click), a numbered chip two steps away (real click), the ArrowRight key. Pass = the new step is the one asked for,
the pager shows page 1 and the first item of the step is visible. Off page 1 / tested:

| | derivation 360×640 | derivation 844×390 | derivation 1366×768 | walkthrough 360×640 | walkthrough 844×390 | walkthrough 1366×768 |
|---|---|---|---|---|---|---|
| "Next" button | 0 / 130 | 0 / 130 | 0 / 5 | 0 / 48 | 0 / 35 | 0 / 1 |
| "← Back" | 0 / 130 | 0 / 130 | 0 / 19 | 0 / 55 | 0 / 41 | 0 / 1 |
| chip | 0 / 144 | 0 / 144 | 0 / 19 | 0 / 57 | 0 / 41 | 0 / 1 |
| ArrowRight | 0 / 130 | 0 / 130 | 0 / 5 | 0 / 48 | 0 / 33 | 0 / 1 |
| **ch14 total** | **0 / 534** | **0 / 534** | 0 / 48 | **0 / 208** | **0 / 150** | 0 / 4 |
| round 3 ("Next" only) | 45 / 115 | 45 / 114 | — | 6 / 48 | 3 / 30 | 0 / 6 |

(More cases than round 3 because "Next →" now moves one step, so every step change is sampled, and because the page
no longer leaks into the probe's own set-up. In two walkthrough steps at 844×390 the second page was not there on
one of the four visits — a live note changes the card's height — so the four rows differ by one or two cases.)
Outside ch14 at 360×640: `rossby_waves` 0 / 86 + 0 / 22, `karman_street_stability` 0 / 69 + 0 / 22,
`point_vortex_lab` 0 / 54 + 0 / 18.

**Looked at with my own eyes** (`reports/viz/ch14/_round4/`, real clicks, `r4_shot.py`):
- `L2r4__starting_vortex__390x844__derive_step4_after_Next_step.png` — the round-3 reproduction (D02 step 3, "›" to
  2 / 2, "Next step →"): step 4 opens at **1 / 2** with its title "The starting vortex carries the opposite
  circulation", then "we had Γ_ABCD = Γ_ABD + Γ_BCD → Name the two halves → now Γ_ABD = Γ_foil, Γ_BCD = Γ_starting →
  Why → In words", all whole. (Round 3: 2 / 2, opening on "your numbers".)
- `L2r4__starting_vortex__360x640__derive_step4_after_Next_step.png` — 1 / 3: title, "we had", the move.
- `L2r4__starting_vortex__390x844__tour_step4_after_Next.png` — walkthrough step 3 at 2 / 2, "Next →": "Step 4 of 6 ·
  Cut the circuit" at 1 / 2 with its text and the cut slider (round 3: opened on the code card).
- `L2r4__zhukhovsky_airfoil_builder__844x390__derive_step4_after_Next_step.png` — landscape, D03 step 3 at 2 / 5 →
  step 4 at 1 / 4 with title and "we had". "← Back" and a chip from a turned page open on page 1 as well (log).

## L3 — FIXED (verified)

**The probe** (`_round4/scripts/r4_dom.py`: `page.mouse.click` at the button's centre, `goStep` wrapped to log its
calls). After each click on "Next →" the step is 1, 2, 3, … n — `starting_vortex` 1–6, `kutta_condition` 1–7,
`downwash_induced_drag` 1–8, `templates/viz_example.html` 1–6, `ch13/rossby_waves` 1–7, `ch01/parcel_stability` 1–8
(round 3: 1, 3, 5, 7) — one `goStep` call per click; the n-th click, on "Explore →", lands on the Explore tab.
36 / 36 file × size runs clean (12 files × 390×844, 1366×768, 844×390). Back, chip, ArrowRight / ArrowLeft, a tap and
Playwright's own `locator.click()` twice in a row each move exactly one step per action.

## L1 — still FIXED

`_round4/scripts/r4_l1.py` (the round-3 `r3_l1.py`, output folder changed): presets with a clipped page / presets

| slug | 360×640 | 390×844 | 844×390 | 768×1024 | 1000×700 | 1366×768 |
|---|---|---|---|---|---|---|
| all nine files | 0 each | 0 each | 0 each | 0 each | 0 each | 0 each |
| **total** | 0/48 | 0/48 | 0/48 | 0/48 | 0/48 | 0/48 |

Transport drag 0 / 42, running animation 0 / 42, entry 0 / 54, `.katex-error` 0, NaN 0, overflow 0, page errors 0.
The Explain page survives a preset that keeps the item count (6 files × 6 sizes, last page stays last page); "›"
while playing reaches 2, 3, 4. Seen: `_round4/L1r3__starting_vortex__1366x768__preset-stop__pg4of4.png` — "Reading
the current setting" whole on 4 / 4 (the file name keeps the script's old prefix; the picture is from this round).

---

## Per explainer — Must fix: none in any file. Should fix (carried forward; files unchanged)

## zhukhovsky_airfoil_builder — PASS
**Must fix** — none. **Should fix**
- 1366×768, Derivation tab: both view titles are cut ("centre μ = −0.100 + 0.155 i,…", "t/c 12.0 % (rule 13.…");
  `medium()` should test the view's width, not the window. (Also at 844×390: "Circle · μ (0.000, 0.14…".)

**What works** — a draggable centre with the guides "up = β", "back = e"; rule-against-measured bars with the gap in
words; 3√3/4 pinned by a parity row; the D03 check a reader can do by hand.

## kutta_condition — PASS
**Must fix** — none. **Should fix**
- `app._fromLoop` is still a private library flag (library: `transport.max` / `holdAt` would remove it).

**What works** — the hypothesis as an experiment (rest point crawls to the tail, "unbounded" until Γ = Γ_K); Γ_K two
ways; the "tail ×9" inset; the ramp that ends in a hold card.

## starting_vortex — PASS
**Must fix** — none. **Should fix**
- Derivation D02 step 4 at 390×844: the label "Γ bound +0.92 Γ_K" runs under the cut B–D and its last characters are
  hidden by the dashed line (still so in `_round4/L2r4__starting_vortex__390x844__derive_step4_after_Next_step.png`).
  Move it left of the cut.

**What works** — the conserved total flat in amber; the cut B–D; the lift panel that says what Wagner's function is a
ratio of; the per-setting 90 % crossing read from one computed value everywhere.

## airfoil_pressure_lift — PASS
**Must fix** — none. **Should fix**
- `desktop__tour-step5.png`: the section uses about 40 % of the view and the lower arrows are 3–4 px.
- 667×375 with negative lift: the first lift row is wider than the 164 px view and loses its first letter.

**What works** — three pictures of one flow; leading-edge suction as "the part of the lift the loop does not show";
the inspector line; the phone lift line drawn where the plot has room.

## lift_curve_stall — PASS
**Must fix** — none. **Should fix**
- Code card `model`: two calls with a leading positional argument (`stall_regime(Re, …)`, `generic_lift_curve(al, …)`).
- 360×640: the Explain tab is 21–22 pages of one or two lines; hide the section view on that tab on small phones.

**What works** — "model" said in the title, plot, status, sketch and notes; the faint line beyond stall; σ written out.

## downwash_induced_drag — PASS
**Must fix** — none. **Should fix**
- Phones, N = 3 and 7: the thin singular downwash lines fill a 94 px plot; keep only the dots at the panel centres.
- Station exactly at the tip with constant loading: status prints "upwash 2.78 U"; show "unbounded at the tip filament".
- Quiz 4 says "try the loading slider" (changes Γ too); point to the ground-height slider.
- Code card `down`: station shown to four decimals, so the pasted card prints 0.03495 under a comment saying 0.03492.
- The transport readout for the station reads "0 s" (span units), which looks like seconds; write "0·s" or "y₁/s = 0".
- Storyboard names `Viz.three`; the view is a canvas projection (better: no WebGL) — record it in the knowledge files.

**What works** — "build the sum" 1 → 3 → 7 → 21 → ∞ with signed term bars; the honest staircase; upwash outboard;
ε drawn ×5 and labelled so.

## lifting_line_solver — PASS
**Must fix** — none. **Should fix**
- Explain §3 and the inspector skip n = 3 at station 3 without a word; add "(n = 3: sin 180° = 0)".
- With "perfect twist" on, code card `solve` uses `T(y)` defined in the next card (NameError when pasted in order).
- Code comments cite numbers alone ("# (14.21): Gamma_1 only"); write the formula.
- D15 step 10, D18 step 9: put "(14.19)" and "(14.22)" beside the line.

**What works** — the matrix row with the reader's numbers; log bars with shares; "make it perfect" by twist; the
N = 1 badge "the truncation speaking".

## elliptic_loading_aspect_ratio — PASS
**Must fix** — none. **Should fix**
- Preset "infinite span" sets Λ = 30 (induced drag still 24 %); rename "Λ = 30 (towards infinite span)".
- "Λ doubled" at unchanged speed is badged "fast"; say "above the new best-glide speed (29.0 m/s)".
- "a published estimate … about 37 % induced": name the source (`reference/ch14/SOURCES.md`).

**What works** — the U-shaped curve with a stacked bar, stall wall and ghost wing; trap T9 resolved by the badge;
coefficient and newtons as the same drag; the hollow diamond.

## thrust_from_lift — PASS
**Must fix** — none. **Should fix** — none open (at 844×345 "Uᵣ", "θ" and "drive" sit close together but are readable).

**What works** — a force polygon that closes; no-go and cannot-hold zones with three badges; the same triangle for
fish and bird; selftest rows that test the drawing; the sailing glossary.

---

## For the orchestrator
- **Open Must-fix: none.** L1, L2 and L3 are verified fixed; the gate is open. The library change and the 121
  re-inlined files are uncommitted — commit them together (`viz_lib: …; re-inlined`), then this report.
- **Library, Should fix:** (a) the Explain page falls to 1 when a preset changes the number of items (E1, E6, E7,
  `rossby_waves`, `karman_street_stability`, `point_vortex_lab`) — keep by nearest section heading instead of by
  count; (b) `queueLive`'s `unseen()` re-fits 3–5 times a second on the Equations tab and on a walkthrough step with a
  live formula while playing (no frame loss measured in round 3: 58–60 fps); re-fit only when the visible size
  changes, and settle once on pause; (c) `transport.max` / `holdAt` for E2; (d) **new, to watch (not ch14)**:
  `ch13/rossby_waves`, 844×390, walkthrough step 2 — the live inspector line on card page 4 / 4 was measured 10 px
  clipped in one `tools/shot.py` run of three, under heavy CPU load, while the animation runs; it passed the two
  re-runs and in round 3. If it shows again in the ch13 audit, the re-pack after a live inspector changes height is
  the place to look.
- **`tools/shot.py` gaps** (carried forward; the file is unchanged): no state change on the Explain tab (the L1
  probe), no `.katex-error` probe, no click pass — **in particular no real click on "Next →" (L3) and no "turn a
  page, change step" check (L2)** — no 844×345 size, walkthrough steps screenshotted only at 390×844 and 1366×768, no
  check that numbers promised for phones are visible there, Code tabs never executed. `r4_dom.py` and `r4_keep.py`
  (`reports/viz/ch14/_round4/scripts/`, git-ignored) are ready to be folded in: one real click on "Next →" must raise
  `VIZ.app.step` by exactly 1; after "›" and a step change the pager must show page 1.
- **Design corrections the builders made, which the notebook should share** (from round 1, still valid). E8 quiz 4:
  halving C_D0 raises the best-glide speed by 2^(1/4). E3 step 5: "a third arrives at once". E7 D15 step 11: each
  term tends to ±n²Γₙ/(2s) at the tips.
- **For `knowledge/viz_patterns.md`:** E3's per-setting number read from one computed value; E4's rule for a number
  that must survive on phones; *a reviewer drives the walkthrough with real clicks at least once* — three rounds of
  `VIZ.goStep` audits hid a button that skipped every second step; and, from this round: *a pager keeps its page only
  when the caller asks for it* (`setItems(list, true)` for a live refresh of the same content; a new step, a new
  derivation page always open on page 1), and *a button built with `h(…, {onclick})` must not also get `.onclick`*.

## Verdict: PASS (all nine explainers PASS) — L1 0 / 288 presets (drag 0 / 42, play 0 / 42, Explain page kept); L2 0 / 1 478 step changes off page 1 in ch14 and 0 / 271 in three older explainers (round 3: 45 / 115 and 45 / 114); L3 36 / 36 real-click runs 1, 2, 3, … n → Explore; `tools/shot.py` 9 / 9 PASS, 2032 views, 4924 pager pages, 0 failures; 328 / 328 selftest rows; 0 `.katex-error` in 700 views; every explainer's diff is the library hunk and nothing else. No Must-fix open.
