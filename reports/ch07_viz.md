# Chapter 7 — explainer review                                  2026-09-24 (round 2, after fix round 1)

**Round 2 (this update).** I re-ran all the gates on the four explainers that failed round 1:
- `tools/viz_lint.py`: 4/4 `ok`.
- Full `tools/shot.py`: 4/4 PASS. `capillary_gravity_waves` 8 sizes / 264 views / 7 steps; `group_velocity_packets`
  8 / 392 / 8 steps (the pond step is new); `wave_rays_refraction` 8 / 336 / 6; `internal_wave_beams` 8 / 432 / 7.
- `tools/eq_refs.py --chapter ch07`: 3 hits, all in `internal_wave_beams`. All three are `ref:` badges: the D36 header
  L3026, and the Code-tab badges L3064 and L3076. The two D33 step titles that relied on the number are renamed
  ("Apply ∇_H² to the second link", "…to the first link"). I accept all three hits.

I then ran `tools/shot.py --chapter ch07 --quick` on all 9 explainers:
- 8 PASS.
- `hydraulic_jump` gave **one intermittent FAIL**: `[phone-tall] tour-step2: viz-step-card page 2/4: viz-terms 'in: ρQu₁ 883 N/m push on fac…' clipped 3px`.
- A full 8-size run of `hydraulic_jump` straight afterwards passed (280 views), and so did two more `--quick` runs. The
  file has not changed since 11:09. The failure is a timing-dependent 3 px packing margin, not a steady failure (see
  hydraulic_jump, Should fix 3).

**Parity: 222/222 selftest rows `ok`, of which 172 are JS ↔ fluidpy `py:` rows.** New rows:
- `wave_rays_refraction`: three island rows (x, y, k_y of an island ray at t = 200 s against `ch07.ray_trace` with the
  smoothed-cone depth, rtol 1e-8/1e-7) and an ω-drift invariant for the island bundle.
- `group_velocity_packets`: a capillary energy invariant (2 × mean surface energy ½σ⟨η_x²⟩ = ½σk²a², midpoint rule,
  rtol 1e-9).

**Rule 9 re-scan.** I grepped the four files and compared them with `tests/book_values_ch07.json`:
- "17.8", "4.4 cm", "4.3 cm", "4.39", "23.1 cm/s at 1.71", "19.4" and "1.20 cm" no longer appear in any hard-coded
  string.
- Both explainers now use σ = 0.07274 N/m and ρ = 998.2 kg/m³. They build c_g,min from `cgmin()` / `CGM` and print 4
  significant figures: **17.76 cm/s at 4.354 cm** (c_min 23.12 cm/s at 1.712 cm).
- The one remaining match is the 3-s.f. live label of the reader's own wave. With λ = λ_m it shows "c = 23.1 cm/s"
  (Explore default: the orange dot, the status and the Explain §3 box). This is a computed value of a numbered equation
  (7.58), not a transcription (Should fix below).
- "Fig. 7.29" survives only as a parenthetical in the D36 *why*, as I recommended. A JS source comment says
  "Fig. 7.31: c − c_g is horizontal". That is a pointer, not reproduced content.

**Round 1 (unchanged, still valid for the five that passed):** 9/9 lint `ok`; 9/9 PASS; 29 D ids match the notebook's
step counts; D30 checked with sympy; D05–D36 checked by hand against the rendered pages.

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| dispersion_relation | ok | PASS 8/528 (quick 4/264 ✓) | 23/23 (20 py) | 8 §, live, interp | D01 5 · D05 11 · D06 7 · D07 5 · D08 7 · D09 7 — ok (one rounding slip) | good | good; phone c(λ) squashed | **PASS** |
| particle_orbits | ok | PASS 8/336 (quick 4/168 ✓) | 31/31 (24 py) | 8 §, live, interp | D10 9 · D11 6 · D27 9 — ok | good | phone-land surface is a coarse polyline | **PASS** |
| capillary_gravity_waves | ok | **PASS 8/264** (r2) | 18/18 (15 py) | 8 §, live, interp | D15 10 · D16 7 — ok | good | branch label now clear; minor phone label overlap | **PASS** (r1 Must-fix resolved) |
| seiche_standing_waves | ok | PASS 8/224 (quick 4/112 ✓) | 23/23 (19 py) | 7 §, live, interp | D17 8 · D18 5 — ok | good | good | **PASS** |
| group_velocity_packets | ok | **PASS 8/392** (r2) | 26/26 (21 py) | 8 §, live, interp | D19 7 · D20 12 · D21 11 — ok | good; pond step added | beats labels fixed | **PASS** (r1 Must-fix resolved) |
| wave_rays_refraction | ok | **PASS 8/336** (r2) | 23/23 (18 py) | 8 §, live, interp | D22 8 · D23 8 · D24 9 — ok | good; island step now synced | phone contour labels overlap | **PASS** (r1 Must-fix resolved) |
| hydraulic_jump | ok | PASS 8/280 (quick: 1 intermittent 3 px clip, 3 re-runs PASS) | 24/24 (16 py) | 7 §, live, interp | D25 10 · D26 8 — ok | very good | step 2 now 4 pages at 390×844 | **PASS** |
| two_layer_modes | ok | PASS 8/520 (quick 4/260 ✓) | 27/27 (19 py) | 8 §, live, interp | D28 10 · D29 12 · D30 12 (sympy ok) · D31 11 — ok | very good | good | **PASS** |
| internal_wave_beams | ok | **PASS 8/432** (r2) | 27/27 (20 py) | 9 §, live, interp | D33 12 · D34 8 · D35 5 · D36 9 — ok | good | phone views now legible | **PASS** (r1 Must-fix resolved) |

---

## dispersion_relation
**Must fix** — none.

**Should fix**
1. D09 step 4 (`dispersion_relation.html` L2991): $\sqrt{\tanh kH/kH}=0.969$ is truncated. The value is 0.96961, so
   write **0.970**. The *why* line's "0.968" for the two-term estimate is right (0.96776).
2. `phone__explore.png`: at 360×640 the c(λ) view is about 80 px tall and carries only the tick "1". The reader sees a
   dot on a line, not the deep branch rising into the √(gH) plateau. Give `curve` more of the portrait stage in Explore
   (or drop the help line under the λ slider on phones) so that at least two decades of c are labelled.
3. `phone-tall__derive-d2p9.png` (and other long steps): page 1 of a step shows *we had*, the move title and about
   180 px of empty card; the *now* line only appears on page 2. The pager makes this acceptable, but the D05 regrouping
   steps would read better if the `aligned` *we had* showed only the line the move acts on.

**What works**
- Three linked views on one clock, with ghost deep and shallow asymptotes and shaded 2 %/3 % regime bands.
- The sea-bed banner carries the hidden profile's number onto phones.
- Presets sit at the genuinely special cases.
- Explain works k → ω → c → c_g → bottom pressure with the reader's numbers and gives a regime-dependent reading.
- The `k from T` parity rows prove the shoaling preset.

## particle_orbits
**Must fix** — none.

**Should fix**
1. `phone-land__explore.png`: η is drawn as a kinked polyline. Sample it at ≥ 2 px spacing in every layout.
2. Phone pages: Explain needs 18 pages and tour step 6 needs 7 pages at 360×640. Trim step 6's extras on phones and
   merge Explain §3–§4 on portrait.
3. `phone-tall__tour-step6.png`: the "Orbit size" view is about 60 px high. Consider `hidePortrait` for it in step 6.

**What works**
- Linear vs exact path-line modes, with the measured and formula drift side by side.
- The dyed line leaning forward.
- Invariant rows that exercise the real JS: focal constant, ellipse = 1, clockwise sense, drift ∝ a².

## capillary_gravity_waves
**Round-1 resolution.** Must fix 1 (rule 9, "17.8 cm/s") is **resolved**:
- Tour step 5 now reads "its own minimum, 17.76 cm/s at 4.354 cm" (`desktop__tour-step5.png`, `phone-tall__tour-step5.png`),
  built from `cgmin(0.07274, 998.2)`.
- Check Q4 uses the same `c4(WG.cg_min)`.
- The c(λ) legend card prints "c_min 23.12 cm/s at λ_m = 1.712 cm · c_g,min 17.76 cm/s at 4.354 cm".
- Explain §6 prints 4 significant figures.

Round-1 Should fixes:
- 1 (hard-coded 23.1 / 1.71 / 19.4) is done: tour steps 4 and 7, check Q1 and the D16 check all compute from `cmin()`
  with 4 s.f. (23.12 → 19.40 cm/s, 1.712 → 1.205 cm).
- 2 is done: the rose "√(2πσ/ρλ)" label is now clear of the orange dot label.
- 3 is done: the water constants are unified with the packets explainer.
- No regression: lint ok, 8-size shot PASS, 18/18 parity.

**Must fix** — none.

**Should fix**
1. The live label of the reader's wave uses 3 s.f. At λ = λ_m (the Explore default and the "λ = λ_m" preset) the orange
   dot, the status ("c = 23.1 cm/s") and the Explain §3 box print the book's rounded (7.58) value. It is computed, so
   it is not a rule-9 breach. Printing c with 4 s.f. when |λ − λ_m|/λ_m < 1 % would make every c_min display read 23.12.
2. `phone-tall__tour-step5.png`: the orange "c = 28 cm/s" dot label covers the start of the blue "gravity only (σ = 0)"
   branch label (only "= 0" shows). On narrow views, drop the branch label, or move it right, when the dot is within
   ~80 px.

**What works**
- The two restoring pushes are drawn at the crest with their pressures.
- The term bars really add up to c²/tanh kH, with a c²_min floor.
- Liquids are offered as modes.
- D16 is derived in front of the valley.
- The "scan never below c_min" invariant.
- The new "our computed clean-water numbers" constants (`WM`, `WG`, `HM`) feed every prose string. This is a good
  pattern for `viz_patterns.md`: never type a number that a function can print.

## seiche_standing_waves
**Must fix** — none.

**Should fix**
1. `phone__explore.png`: the Explore intro says "click the period plot", but that view is hidden on portrait phones.
   Use a phone-specific intro.
2. On phones the η/u view is about 70 px high. Give it more height in steps 3–4.

**What works**
- Right-going, left-going and combined components drawn as ghosts.
- D18 derived with the modes view.
- The period-vs-mode plot with exact, shallow and deep curves.
- Parity rows prove right + left = standing.

## group_velocity_packets
**Round-1 resolution.** Must fix 1 (rule 9, "17.8 cm/s") is **resolved**:
- Water constants are now σ = 0.07274 N/m and ρ = 998.2 kg/m³, the same as the capillary explainer.
- `spd()` prints 4 significant figures, so the pond status reads "c_g,min 17.76 cm/s" (`phone-tall__tour-step7.png`),
  and the c_g(k) label says "c_g,min = 17.76 cm/s".
- Check Q4 uses `spd(CGM.cg_min)` at `lam4(CGM.lam)` (17.76 cm/s at 4.354 cm).
- The source comment now says "(computed)".

Round-1 Should fixes: all four are done.
- 1: the capillary energy *why* gives "½σ⟨η_x²⟩ = ¼σk²a², doubled", and there is a new invariant row.
- 2: the beats corner labels are gone and c, c_g sit in the view title (`phone-tall__tour-step3.png` is clean).
- 3: the Explain §3 time/width line is split.
- 4: the new tour step 7 "A stone in a pond" is photographed at every size. It shows the calm disc r = 33.1 cm at
  t = 1.86 s and the flattest tangent on ω(k).

No regression: lint ok, 8-size shot PASS, 26/26 parity.

**Must fix** — none.

**Should fix**
1. Code tab, pond live value `lam0: '0.0435'` (L3156) is hard-coded. Compute it from `CGM.lam` (0.04354), like every
   other number.
2. `phone-tall__tour-step7.png`: the pond x-axis on phones has no unit (0 / 50 / 100). Add "cm" to the tick labels or
   the view title.

**What works**
- The chord vs the tangent on ω(k).
- The dashed "group moving at c" ghost.
- The energy strip with the conserved ∫E dx.
- The printed-slip toggle with its own parity row.
- The pond step: a live "Right now" callout ("calm out to r = 33.1 cm; the ripples at that edge have λ ≈ 4.354 cm")
  tied to the purple minimum on the c_g(k) view.

## wave_rays_refraction
**Round-1 resolution.** Must fix 1 (the island step highlighted beach code and showed a "…" placeholder) is **resolved**:
- A new Code block `island` ("A ray round the island") defines the smoothed-cone depth
  `H0*np.tanh(s*(np.hypot(X[0], X[1]) - R)/H0)`. This is the same geometry as the JS `geomOf` (source header L2361).
- It launches ray 3 of 9 through `ch07.ray_trace(..., H_fn=H, H_min=1.0)`.
- Tour step 5 lights its lines 11–13 (trace + endpoint).
- `codeRay()` traces the chosen geometry independently of the displayed mode, so every live comment holds a real number
  (`desktop__tour-step5.png`: `[0.0, 569.0]`, `# ~(−3.22, −450) m: the lee side`).
- Three island parity rows against `ray_trace` prove the code.

Round-1 Should fix 1 is done: the status reads "2 of the 5 rays that hit the island land on its lee (sheltered) side",
and "2 rays reach the lee side" on phones. No regression: lint ok, 8-size shot PASS, 23/23 parity.

**Must fix** — none.

**Should fix**
1. The code ray ends at (−3.22, −450) m, which is the island's southern flank, only 3 m past x = 0. The `where` label
   calls this "the lee side" because the test is `xe < 0`. The status `ge.lee` uses the same criterion. Use a stricter
   test (for example, x < −R/2 or the polar angle > 120° from the incoming direction) or say "the flank / lee side", so
   that a ray grazing the side is not counted as sheltered.
2. `phone-tall__tour-step5.png`: the "15 m / 10 m / 5 m" contour labels overprint one another on the phone plan view.
   On views narrower than 420 px, keep only the outermost label.
3. `phone-tall__tour-step5.png`: the along-ray view's rotated y title is clipped ("÷ start valu"). Shorten it to "÷ start"
   on phones.
4. Round-1 item 2 is still open: in island mode the unselected rays remain faint behind the island.

**What works**
- Hamilton's ray equations integrated in the page, with ω drift in the status.
- RK4 vs `ray_trace` parity (beach and now island).
- Snell closed form vs RK4.
- The x–t diagram for D22.
- The inspector's Snell arithmetic.
- The island Code block with a live "where did it land" comment, a good pattern (one code block per geometry, traced
  off-screen so live values never depend on the displayed mode).

## hydraulic_jump
**Must fix** — none.

**Should fix**
1. Explore needs 12 pages and Explain 16 pages at 360×640. Tour step 2 "The box" is now 4 pages at 390×844 (it was 9).
2. The "heights ×2.21" tag has moved clear of the face-2 amber arrow (`phone-tall__tour-step2.png`). Done.
3. **Intermittent audit failure.** One `--quick` run failed:
   `[phone-tall] tour-step2: viz-step-card page 2/4: viz-terms … clipped 3px (cut through a line)`. Three re-runs (one
   full, two quick) passed. The term list sits within a few pixels of the page height when KaTeX finishes late. Give
   step 2's term list 1 line of slack on phones (for example, drop its derivation quote there, since step 3 links it)
   so that the pager never packs to the edge.

**What works**
- The momentum and energy budgets as term bars.
- The forbidden preset.
- Bore and solitary modes with KdV invariants.
- The D26 sign step.

## two_layer_modes
**Must fix** — none.

**Should fix**
1. Tour step 3: write the whole of (7.110), not only its first factor.
2. Phone Explain is 19 pages and Explore 10. Merge §5–§6 on portrait.
3. The (7.96) note is opaque. Give the E_k = E_p = ¼(ρ₂ − ρ₁)ga² sentence instead.

**What works**
- Separate magnifications for the surface and the interface.
- The vortex-sheet halo.
- ω(k) of both roots with their limits.
- The "which g′" note.
- D30 verified symbolically.

## internal_wave_beams
**Round-1 resolution.** Must fix 1 (the key views were unreadable on portrait phones) is **resolved**:
- On portrait phones the tank and the K plane now sit side by side, each at full stage height, with the ω/N view
  hidden.
- The K plane uses a square data window ±1.1|K| (ticks −2…2 around the K = 2.5 circle).
- K, **c**, **c_g**, the right-angle mark and the constant-ω rays are legible at 360×640 (`phone__explore.png`,
  `phone__tour-step1.png`), at 390×844 in the step about c ⟂ c_g (`phone-tall__tour-step5.png`), and at 844×390
  (`phone-land__explore.png`).
- The tank is cropped to |x|, |z| ≤ 5 m, so the cross and its green particles fill the view.

Round-1 Should fixes:
- 1 is done: the preset is now "K up-left (k < 0)", and the book figure number survives only as a parenthetical in the
  D36 *why*.
- 2 is done: the D33 titles are renamed.

No regression: lint ok, 8-size shot PASS (432 views), 27/27 parity. The desktop layout is unchanged
(`desktop__tour-step5.png`).

**Must fix** — none.

**Should fix**
1. Phone K plane: the only ray label shown is "0.25" (the ω/N = 0.25 ray), bottom-left with no "ω/N =". It reads as a
   stray number. Drop it on views narrower than 420 px, or write "ω/N 0.25".
2. Phone K plane has no k-axis title (the m axis has one). Add "k" under the ticks or to the view title ("K plane · k, m
   [rad/m]").
3. `phone-tall__tour-step5.png` page 1 of 3 ends with the top 10 px of the orange derivation-quote box. Start the quote
   on page 2 (a page break before the box).

**What works**
- The sign-safe treatment of k < 0.
- The ω/N-vs-θ view.
- The live w-equation residual.
- E_k = E_p and F = c_g E invariants.
- On phones: two square views side by side instead of two squashed strips. This is a good pattern for any "physical
  space + wavenumber space" pair.

---

## Verdict: PASS

All four round-1 Must-fixes are resolved, and there are no new Must-fixes. Every explainer PASSES.
- **capillary_gravity_waves**: no exercise answer remains. c_g,min is computed and printed as 17.76 cm/s at 4.354 cm,
  and c_min as 23.12 cm/s at 1.712 cm (4 s.f.).
- **group_velocity_packets**: the same computed 17.76 cm/s with unified water constants. A pond step was added.
- **wave_rays_refraction**: the island step has its own island code block with real live values and 3 island parity
  rows. There is no placeholder.
- **internal_wave_beams**: the portrait phone views are legible (tank and square K plane side by side).

Across the 9 explainers: lint 9/9 ok; full shot on the four fixed ones 4/4 PASS; quick shot 9/9 PASS after a re-run
(`hydraulic_jump` had one intermittent 3 px clip, see its Should fix 3); parity 222/222 rows ok (172 JS ↔ fluidpy).
The remaining items are Should fixes.
