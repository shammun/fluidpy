# Chapter 1 — explainer review (round 2)                          2026-09-13

Reviewer: viz-reviewer, round 2 (re-check after the builders' fixes). Everything below was re-run today:
- `tools/viz_lint.py --chapter ch01`.
- `tools/shot.py --chapter ch01` (full run: 8 sizes, every tab, step and derivation page).
- The Playwright click/drag/wheel/tap pass over every visible view canvas: 4 tabs × 4 sizes (desktop 1366×768, notebook
  1000×700, phone-tall 390×844 with touch, phone-land 844×390 with touch), 2 141 actions in total.
- 21 targeted interaction checks, each confirming that the documented gesture changes `app.state`.
- A browser check of the E4 status badge against `ch01.lapse_rate_stability`: 12 lapse rates × 2 conventions × 2
  layouts.
- Derivation step counts read from `app.cfg.derivations` and compared with `notebooks/build_ch01.py` and
  `outputs/ch01/executed.ipynb`.
- Extra screenshots in `reports/viz/ch01/_round2/`: E1 at 844×390 on walkthrough step 6; E4 in ocean mode, isothermal
  and meteorology-unstable at 360×640 and 390×844; E3 in stirring mode.

## Summary

| explainer | lint | shot (sizes · views · steps) | parity rows ok | clicks (page errors) | explain | derivations (steps = notebook) | teaching | verdict |
|---|---|---|---|---|---|---|---|---|
| continuum_averaging_volume (E1) | ok | PASS · 8 · 152 · 7 | 14/14 | 449 actions, 0 errors; resample ✓, dot probe ✓ | 7 sections + interpretation, 18 live | D35 5 = 5 ✓ | good; Kn strip now on landscape phones | **PASS** |
| viscosity_momentum_diffusion (E2) | ok | PASS · 8 · 96 · 7 | 13/13 | 345, 0; plate drag ✓, profile probe ✓ | 7 + interpretation, 9 live | none assigned | good (unchanged) | **PASS** |
| heat_work_paths (E3) | ok | PASS · 8 · 272 · 7 | 13/13 | 449, 0; gas scrub ✓, p–v / T–s probe ✓ (now via `stage.onPointer`) | 8 + interpretation, 7 live | D10 7 = 7 ✓ · D14 11 = 11 ✓ | very good | **PASS** |
| parcel_stability (E4) | ok | PASS · 8 · 512 · 8 | 49/49 | 449, 0; column drag ✓, profile probe ✓, ζ scrub ✓ | 8 + interpretation, 6 live | D18 14 = 14 ✓ · D19 16 ✓ · D20 5 ✓ · D36 8 ✓ | very good; both conventions always on screen | **PASS** |
| buckingham_pi_machine (E5) | ok | PASS · 8 · 224 · 8 | 26/26 | 449, 0 (KaTeX ³ warnings gone); card/header/cell/group ✓ | 7 + interpretation | D28 13 = 13 ✓ | very good | **PASS** |

**Checks that apply to all five explainers:**
- KaTeX renders at every size, with no fallback text.
- The console is clean: no errors and no warnings.
- No page overflows at any of the four click-test sizes.
- The contract and depth features are unchanged from round 1.
- **No regressions found.**

**E4 badge vs fluidpy: 48 of 48 cases match.** In every case the status reads `<emoji> <verdict> · <exact fluidpy text in the chosen convention> · <other convention>: <bare fluidpy text>`:
- **Lapse rates tested:** −12, −9.8, −9.76, Γa, −6.5, 0, +5, +10, −9.25, −0.01, −20 and +3.3 K/km.
- **Where:** Kundu and meteorology conventions, at 1366×768 and 390×844.
- **Verdict:** the verdict word equals `lapse_rate_stability(...)[0]`.
- **Chosen-convention text:** the substring equals `lapse_rate_stability(dTdz/1000, convention=conv)[1]` byte for byte.
- **Other-convention text:** the tail equals the other convention's text with `prefix=False`.

The new selftest row "status = verdict · fluidpy text · other convention (both)" pins this inside the page, next to 13
exact-text parity rows.

**Derivation alignment:**

| D | notebook | explainer | move titles |
|---|---|---|---|
| D35 | 5 | 5 | ✓ |
| D10 | 7 | 7 | ✓ |
| D14 | 11 | 11 | ✓ |
| D18 | 14 | 14 | ✓ |
| D19 | 16 | 16 | ✓ |
| D20 | 5 | 5 | ✓ |
| D36 | 8 | 8 | ✓ |
| D28 | 13 | 13 | ✓ |

The move titles match one to one, including the new D14 steps 3–4, D18 steps 13–14 and D28 step 9. The
walkthrough deep links point at the renumbered steps: E3 → D14 step 10 "Integrate", E5 → D28 step 12, and E4 → D18
step 11 "Name N²".

Algebra re-checked on paper:
- **D14 step 5:** C_p dT / C_v dT = v dp / (−p dv).
- **D18 step 4:** dividing (−ρ_pVg + ρ(z_o+ζ)Vg) by ρ_pV gives −g(ρ_p − ρ(z_o+ζ))/ρ_p.
- **D18 step 14:** the double root gives ζ = A + Bt, and release from rest sets B = 0.
- **D28 step 9, live system for Π₁ = Δp U^a d^b ρ^c:** M: 1 + c = 0, L: −1 + a + b − 3c = 0, T: −2 − a = 0, so
  a = −2, b = 0, c = −1, which gives Δp/(U²ρ).
- **D20's new sentence:** θ at z = 0 is T(1000/1013.25)^0.2857 = 0.9962 T, about 1.08 K below T. The sentence is true.

---

## continuum_averaging_volume (E1) — PASS

**Round-1 items**

| # | item | status |
|---|---|---|
| Must 1 | Kn strip invisible on phones although the text said to turn the phone sideways | **fixed.** At 844×390 the strip is its own row under the graph, showing the continuum/slip/transition bands with "body 0.064 · box 0.064" (`_round2/…phone-land_tour6.png`, `phone-land__explore.png`). Upright phones get the curve title "Kn(body) 0.064 → slip flow" and a Kn (body) readout on step 6 (`phone-tall__tour-step6.png`). The Explain hint matches the layouts |
| Should 1 | step 4 showed two different boxes | **fixed.** One box, L = 21 nm: N̄ = 271, Eq. card 1.47 kg/m³, inspector N = 274 → 1.49 kg/m³ (`desktop__tour-step4.png`) |
| Should 2 | "zero, one or two molecules" | **fixed.** Now "zero, one or a few" |
| Should 3 | 360×640 label clipping | **fixed.** The portrait scene is a compact picture with a text column; the y label is shortened to "÷ ρ(x₀)" (`phone__explore.png`) |
| Should 4 | undefined names in the code snippet | **fixed.** `eps, L_flow, L_body` and `T, p` lines added, with live values |
| Should 5 | check 3 "2.2 decades" | **fixed.** Now "≈ 2.25 decades" |

**Must fix** — none.

**Should fix**
1. **The status noise and the inspector noise disagree on the same screen.** The status badge uses the *mean* count,
   `noiseExpected(L, n)` with n from ρ. The sample text and inspector use the count *at the crest*,
   λ = nL³(1 + ε·sinc), which is 1.2× larger at ε = 0.2.
   - **Step 4** (`desktop__tour-step4.png`): the badge reads "molecular noise ±6.7 %" (1/√226) while the inspector shows
     "expected scatter 0.0607" (1/√271).
   - **Step 1** (`phone-tall__tour-step1.png`): "±111 %" sits next to "expected 0.967", which implies ±102 %.

   **Fix:** use the crest count λ in `model().noise`, or label the badge "at the mean density". This matters on the one
   step whose purpose is tracing the arithmetic.

---

## viscosity_momentum_diffusion (E2) — PASS

The file is unchanged since commit b45cbef. The audit and click results are identical to round 1. The round-1 Should
items remain open and none of them blocks:
- the phone profile view is small;
- the stress y-label shows "τ_w" with an underscore;
- the Eq. card value lags the canvas while playing on step 3;
- the heat-mode q label sits on the AB line.

**Must fix** — none.

---

## heat_work_paths (E3) — PASS

**Round-1 items**

| # | item | status |
|---|---|---|
| Should 1 | D14 had 10 steps; the notebook has 11 | **fixed.** Step 3 is split into "Set ds = 0 in the second Gibbs form" and "Use h = h(T)". The references are renumbered: "Divide step 4 by step 2", "Substitute into step 7". The walkthrough links to step 10 (`desktop__derive-d2p4.png`) |
| Should 2 | D10 step 7 *why* 46 words | **fixed.** Now 25 words, and the replacement argument moved to *in words* |
| Should 3 | D10 step 2 live line is perfect-gas only | **fixed.** The `watch` text now says "perfect gas with constant C_v" |
| Should 4 | phone p–v legend overlapped the route | **fixed.** The legend moved into the view title on narrow views (`phone-tall__derive-d2p10.png`) |
| Should 5 | local capture-phase pointer workaround | **fixed.** `installPointer` is removed and gas scrub and probe run through `stage.onPointer`. Targeted checks pass: drag scrubs `prog`, p–v and T–s clicks set `probe`, 0 page errors |
| Should 6 | moving dot on the "no equilibrium path" in stirring mode | **fixed.** No dot is drawn in the p–v view (`_round2/heat_work_paths__desktop_stir.png`) |

**Must fix** — none.

**Should fix**
1. **The end-of-run card is clipped on portrait phones.** In `phone-tall__derive-d2p10.png` the card "route done: q = 0,
   w = 68.8, Δe = 68.8 kJ/kg · Δs = 0 J/(kg K)" runs past the right edge of the p–v view: the unit is cut to "J/(kg l".
   It also covers the "2" state label.
   - **Fix:** on views narrower than about 420 px, break the card into two lines, or drop the Δs part. The status line
     already carries Δs.

---

## parcel_stability (E4) — PASS

**Round-1 items**

| # | item | status |
|---|---|---|
| Must 1 | meteorology Γ not always on screen | **fixed.** The badge carries both conventions at every size and on every step. For example, "🌊 stable · stable ⇔ dT/dz > Γa: −6.5 > −9.8 K/km · met: 6.5 < 9.8 K/km" (`phone-tall__tour-step1.png`, `phone__tour-step1.png`). The same holds in the meteorology convention: "… · Kundu: −12.0 < −9.8 K/km" (`_round2/…phone_met_unstable.png`). The badge wraps to two lines on phones without truncation. Verified in 48/48 browser cases against fluidpy |
| Must 2 | unstable badge led with "stable" | **fixed.** "🚀 unstable · stable ⇔ dT/dz > Γa: −12.0 < −9.8 K/km · met: 12.0 > 9.8 K/km" (`phone-tall__tour-step6.png`). Neutral reads "⚖️ neutral · … −9.8 = −9.8 K/km" (`desktop__derive-d1p14.png`) |
| Must 3 | wrong tick glyphs on portrait phones | **fixed.** The narrow-view y-title strip is separate from the ticks. The atmosphere shows "1500" (`phone-tall__tour-step1.png`), "2000 / 0" (`phone-tall__tour-step6.png`), and the ocean "−150" with its minus (`phone__explore.png`, `_round2/…phone-tall_ocean_explore.png`) |
| Should 1 | D18 13 vs 14 steps | **fixed.** Step 13 covers release from rest for N² ≠ 0 (thermocline, period 10.7 min). Step 14 treats the double root, moving the picture to Γ = Γa with a flat ζ(t) (`desktop__derive-d1p13.png`, `d1p14.png`) |
| Should 2 | step 1 WATCH refers to the hidden column | **fixed.** It now reads "(on phones, the coloured bracket in T(z))" |
| Should 3 | "cools more slowly" wording for dT/dz ≥ 0 | **fixed.** A separate branch now says "does not cool with height at all" |
| Should 4 | ΔT/Δρ labels drawn across the curve | **fixed.** The labels sit on solid cards (`phone-tall__tour-step6.png`, `desktop__derive-d1p13.png`) |
| Should 5 | D20 p_o = sea level vs 1000 hPa | **fixed.** The sentence is in the D20 goal and in the Equations `p_o` symbol row, and its number is correct (≈ 1.08 K) |
| Should 6 | "switch the medium" hint on phones | **fixed.** It now says "Use the thermocline or salt tank preset (or the Medium chips on larger screens)" |

**Must fix** — none.

**Should fix**
1. **Raw TeX shows in D18 step 13 *why*.** The plain-text *why* prints "e^{±λt}" literally (`desktop__derive-d1p13.png`).
   Write it as "e^(±λt)" or wrap it in `$…$`.
2. **Tick labels collide at 360×640 in the unstable case.** In Explore at 360×640 with the −12 K/km preset, the ζ(t) view
   is about 45 px tall and its "1000" / "−1000" tick labels overlap the "e-fold" marker (`_round2/…phone_met_unstable.png`).
   The two-line badge takes some of the height. The walkthrough and the stable cases read fine.
   - **Fix:** when the plot is under 60 px, draw only the 0 tick, or put the ±ζ_max value in the view title.

---

## buckingham_pi_machine (E5) — PASS

**Round-1 items**

| # | item | status |
|---|---|---|
| Should 1 | D28 12 vs 13 steps | **fixed.** New step 9, "Build one group for each other variable". Its live exponent system for Π₁ checks out (a = −2, b = 0, c = −1). Old steps 9–12 are renumbered 10–13 and the walkthrough links to step 12 (`desktop__derive-d1p9.png`) |
| Should 2 | SI / cgs / imperial columns lost at 1000×700 | **fixed.** All three columns show at notebook size, with the footer "columns: SI · cgs · imperial" (`notebook__equations.png`) |
| Should 3 | "test every 3×3 block" wording | **fixed.** It now reads "the first two 3×3 blocks have determinant 0; the third, (Δp, Δx, U), gives −1 … 21 are nonzero". The first two blocks, (Δp, Δx, d) and (Δp, Δx, ε), are singular and the third has det −1, which I confirmed |
| Should 4 | KaTeX "³" warnings | **fixed.** Units are mapped to TeX powers and the console is clean at all 8 sizes |
| Should 5 | clicking the starred Δp header adds it to a custom set | open (optional). A header click still gives `custom: "U,d,rho,dp"`. The status explains the error, so this does not block |

**Must fix** — none.

**Should fix**
1. **D28 step 9 *why* is too long.** At 53 words it exceeds the 35-word limit (`teaching-style` §5); round 1 proposed that
   text. Keep "A_rep a_j = −A_·j is square; det A_rep ≠ 0 (step 8) makes a_j unique." Move "each group holds a different
   q_j to the first power, so the n − r groups are independent" into *in words* or *watch*.

---

## Candidates for `knowledge/viz_patterns.md` (new this round)
- A status badge that holds two conventions: verdict word first, then the library's exact text, then the other
  convention's bare relation. Pin it with a selftest row that rebuilds the badge from `cfg.status` (E4).
- On narrow views, a separate y-title strip keeps rotated axis titles from overlapping tick glyphs (E4).
- When a view is hidden on portrait phones, put its key number in the title of a visible view, e.g. "Kn(body) 0.064 →
  slip flow" (E1).
- A step that traces one sample's arithmetic must set the global parameter to that sample's own value (E1 step 4).

## Verdict: PASS
All five explainers pass lint and the full shot audit at 8 sizes, and all 115 parity rows are ok. The 2 141 scripted
interactions produced 0 page errors and 0 console warnings. Every round-1 Must fix is fixed and confirmed by eye: 1 in
E1 and 3 in E4. Every derivation now has the same number of steps as the notebook, with matching moves. The E4 badge
matches `ch01.lapse_rate_stability` text in both conventions. The remaining Should-fix items do not block publishing:
- E1 noise basis;
- E3 clipped end card on phones;
- E4 raw "e^{±λt}" and the 360 px ζ ticks;
- E5 step-9 *why* length;
- E2's round-1 cosmetics.
