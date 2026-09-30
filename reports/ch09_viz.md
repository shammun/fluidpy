# Chapter 9 — explainer review                                  2026-09-30

Reviewer: viz-reviewer (phase 7 gate). Re-ran `tools/viz_lint.py --chapter ch09` (9/9 ok) and a **full** `tools/shot.py --chapter ch09`
(8 sizes, every tab, every walkthrough step, every derivation page; 8 PASS, 1 FAIL — wall_jet_invariant at tablet) after the fluidpy review fixes. Screenshots and audits:
`reports/viz/ch09/<slug>/`. Also ran `tools/eq_refs.py --chapter ch09` (13 hits, judged below), a control-character /
lost-backslash scan of every chapter script, rendered pp. 365 (9.7–9.11), 403–407 (9.62–9.85) and cross-checked every
number quoted in text against `fluidpy.ch09` (φ_sep 103.11° / 100.90°, x_sep/L 0.15829 / 0.12569, λ_FS = −0.068148,
C_D,p 0.8838 / 0.5778, street roots 0.2259 / 0.3133 / 0.3428, m_sep −0.090429).

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| bl_scaling_thicknesses | ok | PASS 8/8, 0 failures | 30/30 (45/45 rows) | 8 sections, 5 boxes, live | D01 14/14 · D03 6/6 · D04 9/9 | strong | good | **PASS** |
| blasius_similarity_collapse | ok | PASS 8/8, 0 failures | 24/24 (35/35) | 8 sections, table, live | D05 13/13 · D06 10/10 | strong | good | **PASS** |
| falkner_skan_family | ok | PASS 8/8, 0 failures | 18/18 (22/22) | 8 sections, preset table | D07 12/12 · D12 6/6 (step 1 title garbled) | strong | TeX escapes broken | **FAIL** |
| thwaites_marching | ok | PASS 8/8, 0 failures | 28/28 (36/36) | 8 sections, accuracy table | D08 11/11 · D09 10/10 · D10 9/9 · D11 8/8 | strong | good | **FAIL** (Example 9.2 criterion) |
| cylinder_drag_crisis | ok | PASS 8/8, 0 failures | 25/25 (28/28) | 10 sections, regime table | D13 7/7 | strong | Explain §2 hint garbled | **FAIL** |
| karman_street_stability | ok | PASS 8/8, 0 failures | 9/9 (13/13) | 10 sections, growth table | D14 15/15 (★★★, algebra checked) | strong | minor | **PASS** |
| free_jet_similarity | ok | PASS 8/8, 0 failures | 17/17 (27/27) | 9 sections, cases table | D15 7/7 · D16 14/14 · D17 11/11 · D18 5/5 | strong | minor overlaps | **PASS** |
| wall_jet_invariant | ok | **FAIL** at tablet (768×1024): 2 Explain boxes 1–2 px too wide; 7 other sizes ok | 19/19 (29/29) | 9 sections, 7 boxes, live | D19 13/13 · D20 12/12 · D21 15/15 (result line cites bare (9.83)) | strong | minor | **FAIL** (shot + bare (9.83)) |
| teacup_secondary_flow | ok | PASS 8/8, 0 failures | 5/5 (9/9) | 7 sections, live | D22 7/7 | good | phone cup cluttered | **PASS** |

Physics spot checks done by hand (all correct): (9.8) scaled y-momentum and the sign of ∂p*/∂y* in D01 step 12 match p. 365 and
the corrected `bl_dpdy_scaled`; D05/D07/D16/D20 cancellations (the ηf′f″ terms) and the wall-jet coefficient 4 (D20: 4f‴ + ff″ + 2f′² = 0,
re-derived); D13 closed form sin φ_s(1 − 4/3 sin²φ_s − C_b); D14 steps 8–11 (bracket π²[γd̄_A − iσd̄_B], the 4 × 4 matrix rows);
free-jet (9.62)–(9.76) and Re_x against p. 403–404; wall-jet (9.80)–(9.85) against p. 405–407. Parity rtols are honest (1e-12 for
closed forms, 1e-6…1e-9 for ODE/table routes; the loosest are table interpolations, 1e-3 for an inflection height).

Consistency with the fluidpy review fixes: (9.8) sign — consistent (E1). R5 "trap, not slip" — no explainer calls (9.33) a slip; E2/E5 say
"one side of the plate" (consistent). Falkner–Skan η_max = 16 for m ≤ −0.05 — E3 parity rows pass (`falkner_skan_state`, reversed branch).
LAMBDA_SEP_FS — E4 takes −0.068148 from the embedded closure table and its parity row `thwaites_cylinder_separation(ch09.LAMBDA_SEP_FS)`
passes; **but the Example 9.2 preset reports the exact-FS criterion first** (see E4 Must fix). `separated_pressure_drag` cp_base default —
every E5 parity row passes C_b explicitly (82/−1.2 → 0.8838, 125/−0.6 → 0.5778). `magnus_sign` — not used by any built explainer (B1 not built).

---

## bl_scaling_thicknesses
**Must fix** — none.

**Should fix**
1. Explain §5 (`explainHtml`, line ~2749): `2\,tau_0\,x` is missing a backslash; it renders "2 tau₀ x" in italic letters. Write `2\,\tau_0\,x`.
2. Amber (δ*) and orange (θ) are hard to tell apart in the u/U view (`desktop__tour-step2.png`, `phone-land__explore.png`): both dashed
   height lines and both shaded areas read as orange. Use a yellow-amber (e.g. `#d9a300`) for δ* or hatch one of the areas.
3. Phone title of the u/U view is clipped ("θ 2…", `phone-tall__tour-step3.png`); drop "δ*" from the title on phones (it is in the readouts).

**What works** — the plate view with real vs ideal streamlines lifted by exactly δ* (a static figure cannot show this); the shape chips with a
fixed δ₉₉ make "same edge, different δ* and θ" click; D01 carries the book's y-momentum (9.8) with the correct signs and a live 1/Re; Explain §5
closes D04 numerically (ρU²θ = 2τ₀x).

## blasius_similarity_collapse
**Must fix** — none.

**Should fix**
1. In rescaled mode the three stations are drawn as the same curve `fpv(η)` three times; plot each station's own u(y) samples divided by its
   own δ(x_k) (as dots at different y) so the collapse is shown from data, not assumed.
2. Phone walkthrough step 2 slices the derivation quote card so that "All steps →" is cut by "continued on the next page"
   (`phone-tall__tour-step2.png`); shorten the step text or drop the quote on phones.

**What works** — raw ↔ rescaled mode is the aha; the shooting overlay with "guess too small / too large / converged" and the truncation
preset make the boundary-value problem tangible; the τ₀(x) panel with "area = F_D" ties D06 step 9 to the picture.

## falkner_skan_family
**Must fix**
1. **Broken TeX in D12 step 1** (`desktop__derive-d2p1.png`: the move reads "Evaluate $uu_x+vu_y=-♠rac1 ho♠rac{dp}{dx}+ u u_{yy}$ (9.9) at the
   wall"). Line 2816: the `did` is a plain single-quoted JS string with single backslashes (`\f` = form feed, `\r` = carriage return, `\n`
   in `\nu` = newline). Double every backslash (`-\\frac1\\rho\\frac{dp}{dx}+\\nu u_{yy}`) or use `String.raw`.
2. Same class of bug in Explain §6 (`presetTable`, line 2715): `'$\eta_i = '` renders "eta_i". Write `'$\\eta_i = '`. Run a scan for
   single-backslash TeX in plain `'…'` strings before re-auditing.

**Should fix**
1. Parity row `m_sep (fold)` uses rtol 1e-4 for a number copied from the Python table (`FS_TAB.fold`); it agrees to ~1e-12 — tighten to 1e-9
   so the row can catch a stale table.
2. At 360×640 the two views are ~80 px tall and the profile's y title is clipped (`phone__explore.png`); hide the `curve` view on the
   smallest portrait size (its f″(0) is in the status and readouts) or collapse the preset chips to one row.
3. Explain §1 box writes "−dp/dx = U_eU_e′ (9.35)" (per unit density, as printed) right under the dimensional line; add "per unit density"
   inside the box to avoid a units jump.

**What works** — one dial n with the fold, the dashed second branch and the "reversed" preset; the rose/orange term bars that stay equal
(μu_yy = dp/dx) for every n; the independent RK4 shooting check printed in Explain §2; the wall-curvature parabola on click.

## thwaites_marching
**Must fix**
1. **Example 9.2 preset reports the non-book criterion first.** The preset "diffuser (Example 9.2)" keeps the default closure `fs`, so the
   status, readout and Explain §5 headline give x_sep/L = 0.126 at λ = −0.0681 (`desktop__tour-step5.png`), while the book's example (and
   now `ch09.example_9_2()["criteria"]`, "the book's first") uses λ = −0.090, x/L = 0.158. The −0.09 line is labelled only "Thwaites fit"
   and never "the book's criterion". Change: (a) the Example 9.2 preset sets `closure: 'white'` (or the status for the diffuser shows
   "book criterion λ = −0.090: x/L = 0.158 · exact Falkner–Skan λ = −0.0681: x/L = 0.126", book first); (b) relabel the −0.09 line/chip
   "−0.090 book (Thwaites' table)"; (c) the cylinder regime text and D11 already give both angles — list 103.1° (book criterion) first.

**Should fix**
1. L(λ) view: the label "fold n = −0.0904" is right-aligned against λ = −0.0681 and therefore sits on the λ = −0.09 dotted line
   (`desktop__derive-d3p8.png`) — exactly the n ≈ λ ≈ −0.09 coincidence a reader will confuse. Put it right of the red line and write
   "fold n = −0.0904 (λ = −0.0681)".
2. D10 step 6 says 64 values of n; the design says 60 rows — the text is right for the embedded table (64 rows); update the design note or
   keep the 64 but say "64 values from the fold to n = 50".

**What works** — the teal/rose θ² bar (integral vs inlet memory) and the inspector that re-does (9.50) → λ → l → τ₀ at the cursor; the
diffuser closed form checked against the march in a selftest row; the cylinder body view with φ ticks; the live "how good is it" table.

## cylinder_drag_crisis
**Must fix**
1. **Explain §2 hint is garbled** (line 2809): the source contains a literal form feed and a tab (`$St=\x0crac{…}` and `2\pi\x09imes0.2`),
   so KaTeX gets "rac{Ωd}{U∞}" and "imes". Replace with `\frac` and `\times` (the string is `String.raw`, so single backslashes are right —
   the control characters were introduced when the file was written).

**Should fix**
1. Explain §0 is not mode-aware: in sphere mode it still describes "C_p(φ): pressure around the cylinder" while the view shows the low-Re
   sphere drag laws (`desktop__explain.png`). Branch §0 on `n.sph`.
2. With U fixed at 10 m/s the street preset gives d = 0.15 mm and f = 13 333 Hz (`phone-tall__tour-step2.png`) — correct but odd for a
   first reading; default U = 1 m/s (d = 1.5 mm, f = 133 Hz) makes the Aeolian tone audible and matches the design's "2 mm wire" idea.
3. Notebook frame: the "φ_s = 125°" label overlaps the U arrow in the flow view (`notebook__equations.png`); offset the label below the dot.

**What works** — the log-Re dial that moves cartoon, C_p(φ), C_D and the regime table together; the "set C_b by hand" check that later
separation alone does not lower drag (a genuinely non-obvious lesson); honest "illustrative / rounded" labelling everywhere; the rose
"missing recovery" area.

## karman_street_stability
**Must fix** — none.

**Should fix**
1. Inspector (line 2686): `'lower row, $+\Gamma$'` is a plain string inside the template, so `\G` becomes "G" and the inspector shows
   "+Gamma". Use `'$+\\Gamma$'`.
2. Explain §0 says the orange dots are the maximum "over 31 wavenumbers", §6 says 61 — make them agree with `growthMax` (61).
3. D14 step 1 puts +Γ on the upper row; the picture, Explain §0 and the code put −Γ on the upper row. The why says the sign only reverses
   the travel direction, but a reader comparing colours will stumble — use one convention in both.

**What works** — the ★★★ derivation is complete and correct move by move (checked), with live γ, σ and P, Q; the eigenvalue view whose
dots reach the imaginary axis at 0.2805; kick-and-watch growth; facing-rows preset showing π/4 at every b/a.

## free_jet_similarity
**Must fix** — none.

**Should fix**
1. Laptop/notebook: "inflow v = 59.3 mm/s" collides with "y stretched ×40" in the jet view (`laptop__check.png`); stack them or drop the
   stretch note below 420 px width.
2. Tour step 3: `x^{-1/3}` wraps between base and exponent on desktop (`desktop__tour-step3.png`); wrap the power in `\mathrm{}`-free
   braces `{x^{-1/3}}` or put the pair in one `$…$`.

**What works** — the wrong-exponent toggle that makes ρ∫u²dy visibly fall (the invariant earns its keep); momentum check in the status at
three stations; D18 exposes the printed 5.6152 as the 4 % point with a level slider (a correction the reader can verify by dragging).

## wall_jet_invariant
**Must fix**
1. **`tools/shot.py` FAIL at tablet 768×1024**, quoted: `[tablet] explain: pager page 3/4: viz-work-res 'miss of the closed form: 1.8×10−121−g≈4.' 2px too wide` and `[tablet] explain: pager page 4/4: viz-work-res 'printed: f∞=0.397 (not 1)miss=0.309 (cor' 1px too wide`. The two boxed `aligned` results of Explain §6 and §7 are one-line-too-wide; split each into shorter rows (put the tail formula 1 − g ≈ 4.29 e^{−f∞η/4} and the '(correct: …)' part on their own rows, or use `W.line` instead of a box).
2. **Bare equation number in the D21 result and check.** Result tex (line 2938) reads "… ⇒ (9.83), f″(0) = f∞³/72" and the check text says
   "satisfies (9.83) to 1e-11"; step 12 live says "residual of (9.83)". Write the implicit solution (or a short form such as
   $-\ln(1-g)+\sqrt3\tan^{-1}\frac{2g+1}{\sqrt3}+\ln\sqrt{1+g+g^2}=\frac{f_\infty}4\eta+\frac{\sqrt3\pi}6$) next to the number in the result,
   and "the implicit solution (9.83)" with the formula in the check.

**Should fix**
1. D21 step 15 tex (line 2935) uses `\&` inside `aligned`, so the line renders "f → λf(λη)&f∞ → λf∞" (`desktop__derive-d3p15.png`); use `\\&`.
2. Explain §1: K = (ṁ₁/ρ)²/ν is dimensionally m²/s, but K and C are in m^{3/2}/s; write K = (ṁ₁/ρ)²/(ν x₁^{1/2}) with x₁ = 1 m.
3. Explain §2 calls u₀ = Cx^{−1/2} = 183 m/s the "centre speed" while the peak speed is 14.4 m/s; in wall mode call u₀ the "velocity scale".
4. Profile view: "peak f′ = 0.07875 at η = 8.11" is overprinted by the "f∞ = 1 (dotted)" legend (`desktop__explain.png`); phone bars clip the
   "1.7" label at the top.

**What works** — rose ρ∫u²dy falling 41 % while the teal Ψ bar stays flat to 4×10⁻⁷ is the chapter's best "which quantity is conserved" picture;
the printed-ODE toggle with hollow dots leaving the curves; the gauge preset showing C quartering with an unchanged dimensional profile.

## teacup_secondary_flow
**Must fix** — none.

**Should fix**
1. The blue u(z) curve is drawn on the force axis of the "u(z) and the forces" view without its own scale, and Explain §0 never mentions it;
   add a top axis u/u_e (0 … 1) or a label "blue: u(z)/u_e × ρu_e²/R", and name it in §0.
2. At 360×640 the cup view is ~150 px high and "floor layer, zoomed ×3.5" and "z = 0 mm" overprint the swirl symbols and leaves
   (`phone__explore.png`); hide those labels below ~200 px.
3. The grey dashed curves in the net-force view (other profile shapes?) are unlabelled; label or remove.
4. D22 step 6 invokes plane continuity (6.2) for an axisymmetric loop; say "(in its axisymmetric form (1/R)∂(Ru_R)/∂R + ∂w/∂z = 0)" or keep
   (6.2) but call it the plane analogue.

**What works** — orange never moves with height, rose shrinks toward the floor, teal is the whole effect: the force bars that add
(rose + teal = orange) are exactly the reference-lab idea; river-bend mode generalises the balance; Ekman preview with live scales.

---

`tools/eq_refs.py` hits judged: E3 L2846 and E7 L2375 are code/comments (not reader text); E7 L2643/L2646 and E8 L2701/L2799/L2811 show the
formula in the same line or title (ok); E8 L3007 is an Equations-card ref field and L3043–L3050 are selftest names (not displayed); **E8
L2940 (and the untagged result tex L2938, live L2931) are real bare references → E8 Must fix 2.**

## Round-1 verdict: FAIL
Three explainers need small, mechanical fixes (falkner_skan_family, cylinder_drag_crisis: broken TeX escapes; wall_jet_invariant: a
tablet fit failure in two Explain boxes and a bare (9.83) in the derivation result) and one needs a content fix (thwaites_marching: the Example 9.2 preset must give the book's criterion
λ = −0.090, x/L = 0.158, first). The other five PASS. After the fixes, re-run `tools/viz_lint.py` and `tools/shot.py` on the four files.


---

# Round 2 — re-review of the four fixed explainers                    2026-09-30

Re-ran `tools/viz_lint.py --chapter ch09` (9/9 ok), the control-character / lost-backslash scan (clean apart from one Should-fix in
karman_street_stability), `tools/eq_refs.py` (remaining hits are titles, Equations-card refs and selftest names, all fine) and a **full**
`tools/shot.py` on the four files: **all four PASS at 8 sizes, 0 failures** (falkner_skan_family 272 views / 22 rows,
thwaites_marching 464 / 40, wall_jet_invariant 464 / 29, cylinder_drag_crisis 168 / 28; every parity row ok: 18/18, 31/31, 19/19, 25/25).
I also drove each explainer in a real browser (Edge, Playwright) through its tab buttons and pager, and screenshotted Explain pages beyond page 1.
The five PASS explainers have not changed since round 1 (file times 04:41–05:48), so their round-1 audits (all PASS, 8 sizes) stand.

| slug | round-1 Must fix | round 2 | verdict |
|---|---|---|---|
| falkner_skan_family | D12 step 1 title garbled; Explain §6 "eta_i" | both fixed (`desktop__derive-d2p1.png` now "Evaluate uu_x + vu_y = −(1/ρ)dp/dx + νu_yy (9.9) at the wall"; line 2720 `'$\eta_i = '`); m_sep rtol tightened to 1e-9 | PASS* |
| thwaites_marching | Example 9.2 must use the book's criterion first | fixed: preset `closure: 'white'`; status "book criterion λ = −0.090: x/L = 0.158 · exact Falkner–Skan λ = −0.0681: x/L = 0.126"; Explain §5 "the book's criterion first" with the closed form; chips "book λ = −0.090 / exact FS −0.0681"; cylinder preset 103.1° first; 4 new parity rows against `example_9_2()["criteria"]`; the L(λ) label now "fold n = −0.0904 (λ = −0.0681)" at the red line and "book λ = −0.090" at the amber one | PASS* |
| cylinder_drag_crisis | Explain §2 hint garbled | fixed (renders "St = Ωd/U∞ (4.102) would give St = 2π × 0.2 = 1.26"); default U now 1 m/s (street 266 Hz, "audible range" note) | PASS* |
| wall_jet_invariant | tablet fit failure; bare (9.83) | tablet PASS; D21 result and check now write the implicit solution next to (9.83); step 12 live shows the residual formula; also fixed: `\&` in step 15, K = (ṁ₁/ρ)²/(ν x₁^{1/2}), "velocity scale" wording | PASS* |

\* individually PASS; see the chapter-level Must fix below.

## Chapter-level Must fix (library, confirmed by driving the UI)
**R2-M1. Explain goes stale when only the transport parameter is changed by hand.** Confirmed. In thwaites_marching I pressed → on the cursor
slider until it read 34.8 % (x = 0.087 m, view title λ = −0.101, θ = 0.323 mm); Explain §1–§4 still showed x = 0.0757 m, θ = 0.291 mm
(`scratchpad/stale_thwaites_marching.png`). In falkner_skan_family the slider read n = 3.56 while Explain §1 still said n = 1.89.
The builder's report (λ −0.0611 vs −0.0843) is the same bug.
*Cause:* `assets/viz_lib.js` `app.set` (≈ line 1567): `onlyTime = changed === [transport.param]` → `updateDynamic(true, !onlyTime)` →
`updateCalc(false)`, so only `explain.live` fields refresh and the static Explain HTML is not rebuilt. (During ▶ the loop calls
`updateDynamic(false)`, whose `structural` is undefined → rebuilt every 110 ms, so it looks right while playing and after a click on a
view that changes a second parameter.)
*Affected:* every explainer whose transport parameter is physical and used in Explain HTML — bl_scaling_thicknesses (x),
blasius_similarity_collapse (x), falkner_skan_family (n), thwaites_marching (q), free_jet_similarity (x), wall_jet_invariant (x).
Not affected: cylinder_drag_crisis, karman_street_stability, teacup_secondary_flow (transport is a clock shown through `explain.live`).
The shot audit cannot see this (it sets states, it does not drag).
*Fix (orchestrator / knowledge-keeper, library only):* rebuild the Explain HTML on every user `set()` — e.g. `updateDynamic(true, true)`
from `app.set` when not called from the loop (the loop already rebuilds structurally), or add `transport.clock: true` and treat only
clock parameters as time-only. Then `tools/viz_inline.py --all`, re-run `tools/shot.py --chapter ch09`, and add a machinery test that
moves a transport slider and asserts the Explain text changes. Add a line to `interactive-viz` §9 Lessons.

## Remaining Should-fix (compact)
- **bl_scaling_thicknesses:** `2\,tau_0\,x` typo (Explain §5); amber δ* vs orange θ too similar; phone u/U title clipped.
- **blasius_similarity_collapse:** rescaled collapse drawn from one function, not from each station's data; phone step 2 slices the derivation quote.
- **falkner_skan_family:** 360×640 views ~80 px (hide `curve` on the smallest phones).
- **thwaites_marching:** design note says 60 closure rows, the explainer (correctly) uses 64.
- **cylinder_drag_crisis:** Explain §0 still describes the C_p view in sphere mode; φ_s label vs U arrow overlap in the notebook frame.
- **karman_street_stability:** inspector `'$+\Gamma$'` → "Gamma" (line 2686); 31 vs 61 wavenumbers (Explain §0 / §6); +Γ upper row in D14 vs −Γ in the picture.
- **free_jet_similarity:** "inflow v" and "y stretched" labels collide (laptop); `x^{-1/3}` wraps in tour step 3.
- **wall_jet_invariant:** "peak f′ … at η = 8.11" label overprinted by the f∞ legend; phone bars clip "1.7".
- **teacup_secondary_flow:** blue u(z) curve has no scale and is not named in Explain §0; phone cup labels overprint; unlabelled grey dashed curves; plane (6.2) used for an axisymmetric loop.

## Round-2 verdict: FAIL
All nine explainers meet their own round-1 Must-fix lists and pass lint + a full 8-size shot audit with every parity row ok, but the
chapter cannot pass while R2-M1 makes the Explain tab of six explainers show numbers for a different cursor/station than the picture.
This is a one-line library fix; after `viz_inline.py --all` and a re-run of `shot.py --chapter ch09` plus a manual slider check, the
chapter is PASS without further explainer changes.

## R2-M1 resolution (orchestrator)
- `assets/viz_lib.js` `app.set`: a user-driven change now always calls `updateDynamic(true, true)`, so the Explain HTML is rebuilt
  when only the transport parameter moves (it was treated as "time only"). `updateCalc` still skips DOM work when the HTML is
  unchanged. `tools/viz_inline.py --all` re-inlined 75 files.
- Manual slider check (the reviewer's `stale_check.py`, Edge, keyboard-driven transport slider): Explain static text now changes in
  all six affected explainers (bl_scaling_thicknesses, blasius_similarity_collapse, falkner_skan_family, free_jet_similarity,
  thwaites_marching, wall_jet_invariant); the three clock-driven ones (cylinder_drag_crisis, karman_street_stability,
  teacup_secondary_flow) correctly unchanged.
- `tests/test_machinery.py` 9 passed; `tools/shot.py templates/viz_example.html` PASS; `tools/shot.py --chapter ch09 --quick` 9/9 PASS.
- Remaining Should-fix items above are carried as open items to `knowledge/viz_patterns.md`.

## Verdict: PASS
