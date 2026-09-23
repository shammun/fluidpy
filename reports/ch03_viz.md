# Chapter 3 — explainer review (round 2)                                  2026-09-23

History: round 1 (2026-09-23): 6 PASS, E6 vortex_paddle_wheels FAIL on one Must fix (real-vortex presets
overrode Γ and σ while the sliders showed lab values). Round 2: fixed and verified; all seven PASS.

Reviewer: viz-reviewer (fresh eyes, Phase 7 gate). Re-run for round 2:
- `tools/viz_lint.py --chapter ch03` → all seven `ok`.
- `tools/shot.py` (full: 8 sizes, every tab, step and pager page) on the four changed explainers. All PASS with 0 failures:
  - E6: 392 views, 29/29 rows
  - E5: 456 views, 20/20 rows
  - E4: 376 views, 24/24 rows
  - E7: 344 views, 27/27 rows

  The round-1 full-chapter runs of E1, E2 and E3 (unchanged files) were PASS with 0 failures: 304 + 256 + 184 views and
  16 + 12 + 18 rows. **Chapter total: 145/145 selftest rows ok**, and KaTeX rendered at every size.
- `tools/eq_refs.py --chapter ch03` → the same 4 hits as round 1, each acceptable. E5 L2954 is the move title "Assemble
  (3.19)", and the step body shows (3.19) in TeX. E6 L2941/L2952/L2999 are the `ref:` metadata of Code and Equations
  cards, whose TeX is shown. No bare equation number appears in tour, Explain, Derivation or quiz text.
- Derivations were compared by script with `outputs/ch03/executed.ipynb` in round 1: all 21 assigned D ids match step
  for step (D01 6 · D02 9 · D03 5 · D04 7 · D05 8 · D06 9 · D07 4 · D08 8 · D09 8 · D11 7 · D12 8 · D13 6 · D14 6 · D15 6 ·
  D16 8 · D17 10 · D18 5 · D19 7 · D20 7 · D21 7 · D22 12 · D23 5). Round 2 changed only titles, not steps.
- Book pages read: p098, p101, p105, p107, p109, p110, p112, p115. Every Equations card matches its page, and every
  convention the project pins is honoured:
  - R = G − Gᵀ.
  - ω = 2 × spin.
  - γ = 2S₁₂ (distinct from ch02's Γ ≡ S₁₂).
  - Circulation Γ, with the line-vortex strength B (the book's Fig. 3.16 calls it C).
  - The frame-velocity sign in (3.9).
  - The signed RTT wall term.
  - (3.6) written with F.

## Summary

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| flow_lines_unsteady (E1) | ok | PASS · 8 · 304 | 16/16 (1e-12 closed forms; 1e-7 RK4 vs solve_ivp) | 8 sections + interpretation, 5 live | D03 5 ✓ · D04 7 ✓ · D05 8 ✓ | very good | minor (D04 step 3 live numbers use t′ = 0) | **PASS** |
| material_derivative_probe (E2) | ok | PASS · 8 · 256 | 12/12 (1e-10…1e-12; stencils 5e-5, justified) | 9 sections + interpretation, 19 live | D01 6 ✓ · D02 9 ✓ | very good | minor (no direction cue on the negative bar; label collision) | **PASS** |
| galilean_frames_cylinder (E3) | ok | PASS · 8 · 184 | 18/18 (1e-12; 1e-6 vs fluidpy FD, honest) | 9 sections + interpretation, 6 live | D06 9 ✓ | very good | minor (legend strip over the particle on phones) | **PASS** |
| fluid_element_deformation (E4) | ok | PASS · 8 · 376 | 24/24 (1e-12 / 1e-10; 1e-8 measured) | 8 sections + interpretation, 4 live | D07 4 ✓ · D08 8 ✓ · D09 8 ✓ · D11 7 ✓ | very good | step 1 now shows its question on page 1 at 360×640 ✓ | **PASS** |
| spin_and_principal_axes (E5) | ok | PASS · 8 · 456 | 20/20 (1e-12 / 1e-10) | 9 sections + interpretation, 7 live | D12 8 ✓ · D13 6 ✓ · D14 6 ✓ · D15 6 ✓ · D16 8 ✓ | excellent | "Rate k" with per-flow meaning ✓; D13 title ✓; 360 px labels ✓ | **PASS** |
| vortex_paddle_wheels (E6) | ok | PASS · 8 · 392 | 29/29 (1e-12; brentq/Lambert W 1e-10) | 10 sections + interpretation, 3 live | D17 10 ✓ · D18 5 ✓ · D19 7 ✓ · D20 7 ✓ | very good | Must fix closed ✓; empty phone band gone ✓ | **PASS** (round 2) |
| reynolds_transport_cv (E7) | ok | PASS · 8 · 344 | 27/27 (1e-9…1e-12 quadrature; FD 1e-6/1e-7) | 9 sections + interpretation, 4 live | D21 7 ✓ · D22 12 ✓ · D23 5 ✓ | excellent | short titles ✓; "now" box on page 1 ✓; Explore 15 → 11 pages ✓; display-style integrals on desktop ✓ | **PASS** |

**Checks that apply to all seven.**
- Tabs present: Walkthrough · Explore · Explain · Derivation · Equations · Code · Check yourself.
- Each has ≥ 6 depth features: all have transport, presets, status, inspector, linked views and modes; E2, E3 and E7
  add term bars; E2, E5 and E6 add notes.
- 4 check questions each, all answerable with a preset.
- Step 1 asks a plain-words question in every explainer.
- Audit: no overflow, no text below 12 px and no view below 60 px at any size.
- Every `py:` row calls `ch03.*` at non-trivial inputs against the explainer's own JS function.
- No book prose, figures or tables are reproduced.

## flow_lines_unsteady (E1)
**Must fix** — none.

**Should fix**
1. `desktop__derive-d2p3.png` (D04 step 3): the live line computes the path-circle centre for the orange particle
   (t′ = 0), while the answer-sheet view beside it shows t′ = t = 1.5 s. Label the live line "for the orange particle
   (it passed the port at t′ = 0)".
2. `notebook__equations.png`: the ~70 px velocity window has touching y ticks. Draw only −1, 0 and 1 below ~90 px.

**What works**: the answer-sheet view beside the experiment; the ringed dye dot, whose position (1.37, −1.37) m matches
the inspector arithmetic; the pre-filled streak circle; loops/no-loops presets at the U₀ = ωξ_o threshold.

## material_derivative_probe (E2)
**Must fix** — none.

**Should fix**
1. `desktop__tour-step4.png`: the −0.36 advective bar is drawn as a waterfall segment pointing the same way as the
   +0.36 local bar. Add a left-pointing arrowhead or a hatch for negative segments (as the E7 "swept out" bar has).
2. `desktop__derive-d2p6.png`: the probe label is cut by the float marker. Offset it when the two markers are closer
   than ~60 px.

**What works**: the float's independently measured rate (◇) lands on the DT/Dt bar; the stalled front reverses the
roles; the stretching-map t–x diagram; the (3.6) note about the missing F.

## galilean_frames_cylinder (E3)
**Must fix** — none.

**Should fix**
1. `phone-tall__tour-step3.png`: the in-canvas legend strip covers the top ~40 px of the flow view, where the tagged
   particle sits. On portrait, move the legend into the view title.
2. At 360×640 Explain runs to 21 pages. Merge sections 5 and 6 on phones.

**What works**: the local and advective bars trade places while the purple total and the dashed body-frame line stay
put; the u = U + u′ vector view; the lake, body, halfway and overtake presets.

## fluid_element_deformation (E4)
**Must fix** — none. (Round-1 Should fix closed: step 1 is 24 words, and its question is on page 1 of 3 at 360×640,
`phone__tour-step1.png`.)

**Should fix** — none new.

**What works**: the measured stretching and closing dots land on the n·S·n and 2n₁·S·n₂ curves; the shear corner
reaches 68.2° after 0.4 s; the area graph against e^{t tr G}; the S + ½R = G split in the R = G − Gᵀ convention.

## spin_and_principal_axes (E5)
**Must fix** — none. The round-1 Should-fixes are closed:
- The slider is now "Rate k". `kMeaning()` agrees with `Gof()`: shear gives γ = k, solid body ω₀ = k, and pure strain
  S₁₁ = −S₂₂ = k/2 (`desktop__explain.png` §1).
- The D13 title reads "ω′ = ω − 2Ω (z-components)".
- The 360 px labels no longer collide ("1: −0.433 / 2: +0.433 / avg 0.000", `phone__explore.png`).

**Should fix** — none new.

**What works** (candidate for viz_patterns): single-thread rates against the flat pair-average line (the "why ½ω"
aha); the co-rotating observer that stops the wheel; the split view with a three-row ellipse table.

## vortex_paddle_wheels (E6)
**Must fix** — none. **Round-1 Must fix closed**:
- `desktop__explore.png`: in cyclone mode the log sliders read Γ = 1.8×10⁷ m²/s and σ = 36000 m, the values the
  picture uses.
- `phone__explore.png`, `phone-tall__tour-step7.png`: the status line leads with "🌀 cyclone Γ 1.8×10⁷ m²/s, σ 36000 m"
  at every size.
- 5 new selftest rows, all ok.
- `phone-tall__derive-d3p2.png`: the empty band is gone, and the profiles view now also carries Γ(r)/Γ for D19.

**Should fix**
1. The status line and the σ chip print "36000 m". Use `Viz.fmt` with a km switch ("36 km"), as the walkthrough text
   already does.

**What works**: wheels that orbit without turning in the line vortex while the loop still reads Γ = 2πB; the log–log
mean-vorticity plot with slope −2; the D17 sector drawn on the vortex; the real-vortex table with the current row.

## reynolds_transport_cv (E7)
**Must fix** — none. The round-1 Should-fixes are closed:
- The D21/D22 titles are short, and (3.30)/(3.35) are written out on the goal pages.
- `phone-tall__derive-d2p3.png`: the "now" box is on page 1.
- `phone__explore.png`: Explore is 11 pages, down from 15, with per-mode control hiding.
- `desktop__derive-d2p6.png`: integrals render in display style on desktop.

**Should fix**
1. `phone-tall__tour-step3.png`: the precomposed "ḃ" loses its dot at phone size. Write `$\dot b F(b)$` in TeX.
2. `notebook__equations.png`: the budget labels are cut inside ~45 px bars. Use 2 significant figures, or put the
   labels above the bars when a bar is under 60 px wide.
3. The D22 result page is still 19 pages at 360×640 (the whole 12-line chain). Consider showing the chain collapsed to
   its first and last lines on phones.

**What works** (candidate for viz_patterns): swept bands coloured by the sign of b·n; the budget waterfall landing on
the measured ◇; the dropped-term log–log panel (T4 at slope 2); three geometries on one engine.

## Library bugs (for the orchestrator / knowledge-keeper; `assets/viz_lib.js`)
Still open. The chapter explainers now work around them.
1. **Pager has no keep-with-next.** A heading, or a derivation move line, can end a page alone.
2. **Packed height ≠ shown height.** A page packed at 145 px is shown at 135 px because layout runs before the final
   `fit()`. Re-pack after `fit()`.
3. **Tall display integrals clipped at a slice break.** The slice boundary ignores the KaTeX `.vlist` overhang.
4. **`optional` controls hidden on phones even when their value drives the picture.** E6 now works around this through
   the status line.

## Best screenshots of the chapter
- `reports/viz/ch03/flow_lines_unsteady/desktop__tour-step1.png`
- `reports/viz/ch03/reynolds_transport_cv/desktop__tour-step5.png`

## Verdict: PASS
All seven explainers PASS. The round-1 Must fix (E6) is closed and nothing regressed. The Should-fix items and the
library bugs do not block the gate.
