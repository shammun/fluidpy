# Chapter 2 — explainer review (round 1)                          2026-09-16

Reviewer: viz-reviewer (fresh eyes, Phase 7 gate). Everything below was re-run today:
- `.venv/Scripts/python.exe tools/viz_lint.py --chapter ch02` → all five `ok`.
- `.venv/Scripts/python.exe tools/shot.py --chapter ch02` (full run, 8 sizes, every tab, walkthrough step and derivation
  page; 1 496 views in total) → all five PASS, 0 failures, 127/127 parity rows ok, KaTeX rendered at every size.
- Screenshots read for every explainer: `desktop__tour-step1…8`, `phone-tall__tour-step*`, `phone__explore`,
  `phone-land__explore`, `notebook__equations`, `laptop__check`, `desktop__explain`, `phone-tall__explain`, the Code tab,
  and the derivation pages named below at desktop and phone-tall.
- Physics read against the rendered book pages `chapters/pages/ch02/p069, p070, p073, p075, p076, p081, p082, p084,
  p085, p086, p088, p089` (Eqs. 2.4–2.8, 2.12, 2.15, 2.24–2.29, 2.30–2.35; Examples 2.1–2.6).
- Derivation step counts and move titles compared programmatically with `analysis/ch02_design.md` Part F: all 12
  derivations match (D01 6 · D02 8 · D05 9 · D06 8 · D12 5 · D14 6 · D15 7 · D17 15 · D21 7 · D22 9 · D25 9 · D26 10).

## Summary

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| rotation_of_axes (E1) | ok | PASS · 8 · 328 | 20/20 (rtol 1e-12; 4 invariants) | 10 sections + interpretation, 2 live | D01 6 ✓ · D02 8 ✓ · D06 8 ✓ | good story (arrow vs ruler, C twice for a tensor, polar = rotated Cartesian) | transport autoplays under the trace steps; labels collide in the plane view; matrix assembly line clipped on short views | **FAIL** (2 Must) |
| cauchy_traction_principal_axes (E2) | ok | PASS · 8 · 328 | 27/27 (rtol 1e-12; 3 invariants) | 9 sections + interpretation, 4 live | D05 9 ✓ (2-D tetrahedron sketch, shrinks) · D17 15 ✓ (live numbers on 6 steps) | very good — the reference-depth panel of the chapter | 4 *why* lines > 35 words; strain-mode axis labels still say σn/τs | **PASS** |
| strain_vs_rotation_split (E3) | ok | PASS · 8 · 232 | 25/25 (rtol ≤ 1e-9; 2 invariants) | 8 sections + interpretation, 10 live | D14 6 ✓ · D15 7 ✓ | very good; conventions pinned exactly (A = ½R, ½ω, never "ω") | motion view small on portrait phones in "all" mode; end card covers the G panel's axis | **PASS** |
| gauss_flux_box (E4) | ok | PASS · 8 · 352 | 25/25 (rtol ≤ 1e-9 except the singular-source row at 2e-4, justified) | 8 sections + interpretation, 5 live | D25 9 ✓ · D22 9 ✓ · D21 7 ✓ | very good (tiles, Taylor faces, shrink, singular source) | bars gap label wrong in the singular regime; D21's key number hidden on portrait phones | **FAIL** (2 Must, both one-line) |
| stokes_circulation_loop (E5) | ok | PASS · 8 · 256 | 30/30 (rtol ≤ 1e-8; 5 invariants) | 8 sections + interpretation, 5 live | D26 10 ✓ (tiles drawn, staircase → circle) · D12 5 ✓ | very good ("curl without curves", Stokes failing on purpose) | loop small on portrait phones (x-range grows to ±5); step 6 points at a hidden view | **PASS** |

**Checks that apply to all five:** KaTeX renders everywhere (no fallback text); no overflow, no text below 12 px and no
views below 60 px at any of the 8 sizes (audit); tabs Walkthrough · Explore · Explain · Derivation · Equations · Code ·
Check yourself present; ≥ 2 depth features each (all five have presets, status, inspector, notes, transport; E1/E3 add
modes, E2/E3/E4/E5 add term bars); every `py:` expression exercises the explainer's own JS function against
`ch02.*` at non-trivial inputs (30° cuts, 0.7 rad rotations, off-centre loops, the smooth test field); tolerances are
honest (1e-12 for closed forms, 1e-10 for `expm`/quadrature, 2e-4 only for the point-source flux with 32 midpoints).
Equation numbers cited match the rendered pages: (2.5)/(2.6)/(2.7)/(2.8) p069–070, (2.12) p073, (2.15) p075,
(2.24)–(2.25) p081, (2.26)–(2.29) p082 (the "i−1" misprint is flagged), (2.30)–(2.33) p085–086, (2.34)–(2.35) p088;
orthogonality is credited to Exercise 2.8, not (2.6); Ex. 2.4 uses Γ ≡ S₁₂ with the book's "2S₁₂ = Γ" line flagged;
Ex. 2.6's second bracket is read as u_y. No book prose, figures or tables are reproduced.

**Conventions pinned by the chapter — all honoured:** C_ij = e_i·e'_j passive, x' = Cᵀx (E1 `transformVector` sums the
first index; the active toggle is labelled "Wikipedia R" and computes Cx); f_i = τ_ji n_j contracts the first index (E2
`traction` uses `tau[j][i] * n[j]`; the non-symmetric demo shows τ·n differing by √3); the book's R = G − Gᵀ carries
ω = ∇×u while A = ½R carries ½ω (E3 `omegaOf(A) = A₂₁`, readout "Spin ½ω₃ = A₂₁", matrix line "½ω₃ = A₂₁ = −A₁₂ … (ω = ∇×u)");
presets pure_strain = diag(1, −1), uniaxial_extension = diag(1, 0), irrotational_strain = [[0, 1], [1, 0]]; Stokes
n_c into A and t = n_c × n counterclockwise about n (E5 `tangentOf([-1,0,0],[0,0,1]) = (0, 1, 0)`, flip reverses both
sides together).

---

## rotation_of_axes (E1)

**Must fix**
1. **The transport autoplays under the trace steps, so the numbers the text quotes are not the numbers on screen**
   (`desktop__tour-step1…4`, `phone-tall__tour-step3`). The config has no `autoplay: false` (E2 and E4 set it) and
   steps 1–3 have no `play: false`, so θ sweeps from the moment the page opens: step 1 shows θ = 35.5°, step 3 says
   "x′₁ = 1×0.866 + 2×0.5 = 1.866" while the plane reads x′₁ = 2.215 at θ = 55.5° and the Eq. (2.5) card reads
   (2.175, 0.520); step 4 says "length 2.236 is the same" while the frame keeps moving; in step 2 the inspector
   (rendered once) says "cos(128°) = −0.623" while the matrix cell shows −0.719 at θ = 46°. This is the ch01 lesson
   "a step that traces one sample's arithmetic sets the global parameter to that sample's value" — the value is set
   but immediately overwritten. Fix: add `autoplay: false` to the app config and `play: false` to steps 1, 2, 3, 6, 7
   and 8 (step 4 keeps `play: true`, it says "Press ▶ to sweep"); also have the inspector re-render while playing or
   pause on cell click (it already calls `A.play(false)` on plane drags — do the same for the matrix click).
2. **Labels collide in the plane view at the default angles** (`desktop__tour-step3`, `desktop__tour-step4`,
   `desktop__tour-step6`, `phone-tall__tour-step3`). Vector mode: the cell-inspector arc label "cos 55.5° = 0.567"
   is drawn at radius 0.7 where the shadow feet of x′₂ and x₁ sit, so it covers "x′₂ = 0.309" and "x₁ = 1.000"
   (phone: "cos 60° = 0.493" hides x′₂ completely); the θ arc label "θ = 35.5°" (radius 0.32 R) sits on "x₁ = 1.000";
   in step 4 "x′₁ = 2.236" overprints the "x" tip label and "x₂ = 2.000". Tensor mode (step 6): "τ′₁₂ = 0.500" and
   "τ′₁₁ = 0.866" overlap each other on the 1′ face and the orange arc label "cos 120° = −0.500" sits on them; in
   step 7 "τ₁₂ = 1.000" is overdrawn by the rose arrow. Fix: put the inspector arc label outside the unit-vector pair
   (radius ≥ 1.3 in vector mode, on the far side from the shadows) and the θ label beyond the 1′ arrowhead; offset
   shadow labels perpendicular to their axis (teal below/left, orange above/right) and, in tensor mode, place the two
   primed labels on different faces (τ′₁₁ on the +1′ face, τ′₁₂ on the −1′ face) as the unprimed ones already do.

**Should fix**
- The matrix window's third assembly line is clipped by the canvas bottom on short views: `phone__explore` (360×640)
  and `phone-land__explore` show "u_θ = u₁C₁₂ + u₂C₂₂" cut in half. `drawMatrix` sizes the cells from
  `(H − textH − hh − 18)/n` but then draws three lines of `lh` each below `y0 + n·cw + 10`; when `H < ~180` drop the
  third (orthogonality) line — it is already in the status badge — or reduce `lh`.
- *Why* over 35 words: D01 step 5 (36), D02 step 8 (40). Move the consequence into *in words*.
- Step 7 ("Drag θ until τ′₁₂ = 0 — the rose curve crosses zero") talks about the curves view, hidden on portrait
  phones; the `tau` readout is shown, so say "…or watch the τ′₁₂ readout reach 0" and keep "τ′₁₂ = 0 at 45°" in the
  matrix title (it is there).
- The vector-mode plane at desktop is a 3.6-unit square in a wide panel: the arrow x = (1, 2) uses a third of it.
  Consider `planeRange` 3.0 for |x| ≤ 2.5.
- Equations tab `x8` live line prints four residuals in a two-row `aligned` block that wraps on the notebook size;
  fine, but a one-column table would read better.

**What works** (keep; candidates for `knowledge/viz_patterns.md`)
- The passive/active table with the current row highlighted and the toggle "what rotates" — the confusion is settled
  by a picture, not a sentence; `activeRotate(θ) = transformVector(−θ)` is pinned by a selftest row.
- The "Is it a vector? test of Eq. (2.8)" table (four candidate triples, residuals from the same C) and the
  "fixed array (not a tensor)" preset with its (2.12) residual in the status badge — the chapter's CORE idea taught by
  a counter-example.
- Direction-cosine matrix as clickable cells with the inspector drawing the two unit vectors and their angle.
- Mirror toggle wired to D02 step 8 (det C = −1 while CᵀC − I stays 0, status turns amber).

---

## cauchy_traction_principal_axes (E2)

**Must fix** — none.

**Should fix**
- *Why* over 35 words in D17 steps 5 (39), 8 (42), 9 (39), 15 (42); e.g. move "A real λ makes the linear system real,
  so b can be chosen real" (step 5) and the Gram–Schmidt remark (step 8) into *in words*.
- Strain-rate mode still labels the Mohr axes "σn [1/s]" / "τs [1/s]" and the readouts "Normal part σ_n", "Shear
  part τ_s" (`desktop__tour-step8`); use "n·S·n" / "s·S·n" (or "stretch rate along n" / "shear rate") when
  `s.mode === 'strain'`, as the element labels already switch to S₁₂/S₂₁.
- `desktop__tour-step1`: the element is a 300 px square in a 300 × 540 px panel — half the panel is empty above and
  below. Draw the element at the panel's full height (equal aspect, `pad` computed from `v.h`) or move the
  "grey half…" caption into that space at a larger size.
- The D17 result page inherits `play: true` from step 14 (`desktop__derive-d2p16` shows φ sweeping while "your
  numbers" are read); set `play: false` on the result as step 15 does.
- Explore at 844×390 hides the step/speed buttons (engine), fine; the status line "↗ f = (0.500, 0.866) 1/s …" wraps to
  two lines at 1000×700 — the mid-length variant could trigger at `Wd < 1100`.

**What works**
- The element view with the grey removed half, n, f and its blue/rose decomposition, plus σn = n·τ·n term bars
  (τ₁₁n₁², (τ₁₂+τ₂₁)n₁n₂, τ₂₂n₂²) that add exactly — a stress tensor made visible.
- Three linked views on one φ (element, σn/τs curves with λ bounds and the special-cut ticks, Mohr's circle with the
  double angle) and a status verdict with correct thresholds (principal within 0.5°, max-shear within 0.5°,
  hydrostatic when radius < 1e-9, non-symmetric warning).
- D05's 2-D tetrahedron sketch that shrinks with the step (h = 1 → 0.5 → 0.25) while the face arrows stay and the
  dashed volume term shrinks — the orders-of-smallness argument seen, not stated.
- D17 with live numbers on six of fifteen steps (discriminant, b_max·b_min = 0, τ′ = diag, c_k, σn = Σλc², shear bound).
- The non-symmetric toggle drawing τ·n in amber next to n·τ and the Mohr view saying "not a circle" — the first-index
  convention taught by contrast.

---

## strain_vs_rotation_split (E3)

**Must fix** — none.

**Should fix**
- On portrait phones the motion view is small when all three squares are overlaid: `phone-tall__tour-step5` gives the
  square a ~150 px plot (squares ~40 px), `phone-tall__explain` ~110 px (squares ~25 px, the lean is not readable).
  The single-panel layout also caps the plot width at `v.h·1.2 + 44` and centres it, leaving wide margins. On
  portrait give the motion row more height (`rows: [2, 1]`) and let the panel use the full width; or hide the matrix
  view's fact lines below ~120 px so the row can be shorter.
- The end-of-run card sits over the G panel's lower-left quadrant and covers the x₁ axis labels
  (`desktop__derive-d1p6`, `d1p7`). Place it in the empty band above the panels (the title row) or in the A panel's
  free corner.
- In "G only" mode at desktop (`desktop__tour-step1`) the single panel is ~340 px wide in a 930 px view; the
  streamlines fill it but two thirds of the view is blank. Consider drawing the S and A ghosts' panels faintly beside
  it, or widening the single panel (`pw = Math.min(v.w, v.h * 1.6)`).
- The `S:A` term is labelled "2 S:A (2.29) — always 0" so that the bars add to G:G; say so in the terms title
  ("S:S + 2 S:A + A:A = G:G") so a reader does not look for a missing factor.

**What works**
- The whole explainer answers its title with one clock: black square leans, teal stretches along ±45°, orange spins;
  the walkthrough (G alone → split → S alone → A alone → all → blind → your turn) is the cleanest story of the five.
- Conventions exactly as pinned: ½ω₃ = A₂₁ everywhere (readout, matrix line, status, Explain §3, D15 goal/result),
  R = G − Gᵀ named as the book's tensor with ω = ∇×u, and selftest rows "vorticity ω3 = vector(R)" and "R = 2A" pin it.
- `expm2` closed form (cosh/cos/shear-limit branches) matched to `scipy.linalg.expm` to 1e-10 at three regimes, and
  the tracer inspector that shows e^{Mt} and x(t) under G, S and A for the clicked dot.
- Frobenius term bars S:S, A:A, 2S:A → G:G with the S:A bar pinned at 0 for every slider setting (Eq. 2.29 as a
  live invariant); the material-line angle view with the "initial rate" ghost.
- The Explain §5 "where the square is at t" section: e^{λ₁t}, e^{λ₂t}, the area factor e^{(∇·u)t}, the turn ½ω₃t and
  the shear displacement Γt·x₂ — every number on the stage accounted for.

---

## gauss_flux_box (E4)

**Must fix**
1. **The bars-view gap label gives the wrong reason in the singular regime** (`desktop__tour-step7`,
   `phone-tall__tour-step7`): with the point source enclosed it reads "|∮ − ∬| = 1 m²/s (midpoint rule, n = 8:
   ∝ 1/n²)", i.e. it tells the reader the 1 m²/s discrepancy is quadrature error — the opposite of the step's lesson
   (the status badge and the notes say the right thing). In `drawBars`, when `p.singular` print
   "|∮ − ∬| = 1.002 m²/s: the delta at the origin carries the flux m — Gauss needs a smooth Q" (colour `C.interior`),
   and keep the 1/n² text for the non-polynomial smooth field only.
2. **D21 step 4's key number is not visible on portrait phones** (`phone-tall__derive-d3p4`, also `phone__explore`
   and `phone-tall__explain`): the WATCH line says "At h = 0.05 m the ratio (1/A)∮ in the bars title has converged to
   ∇·Q(x₀)", but on phones the bars title is truncated to "… ∮ = 0.004 = ∬ ∇·Q dA = 0.004…" — the ratio and ∇·Q(x₀)
   at the end are cut off, the limit view is hidden and the status badge is not shown on the Derivation tab. Same
   ch01 lesson (put a hidden view's key number in a visible place). Fix: on narrow views lead the bars title with
   the ratio ("(1/A)∮ = 1.540 → ∇·Q(x₀) = 1.540 1/s · ∮ = 0.004 = ∬ …"), and add a `live` line to D21 step 4
   (`(1/A)∮ = … at h = … vs ∇·Q(x₀) = …`, as steps 3 and 6 already have).

**Should fix**
- The field view is small on portrait phones: `phone__explore` (360×640) shows a ~280 × 140 px plot with a ~70 px box,
  `phone-tall__explain` ~330 × 190 px; `fieldPlot` caps the plot to `min(v.w, v.h + 80)` and centres it. Give the
  field row more height on portrait (`rows: [2, 1]`) or let the heatmap fill the panel width with the box drawn at
  equal aspect inside it.
- The status text for the singular source is three lines long on phones (`phone-tall__tour-step7`); a short variant
  below 520 px ("⚠️ singular inside: ∮ = 1.002 ≈ m whatever h") as E2 does.
- `desktop__tour-step5` starts the shrink from h = 2.5 m, so for the first second the box grows (h 2.256 → 2.09 …)
  while the text says "the box shrinks"; start the sweep at `shrink: 0` (h = 1) or say "the box first grows, then
  shrinks".
- Explain §2 prints "flux = Q̄₁ h = 0.739 × 2.000" for the smooth field — good — but for the linear field the
  "constant along this face" reason and the `F.poly` test would also flag (x², 0) as exact along vertical faces
  (it is), fine; only the `Explain §3` special-case `F.chip !== '(x², 0)'` reads oddly. Minor.

**What works**
- Face arrows coloured by sign (orange out, blue in) on the box and the waterfall bars left → bottom → right → top →
  ∮ next to ∬, with the interior-faces bar appearing when tiled — Gauss' theorem as bookkeeping you can watch.
- Tiling with the shared faces drawn as opposite pairs and a shared-face inspector that prints "+0.250 × 0.250 … from
  the left tile, −… from the right tile, sum = 0 — D25 step 8".
- Ex. 2.5's Taylor face values printed on the faces ("Q₁ = 0.500 + 0.500 = 1.000") wired to D22 with live numbers on
  five steps.
- The limit view: (1/A)∮ vs h on a log axis with the exact value dashed, "exact at every h" for polynomial fields,
  the slope-2 log–log inset for the smooth field and the 1/h² blow-up for the singular source — three regimes, one
  picture.
- `Viz.CMAPS.bwo` diverging map with a white zero (promotion candidate), `labIn` (labels clamped inside the plot).

---

## stokes_circulation_loop (E5)

**Must fix** — none.

**Should fix**
- On portrait phones the field's x-range expands to ±5 m to keep equal aspect in a short panel (`phone-tall__explain`,
  `phone-tall__tour-step5`, `phone-tall__derive-d1p8`): the loop is ~80 px and the unit rectangle ~50 px. Give the
  field row more height on portrait, or clamp the view to the box [−2, 2]² and let the heat image fill the panel.
- Walkthrough step 6 ("watch Γ/A hold at the point value — the dot on the flat curve") refers to the bars view, hidden
  on portrait phones; add `readouts: ['ratio', 'curl']` to the step (the field title also carries Γ/A).
- The u·t view on `phone-tall__explain` is a ~60 px strip and its "shaded area = Γ" label sits on the curve; below
  ~80 px drop the label (the title already has Γ).
- D26 step 9 sets `tiles: 6` on the circle: "32 tiles: Σ Γₖ = −0.889 · staircase −0.889" vs the circle's Γ = −0.785
  (`desktop__derive-d1p9`). The WATCH says they converge as the tiles shrink — true, but a reader sees a 13 % gap;
  use `tiles: 8` (the slider max) or state the gap in the live line ("staircase −0.889 → circle −0.785 as k → ∞").
- The "⚠ singular core" label at the origin and the wheel-rate label at the bottom right are fine; the orientation
  label "n = +e₃ (out of screen) · t counterclockwise" and the tiles label overlap the top streamlines slightly —
  cosmetic.

**What works**
- Curl heat image (purple/teal with a white zero, `curlImage`), streamlines, RK4 tracers and the paddle wheel at
  ½(∇×u)₃ on one stage; the shear preset makes "curl without curves" visible in one glance.
- The u·t "unrolled loop" view whose shaded area is Γ, split into Ex. 2.6's four sides in rectangle mode with their
  sums as term bars (hatched for negative values), and the loop-point inspector giving one term of the 400.
- Stokes failing on purpose: the vortex preset with the core inside shows Γ = 2πK against curl flux 0, the bars
  "core inside: sides differ by 6.283", the Γ/A ∝ 1/A curve and `hypothesis_ok = False` pinned by a parity row.
- D26 with the tiles drawn (shared edges walked twice, blue/orange arrows), the staircase-to-circle refinement and
  the ∇φ corollary that turns the field white.
- Flip-n toggle changing both sides of (2.34) together, with the status saying so.

---

## Library-promotion candidates reported by the builders (for the knowledge-keeper)

| Helper | From | Note |
|---|---|---|
| `tipLabel(P, from, to, str, color)` — label pushed past an arrow tip | E1 | useful, but see E1 Must 2: needs a collision-avoiding offset |
| `drawArc(P, a0, a1, r, color, label)` — angle arc with label | E1 | promote with a `labelRadius` option |
| `fitText` (try progressively shorter alternatives until one fits) | E1 `drawMatrix` | promote as `Viz.fitText(ctx, alts, maxW)` |
| `drawMatrix` cells with row/column highlight | E1, E3 (`matrixLayout`) | second use in ch02 → promote `Viz.matrix(ctx, rect, M, {hot, colors})` |
| `hideLandscape` / collapse-row engine flag | E1 | engine feature request |
| `fx/tx` fixed-decimal formatters (HTML minus vs TeX minus) | E2, also E3 `f2z/tz`, E4 `fz/tz`, E5 `fz/tz` | four copies → promote `Viz.fixed(v, d, {tex})` |
| `halfSquare` / `fillPoly` | E2 | `Viz.fillPoly(P, pts, color, alpha)`; clipping stays local |
| adaptive status (short / mid / full by `window.innerWidth`) | E2 | promote as `status: s => ({short, mid, text})` in the engine |
| `terms.unit` getter | E2 | engine: allow `unit: s => …` |
| readout wrap CSS | E2 | engine |
| `expm2`, `principal2d`, `eps3` | E3 (`principal2d` also in E2) | promote `Viz.num.expm2`, `Viz.num.principal2d`, `Viz.num.eps3` |
| `Viz.card` `align: 'left'` | E3 | promote |
| `Viz.CMAPS.bwo` diverging blue–white–orange | E4 | promote; also make `Viz.field.heatmap` accept a diverging map with `vmin = −vmax` |
| `labIn` (text clamped inside the plot rect) | E4 | promote `Viz.textIn(P, str, x, y, opt)` |
| xlog axis (manual log ticks) | E4 `drawLimit`, E5 `drawBars` | second use → `P.axes({xlog: true})` |
| shrink-transport pattern (transport on log h) | E4 | pattern note |
| `curlImage` diverging heat image with white zero | E5 | merge with `bwo` above |
| `.viz-seg` wrap | E5 | engine CSS |
| `tilesOf` / `tileSums` | E5 (E4 has `tiled`) | promote `Viz.num.tiles(rect|disc, k)` |
| hatched bars | E5 (ch01 E3 `hatchPattern`) | second use → `Viz.hatch(ctx, color)` |

## Verdict: FAIL

E2, E3 and E5 PASS. E1 and E4 each have two Must-fix items (E1: `autoplay: false` + `play: false` on the trace steps
and label placement in the plane view; E4: the singular-regime gap label and D21 step 4's number on portrait phones).
All four fixes are local edits; re-run `tools/shot.py --chapter ch02` afterwards and re-check
`desktop__tour-step1…4`, `phone-tall__tour-step3`, `desktop__tour-step6` (E1) and `desktop__tour-step7`,
`phone-tall__derive-d3p4` (E4).
