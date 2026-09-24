# Chapter 6 — explainer review                                  2026-09-24

Phase 7 gate, fresh-eyes review of the 9 explainers against `analysis/ch06_design.md` Part B (E1–E9), the errata in
Part G (G1–G8) and curation §5. Re-run today, after the fluidpy changes (ρ = 1.2 kg/m³ default for the 2-D force
helpers, `superposition_state` net-sink status, new `blasius_state` keys, the odd-N branch of the axial solver):
`tools/viz_lint.py --chapter ch06` (9 × ok), `tools/shot.py` on all 9 files (8 sizes each, every tab, step and pager
page), `tools/eq_refs.py --chapter ch06`, a scan for single-backslash TeX in the chapter scripts, and a look at about
70 screenshots, the rendered book pages (6.3, 6.10, 6.24/6.27, 6.46, 6.57–6.62, 6.68–6.69, 6.72–6.73, 6.93–6.95,
6.103–6.109) and fluidpy spot checks (relaxation sweep counts, axial_state numbers).

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| superposition_sandbox (E1) | ok | PASS 8 / 336 | 37/37 | 8 §§ + interpret | D05 5, D06 9, D07 10 ✓ | strong | 2 small overlaps, KaTeX glyph warnings | **PASS** |
| cylinder_circulation_lift (E2) | ok | PASS 8 / 328 | 32/32 | 8 §§ + interpret | D09 6, D10 9, D11 9: **2 *why* texts render as TeX errors** | strong | — | **FAIL** (r1) → **PASS** (r2) |
| vortex_wall_images (E3) | ok | PASS 8 / 256 | 28/28 | 8 §§ + interpret | D12 5, D13 11 ✓ | strong | clean | **PASS** |
| complex_potential_corners (E4) | ok | PASS 8 / 240 | 30/30 | 7 §§ + interpret | D14 8, D15 7 ✓ | good (step 2 could show convergence) | phone CR label overlap | **PASS** |
| blasius_kutta_contour (E5) | ok | PASS 8 / 384 | 30/30 | 8 §§ + interpret | D16 8, D17 12, D18 11 ✓ | excellent | clean | **PASS** |
| conformal_joukowski (E6) | ok | PASS 8 / 344 | 30/30 | 8 §§ + interpret | D19 6, D20 8, D21 11 ✓ | strong | tracers not clipped; **eq_refs hit** | **FAIL** (r1) → **PASS** (r2) |
| laplace_relaxation (E7) | ok | PASS 8 / 248 | 36/36 | 7 §§ + interpret | D22 6, D23 9 ✓ (one tolerance inconsistency) | strong | phone label overlap | **PASS** |
| axial_singularity_bodies (E8) | ok | PASS 8 / 264 | 33/33 | 8 §§ + interpret | D26 10, D27 6 ✓ | strong, honest (G3 applied) | stray contour spurs | **PASS** |
| added_mass_sphere (E9) | ok | PASS 8 / 464 | 35/35 | 8 §§ + interpret | D28 5, D29 14, D30 9, D31 10 ✓ | strong | phone label clipped | **PASS** |

Parity: every row is `ok`. Closed-form rows use rtol 1e-10…1e-14. The only looser tolerances are honest convergence
invariants (E1 pair-error ratio ¼ at rtol 0.02, E8 moment → exact at rtol 1e-3). The rows call the explainers' own
drawing functions at non-trivial inputs, and they cover the changed fluidpy paths: E1 has an exact-text row for the
net-sink status, E5 rows use the new `blasius_state` keys (F_perp, F_par, c0, D/L_pressure, offsets), and E2, E5 and E6
pass ρ = 1.2 explicitly. E8 only uses even N (slider step 2), so the new odd-N least-squares branch is never reached from
JS.

Sign conventions: all checked. The book's Γ is clockwise in (6.36)–(6.40), (6.52), (6.62) and (6.68), and
L = +ρUΓ_cw (E2, E5, E6). E3's eddy is −Γ with a +Γ image. E1 uses counterclockwise Γ and labels it so. Real slips are
shown in corrected form: (6.61) in E5 (the 1/z² coefficient −(Ud/π + Γ²/4π²), with the printed form noted as a slip
that changes nothing), (6.104) in E9 (the printed −u_s/a³ bracket, with a live "printed bracket" number) and (6.108) in
E9 (the leftover dφ). (6.82) is not presented as a slip anywhere (G6). G1, G3, G5 and G7 are applied: E6's principal
root fails at every left-half point off the slit, E8's alternating strengths / airship bars / panel 1/N² only on the
ellipse / N ≤ 40, E2's separated band is marked "qualitative sketch", and E5's text follows the tilted-stream mode.
None of the explainers copies book prose, figures or tables. Example 6.2 uses Q = 1 (not the book's value), and
Example 6.1 uses Γ = h = 1.

---

## superposition_sandbox (E1) — PASS
**Must fix**: none.

**Should fix**
1. `desktop__tour-step2.png`, C_p view: the orange "S: Cₚ = 1" label runs into "Cₚ = 0 at θ = 113.2°" at the top
   left. Move the zero-angle label down (e.g. `P.px.y + 26`) or right-align it past the marker.
2. KaTeX console warnings at every size ("No character metrics for '✚' / '²' / '³'"). These come from Unicode inside
   `\text{}`: EQS `src`/`vtx`/`pair`/`wall` live strings (`\text{m²/s}`, `\text{m³/s}`) and D06 step 9 live
   (`\text{at ✚: }`). Use `\text{m}^{2}\text{/s}` and write "at the probe" (knowledge lesson ch01).
3. The Explore callout says "move the sink twice as far away (drag it)", but dragging a marker sets ε, which moves source
   *and* sink. Reword: "drag the sink to double the gap 2ε".
4. The lifting-cylinder kit uses counterclockwise Γ, but its status shows no conversion. The Part B common rule asks for
   both conventions: add "· Γ_ccw = −6.28 ⇔ Γ_cw = 6.28 m²/s" to the lift-kit status.
5. `phone__explore.png` has 16 Explore pages of one control each at 360×640. Mark `ms`, `d` and `G` `optional` when they
   are not in the current kit (they are already hidden per kit), or accept this.

**What works** (keep): the ψ-bars that add tip-to-tail to the black arrow, with a click that isolates one element. Teal
source-fluid tracers fill the half-body and never leave it, which is the mass-balance argument of D07 step 9 made
visible. The "Squeeze ▶ (2mε fixed)" button animates the doublet limit onto grey ghost circles. The D05 step 3 tangent
and normal arrows on the body carry the live n·u ≈ 0.

## cylinder_circulation_lift (E2) — FAIL
**Must fix**
1. **Derivation text renders as TeX errors** (`desktop__derive-d1p4.png`: red "u_r=♀rac1rpartial_ hetapsi";
   `desktop__derive-d3p8.png`: "L = hoUGamma (6.40)"). Two JS strings have single backslashes, so `\f`, `\t` and `\r`
   become form-feed, tab and carriage-return characters:
   - line 2764 (D09 step 4, *why*): `'The polar components $u_r=\frac1r\partial_\theta\psi$ and $u_\theta=-\partial_r\psi$ …'`
     → `$u_r=\\frac1r\\partial_\\theta\\psi$ and $u_\\theta=-\\partial_r\\psi$`.
   - line 2851 (D11 step 8, *why*): `'… the lift theorem $L=\rho U\Gamma$ (6.40).'` → `$L=\\rho U\\Gamma$`.

   The second one is also a bare equation reference in effect, because the equation next to "(6.40)" does not render.
   Suggest adding a lint rule for `(?<!\\)\\(frac|theta|rho|Gamma|partial|psi|…)` inside `'…'` strings (my scan found
   only these two cases in ch06).

**Should fix**
1. The `forces` view is not `hidePortrait` (it shows on phones as the D11 view, which is fine), but Explain §0 says
   "(hidden on phones)". Align the hint text with the layout.

**What works** (keep): term bars for the four pieces of the expanded (6.39), with ◇ ρUΓ on the total and the cross term
highlighted in amber in C_p(θ) (shaded against a dotted "without it" curve). The Γ-keeping logic when U or a is dragged
makes "L stays put when a changes" discoverable. The status carries both conventions at every size. The merged and free
stagnation regimes use exact-text selftest rows.

## vortex_wall_images (E3) — PASS
**Must fix**: none.

**Should fix**
1. Tour step 1 has `play: true` on the question step, so the trace is already running before the reader has read the
   question (`phone__tour-step1.png` at t = −38 s). Consider `play: false` and let step 3 start the motion.
2. The phone status drops the units: "Γ_cw 1 ⇔ Γ_ccw −1". Add m²/s if it fits.

**What works** (keep): three linked views on one clock. The sensor trace shows the faint whole curve, the bold part so
far, and both parts dashed in their colours with markers at 0, ±4πh²/Γ and ±4√3πh²/Γ. The wall and fluid are washed by
p − p∞. The "same-sign image" toggle shows the leaking wall in rose (the ch05 "deliberate miss" pattern). The source
mode checks (6.41)'s hyperbola ψ = m/4 with a number.

## complex_potential_corners (E4) — PASS
**Must fix**: none.

**Should fix**
1. **Teaching, tour step 2** (`desktop__tour-step2.png`): with w = z² the central quotients along x and along iy are
   exact at every h, so the "two quotients" view shows two coincident flat lines and nothing converges. The idea the
   explainer is named for (two directions converge to one number only for analytic w) is only visible in the D14 pages,
   which use z³. Set step 2 to `n: 3, sig: 0.4` and say "they close in from above and below as h shrinks", or add a
   `sig` animation.
2. `phone-tall__tour-step4.png`, CR bars: the "∂ψ/∂y" label is overlapped by its bar ("∂ψ/∂ʏ −1"). Put the value on the
   empty side or shorten the bar area on narrow views.

**What works** (keep): purple dw/dz and teal u + iv arrows as mirror images. The branch cut is drawn and the status
reports when the step crosses it ("⚠️ the step crosses the branch cut"). The log–log speed-near-the-tip view has a
slope readout, with the 270° corner of Example 6.2 as a preset. z* serves as the non-flow contrast.

## blasius_kutta_contour (E5) — PASS
**Must fix**: none.

**Should fix**
1. The status never states the Γ convention, although the G slider is "clockwise". Add
   "· Γ_cw 2 ⇔ Γ_ccw −2 m²/s" (Part B common rule, as in E2, E3 and E6).
2. Quiz 4 answer: "n = 32 gives 0.003 N/m and n = 64 0.00002 N/m" reads as if the lift were 0.003 N/m. Write "off by
   0.003 N/m … 0.00002 N/m".
3. `phone__explore.png`: 15 Explore pages at 360×640. Mark `nq` and the offsets hidden on phones, or accept this.

**What works** (keep; best-in-chapter candidate for `knowledge/viz_patterns.md`): a "faint = size on the circle,
solid = what survives the trip round" bar per power of z, which makes the residue argument visible. The force-vs-R view
is flat wherever the contour encloses the body, with a rose band where it cuts it. Teal per-arc shares change while the
amber total stays pinned. The tilted-stream mode handles G7 correctly (x/y parts vs F⊥/F∥). The (6.61) correction is
shown with live numbers.

## conformal_joukowski (E6) — FAIL
**Must fix**
1. **`tools/eq_refs.py` flags 1 string** (the only hit in the chapter): `<meta name="viz:concept">` (line 12) reads
   "angle preservation (6.63)–(6.64), the Zhukhovsky map (6.65), … (6.68) … (6.69) (C11)". Numbers are given without
   their equations, and this text is published on the gallery card. Drop the numbers ("Conformal mapping: angles kept,
   the Zhukhovsky map z = ζ + b²/ζ, circle → slit and ellipse, the flow round an elliptic cylinder and the outside-root
   inverse (C11)") or write each equation next to its number. Convention 10 requires `eq_refs → 0`.

**Should fix**
1. Tracer dots are not clipped to the plot rectangle. In `desktop__tour-step5.png`, dots sit below the ζ-plane axes
   among the tick labels and above the plot top, and in `phone-tall__tour-step3.png` they spill into the left and right
   gutters. Draw the particles inside `P.clip(...)`.
2. The D21 step 5 *why* ends with "Our check, not the book's". Also add that numpy's root is fine for Re z > 0, so the
   reader does not overgeneralise (the Explain §4 line already says this).

**What works** (keep): the cross and its image carry a live β = α readout that doubles at the critical point. A ghost
second root b²/ζ₁ appears inside the circle. Rose streamlines computed through numpy's root show the branch problem
directly. The same particles move in both planes. The "Stretch round the circle" plot shows |dz/dζ| together with the
surface q/U.

## laplace_relaxation (E7) — PASS
**Must fix**: none.

**Should fix**
1. D23 step 7 uses a different tolerance from every other number on the page: "k ≈ ln(10⁻⁸)/ln ρ: 27 (Jacobi), 14
   (Gauss–Seidel)", while the tour, the status and the check use 10⁻¹⁰ from r₀ = 3 (34 and 19 sweeps). Use
   ln(10⁻¹⁰/r₀)/ln ρ ≈ 35 / 18 to match the bars' ▼ markers, or say explicitly that 10⁻⁸ is a reduction factor.
2. `phone-tall__tour-step4.png`, residual view: the "10⁻¹⁰ stop" label collides with "SOR 12" at the right edge. Move
   one of them (e.g. put the stop label at the left end of the dashed line).
3. Contraction grid, mid-sweep (`desktop__tour-step6.png`): the unknowns at (3, 1) and (4, 1) show no number while
   their neighbours do. Print "0" (or "·") for not-yet-updated nodes so the sweep front reads cleanly.

**What works** (keep): node-by-node stepping with the purple stencil, and the four-equation matrix with columns
"ψ now | b | exact" and b − Aψ. Faint ghosts of the other two methods sit on the residual plot, with a
predicted-vs-measured ρ table. The refinement study honestly reports order 4/3 and gives the r^{2/3} corner reason
(G2).

## axial_singularity_bodies (E8) — PASS
**Must fix**: none.

**Should fix**
1. `desktop__tour-step4.png`, body view: short black spurs of the computed ψ = 0 contour poke below the axis at the
   nose and tail (around z ≈ ±1.3 m). Clip the computed contour to R ≥ 0, or stop the marching squares at the axis
   row.
2. The Explain §5 hint says odd N makes A "exactly singular … the fit fails". fluidpy now detects this case and falls
   back to a minimum-norm least-squares solve with a RuntimeWarning. Say "LU fails; fluidpy warns and uses least
   squares — use even N" so the text matches the library.

**What works** (keep): sums converge while the bars do not (Σkₙ∆ξ and m₁ against the exact moment of each target, with
the alternating bars shown). The fit-vs-N view has twin axes with a "cond > 10¹²: noise" band. The inspector traces
one row of A with its two largest terms. The panel mode shows the circle's round-off exactness beside the ellipse's
slope −2 (G3c).

## added_mass_sphere (E9) — PASS
**Must fix**: none.

**Should fix**
1. `phone__explore.png` (steel preset): the "duₛ/dt" label at the bottom of the sphere view loses its subscript to the
   view edge ("du./dt"). Keep labels ≥ 12 px inside the view, or place the label beside the arrow tip.
2. Tour step 6 quotes numbers (2g) with `play: true`. Consider starting paused, with a "press ▶" line.
3. Explain §0 says "rim coloured … (red above, blue below)". The chapter's pressure colours elsewhere are orange/blue.
   Fine as is, but name the colours exactly as drawn.

**What works** (keep): the teal speed part and orange acceleration part add up to the black total in the θ view, and
the force bars show the speed part at 0. Added mass comes out three ways (pressure, energy, formula). Release mode has a
with/without-M comparison and a density table for bubble, wood, water and steel. A "sideways" preset shows the force
following the acceleration, not the velocity. The (6.104) printed-bracket number is live in D29 step 8.

---

## Round 2 — 2026-09-24 (E2, E6, E8 re-audited)

Re-run on the fixed files: `tools/viz_lint.py --chapter ch06` → 9 × ok; `tools/shot.py` on the three files → PASS
(E2 8 sizes / 328 views / 32 selftest rows, E6 8 / 344 / 30, E8 8 / 264 / 33; every parity row ok, no overflow, no
small text, KaTeX present); `tools/eq_refs.py --chapter ch06` → **0** strings; the single-backslash TeX scan over all
nine chapter scripts → 0 hits.

| slug | round-1 Must fix | round 2 | verdict |
|---|---|---|---|
| cylinder_circulation_lift (E2) | 2 broken *why* strings | `desktop__derive-d1p4.png` now renders $u_r=rac1r\partial_	heta\psi$, $u_	heta=-\partial_r\psi$; `desktop__derive-d3p8.png` renders $L=ho U\Gamma$ (6.40). Fixed. | **PASS** |
| conformal_joukowski (E6) | bare numbers in `viz:concept` | meta now "Conformal mapping: angles kept, the Zhukhovsky map z = ζ + b²/ζ, …" — eq_refs 0. Also done: tracers clipped to the plot (`desktop__tour-step5.png`: no dots in the tick area), D21 step 5 *why* now says numpy's root is right for Re z > 0. | **PASS** |
| axial_singularity_bodies (E8) | none | Explain §5 hint and the N help now say fluidpy falls back to the minimum-norm least-squares solve for odd N (matches `core.panels`). Spurs: `desktop__tour-step4.png` still shows the computed ψ = 0 line curling onto the axis at both ends (a very short hook near z ≈ ±1.3 m) — barely visible now; stays a Should fix. | **PASS** |

The remaining Should-fix items of round 1 (Γ conventions in E1/E5 status, E4 step 2 with z³, E7 D23 step 7 tolerance,
label overlaps in E1/E4/E7/E9, E1 KaTeX glyph warnings) are open but none blocks the gate.

## Verdict: PASS (all explainers PASS — round 2)

Round 1 was FAIL (E2 and E6).


Round-1 note: 7 of 9 explainers PASSED. **E2 `cylinder_circulation_lift`** fails because two derivation *why* strings render as TeX
errors (single backslashes at lines 2764 and 2851). **E6 `conformal_joukowski`** fails on its one `eq_refs` hit
(bare equation numbers in `viz:concept`). Both are one-line fixes. After them, re-run
`tools/viz_lint.py --chapter ch06`, `tools/shot.py viz/ch06/cylinder_circulation_lift.html viz/ch06/conformal_joukowski.html`
and `tools/eq_refs.py --chapter ch06` (expect 0). The Should-fix items can go in the same round.
