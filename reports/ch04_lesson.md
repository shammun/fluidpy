# Chapter 4 — lesson review (round 1)                              2026-09-23

Reviewed: `outputs/ch04/executed.ipynb` (599 cells, 0 execution errors), built by `notebooks/build_ch04.py`. I read the
whole text dump top to bottom (`outputs/ch04/nb_dump.txt`), looked at every saved figure, rendered p. 117 (Eqs. 4.43–4.45),
and re-ran the doubtful numbers in scratch files (`outputs/ch04/check_lo.py`, `outputs/ch04/ngram_check.py`).

coverage_check: **OK** (0 errors, 7 warnings: all are Derivation tabs of explainers not built yet) · sections covered
**11/11** · CORE blocks **15/15** with code + visual · derivations **30/30** (235 steps) · 12 recaps, 23 new primers
(P111–P133) · the paraphrase check found 9-word overlaps with the book in only 2 cells (one is a list of symbols, the other
the stock phrase "can be written as the gradient of a scalar"): no copied prose.

Checks the orchestrator asked for:
- **Coriolis signs.** (4.43) has +2Ω × u′ on the acceleration side and (4.45) has −2Ω × u′ on the force side. Both match
  the rendered page, and cell 314 explains the difference between the two.
- **Rossby number.** Ro = U/(2Ωl) is correct wherever it appears (cells 593, 595).
- **Jet split.** Q(1 ± cos θ)/2 is not in the notebook; it lives only in the explainer, so there is nothing to check here.
- **Book typos.** All 14 slips listed in analysis §9 are taught in corrected form, and so are the ambiguities in that
  section: γ = μ, "deviatoric", p versus p̄, the pitot tube's "irrotational", and the two meanings of (4.9).
- **Gravity value.** g is inconsistent inside the notebook (Should-fix 3). Cell 457 uses 9.80665 without saying so.
  Every other hydrostatic number (4905 N, 9810 N/m³, 981 Pa, 11.8 km) uses 9.81.

## Must fix

1. **Cells 288–289 · the Lamb–Oseen reading notes are physically false.** "Arrows pointing backwards inside the core and
   forwards outside" and "Viscosity slows the fast core and speeds the slow outskirts" are both wrong. For u_θ =
   Γ/(2πr)(1 − e^{−r²/σ²}), the viscous force is μ ω′(r) e_θ, and it opposes the swirl at every radius. I checked with
   `viscous_force_forms`: the θ-component is −7.7e-3, −1.56e-2, −1.47e-2, −6.3e-3, −1.5e-3 N/m³ at r = 0.2, 0.5, 1, 1.5,
   2 m.
   - **Change:** "Every arrow points against the swirl, strongest near r ≈ σ. Viscosity slows the vortex at every
     radius, while its vorticity spreads outward (ω rises outside the core, σ² = 4ν(t + t₀))."
   - Keep the contrast with solid-body rotation, where no force acts.
2. **Cell 255 · wrong scaling for body couples.** The note says "a torque per unit *mass* scales like h⁵". A couple per
   unit mass times the cube's mass ρh³ scales like **h³**. That is why it can balance (τ₁₂ − τ₂₁)h³ at every size, and it
   is what D08 step 8 itself says.
   - **Change:** h⁵ → h³.
3. **Cells 531–532 · the heat-flux arrows in the composite-wall figure point the wrong way.** In panel (a),
   `annotate("", (24.5, yv − 3), (24.5, yv + 3))` draws downward arrows, from the 20 °C face toward the 30 °C face. Heat
   flows from the 30 °C wall at y = 0 **upward**.
   - **Change:** swap the arrow's head and tail.
   - Add to the reading notes: "heat flows up, from hot to cold, at the same 9.09 W/m² in both layers".
4. **Derivation checks cite computations that the notebook never runs.** The ch03 round-1 review treated this as a
   Must-fix. Either add the cell, or reword the check to name a real cell:
   - **D13 (cell 274).** "the three forms agree to stencil accuracy (`viscous_force_forms`)": cell 288 discards two of
     the three forms (`lap, _, _ =`). Add `assert np.allclose(lap, f2) and np.allclose(lap, f3)`; I checked that they
     agree.
   - **D21 (cell 370).** "`ch04.kinetic_energy_budget` residual ≈ 0 on a test field (C10)": this is never called.
   - **D24 (cell 408).** "`ch04.lamb_identity_sym` returns zero vectors": this is never called (cell 425 checks the
     identity numerically only).
   - **D28 (cell 474).** "sympy (C13 note N110): e = C_vT, p = ρRT, α = 1/T ⇒ ρC_v DT/Dt + p∇·u = ρC_p DT/Dt": N110
     (cell 476) is markdown only, and no sympy cell exists.
   - **D12 (cell 269).** "The toy field u = (ax², 0, 0) has … a nonzero (μ_v + ⅓μ)·2a term": this is not computed
     anywhere.
5. **Cell 54 · the text disagrees with the printed number.** The text says "local changes of ±0.2 kg/(m³ s)", but cell
   53 prints min −0.153 and max +0.153.
   - **Change:** "±0.15 kg/(m³ s)".

## Should fix

1. **All nine 🎮 cells lack a "why interactive" paragraph.** Add one sentence to each on what the interaction shows that
   a static figure cannot. Each cell already has good **What to try** bullets.
2. **Explainer IDs used as jargon.** "E1" (cell 121, D05 check), "E3 parity row" (cell 236, D10 check), "E4's wrapper"
   (cell 281) and "E7's oil bearing" (cell 475) are never explained. Name the explainer instead, e.g. `control_volume_budgets`.
3. **Inconsistent g.**
   - Cell 457 calls `statics.integrate_hydrostatic` with the Ch. 1 default G0 = 9.80665. The printed 352 595 Pa at 25 m
     would be 352 681 Pa with 9.81. Pass `g=9.81` or add a comment.
   - Cells 332–333 and 340 use a spherical-Earth toy with g_n = 9.80, which gives |g_e| = 9.783 m/s² at 45°. A reader who
     knows standard gravity (9.80665 at 45°) will be puzzled. Say that it is a toy with uniform g_n.
4. **Three different "lakes".**
   - The C13 question (cell 458) has 2 °C, less than 10⁻³.
   - The tiny example (cell 481) has αδT = 2×10⁻³, g′ = 0.0196.
   - The code's `"lake"` scenario (cell 483) prints αδT = 7.5×10⁻⁴, g′ = 0.00736, heating 2.6×10⁻¹³, versus the tiny
     example's 2×10⁻¹¹.
   - Either align them, or say that the tiny example's lake is the code's `"thermocline"` row.
5. **Terms used before they are introduced.** Add a forward pointer or a gloss at the first use:
   - Poiseuille numbers in the D11 check (cell 263); the flow is introduced at cell 280.
   - The Rankine core in the D24 check (cell 408); it is worked in cell 422.
   - ζ and f in D14 "What it means" (cell 304); f is primed only at cell 317, and ζ is never defined.
   - g_n and Φ_n in D18 step 5 (cell 325); they are explained only in N64 (cell 327).
6. **D06 check (cell 166).** It says the cell prints `U*dU + g*dz + dp/rho`. It actually prints
   `U*U_s + g*sin(theta) + p_s/rho` and a residual of 0.
7. **D20 check (cell 362).** It says `stress_work_split` runs "on a random field". Cell 365 uses the Taylor–Green field.
8. **D04 (cell 82).**
   - The Tools line says the curl-of-a-product gloss has "the sympy check in N14's cell". That cell checks
     ∇·(∇χ×∇ψ) = 0, not the product rule.
   - Step 4's *why* ("the minus sign is the usual convention") should say what the sign buys: it gives u = +∂ψ/∂y.
9. **Cell 35.** `MassBudget(... local=nan ...)` is printed with no explanation. Say that `local` is not computed for a
   moving box, or hide it.
10. **Cell 240 formatting.** "…flow `N47` — > ⚠️ **Mean vs thermodynamic pressure.**" renders a literal ">". Start the
    callout on a new line.
11. **Cells 551 and 552 repeat the same book-slip callout** ((4.100) → (4.101)). Keep one.
12. **D26 step 3 heading (cell 435).** "Insert into the steady-less (4.69), $…$" puts a full equation in a bold heading
    and uses an odd word. Make it "Insert into (4.69) with ω = 0" and move the equation to the *why* line.
13. **Cell 253.** "every plane feels the same extra pull 3μ_v" is missing a unit. Write μ_v S_mm = μ_v × 3 s⁻¹ (Pa).
14. **Cube-spin figure (cells 254–255).** The grey band "fastest spin-up plausible for a fluid element ~10⁴ rad/s²" is an
    unexplained, arbitrary number. Remove it, or justify it in the reading notes.
15. **Wake figure (cell 136).** The rose "−F_D/l on the fluid" arrow is drawn diagonally, but the force is along −x. Make
    the arrow horizontal.
16. **Cylinder ψ figure (cell 90).** White pixel notches show at the body edge, where the pcolormesh mask and the drawn
    circle differ. Draw the circle slightly larger, or mask at r < 1.02.
17. **Boussinesq validity figure (cell 488).**
    - The legend covers the "lake" bars.
    - Colouring every above-threshold bar amber hides which quantity crossed. Keep the colours and mark the crossing
      with an outline or a hatch instead.
18. **Cells 526–527.** The cap forces are computed with Δp = 146 Pa, not the exact 145.6 Pa. Most of the printed
    mismatch comes from that rounding, not from "higher order in ζ". Pass `2*0.0728/1e-3`.
19. **N156 (cell 580) and cell 581.** The wave-drag factor is written "λ³", but λ = 1/25 in the tiny example. It is
    λ⁻³ (or define λ = l_p/l_m).
20. **Equation numbers wrapped in italics** around the LaTeX (cells 0, 398, 598), e.g. `*(4.65), $…$*`. Cell 398 has
    nested asterisks that break the italics. Use the plain "(N.M), $…$" form used elsewhere.
21. **Several from-scratch cells have no "What does the code above do?"** (cells 67, 202, 283, 425, 485). One line each
    would do.
22. **Stray empty stream outputs** (cells 167, 194, 334). This is cosmetic.

## Unexplained-first-use list

Rows explained at or before first use: all 23 new primers (P111–P133) and all 12 recaps sit before their first use. The
"Tools from earlier chapters" blocks name the chapter and primer for every tool the chapter reuses. Python idioms that
appear for the first time here (`sp.lambdify`, `np.ma.masked_where`, `contourpy`, `make_subplots`,
`sp.KroneckerDelta`, `sp.solve`/`.coeff`, `np.outer`) are each explained by their line comment. That is adequate, though
there is no ch03-style idiom gloss.

These rows are missing or late:

| term | first cell | explained at |
|---|---|---|
| ζ (relative vorticity) | 304 (D14 "What it means") | **missing** |
| f = 2Ω sin φ | 304 | 317 (P126): **late** |
| g_n, Φ_n (gravitation without the centrifugal part) | 325 (D18 step 5) | 327 (N64): **late** |
| Poiseuille flow and its numbers | 263 (D11 check) | 280: **late** |
| Rankine-core numbers (B = r² − 1) | 408 (D24 check) | 422: **late** (the vortex itself is Ch. 3) |
| "E1", "E3", "E4", "E7" | 121, 236, 281, 475 | **missing** |
| Stokes' first problem | 262 (idea sketch) | 278: **late**, but only named in a sketch (acceptable) |

## Derivation audit

| D | steps | every step follows? | why / in-words ok | check ok | verdict |
|---|---|---|---|---|---|
| D01 | 9 | yes | yes | yes (3-box numbers run) | pass |
| D02 | 8 | yes (localisation argument written out) | yes | yes | pass |
| D03 | 5 | yes | yes | yes | pass |
| D04 | 8 | yes (normal to the right, sign checked) | step 4's why is thin; Tools misattributes a sympy check | yes | pass (Should 8) |
| D05 | 8 | yes | yes | "E1" jargon | pass (Should 2) |
| D06 | 11 | yes (steps 5–8 expanded by hand and by sympy) | yes | printed text misquoted | pass (Should 6) |
| D07 | 10 | yes | yes | yes | pass |
| D08 | 8 | yes (torques, I = ρh⁵/6, h⁻² checked) | yes | yes | pass; the linked cell 255 fails (Must 2) |
| D09 ★★★ | 13 | yes | yes | sympy rebuilds K, contracts it, rotates it | pass |
| D10 | 4 | yes | yes | "E3 parity row" jargon | pass (Should 2) |
| D11 | 7 | yes | yes | Poiseuille used before it is introduced | pass (Should 5) |
| D12 | 7 | yes (μ + μ_v − ⅔μ = μ_v + ⅓μ) | yes | toy-field claim not computed | **fix (Must 4)** |
| D13 | 9 | yes (ε_jik = ε_kji, δ substitution checked) | yes | "three forms agree" not shown | **fix (Must 4)** |
| D14 | 7 | yes | yes, but ζ and f undefined | yes | pass (Should 5) |
| D15 ★★★ | 12 | yes | yes | sympy rebuilds x = X + Rx′ and tests the factor 2 | pass |
| D16 | 6 | yes; matches p. 117 | yes | yes | pass |
| D17 | 5 | yes (9451 m against the exact 9342 m) | yes | yes | pass |
| D18 | 5 | yes | g_n undefined | yes (0.0339 m/s²) | pass (Should 5) |
| D19 | 8 | yes | yes | yes | pass |
| D20 | 7 | yes | yes | "random field" is really Taylor–Green | pass (Should 7) |
| D21 | 6 | yes | yes | `kinetic_energy_budget` not run | **fix (Must 4)** |
| D22 | 8 | yes (force work cancels; the ½R slip noted) | yes | yes | pass |
| D23 | 8 | yes (δ_ijδ_ij = 3) | yes | yes (1000 random gradients) | pass |
| D24 | 10 | yes (ε_jkl = ε_ljk, δ substitution) | yes | `lamb_identity_sym` not run | **fix (Must 4)** |
| D25 | 5 | yes | yes | yes | pass |
| D26 | 8 | yes (sign of the absorbed function fixed) | step 3 heading | yes (sympy: 0 and 2B) | pass (Should 12) |
| D27 | 9 | yes | yes | yes | pass |
| D28 | 8 | yes (−pα DT/Dt = −ρR DT/Dt) | yes | cited sympy cell does not exist | **fix (Must 4)** |
| D29 | 8 | yes | yes | yes | pass |
| D30 | 8 | yes | yes | yes (sympy brackets) | pass |

## What works well

Keep these, and consider them for `knowledge/` and the teaching-style skill:

- **Every book slip is taught as a correction with a reason.** Examples: the "= 0" of (4.15), the φ-gauge sign of
  (4.74), the index slip after (4.24), and "κ" meaning μ_v with ≥ rather than >. Each gets a code or sympy
  demonstration (`gauge_absorbed_bracket` gives 0 against 2B).
- **The equation ledger** (pandas table 0 → 6/13 → 4/5 → 5/5 → 7/7) threads §4.1 through C06, C08 and C10. Candidate
  pattern: "running count of equations against unknowns".
- **The "four primes" and "reused symbols" callouts** up front, plus the code rule "one name per meaning" (`u_rot`,
  `p_pert`), deal with Kundu's reuse of symbols head-on.
- **Coriolis term versus force** (cell 314), with `frame_acceleration_terms` against `apparent_body_forces`: exactly the
  sign language Ch. 13 needs.
- **Both ★★★ sympy checks rebuild the derivation's own construction.** D09 builds K from Kronecker deltas and rotates it;
  D15 builds x = X + Rx′ and shows the factor 2 is necessary. This follows the ch01 lesson.
- **Two routes to every key number:** the budget by hand against the library, the D06 stream-tube expansion against the
  library, and the stress law by explicit loops, by `einsum` and by the library.
- **The lapse-rate convention is shown both ways** in N137 and N114 (Kundu's Γ and Γ_met), following the ch01 rule.
- **The climate hooks are real:** geostrophy and Rossby number, geopotential, Boussinesq ocean and atmosphere,
  thermal-wind baroclinicity, and highs and lows.

## Verdict: FAIL

Five Must-fix items, all small edits to `notebooks/build_ch04.py`:
- rewrite one false reading note (Lamb–Oseen);
- correct one scaling (h⁵ → h³);
- flip one pair of figure arrows (composite wall);
- make five derivation checks point at computations that exist;
- fix one number (±0.2 → ±0.15).

After those, and a grep for any remaining "E[0-9]" labels, the lesson should pass. The only re-check needed is on the
changed cells.

## Round 2 (2026-09-23)

Re-reviewed: the regenerated `outputs/ch04/executed.ipynb` (611 cells; 113 s run with 0 execution errors). I dumped it
again and diffed it against round 1 line by line (`outputs/ch04/nb_diff.txt`, 389 diff lines). I re-read every changed
cell and every new cell, and looked again at the regenerated figures (interfaces, wake budget, Boussinesq validity).

coverage_check: **OK** (0 errors, 2 warnings). The warnings are the Derivation tabs of `kinematic_free_surface` (the
backup explainer) and `dynamic_similarity_models` (still being built). Sections 11/11 · CORE 15/15 with code + visual ·
derivations 30/30.

### Round-1 Must-fix: all resolved

| # | item | status |
|---|---|---|
| 1 | Lamb–Oseen reading notes | Fixed. A new printed line shows the θ-component of the viscous force at r = 0.2–2 m: −0.0077, −0.0156, −0.0147, −0.0063, −0.0015 N/m³, the same as my check. The notes now say every arrow opposes the swirl, largest near r ≈ σ, while the vorticity spreads outward. |
| 2 | Body-couple scaling | Fixed. The text now says "couple per unit mass × ρh³ scales like h³, like (τ₁₂ − τ₂₁)h³", pointing to D08 step 8. |
| 3 | Heat-flux arrows | Fixed. The arrows in the figure now point upward, from the 30 °C wall to the 20 °C wall. The reading notes say "heat flows up, from hot to cold, 9.09 W/m² in both layers". |
| 4 | Checks that cited missing computations | Fixed. Each check is now a real cell (details below). |
| 5 | ±0.2 against a printed ±0.153 | Fixed. The text now says "±0.15 kg/(m³ s) (the printed min and max)". |

The checks from Must-fix 4:
- **D12:** a new sympy cell prints ∇·u = 0 and [0, 0, 0] for the Poiseuille-like field, and 2a(μ + 3μ_v)/3 = 2a(μ_v + ⅓μ) for the toy field. Both are correct.
- **D13:** the viscous-paradox cell now asserts that the three forms of (4.40) agree.
- **D21:** a new cell runs `kinetic_energy_budget` on Taylor–Green. It prints −240.0047 = −239.1457 − 0.8590 with residual 0, which checks by hand.
- **D24:** a new cell runs `lamb_identity_sym` on a general 3-D field and prints [0, 0, 0].
- **D28:** a new sympy cell prints 0. I checked its construction: p∇·u = ρRT · (1/T) DT/Dt = ρR DT/Dt, which turns ρC_v into ρC_p.

### Round-1 Should-fix: 21 of 22 done

The builder did all of these, and I checked each one in the diff:
- There is now a "Why interactive" paragraph in all 9 🎮 cells.
- The E1/E3/E4/E7 labels are replaced by explainer names; a grep finds none left.
- g is consistent: `integrate_hydrostatic(..., g=9.81)` now prints 352 681 Pa at 25 m. The toy g_n = 9.80 is now explained in two places, next to standard gravity (9.80665 at 45°).
- The lake numbers are reconciled: the tiny example is the code's `"thermocline"` row, and the question's lake and the code's `"lake"` row are each explained.
- Forward pointers or glosses were added for Poiseuille (D11), Rankine (D24), ζ and f (D14), and g_n and Φ_n (D18). The three example flows of C08 are named at their first mention.
- The D06 check now quotes the printed expression. The D20 check now says the Taylor–Green field. D04's Tools line and step 4's *why* are corrected.
- `local = nan` is explained, and the N47 callout renders correctly.
- The duplicate slip callout is removed, and the heading of D26 step 3 is fixed.
- The μ_v S_mm unit is added, and the cube-spin band is removed.
- The wake arrow is now horizontal and the cylinder no longer shows notches.
- In the Boussinesq figure, bars over the threshold are now hatched with their colour kept, and the legend sits outside the bars.
- The cap forces use the exact Δp.
- The ship scaling now reads λ⁻³ = 25³ and Re_p/Re_m = (l_p/l_m)^{3/2}.
- The nested italics around equation numbers are removed.
- "What does the code above do?" was added to the from-scratch cells.

Left open: #22, the stray empty stream outputs. There are now about 4 of them (cosmetic, from the kernel splitting output).

### New cells checked for correctness and first-use order

The four new check cells and their explanations use only tools primed or recapped earlier: sympy (P40), and
`viscous_force_forms`, `kinetic_energy_budget` and `lamb_identity_sym`, each with a line comment. Every number they print
matches their prose. The new "Why interactive" paragraphs and pointers introduce no new terms.

The unexplained-first-use rows from round 1 are all resolved:
- ζ and f: glossed at cell 304, with a pointer to P126.
- g_n and Φ_n: defined in D18 step 5.
- Poiseuille and Rankine: forward pointers added.
- The E-labels: replaced by explainer names.

No new missing rows.

Cosmetic, no action needed: in the Boussinesq figure the "threshold 0.1" label touches the legend.

**Verdict: PASS**
