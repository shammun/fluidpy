# Chapter 5 — lesson review (round 1)                              2026-09-23

Reviewed: `outputs/ch05/executed.ipynb` (487 cells, 0 execution errors), built by `notebooks/build_ch05.py`. I read the
whole text dump top to bottom (`outputs/ch05/review_dump.txt`) and looked at every saved figure
(`outputs/ch05/review_imgs/`). I also re-ran the doubtful numbers: the planetary term of cell 394, and R·T·ln 2 in cell 102.

coverage_check: **OK** (0 errors, 7 warnings; all are Derivation tabs of the three explainers still being built) ·
sections covered **8/8** (§5.1–5.8) · CORE blocks **14/14** with code + visual · derivations **23/23** (194 steps) ·
18 recaps, 15 new primers (P134–P148). The embeds `vorticity_equation_rotating`, `point_vortex_lab` and
`vortex_sheet_rollup` (cells 405, 458, 482) return stub HTML because the explainers are still being built. That is not
a notebook defect.

Checks the orchestrator asked for:
- **Book slips.** All of them are taught in corrected form:
  - (5.14) +1/(4π): D10 step 11, with both signs computed in cells 249 and 255.
  - The integrand's second sign slip: D11 step 5, sympy cell 258.
  - Fig. 5.2 "2ω" → ω/2: R02 and D02.
  - Kelvin needs ρ = ρ(p), not "single valued": D05 step 6.
  - Fig. 5.11: G is the centre of vorticity, not a stagnation point; the fluid moves at −1.70 m/s there (cell 420).
  - The sheet caption's u₁ − u₂: D23 and cell 467.
  - Exercise references 5.8 → 5.9 and 5.11 → 5.10: D10, D19 and N25.
  - Also Π → Φ, the silent u_{j,j} drop, the ½ in the Lamb step, and "hyperboloids".
- **Baroclinic sense.**
  - Correct everywhere. With ∇p = (0, −P) and ∇ρ = |∇ρ|(sin θ, −cos θ), (∇ρ×∇p)_z = −|∇ρ|P sin θ. So θ = 0 gives zero
    torque, and 0 < θ < 180° spins clockwise (cells 153, 155, 161–162).
  - In the disc figure (cell 159), G sits lower-right and the arrow turns clockwise.
  - Lock exchange with the heavy fluid on the left: +2.42 s⁻², counterclockwise (D07; the arrows in cell 169 turn
    counterclockwise).
- **Planetary stretching and columns.** 2Ω∂w/∂z = 1.458×10⁻⁹ s⁻² is correct. The ridge makes ζ = −0.2f (anticyclonic)
  and the trough +0.2f. The poleward ring gives −4.19×10⁷ m²/s. The Southern-Hemisphere reading is correct.
- **Images.** Correct. The image of Γ > 0 above the wall is clockwise and drives the vortex toward +x at Γ/4πh; the sum
  arrows under the vortex point +x, 0.637 m/s (cells 437, 440, 441).
- **Qualitative models are labelled.**
  - Ring near a wall: `valid` samples only, stop time 3.117 s at a gap of 3 core radii (cells 446, 449–452).
  - Leap-frogging (N43) and sheet roll-up (cell 478) are both called qualitative.

## Must fix

1. **The text quotes numbers the code does not print.**
   - **Cell 128, item 2.** The text says "spread ~10⁻¹⁵". Cell 127 prints `spread 3.0e-13 m²/s`.
     - **Change:** "spread 3×10⁻¹³".
   - **Cell 211, item 1.** The text quotes "advective +15.4937, stretching +61.9750, diffusion −46.4813". Cell 210
     prints `+15.49 +61.97 −46.48` for the preset and `15.4931, 61.9744, −46.48` for the general call. None of the quoted
     digits are printed, and the 4th decimal disagrees with the second call.
     - **Change:** print the preset with `:.4f`, then quote those numbers. Or quote 15.49 / 61.97 / −46.48.
   - **Cell 219, item 2.** The text says "agrees … to a few times 10⁻⁵ of the peak". Cell 218 prints `7.9e-06`.
     - **Change:** "to about 10⁻⁵ (8×10⁻⁶) of the peak".
2. **Cell 394 · the output looks wrong.**
   - `print(ch05.planetary_vorticity_terms(...))` shows `[0. 0. 0.]`, because the notebook's numpy print options round
     2.9e-9 to zero. The comment claims 2.9168e-9, and the true value is `2.916846e-09` (I re-ran it).
   - **Change:** `print([f"{v:.4e}" for v in ...])` and say "= 2Ω·2×10⁻⁵ = 2.9168×10⁻⁹ s⁻²" in the text.
3. **Cell 99 · false physics claim.** "…the laboratory analogue of geostrophic balance (Ch. 13)" is wrong.
   - Water at rest in the rotating frame has **no** Coriolis force, because u′ = 0. Its paraboloid comes from gravity
     plus the centrifugal force: an effective-gravity geopotential, the same reason the Earth is oblate.
   - Geostrophic balance is Coriolis against pressure gradient, and it needs motion relative to the frame.
   - **Change:** "…the centrifugal force folds into an effective gravity (Ch. 4 §4.7); the free surface is a surface of
     constant geopotential — the reason the Earth's equator bulges."
4. **Cell 135 · wrong prediction.** It says: "…ν switched on in the cellular flow: the teal curve would start to fall …
   fastest where the spiral is thinnest".
   - ψ = sin x sin y has ∇²u = −2u, so it is an exact Navier–Stokes solution that decays as e^{−2νt}. For any material
     loop, DΓ/Dt = ν∮∇²u·dx = −2νΓ.
   - So Γ(t) = Γ(0)e^{−2νt} exactly. The rate does not depend on the loop's shape, and "fastest where…" has no meaning
     for the single number Γ.
   - **Change:** "…the teal curve would fall as e^{−2νt} (∇²u = −2u for this flow, so the viscous term of (5.11) is
     −2νΓ), whatever the loop's shape."
5. **The checks of two ★★★ derivations describe tests that the cells do not run.** The ch03 and ch04 reviews treated
   this as a Must-fix.
   - **D09 (cell 206).** The check says: "substitute a divergence-free polynomial field (u = (y z², x z, −x y)·t-free
     plus a Taylor–Green piece)". Cell 207 uses only `[y*z**2, x*z, -x*y]`, with no Taylor–Green piece. It also
     computes no separate curl for each term of (4.39b) beyond pressure and gravity.
     - **Change:** add the Taylor–Green piece, a time-dependent divergence-free field, to `u0`. Or delete that clause
       and describe what the cell prints.
   - **D15 (cell 324).** The check says: "then that the extra terms vanish for a divergence-free polynomial u (steps 7,
     9) … the viscous term equals νω_{n,jj} for that u". Cell 325 uses a generic u throughout and carries u_{j,j}
     symbolically. That test is stronger, but it is not the one described.
     - **Change:** reword the check so it matches the cell: generic u, and the last line isolates u_{j,j}(ω_n + 2Ω_n).
6. **Cells 210, 211, 214, 215 · "Taylor–Green" is used but never explained.** It appears in the budget code, in the
   list "Taylor–Green (2-D)" and as a figure panel title. No earlier chapter defines it (it is not in
   `knowledge/primers.md` or `concept_map.md`).
   - **Change:** add a one-line gloss before cell 210. It is the decaying cellular flow u = (sin x cos y, −cos x sin y)e^{−2νt}
     of C03, an exact 2-D Navier–Stokes solution whose vorticity is only diffused. Name it the viscous version of the
     C03 cellular flow.
7. **Cell 35 · `np.meshgrid(..., indexing="ij")` is used before its recap.** The recap of P76 comes only in cell 243.
   - **Change:** move the P76 recap to the "Tools from earlier chapters" block of cell 31, and mention `indexing="ij"`
     (rows run along the first argument).
8. **Cell 137 · `itg.gauss_legendre_nodes` (40 × 40 Gauss nodes) is used before its primer.** P143 comes only in
   cell 252.
   - **Change:** move the 1-D part of P143 (with its 4-node demo, cell 253) before cell 137. Or add a recap line there
     that points forward to P143.
9. **Cell 161 · `np.deg2rad` is used before its recap.** The recap of P66 comes only in cell 380.
   - **Change:** move that recap line into cell 154's tools block, or add it just before cell 161.
10. **D03 step 8 (cell 80) · a tool is used without a primer.** The step uses the polar vector Laplacian,
    (∇²u)_θ = u_θ″ + u_θ′/r − u_θ/r², citing only "Appendix B's vector Laplacian".
    - The gloss in cell 79 covers the strain rate and the stress divergence, not this.
    - Only the Cartesian vector Laplacian has been primed before (`primers.md`).
    - **Change:** add this formula to the cell-79 gloss, with one sentence on the extra −u_θ/r² (e_θ turns). A 2-line
      sympy demo would help.

## Should fix

1. **Cell 8 · the book-slip table is garbled by auto-inserted equations.**
   - The "(5.14)" row now shows *our* corrected equation in the "The book prints" column, then "as −1/4π…".
   - The exercise row has "(5.33), DΓ_a/Dt = 0 …" wedged inside "after (5.33)".
   - **Change:** write these two rows by hand. Show the printed −1/(4π) form alone in column 1, and put the exercise
     pointers without the equation insert (or with a compact inline form).
2. **Cell 102 · "the same by hand" is not the same.** The library gives 57 334.197 and the hand line gives 57 332.600,
   because it uses R = 287.05 while the library uses 287.058.
   - **Change:** use 287.058, which gives 57 334.198, or say "≈" and give the reason.
3. **Cell 168 · numbers the output cannot show.** The text quotes "2.422219" and "the 3×10⁻⁶ gap", but cell 167 prints
   `[0. 0. 2.4222]`.
   - **Change:** print with `:.6f`. (Cell 337 later prints 2.4222189926.)
4. **D13 result (cell 286) · dangling "=".** The line "u_φ = …; infinite line Γ/2πd =" ends in "=".
   - **Change:** "… infinite line u_φ = Γ/2πd, the line vortex (5.2)".
5. **Cell 280 gloss and D12 step 3 · the variation of the frozen kernel is overstated.** The kernel
   (x − x′)/|x − x′|³ has magnitude 1/d², so moving across a core of radius a changes it by about **2a/d**, not 3a/d:
   20 % at 10 radii and 2 % at 100.
6. **Checks that quote numbers no cell prints.**
   - D13: "semi-infinite 0.079577" is not printed anywhere.
   - D21: "ΣΓx … flat" is not printed.
   - D22: "0 at h = H/2" is not printed, and the check says "200 wall points" where cell 437 uses 201.
   - **Change:** add the one-liners or trim the claims.
7. **Cell 27 · "spectrally accurate on a circle" comes before the reader knows the idea.** The periodic trapezoid
   primer P136 appears only in cell 110.
   - **Change:** add "(the periodic trapezoid rule, primer P136 in C03)".
8. **Earlier-chapter Python tools used with no recap line:**
   - `np.cross`: cell 184 (Ch. 4).
   - `np.linalg.norm`: cell 186 (P67).
   - `np.gradient`: cell 449 (P22).
   - `@` for matrix–vector products: cell 186.
   - **Change:** add them to the nearest tools blocks.
9. **"convention 3" and "convention 4" are undefined** (D07 step 7, D21 assumptions, D23 step 5). The notation table in
   cell 6 is not numbered.
   - **Change:** number the conventions in cell 6, or write the convention out ("ω_z > 0 is counterclockwise").
10. **Derivation tools used before they are stated.**
    - D06 step 6 moves the torque with M_G = M_O + (r_O − r_G) × F, citing P116. P116 defines torque, not the change of
      reference point. Add it to the cell-149 gloss.
    - D05 step 9 uses −∮dp/ρ = ∫(∇ρ×∇p/ρ²)·n dA, which needs the quotient rule ∇(1/ρ) = −∇ρ/ρ². That gloss appears only
      in cell 316. Add a half-line or a forward pointer.
11. **D19 tools (cell 388) · wrong pointer.** It says "vector area of a closed loop (primer in C11)". The primer is P138
    in C05 (cell 183); C11 has only a reminder.
12. **Cell 267 figure · the arrows cannot be read.** The orange du and the purple total are drawn overlapping at the
    field point, and the orange one looks as long as the purple.
    - **Change:** draw them from the point in different directions or lengths, or offset one. Also move the "About 7 %"
      in the notes into the title as a number.
13. **Cell 373 · wrong word.** "the Gaussian's inflexion": the 1-D inflection of e^{−R²/a²} is at a/√2. At R = a it is
    the zero of the *cylindrical* Laplacian.
    - **Change:** "where the cylindrical Laplacian (1/R)(Rω′)′ changes sign".
14. **N43 (cell 453) · overclaim.** "for ever in an ideal fluid" should be "for ever in this thin-core model". Real
    leap-frogging breaks down after a few passes.
15. **Cell 486 · wrong chapter.** "Ch. 11: … baroclinic instability" should be Ch. 13 (§13 baroclinic instability).
    Ch. 11 keeps Kelvin–Helmholtz.
16. **Cell 40 · "on a third of the area".** The ratio is e⁻¹ ≈ 0.37 (0.011557 / 0.031416).
    - **Change:** "on 37 % (e⁻¹) of the area".
17. **Cell 175 (C05) · missing heading.** The *In one line* statement is missing. Add "vortex lines are material lines
    when Kelvin's four restrictions hold".
18. **Cell 215 · "round-off floor" is wrong for Burgers.** Its residual, 2.6×10⁻⁵, is stencil truncation error, not
    round-off.
    - **Change:** "at the stencil-error floor".
19. **D04 check (cell 117) · sign dropped.** "contour term ≈ 3e-14" should be ≈ −3×10⁻¹⁴, as printed and as cell 128
    says.
20. **Cell 221 · "ν ten times larger (air instead of water)".** Air's ν is ≈ 1.5×10⁻⁵ m²/s, about 15 times water's.
    - **Change:** say "15 times … one fifteenth of the time", or drop "(air…)".

## Unexplained-first-use list

| Term / symbol / tool | First used (cell) | Explained at | Status |
|---|---|---|---|
| Taylor–Green flow | 210 (code), 211, 214 | — | **missing** (Must 6) |
| `np.meshgrid(..., indexing="ij")` | 35 | recap P76 at 243 | **too late** (Must 7) |
| `itg.gauss_legendre_nodes` / Gauss–Legendre | 137 | primer P143 at 252 | **too late** (Must 8) |
| `np.deg2rad` | 161 | recap P66 at 380 | **too late** (Must 9) |
| polar vector Laplacian (∇²u)_θ = u″ + u′/r − u/r² | 80 (D03 step 8) | — | **missing** (Must 10) |
| "spectrally accurate" (periodic trapezoid) | 27 | P136 at 110 | late (Should 7) |
| `np.cross` | 184 | Ch. 4 (no recap here) | recap missing (Should 8) |
| `np.linalg.norm`, `@` | 186 | Ch. 2 P67 (no recap here) | recap missing (Should 8) |
| `np.gradient` | 449 | Ch. 1 P22 (no recap here) | recap missing (Should 8) |
| "convention 3", "convention 4" | 165, 415, 465 | — | undefined (Should 9) |
| moment transfer M_G = M_O + (r_O − r_G)×F | 151 (D06 step 6) | — (P116 is torque only) | missing (Should 10) |
| quotient rule ∇(1/ρ) = −∇ρ/ρ² | 122 (D05 step 9) | gloss at 316 | late (Should 10) |
| "cat's eye" | 480 | — | name only (minor) |

Everything else is explained at or before its first use, including these:
- Primers P134–P148, each placed before the derivation that needs it.
- The ε–δ, Schwarz and Lamb recaps.
- `solve_ivp` with many points (P137); FFT differentiation (P136/P142); elliptic integrals with m = k² (P148).
- `itertools.product`, `np.ndindex`, `np.roll`, `np.ptp`, `np.hypot` (explained inline).
- symlog axes (cell 215) and plotly `updatemenus` (cell 363).

## Derivation audit

| D | ★ | Steps | Every step follows? | Why / in words | Check (cell, printed numbers) | Verdict |
|---|---|---|---|---|---|---|
| D01 (5.4) | ★ | 7 | yes | ok | cell 32: flux 0.632121 at z = 0 and 1 m; broken field −0.632 / +0.233 / −0.400 ✓ | pass |
| D02 (5.6) | ★★ | 10 | yes (steps 3–4 fill the book's skipped centripetal move) | ok | cell 71: 2500 Pa/m, 12.742 mm, 6.371 mm, sympy 0 0 ✓ | pass |
| D03 (5.7) | ★★ | 9 | yes; I checked step 8 algebra (2 − 1 − 1 = 0) | step 8 uses the polar vector Laplacian with no primer | cells 83, 87: −0.2342 three ways; torque −0.002 at every r ✓ | **Must 10** |
| D04 (5.9) | ★★ | 8 | yes | ok | cell 127 ✓ (text drops a sign; Should 19) | pass |
| D05 (5.10–5.11, 5.8) | ★★ | 9 | yes; Stokes form of −∮dp/ρ verified | step 9 uses the quotient rule early | cells 127, 137: −3.345e-4; 0.046708396673 twice ✓ | pass (Should 10) |
| D06 (5.28) | ★★ | 10 | yes; I re-derived x_G = R²∇ρ/4ρ₀, M_G = πR⁴∇ρ×∇p/4ρ₀, I = ρ₀πR⁴/2 | step 6 transfer rule not primed | cells 155, 158: −0.0981 both routes; order 2.00 ✓ | pass (Should 10) |
| D07 lock exchange | ★ | 7 | yes; sign ✓ (heavy left → +, counterclockwise) | "convention 3" undefined | cell 167: 2.422222 / 2.4222 ✓ | pass (Should 9) |
| D08 Helmholtz | ★★ | 9 | yes | ok | cells 186, 226: Γ 8.0159e-10 fixed; \|cos\| 1.3e-5; ABC angle 6.7e-13 ✓ | pass |
| D09 (5.13) | ★★★ | 12 | yes | ok | cell 207 prints zeros, but the check text describes a Taylor–Green piece the cell does not have | **Must 5** |
| D10 (5.14 corrected) | ★★★ | 12 | yes; the source strength −4π and G = −1/4πr ✓ | ok | cell 249: −4π, 1, ±0.3088 vs 0.308838 ✓ | pass |
| D11 (5.16) | ★★★ | 12 | yes; both book sign slips shown | ok | cells 258, 262, 264, 270 ✓ | pass |
| D12 (5.17) | ★★ | 6 | yes | the kernel-variation factor is overstated | cell 290: 1.054786, 1.013052, 1.000804 ✓ | pass (Should 5) |
| D13 segment | ★★ | 8 | yes; substitution l = −d cot θ ✓ | result line has a dangling "=" | cell 290: 0.112540, 0.159155, square 0.900316 ✓; semi-infinite not printed | pass (Should 4, 6) |
| D14 (5.25) | ★★ | 8 | yes | ok | cells 320, 322: residuals equal (1.40e-10, 9.00e-11) ✓ | pass |
| D15 (5.30) | ★★★ | 15 | yes; the dropped u_{j,j} term is made explicit | ok | cell 325: every line [0, 0, 0] ✓ — but the text describes a polynomial-u test the cell does not run | **Must 5** |
| D16 (5.31) | ★ | 6 | yes | ok | cell 358 ✓ | pass |
| D17 (5.32) | ★ | 5 | yes | ok | cell 358: 7.3891, (1, 0, 2) ✓ | pass |
| D18 Burgers | ★★ | 9 | yes; I checked the total derivative and C = αΓ/4πν | ok | cells 210, 372: 2 mm, 79.5775, Γ 0.001 ✓ | pass |
| D19 (5.33) | ★★ | 10 | yes; I checked the by-parts step 7 | tools pointer says C11 (should be C05) | cell 396: Γ 0 → 1.985865, Γ_a = π; (−π, 0) ✓ | pass (Should 11) |
| D20 column | ★ | 6 | yes | ok | cell 396: +1.00e-05 ✓ | pass |
| D21 pair | ★★ | 6 | yes; G at −1.70 m/s re-derived | "convention 3" | cell 420: 19.739209 s; I spread 5e-11 ✓; ΣΓx not printed | pass (Should 6, 9) |
| D22 image | ★ | 5 | yes; direction +x ✓ | ok | cells 437, 440 ✓ | pass |
| D23 sheet | ★ | 5 | yes; counterclockwise walk ✓ | "convention 4" | cell 472: L1 0.043968 → order 1.00 ✓ | pass (Should 9) |

## What works well

Keep these; they are candidates for `knowledge/` and the teaching-style skill:
- **Book slips are taught by computing both versions.** Examples:
  - The printed −1/(4π) is run next to the corrected sign on a real tube (cell 255, ±0.308837) and appears as a dashed
    "printed sign" trace in the slider (cell 274).
  - The unequal pair's G is shown to have a fluid velocity of −1.70 m/s.
  - The reader *sees* each slip fail.
- **Three routes to one number.** Stress divergence, Laplacian and −μ∇×ω for the net viscous force (cell 83). Torque
  route against formula for the baroclinic spin-up, with the gap shrinking as R² (cell 155). Loop against area for the
  baroclinic circulation, agreeing to 12 digits (cell 137).
- **The sympy check of D15 keeps the dropped term.** It prints curl(5.25) − (5.30) − u_{j,j}(ω_n + 2Ω_n) ≡ 0 for a
  *generic* u. Pattern: "check what the book dropped by keeping it".
- **The Kelvin decision table** over all 16 hypothesis combinations (cell 140), with the surviving term named in each
  case.
- **Qualitative models are fenced in.** The ring-near-wall model stops at a gap of 3 core radii, with a `valid` mask
  and a `stop_time`, and the dot in the figure marks where it stops.
- **Convention checks put numbers on the screen.** The baroclinic sine curve with the formula ghost (cell 161), the
  hemisphere slider for the column (cell 402), and both sheet conventions computed (cell 472).
- **Every derivation check names a real cell and quotes its output**, with only the few exceptions listed above.

## Round-1 verdict: FAIL

There are 10 Must-fix items. The notebook's structure, coverage and physics are otherwise strong: 8/8 sections, 14/14
CORE blocks, 23 derivations that I checked step by step, and every book slip is corrected. The fixes are small:
- four text or number corrections (Must 1–2);
- two wrong physics sentences (Must 3–4);
- two check descriptions (Must 5);
- four first-use-order and primer gaps (Must 6–10).


---

# Round 2 — 2026-09-23

Reviewed: the new `outputs/ch05/executed.ipynb` (496 cells, 0 execution errors).
- I re-dumped it to `outputs/ch05/review_dump.txt` and diffed it against round 1 (`review_dump_r1.txt`).
- I looked again at the changed figures.
- I re-ran my scan for equation numbers cited without maths.

coverage_check: **OK** (0 errors, 3 warnings: the Derivation tabs of `point_vortex_lab` and `vortex_sheet_rollup`,
not built yet). `vorticity_equation_rotating` is now built and embedded (cell 413). Sections 8/8 · CORE 14/14 ·
derivations 23/23.

## Round-1 Must-fix items: all 10 resolved

| # | Status | Evidence |
|---|---|---|
| 1 | fixed | Cell 128: "spread 3×10⁻¹³". The budget prints `:.4f` (+15.4937, +61.9750, −46.4813), and the text quotes those digits. FTCS text: "about 10⁻⁵ (8×10⁻⁶)", printed 7.9e-06 ✓ |
| 2 | fixed | Planetary term prints `2.9168e-09`; the N39 text gives 2Ω × 2×10⁻⁵ = 2.9168×10⁻⁹ ✓ |
| 3 | fixed | Rotating-tank what-if: no Coriolis force, centrifugal force folded into effective gravity, the Earth's equatorial bulge; geostrophy contrasted correctly ✓ |
| 4 | fixed | Now e^{−2νt} whatever the shape, with a **new demo cell** (circle and square, ν = 0.01 m²/s): both print 0.8187, 0.6703, 0.4705 = e^{−2νt} at 10, 20, 37.7 s. The physics is right: u = e^{−2νt}(sin x cos y, −cos x sin y) is exact Navier–Stokes with ∇²u = −2u, so DΓ/Dt = −2νΓ ✓ |
| 5 | fixed | **D09:** `u0 = (yz² + sin x cos y, xz − cos x sin y, −xy)`. I checked ∇·u0 = cos x cos y − cos x cos y = 0, and the cell prints `its divergence: 0` and a zero residual. The text now describes exactly this field; its remark that a time factor is unnecessary is correct, because the tested identity is instantaneous. **D15:** the text now describes the generic-u test the cell runs ✓ |
| 6 | fixed | Taylor–Green is explained where it first appears (the new what-if demo, "the viscous version of our cellular flow"). A full gloss (with Lamb–Oseen, (3.29)) sits before the budget cell ✓ |
| 7 | fixed | The P76 recap (with `indexing="ij"`) now comes before the from-scratch lid sums ✓ |
| 8 | fixed | A Gauss–Legendre gloss comes right before the lock-exchange Stokes cell, pointing to P143 ✓ |
| 9 | fixed | The `np.deg2rad` recap comes right before the tilt slider ✓ |
| 10 | fixed | The polar vector Laplacian is glossed before D03, with the reason for −u_θ/r² (e_θ turns, Ch. 4 P124). A sympy cell prints `0 0` for the line vortex and the solid body. I checked it by hand: for Γ/2πr the terms are 2 − 1 − 1 = 0; for ωr/2 they are 0 + ω/2r − ω/2r = 0 ✓ |

## Round-1 Should-fix items: 20 of 20 resolved

- **Slip table.** The printed −1/(4π) form now stands alone in column 1.
- **R = 287.058.** Hand value 57 334.198 against the library's 57 334.197.
- **2.422219 printed.**
- **D13 result line** completed with "the line vortex (5.2)".
- **Semi-infinite value, ΣΓx spread and mid-channel value** are all printed now: 0.079577; P spread 0.0; 0.0e+00.
- **Periodic trapezoid.** "Spectrally accurate" now points to P136.
- **Missing recaps added:**
  - `np.cross`, `np.linalg.norm`, `@`, before the loop vector-area cell;
  - `np.gradient`, before the ring cell.
- **Conventions written out** in D07, D21 and D23.
- **Glosses:** the moment-transfer gloss before D06; the quotient rule written into D05 step 9.
- **D19 pointer** now "P138 in C05".
- **Filament figure is readable.**
  - One scale for both arrows; du is ×5 and offset, and the notes say so.
  - The title says "7 % of |u|" (0.0092/0.141 ✓).
  - Only a little of the purple arrow runs past the axes box, which does not matter.
- **Wording:**
  - "cylindrical Laplacian changes sign";
  - leap-frogging "for ever in this thin-core model";
  - baroclinic instability moved to the Ch. 13 bullet;
  - "37 % (e⁻¹)";
  - C05 now has its *In one line*;
  - "stencil-error floor (2.6×10⁻⁵)";
  - D04 check "≈ −3e-14";
  - ν "fifteen times larger";
  - "cat's eye" glossed.

**The point to judge: the frozen-kernel variation 2a/|x − x′| is correct.** The kernel K = (x − x′)/|x − x′|³ has size
1/d².
- Moving the source point x′ by a along x − x′ changes |K| by d(d⁻²) = 2a/d³, a relative change of **2a/d**.
- Moving it sideways by a only turns K, by a relative amount ≈ a/d.
- So 2a/d is the largest first-order variation: 20 % at 10 core radii and 2 % at 100, as the gloss and D12 step 3 now
  say. Averaged over a symmetric core the first-order term cancels, so the true error is even smaller. The gloss says
  "about", which is fine.

## New cells: checked

- **The e^{−2νt} demo.**
  - Its markdown explains the reasoning: (5.11) shown, then ∇²u = −2u, DΓ/Dt = −2νΓ.
  - The code is commented.
  - Its "What does the code above do?" quotes the printed numbers.
  - ν = 0.01 m²/s is labelled as deliberately large.
  - First-use order is fine: `ch05.cellular_flow`, `square_loop_points` and `material_circulation` all come after P136/P137 and the C03 scenario.
- **The polar vector Laplacian gloss and its sympy cell.** They come before D03; D03 step 8 cites them.
- **The Gauss–Legendre gloss, the meshgrid, deg2rad, np.cross/norm/@ and np.gradient recaps.** Each is placed before the
  first use it serves; I checked this against the dump line order.
- **The Taylor–Green gloss** is correct: (cos x sin y, −sin x cos y) is the same flow shifted by π/2, and the gloss
  says so.

## Must fix (round 2)

1. **Equation numbers cited with no equation next to them.** `coverage_check` passes these because each cell contains
   some other maths, but the house rule requires the equation itself.
   - **Cell 8, slip table.** The row "Exercise 5.8 after Eq. 5.14, Exercise 5.11 after Eq. 5.33" is new in round 2,
     from the table rewrite. Eq. 5.33 is not written anywhere in that cell.
     - **Change:** "after (5.14), $\mathbf u=\frac1{4\pi}\int\cdots$" and "after (5.33), $D\Gamma_a/Dt=0$". Compact
       inline forms are enough.
   - **Cell 495, "What later chapters build on".** I missed these in round 1.
     - "the rotating cylinder (8.11)": write the equation, $u_\theta=\omega a^2/2r$ (or Ch. 8's form), or drop the
       number.
     - "the rotating vorticity equation *(5.30)*, absolute circulation *(5.33)*": add
       $\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega+2\boldsymbol\Omega)\cdot\nabla\mathbf u+\dots$ and
       $D\Gamma_a/Dt=0$.

My line scan found no other bare equation numbers. Every other hit is a derivation heading or an ASCII sketch followed
by a "*The numbered equations in the sketch:*" line, or a line whose own cell shows the equation.

## Should fix (round 2)

None.

## Verdict: FAIL (round 2) — 1 Must-fix: write out the cited equations in cells 8 and 495. Once they are written out, the notebook passes; no other open items.


## Orchestrator follow-up (after round 2)
The one remaining round-2 must-fix, equation numbers cited without their equations, was fixed by the orchestrator in `notebooks/build_ch05.py`:
- The slip-table row now writes the sign-corrected (5.14) and $D\Gamma_a/Dt=0$ (5.33).
- The feeds-forward line writes (5.30) and (5.33) in full.
- The Ch. 8 pointer no longer cites a bare (8.11).

The notebook is re-run at the merge gate.

**Verdict: PASS**
