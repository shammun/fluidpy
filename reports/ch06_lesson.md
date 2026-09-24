# Chapter 6 — lesson review                                   2026-09-23

Reviewed: `outputs/ch06/executed.ipynb` (562 cells, 0 errors, 20 static figures, 7 plotly slider figures, 3 animations,
1 live widget, 9 `show_viz` embeds), against `analysis/ch06_curation.md`, `analysis/ch06_design.md` (Parts A, E, F and
the Part G errata G1–G8), `knowledge/primers.md` and `knowledge/concept_map.md`. Every figure image was examined. Every
derivation step was checked by hand. The sympy check cells were read against the steps they claim to test.

coverage_check: **OK** (0 errors, 8 warnings; all 8 say the explainers `laplace_relaxation`, `axial_singularity_bodies`
and `added_mass_sphere` are not built yet. That is expected and is not a notebook defect.) · sections covered **10/10**
(§6.1–§6.10, each with a "What is this section about?") · CORE blocks **15/15 with code + visual** (each has plain words,
idea, maths, tiny example, commented code, from-scratch parity and What you see / How to read it / What would change
if…) · derivations **31/31** present (D01–D31). The five ★★★ derivations (D17, D18, D21, D29, D31) all have sympy cells.

Errata and conventions checked (all respected):
- **Circulation sign.** The book's Γ is clockwise in (6.36)–(6.40), (6.52), (6.62) and Example 6.1; the code uses
  `Gamma_cw=`/`Gamma_ccw=`. A convention table and a numbers callout are in place (cells 6 and 204).
- **(6.82)** is taught as correct: r × the App. B divergence (cell 433, G6).
- **Real slips in corrected form:**
  - (6.61)'s 1/z² coefficient (D18 step 6, with a sympy cell)
  - (6.104)'s middle bracket (D29 step 8, with a sympy cell)
  - the stray dφ in (6.108) (D30 step 7)
- **Example 6.2** converges at order 4/3 because of the r^{2/3} corner singularity (cells 411–413, 561; G2).
- **G1:** the slit is excluded, so the failure fraction is 1.0 (cell 369).
- **G3:** k_n alternate in sign, N is capped at 40, panels are exact on the circle and O(1/N²) on the ellipse
  (cells 484–496, 491).
- **G4:** off-body panel velocity is first order (N90).
- **G5:** the cylinder pressure band is labelled qualitative.
- **G7:** tilted-ellipse lift is perpendicular to the stream (cell 347).
- **G8:** ρ is passed explicitly.
- **No book text or figures copied** (§6.1 and §6.10 spot-checked against `chapters/ch06.txt`: our own wording).
  Example 6.2 is run with our Q = 1 and ψ/Q.

## Must fix

1. **Literal TeX commands in markdown prose (98 hits in 17 cells). They render as raw text such as "\lvert d\rvert".**
   - Cells: 79, 132, 158, 183, 207, 283, 320, 335, 355, 361, 368, 451, 510, 515, 517, 530, 532.
   - The hits sit in derivation *why* / *in words* / plan / trap lines outside `$…$`. Examples:
     - D04 step 4: "\lvert∇φ\rvert = \lvert∇ψ\rvert = speed"
     - D06 plan: "2mε = \lvert d\rvert fixed"
     - D09 step 5: "2U\lvert sin θ\rvert"
     - D15 step 6: "\lvert dw/dz\rvert"
     - D29 steps 5, 7, 11, 13
     - D31 title "½ρ∫\lvert∇φ\rvert² dV"
   - Cell 361 (D20 Result) is worse. Inside `\text{}` it contains `\text{for \backslash lvert ζ\backslash rvert = b …}`,
     which renders as "\ lvert".
   - ch05 had 0 such hits, so this is a regression in `notebooks/build_ch06.py`: Part F's `\lvert … \rvert` convention
     was copied into non-maths text.
   - **Change:** in the builder, convert `\lvert x\rvert` to `|x|` in every non-maths string passed to `nb.derivation` /
     `nb.trap` / `nb.note`, or wrap those fragments in `$…$`. Rebuild, then check with a scan for `\\(lvert|rvert|backslash)`
     outside `$…$`.

2. **Markdown asterisks inside `\text{}` in 43 display equations. The italic markers render as literal "\*(6.1)\*".**
   - Every derivation's "We start from" and "Result" line that cites a number is affected. Examples:
     - cell 17: `\text{*(4.7)* and}` and `\text{*(6.1)* with only}`
     - cell 54: `\text{*(6.3)*}`
     - cells 60, 152, 207, 215, 241–248, 276, 317, 320, 335, 361, 368, 404, 441, 451, 476, 510, 515, 523, 530 …
   - ch05 has 0 such patterns.
   - **Change:** in `build_ch06.py`, strip the `*…*` emphasis from equation labels before they are placed inside
     `\text{}` (write `\text{(6.1)}`).

3. **Bare equation numbers inside derivation bodies and notes, with the equation not written.** This breaks the "show the
   equation, not just its number" rule. Each of these needs the LaTeX next to the number:
   - cell 7 (slip table): "(6.104)'s middle bracket" → write $\mathbf u_a=\tfrac32(\mathbf u_s\cdot\mathbf e_\xi)\mathbf e_\xi-\tfrac12\mathbf u_s$ (6.104).
   - cell 183, D09:
     - plan item 3 "Velocities by (6.21)–(6.22)"
     - step 4 *why* "The polar forms (6.21)–(6.22)"
     - "We start from" "(from Im of (6.49))"
   - cell 335, D18:
     - Tools "the element potentials (6.47)–(6.49)"
     - step 3 "((6.47) with −Γ)"
   - cell 510, D28:
     - step 3 "(6.89, 6.92)"
     - step 5 "(6.16 with U_s = u_s)"
   - cell 515, D29 plan: "(6.101)–(6.102)" and "(6.103)–(6.104)".
   - cell 523, D30 plan: "(6.107)–(6.108)".
   - cell 518 (N96): "The book's middle expression of (6.104) prints…".
   - cell 525 (N103): "The book's (6.108) still shows dφ…".
   - cell 94 "(Eq. 6.18, from 4.72)", cell 306 "(Eq. 6.54, from 4.17)" and cell 502 "(Eq. 6.99, from 4.75)": either show
     (4.72)/(4.17) or drop the provenance number.

4. **Cell 320, D17 step 8: the *why* is wrong for part of the body.** It says the velocity and the step dz "have the same
   direction angle ϑ". On a counterclockwise contour the velocity is **antiparallel** to dz wherever the flow runs
   clockwise round the body, e.g. the top of a cylinder in a stream from the left.
   - **Change:** write $u+iv=\pm\lvert q\rvert e^{i\vartheta}$ (tangent: parallel or antiparallel). The sign appears on
     both sides of (6.59), $(u+iv)\,dz^*=(u-iv)\,dz$, so it cancels. Also make step 9's line read "$=\pm\lvert q\rvert\lvert dz\rvert$".

5. **Cells 475 and 477 (D26 tools gloss and traps) claim the limits "turn around" / are reversed by d(cot α). That is
   false.**
   - With $z-\xi=R\cot\alpha$ at fixed P, the two minus signs cancel: $d\xi=R\,d\alpha/\sin^2\alpha$ (step 4).
   - ξ = 0 → a maps to α = θ → α₁ with α₁ > θ (step 5), so no reversal.
   - The derivation steps themselves are right. Only the gloss and the trap misstate them.
   - **Change:** "d(cot α)/dα = −1/sin²α, and the minus cancels the one from −dξ, so the limits keep their order
     (ξ = 0 → a becomes α = θ → α₁)." Replace the trap with "the minus sign of d(cot α) is cancelled by −dξ; do not flip
     the limits."

6. **Cell 259 ("What you see", two-sources figure c258) makes a claim the figure contradicts.** The claim is that the
   dashed algebraic curves "lie exactly on three of them [the teal streamlines]".
   - The dashed curves are ψ = π/4, π/2, 3π/4.
   - The teal contours use `levels=np.linspace(0.2, 2π − 0.2, 16)`, which contains none of those values.
   - In the image the dashed curves run *between* teal lines.
   - **Change:** in cell 258, add π/4, π/2, 3π/4 to the teal levels, or draw those three levels in teal underneath. The
     text then becomes true. Also start the upper dashed curve at the axis (it begins at y ≈ 2).

## Should fix

1. **Derivation headings and ASCII idea sketches cite numbers without the equation.**
   - Headings in cells 17, 54, 100, 132, 152, 183, 207, 215, 276, 283, 317, 335, 355, 361, 368, 397, 404, 441, 451,
     476, 510, 515, 523. Example: "(6.4) → (6.5), (6.6)".
   - Sketches in cells 15, 98, 311, 353, 394, 473.
   - The equations do follow in the body. Still, put the key equation in the heading (as D17 does) or drop the numbers
     from the heading.
   - The cell 6 convention table's "(6.36)–(6.40)", "(6.61)–(6.62)" could show one representative equation each.
2. **Python idioms used before their reminder.**
   - Cell 26 uses `lambda` and `assert np.allclose`, whose 🔁 reminders come in cell 31.
   - Cell 23 uses a dict comprehension `{k: … for k, v in r.items()}` that nothing recalls.
   - **Change:** move cell 31's lambda / allclose bullets into cell 22, and add a one-line gloss for the dict
     comprehension.
3. **Internal design jargon in learner-facing cells.** "Part G1" (cell 369) and "Part G3" (cells 486, 489, 494 prose, 495)
   mean nothing to a reader. Replace them with the reason itself ("both roots have |ζ| = b on the slit", "the strengths
   alternate; only their moments converge").
4. **Cell 33 "What you see" vs the c032 image.** The text says the background is log₁₀ ≈ −12 … −8 (round-off). The image
   shows a mid-grey halo of ≈ −6 … −4 next to the cylinder: stencil truncation, not round-off. Say "round-off far away,
   ~10⁻⁵ stencil error next to the body — still 10⁸ below the inertia term".
5. **Cell 226 / c225, Γ = 6πaU panel.** The text says "a free point below the body where two streamlines cross". No
   contour passes through the orange dot. Add `levels=[ψ(r₊)]` so the separatrix is drawn.
6. **Cell 294 / c293.** The claimed "dark orange spot at the tip of the 270° and 360° corners" is not visible: the shade is
   uniformly light at `alpha=0.5, vmin=-1.5, vmax=1`. Tighten the colour range or the grid near the tip, or reword.
7. **Figure hygiene.**
   - c466: the inset hides the right end of the sphere C_p curve. Move it to lower right.
   - c378: the "× = ±b" markers are drawn under the ζ-plane disc and cannot be seen. Plot them with a higher zorder.
   - c380: w = ln z shows the cut as a heavy black line on the negative x-axis although the code says "cuts masked".
     Mask it or say what it is.
   - c090: white axis-grid lines cross the Rankine core. Turn off the grid.
   - c256: the legend at fontsize 6 overlaps the streamlines.
8. **c194 band.** The qualitative band leaves the ideal curve at ≈ 45° and plateaus at −1, while the slider/label says
   separation at 80°. Add one sentence saying the band is a sketch: it falls to a shallower minimum before the
   separation angle, then stays flat. Or reshape it so the plateau starts at `sep_deg`.
9. **D21 (★★★) sympy cell 369** tests the roots and Vieta (steps 2–3) and the branch numerically. It does not test the
   derived velocity (step 11). Add: with $z=\zeta+b^2/\zeta$, check symbolically that
   $\frac{d}{dz}W(\zeta(z))=W'(\zeta)/(1-b^2/\zeta^2)$ (implicit differentiation). That way the check covers the result
   actually derived.
10. **Cell 207, D10 step 8.** "multiply by r²/1" should read "multiply by r²".
11. **Cell 29 prints "Re = 50 < 1e+03".** Design rule 13 lists the §6.1 Re threshold as a private book value. Either
    confirm that an order-of-magnitude threshold is acceptable in public, or have `ideal_flow_applicability` word it
    without the book's number.
12. **Recap labels.** Cells 43, 45, 46, 430 and 433–435 read "(Ch. 4 (App. B operators))" / "(Ch. 4 (App. B))". Use
    "(App. B; operators recapped in Ch. 4)".
13. **Missing interpretation.** Cell 224 (`circulation_family_check`) has no "What does this show?". Say why the loop
    circulation is −Γ_cw (counterclockwise loops) and that the far-field error is exactly Γ/2πR.
14. **Cell 407 comment.** "2/(1 + sin(π/N))" returns 1.0718 only with N = 3 grid *intervals*. Say so.
15. **Missing reminders for idioms primed in earlier chapters.**
    - `A @ x` (cell 403; Ch. 2 P63)
    - `np.linalg.norm` (cell 32; Ch. 2 P67)
    - `ax.pcolormesh` (cell 32)
    - `np.meshgrid(..., indexing="ij")` (cell 538)
    - One reminder line each where first used.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / tool | First used | Explained at |
|---|---|---|
| ideal flow, inviscid vs irrotational | 13 | 14–15, D01 (17) |
| Kelvin's theorem, Mach number, baroclinic torque | 11 | 11 (recap) |
| Reynolds number (4.103) | 28 | 28 |
| plane Poiseuille flow | 23 (output) | 24 (right after; Ch. 4 §4.6) |
| `lambda`, `assert np.allclose` | 26 | 31 (after), named in 9 — **move earlier (Should 2)** |
| dict comprehension | 23 | **missing** (Should 2) |
| `np.linalg.norm`, `ax.pcolormesh` | 32 | **missing reminder** (Should 15) |
| ψ, polar continuity, polar Laplacian | 40–46 | recaps 39–46 |
| harmonic | 49 | 49 |
| Dirac delta in the plane, 2-D divergence theorem | 57 (D02 step 5 names it) | primer 58 (before D03) |
| simply connected | 42 | 42 |
| level sets, directional derivative | 79, 100 | 78, 99 |
| complex numbers as points (`-1+0j`) | 107 | gloss 99; primer P153 at 270 |
| Newton's method (via `stagnation_points`) | 107 | pointer 108, primer 211 |
| limit with product fixed, O(ε²) | 131–132 | 122, 130 |
| `np.arctan2` range / branch of θ | 150 | 150 |
| `brentq` | 164 | 150 |
| trapezoid / periodic sums | 187 | 182 |
| integrals of sin/cos over a period | 215 (D11) | primer 186 |
| Vieta | 207 (D10) | 206 |
| images, drift | 233–236 | recaps 233–235 |
| arctan derivative with moving argument, unsteady Bernoulli (4.83) | 248 (D13) | gloss 241 |
| Milne-Thomson circle theorem, w(z) | 233, 177 | forward pointer to C09 |
| complex plane, analytic, Cauchy–Riemann | 270–276 | primers 270, 273 |
| branch cuts, complex log and powers | 282–283 | primer 281 |
| Cauchy's theorem | 311 | primer 312 |
| residue | 314 | pointer to primer 326 |
| Laurent series, FFT on a circle | 326–331 | primer 326, recap 325 |
| sympy `.coeff`, `residue` | 334 | primer 333 |
| complex square roots / principal branch | 367 | primer 366 |
| grid differences (6.70), (6.71) | 388–391 | recaps 388–390 |
| spectral radius, `np.linalg.solve` | 404 (D23) | 395 |
| Jacobi / Gauss–Seidel / SOR, residual | 403 | primer 402 |
| `A @ x` matrix product | 403 | **missing reminder** (Ch. 2 P63) |
| boolean masks, `scipy.sparse` | 408–412 | primer 409 |
| Stokes stream function, spherical operators | 426–436 | recaps 426–435 |
| cylindrical cross products, 3-D source flux | 441 (D24), 451 (D25) | gloss 440 |
| substitution in an integral, cot substitution | 476 (D26) | 475 (gloss **misstates** limits — Must 5) |
| collocation, condition number, `np.vander` | 483–486 | primer 482 |
| chain rule along a path, ∂/∂x_s = −∇ | 515 (D29) | 509 |
| surface integrals on a sphere, `dblquad` | 522–523 | primer 521, recap 475 |
| Green's first identity | 530 (D31) | primer 527 |
| Gauss–Legendre, `solve_ivp` | 535 | 534 |
| `np.meshgrid(indexing="ij")` | 538 | **missing** (Should 15) |
| live widgets | 548 | 546 |

## Derivation audit

| D | Steps | Every step follows? | why / in-words ok | Check ok | Verdict |
|---|---|---|---|---|---|
| D01 ideal-flow eqs (6.1) | 8 | yes | yes | units, μ = 0 limit, numbers (cell 23) | pass |
| D02 ω = −∇²ψ (6.4)–(6.6) | 5 | yes | yes | units, Rankine core 31.8 s⁻¹ | pass |
| D03 ln r and the delta | 9 | yes | yes | flux on 0.01…100 m, off-centre 0 | pass |
| D04 orthogonality (6.10) | 5 | yes | yes (literal `\lvert`, Must 1) | corner at (1,1), 1.5e-10 | pass after M1 |
| D05 wall = streamline (6.16) | 5 | yes | yes | half-body, cylinder | pass |
| D06 doublet limit (6.29) | 9 | yes (ε² terms cancel, checked) | yes (Must 1) | ε²/3 at (1,0), order 2.01 | pass after M1 |
| D07 half-body (6.31) | 10 | yes | yes | mass balance = geometry | pass |
| D08 half-body C_p zero | 8 | yes (tan θ = −2(π − θ), 113.2°) | yes (Must 1) | θ = 90°: −0.405 | pass after M1 |
| D09 cylinder (6.33)–(6.35) | 6 | yes | bare (6.21)–(6.22) (Must 3) | 60/−180 Pa | pass after M1, M3 |
| D10 circulation stagnation (6.38) | 9 | yes (quadratic, r₊ = 0.2618) | "r²/1" typo (Should 10) | Vieta product a² | pass after M1 |
| D11 L = ρUΓ (6.40) | 9 | yes | yes | 64-pt sum 24.000 | pass |
| D12 image rules | 5 | yes | yes | wall v ≤ 1e-16 | pass |
| D13 Example 6.1 wall pressure | 11 | yes (re-derived: −Γ²/4π²(h²+s²), v = Γh/π(h²+s²), max at s = √3h) | yes | numeric route agrees 7 digits | pass |
| D14 Cauchy–Riemann, u − iv | 8 | yes | yes | z², z* control | pass |
| D15 corner w = Azⁿ (6.46) | 7 | yes | yes (Must 1) | slopes n − 1 | pass after M1 |
| D16 force from pressure (6.56) | 8 | yes | yes | circle → D11 step 4 | pass |
| D17 Blasius (6.60) ★★★ | 12 | step 8 claims same direction; antiparallel on part of the body | **wrong why** (Must 4) | sympy tests steps 6, 9, result, contour move | **fail → fix M4** |
| D18 Kutta–Zhukhovsky (6.62) ★★★ | 11 | yes (1/z² coefficient re-derived) | bare (6.47)–(6.49) (Must 3) | sympy: residue, true vs printed coefficient, D, L | pass after M3 |
| D19 conformal angles (6.64) | 6 | yes | yes | 90°→90°, 45°→90° | pass |
| D20 Zhukhovsky ellipse (6.67) | 8 | yes | Result renders `\backslash lvert` (Must 1) | c² = 4 | pass after M1 |
| D21 inverse map (6.69) ★★★ | 11 | yes (cuts of √(z∓2b) checked) | yes | sympy covers steps 2–3 only (Should 9) | pass (weak check) |
| D22 five-point rule (6.72) | 6 | yes | yes | orders 1.94/1.97, xy exact | pass |
| D23 Aψ = b, Gauss–Seidel | 9 | yes (ρ_J = ½, ρ_GS = ¼ checked) | yes | 34/19/12 sweeps | pass |
| D24 Stokes ψ (6.77) | 8 | yes (cross products, sign) | yes | ½UR², Rz control −z/R² | pass |
| D25 sphere (6.86)–(6.91) | 10 | yes (dipole −2εQ e_z, u_θ checked) | yes (Must 1) | 1.5U, −1.25 | pass after M1 |
| D26 line sink / airship (6.95) | 10 | yes | steps right; gloss 475 and trap 477 say limits reverse (Must 5) | quad = closed form, far field | **fix M5** |
| D27 axial singularity method | 6 | yes | yes (G3 honest) | error/cond table | pass |
| D28 moving-sphere φ (6.97) | 5 | yes | bare (6.89, 6.92), (6.16) (Must 3) | n·∇φ = n·u_s | pass after M3 |
| D29 moving-sphere pressure (6.105) ★★★ | 14 | yes (\|A − B\|² expansion checked) | yes (Must 1) | sympy: chain rule, (6.104) vs printed, Bernoulli = (6.105) | pass after M1 |
| D30 force, added mass (6.108)–(6.109) | 9 | yes (πρa³ × 2/3) | yes | 2.094 kg, 2g bubble, surface sum | pass |
| D31 added mass by energy ★★★ | 10 | yes (inner-normal sign, R⁻³ far term) | yes (Must 1 in title) | sympy vol = surf, M, far term | pass after M1 |

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **Convention discipline.** A conventions table with numbers (cell 6) and a callout at first bite (cell 204). Code keys
  are named for the convention (`Gamma_cw=`/`Gamma_ccw=`). The same numeric case (stagnation points −9.16°/−170.84°,
  L = 24 N/m) appears in the table, the tiny example, the code and the explainer. Good pattern for knowledge/.
- **Book slips handled as teaching moments.** Each slip is corrected where it is used, and the sympy cell prints both the
  correct and the printed form ((6.61) coefficient, (6.104) bracket). (6.82) is correctly *not* flagged. Reuse for later
  chapters.
- **Honest numerics.**
  - Example 6.2's order 4/3 is explained by the corner exponent of (6.46), with pairwise orders printed.
  - The axial method shows that the shape converges while the strengths do not.
  - Circle panels are exact and ellipse panels O(1/N²).
  - Candidate rule for knowledge/: "state what converges and what doesn't".
- **Two independent routes for every load-bearing number.** Lift by the surface integral, the far-field momentum CV,
  Blasius on many contours and residues. Added mass by the force and by the energy (surface and volume). Example 6.1 in
  closed form and by moving vortices.
- **From-scratch parity cells** after every library call, each with an assert and a one-line conclusion.
- **Checks re-run the construction.** The ★★★ sympy cells rebuild the chain rule step by step (D29 cell 516 compares
  direct ∂φ/∂t with the chain-rule form) and do not only test the final formula.
- **Climate links where real:** ψ inversion as Poisson (C02, C12), deformation flow and fronts (C04), bottom-pressure
  recorders (C08), the Stokes geometry of the zonal-mean overturning streamfunction (D24).

## Verdict: FAIL

Six Must-fix items. Two are systematic rendering defects from `notebooks/build_ch06.py` (M1: literal `\lvert` / `\backslash`
in prose, 98 hits; M2: `*…*` inside `\text{}`, 43 displays). One is the show-the-equation rule (M3, about 20 bare
references in derivation bodies and notes). Two are wrong derivation explanations (M4 D17 step 8; M5 the D26
limits gloss/trap). One is a false figure claim (M6, cell 259). All are local builder or text fixes. The physics, the
signs and the errata handling are correct, and the teaching structure is complete. After these fixes and a rebuild, I
expect this chapter to pass.

---

# Round 2 — re-review                                          2026-09-24

Re-read `outputs/ch06/executed.ipynb` (564 cells, 0 errors, 20 figures re-extracted and inspected). coverage_check:
**OK** (0 errors, 4 warnings — explainers not built yet, not a notebook defect). Scans re-run: literal TeX commands
outside maths **0** (was 98); `*…*` inside `\text{}` **0** (was 43); "Part G" in learner text **0**. The 22 bare-number hits
left by my heuristic are all acceptable: ASCII idea sketches followed by their "numbered equations in the sketch" line
(cells 15, 98, 312, 354, 395, 475), repeat mentions in a cell that already displays the equation (126, 280, 324, 520, 522,
527), or non-equation numbers (91, 173, 398, 486).

## Round-1 items — status
| Item | Status | Evidence |
|---|---|---|
| M1 literal `\lvert`/`\backslash` in prose | fixed | scan 0; D20 Result renders |
| M2 asterisks inside `\text{}` | fixed | 0 occurrences |
| M3 bare numbers in derivation bodies/notes | fixed | D09, D18, D28, D29, D30 plans/steps, N96, N103, slip table now show the equations; provenance "from 4.72/4.17/4.75" handled |
| M4 D17 step 8 direction | fixed | steps 8–9 now $u+iv=\pm\lvert q\rvert e^{i\vartheta}$ and the ± cancels; step 12 adds a correct note (a circle cutting the cylinder still encloses its only singularity z = 0 → ρUΓ; for the ellipse it crosses the map's cut) |
| M5 D26 limits | fixed | gloss (cell 475 area) and trap: "limits keep their order (ξ = 0 → a becomes α = θ → α₁)" |
| M6 two-sources caption | fixed | levels kπ/8 include π/4, π/2, 3π/4; c259 shows the dashed curves on teal lines |
| S1 headings | fixed | every derivation heading now carries its equation in LaTeX |
| S2 lambda/allclose/dict comprehension | fixed | reminders and a dict-comprehension gloss before first use (cell 22) |
| S3 "Part G" jargon | fixed | none left |
| S4 c032 caption | fixed | caption states round-off far away, ≈ −6…−4 next to the body |
| S5 separatrix at Γ = 6πaU | fixed | c226 draws the crossing ψ(r₊) contour in orange |
| S6 corner tip shading | fixed | c294 shows the bright tip at 270° and 360° |
| S7 figure hygiene | fixed (c379 markers visible, c381 cut masked, c090 grid off, c257 legend moved) — c468 see new item N2 |
| S8 band caption | fixed | caption describes the shallower minimum near −1, plateau until separation (intended qualitative sketch) |
| S9 D21 sympy covers step 11 | fixed | "step 11 − chain rule: 0" |
| S10 "r²/1" | fixed | "multiply by r²" |
| S11 Re threshold | fixed | verdict text no longer prints the book's number |
| S12–S15 | fixed | recap labels, cell 224 interpretation, ω_opt comment, idiom reminders present |

## New in round 2 (Should fix; none blocks)
- **N1 · Glued absolute-value bars.** The `\rvert`→`|` conversion swallowed the following space in ~15 prose places:
  "2mε = |d|fixed" (D06 plan), "−|d|e_x … |d|x" (D06 step 8), "|u|= 2U|sin θ|" (D09 step 6), "|sin θ|would exceed"
  (D10 step 8), "|dw/dz|= |u − iv|= speed" and "|dw/dz|vs r" (D15), "|f′|and" (D19), "|ζ|grows" (D20), "|u|= (3/2)U"
  (D25), "|ξ|u_s·e_ξ" (D28), "|u_s|cos θ_s" (D29). Fix the converter to keep the trailing space (`\rvert ` → `| `).
- **N2 · c468 inset** was moved to lower right and now covers the cylinder's dashed C_p curve between θ ≈ 120° and 140°.
  Put the inset in the upper middle (C_p ≈ 0.5…1, θ ≈ 60…120°, empty there) or shrink it.

## Verdict: PASS
All 6 Must-fix and 15 Should-fix items of round 1 are fixed and verified in the executed notebook; no new Must-fix.
Two cosmetic Should-fix items (N1, N2) remain.
