# Chapter 10 — lesson review                                   2026-09-30

coverage_check: **OK** (0 errors, 2 warnings — D22/D23 explainers `mixed_fe_lbb`, `lid_driven_cavity` not built yet, ignored as instructed)
· sections covered **6/6** (10.1–10.6) · CORE blocks **15/15** with code + visual · derivations 23/23 written out
(D18 ★★★ with a sympy cell that re-runs the construction) · 27 static figures looked at · 0 execution errors
(1 stray `UserWarning` in the output).

Reviewed: `outputs/ch10/executed.ipynb` (589 cells) read in order; dump and figures in `outputs/ch10/review/`.
Pages p469 (10.110) and p481 (10.147–10.152) rendered and compared: they match. A sample of prose numbers was recomputed
(FTCS tiny example, |G| = 1.04 and 939 steps, CFL time steps 126 s / 78 s, D16/D17 tiny examples, D03 E = 0.0406, element
blocks, (10.110) = 2.935e-4, 28 000-step estimate, Richardson): all agree with the cells, apart from the items below.
The brief's required statements are present: pre-asymptotic default-grid orders (BTCS 1.80, upwind 0.76, Θ = ¼ 0.86),
MAC 64² 0.26 %, MacCormack 1.38 % / 0.44 %, ψ_min p = 2.00, block C_D p = 0.43 non-asymptotic, MacCormack mass drift
and worse at lower Ma on a fixed grid, confined St/C_D qualitative. **Missing: the MacCormack default-Δt DEVIATION (M2).**
Public-repo check: `tools/check_public.py` passes on the notebook and builder; none of the `_forbidden_public` strings
appear; an 8-word shingle scan against `chapters/ch10.txt` finds no copied prose (only stencil subscripts and one
short technical phrase).

## Must fix

1. **Cells 15–16, 491–492 (and 406–407): the text describes things the figures do not show.**
   `spacetime_grid` (scripts/ch10_drawings.py:18) draws no highlighted stencil, yet the cell-15 comment says "the FTCS
   stencil highlighted" and cell 16 says "four dots highlighted in orange — three at tₙ feeding one at tₙ₊₁" and "the orange
   stencil is FTCS". `cavity_sketch` (:169) draws only box, lid, U and L, yet the cell-491 comment and cell 492 describe
   "the dashed vertical centreline x = ½ and a sketch of the primary eddy with two small corner eddies". `staggered_grid`
   draws no thick boundary faces, yet cell 407 says "the boundary faces drawn thick".
   *Change:* add the orange FTCS stencil (three tₙ nodes → one tₙ₊₁ node, with connecting lines) to `spacetime_grid`, and the
   centreline plus eddy cartoon to `cavity_sketch`, thicken the outer faces in `staggered_grid` — or rewrite the three
   comments and reading notes to describe what is actually drawn.

2. **Cells 359–361, 498, 504–506, 526, 543: the MacCormack step and boundary DEVIATIONs are not stated.**
   `MCK.cavity_maccormack` defaults to `dt_rule="additive"`, i.e. Δt = σ/[|u|/Δx + |v|/Δy + c√(1/Δx²+1/Δy²) + 2ν(1/Δx²+1/Δy²)],
   not (10.110) $\Delta t\le\frac{\sigma}{1+2/Re_\Delta}\big[\frac{|u|}{\Delta x}+\frac{|v|}{\Delta y}+c\sqrt{\frac1{\Delta x^2}+\frac1{\Delta y^2}}\big]^{-1}$
   (code marks it `# DEVIATION`, reports/ch10_review.md M3). The notebook teaches only (10.110) (N58) and (10.155) (N91).
   The tiny example of cell 504 even estimates "about 28 000 steps" with (10.155), while the cached 64² run uses the
   additive rule (about 5.6e-4 → ≈ 36 000 steps). The top-corner density extrapolation (`corner="extrapolate"`) and the
   block's zero-gradient inflow density (`block_channel`) are also DEVIATIONs from §10.5 that the reader never hears about.
   *Change:* after N58 (cell 359) add a `> ⚠️ Our deviation` box: the formula of the additive rule, why (the book itself calls
   (10.110) "quite conservative" at small mesh Reynolds numbers; the inviscid bracket alone blows up at low Re), and a one-line
   cell printing both steps for the 64² cavity. In C14 (after cell 505) and C15 (after cell 519/543) add one sentence each
   on the corner and inflow density treatments. Fix cell 504 step 2 to use the Δt the solver uses (or say that (10.155) is
   only an estimate and the run uses the additive rule).

3. **Cell 405: a code comment contradicts the printed output.** The comment says the shapes are `(2, 3) (2, 4) (3, 2)`;
   the cell prints `(2, 3) (2, 4) (3, 3)`. It also says "one more face than cells in each direction", which is true only
   along the staggering direction of each array. *Change:* comment `# (2, 3) (2, 4) (3, 3): u has one extra column (x-faces), v one extra row (y-faces)`.

4. **Cells 336, 339–340, 439–440: "Taylor–Green vortex" is used before it is explained** (not in knowledge/primers.md or
   concept_map.md). It is the verification tool of C10 and C12. *Change:* a 📎 primer before cell 339: the exact periodic
   solution $u=-\cos x\sin y\,e^{-2t/Re}$, $v=\sin x\cos y\,e^{-2t/Re}$ (state the form `MAC.taylor_green` uses), why it is
   exact (nonlinear term balanced by pressure, pure viscous decay), plus a 2-line check.

5. **Equations referred to by number without the equation.**
   - cell 170: "…u < 0 with the printed stencil (10.29): …" → add $T^{n+1}_i=T^n_i-2\alpha(T^n_i-T^n_{i-1})$.
   - cell 248 (D13 step 7 *why*): "The book prints (10.67)–(10.68) with the labels …" → show the printed pair (as cell 251 does).
   - cell 385 (D19 goal): "the discrete (10.124)" → show $\nabla^2_dp^{n+1}=\nabla_d\cdot\mathbf u^{n+1/2}/\Delta t$.
   - cell 408 (N64): "the same for v on its faces *(10.120)*" → write the v-predictor.
   - cell 564 (N95): "The Galerkin statements *(10.160)–(10.162)* are the same …" → show at least (10.160).
   - cell 570 (N98…): "(10.168)–(10.169)" → show the expansions $u'=\sum_Au_AN^u_A$, $p'=\sum_Bp_BN^p_B$.
   - cell 574 (N104…): "the element matrices (10.188)–(10.197)" → show one representative entry (e.g. the $B^{up}$ element entry).
   - cell 6 (slip table R4): "the book prints '(10.120)'" → add the v-predictor equation in the cell.

6. **False or unsupported factual claims in prose.**
   - cell 504: "A honey-like liquid … ν = 10⁻⁴ m²/s" — honey has ν ≈ 10⁻³–10⁻² m²/s; 10⁻⁴ m²/s is a light oil / glycerol–water
     mix. *Change:* "a light oil (ν = 10⁻⁴ m²/s, about 100 × water)".
   - cell 306 (and N49, cell 304): "uΔx/2 = 500 m²/s — larger than typical eddy diffusivities used in coarse ocean models" —
     coarse (≈ 1°) models use O(1000) m²/s for isopycnal/GM mixing. *Change:* compare with what is used at the stated 10 km
     resolution (typically tens to ~100 m²/s), and say "at this resolution".

## Should fix

- **Symbol conventions announced in cell 6 but broken later.** (a) Cells 398, 569, 570 show (10.129)–(10.133), (10.163),
  (10.165) with bare α, β, θ, although cell 6/N74/N97 announce $\alpha_\Theta,\beta_\Theta,\Theta$ and $\alpha_t,\beta_t$ — and θ is
  the Fourier angle. (b) Cell 433 shows (10.119)–(10.120) with the book's $f$, $g$ body force; N64 says we write $(g_x,g_y)$.
  (c) (10.151)–(10.155) are displayed with $M$ (cells 360, 361, 522, 528–531) while the table says Mach is $Ma$ — add
  "$M\equiv Ma$ in the book's displays" once. (d) The table says "N cells, n_el elements", but the notebook uses $n$ for
  cells/elements throughout (D16, C07, C13 table, (10.92)) — change the table row. (e) Not announced: $\delta$ (layer
  thickness) vs Kronecker $\delta_{AB}$, $\delta_{e\,n_{el}}$; bold $\mathbf F$ as FE force vector (C07) and as y-flux (C10);
  $G$ (amplification) vs $G_A$ (D11); $p$ pressure vs order $p$ (D23); $\mathbf D[\mathbf u]$ — D22 step 5 "the bracket is 2D" reads as
  "2 × diffusivity": write $2\mathbf D[\mathbf u]$.
- **Label collisions.** Recap IDs "R02, R04, R05, R10" are cited (D01, D04, D06, D08, D17, D22) but no recap cell shows an ID,
  and the printed slips are "R1…R12" — "slip R4" vs "R04" is confusing. Rename slips (e.g. "misprint 1–12") or show the
  recap IDs in the recap headings. Explainers are cited as "E2", "E3", "E6" (cells 109, 345, 420) while cell 1 numbers them 1–8.
- **Boilerplate comments.** Many code comments are machine fillers that restate syntax or are wrong: "# print the result",
  "# loop over the cases" on time-step loops (cells 159, 526), "# a coloured map of the field" on contour lines (cell 511),
  "# axis labels (with units), limits and the panel's message" on `set_title`. Replace with meaning + unit comments
  (teaching-style §3).
- **Cell 159**: the output shows a matplotlib `UserWarning` with a local path (`C:\Users\...\ipykernel_…`). Drop
  `fig.subplots_adjust` (the house style uses a layout engine) or set `layout=None`.
- **Figures.** Cell 87: x-tick labels of (a) overlap ("2×10⁻² 3×10⁻² 4×10⁻²") — use `NullFormatter` for minor ticks as in cell 541.
  Cell 128 (c): the "2β = 1" label sits on the (0, 0.5)/(0, 0.51) dots and the (0.2, 0) dot is clipped; cell 129 calls the
  θ ∈ [0, π] curves "closed" (they are half-curves; the amber one is a segment on the real axis). Cell 222: purple half-hats
  $N_0$, $N_n$ and stray teal/purple zero segments are not mentioned in cell 223. Cell 392: titles print raw "u^{n+½}" —
  use mathtext. Cell 475: both pressure panels are labelled "(b)" and have no colour bar. Cell 509 (b): the "MAC 32²" curve is
  hidden under the 64² curve, yet cell 510 compares it with MacCormack 32² — change style/offset or drop the claim.
- **Cell 41** "float32: the floor rises to about 10⁻³ and moves to h ≈ 10⁻³" — by P222 the forward floor is ≈ 10⁻⁴ at
  h ≈ 3 × 10⁻⁴ and the centred one ≈ 10⁻⁵ at h ≈ 5 × 10⁻³; compute it in a float32 cell or quote those.
- **Cells 192 (D10 check) and 198**: "natural slope approached at order 1" / "roughly like h" — the legend prints 0.79; say
  "≈ 0.8 on these meshes". "Nodal values are second-order accurate" is asserted but never computed — add the number or drop.
- **Cell 546**: "Beyond 64² the curves differ by less than Ghia's own grid error" is unsupported; say instead that the
  64²→128² change is below the 0.26 % deviation (compute it).
- **Cell 519 (N86)** says the block is run at Re = 20 **and 100**, but only Re = 20 appears; `reference/ch10/block_forces_re100_dx16.csv`
  exists — plot its lift history (shedding, confined St) or drop "and 100".
- **Cell 363/364**: the numerical phase speed of panel (c) is used but never defined — add one line: $c_{num}/u=-\arg G/(C\theta)$.
- **Cell 220 (N29)** "finite differences in time (the θ-scheme)" — remind it is Ch. 8's θ-method (θ = 0, ½, 1), not
  Glowinski's Θ-scheme of N74.
- **Cell 285 (D15 "What it means")** "with R in the role of √Re's square" — write "δ/L ~ 1/R here, ~ 1/√Re in Ch. 9".
- D15 step 1 *why* "Move uT′ to the right-hand side's side" — garbled; "subtract uT′ from both sides".
- D17 tiny example (cell 308) "exact values 0.00001" — print 6.0 × 10⁻⁶.

## Unexplained-first-use list (term · first cell · explained at cell / missing)

| Term / symbol / tool | First used | Explained at |
|---|---|---|
| big-O, order of accuracy | 0 (road map), 19 | 23 (P221) ✓ |
| FTCS, α, β | 0, 12 | 47–51 (C02, D02) ✓ |
| Noye's region, CFL, Lax theorem, weak form, LBB | 0 (road map only) | 114, 146, 158, 185, 465 ✓ (road-map preview acceptable) |
| machine epsilon, round-off floor | 38 | 38 (P222) ✓ |
| `np.roll` | 56 | 57 ✓ |
| ghost node | 58 | 58 (P223) ✓ |
| error norms | 55 (pointer) | 74 (P224) ✓ |
| half-angle identities | 104 | 104 (P225) ✓ |
| domain of dependence | 141 | 142 (P227) ✓ |
| well-posed | 156 | 156 (P228) ✓ |
| Péclet R, R_cell | 0, 164 | 164 (P229) ✓ |
| Lamb-wave speed | 154 | **missing** (one clause: the fastest external acoustic–gravity mode of the atmosphere, ≈ 310–320 m/s) |
| H¹, 𝒮, V | 182 | 182 (P230) ✓ |
| fundamental lemma | 190 | 190 (P231) ✓ |
| bilinear form | 206 | 206 (P232) ✓ |
| θ-scheme (time) | 220 | **missing** reminder of Ch. 8's θ-method (clashes with Glowinski's Θ) |
| scatter-add, `np.add.at` | 255 | 255 (P235) ✓ |
| linear recurrence | 291 | 291 (P236) ✓ |
| GLS (Galerkin/least-squares) | 307, 462 | **missing** acronym never expanded (465 names it only) |
| Ma (Mach) | 328 | 336 (N53) — Ch. 1/4 term, acceptable |
| Taylor–Green vortex | 336, 339 | **missing** (Must fix 4) |
| Lagrange multiplier | 330 | 331 (P239) ✓ |
| conservation form, predictor–corrector | 341, 343 | P237, P238 ✓ |
| Lax–Wendroff | 345 | 345 (D18 step 5) ✓ |
| numerical phase speed | 363 | **missing** (Should fix) |
| Gibbs-like overshoot | 345 (D18) | **missing** (one clause) |
| commutator, splitting | 373 | 373 (P241) ✓ |
| Helmholtz–Hodge | 383 | 383 (P240) ✓ |
| half-index arrays | 404 | 404 (P242) ✓ |
| compatibility / singular system, `np.linalg.lstsq` | 411–412 | P243 ✓ for the idea; **`lstsq` missing** |
| `diags`, `kron`, `splu` | 413 | 413 (P244) ✓ |
| SVD null space | 424 | 424 (P245) ✓ |
| saddle-point matrix, `np.block`, `eigvalsh` | 456–457 | P246 ✓ for the idea; **`np.block`, `eigvalsh` missing** |
| GMRES | 458 | 458 (P250) ✓ |
| `eigh(A, B)` | 463 | 463 (P247) ✓ |
| barycentric coordinates, Jacobian | 466, 468 | P248, P249 ✓ |
| Euler's formula V − E + T = 1 / 0 | 470, 559 | **missing** (one line: vertices − edges + faces = 1 for a disc, 0 with one hole) |
| `tripcolor`, `tricontourf`, `triplot` | 475, 553, 562 | **missing** (one line: plotting on a triangle mesh) |
| verification vs validation | 488 | 488 (P253) ✓ |
| `np.interp`, cached runs | 500–503 | P251, P252 ✓ |
| Kovasznay flow | 572 | named as "exact steady NS flow" — acceptable |
| `re.findall` | 583 | comment only — give one line |
| Hopf bifurcation | 578 | 578 ✓ |

## Derivation audit

| D | steps | every step follows? | why / in-words ok | check ok | verdict |
|---|---|---|---|---|---|
| D01 ★ | 9 | yes | yes | units + numbers on sin(1) ✓ | OK |
| D02 ★ | 6 | yes | yes | units, tiny example, u = 0 → Ch. 1 ✓ | OK |
| D03 ★★ | 11 | yes (time series added) | yes | units, tiny example 0.0406, sympy, measured residual within 1.6 % ✓ | OK |
| D04 ★★ | 11 | yes; index rename announced | yes | θ = 0 gives 1; tiny G = 0.6 − 0.2i ✓ | OK |
| D05 ★★ | 8 | yes | yes | numbers ✓ | OK |
| D06 ★★ | 13 | yes (checked with sympy: f(0) = 8(2α² − β), f(1) = 8β(2β − 1)) | yes | three cases + brute-force scan ✓ | OK |
| D07 ★★ | 7 | yes | yes | β = 100 → 1/401 ✓ | OK |
| D08 ★★ | 9 | yes | yes | C = ½, θ = π ✓; C = 1.1 → 1.2 ✓ | OK |
| D09 ★★ | 9 | yes | yes | units, x²/x tiny example, sympy ✓ | OK |
| D10 ★★ | 9 | yes | yes | sympy ✓; "order 1" vs printed 0.79 (Should fix) | OK |
| D11 ★★ | 13 | yes | yes | sympy for 3 hats, units ✓ | OK |
| D12 ★★ | 10 | yes | yes | row sums, (−1.5, 2, −0.5) ✓ | OK |
| D13 ★ | 7 | yes; step 7 cites (10.67)–(10.68) without showing them (Must fix 5) | yes | slopes ±4 ✓ | fix M5 |
| D14 ★★ | 11 | yes | yes | numbers, Gauss–Legendre parity ✓ | OK |
| D15 ★ | 8 | yes | step 1 wording garbled; "√Re's square" confusing (Should fix) | limits, 0.3561 ✓ | OK |
| D16 ★★ | 12 | yes (sympy re-run in cell 294) | yes | n = 4 tiny example, P → 0 ✓ | OK |
| D17 ★★ | 8 | yes (R_cell held fixed stated at step 7) | yes | numbers 0.0667… vs 0.0708… recomputed ✓ | OK |
| D18 ★★★ | 14 | yes (sympy re-runs steps 1–5, 6–12, 13–14) | yes | tiny example, C = 1, measured order 1.96 ✓ | OK |
| D19 ★★ | 9 | yes | yes | units, tiny projection ✓; goal cites (10.124) bare (Must fix 5) | fix M5 |
| D20 ★★ | 12 | yes | yes | pipe example, rank 63 ✓ | OK |
| D21 ★★ | 7 | yes | yes | null-space counts 4 / 1 ✓ | OK |
| D22 ★★ | 11 | yes | step 5 "2D" ambiguous (Should fix) | sympy identity ✓ | OK |
| D23 ★★ | 8 | yes | yes | f = 1 + h² tiny example, parity with tools ✓ | OK |

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **Pre-asymptotic honesty as a lesson**: C03 (BTCS 1.80 → 1.98, upwind 0.76 → 0.95 on finer grids), C10 (upwind 0.70),
  N74 (Θ = ¼ 0.86 → 0.995) and the block drag labelled a *failed* convergence study (p = 0.43) — the reader learns that an
  observed order is only meaningful in the asymptotic range. Candidate lesson for `teaching-style` §6.
- **Printed slips as executable counter-examples**: `upwind_printed` explodes to 1.4e16, the printed Step 5 loses the density in
  35 steps, the printed (10.172) gives a 1.7e55 Poiseuille error — each slip is a code option plus a test.
- **Conventions table up front** for heavily overloaded letters (i/√−1, α, β, θ, g, D, M, K, S, n, R) — keep the pattern;
  just enforce it everywhere (Should fix).
- **Closed-form discrete solution** $r^j$ for the wiggles (D16) with the sympy re-run and a two-line from-scratch check against
  the tridiagonal solver — a model for "derive the discrete solution, not only the continuous one".
- **Every stability derivation has a brute-force scan** (3592 (α, β) points against Noye; 2-D von Neumann at 0.95/1.05 × limit).
- **Climate hooks with numbers** (126 s ocean / 78 s atmosphere CFL steps, the C-grid, dynamics/physics splitting).
- **Figures of the machinery itself**: assembly heat maps after each element, the projection stages as a frame player,
  the saddle-point sparsity pattern.

## Verdict: FAIL

Six Must-fix items, all local edits (two drawing helpers, one DEVIATION box and three sentences, one comment, one
primer, eight equation displays, two factual corrections). No derivation is wrong; every CORE block is complete with
code and a figure; the brief's required numbers are present and correct.

---

## Round 2 — re-review of the rebuilt notebook                    2026-09-30

coverage_check: **OK** (0 errors, 0 warnings) · 602 cells · 0 execution errors · **no `UserWarning` / local path in any output**
(cell 159's warning is gone) · 28 figures re-extracted to `outputs/ch10/review/`; the changed ones looked at by eye.

### Must-fix status (round 1 → now)

| # | Round-1 item | Status | Evidence (new cell numbers) |
|---|---|---|---|
| 1 | drawings contradict text | **fixed** | cell 15: orange FTCS stencil (three tₙ nodes → square Tᵢⁿ⁺¹); cell 498: dashed x = ½ centreline, clockwise primary-eddy loop, two corner eddies; cell 412: thick outer faces — all three match the reading notes (16, 499, 413) |
| 2 | MacCormack DEVIATIONs unstated | **fixed** | cell 364 box shows (10.110) and the additive rule in full with the reason; cell 365 computes both (Re = 100: 2.79e-4 vs 5.62e-4, ratio 2.0; Re = 10: 6.0; Re = 1: 9.2 — matches cell 366's "six to nine"); corner box after N84 (cell ~497), inflow box in N86 (cell 528, ρ₀ = (4ρ₁ − ρ₂)/3 with its one-sided-derivative reason); cell 511/512 step counts computed: MAC 4 096, MacCormack 35 580 (5.62e-4), ratio 8.7 — equal to my round-1 recomputation |
| 3 | shape comment | **fixed** | cell 405 comment now `(2, 3) (2, 4) (3, 3)` with the right per-array explanation |
| 4 | Taylor–Green unexplained | **fixed** | primer P254 (cell 336) before first use: $u=\sin x\cos y\,e^{-2t/Re}$, $v=-\cos x\sin y\,e^{-2t/Re}$, $p=\tfrac14(\cos2x+\cos2y)e^{-4t/Re}$. Checked by hand: divergence 0, $uu_x+vu_y=\tfrac12\sin2x=-p_x$, eigenvalue −2 ⇒ $e^{-2t/Re}$; same form as `MAC.taylor_green` docstring |
| 5 | bare equation numbers | **fixed** | re-scan of all markdown: every previously flagged citation now has its equation beside it ((10.29), (10.67)–(10.68), (10.120), (10.124), (10.160)–(10.162), (10.168)–(10.169), a representative (10.188)–(10.197) entry, slip #4) |
| 6 | honey / ocean claims | **fixed** | cell 511 "a light oil (ν = 10⁻⁴ m²/s, about 100 × water)"; cells 304/306 compare 500 m²/s with "tens to about a hundred m²/s at 10 km" — correct |

### Round-1 Should-fix: done
Symbols (Θ, α_Θ, β_Θ; α_t, β_t; Ma in (10.151)–(10.155); $2\mathbf D[\mathbf u]$ in D22 step 5), slips renamed "#1…#12", explainer
numbering, phase speed defined (cell 368, $c_{num}/u=-\arg G/(C\theta)$), θ-method reminder (N29), GLS expanded, Lamb wave,
`lstsq`, `np.block`/`eigvalsh`, Euler's V − E + T, `tripcolor`, `re.findall`, float32 floor numbers, D10 order ≈ 0.8, figure fixes
(cell 87 ticks, cell 398 mathtext titles, cell 481 panels (b)/(c) with colour bars, cell 518 MAC 32² now dotted and visible),
N86 Re = 100 lift history plotted (cell 531, from `block_forces_re100_dx16.csv`), cavity-convergence claim now computed
(cell 558: 64²→128² change 0.204 % vs 64² deviation 0.258 % of U).

### New or remaining (all Should fix; none blocks)
1. **Cell 557: a broken LaTeX escape.** The source contains a literal TAB: `$64^2<TAB>o128^2$` (`\to` was read as `\t` + `o`),
   so the page shows "64² o128²". Write the builder string as a raw string (`r"""…"""`) or `\to`; a builder self-check for
   control characters in markdown would catch every future case.
2. Cell 499 "the centreline cuts the primary eddy through its core" — the eddy centre is at x ≈ 0.62 (and so drawn); say
   "just left of its core".
3. Cell 532 "the lift oscillates about zero" — the Re = 100 block lift in cell 531 swings about ≈ +0.1 (confined, coarse grid);
   say "about a small positive mean" or quote the mean.
4. Templated code comments ("# show the numbers computed above", "# the panel's title: what this panel shows") remain in many
   cells — acknowledged by the builder; keep on the list for the next chapter's builder template.

## Verdict (round 2): PASS
All six Must-fix items are verified fixed in the executed notebook; no new Must-fix found. Item 1 above is a one-character
rendering fix that should go in before publishing.
