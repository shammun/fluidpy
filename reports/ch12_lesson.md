# Chapter 12 — lesson review, round 3 (final)                      2026-10-07

## Verdict: PASS — 0 Must-fix · 9 Should-fix carried

coverage_check: **OK** (0 errors, 0 warnings) · sections covered **13/13** (12.1–12.13) · CORE blocks **16/16** with
plain words + idea + derivation + tiny example + commented code + from-scratch assert + visual · notebook: 683 cells
(520 markdown, 163 code), 0 error outputs, 0 stderr outputs, 10 explainers embedded (one `show_viz` each) ·
derivations **28/28** (round-1 audit; the one step that was open, D13 step 5, re-read and compared with the page) ·
`tools/check_public.py`: OK (647 files).

Reviewed: `outputs/ch12/executed.ipynb`. Scratch files of this round, all in `outputs/ch12/`:
`_rev3_scan.py` → `_rev3_scan.txt` and `_rev3_scan2.py` (control characters, TeX commands without a backslash,
unbalanced `$` and braces), `_rev3_teeth.py` (the scanner run on planted damage), `_rev3_chk.py`, `_rev3_chk2.py` (the
spectrum numbers below), `_rev3_diff.py` → `_rev3_diff.txt` (every cell that differs from round 2), `_rev3_wdc.py`
(code cells and their explanation blocks), `_rev3_sample.py` → `_rev3_sample.txt` (the ten sampled cells);
`review_dump.txt`, `review_img/` regenerated (round-2 copies kept as `review_dump_r2.txt`, `review_img_r2/`).

What changed since round 2 (`_rev3_diff.py`): 14 markdown cells and 15 code cells differ in source; 28 more differ in
output only (10 explainer frames now embedded, 14 `Legend at 0x…` addresses, 4 empty lines). Every one of the 29
source changes was read; none introduces a defect.

## Round-3 rulings on the two Must-fix items of round 2

### 1 · D13 step 5, the broken equation — **fixed**

Cell 359 now reads: "*Why we can do this:* Viscous terms are ~ 1/Re of the inertia terms; pressure and
$\overline{u^2}$ gradients were shown small in step 4. This is the first of (12.61),
$U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}\cong-\frac{\partial\overline{uv}}{\partial y}$."
Page 575 (rendered, `chapters/pages/ch12/p602.png`) prints (12.61) as
$U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}\cong-\frac{\partial\overline{uv}}{\partial y}$, and
$0\cong-\frac1\rho\frac{\partial}{\partial y}(P+\rho\overline{v^2})$ — the first member is what the cell shows.

Scans (counts):

| what was scanned | control chars | TeX name without `\` | escape tails (`rac{`, `ext{`, `heta`, `abla`, `ightarrow`, …) | odd `$` | unbalanced braces in a maths span | maths spans |
|---|---|---|---|---|---|---|
| `_rev2_ctrl.py` re-run, 520 markdown cells | 0 | 0 | 0 | — | — | — |
| mine: 520 markdown cells | 0 | 0 | 0 | 0 | 0 | 2824 |
| mine: 214 text outputs (streams, `text/plain`, plotly titles and trace names, small HTML) | 0 | 0 | 0 | 0 | 0 | 4 |
| mine: `notebooks/ch12_turbulence.ipynb`, 520 markdown cells | 0 | 0 | 0 | 0 | 0 | 2824 |
| mine: `analysis/ch12_design.md` (502 982 characters, line by line) | 0 | 0 | 0 | 0 | 0 | 1800 |
| mine: `analysis/ch12_curation.md` | 0 | 0 | 0 | 0 | 0 | 69 |
| mine: `notebooks/build_ch12.py` | 0 | — | 5, all inside the builder's own guard pattern `_EATEN` (lines 5332, 5337) | — | — | — |
| supplement: tabs and carriage returns in markdown, outputs and the design file; short tails (`au`, `ho`, `imes`, `ilde`, `ar`, `eq`, …) in 4624 maths spans | 0 tabs, 0 CR | — | 2 hits, both the legitimate product $a\,u_*$ in $\ln(au_*/\nu)$ (cell 493 and its design line) | — | — | — |

"Control chars" covers U+0000–0008, 000B, 000C (form feed), 000E–001F, 007F, 0080–009F, FFFD, 2028, 2029. The scanner
was tested on planted damage before its zeros were trusted (`_rev3_teeth.py`): the round-2 line with three form feeds
gives 3 control + 3 tail hits, the same line with the form feeds turned into spaces 3 tail hits, `\text` → ` ext{` 1,
`\theta` → ` heta` 1, `\nabla` → ` abla` 1, `\rightarrow` → ` ightarrow` 1, bare `overline{`, `partial`, `mathrm{`,
`left(`/`right)` 5, an unclosed `$` 1, and the correct line 0. One planted case it does not catch by pattern — a tab
that ate the `t` of `\tau` — is covered by the tab count (0). The builder now carries its own guard
(`self_check_eaten`), which is the right place for it.

### 2 · The two routes to $S_{11}(k_1)$ (C05, N77; cells 255–257) — **fixed**

The cell (256) now computes the periodic correlation,
`np.fft.ifft(np.abs(np.fft.fft(row))**2).real/n3`, transforms it with a plain cosine sum over one period,
`np.sum(R_circ*np.cos(k_*r_circ))*dx3/(2*np.pi)`, and draws one bold line over the whole range; a dotted curve shows
the biased estimate of a non-periodic record and a thin dashed guide $R_{11}(0)/(\pi Lk_1^2)$. Printed: "periodic
correlation route / periodogram over the plotted range: between 1.0000 and 1.0000" and "biased estimate at k1 = 10
and 30 rad/m: 1.3e-03, 6.4e-05; R(0)/(π L k²) there: 4.5e-04, 5.0e-05 m³/s²". Title: "The two routes to $S_{11}(k_1)$
give the same spectrum". The figure shows exactly that: circles and line together from $10^{-1}$ down to $10^{-12}$,
the dotted curve leaving them near $k_1\approx8$ and running just above the dashed guide.

My own numbers (`_rev3_chk.py`; nothing taken from the notebook except the definition of the field, 586 lines of 64
points, $L=2\pi$ m):

| check | result |
|---|---|
| circular correlation by direct wrapped lag products (no FFT) = `ifft(|fft|²)/n` | equal to 1e-12; $R_{11}(0)$ = 0.8856 m²/s² |
| my periodogram $\lvert\sum u\,e^{-ikx}\Delta x\rvert^2/(2\pi L)$ (slow sum = `np.fft`) against the cosine sum of the correlation, $k_1$ = 1 … 21 rad/m (the plotted range, $S>10^{-12}$) | ratio between 1.000000 and 1.000032; the 3 × 10⁻⁵ is at $S$ = 2.8 × 10⁻¹², an absolute difference of 10⁻¹⁶ — round-off, as the text says |
| library `TS.periodogram` against mine, same range | equal to 2 × 10⁻¹³ |
| normalisation: $\sum S\,\Delta k$ over both signs of $k$ | 0.885595 = the variance of the lines (0.885595): two-sided, $\int S_{11}\,dk_1=\overline{u_1^2}$, with the $1/2\pi$ in the forward transform — the convention of $S_{11}(k_1)=\frac1{2\pi}\int R_{11}(r_1)e^{-ik_1r_1}dr_1$ (12.45) in N76 (compared with p. 563), of `TS.spectrum_from_correlation` and of C03 |
| biased estimate ($n-m$ products divided by $n$): mine = library | equal; and it differs from $(1-r/L)\times$ the circular correlation by at most 0.008 (of 0.886) |
| the tail | see the table below |

| $k_1$ [rad/m] | true $S_{11}$ | biased estimate | guide $R_{11}(0)/(\pi Lk_1^2)$ | biased / guide | (biased − true) / guide |
|---|---|---|---|---|---|
| 8 | 3.95e-3 | 4.79e-3 | 7.01e-4 | 6.83 | 1.19 |
| 10 | 5.52e-4 | 1.29e-3 | 4.49e-4 | 2.88 | 1.65 |
| 14 | 2.51e-6 | 3.31e-4 | 2.29e-4 | 1.44 | 1.43 |
| 20 | 3.01e-11 | 1.48e-4 | 1.12e-4 | 1.32 | 1.32 |
| 30 | 2.54e-23 | 6.42e-5 | 4.98e-5 | 1.29 | 1.29 |

Log–log slope of the biased estimate for $k_1\ge14$: −2.2. For a smooth model correlation times the triangle,
integrated exactly, the transform over $R(0)/(\pi Lk^2)$ is 1.03 at $k=30$ and 1.002 at $k=100$: the law and its
coefficient are right; the 30 % excess in the notebook's curve is the 64-point grid. The printed numbers (1.3e-03,
6.4e-05, 4.5e-04, 5.0e-05) are reproduced.

Answers to the four questions:

- **First-use rule.** Met. The Wiener–Khinchin statement is N38 in C03 (cell 139, "the spectrum is also the squared
  magnitude of the Fourier transform of the record itself, $S_e(\omega)=\lim_{T\to\infty}\frac1{2\pi T}\lvert\int u\,e^{-i\omega t}dt\rvert^2$",
  with the triangular weight named). The FFT form is primer P289 (cell 99–100: "transform the record, multiply by its
  complex conjugate, transform back. The record must first be padded with zeros …, otherwise the end wraps round onto
  the beginning"), with a five-line demo. Cell 256's comment says what is new here: "the PERIODIC (circular)
  correlation — Wiener–Khinchin: the inverse FFT of |FFT|² is the correlation of a periodic record, exactly"; cell 257
  names "the Wiener–Khinchin theorem of C03". What is missing is a sentence that links back to P289 (Should-fix 2).
- **Reason for the tail.** Now correct. Cell 257: "Its usual 'biased' correlation estimate divides the sum of the
  $n-m$ available products by $n$, which multiplies the true correlation by the triangle $1-r/L$. A triangle has a
  kink at $r=0$, and a kink transforms into a tail $R_{11}(0)/(\pi Lk_1^2)$ that falls only like $k_1^{-2}$". A slope
  jump of $-2R(0)/L$ at zero lag does transform to $R(0)/(\pi Lk^2)$ with this normalisation (checked by integration,
  above). "A longer record ($L$ larger) lowers it" follows from the $1/L$.
- **Honesty of the size.** Acceptable, with one loose number. The text says the dotted curve "levels off along the
  thin dashed line" and is "still about $5\times10^{-5}$ at $k_1=30$ rad/m … (the cell prints both)". The cell prints
  6.4e-05 for the curve and 5.0e-05 for the guide: "about 5 × 10⁻⁵" is the guide's value given to the curve, 28 % low.
  The text makes no claim at $k_1=10$, where the curve (1.3e-3) is the true spectrum (5.5e-4) plus the tail and is 2.9
  times the guide. Nothing false is asserted and both numbers are on screen; the wording should say "within about
  30 %" (Should-fix 1).
- **Normalisation.** Consistent (table above).

## Must fix

None.

## Regression check on the other changes of this round

| What (cell) | Ruling |
|---|---|
| C02 ASCII sketch (81): "signal u(t)", "copy u(t + τ)", "lag τ = 0", "lag τ large" | ok — no line starts with a bare symbol; the columns still line up |
| Surface-layer sketch (589) and its notes (590) | ok — the lower label now reads "forced convection: shear-made turbulence, near-logarithmic wind", clear of the line "z = ∣L_M∣: buoyancy and shear production equal in size" (its first letters touch the wind curve; readable). The sketch is titled "unstable surface layer, L_M = −20 m", for which "forced convection" is the right name. For either sign, cell 590 says: "it says 'forced convection' for any height well below $\lvert L_M\rvert$ **whatever the sign** of $L_M$ — the layer there is shear-driven; whether the stratification is stable or unstable must be read from the sign of $L_M$". "Shear-made turbulence, near-logarithmic wind" is true for $z\ll\lvert L_M\rvert$ of either sign, and the printed list shows it ((10 m, +100 m) → "forced convection") |
| N95 (322) | ok — "the energy spectrum $E(K)$ we have used since §12.1 (the book's letter for it is $S(K)$)", then $\bar e=\int_0^\infty E(K)\,dK$ and $E(K)=C\bar\varepsilon^{2/3}K^{-5/3}$ |
| R01 (161), the two forms of (4.86) | ok — p. 136 (rendered) prints $\frac{D\mathbf u}{Dt}=-\frac1{\rho_0}\nabla p'+\frac{\rho'}{\rho_0}\mathbf g+\nu\nabla^2\mathbf u$ (4.86); §12.5 reprints the index form with full pressure under the same number. The new sentence: "put $\tilde p=p_0(z)+p'$ with $dp_0/dz=-\rho_0g$ and $\rho'=-\rho_0\alpha(\tilde T-T_0)$ into the form above and the two are the same equation". By hand: $-\frac1{\rho_0}\partial_i\tilde p-g[1-\alpha(\tilde T-T_0)]\delta_{i3}=-\frac1{\rho_0}\partial_ip'+g\delta_{i3}-g\delta_{i3}+g\alpha(\tilde T-T_0)\delta_{i3}$, and $(\rho'/\rho_0)\mathbf g=+g\alpha(\tilde T-T_0)\hat{\mathbf z}$ ✓ |
| "Γ < Γa" explained beside each printed string | ok — R05 (565): "read it as $\Gamma_{met}<\Gamma_d$ — its 'Γ' is our $\Gamma_{met}\equiv-dT/dz$ and its 'Γa' the positive dry-adiabatic rate $\Gamma_d\approx+9.8$ K/km"; C14 (579): "means $\Gamma_{met}<\Gamma_d$: here $-10<9.8$ K/km" (the output line is "−10.0 < 9.8"); C15 (598): "read the printed 'Γ < Γa' as $\Gamma_{met}<\Gamma_d$" |
| Both lapse-rate conventions together | ok — every stratification mention has both: the title-cell table (Kundu $\Gamma\equiv dT/dz$, $\Gamma_a\approx-9.8$ K/km, stable when $dT/dz>\Gamma_a$; meteorology $\Gamma_{met}\equiv-dT/dz$, $\Gamma_d\approx+9.8$ K/km, stable when $\Gamma_{met}<\Gamma_d$), C06's mixing example and its legend, R05's table and output, D24's check and trap, C14's tiny example and outputs, the three bar-group labels, C15's verdicts, the slider trace names and the widget title. The code computes with Kundu's sign throughout (`Gamma_a = STRAT.adiabatic_lapse_rate()` = −9.76 K/km) |
| C15, κ = 0.4 (596) | ok — "$\kappa=0.4$ (the round value customary in micrometeorology; 0.41 elsewhere in this chapter — a 2.5 % difference in every wind speed here)"; 0.41/0.40 = 1.025 ✓ |
| "the printed B is our evaluation" (465–467; the cell was later moved to an illustrative κ = 0.40 so its output coincides with no tabulated pair) | ok — the equation stands directly above the cell, $\kappa B=1.6[\exp(0.1663B)-1]$ (12.92), and again in item 3 of the explanation; it matches p. 590. The value is computed (`B_nc = WT.nagib_chauhan_B(0.384)`, a `brentq` root) and printed with an f-string; by hand both sides of the relation equal 1.601 ✓. "The printed $B$ is **our evaluation** of that formula at the DNS value $\kappa=0.384$, not a number taken from a table" |
| N168 (514) | ok — "(12.99) (the Start and step 1)" |
| New code comments | ok — 26 outside cell 256 (cells 374, 376, 395, 403, 469, 472, 474, 521, 557, 582, 602, 645, 647, 657) and 13 in cell 256. Each says what the line is for, e.g. "the rows of this DNS case", "D21's profile; van Driest damping only when A⁺ > 0", "an all-NaN curve: plotly draws nothing for it", "three times: early, at the memory time, late", "squared end-to-end distance of the library's walkers". None restates syntax |
| Frozen-field sentence (151), slider note (156) | ok — "the error is 0.07 over one eddy (0.125 m) and 0.35 over the whole 1 m line" equals the output; "corner frequency × $\Lambda_t$ = 1 at every slider value" |
| Animations A2, A3, A5 | unchanged from round 2 (not re-judged) |

## Spot check — ten markdown cells drawn at random (seed 20261007, from 171 cells longer than 350 characters and at least three cells away from anything quoted in rounds 1–2)

| cell | block | what it is | equation beside every number | nothing used before explained | own words, no book table/example numbers |
|---|---|---|---|---|---|
| 81 | C02 | the idea, ASCII sketch, symbols with units | no number cited | ok | ok |
| 102 | C02 | N32, correlation time and independent samples | no number cited | ok — P281 named | ok |
| 126 | C03 | tools from earlier chapters | no number cited | ok — each tool with chapter and primer id | ok |
| 207 | C04 | N59, the Reynolds-stress transport equation (12.35) with a term table | ok — (12.35) in full, compared with p. 560: every sign and index agrees; (12.47) written out | ok | ok |
| 243 | C05 | N73, $\bar\varepsilon=2\nu\overline{S'_{ij}S'_{ij}}=\frac\nu2\overline{(\partial u_i/\partial x_j+\partial u_j/\partial x_i)^2}$ (12.42) | ok — matches p. 563; $\varepsilon=2\nu S_{ij}S_{ij}$ (4.58) written | ok | ok |
| 247 | C05 | traps in D08 | ok — the bracket form of (12.43) written | ok | ok |
| 314 | C08 | the idea, table of the three ranges | no number cited | ok — units given for each quantity | ok |
| 544 | C13 | N177, the modelled dissipation equation (12.105) | ok — in full | ok — term-by-term reading | ok |
| 626 | C16 | the idea, two limits, symbols, warning about α | no number cited | ok | ok |
| 654 | C16 | N209, the random walk, $\overline{R_n^2}=\overline{R_{n-1}^2}+L^2+2\overline{\mathbf R_{n-1}\cdot\mathbf L}$ (12.124) | ok — matches p. 605 | ok — "the product rule of C01" | ok |

Ten of ten pass. The round-2 state held.

## Should fix (carried; none blocks)

1. **Cell 257, the size of the tail.** Replace "levels off along the thin dashed line" by "bends onto a line of slope
   −2, within about 30 % of the dashed guide for $k_1\ge14$ rad/m", and "still about $5\times10^{-5}$ at $k_1=30$" by
   "still $6.4\times10^{-5}$ at $k_1=30$ rad/m (guide $5.0\times10^{-5}$)"; add that at $k_1=10$ the curve,
   $1.3\times10^{-3}$, is the true spectrum $5.5\times10^{-4}$ plus the tail. In cell 256 the comment "the tail it must
   have" → "the tail the kink predicts (leading order)".
2. **Cell 256 has no "What does the code above do?" block.** It is the longest non-derivation cell of the chapter
   (26 statements, three estimators). Add four items: the lines and their means; the circular correlation, with the
   link to primer P289 ("there we padded with zeros to stop the wrap-round; here the record is periodic, so the
   wrap-round is the correct lag product"); the two transforms; the biased estimate and the guide. Break the 230-
   character comment on the `R_circ` line. Correction to round 2: that report said every non-trivial cell is followed
   by the block. The count (`_rev3_wdc.py`) is 90 of 163; of the 73 without it, 10 are `show_viz` one-liners, 32 are
   figure cells followed by *What you see / How to read it / What would change if* (cell 256 is the only one of these
   whose mechanics are not obvious from its comments), and 31 are primer or note demos of 1–9 lines and four sympy
   check cells with every statement commented. Every CORE block has the block after its library cell and after its
   from-scratch cell.
3. **Sketch labels.** Jet sketch (cell 348): caption texts crossed by the profile lines; parcel sketch (cell 187): the
   two stress labels touch. Both images are byte-identical to round 2; both readable.
4. **`WT.boundary_layer_stress_from_profile`** (cell 427–429) is named and not run; one printed number (total
   stress/τ₀ at y⁺ = 30 and at y/δ = 0.1) would back "nearly constant across the inner layer".
5. **`ch12.book_slips()`** (shown in cell 7): row 5 "the power family … is the valid one" → "a valid one; so is
   δ ~ e^{ax}, U_CL ~ e^{−ax/2}, Ψ = const" (as the box in cell 394 says); add rows for (12.16) printed without the
   absolute value and for "V → 0 at the jet edge".
6. **The library string "Γ < Γa"** in the meteorological line of `STRAT.lapse_rate_stability` should read
   "Γ_met < Γ_d" (a core function — for the implementer). The notebook explains the string at each of its three
   appearances, so the rule is met.
7. **Generic comment phrases** remain on 113 lines: "a new figure and its panel(s); figsize in inches", "the legend
   names every curve", "hand the result back to the caller", "a thin reference line" (say which line), "move the
   drawn objects to this frame's values".
8. **Stray output lines.** 18 `<matplotlib.legend.Legend at 0x…>` (or similar) lines under figures and a few empty
   stream outputs: end those cells with a semicolon or assign the last call.
9. **`notebooks/ch12_turbulence.ipynb` carries no outputs** (0 outputs in 163 code cells; its markdown is identical to
   the executed copy). For the publish phase: the committed notebook is meant to be the executed one.

## Unexplained-first-use list

Round-1 table re-walked in round 2: no row missing. New first uses in this round:

| term | first cell | explained | status |
|---|---|---|---|
| periodic (circular) correlation | 256 | 256 (comment), 257; wrap-round in P289, cell 99 | ok (link sentence: Should-fix 2) |
| `np.fft.fft` / `np.fft.ifft` (complex, full) | 256 | P286 (`ifft2`, `ifftn`), P289 and P291 (`rfft`, `irfft`), cells 13, 99, 136 | ok |
| biased / unbiased correlation estimate | 256 | P289, cell 99 ("'Unbiased' divides lag $m$ by $N-m$"); `TS.autocorrelation` used since C02 | ok |
| triangle weight $1-r/L$ | 257 | 257; the same triangle in N38 (cell 139) and C16 | ok |

No row is missing an explanation.

## Derivation audit

All 28 derivations: every step follows (round 1), why / in-words ok, checks ok, ★★★ with sympy checks that test the
derived result (D07 both directions since round 2). The single open item of round 2 — D13 step 5, an unreadable
equation — is closed (ruling 1). No derivation text other than that line changed in this round.

| D | verdict |
|---|---|
| D01–D12, D14–D28 | ok (unchanged since round 2) |
| D13 | ok — step 2 cites the second of (12.61), $0\cong-\frac1\rho\frac{\partial}{\partial y}(P+\rho\overline{v^2})$; step 5 the first, $U\frac{\partial U}{\partial x}+V\frac{\partial U}{\partial y}\cong-\frac{\partial\overline{uv}}{\partial y}$; both as on p. 575 |

## Three best blocks, three weakest

Best:
1. **C06, the energy budget** — the order-of-magnitude estimate, a "Where this holds" box and four computed shares
   (0.56, 0.45, 0.37, 0.33) that show where the estimate fails; then the stratified mixing example in both lapse-rate
   conventions.
2. **C11, the log law against DNS** — the gap is printed before it is described, κ is fitted in a stated window, and
   the text says which constant the data favour and by how much.
3. **C09, the plane jet** — similarity derivation, a three-station check of the invariant, and a slip box that states
   what is printed, why it fails, what works, with two sympy lines for both.

Weakest (all passing):
1. **C05, the spectrum cell (256–257)** — now right, but the densest cell of the chapter has no numbered explanation
   and its reading note rounds the tail to the guide's value (Should-fix 1, 2).
2. **C10, the stress across the wall layer (427–430)** — "nearly constant across the inner layer" rests on the
   equation alone; the function that would show it is named and not run (Should-fix 4).
3. **C09's opening sketch (348)** and **C04's parcel sketch (187)** — labels crossed or touching; the two figures a
   reader meets first in those blocks are the least tidy (Should-fix 3).

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **An artefact shown beside the clean result, with its predicted size**: the exact route as the bold line, the
  everyday estimator as a dotted curve, a one-term guide for its tail, and both numbers printed.
- **"Where this holds" under an order-of-magnitude estimate** (N80).
- **A companion slider on the right axes** (F2, F8) and printed numbers for what a log–log plot cannot show.
- **A printed gap before a comparative sentence** (cell 469).
- **A slip box that says what is wrong, what would be right, and checks both** (slip #5).
- **A library string that uses another notation is translated at every appearance** ("read 'Γ < Γa' as
  $\Gamma_{met}<\Gamma_d$").
- **For the skill (lessons of rounds 2–3):** never write LaTeX through a shell heredoc or a non-raw string; let the
  builder scan its own output for command tails without a backslash and for control characters, and test that scan on
  planted damage; before explaining a numerical artefact, switch the suspected cause off and look; when a guide line
  is drawn, state how closely the curve follows it.

## History

| round | cells | verdict | Must-fix | Should-fix | notes |
|---|---|---|---|---|---|
| 1 | 681 | FAIL | 12 | 17 | no derivation step failed; no book equation with a wrong sign or exponent |
| 2 | 683 | FAIL | 2 (both new; 12 of 12 from round 1 closed) | 12 | one equation broken by form feeds; one wrong explanation of a spectrum tail |
| 3 | 683 | **PASS** | 0 (2 of 2 closed) | 9 | 8 of the 12 round-2 Should-fix items done (surface-layer label, κ = 0.4 clause, N95, the (4.86) bridge, bare code lines, two numbers in prose, "our evaluation" of (12.92), N168); 5 carried, 4 new and minor |

## Verdict: PASS
