# Chapter 16 — explainer review                                  2026-10-14

**Resolution (orchestrator, 2026-10-14): PASS, 10/10.** The four Must-fix strings (single backslashes in ordinary JS
strings eating TeX) are fixed by doubling the backslashes: `pulse_wave_artery.html:2932` (16.63) note,
`womersley_profiles.html:2737` and `:2747` (D13 steps 1 and 4, (16.72) and (16.77)), `junction_reflection.html:2727`
(Explain §1, (16.151)). A scan of every ch15 and ch16 explainer's own script for lone-backslash TeX commands in untagged
strings finds no other case. viz_lint ok and shot.py PASS on all three (216 / 304 / 416 views); the D13 step 4 page
(`desktop__derive-d2p4.png`) now renders (16.77) with no KaTeX error. Should-fix items stay open (non-blocking).

Fresh, read-only review of the ten explainers in `viz/ch16/` against `analysis/ch16_design.md` (E1–E10), the curation, the
corrections in `reports/ch16_review.md` (S34 Dean sign, collapsible law α ≤ 1, throat α_t = 0.35, Q* = 47.646 mL/s, Starling
Q_choke = 19.99 mL/s at −946.1 Pa) and `notebooks/build_ch16.py --dump`.

**Runs (real numbers).**
- `tools/viz_lint.py --chapter ch16`: 10/10 `ok`.
- `tools/shot.py --chapter ch16` (full, 8 sizes, every tab, every tour step, every derivation page): **10/10 PASS**, no time-outs,
  KaTeX loaded at every size: casson 408 views / 22 rows · collapsible 368 / 33 · dean 224 / 25 · junction 216 / 20 · murray 320 / 26 ·
  pulse 304 / 19 · windkessel 296 / 17 · womersley_flow_phase 232 / 19 · womersley_profiles 416 / 22 · xylem 144 / 26.
  All 229 selftest rows `ok` (177 fluidpy parity rows + 52 invariants).
- `tools/eq_refs.py --chapter ch16`: 34 hits, judged one by one: 25 are `ref:` labels of derivations/code/equations or derivation
  titles that the notebook uses verbatim (accepted, as in ch14/ch15); 6 are *did* titles copied from the notebook ("Substitute
  into (16.55)", "Average (16.49) over one beat", …) whose *why* or "we had" line writes the equation out (accepted); 3 are in
  check/answer strings that carry the TeX via a constant (`EQ178`, `EQ182`) — false positives. No bare number in tour, Explain or quiz
  text. Three equations that *are* written out are rendered broken (Must fix below).
- Reviewer probes (scratchpad, not committed): real clicks on every preset (exactly one chip pressed after each, 43/43), the
  walkthrough driven by real "Next →" clicks (every step visited, 0 skipped, 10/10), 9 clicks + 1 drag on every visible view and
  every preset applied with the Explain tab open — **0 page errors** in all ten files; derivation step counts and `did` titles
  compared with `build_ch16.py --dump`: **all 34 counts equal**; title differences only cosmetic (casson D31 s2 / D32 s3 and dean D27 s1/s9 use
  TeX `$\Delta x$`/`$\kappa^0$` where the notebook has Unicode; collapsible D30 s5/s8/s9/s11 append the equation itself to the
  notebook's title — an improvement); a deep scan of every config string and live function output for broken TeX escapes and a
  `.katex-error` count over all tabs, steps and derivation pages (found the three Must-fix strings; 1 visible KaTeX error).
- sympy, on the derivation lines: Dean ψ₀ (∇⁴ψ₀ = 8r(1 − r²) sin α, f(1) = f′(1) = 0), w₁ residual 0 with w₁(1) = 0, ∫r f g dr =
  1541/1045094400 → 1541/50400; Casson ring integral 4I and the plug polynomial reproduce D33 steps 7–9 exactly; D30 steps 8–13
  and Murray D22 steps 4–8 checked by hand. All correct.

| slug | lint | shot (sizes/views) | parity rows ok | explain | derivations (steps ok) | teaching | polish | verdict |
|---|---|---|---|---|---|---|---|---|
| casson_blood_flow | ok | PASS 8 / 408 | 18/18 py + 4 inv | 9 sections, live, regime by ξ | D02 5/5 · D31 5/5 · D32 11/11 · D33 10/10 | good | good | **PASS** |
| windkessel_heartbeat | ok | PASS 8 / 296 | 14/14 + 3 | 8 sections, 7 live, RK regimes | D01 6/6 · D09 4/4 · D10 9/9 | good | good | **PASS** |
| pulse_wave_artery | ok | PASS 8 / 304 | 15/15 + 4 | 7 sections, live, E regimes + linearity | D08 6/6 · D11 7/7 · D12 7/7 | good | Equations note broken | **FAIL** |
| womersley_profiles | ok | PASS 8 / 416 | 16/16 + 6 | 8 sections, table, regime | D03 4/4 · D13 13/13 (2 *why* broken) · D14 6/6 · D15 8/8 | good | phone cycle plot flat | **FAIL** |
| womersley_flow_phase | ok | PASS 8 / 232 | 14/14 + 5 | 8 sections + regime table | D16 7/7 · D17 7/7 | very good | minor label overlap | **PASS** |
| murray_branching | ok | PASS 8 / 320 | 19/19 + 7 | 7 sections, regime | D22 8/8 · D23 3/3 · D24 11/11 | good | pumping/total colours close | **PASS** |
| junction_reflection | ok | PASS 8 / 216 | 15/15 + 5 | 6 sections, broken (16.151) | D25 6/6 · D26 6/6 | good | good | **FAIL** |
| dean_curved_pipe | ok | PASS 8 / 224 | 19/19 + 6 | 8 sections incl. S34 | D27 14/14 (sympy) | good | slice small on phones | **PASS** |
| collapsible_tube_choking | ok | PASS 8 / 368 | 27/27 + 6 | 7 sections, 3 modes | D28 4/4 · D29 9/9 · D30 15/15 | good | 17 Explore pages at 360×640 | **PASS** |
| xylem_water_potential | ok | PASS 8 / 144 | 20/20 + 6 | 8 sections, regime | D34 5/5 | good | good | **PASS** |

Tolerances are honest: closed forms at 1e-12…1e-10; 1e-9…1e-6 only where the two sides use different numerics (Starling
α-integration vs fluidpy bisection 1e-9, Weiszfeld vs Nelder–Mead 1e-6, grid echo `network_pulse_run` vs closed-form R 1e-6, fluidpy's
numerical dA₀/dx 1e-7); quadrature invariants 1e-6…1e-4. The corrections are all live and on screen: Dean (16.174) shown with
"+" and c₄ = 0.0119312 vs 0.0121130 (tour step 6, Explain §5, `flux(576, S34)` parity row); collapsible throat default α_t = 0.35,
Q* = 47.65 mL/s, H = −41.47 Pa, warning "outside the tube-law fit" for α_t > (1 − n/2)^{1/n}; Starling plateau 19.99 mL/s at
p₂* = −7.10 mmHg (parity row against `flow_limitation_curve` at 1e-9).

## casson_blood_flow
**Must fix** — none.

**Should fix**
1. Explain §6 in the no-flow branch says "the bracket of (16.220) are both zero" — write the bracket ($1-\tfrac{16}{7}\xi^{1/2}+\dots$) or
   reuse `EQ.e220i`.
2. `phone__explore.png`: at 360×640 the material view is 60 px tall and the pipe view has no title — acceptable, but the τ(γ̇) curve
   is barely readable; consider `hidePortrait` on `law` during Explore or a tighter y range.

**What works** — dye line released at t = 0 that shows plug vs Poiseuille ghost in one picture; stress panel sharing the r axis so the
crossing τ = τ_y *is* the plug edge; inspector traces (16.209) → (16.211) → (16.216) with the reader's radius; D33 fills the book's
"considerable algebra" with polynomials in t = √ξ (sympy-exact); trapezoid-flux invariant.

## windkessel_heartbeat
**Must fix** — none.

**Should fix**
1. D10 `result.tex` reads "p(t) from (16.52)" and D09 `start.plain` "The balance (16.47) in words" — write (16.52) and (16.47)
   out (the step pages already have them; the goal/result pages do not).
2. D01 `view: 'system'`; step 3's *watch* points at "the two vertical legs of the loop", which phones cannot see in the Derivation
   tab. Use `view: 'loop'` (storyboard) or reword the watch for phones.
3. `desktop__tour-step2.png`: the rigid-aorta ghost (p = QR, peaks 403 mmHg) is drawn as clipped vertical dashes that clutter the
   trace; draw it faint or only where it is on-axis.

**What works** — beat-by-beat closed form (no ODE solver) with parity at 1e-10; terms bars that add to the inflow every instant;
true mean R·Q̄ and the ⅓ rule as two lines; schematic PV loop built from the model's own aortic pressures (no book table).

## pulse_wave_artery
**Must fix**
1. `viz/ch16/pulse_wave_artery.html:2932` (Equations tab, `mk` note): the string is single-quoted with single backslashes —
   `'… it reads $\sqrt{Eh/(\rho d)}$; thick wall: $c=\sqrt{Eh/(2\rho(a_0+h/2))}$ (16.63)'`. In JS `\s` → "s" and `\r` → a carriage
   return, so the page shows "sqrt{Eh/( ho d)}" and the thick-wall (16.63) is printed as garbage text next to its number. Double
   the backslashes (`\\sqrt`, `\\rho`) or use the file's own `R\`…\`` raw strings.

**Should fix**
1. Preset "soft vein-like" (E 0.1 MPa, h/a₀ 0.05, a₀ 5 mm) with the default 25 mmHg pulse lands at p̂/(Eh/a₀) = 0.67, u/c = 1.3
   (status ⛔ "too strong for linear theory"). Either lower P in that preset (e.g. 5 mmHg) or say so in its title.

**What works** — tracer shifts integrated in closed form (erf) with a path-integral invariant; one exaggeration factor for bulge and
shifts, labelled with the true sizes; linearity status u/c = 2p̂/(Eh/a₀) shown as one condition; D11 shows the printed (16.54) slip.

## womersley_profiles
**Must fix**
1. `viz/ch16/womersley_profiles.html:2737` (D13 step 1 *why*): `'… $\mu[u_{rr}+u_r/r]-\rho u_t=p_x$ (16.72) …'` — `\m` → "m",
   `\r` → carriage return: the equation (16.72) renders as "mu[u_rr + u_r/r] − ho u_t = p_x".
2. `viz/ch16/womersley_profiles.html:2747` (D13 step 4 *why*): `'… $U=C_1J_0+C_2Y_0+\frac{G}{\omega\rho i}$ (16.77) …'` — `\f` is a
   form feed: KaTeX throws and the line is shown in red ("Unexpected character '\x0c'"), see `desktop__derive-d2p4.png`. Double every
   backslash in both strings (the rest of the file already does).

**Should fix**
1. `phone-tall__tour-step5.png` ("Reversal starts at the wall"): the cycle plot is ~120 px high and the push is scaled to the largest
   speed, so the purple (r = 0.9a) and blue (centre) traces are a few pixels tall — the step's whole point. On portrait use a y range
   fitted to the two speeds (the reversal dots are there; the curves are not readable).
2. D14 step 6 "Substitute into (16.82)": the *why* shows only the prefactor product; add $u=\frac{Ga^2}{i\mu\alpha^2}(1-F_1)e^{i\omega t}$ (16.82) (EQ82i).
3. `EQ88` (Explain §5, check Q1) mixes $G_n$, $\alpha_n$ with $\sqrt{\omega/2\nu}$; the Equations tab notes "use nω (slip S7)". Write
   the n = 1 form here or $\sqrt{n\omega/2\nu}$.

**What works** — full-diameter profile with the last half-cycle faint, both limits as ghosts and the Stokes band; the Hankel-scaled J₀
(no overflow to α = 50, parity rows on that branch); "reverses first" readout; vessel axis built from `vessel_scales_table`.

## womersley_flow_phase
**Must fix** — none.

**Should fix**
1. `phone-tall__tour-step2.png`: the in-plot x label "α (log scale)" sits on the purple φ curve; move it below the ticks or drop it on
   narrow views.

**What works** — the APS pattern done right: time traces, Bode pair and the Z = X + iY locus on one α; X/Y term bars with |Z| as total;
printed (16.96) vs corrected vs leading term with live errors (the printed one is worse than the leading term); wall-shear lead
36.7° read off arg F₂; ring-sum invariant 2∫F₁s ds = F₂.

## murray_branching
**Must fix** — none.

**Should fix**
1. Cost view: pumping (orange) and total (amber) are nearly the same colour; on phones (`phone-tall__tour-step2.png`) the total has
   no label. Make total clearly distinct (darker/thicker, or label it on small views).
2. Phone status "a −57 %: cost +887 % · B +2.0 %, gap 27 %" is cryptic; say "radius 57 % below a*" or drop the B part on phones.

**What works** — the 1 : 2 pumping : upkeep signature in the term bars; corrected vs printed (16.140) as a mode with units argument;
draggable branch point with a live "three pulls" triangle that closes at the Weiszfeld point (parity with Nelder–Mead 1e-6);
D24 fills the book's "making displacements" with a directional-derivative argument.

## junction_reflection
**Must fix**
1. `viz/ch16/junction_reflection.html:2727` (Explain §1, a template literal): `` `… $Q=Au=\pm Yp=\pm p/Z$ (16.151) …` `` — `\p` → "p",
   so (16.151) reads "Q=Au=pm Yp=pm p/Z". Write `\\pm` (or use `Y_TEX`, which is already correct).

**What works** — exact closed-form pulses (no grid) with a grid-echo parity row; energy bars that stack to 1 with R and T as signed
bars; inspector arithmetic of (16.159)–(16.160) at the probe; slips S16/S17 named where the equations are shown.

## dean_curved_pipe
**Must fix** — none.

**Should fix**
1. `phone-tall__tour-step2.png`: the slice is drawn ~160 px across inside a 260 px plot whose x range is ±2a, so the outward shift
   (0.11a ≈ 9 px) that step 2 is about is hard to see. On portrait fit the disc to the view height and put the inner/outer labels in the
   title.
2. Step 1 text is 37 words (storyboard cap 24); same for collapsible (33) and xylem (30).
3. Explain §0 says "Left / Right / Below" also on phones, where the views are stacked; use "Top / Middle" in portrait or neutral names.

**What works** — S34 made visible (tour step 6, Explain §5, both c₄ values, parity row); "which Dean number?" chips κ / K / D on one
slider; our reliability criterion stated as ours with its threshold K ≈ 460; true-scale tracers in a 3-D arch with a 2-D fallback;
D27 derivation is complete and sympy-exact.

## collapsible_tube_choking
**Must fix** — none.

**Should fix**
1. At 360×640 the Explore pager has 17 pages (one control or a third of a list per page). Hide `along` in Explore on small portrait or
   mark `xc` and the mode-irrelevant controls so they drop out.
2. `markMode()` reads the private flag `a._ready` (ch14 lesson: never use `app._…`); call `a.fit()` from `onChange` only, or test
   `a.root.isConnected`.
3. D29 `result.plain` names "(16.201) corrected" without the relation; write $2\tfrac{dc}c=\dots$ as in D30 step 8.

**What works** — three systems (Starling, smooth throat, giraffe jugular) on one area equation with term bars formed before dividing by
1 − S² (trap T10); the closed-form x(α) for the Starling tube (no nested root-finding) with parity to fluidpy at 1e-9; amber "news"
arcs that stop at S = 1 — choking shown, not told; the α ≤ 1 fit limit enforced and explained in the throat mode.

## xylem_water_potential
**Must fix** — none.

**Should fix**
1. `phone-tall__tour-step2.png`: the label "top p −0.78 MPa" is drawn half-way up the trunk; anchor it at the top station.

**What works** — ψ ledger as term bars at four stations; still vs moving column with the viscous drop shaded; tallest-tree H_max solved
and pinned by an invariant (p(H_max) = p_cav); S30 (book units) and S28 (ν) shown where used; phloem mode for (16.222).

## Candidate lessons for knowledge/viz_patterns.md
- (ch16) **A single backslash in a normal JS string silently eats TeX**: `\m`, `\p`, `\s` become letters, `\r`, `\f`, `\t`, `\b` become
  control characters (`\frac` → form feed). Lint and shot passed three files with broken equations. Use `String.raw` (`R\`…\``) for
  every TeX-bearing string, and let `shot.py` fail on any `.katex-error` and on control characters in config strings (scan
  `app.cfg` and the live text functions, as the reviewer's probe did).
- (ch16) Replace nested root-finding by a closed form whenever the physics allows it (Windkessel beat-by-beat, Starling tube x(α));
  parity rows then hold at 1e-9…1e-12 and the page stays fast.
- (ch16) Show a printed slip as a selectable mode with both numbers live (Murray (16.140), Womersley (16.96), Dean c₄/S34) and a parity
  row for each form.
- (ch16) A model's fit range belongs in the status line ("αₜ > 0.397: upstream α > 1, outside the tube-law fit") with the
  threshold computed, not typed.
- (ch16) Every preset should land inside the explainer's own validity range or say in its title that it does not (pulse "soft
  vein-like").

## Verdict: FAIL at review → PASS after the fix (see the resolution at the top)
Seven explainers PASS; **pulse_wave_artery, womersley_profiles and junction_reflection FAIL** on four broken-escape TeX strings
(one-line fixes each: double the backslashes at `pulse_wave_artery.html:2932`, `womersley_profiles.html:2737` and `:2747`,
`junction_reflection.html:2727`). After the fix, re-run `tools/shot.py` on the three files and confirm 0 `.katex-error` on
`desktop__derive-d2p4` of womersley_profiles.
