# ch01 + ch02 explainers: review of the "show every equation" text retrofit (commit a54f4c3), 2026-09-23

Scope: the 166 changed lines in `viz/ch01/*.html` and `viz/ch02/*.html` (`git diff a54f4c3~1 a54f4c3 -- viz/`). The reviewer did not edit any explainer.

## How it was checked
- **Mechanical.** `tools/viz_lint.py --chapter ch01/ch02`: 10/10 ok. `tools/shot.py --chapter ch01/ch02` (re-run today): 10/10 PASS, 8 sizes, every parity row ok (26/14/13/49/13 and 27/25/20/30/25 rows).
- **Sub-pages that shot.py never pages through.** A scratch Playwright driver (`<scratchpad>/subpages.py`, not in `tools/`) ran 3,891 views: 4 sizes (360×640, 390×844, 1920×1080, 1000×700) × every tab × every pager page × modes and presets:
  - rotation `vector / tensor+fixed / polar`
  - viscosity `momentum / heat / species`
  - heat_work `rev / stir / free`
  - parcel `atm / ocean / lab`
  - gauss `default / → point / source`
  - stokes `default / enclose the core`
  - cauchy `default / nonsym`
  - every walkthrough step and every Derivation page, with Check answers revealed.

  For each view it checks:
  - KaTeX loaded;
  - `.katex-error` elements;
  - visible text holding `$`, `\cmd`, `^{` or `_{` outside KaTeX (raw TeX);
  - KaTeX boxes wider than their card;
  - `VIZ.audit()` overflow.

  Screenshots are in `reports/viz/eq_retrofit/` (+ `results.json`). Targeted shots of every changed string that sits on a sub-page are in `reports/viz/eq_retrofit/located/`.
- **Before/after.** The same driver ran on the pre-retrofit files (`git show a54f4c3~1:…`) at 360×640, to separate retrofit effects from pre-existing defects.
- **TeX against the book.** Each changed equation was compared with `analysis/ch01.md` / `analysis/ch02.md`. Pages were rendered for (2.28)–(2.29) and (1.21).

## Results in one table
| explainer | TeX correct | renders (KaTeX, 0 errors, no raw $) | changed strings fit | retrofit verdict |
|---|---|---|---|---|
| ch01 buckingham_pi_machine | yes | yes | yes (Explain §6 p.11/12 phone; D28 result p.10/10 phone) | **PASS** |
| ch01 continuum_averaging_volume | yes (symbol renamed, S1) | yes | note fits; the live lines on the same page are clipped (P2) | **PASS** (P2 is pre-existing) |
| ch01 heat_work_paths | yes | yes | yes (stir: p.5–6/9 phone) | **PASS** |
| ch01 parcel_stability | yes | yes | tour quote clipping +23 px vs before (L1) | **PASS** (L1) |
| ch01 viscosity_momentum_diffusion | yes | yes | yes (heat/species: Explain §0 p.1; §6 p.8/11 phone) | **PASS** |
| ch02 cauchy_traction_principal_axes | yes | yes | yes | **PASS** |
| ch02 gauss_flux_box | yes (S2) | yes | **no**: Explore "Right now" +43 px clipped at 360 | **FAIL** (R1) |
| ch02 rotation_of_axes | yes | yes | yes; the fixed-array status wraps to 2 lines at 360, no overflow | **PASS** |
| ch02 stokes_circulation_loop | yes (S3) | yes | yes | **PASS** |
| ch02 strain_vs_rotation_split | **no**: (2.28)/(2.29) mislabelled | yes | yes | **FAIL** (R2) |

Across 3,891 views: KaTeX loaded everywhere, 0 `.katex-error`, and none of the retrofit strings left raw `$…$` or halved backslashes. All "raw" hits are pre-existing plain-text braces (S5).

## Must fix: caused by the retrofit
**R1 · gauss_flux_box, line 2237 (`regimeName.small`).**
- **Problem.** The label now carries the whole of $\nabla\cdot\mathbf Q=\lim_{V\to0}\frac1V\iint_A\mathbf n\cdot\mathbf Q\,dA$ (2.32). It is used in the Explore/notes "Right now (…)" heading, and the `meaning` paragraph right below already shows the same equation, so it appears twice.
- **Effect.** At 360×640 the "Right now" item overflows its page by 163 px (was 120 px), and the face table is cut off. Screenshot: `located/phone__10__gauss_flux_box__explore__pg7.png`.
- **Fix.** Revert to `small: 'small box, close to the limit V → 0'`, with no equation number. The equation stays visible in the `meaning` line at 2628 and in Explain §6.

**R2 · strain_vs_rotation_split, lines 2625 (Explain §6 title) and 2556 (term label `SA`).**
- **Problem.** The book (p. 55) has:
  - (2.28): $P=\tau_{kl}B_{kl}=\tau_{kl}S_{kl}+\tau_{kl}A_{kl}$
  - (2.29): $P=\tau_{kl}S_{kl}-\tau_{kl}A_{kl}$

  $\tau_{ij}A_{ij}=0$ and $\tau_{ij}B_{ij}=\tau_{ij}S_{ij}$ are *consequences* of the two together. The retrofit labels $\tau_{ij}B_{ij}=\tau_{ij}S_{ij}$ as "(2.28)" and $\tau_{ij}A_{ij}=0$ as "(2.29)", so the TeX does not match the numbered equations.
- **Fix, §6 title.** `'A symmetric tensor cannot see A: $\\tau_{ij}B_{ij} = \\tau_{ij}S_{ij} \\pm \\tau_{ij}A_{ij}$ (2.28)–(2.29)'`. If that is too long, use `'… — (2.28) minus (2.29) gives $\\tau_{ij}A_{ij} = 0$'` and write both equations out in the section's first `W.say`.
- **Fix, term label.** `'2 S:A — $\\tau_{ij}A_{ij} = 0$ by (2.28)–(2.29)'`.
- Also re-word the selftest names `τ:A = 0 (2.29)` and `τ:G = τ:S (2.28)` (lines 2819–2820). They are exempt from the rule but carry the same mislabel.

## Must fix: pre-existing, found on the sub-pages shot.py does not page through
**L1 · all 10 explainers (library `Pager` + `tools/shot.py`).**
- **Problem.** At 360×640, and in a few cases at 390×844, one paged item can be taller than the pager body. `Pager.layout` cannot split it, so it is clipped. The worst cases:
  - the walkthrough "Derivation · step k" quote cards: the *why* box is cut;
  - Explain §0 paragraphs;
  - Explore "Right now" tables;
  - Code blocks.
- **How bad.** Maximum clip per explainer (px, pre → post retrofit):
  - buckingham tour: 84 → 84
  - cauchy tour: 151 → 151
  - continuum tour: 102 → 102
  - gauss tour: 206 → 209
  - gauss Explore: 139 → 139
  - heat tour: 95 → 95
  - parcel tour: 102 → 125
  - rotation tour: 137 → 139
  - stokes tour: 198 → 200
  - strain tour: 142 → 142
  - viscosity: 24–32 → unchanged

  The retrofit made it worse only in parcel (+23 px, lengthened D19/D21 *why* strings quoted in the tour) and in gauss (R1).
- **Examples.**
  - `ch02__rotation_of_axes__phone__vector__tour-step3__pg3.png` (D01 step 5 quote cut)
  - `ch01__buckingham_pi_machine__phone__default__tour-step6__pg2.png`
  - `ch02__gauss_flux_box__phone__default__tour-step7__pg5.png`
  - `ch02__cauchy_traction_principal_axes__phone__default__explain__pg2.png`
- **Fix (orchestrator; library + tool).**
  - In `Pager.layout`, when a single item exceeds `avail`, split `.viz-work-*` / derive-quote items into sub-items (why / in words / live as separate pageable children), or cap the item with an internal "more ›".
  - Make `shot.py` page through every pager page (click `.viz-pager-btn[aria-label="Next page"]` until disabled) and audit each one. Without that, this class of defect passes the audit.

**P2 · continuum_averaging_volume, Equations page 2 (`pkin` line 2794, `mfp` line 2790), the page flagged as unchecked.**
- **Problem.** The live lines are clipped:
  - `p = 2.55×10^25 × 1.381×10^-23 × 288.15 = 1.013×10^5 Pa` is 428 px in a 332 px card at 360×640 (`equation-too-wide`, 96 px);
  - at 390×844 the "Pa" / "Pa)" ends are visibly cut, but the audit does not flag it (`located/phone-tall__00__continuum_averaging_volume__equations__pg2.png`).
- **Fix.** Split each `live` into `\begin{aligned}` two lines (product on line 1, `= 1.013\times10^5` Pa on line 2), or drop the operands and show only $p = n k_B T = 1.013\times10^5$ Pa.

**P3 · parcel_stability, Explain last page at 360×640.**
- **Problem.** The "Reading the current setting" paragraph (line 2824) runs past the card, 13 px clipped ("troposphere more than the surface…"). The profile view shrinks to 56 px and draws nothing (`view-too-small:profile`); ζ(t) shrinks to 43 px. Screenshot: `located/phone__06__parcel_stability__explain__pg12.png`.
- **Fix.** Split the interpretation into two `W.interpret` blocks (regime text · climate note), or move the climate note into its own `W.say`.

## Should fix
- **S1 · continuum line 2796.** The book's (1.21) is $pV = n k_B T$ with $n$ = number of molecules (p. 16). The note writes $pV = N k_B T$. Add "(the book writes $n$ for this count)" so the rename is visible. `knowledge/notation.md` already flags n ⚠️.
- **S2 · gauss lines 2651 and 2655.** The TeX $\iint_A \nabla\cdot\mathbf Q\,dA = \oint \mathbf n\cdot\mathbf Q\,ds$ is the planar form but is labelled "(2.30)". Write "Gauss in 2-D, the planar form of (2.30)", or show the book's $\iiint_V \partial Q/\partial x_i\,dV = \iint_A n_iQ\,dA$.
- **S3 · stokes line 2743.** `Stokes $\\oint_C = \\iint_A$ (2.34)` is not an equation. Write the integrands: `$\\oint_C \\mathbf u\\cdot\\mathbf t\\,ds = \\iint_A (\\nabla\\times\\mathbf u)\\cdot\\mathbf n\\,dA$`.
- **S4 · stokes line 2798 (Code title).** The TeX duplicates the `ref` "Eq. (2.25)". At 360×640 the title is orphaned at the bottom of page 5/7 with its code on the next page (`located/phone__16__stokes_circulation_loop__code__pg5.png`). Revert to "The curl and the orientation of C", or keep the title with its code in the pager.
- **S4b · rotation line 2574.** The status wraps to two lines at 360 px ("…by 0.866" / "Pa"). Shorten to `'❌ not a tensor: off $\\boldsymbol\\tau\' = \\mathbf C^{\\rm T}\\boldsymbol\\tau\\mathbf C$ (2.12) by ' + … + ' Pa'`, or drop "misses".
- **S5 · pre-existing plain-text braces.**
  - parcel (D18 *why* "e^{λt}", "e^{±λt}"; D20 check "2^{2/7}")
  - strain (readouts "× e^{λ₁t} = e^{0.00}", "det e^{Gt} = e^{(∇·u)t}"; Equations note of the flow-map card)

  Use e^(λt) / 2^(2/7) (ch01 lesson), or put them in TeX.
- **S6 · canvas labels in rotation tensor mode (pre-existing).** On phones the left view clips "cos 30° = 0.8…" and "…eclared "the same in every fram…". The matrix view prints raw "C_i₁ C_j₁ τ_ij" (`ch02__rotation_of_axes__phone__fixed__explore__pg1.png`). Use `Viz.text` with `maxWidth` and Unicode subscripts.

## eq_refs.py: remaining hits judged
- **ch01, 19 hits.** All exempt or false positives:
  - `ref:` labels next to shown TeX (bpm 3204, 3247; heat 2938; parcel 2937, 2973);
  - derivation `short:` (heat 2795, 2833);
  - meta (continuum 14);
  - selftest names (bpm 3319; parcel 3094, 3108, 3112; viscosity 2733, 2737);
  - code/py (continuum 2821);
  - numbers in parentheses (parcel 2916; viscosity 2722, 2726);
  - parcel 2921: "the isentropic relation above, (1.26)" sits right under its start TeX, so it is shown.
- **ch02, 47 hits.** All exempt or false positives:
  - meta (cauchy 12, gauss 12, stokes 11);
  - CSS comment (cauchy 458);
  - `ref:` labels (cauchy 2626, 2724, 2756; gauss 2730, 2774; stokes 2798; strain 2687, 2726);
  - `short:` (stokes 2697, 2751);
  - selftest names (the rest);
  - coordinates (rotation 2736, 2875; stokes 2676, 2845);
  - stokes 2617: a preset tooltip "Γ/A → (∇×u)₃ at the centre (2.35)", a compact written form. Acceptable.

## Screenshots examined
- `located/`:
  - `phone__00`, `phone-tall__00` (continuum Equations p.2)
  - `phone__02` (viscosity species §0)
  - `phone__04` (heat stir paddle box)
  - `phone__06` (parcel Explain last page)
  - `phone__10` (gauss Explore Right now)
  - `phone__11` (gauss §6)
  - `phone__16` (stokes Code title)
  - `phone__17` (strain §6)
  - `desktop__15` (stokes Explain + subtitle)
- `eq_retrofit/`:
  - `ch02__rotation_of_axes__phone__fixed__explore__pg1` (status on phone)
  - `ch02__rotation_of_axes__phone__vector__tour-step3__pg3`
  - `ch01__buckingham_pi_machine__phone__default__derive-d1p14__pg8`, `__pg10` (D28 result and check)
  - `ch01__buckingham_pi_machine__phone__default__explain__pg10`, `__pg11` (§5–§6)
  - `ch01__buckingham_pi_machine__phone__default__tour-step6__pg2`
  - `ch02__gauss_flux_box__phone__default__tour-step7__pg5`
  - `ch02__cauchy_traction_principal_axes__phone__default__explain__pg2`

## Verdict: FAIL
- **R1 and R2 are retrofit regressions.** Both are small text edits.
- **L1, P2 and P3 are pre-existing** fit violations on sub-pages that `shot.py` does not audit. L1 needs a library and `shot.py` change by the orchestrator.
- The other eight explainers' retrofit strings are correct, render in KaTeX and fit.
