# Chapter 7 — lesson review (round 2, after fix round 1)          2026-09-24

Reviewer: `lesson-reviewer` (first-time reader: strong maths/statistics/Python, no fluid mechanics).
Input: `outputs/ch07/executed.ipynb` (611 cells, executed 11:50, newer than every `viz/ch07/*.html` (latest 11:39) and
than `notebooks/build_ch07.py` / `ch07_gravity_waves.ipynb` (11:47); cell sources identical to the committed notebook,
0 error outputs, **0 stderr outputs**). Scratch: `outputs/ch07/lesson_r2_{dump,checks,latex,img}.py`, images in
`outputs/ch07/r2img/`.

coverage_check: **OK** (0 errors, 0 warnings) · sections covered **8/8** (§7.1–§7.8 + S01 pointer) · CORE blocks
**16/16** with code + visual + tiny example · derivations **37/37**, 321 steps, 4 ★★★ with sympy cells · primers
P165–P184 (20; P179–P184 new) all placed before first use · **9/9** explainers embedded with real content (246–281 kB each,
no `fluidpy-viz-missing` placeholder), each with "why interactive" + "What to try".

Re-verified against the rendered pages (p317, p318, p321, p322, p324): the equations newly written in *We start from*
of D30 ((7.101), (7.106)–(7.109)), D32 ((4.9), (4.10), (7.120)–(7.124)), D35 ((7.140), (4.10)) and D37 ((7.128), (7.131),
(7.132)) match the book in every sign and factor. The new energy note (cell 595) was re-derived:
⟨u²+w²⟩ = ½ŵ²(m²/k²+1) → E_k = ¼ρ₀(m²/k²+1)ŵ² (7.154); E_p = ½g²⟨ρ′²⟩/(ρ₀N²) = N²ρ₀ŵ²/(4ω²) (7.155); equal because
N²/ω² = K²/k² — correct. The KdV crossover 4π²/9 (cell 413) was re-derived: (3/2)(a/H)k ÷ (1/6)H²k² = 9a/(k²H³) — correct.

## Round 1 resolution

| # | Round-1 item | Now (cell) | Status |
|---|---|---|---|
| M1 | `internal_wave_beams` placeholder | 607: embedded, 270 kB; all 9 explainers older than the run | **resolved** |
| M2 | aω ≈ 0.14 m/s | 56–57: code prints aω = 0.154, ka·aω = 7.7e-3; prose "≈ 0.15 / ≈ 8×10⁻³" | **resolved** |
| M3 | "0.73" on the 30 m shelf | 220 prints F/(Ec) = 0.68; 221 says 0.68 | **resolved** |
| M4 | island rays never reach the lee | 373–374: rays launched to ±4.6 km, 2 rays counted landing on the east half (thick orange, printed); caption rewritten; explainer bullet (376) still achievable | **resolved** (see S5) |
| M5 | E (7.157) used in D37 before explained | 595 (energy note with E_k, E_p, E_k = E_p, E, and the ½Re(AB*) rule) now precedes D37 (596) | **resolved** |
| M6 | bare numbers in D30/D32/D35/D37 | 493, 537, 571, 596: every cited equation written in *We start from*; D32 steps 9–10 written out | **resolved** |
| M7 | 17.8 cm/s (exercise answer) | 238, 254, 337: "0.1776 m/s at λ ≈ 4.35 cm, computed with `min_group_velocity`"; printed by 336 | **resolved** |
| M8 | unglossed idioms / terms | P179 np.stack/np.c_ (38), P180 np.sign/np.nonzero (296), P181 rfft/rfftfreq (307), P182 np.interp (330), P183 catch_warnings (471), P184 np.errstate (559), WKB defined (348), Hamilton's equations glossed (365) — each before first use, each with a runnable demo whose output is right | **resolved** (see S3 for a notation nit) |
| S1 | findfont / subplots_adjust warnings | no stderr output anywhere | resolved |
| S2 | duplicated auto-inserted equations | mostly gone (1 regex hit left), but see S1 below: 319, 466, 487, 497, 588 | **partly** |
| S3 | D14 step 5 "first form" | 212: shows $p'=-\rho\phi_t=\rho\frac{a\omega^2}{k}\frac{\cosh}{\sinh}\cos$ | resolved |
| S4 | broken LaTeX in `\text{}` results | D08, D28, D31 results now in maths | resolved (raw `e^{…}` in *why* lines remains, S2 below) |
| S5 | forward uses in §7.1/C02 | 32, 43 "take it on trust for now"; D03/D04 size from η alone | resolved |
| S6 | N20–N23 step pointer | 105: "(step 6)" | resolved |
| S7 | D07 constant factor | 131: "c² up to the positive constant g/2π" | resolved |
| S8 | "extracted text" jargon | gone | resolved |
| S9 | depth rule three ways | 126, 153, 610 and the plotly title: "about three times the depth (H > 0.32λ)" | resolved |
| S10 | a = 1 m tiny example | 179: 10 cm radii | resolved |
| S11 | plotly figures lack "What would change if" | 152, 154, 250, 280, 335, 375, 514, 563 all have it | resolved |
| S12 | figure readability | 358 x–t: every 2nd crest + one traced crest; 510: baroclinic dashed with markers + label; 479: γ labels above the arrows | resolved |
| S13 | uncomputed book numbers | Ursell "≈ 16" and "10–20" gone, crossover computed (414); a_max computed (434); "0.71N photograph" gone | resolved |
| S14 | 0.1411 attribution | 433: Schwartz & Fenton "about one-seventh"; digits from Dyachenko, Lushnikov & Korotkevich (2016) | resolved |
| S15 | Stokes trap paraphrased an exercise | 423: "a tempting shortcut", no exercise description | resolved |
| S16 | simple-wave speed unsourced | 386: "Riemann invariant … stated, not derived" | resolved |
| S17 | cell 0 wording | "every equation is shown in full together with its number" | resolved |
| S18 | small "What does the code above do?" gaps | 189, 219, 277, 405, 506 still have no line after them | **open** (S6 below) |

## Must fix

None.

## Should fix

1. **Leftover auto-inserted duplicate equations** (the same equation twice, the second copy inside italics):
   - cell 319: "*(7.70), $c_g=c/2\ …$*" after the text form;
   - cell 466: "… at z = 0 *(7.93), $\frac{\partial\phi_1}{\partial z}=…$*";
   - cell 487: (7.98), (7.100), (7.101) each twice;
   - cell 497 step: "In (7.109), $b=…$ gk/ω² = 1";
   - cell 588: the bold bullet title "**density jump as a delta function (7.151), $E_p=…$, (7.152), $N^2=…$**".

   Keep one copy, outside the italics.
2. **Plain-text TeX in *why* / *in words* lines:** `e^{kz}`, `Ae^{−2kH}`, `e^{iθ}`, `(4gσ/ρ)^{1/4}` appear as literal braces in
   cells 102, 103, 134, 135, 238, 312, 439, 440, 468, 489, 493, 497, 508, 571, 596 (reader sees "e^{kz₀}"). Wrap them in
   `$…$` (e.g. $e^{kz_0}$) or write them in Unicode (e^(kz₀)). Also cell 171 step 1: `\text{and the same for z\_p …}`
   shows a literal "z_p". Write $z_p$ in maths.
3. **Cell 365 (D24 step 1, Hamilton gloss): letter clash.** In the same sentence, H is the depth H(x, y) and the Hamiltonian
   ("dx/dt = ∂H/∂p … ω as the energy H"). Use $\mathcal H$ (or "the Hamiltonian") for the mechanics analogy.
4. **Exercise labels / book quote left in prose:**
   - cell 317 "(Exercise 7.9, `group_velocity(g=0)`)";
   - cell 493 "(The book: "after some algebraic manipulations", Exercise 7.19; written out here.)".

   Drop the labels and the quoted phrase. Say "the book skips this algebra; we write it out". The S01 pointer (609) is the
   one agreed home for exercise numbers.
5. **Cell 374 island caption:**
   - "Rays that start beyond the shelf edge (dashed circle) pass by almost straight" is not shown: no ray is launched
     beyond |y| = 4.6 km, inside the 5 km shelf edge. Launch one or two rays at |y| ≈ 5.5–6 km, or drop the sentence.
   - The outermost purple rays (±4.6 km) swing furthest behind the island and cross each other, but are not orange because
     they end offshore. Add one clause saying so, so the reader does not wonder why the most-wrapped rays are not
     highlighted.
6. **Cells 189, 219, 277, 405, 506:** checks that print numbers are followed directly by a figure, with no one-line reading.
   Add a short "**What does this show?**", e.g. 189 "ψ-derived u, w equal (7.27) to 1e-10 m/s", 405 "the negative root
   is unphysical, 3.772 = (7.81)". Same for 603 ("our finite-difference gradient = sign-safe c_g, and c·c_g = 0").
7. **Heading mismatch on plotly cells** 152, 154, 250, 280, 335, 375, 514, 563: the text under "**What does the code above
   do?**" is really *What you see / How to read it*. Rename the heading (the content is good).
8. **Cell 595 wording:** "the amplitudes are the plane-wave ones D37 finds in steps 3–4". Add "(derived line by line just
   below)" so the reader knows the forward pointer is deliberate. Notes 595 and 598 also share the ID `N159–N166`. Split them
   (e.g. N159–N162 energies, N163–N166 polarization/flux).

## Unexplained-first-use list (term · first cell · explained at cell / missing)

Only the rows that were open in round 1, plus the new material. All the other round-1 rows stay "ok".

| Term / symbol / tool | First used | Explained at | Status |
|---|---|---|---|
| `np.stack(axis=-1)`, `np.c_` | 40, 42 | P179 (38–39) | ok |
| `np.sign`, `np.nonzero` sign-change | 298 | P180 (296–297) | ok |
| `np.fft.rfft`, `rfftfreq` | 309 | P181 (307–308) | ok |
| `np.interp` | 332 | P182 (330–331) | ok |
| `warnings.catch_warnings(record=True)` | 473 | P183 (471–472) | ok |
| `np.errstate(invalid="ignore")` | 561 | P184 (559–560) | ok |
| WKB | 348 | 348 (defined in place) | ok |
| Hamilton's equations (rays) | 365 | 365 (gloss; letter clash, S3) | ok |
| internal-wave energy E (7.157), E_k, E_p | 596 (D37 step 11) | 595 | ok (was late) |
| (7.128), (7.131), (7.132) in D37 | 596 | 596 *We start from* | ok |
| deep-water ω = √(gk) in §7.1 | 32, 43 | "take it on trust", derived in C03–C04 | ok (declared forward use) |
| 1-D Dirac delta | 589 | 588 (reminder of Ch. 6 P149) | ok |
| Ursell number, KdV crossover | 413 | 413–415 | ok |

No missing rows.

## Derivation audit (D id · steps · every step follows? · why/in-words ok · check ok · verdict)

Round-1 verdicts stand for D01–D36 (all pass). Changes since round 1, re-read line by line:

| D | Steps | Every step follows? | Why / in words | Check | Verdict |
|---|---|---|---|---|---|
| D03 linearise + transfer | 8 | yes; sizing now from η alone | ok | residual slopes ✓ | pass |
| D04 dynamic BC | 8 | yes; φ ~ aω/k argued without C03 | ok | ✓ | pass |
| D07 c(λ) increasing | 5 | yes; constant g/2π stated | ok | ✓ | pass |
| D10 linear orbits | 9 | yes | step 1 `\text{z\_p}` (S2) | ✓ | pass |
| D14 energy flux | 10 | yes; correct form of p′ in step 5 | ok | ✓ | pass |
| D24 Snell | 9 | yes | Hamilton gloss present; H clash (S3) | RK4 = Snell 1e-6 ✓ | pass |
| D30 ★★★ two-layer dispersion | 12 | yes | *We start from* has all equations (checked with p317/p318) | sympy ✓ | pass |
| D32 linear equations | 10 | yes | *We start from* and steps 9–10 written (checked with p321/p322) | ✓ | pass |
| D35 K·u = 0 | 5 | yes | (7.140), (4.10) written | ✓ | pass |
| D37 polarization, F = c_gE | 12 | yes (steps 2–12 re-done: p̂ = −ωmρ₀ŵ/k², ρ̂ = iN²ρ₀ŵ/(ωg), û = −mŵ/k, F and c_gE agree via N/K = ω/k) | E now explained before step 11; (7.128)/(7.131)/(7.132) shown | units ✓, `internal_wave_energy` F = c_gE, E_k = E_p to 1e-12 ✓ | **pass** (was fail) |

## What works well (keep; candidates for knowledge/ and the teaching-style skill)

- **Idiom primers with a 2-line demo whose printed output is the explanation** (P179–P184): e.g. `np.c_[[0,0],[3,4]]` prints
  the x-row and the y-row that `ax.plot(*…)` needs; `np.errstate` shows `[nan 0.5]`. A good pattern for `knowledge/primers.md`.
- **"Take it on trust for now (derived in C03–C04)"**: an honest, explicit forward use that keeps §7.1 concrete without
  hiding the dependency. Worth adding to teaching-style §1b.
- **Count what the caption claims**: the island figure now prints the number of rays that reach the lee side, so the caption
  is checked by the cell itself.
- **Scaling derived from the forcing, not from a later result** (D03/D04: "the surface rises at aω, so velocities are aω;
  a potential varying over 1/k is aω/k"): a reusable derivation move.
- **Crossover computed rather than quoted** (Ursell 4π²/9 from the two KdV terms) is the rule-9 fix pattern for
  "book threshold" numbers.
- The round-1 strengths remain: the C02 residual scan, ghosts of printed slips, the Stokes γ table, the c ⟂ c_g triangle,
  and cross-chapter parity checks.

## Verdict: PASS

All 8 round-1 Must-fixes are resolved and verified cell by cell (M1 re-execution, three numeric/claim fixes, D37 ordering,
equations written beside every number, the exercise answer replaced by our computed 0.1776 m/s, eight idioms/terms
glossed before first use). 16 of 18 round-1 Should-fixes are resolved, 1 partly (S2), 1 open (S18). The fixes created no
new forward references, use-before-explanation or book-quoted numbers. What remains are the 8 Should-fixes above:
formatting (duplicate/raw TeX), one notation clash, two leftover exercise labels, a caption sentence and small reading
notes.

## Round 2 (should-fix clean-up, orchestrator note, 2026-09-24)
All 8 remaining Should-fix items fixed by the notebook-builder (duplicated equations, raw TeX pass `tidy_raw_tex` + build-time guard, ℋ vs H, exercise labels and the verbatim book phrase removed, island caption matches the traced rays, 6 new reading cells, plotly "What you see" headings, N159–N166 IDs split). Re-executed: 617 cells, 0 errors, 189 s; coverage_check OK (0/0). Known cosmetic leftover: repeated "at z = 0" in two D28/D29 step lines.

**Verdict stays PASS.**
