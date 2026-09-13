---
name: teaching-style
description: Shammunul's educational style for fluidpy notebooks and explainers - exhaustive coverage with tiered depth and nothing used unexplained (A/CORE blocks in full, B/C notes, recaps, primers), derivations written out one small move per step with what-we-did / why-allowed / in-words and a check + interpretation, plain words before symbols, the problem → idea → maths → code → output chain, tiny worked examples with easy numbers, every code line commented for a novice, "What does the code above do?" and "What you see / How to read it / What would change if" blocks, tables and ASCII sketches, emoji section markers, and the rules for walkthrough text. Load when writing or reviewing any notebook cell, explainer text, curation or design document.
---

# teaching-style — make it click for a first-time reader

Distilled from Shammunul's own annotated notebooks (`fast.ai/shammunul-fastai-notes`, e.g. `18_fid_explained.ipynb`,
`14_augment_explained.ipynb`). The reader is smart, busy, strong in Python and statistics, and new to this physics.

## 1. The logic chain (every CORE idea, in this order)
1. **The problem in plain words** — why anyone cares. A concrete scene first (a wave arriving at a beach, air over a
   wing, a pollutant in a river, the atmosphere over a mountain). No symbols yet.
2. **The idea** — the one mental picture that makes the result feel inevitable. An ASCII sketch or a two-column table
   often does it:
   ```
   deep water:  particles move in circles ──► shrink with depth like e^{kz}
   shallow:     particles move in flat ellipses ──► the bottom squashes them
   ```
3. **The maths, step by step** — one algebraic move per line, each with a short reason ("substitute (7.22) into the
   Laplace equation", "the sinh cancels"). Define every symbol with its unit the first time it appears. Cite the
   book's equation number: *(Eq. 7.27)*. Anything longer than a line or two is a full **derivation** (§1c).
4. **A tiny example with easy numbers** — traced by hand before any code: H = 1 m, k = 1 m⁻¹, g ≈ 10 m s⁻². Show the
   intermediate numbers. The reader should be able to check it on paper.
5. **Code** — short cells that call the tested `fluidpy` function; every line commented for a novice.
6. **What the output shows** — never leave a figure or number uninterpreted.

## 1b. Nothing unexplained (the self-contained rule)
Every A-depth (load-bearing) idea of the chapter is a CORE block with code and a visualization; B ideas are stated
and explained inside the nearest CORE block, C ideas named with a pointer (`book.yaml → policy`). Beyond that, **every concept, term, symbol,
mathematical tool and Python function/idiom is explained no later than the cell where it is first used** — by its own
CORE block, by a 🔁 recap of an earlier chapter, or by a 📎 primer right there. Typical primer candidates:
- maths: partial derivatives, gradient/divergence/curl, the chain rule in several variables, Taylor expansion, complex
  exponentials, Fourier modes, eigenvalues, integrals over surfaces and volumes, order-of-magnitude estimates, logarithmic axes;
- physics vocabulary: control volume, steady vs unsteady, inviscid, incompressible, boundary condition, dimensionless group;
- Python: numpy broadcasting, `np.meshgrid`, vectorised `np.where`, `scipy.integrate.solve_ivp`, `scipy.optimize.brentq`,
  `matplotlib` contour/quiver/streamplot, f-strings, `lambda`, `assert np.allclose`.
Primer format (`nb.primer(term, text, code)`): one or two plain sentences, what it means *here*, and for maths or Python a
2–4-line runnable demo with easy numbers. If an earlier chapter already primed it (`knowledge/primers.md`), write one
reminder sentence and name the chapter instead of repeating the primer. When unsure whether the reader knows something,
prime it — a two-line primer costs less than a lost reader.

## 1c. Derivations — every move shown, every move explained (hard maths made easy)
Whenever a result is obtained by manipulating equations (more than one line of algebra or calculus from things the
reader already knows), it is a **derivation** (`D01…` in the curation) and is written out in full — in the notebook
(`nb.derivation`, so also in Colab and on the web page) and, when an explainer covers that idea, in the explainer's
**Derivation** tab with the same steps.
- **Say the goal first, in plain words**, and why we want it ("show that friction can only remove energy"). Then **the
  plan** in 2–4 bullets ("multiply by the velocity; spot two exact derivatives; read off the sign"), then **the tools**
  it uses (each already taught, recapped or primed — a tool first used inside a derivation gets its primer *before* it).
- **One small move per step.** Never "it can be shown", "after some algebra", "similarly", "clearly". If a step needs
  two moves, it is two steps. Hard derivations (★★★) may have 10–15 steps; that is fine.
- **Every step has four parts**:
  1. *What we did* — the move in a few words ("divide both sides by $m$", "integrate from the wall to $y$").
  2. *The new line* — the equation after the move; keep it short enough for a phone (split long lines).
  3. *Why we can do this* — why the move is allowed (name the rule: chain rule, product rule, divergence theorem,
     $e^{rt} \neq 0$, "the pressure does not depend on $x$ here") **and** why we make it (what it is heading toward).
  4. *In words* — what the new line says physically, in one sentence.
- **Say every assumption the moment it is used** ("steady, so $\partial/\partial t = 0$"; "inviscid, so the viscous
  term drops") and mark approximations with ≈ and the reason ("$kH \ll 1$, so $\tanh kH \approx kH$").
- **Colour-code terms consistently**: the same term has the same colour in the derivation, the explainer's bars and
  curves, and the figure legend (e.g. pressure orange, viscous rose, inertia teal).
- **After the result: check it, then interpret it.** Units of both sides; one or two limits or special cases that the
  reader knows ("$\zeta = 0$ gives energy conservation"); a sympy check cell for ★★★ derivations (`check_src`, every
  line commented); a number with easy values; then **what it means** and when it fails.
- **Link the derivation to the picture**: which curve, bar or readout shows the result, and one "try this" in the
  explainer (a derivation step may `set` the explainer to the case being derived).
- Avoid the classic traps: sign flips when moving terms, dropping a factor 2 from a square, differentiating a product
  without the product rule, swapping the order of integration and differentiation without saying why it is allowed.

## 2. Recurring blocks (use these exact headings — readers learn to scan for them)
- `**What is this section about?**` — 2–4 sentences at the top of a big section.
- `**What does the code above do?**` — after any non-trivial code cell: a numbered list of what each part does, in the
  order it runs, naming the variables.
- `**What does this show?**` — after a printed result: what the number/shape means, with a quick sanity check
  ("512 = one mean per feature").
- After every figure: `**What you see.**` · `**How to read it.**` · `**What would change if…**` (one concrete change
  and its predicted effect — ideally something the reader can try).
- `#### 🧮 Derivation — <result> (Eq. N.M)` (written by `nb.derivation`): What we want to show · The plan · Tools we
  use · We start from · Step k of n — <move> (line · *Why we can do this* · *In words*) · Result · What it means · Check it.
- `### ✏️ Tiny example: …` · `### 🎮 Interactive: <what it makes clear>` followed by `**What to try:**` bullets, each an
  action + what to watch for ("Drag H down to 1 m and watch the orbits flatten into ellipses").
- `> ⚠️ Common confusion:` callouts for misconceptions (phase vs group velocity, ν vs μ, streamline vs pathline).
- `### 🧩 <CORE idea>` opens every A-depth idea (B/C items sit inside it as 📝 notes); `> 📎 **Primer — term.**` and `> 🔁 **Recap — idea** (Ch. N §N.M).` blocks explain prerequisites where they are needed.
- `## ✅ What should have clicked` at the end: one bullet per CORE idea, each a sentence the reader could now explain to a friend.

## 3. Code comments (novice-grade, never noise)
```python
k = 2 * np.pi / wavelength          # wavenumber k [1/m]: how many radians of phase per metre
omega = np.sqrt(g * k * np.tanh(k * H))   # Eq. (7.36): the dispersion relation — frequency set by k and depth H
c = omega / k                       # phase speed [m/s]: how fast a crest moves
```
Comment the *meaning* and the unit, cite the equation; don't restate Python syntax ("# assign x").
From-scratch versions: a transparent loop or formula next to the library call, then
`assert np.allclose(mine, fluidpy_result)  # same numbers → the library does exactly this`.

## 4. Tone and format
- Second person, present tense, friendly and exact. Short paragraphs (≤ 4 sentences). Bold the one phrase to remember.
- Tables for comparisons (regimes, limits, methods: | Method | Problem |). ASCII diagrams for pipelines and geometry.
- Emoji as section markers only (⚙️ setup, ✏️ example, 🎮 interactive, ⚠️ confusion, ✅ summary) — never decoration in prose.
- Units always (SI); numbers with sensible significant figures; say "≈" when rounding.
- Our own words only: never paraphrase the book sentence by sentence, never copy its figures, tables or exercises.
  Equations are fine, with their numbers.
- Link backwards and forwards: "this is the same Laplacian you met in Ch. 4", "Ch. 13 will add rotation to exactly this".
  For Shammunul, name climate/ocean/atmosphere connections when they are real (geostrophy, Ekman, Rossby, stratification).

## 5. Walkthrough text inside explainers
- 4–8 steps; each step ≤ 45 words, one idea, and it *does something* to the picture (`set:` values, an animation, a
  highlight). Step 1 asks the question in plain words; the middle steps build the idea and show the equation with live
  numbers; the last step hands the controls over with a prediction challenge ("Predict first, then drag…").
- Use "you" and imperatives: "Drag the depth slider…", "Watch the red dot…".
- Check-yourself questions must be answerable by experimenting with the explainer, and the answer explains *why*.
- **Explain tab** ("Explanation & interpretation", modelled on `forced_damped_vibrations.html`): numbered sections that
  compute every number on the screen with the reader's settings (formula → numbers substituted → result, with a short
  why), coloured words that match the curves, a hint for any view hidden on phones, the values at the current time, and
  a regime-dependent **Reading the current setting** paragraph.
- **Derivation tab**: the same steps as the notebook's derivation (§1c), each step ≤ 35 words of *why* and one sentence
  *in words*; add `live` (the line with the reader's numbers) on the step where numbers help most.

## 6. Lessons (the knowledge-keeper appends; one line each, dated by chapter)
- (ch01) **Conventions, both on screen.** Lapse rate: compute with Kundu's Γ ≡ dT/dz (Γ_a ≈ −9.8 K/km, stable when
  Γ > Γ_a) and always show the meteorology Γ ≡ −dT/dz (+9.8 K/km, stable when Γ < Γ_a) alongside — a two-row table, a
  worked conversion (negate the number, flip the inequality, P48) and both forms in every slider legend and badge.
- (ch01) **"Which p_o?" pattern** for any book-vs-field-vs-code convention difference: state each convention, give the
  size of the difference in numbers (θ at the ground 287.07 vs 288.15 K), say which one the code uses and why — as a
  `> ⚠️ Common confusion:` callout next to the first use.
- (ch01) ★★★ sympy check cells re-run the derivation's own construction (build what step k builds, then test it), not
  only the final result — the reader sees the moves verified, not just the answer.
- (ch01) Make an approximation measurable: a small table of "linearised vs full" with the gap growing with amplitude
  tells the reader when the assumption fails.
- (ch01) Every number in prose is computed by a cell (or checked against one); after renumbering derivation steps,
  check every "Dxx step N" pointer by script.
- (ch01) A tool first used inside a derivation (∂, ∇, complex roots and cosh, the double root, Maxwell relations) gets
  its 📎 primer *before* that derivation; forward pointers to a later block say "(worked with numbers in C70 below)".
- (ch01) Figure annotations that the reading notes never mention (stray symbols, point letters) are removed or explained;
  a comment's claimed range must match the plotted data.
- (ch01) Reusable derivation moves (name them in *why*; recap by pointing to the chapter's D row):
  divide Newton's second law by the body's **own** mass/density, then Taylor-expand the surroundings to first order
  (D18) · eliminate a path function (dq) between two relations to leave state-function relations (D10) · set ds = 0 in
  both Gibbs forms and divide them so the unknown cancels (D14) · introduce g = h − Ts to get a Maxwell relation that
  trades an unmeasurable derivative for a measurable one (D19) · log-differentiate a product of powers into a sum of
  relative rates (D36) · try e^{λt} to turn a linear ODE into algebra, and treat the double root separately (D18) ·
  integrate pressure over a box's faces: sides cancel, top minus bottom is the net force (D05, D37) · "dimensionless ⇔
  exponent vector in the null space", then rank–nullity counts the groups (D28) · Poisson variance = mean gives relative
  scatter 1/√N (D35).
