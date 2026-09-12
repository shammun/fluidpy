---
name: teaching-style
description: Shammunul's educational style for fluidpy notebooks and explainers - every new idea taught in full and nothing used unexplained (CORE blocks, recaps, primers), plain words before symbols, the problem → idea → maths → code → output chain, tiny worked examples with easy numbers, every code line commented for a novice, "What does the code above do?" and "What you see / How to read it / What would change if" blocks, tables and ASCII sketches, emoji section markers, and the rules for walkthrough text. Load when writing or reviewing any notebook cell, explainer text, curation or design document.
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
   book's equation number: *(Eq. 7.27)*.
4. **A tiny example with easy numbers** — traced by hand before any code: H = 1 m, k = 1 m⁻¹, g ≈ 10 m s⁻². Show the
   intermediate numbers. The reader should be able to check it on paper.
5. **Code** — short cells that call the tested `fluidpy` function; every line commented for a novice.
6. **What the output shows** — never leave a figure or number uninterpreted.

## 1b. Nothing unexplained (the self-contained rule)
Every new idea of the chapter is a CORE block with code and a visualization. Beyond that, **every concept, term, symbol,
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

## 2. Recurring blocks (use these exact headings — readers learn to scan for them)
- `**What is this section about?**` — 2–4 sentences at the top of a big section.
- `**What does the code above do?**` — after any non-trivial code cell: a numbered list of what each part does, in the
  order it runs, naming the variables.
- `**What does this show?**` — after a printed result: what the number/shape means, with a quick sanity check
  ("512 = one mean per feature").
- After every figure: `**What you see.**` · `**How to read it.**` · `**What would change if…**` (one concrete change
  and its predicted effect — ideally something the reader can try).
- `### ✏️ Tiny example: …` · `### 🎮 Interactive: <what it makes clear>` followed by `**What to try:**` bullets, each an
  action + what to watch for ("Drag H down to 1 m and watch the orbits flatten into ellipses").
- `> ⚠️ Common confusion:` callouts for misconceptions (phase vs group velocity, ν vs μ, streamline vs pathline).
- `### 🧩 <CORE idea>` opens every new idea; `> 📎 **Primer — term.**` and `> 🔁 **Recap — idea** (Ch. N §N.M).` blocks explain prerequisites where they are needed.
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

## 6. Lessons (the knowledge-keeper appends; one line each, dated by chapter)
- (none yet)
