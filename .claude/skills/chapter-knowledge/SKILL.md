---
name: chapter-knowledge
description: How knowledge is captured and carried forward between chapters in fluidpy - knowledge/chNN.md, CUMULATIVE.md, concept_map.md, notation.md, viz_patterns.md, progress.json, and promotion of lessons into skills and of shared helpers into fluidpy/core and assets/viz_lib.js. Load at the start of every chapter (to read) and at the end (to write).
---

# chapter-knowledge — make chapter N+1 smarter than chapter N

## Read at chapter start (every agent, in this order)
1. `knowledge/CUMULATIVE.md` — physics pipeline so far, primitives, validation summary, pitfalls, explainer inventory.
2. `knowledge/notation.md` — symbol → meaning → SI unit → sign convention → chapters (our words, not the book's table); `knowledge/primers.md` — concepts and tools already explained (reuse, do not repeat).
3. `knowledge/viz_patterns.md` — explainer stage ideas, controls and walkthrough moves that worked; ideas that failed.
4. `knowledge/concept_map.md` — every implemented concept with its function and validation label.
5. `knowledge/ch<N-1>.md` → "Feeds forward".

## Write at chapter end: `knowledge/chNN.md`
```
# Chapter N — <title>   (completed <date>, commit <hash>)
## 1. What clicked — the teaching spine (in our words, each with its equation numbers)
## 2. Implemented items (concept → fluidpy target → book §/Eq. → validation label → tier)
## 3. New reusable primitives (fluidpy/core and assets/viz_lib.js; signature + one line + who will reuse)
## 4. Mathematics → Python lessons (only NEW pitfalls)
## 5. Conventions used (reference scales, signs, coordinates, dimensional vs non-dimensional) and changes vs earlier chapters
## 6. Results vs benchmarks (and vs the book, as relative errors only)
## 7. Explainers and figures (slug · CORE idea · what worked · what the reviewer flagged) — and animations/plotly figures
## 8. Feeds forward — what later chapters need from this one
## 9. Open questions / revisit (every qualitative or unverified item)
```

## Update `knowledge/CUMULATIVE.md` (rewrite; ≤ ~400 lines)
*Physics pipeline so far* (text diagram) · *Available primitives* (`fluidpy/core` table with labels) · *Explainer
inventory* (chapter · slug · idea · reusable stage pattern) · *Global pitfalls* (deduplicated) · *Benchmark inventory* ·
*Validation summary per chapter* (counts per label + verdict).

## Update `knowledge/concept_map.md`
`| Concept / theorem | Book § / Eq. | fluidpy target | Validation label | Tier | Taught in (notebook § / explainer) | Verified in |`
Never delete a verified row; superseded rows say "superseded chNN".

## Update `knowledge/notation.md`
`| Symbol | Meaning | SI unit | Convention / sign | Chapters | Code name |` — add every symbol the chapter's code uses; flag
symbols whose meaning changes between chapters (e.g. ω angular frequency in Ch. 7, vorticity magnitude in Ch. 5).

## Update `knowledge/primers.md`
One row per primer the chapter's notebook wrote (`metadata.fluidpy.primers`, design Part E):
`| Term | Explained in (chNN · § · CORE id) | One-line gist |`. Later chapters point here instead of re-priming.

## Update `knowledge/viz_patterns.md`
Sections: *Patterns that worked* (stage, control, walkthrough move — with the explainer that proved it) · *Failures and
fixes* (fit problems, confusing controls, parity surprises) · *Promotion candidates* (helpers duplicated across
explainers, with file references) · *Ideas for later chapters*.

## Promotion (only when the orchestrator's brief allows it — nothing else may be reading `viz/` at the time)
- JS helper duplicated in ≥ 2 explainers → `assets/viz_lib.js` (documented in its header, backwards compatible) →
  `tools/viz_inline.py --all` → `tools/shot.py templates/viz_example.html --quick` and `--chapter <each affected> --quick`
  must PASS, else revert.
- Python helper used by ≥ 2 chapters → `fluidpy/core/<topic>.py` with re-exports from the old location → `pytest -q`.
- A lesson that applies to any chapter → the relevant skill's "Lessons" section (one line, tagged with the chapter).

## progress.json
The knowledge-keeper sets `knowledge: "done"` + notes for its chapter; the orchestrator owns every other phase key,
`current_chapter`, and all commits.
