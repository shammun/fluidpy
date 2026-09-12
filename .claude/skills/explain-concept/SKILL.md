---
name: explain-concept
description: Explain one concept, equation or result from a chapter to the user right now, in Shammunul's teaching style, grounded in the book's rendered pages and the project's tested code - plain words, the idea, the maths step by step, a tiny worked example with easy numbers, runnable Python using fluidpy (executed, real output), and links to the matching notebook section and explainer; optionally saves it as a study note. Usage - /explain-concept 7 "why group velocity is half the phase velocity in deep water"
---

# /explain-concept N "<concept>"
For learning, not for the pipeline. Work directly (no subagents unless the chapter text search is large).

1. Locate: grep `analysis/chNN.md`, `analysis/chNN_curation.md`, `knowledge/chNN.md`, `knowledge/concept_map.md` and
   `chapters/chNN.txt` for the concept; render the relevant page(s) with
   `.venv/Scripts/python.exe tools/render_pages.py chNN --eq N.M` (or `--find "<phrase>"`) and read them.
2. Explain in the `teaching-style` chain: problem in plain words → the idea (ASCII sketch if useful) → maths step by step
   with equation numbers → tiny example with easy numbers → Python (call the tested `fluidpy` function if it exists,
   plus a transparent from-scratch version; run it with `.venv/Scripts/python.exe` from a scratch file and show the real
   output) → what the output shows → a "⚠️ Common confusion" and "✅ What should click" line.
3. Point to where it lives: notebook section (`notebooks/chNN_<slug>.ipynb`), explainer (`viz/chNN/<slug>.html` +
   the tab/step), function and test names.
4. If the user asks to keep it: write `notes/chNN_<short_slug>.md` (our words, no book prose; equations with numbers).
   Never commit book text or page images.
