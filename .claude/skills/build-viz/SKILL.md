---
name: build-viz
description: Build, rebuild or improve ONE interactive explainer for a chapter outside the full pipeline - storyboard it if needed (lesson-designer), build it (viz-builder), review it (viz-reviewer), embed it in the notebook if new, and commit. Usage - /build-viz 7 dispersion_relation   ·   /build-viz 7 dispersion_relation "make the particle orbits bigger on phones"   ·   /build-viz 7 new "Kelvin wave along a coast"
disable-model-invocation: true
---

# /build-viz N <slug | new> ["request"]
1. Read `progress.json` and `analysis/chNN_curation.md` §5. Count existing `viz/chNN/*.html` (a chapter always keeps 4–5); a new one must not exceed
   `book.yaml → project.max_explainers_per_chapter` — if it would, ask the user which one to replace.
2. **New explainer** (`new "<idea>"`): brief `concept-curator` to add a `### E<k> · <slug>` entry (CORE idea, why
   interactive, mirrored function) to the curation, then `lesson-designer` to add its storyboard to
   `analysis/chNN_design.md` Part B. If the mirrored fluidpy function does not exist, brief `concept-implementer` and then
   `math-verifier` for it (tests first) before building.
3. **Build**: brief `viz-builder` (brief header from `/do-chapter`) with the slug, storyboard heading, mirrored function,
   and the user's request verbatim. It iterates on `tools/viz_lint.py` + `tools/shot.py` until PASS.
4. **Review**: brief `viz-reviewer` limited to this slug (it appends/updates the slug's section in `reports/chNN_viz.md`).
   Must-fix → back to the builder (max 2 rounds).
5. **Embed** (new explainers only): brief `notebook-builder` to add the `nb.explainer(...)` block where the storyboard
   places it; `tools/run_notebook.py chNN`; `tools/embed_check.py chNN`.
6. Commit `chNN: viz — <slug>: <one line>`. If the chapter is published, run `/publish-chapter N` (a rebuilt existing
   explainer only needs the push — pages load the file directly — but re-run `tools/publish_notebook.py --site-only` so
   the gallery text is current).
Print the shot verdict, the parity table and the two best screenshot paths.
