---
name: notebook-chapter
description: Rebuild and headlessly execute one chapter's teaching notebook with the notebook-builder agent (e.g. after physics fixes, a new explainer, or teaching feedback), optionally applying the user's requested changes to the storyboard first. Usage - /notebook-chapter 7   ·   /notebook-chapter 7 "add a tiny example for group velocity"
disable-model-invocation: true
---

# /notebook-chapter N ["change request"]
1. If a change request is given, append it to `analysis/chNN_design.md` Part A under "## Revisions (<date>)" in one or
   two lines (you, not an agent).
2. Brief `notebook-builder` exactly as phase 8 of `/do-chapter` (brief header from that skill) plus the change request.
3. Gate: `.venv/Scripts/python.exe tools/run_notebook.py chNN` prints OK; `.venv/Scripts/python.exe tools/embed_check.py
   chNN` prints OK (if all explainers exist).
4. `progress.json → notebook: pass`, commit `chNN: notebook — …`.
5. If the chapter was already published, run `/publish-chapter N` so the page, the executed `.ipynb` and the Colab twin
   match. Print notebook path, runtime, cells by kind, figures.
