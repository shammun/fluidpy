---
name: status
description: Print fluidpy progress per chapter and phase (analyze, curate, design, implement, verify, review, viz, notebook, knowledge, publish), the current chapter's open items and explainer verdicts, environment and git state, and the exact next command to run. Usage - /status
disable-model-invocation: true
---

# /status

1. Read `progress.json`. Print one compact table:
   `ch | ana cur des imp ver rev viz nb know pub | explainers | notes` using ✅ pass/done · ⏳ running · ❌ fail ·
   ⛔ blocked · · todo. A chapter is finished only when `publish` is ✅.
2. For the current (first unfinished) chapter, if they exist: the Verdict line and "Open items" of
   `reports/chNN_verification.md`; the summary table of `reports/chNN_viz.md`; the Must-fix list of
   `reports/chNN_review.md` — each trimmed to ≤ 8 lines.
3. Environment: `.venv/Scripts/python.exe --version`, numpy/scipy/sympy/plotly/playwright versions (one line), whether
   `chapters/` is split (count of `chNN.txt`), the `fluidpy-venv` kernel exists
   (`.venv/Scripts/python.exe -m jupyter kernelspec list` contains it — informational).
4. Git: `git log --oneline -5`, working tree clean or not, `origin` configured or not, ahead/behind `origin/main`.
5. Site: if `origin` exists, `curl -s -o /dev/null -w "%{http_code}"` on the Pages root (from `book.yaml → project.site_url`).
6. Next command, one of:
   - `/setup-project` — if `.venv` or `chapters/` or git is missing;
   - `/do-chapter N --from <first phase not pass>` — for an interrupted chapter;
   - `/clear` then `/do-chapter N+1` — when the last chapter is fully published;
   - the exact `! gh …` commands from GUIDE.md §1.6 — if chapter 1 is ready to publish but there is no `origin`.
