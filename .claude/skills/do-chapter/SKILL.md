---
name: do-chapter
description: Run the whole fluidpy pipeline for one book chapter autonomously - analyze, curate (every new idea is CORE with code and a visual), design, implement, verify, review, 4–5 interactive explainers, notebook with lesson review, knowledge, publish - with parallel subagents where safe, gates, fix loops, progress tracking and a commit per phase. Usage - /do-chapter 7   ·   /do-chapter 7 --from viz   ·   /do-chapter 7 --only notebook   ·   /do-chapter 7 --consult
disable-model-invocation: true
---

# /do-chapter N [--from PHASE] [--only PHASE] [--consult] [--no-push]

You are the **orchestrator**. You delegate; you do not write physics, tests, notebooks or explainers yourself. Keep your
context small: read `progress.json`, `book.yaml`, agent *replies* and the headers/verdict lines of artifacts — not whole
files. Every brief must be self-contained (a subagent cannot see this conversation) and must end with:
*"Do not return while a background run is still in progress; wait for it and report the real numbers."*

Phases (keys in `progress.json`): `analyze → curate → design ∥ implement → verify → review ∥ viz ∥ notebook →
knowledge ∥ publish`. Status values: `todo`, `running`, `pass`, `fail` (+ `notes`), `blocked` (waiting on the user);
`knowledge` uses `done`.

---------------------------------------------------------------------------------------------------------------------
## 0. Preflight (no agents)
1. Parse N → `chNN`; read its `book.yaml` row: `number, title, slug, pdf_start/end, printed_start, sections, blurb,
   viz_seeds`. Read `progress.json`.
2. **Sequential rule**: the previous `do_chapter: true` chapter must have `publish: pass` and `knowledge: done`, unless
   the user said to override. Otherwise stop and say what is missing.
3. **Environment**: `.venv/Scripts/python.exe` exists (else tell the user to run `/setup-project`); `git rev-parse` works;
   `chapters/chNN.txt` exists (else run `.venv/Scripts/python.exe tools/split_pdf.py chNN` yourself — it only writes
   git-ignored local files).
4. **Resume**: with `--from P`, every phase before P must be `pass`/`done`; with `--only P`, run just that phase and its
   gate. Without flags, start at the first phase that is not `pass`/`done` (a `running` phase from a crashed session is
   re-run from its start — its artifacts on disk are inputs, not results).
5. Print the plan (phases to run, parallel groups, expected artifacts) in ≤ 12 lines. Continue without waiting unless
   `--consult`.

**Brief header** (paste at the top of every brief, filled in):
```
Project: fluidpy (Kundu, Cohen & Dowling, Fluid Mechanics 5e → tested Python, teaching notebook, interactive explainers).
Repo root: C:\Users\sislam27\Work\Climate Dynamics PHD\Fluid Dynamics  (Windows, path has spaces; run Python as .venv/Scripts/python.exe;
write code files with Write/Edit, not shell heredocs). Rules: CLAUDE.md. Do not commit.
Chapter: chNN — N. <title>  · slug <slug> · sections <…> · PDF pages <s>-<e> (printed = pdf - 27)
Chapter text: chapters/chNN.txt (maths garbled — render equation pages: tools/render_pages.py chNN --eq N.M, then Read the PNG)
Knowledge to read first: knowledge/CUMULATIVE.md, knowledge/notation.md, knowledge/viz_patterns.md, knowledge/concept_map.md
```

After each phase: update `progress.json` (status + one-line `notes` + `updated` date), then commit:
`git add -A && git commit -m "chNN: <phase> — <one line>"` (`.githooks/pre-push` guards pushes; commits are local).

---------------------------------------------------------------------------------------------------------------------
## 1. ANALYZE — `concept-analyst` → `analysis/chNN.md`
Brief: header + the output format lives in the agent definition; add the chapter's `blurb`.
**Gate**: every numbered equation appears once — compare the count of lines that are exactly `(N.M)` in
`chapters/chNN.txt` (`grep -cE "^\(N\.[0-9]+\)$"`) with the inventory rows mentioning `Eq. (N.`; the dependency graph
(§3) exists; CODE rows have validation plans. One send-back allowed. → `analyze: pass`, commit.

## 2. CURATE — `concept-curator` → `analysis/chNN_curation.md`
Brief: header + `book.yaml` `viz_seeds` for the chapter (as priors) + "4–5 explainers plus one backup".
**Gate** (send back once on any failure):
- every inventory row has an ID and a tier (CORE `C..` / RECAP `R..` / NOTE `N..` / SKIP `S..`);
- **every row the analyst marked NEW is CORE** (NOTE only when a parent CORE ID is named; SKIP only for history,
  exercises, bibliography or explicitly deferred material) — count NEW rows in `analysis/chNN.md` and compare;
- every SEEN row is RECAP (or NOTE inside a CORE block);
- the §3 section-coverage table has one row for every section in `book.yaml → sections`, none empty;
- §4 lists the prerequisites needing primers; §5 has 4–5 explainers + 1 backup, each with CORE IDs, why-interactive,
  mirrored fluidpy function, depth features (calc + code + ≥ 2 more) and a named reference explainer.
Print the counts per tier, the teaching order (one line per CORE item, grouped by section) and the explainer list
(slug + aha) to the user.
With `--consult`: stop here and ask the user to approve/edit the shortlist (`blocked`). → `curate: pass`, commit.

## 3 ∥ 4. DESIGN ∥ IMPLEMENT (launch both agents in ONE message)
- `lesson-designer` → `analysis/chNN_design.md` (brief: header; a storyboard block for EVERY CORE ID with code and a
  visual; storyboards for every explainer with headings `### E<k> · <slug>` and the backup `### B1 · <slug>`; Part C =
  function contract; Part E = prerequisite ledger).
  **Design gate**: the number of CORE blocks in Part A equals the curation's CORE count; Part E exists and every row has
  an "Explained by"; every explainer storyboard lists calc, code and ≥ 2 depth features. One send-back allowed.
- `concept-implementer` → `fluidpy/`, `scripts/` (brief: header; implement analysis §4 + curation §9; "the lesson-designer
  is writing analysis/chNN_design.md in parallel — if Part C appears before you finish, honour its names and signatures").
**Merge gate** (after both return): every Part C function exists with a compatible signature
(`grep -n "def <name>" fluidpy/`). Missing → one short follow-up brief to the implementer with just those functions.
All scripts exit 0 (from the implementer's reply; spot-check one yourself). → `design: pass`, `implement: pass`, commit.

## 5. VERIFY — `math-verifier` → `tests/test_chNN.py`, `reference/chNN/`, `reports/chNN_verification.md`
Brief: header + tiers from curation §2 + design Part C list.
**FAIL loop** (max 3): send the failing items + report path to `concept-implementer` ("fix physics, never tolerances,
keep public signatures"), then re-run `math-verifier`. Still failing → `verify: fail`, notes, **stop and report**.
Never accept a relaxed tolerance as a fix. → `verify: pass`, commit.

## 6 ∥ 7 ∥ 8. REVIEW ∥ VIZ ∥ NOTEBOOK (launch together in ONE message)
- `derivation-reviewer` (read-only; brief: header + file list + verification report path).
- `viz-builder` — **one agent per explainer**, each brief naming its single `viz/chNN/<slug>.html`, its storyboard
  heading, its mirrored function and the order number. On a laptop launch at most `project.max_parallel_viz`
  (default 3) at a time; start the next as one finishes.
- `notebook-builder` (brief: header + design Parts A/C/E + list of explainer slugs + "explainers are being built in
  parallel; embed by slug; run tools/run_notebook.py and tools/coverage_check.py until clean").
Folders never overlap: reviewer writes nothing, each viz-builder one file, notebook-builder `notebooks/`.

**When the notebook-builder returns**: `lesson-reviewer` → `reports/chNN_lesson.md` (reads the executed notebook as a
first-time learner: coverage, CORE blocks complete with code + visual, nothing used before it is explained, style,
correctness). Must-fix → re-brief `notebook-builder` with the report (max 2 rounds), then `lesson-reviewer` again.

**When the reviewer returns**: write `reports/chNN_review.md` yourself from its reply (Must fix / Should fix / Verified).
Must-fix code items → `concept-implementer` (fluidpy only; signatures stable), then `.venv/Scripts/python.exe -m pytest
tests/test_chNN.py -q` and, if tests had to change, `math-verifier` for those tests. Doc findings you fix yourself.
**When all viz-builders return**: `viz-reviewer` → `reports/chNN_viz.md`. Must-fix items → re-brief the specific
`viz-builder` with its list (max 2 rounds per explainer), then `viz-reviewer` again. An explainer that still fails is
removed (delete the file, note it in `knowledge/viz_patterns.md` as a failed idea) and — if fewer than 4 remain — the
**backup** (`### B1`) is built by a `viz-builder` and reviewed the same way; the notebook-builder swaps the embed. Never
publish a failing explainer and never publish fewer than 4.
**Merge gate** (all done):
`.venv/Scripts/python.exe -m pytest -q` (full) · `.venv/Scripts/python.exe tools/shot.py --chapter chNN --quick`
(parity after any physics fix) · `.venv/Scripts/python.exe tools/run_notebook.py chNN --save` ·
`.venv/Scripts/python.exe tools/coverage_check.py chNN --nb outputs/chNN/executed.ipynb` · `.venv/Scripts/python.exe
tools/embed_check.py chNN` (4–5 explainers, each embedded once) · `reports/chNN_lesson.md` and `reports/chNN_viz.md`
Verdict PASS. All clean → `review: pass`, `viz: pass` (+ `explainers: [slugs]`), `notebook: pass`, commit.

## 9 ∥ 10. KNOWLEDGE ∥ PUBLISH (launch together)
- `knowledge-keeper` (brief: header + report paths + "library/skill promotion is NOT allowed in this run; list
  promotion candidates in knowledge/viz_patterns.md instead" — the publisher is reading viz/ at the same time).
- `site-publisher` (brief: header + "run tools/embed_check.py, tools/publish_notebook.py chNN, tools/shot.py --page,
  tools/check_public.py; report the JSON").
Then, if promotion candidates were listed and are clearly worth it: brief `knowledge-keeper` again with "promotion
allowed" (serial; it re-inlines and re-audits).
**Publish gate** (yours):
1. `git add -A`, then `.venv/Scripts/python.exe tools/check_public.py` → OK. Never stage `chapters/`, `outputs/`,
   `tests/book_values_*`, `reports/viz/`.
2. Commit `chNN: publish — <n> explainers, notebook <runtime>s`.
3. Remote: if no `origin`, stop and ask the user to create the repo (exact commands in GUIDE.md §1.6) — do not create it.
   Otherwise `git fetch origin`; if `origin/main` has commits you lack, inspect them (a "Created using Colab" commit
   means outputs were saved from Colab: rebase, re-run `tools/publish_notebook.py chNN`, commit
   `chNN: restore outputs-stripped _colab.ipynb`, tell the user). Never force-push.
4. Unless `--no-push`: `git push` (if the classifier blocks it, ask the user to run `! git push`).
5. Poll `curl -s -o /dev/null -w "%{http_code}" <page_url>` until 200 (Pages takes ~1 min; give up after 5 min and say
   so). Check one explainer URL and the Colab URL (`https://colab.research.google.com/github/…` — a 200 from the GitHub
   blob URL of the `_colab.ipynb` is sufficient evidence).
→ `knowledge: done`, `publish: pass`, commit progress.json (and push it).

---------------------------------------------------------------------------------------------------------------------
## Final report to the user (≤ 18 lines)
- Verdict + validation counts per label; convergence orders; benchmarks used.
- Counts per tier (CORE / RECAP / NOTE / SKIP) and primers; the CORE items (one line each, grouped by section); the explainers (slug · aha · shot PASS · depth features).
- Notebook runtime; page size; animations and interactive figures count.
- URLs: chapter page, explainer gallery, Colab.
- Open items / things the user should look at first (the 2 best screenshots under `reports/viz/`).
- "Next: `/clear` then `/do-chapter N+1`."

## Orchestration rules (from experience)
- One chapter per session; `/clear` between chapters — everything needed is on disk.
- Parallel only inside the ∥ groups above. Never two agents writing the same folder.
- A failed gate leaves `fail` + notes and stops the pipeline; nothing downstream runs on unverified code.
- An agent that returns placeholders or "still running" is re-briefed to wait and report real numbers.
- If a machinery tool itself misbehaves (shot.py, publish, embed), run `tools/pipeline_selftest.py`; fix the tool (small,
  careful), re-run the self-test, commit `tools: …`, then continue.
