---
name: setup-project
description: One-time (idempotent) bootstrap for fluidpy - verify the book PDF and the chapter map in book.yaml, create or repair the Python venv with every dependency (incl. Playwright using the installed Edge/Chrome and imageio-ffmpeg), register the fluidpy-venv kernel, split the book into chapters/, run the machinery self-tests, initialise git with the public-repo guard, and print a checklist. Usage - /setup-project
disable-model-invocation: true
---

# /setup-project

Idempotent: every step checks before acting. Print a ✅/❌ checklist at the end. Do not start chapter work here.

## 1. Layout
- `pwd` is the project root (it contains `CLAUDE.md`, `book.yaml`, `.claude/`, exactly one large `*.pdf`).
- `.venv/Scripts/python.exe tools/split_pdf.py --check` → "check OK: every chapter starts on its CHAPTER opener page".
  (`book.yaml` page ranges come from the PDF's own outline; if the check fails, stop and show the user the mismatch.)

## 2. Python environment
- If `.venv` is missing: `python -m venv .venv` (Python ≥ 3.11).
- `.venv/Scripts/python.exe -m pip install -U pip -r requirements.txt` (includes plotly, playwright, imageio-ffmpeg).
- Browser for audits: `tools/shot.py` uses the installed Microsoft Edge or Google Chrome through Playwright (no
  download). Only if neither exists: `.venv/Scripts/python.exe -m playwright install chromium`.
- Kernel: `.venv/Scripts/python.exe -m ipykernel install --sys-prefix --name fluidpy-venv --display-name "Python (fluidpy .venv)"`.
- Record versions into `progress.json → environment` (python, numpy, scipy, sympy, plotly, nbconvert, playwright,
  browser channel found, checked_on) and `notes`: "equations must be read from rendered pages (tools/render_pages.py)".

## 3. Split the book (local, git-ignored)
`.venv/Scripts/python.exe tools/split_pdf.py --appendices` → `chapters/chNN.{pdf,txt}`, `chapters/_page_map.json`.
Spot-check: the first lines of `chapters/ch01.txt` show "C H A P T E R / 1 / Introduction" and the last page of
`chapters/ch01.txt` is chapter 1's literature list, not chapter 2.

## 4. Machinery self-tests (all must pass)
- `.venv/Scripts/python.exe -m pytest -q tests/test_machinery.py`
- `.venv/Scripts/python.exe tools/viz_inline.py --all --check`
- `.venv/Scripts/python.exe tools/shot.py templates/viz_example.html --quick`
- `.venv/Scripts/python.exe tools/pipeline_selftest.py`

## 5. Git and the public-repo guard
- `git init -b main` if needed; `git config core.hooksPath .githooks`.
- `git add -A`, then `.venv/Scripts/python.exe tools/check_public.py` → OK, then
  `git commit -m "setup: fluidpy kit — agents, skills, explainer engine, publishing machinery"`.
- **Do not** create the GitHub repository or push. Tell the user the one-time commands (GUIDE.md §1.6), to run with a
  leading `!` when they are ready (normally after chapter 1 is published locally):
  `! gh repo create shammun/fluidpy --public --source . --remote origin --push` and
  `! gh api -X POST repos/shammun/fluidpy/pages -f "source[branch]=main" -f "source[path]=/"`.

## 6. Checklist + next step
Print the checklist, then: "Next: `/do-chapter 1`".
