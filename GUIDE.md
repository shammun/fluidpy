# Guide: Kundu's *Fluid Mechanics* → notebooks, Colab, web pages and interactive explainers, with Claude Code

**Project folder (Claude Code runs here):** `C:\Users\sislam27\Work\Climate Dynamics PHD\Fluid Dynamics`
**Public results:** repository `https://github.com/shammun/fluidpy` · site `https://shammun.github.io/fluidpy/` ·
explainer gallery `https://shammun.github.io/fluidpy/viz/` · one *Open in Colab* button per chapter page.

Part 0 explains the machine. Part 1 is the one-time setup. Part 2 is the everyday loop (one command per chapter).
Part 3 walks through the ten phases and what to look at. Part 4 is the prompt playbook. Part 5 covers the explainers
and the three deliverables in detail. Part 6 is GitHub/Colab, troubleshooting and reference tables.
Companion: **LESSONS-FROM-SEAICE.md** (what went wrong in the first book project and what now prevents it).

---

## 0. How the system works

For every chapter you get **one learning package**:

| Deliverable | Where | What is in it |
|---|---|---|
| Notebook | `notebooks/chNN_<slug>.ipynb` (executed) | plain-words explanations → maths step by step → tiny worked example → commented code → figures with "how to read it" → animations, plotly slider figures, live widgets → the chapter's explainers embedded |
| Colab notebook | `notebooks/chNN_<slug>_colab.ipynb` | the same, outputs stripped; its first cell clones the repo; explainers load from GitHub Pages |
| Web page | `notebooks/chNN_<slug>.html` | the executed notebook in the site theme; each explainer fills 100 % of the window's width and height; Open-in-Colab; prev/next; contents |
| Explainers (5–10) | `viz/chNN/<slug>.html` | walkthrough · explore (sliders, presets) · step-by-step working · equations (live numbers) · code · check yourself; linked views; fits any screen with no scrolling |

The **kit** that produces them:

| Component | Role |
|---|---|
| `CLAUDE.md` | the rules Claude Code reads every session (sequential chapters, equations from page images, teaching style, exhaustive coverage with tiered depth, 5–10 explainers that fit the window, evidence, public-repo safety, Windows specifics) |
| `book.yaml` | chapter map with exact PDF pages (from the PDF's own outline), sections, blurbs, explainer seeds; project URLs |
| `.claude/agents/` (12) | analyst · curator · lesson-designer · implementer · verifier · reviewer · viz-builder · viz-reviewer · notebook-builder · lesson-reviewer · knowledge-keeper · site-publisher |
| `.claude/skills/` (17) | 8 commands (`/do-chapter`, `/status`, …) + 9 background skills (teaching style, explainer engine, Python viz, maths→Python, verification, …) |
| `assets/viz_lib.js`, `viz_base.css` | the explainer engine: layout that fits any window, tabs, walkthrough, equations, plotting and flow-field tools, audit hooks |
| `fluidpy/core/` | `show_viz` (embedding), `anim` (animations), `interact` (plotly sliders), `style`, `project`, `units`, `refdata` |
| `tools/` | split/render the PDF · notebook kit · explainer scaffold/inline/lint · browser audit (`shot.py`) · publish · site · public-repo check · self-tests |
| `knowledge/`, `progress.json` | what each chapter taught and built (read by the next chapter) · phase tracker |

**Why agents.** Each agent works in a fresh context with a written brief and returns a short report; the main session
(the *orchestrator*) only coordinates, so a 70-page chapter never floods one context. Everything important is on disk,
so any session can resume.

**Coverage exhaustive, depth tiered; nothing used unexplained.** The analyst lists *everything* in the chapter (every
definition, theorem, numbered equation, example and figure), marks each as NEW or SEEN (taught in an earlier chapter) and
lists what it needs to be understood. The curator gives every item an ID and a depth in an auditable chapter map
(`book.yaml → policy`: `tier_a_full_treatment: [12, 18]`, `coverage: exhaustive`) — nothing is silently dropped and the
notebook covers **every book section**, but only the load-bearing ideas get the full treatment:

| Depth · tier | Who gets it | Treatment in the notebook |
|---|---|---|
| **A** · **CORE** `C01…` | **the 12–18 ideas that drag the rest along** (hard max 18; always the minority) | 🧩 block: picture → question → derivation step by step → worked number → commented code → **at least one visualization** (figure, animation, plotly slider or explainer) → how to read it; ≥ 2 kinds of test evidence |
| **B** · **NOTE** `N01…` | other new ideas, stated and explained | 📝 inside the nearest A block: a paragraph, the equation, a number where it helps; no separate derivation (the result is given) |
| **C** · **NOTE** `N01…` | new ideas that only need naming | one sentence inside the nearest A block, with a pointer to where the idea is used later |
| **RECAP** `R01…` (B or C) | ideas taught in an earlier chapter that this chapter uses | 🔁 a short reminder where it is needed, with a pointer to where it was taught |
| **SKIP** `S01…` (C) | only exercises, bibliography, material the book defers to a later chapter | one line in its section saying where the idea is covered |

**Hard maths is derived step by step, and explained.** Every result the chapter gets by manipulating equations is a
**derivation** `D01…` with a difficulty (★ short · ★★ several ideas · ★★★ long or conceptually hard). The analyst notes
which moves the book skips; the designer writes every derivation out in full; the notebook shows it as a 🧮 block:

```
What we want to show (plain words) · The plan (2–4 bullets) · Tools we use (each already explained) · We start from
Step k of n — <the move, e.g. "multiply every term by the velocity">
    the new line (one relation)
    • Why we can do this: the rule or assumption that allows it, and why we make this move
    • In words: what the new line says physically
Result · What it means · Check it (units, a limit, and a sympy cell for ★★★)
```

The explainer for that idea has a **Derivation** tab with the same steps: the line before (faded), the move, the new
line (highlighted), *Why*, *In words*, the line with **your numbers**, and each step moves the picture to the case it is
about; the last page shows the whole chain, the result and what it means right now. A walkthrough step can quote the key
step with an "All steps →" button. Every explainer also has a live **Explain** tab ("Explanation & interpretation"),
modelled on your `forced_damped_vibrations.html`: what each colour on screen is, every displayed number worked out from
your current settings in numbered sections with boxed results, the values at the current time, and a paragraph
interpreting the current regime (on a big screen, ↗ opens it in its own live browser tab).

The notebook kit refuses to save a notebook where a planned derivation is missing, a step lacks its *why* or *in words*,
or a ★★★ derivation has no sympy check; the browser audit visits every page of every derivation at 8 screen sizes; the
`math-verifier` re-derives ★★/★★★ results symbolically; the `lesson-reviewer` and `viz-reviewer` work through every
step and fail any line that does not follow from the one above by the stated move.

**Nothing is used unexplained:** every concept, symbol, maths tool (partial derivatives, Taylor series, complex
exponentials, …) and Python function (broadcasting, `np.meshgrid`, `solve_ivp`, …) the notebook uses is explained no later
than where it is first used — by its own CORE block, a recap, or a 📎 **primer** (plain words + a tiny runnable demo).
The designer keeps a *prerequisite ledger*; `knowledge/primers.md` stops later chapters repeating primers.

This is enforced mechanically: the notebook kit refuses to save a notebook with a missing section, a missing CORE/RECAP
block, or a CORE block without code or a visual; `tools/coverage_check.py` repeats the checks on the executed notebook
and checks the ledger; and the `lesson-reviewer` reads the notebook as a first-time learner and fails anything used
before it is explained. The **5–10 explainers** are chosen from the CORE ideas where manipulation or motion teaches what
a static figure cannot; every other CORE idea gets its figures, animations and plotly sliders.

**Correctness without an answer key.** The book has no code, so each function proves itself: V1 analytic solutions ·
V2 symbolic re-derivation and units · V3 convergence order · V4 conservation · V5 cited benchmarks · V6 the book's own
numbers (kept private) · V7 limits and symmetry. CORE items need two independent levels. Explainer JavaScript is
checked against the Python it mirrors (`selftest()` parity rows).

---

## PART 1 — ONE-TIME SETUP

### 1.1 What is already done in this folder
The kit is installed and was tested on this machine: `.venv` (Python 3.11.5 with numpy, scipy, sympy, matplotlib,
plotly, nbconvert, Playwright, imageio-ffmpeg), the `fluidpy-venv` Jupyter kernel, the book split into `chapters/`
(16 chapters + 4 appendices, every chapter verified to start on its opener page), the machinery self-tests passing
(unit tests, explainer audit at 8 screen sizes, end-to-end publish of a test notebook), and a local git repository.

### 1.2 Start Claude Code in the folder
```powershell
cd "C:\Users\sislam27\Work\Climate Dynamics PHD\Fluid Dynamics"
claude
```
Check the commands are visible:
```prompt
/status
```

### 1.3 If you ever need to rebuild the environment
```prompt
/setup-project
```
It is idempotent: it repairs the venv, kernel, split, hooks and runs the self-tests, then prints a checklist.

### 1.4 GitHub CLI (needed once, when chapter 1 is ready to go online)
```powershell
gh auth status
```

### 1.5 Google Drive — optional
Colab gets everything from GitHub. Drive is only needed if a chapter asks you to download data by hand.

### 1.6 Create the public repository and switch on GitHub Pages (once)
Do this after `/do-chapter 1` finishes locally (or before, if you prefer an empty site first). Claude Code will stop and
ask; run these yourself with a leading `!` so the output lands in the conversation:
```prompt
! gh repo create shammun/fluidpy --public --source . --remote origin --push --description "Kundu's Fluid Mechanics learned through Python: teaching notebooks, Colab, interactive explainers"
! gh api -X POST repos/shammun/fluidpy/pages -f "source[branch]=main" -f "source[path]=/"
```
Then tell Claude Code: `Continue publishing chapter 1.` Pages is live about a minute later.
If you choose another repository name, change `project.repo` and `project.site_url` in `book.yaml` first — every tool
reads the URLs from there.

---

## PART 2 — THE EVERYDAY LOOP

```prompt
/clear
/do-chapter 1
```
That is the whole workflow. The orchestrator runs ten phases (parallel where safe), commits after each, stops only on
a failed gate or an outward-facing action (first push), and ends with a ≤ 18-line report: the verdict, the teaching
spine, the explainers, runtime, URLs, and the two screenshots you should look at first.

Then study the chapter (page, Colab, explainers), and when you are ready:
```prompt
/clear
/do-chapter 2
```
Chapters are sequential: chapter N starts only when N−1 is published and its knowledge captured.

Useful options:

| You want | Type |
|---|---|
| to approve the explainer shortlist before building | `/do-chapter 5 --consult` |
| to resume after an interruption | `/do-chapter 5` (it resumes at the first unfinished phase) or `/do-chapter 5 --from viz` |
| to redo one phase only | `/do-chapter 5 --only notebook` |
| to publish locally without pushing | `/do-chapter 5 --no-push` |
| to see where you are | `/status` |

---

## PART 3 — THE TEN PHASES (what happens, what to look at)

| # | Phase (agent) | Reads | Writes | Look at |
|---|---|---|---|---|
| 1 | **analyze** (`concept-analyst`) | chapter text, rendered equation pages, knowledge | `analysis/chNN.md` | the dependency graph (§3), the derivations and the moves the book skips (§2b), the risks (§9) |
| 2 | **curate** (`concept-curator`) | analysis, `viz_seeds`, viz patterns | `analysis/chNN_curation.md` | the teaching spine, the count per tier, the derivations table (difficulty, where shown), the section-coverage table, the explainer shortlist with each "aha", the SKIP list with pointers |
| 3 | **design** (`lesson-designer`) ∥ | curation | `analysis/chNN_design.md` | explainer storyboards (steps, controls, Explain sections, equations) and Part F: every derivation written out step by step |
| 4 | **implement** (`concept-implementer`) ∥ | analysis, curation, rendered pages | `fluidpy/chNN_<slug>.py`, `scripts/` | docstrings cite section + equation + units |
| 5 | **verify** (`math-verifier` ⟲ implementer) | code, tiers | `tests/test_chNN.py`, `reports/chNN_verification.md` | Verdict, labels, convergence orders, Open items |
| 6 | **review** (`derivation-reviewer`) ∥ | code vs page images | `reports/chNN_review.md` | Must-fix list and that each was fixed |
| 7 | **viz** (`viz-builder` ×5–10 → `viz-reviewer`) ∥ | storyboards, fluidpy | `viz/chNN/*.html`, `reports/chNN_viz.md` | screenshots in `reports/viz/chNN/<slug>/` |
| 8 | **notebook** (`notebook-builder` → `lesson-reviewer`) ∥ | design Parts A/C/E/F, fluidpy | `notebooks/build_chNN.py`, `.ipynb`, `reports/chNN_lesson.md` | every CORE idea with code + visual, every derivation step-by-step and audited, nothing unexplained, runtime |
| 9 | **knowledge** (`knowledge-keeper`) ∥ | everything | `knowledge/*` | "Feeds forward" |
| 10 | **publish** (`site-publisher` + orchestrator) | notebook, explainers | executed `.ipynb`, `_colab.ipynb`, `.html`, index, gallery | the live page on your phone |

Gates: analysis covers every numbered equation · every item has an ID and a depth (A ≤ 18), every section covered, every derivation is a D row (A item or never written out by the book) or a stated result · 5–10 explainers + backup at reference depth · every CORE block has code + visual, every derivation is written out and checked, and nothing is used unexplained (coverage_check + lesson review) · every design function exists · tests
PASS (≤ 3 fix loops, never by loosening a tolerance) · no open Must-fix · every explainer passes lint + browser audit at
8 sizes + parity + visual review (≤ 2 rebuild rounds, else it is dropped, never published broken) · notebook executes
with 0 errors · embed check (every explainer embedded exactly once) · public-repo check · page audit (explainers
full-window).

---

## PART 4 — PROMPT PLAYBOOK

### 4.1 Decision table
| Situation | Prompt |
|---|---|
| start / continue the book | `/clear` then `/do-chapter N` |
| where am I? | `/status` |
| an explainer could be better | `/build-viz N <slug> "what to improve"` |
| you want another or a better explainer (5–10 per chapter) | `/build-viz N new "the idea"` |
| a notebook section is unclear | `/notebook-chapter N "add a tiny example for …"` |
| you changed physics code | `/verify-chapter N` then `/publish-chapter N` |
| just re-publish | `/publish-chapter N` |
| understand something right now | `/explain-concept N "…"` |
| verification stopped with FAIL | the root-cause prompt in 4.3 |
| Claude says a command was blocked | run it yourself with a leading `!` |

### 4.2 Rules of thumb
- **One chapter per session**; `/clear` before `/do-chapter`. Never paste book text into the chat.
- **Don't accept loosened tolerances** or "skip verification for now":
  ```prompt
  Do not loosen tolerances or skip a level. Report the discrepancy in Open items with your best hypothesis (sign
  convention, factor of 2, nu vs mu, reference scale, half-cell boundary, degrees vs radians) and stop.
  ```
- **Never use Colab's File → Save a copy in GitHub.** Save to Drive.

### 4.3 Follow-up prompts
**Verification failed:**
```prompt
Read reports/ch05_verification.md. For each failing item state the most likely cause from the standard list, test the
hypothesis with the smallest experiment, fix the code (not the tolerance) via the implementer, and re-run /verify-chapter 5.
```
**An explainer does not click:**
```prompt
/build-viz 7 dispersion_relation "Step 3 is confusing: show the crest and the wave packet moving at the same time, and
make the group-velocity dot bigger on phones. Keep the walkthrough to 6 steps."
```
**Teaching depth:**
```prompt
/notebook-chapter 4 "Section 4.6: derive the Navier-Stokes viscous term one index step at a time, with a 2-D
Couette example using easy numbers, and add a plotly slider for the pressure gradient."
```
**A derivation is still too hard to follow:**
```prompt
/notebook-chapter 7 "D04 (the dispersion relation): step 5 jumps from the surface condition to the tanh. Split it into
one move per step, say why each move is allowed and what each line means in words, add a sympy check of the
intermediate line, and copy the same steps into the Derivation tab of dispersion_relation (/build-viz 7 dispersion_relation)."
```
**Weak label:**
```prompt
In reports/ch06_verification.md the item "<name>" is labelled qualitative. Add a real assertion (exact solution,
manufactured-solution convergence study, or a cited benchmark). If none is possible, say precisely why.
```
**Study sheet after a chapter:**
```prompt
Write notes/ch07_studysheet.md: the chapter's argument on one page — assumptions, the load-bearing equations, for each
the fluidpy function, the test that proves it and the explainer step that shows it — then five questions I can answer by
running code or playing with the explainers.
```
**Keep the workflow improving:**
```prompt
Save what we just learned so the next chapter already knows it: update the relevant agent/skill file (Lessons section),
knowledge/viz_patterns.md or knowledge/CUMULATIVE.md, and show me the diff.
```

---

## PART 5 — THE DELIVERABLES IN DETAIL

### 5.1 How an explainer is built
1. The curator picks 5–10 CORE ideas (plus a backup), says why interaction beats a static figure, which depth features it
   needs, which derivations it steps through, and which **reference explainer** it follows — your own, with the three
   Unit 4 MIT-mathlet re-implementations preferred: `forced_damped_vibrations.html` (its live explanation panel is the
   model for every Explain tab), `amplitude_phase_second_order_II_3.html` (a numbered derivation with live numbers) and
   `angular_frequency_explorer_1.html` (linked views and modes); then the fast.ai labs `fid_formula_lab`,
   `forward_noising_lab`, `pixels_as_parameters`, `random_copy_lab`, `overfitting_curves`, `stride_padding_playground`,
   `ddpm3_unet_3d`, `np_resnet_3d` (paths in skill `interactive-viz` §5).
2. The lesson-designer storyboards it: linked views, controls, presets, status verdict, the Explain sections and their
   interpretation texts, the Derivation tab (what the picture shows at each step), code, 4–8 walkthrough steps,
   equations with book numbers and live substitutions, ≥ 3 questions, parity rows.
3. `tools/new_viz.py` scaffolds the file from `templates/viz_template.html`; the builder writes the physics and the
   `Viz.app({...})` configuration, copying patterns from the three passing references in `templates/`
   (`viz_example.html` — linked views, transport, presets, status, energy-term bars, inspector, notes, Explain tab,
   two derivations (why energy can only fall; where ω_d comes from), synced code; `viz_example_field.html` — flow fields and particles; `viz_example_3d.html` — a three.js scene).
4. `tools/viz_lint.py` (static rules) and `tools/shot.py` (headless Edge at 360×640, 390×844, 844×390, 768×1024,
   1280×720, 1000×700, 1366×768, 1920×1080; every tab, every step): fails on any overflow, a window not filled, text
   below 12 px, small tap targets, JS errors, missing equations, parity mismatch with the Python function, and on the
   **quality floor** — an Explain tab with ≥ 3 numbered sections and an interpretation, a Derivation tab for every
   derivation the storyboard lists (each step with its move, line, why and in-words), a Code tab, ≥ 2 depth features,
   ≥ 3 questions.
5. The builder and then the `viz-reviewer` *look at* the screenshots and compare with the reference. Nothing is
   published that fails; if one fails twice, the backup idea is built so the chapter still has at least 4.

**What every explainer offers the reader** (tabs): Walkthrough · Explore (controls, presets, live numbers, term bars,
inspector, "right now" interpretation) · Explain (every number on screen worked out with your current settings, then
what it means) · Derivation (the key formula built one move at a time, when the idea has one) · Equations · Code (the
Python, with the walkthrough lighting up the lines) · Check yourself. Above the picture: presets and a live
status verdict; below it: ▶ ⏮ ⏭ ↺, a time scrubber and speed.

**Why it never needs scrolling:** the engine measures itself. It chooses a layout for the window's shape, tightens
spacing in up to three steps, hides optional items, and splits long lists (controls, equations, questions) into pages
with ‹ 1/2 › buttons; the audit proves the result at every size.

### 5.2 How explainers appear in the three places
| Place | Mechanism | Size |
|---|---|---|
| Jupyter / VS Code | `show_viz("chNN", "<slug>")` → `<iframe srcdoc>` with the file from disk (works offline) | output width × visible window height; ⤢ Full screen; Open in new tab |
| Google Colab | `show_viz` → `<iframe src>` of the GitHub Pages copy (srcdoc from the clone if not yet published) | output width × screen height − 170 px; Full screen / new tab |
| Web page | the publish tool replaces the output with a block the page script moves out of the cell | exactly 100 % of the window width and height; the page's Contents button hides while it is on screen |

### 5.3 Python animations and interactive figures
- `show_animation(animate(update, frames=60, fig=fig))` → an MP4 video (compact) or `player="frames"` for a
  step-through player. Both play in Jupyter, Colab and the web page.
- `slider_figure(fn, "H", values, …)` and `animate_figure(frame_fn, times, …)` → plotly figures whose sliders work on
  the web page without Python (every slider position is precomputed).
- `live(fn, a=(0, 1, 0.05))` → ipywidgets for free exploration while a kernel runs; the web page shows a "run it in
  Colab" note in its place.

### 5.4 The web site
`index.html` lists all 16 chapters as cards (published ones with Read · Open in Colab · .ipynb and their explainers;
the rest "coming"). `viz/index.html` is the gallery of every explainer, each opening full-window or at its place in
the chapter. Everything is mobile-first.

---

## PART 6 — GITHUB, COLAB, TROUBLESHOOTING, REFERENCE

### 6.1 Colab
Open a chapter page → **Open in Colab** → run the ⚙️ Setup cell (clones the repo, installs `pint`) → run cells top to
bottom. Explainer cells show the published explainers. Do not use *Save a copy in GitHub*: it commits outputs into the
public repository. If it happens anyway:
```prompt
Check origin/main for a "Created using Colab" commit. If one exists, rebase onto it, restore the outputs-stripped
_colab.ipynb by re-running tools/publish_notebook.py, commit, run tools/check_public.py, and tell me what remains in history.
```

### 6.2 Troubleshooting
| # | Symptom | Fix |
|---|---|---|
| T1 | `/do-chapter` not found | start `claude` inside the project folder; `.claude/` must exist |
| T2 | an equation in the analysis looks wrong | `/explain-concept N "(N.M)"`; agents must read `tools/render_pages.py chNN --eq N.M` images, never the text |
| T3 | verification loops and stops | prompt 4.3; typical causes: sign convention, factor 2, ν vs μ, reference scale, half-cell BC, degrees vs radians |
| T4 | an explainer fails the audit repeatedly | `/build-viz N <slug> "…"` with a simpler stage or a custom tab; the engine's rules are in skill `interactive-viz` |
| T5 | notebook execution error only in the publish step | `.venv/Scripts/python.exe tools/run_notebook.py chNN --save` and read the error |
| T6 | `jupyter nbconvert` fails with an unrelated binary error | it ran Anaconda's; always use the tools (they use the `fluidpy-venv` kernel) |
| T7 | explainers look wrong only on the web page | `.venv/Scripts/python.exe tools/shot.py --page notebooks/chNN_<slug>.html` and look at `reports/viz/_pages/` |
| T8 | machinery doubts after editing a tool | `pytest -q tests/test_machinery.py`, `tools/shot.py templates/viz_example.html`, `tools/pipeline_selftest.py` |
| T9 | `git push` rejected after using Colab | 6.1 |
| T10 | a command is "blocked by the classifier" | run it yourself with a leading `!` |
| T11 | page not updated after push | Pages takes ~1 min; hard-refresh; check the Actions tab on GitHub |
| T12 | garbled Unicode in the terminal | the project sets `PYTHONUTF8=1` in `.claude/settings.json`; restart Claude Code |

### 6.3 Reference — agents
| Agent | Job | Writes |
|---|---|---|
| concept-analyst | full inventory from text + page images, dependency graph, validation plan | `analysis/chNN.md` |
| concept-curator | IDs + depth (A 12–18 CORE / B / C, nothing dropped), derivations (D rows), teaching order, section coverage, primers needed, 5–10 explainers + backup, animations, interactives | `analysis/chNN_curation.md` |
| lesson-designer | notebook storyboard + explainer storyboards + function contract + prerequisite ledger + every derivation step by step | `analysis/chNN_design.md` |
| concept-implementer | physics as documented functions; fixes | `fluidpy/`, `scripts/` |
| math-verifier | tests on the evidence ladder; report | `tests/`, `reference/`, `reports/chNN_verification.md` |
| derivation-reviewer | fresh-eyes code-vs-book review (read-only) | → `reports/chNN_review.md` |
| viz-builder | one explainer, iterated to PASS | `viz/chNN/<slug>.html` |
| viz-reviewer | audits, screenshots, parity, teaching quality | `reports/chNN_viz.md` |
| notebook-builder | teaching notebook via nbkit (CORE blocks, recaps, primers, derivations with sympy checks), executed, coverage-checked | `notebooks/build_chNN.py`, `.ipynb` |
| lesson-reviewer | reads the executed notebook as a first-time learner: coverage, nothing unexplained, every derivation step follows, correctness, style | `reports/chNN_lesson.md` |
| knowledge-keeper | memory + promotion of patterns and helpers | `knowledge/`, skill Lessons |
| site-publisher | embed check, publish, page audit, public check | pages, Colab twin, index, gallery |

### 6.4 Reference — commands and tools
| Command / tool | Purpose |
|---|---|
| `/setup-project` · `/status` | environment + checklist · where you are and what next |
| `/do-chapter N` | the full pipeline |
| `/build-viz` · `/notebook-chapter` · `/verify-chapter` · `/publish-chapter` | redo one part |
| `/explain-concept N "…"` | learn one thing now |
| `tools/split_pdf.py`, `tools/render_pages.py`, `tools/detect_chapters.py` | book → chapters; equation page images; outline check |
| `tools/nbkit.py`, `tools/run_notebook.py` | notebook building blocks; headless execution report |
| `tools/new_viz.py`, `tools/viz_inline.py`, `tools/viz_lint.py`, `tools/shot.py` | explainer scaffold, library inlining, static lint, browser audit |
| `tools/coverage_check.py`, `tools/embed_check.py`, `tools/publish_notebook.py`, `tools/build_site.py` | teaching-coverage gate, integration gate, publishing, site pages |
| `tools/check_public.py` (+ pre-push hook) | nothing book-derived goes public |
| `tools/pipeline_selftest.py`, `tests/test_machinery.py` | prove the machinery works |
| `tools/convergence.py`, `tools/benchmarks.py`, `tools/compare_fields.py` | verification helpers |
| `tools/build_guide_docx.py` | regenerates `GUIDE.docx` from this file |

### 6.5 Reusing the kit for another book
Copy the folder without `chapters/`, `analysis/`, `fluidpy/ch*`, `viz/ch*`, `notebooks/`, `tests/test_ch*`, `reports/`,
`knowledge/ch*`; replace the PDF; regenerate `book.yaml` from the outline (`tools/detect_chapters.py` shows it); edit
`project` in `book.yaml`, the "book in one paragraph" in skill `fluids-book`, and the benchmark catalogue in
`verify-implementation`. Everything else — agents, explainer engine, publishing, audits — is book-independent.
