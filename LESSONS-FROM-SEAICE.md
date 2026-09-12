# What went wrong the first time, and what in this kit prevents it

The kit is the sea-ice project's workflow after nine chapters of use. This file records the actual failures of that
project — what happened, what the diagnosis was, what fixed it — and points at the file in **this** kit that carries the
fix forward. Where a lesson does not transfer (because that book shipped MATLAB code and this one ships none), it says
what replaces it.

Source project: *Sea Ice Image Processing with MATLAB* (Zhang & Skjetne, CRC Press 2018) → `github.com/shammun/seaice-py`,
chapters 2–10 ported, verified against MATLAB R2025a, published to GitHub Pages between 8 and 11 September 2026.

---

## A. Incidents, in the order they hurt

### A1 — The automatic PDF page offset was wrong, and looked right
**What happened.** The splitter auto-detected a front-matter offset of −3 by matching chapter titles. The true offset
was 31. With −3, `chapters/ch02.txt` contained the table of contents, and the analyst would have "analysed" the wrong
pages.
**Diagnosis.** Title matching also hits the table of contents and cross-references. The only reliable signal is the
page number *printed in the running head*.
**Fix.** Score candidate offsets by whether the page carries the chapter's printed start number as a standalone line at
the top or bottom, and always verify by eye.
**In this kit.** `tools/detect_chapters.py` measures the offset from printed page numbers across the whole book (modal
vote) and writes it to `book.yaml → book.pdf_offset`; `tools/split_pdf.py` prefers that value over re-guessing; the
`/setup-project` skill makes "open two chapters and check the first and last lines" a mandatory step; troubleshooting
T2 in GUIDE.md. **This is the single most valuable 90 seconds of the whole setup** — every downstream artifact is built
on the split.

### A2 — A "pre-check" that matched real data at zero error was still wrong
**What happened (chapter 4).** The analyst probed how MATLAB's edge thinning pads the image, concluded "replicate
padding", and confirmed it on a 12-megapixel photograph with **0 differing pixels**. It was wrong: MATLAB zero-pads.
Only constructed fixtures (ties, borders, extreme dtypes) exposed it.
**Diagnosis.** Real data does not exercise the branches where implementations differ. Agreement on a natural image is
weak evidence; it means the code paths were never visited.
**Fix.** Label pre-checks as hypotheses, and verify with constructed inputs that deliberately hit the boundary, the tie
and the extreme value.
**In this kit — and this is the deepest transfer.** The equivalent of "a real photograph" here is "a plausible-looking
flow field": a solution can look perfect and still be first-order accurate, off by a factor of 2, or right only in the
interior. That is exactly why the `verify-implementation` skill demands **two independent evidence levels** and why
`tools/convergence.py` exists: a manufactured solution and an order-of-accuracy measurement visit the branches that a
pretty contour plot never does. The `derivation-reviewer` agent is told in as many words that a passing test suite is
not evidence when the tests and the code share the same misreading.

### A3 — Colab's "Save a copy in GitHub" published copyrighted figures
**What happened (9 September 2026).** A notebook was opened in Colab, run against the *private* book image, and saved
to GitHub from Colab's File menu. That writes the **cell outputs** into `notebooks/chNN_<slug>_colab.ipynb` on `main` —
so figures derived from the book's image landed in a public repository, in a commit titled "Created using Colab".
**Diagnosis.** The publish flow assumed the only writer to `main` was the local machine.
**Fix.** (a) a guard that refuses any tracked notebook with outputs; (b) "fetch and inspect `origin/main` before every
push"; (c) restore the stripped notebook by re-running the publish tool; (d) tell the user the orphaned commit remains
fetchable by SHA until GitHub garbage-collects it — and never force-push for them.
**In this kit.** `tools/check_public.py` + `.githooks/pre-push` (install with `git config core.hooksPath .githooks`);
step 4 of the PUBLISH phase in `do-chapter`; the hazard section of the `colab-notebook` skill; §4.2 of GUIDE.md.
Here the private material is not images but the book's own numbers: they live in `tests/book_values_chNN.json`
(git-ignored), the notebook prints the banner `book values active` when it runs with them, and `check_public.py` fails
on any tracked notebook or HTML page carrying that banner.

### A4 — The permission classifier blocked the outward-facing commands
**What happened.** `gh repo create … --public --push` and some `git push` invocations were refused by Claude Code's
permission classifier, stalling the publish phase.
**Fix.** Have the *user* run the exact command by typing it with a leading `!` in the Claude Code prompt; the output
lands in the conversation and the session continues.
**In this kit.** GUIDE.md §2.2 and §4.1 give the exact commands to paste; the `do-chapter` publish phase says to ask
the user rather than retry; troubleshooting T14.

### A5 — `python -m jupyter nbconvert` ran somebody else's Python
**What happened.** Notebook execution failed with a scikit-image binary error that had nothing to do with the notebook.
`python -m jupyter nbconvert` dispatched to Anaconda's `jupyter-nbconvert.exe`, which was earlier on `PATH`.
**Fix.** Always `.venv\Scripts\python.exe -m nbconvert …`.
**In this kit.** Spelled out in the `colab-notebook` skill and the `notebook-builder` agent; `publish_notebook.py` calls
the `nbconvert` API in-process; troubleshooting T6.

### A6 — An agent returned a report with placeholder numbers
**What happened.** A verifier started pytest in the background and returned its report before the run finished; the
report contained placeholders and had to be resumed.
**Fix.** Every agent brief ends with "do not return while a background run is in progress; wait for it and report the
real numbers."
**In this kit.** That sentence is mandated in the `do-chapter` skill (top of the file) and repeated in the
`concept-implementer` and `math-verifier` agents.

### A7 — Two agents editing the same folder
**What happened.** Parallelising phases to save time caused conflicting edits when two agents both wrote `seaice/`.
**Fix.** A short list of *safe* parallel pairs: reviewer ∥ notebook-builder (both read-only w.r.t. the package),
knowledge-keeper ∥ the publish runs. Never two writers of the same folder.
**In this kit.** The "Orchestration lessons" section of the `do-chapter` skill.

### A8 — The reviewer cannot write, so its findings evaporated
**What happened.** `port-reviewer` is read-only by design; its findings lived only in a reply and had to be re-derived.
**Fix.** The orchestrator writes `reports/chNN_review.md` from the reply, then dispatches: code findings to the
implementer, test findings to the verifier, documentation findings it fixes itself.
**In this kit.** Phase 4 of `do-chapter` says exactly that.

### A9 — Hunting for substitute data took four attempts
**What happened.** Registering a public-domain image to replace each book image needed four tries (cloudy or
featureless satellite scenes) — and each candidate had to be *looked at* before being accepted.
**Fix.** Verify data with a real request and human inspection before relying on it; record provenance, licence and the
verification date.
**In this kit.** The analogue is a benchmark number, and the rule is stronger: `data-and-benchmarks` and
`verify-implementation` both say **never quote a benchmark from memory** — fetch or cite the source, record it in
`reference/chNN/SOURCES.md` with the date, and mark "digitised from a plot" where that is what happened.
`tools/benchmarks.py` ships the usual textbook constants **flagged `VERIFY`** so nobody mistakes them for checked values.

### A10 — The same code written twice
**What happened.** Chapters 6 and 7 shipped byte-identical MATLAB folders, and chapter 5 duplicated chapter 2's chain
codes. Without a check, the same algorithm would have been ported three times, three ways.
**Fix.** "Reuse before re-implementing": check `core/` and the function map first; anything used by two chapters is
promoted to `core/`.
**In this kit.** CLAUDE.md rule 9, the `fluids-book` skill's `fluidpy/core/` layout, and `knowledge/concept_map.md`,
which the analyst must read *before* planning. In a maths book the duplication is subtler — the same operator appears
as "the Laplacian" in one chapter and "the diffusion term" in the next — which is why the concept map is keyed by
concept, not by file name.

### A11 — Pressure to loosen a tolerance
**What happened.** When a parity test failed, the tempting fix was a looser tolerance.
**Fix.** A standing rule: fix the code, or report the discrepancy in Open items with a hypothesis; never relax the
test. Max three fix loops, then stop and tell the user.
**In this kit.** CLAUDE.md rule 5, `verify-implementation` ("Failure loop"), the `do-chapter` verify gate, and the
ready-made refusal prompt in GUIDE.md §2.2.

### A12 — Assumptions about the machine
**What happened.** The plan mirrored the repo to `G:\My Drive\seaice-py`. On that PC, `G:` is an external disk and
Google Drive for Desktop was not installed; the sync script was written and never used. Separately, the Chrome
extension was assumed for screenshots and was not connected (headless Chrome was used instead).
**Fix.** Check the environment; do not build steps on an unverified assumption.
**In this kit.** Drive is optional and only for Tier-4 manual data (GUIDE.md §1.7); `/setup-project` records what it
actually found in `progress.json → environment`.

---

## B. What does **not** transfer, and what replaces it

| Sea-ice project | Why it does not apply here | Replacement in this kit |
|---|---|---|
| MATLAB R2025a as the reference engine (`tools/run_matlab_ref.py`) | The fluid-dynamics book ships no code — there is no oracle to diff against | The V1–V7 evidence ladder in `verify-implementation`; `tools/convergence.py`, `tools/benchmarks.py`, `tools/compare_fields.py` |
| Parity labels `exact / near / approx / reimplemented` | Nothing to be exact *to* | Validation labels `analytic / symbolic / converged / conserved / benchmark / book-value / qualitative / unverified` |
| `analysis/_matlab_inventory.md`, "every `.m` file has one Python counterpart" | No `.m` files | The **concept inventory**: every definition, theorem, numbered equation, example and quoted number classified CODE/DEMO/SKIP — the analyst's list is the specification |
| Image-data tiers (book images → NASA substitutes → synthetic) | Most chapters need no data at all | Benchmark tiers: analytic/synthetic → cited published benchmarks → open datasets → manual download |
| `seaice/core/public_images.py`, substitute-image registry | Not needed | `fluidpy/core/refdata.py` (reference-data resolution + registry) and `reference/chNN/SOURCES.md` |
| Chapter ↔ code-folder mapping confirmed by reading the code | The book's own table of contents is the only structure | `tools/detect_chapters.py` + a mandatory human check of `book.yaml` |
| MATLAB semantics pitfalls (1-based indexing, column-major, `imfilter` = correlation) | No MATLAB | Numerical pitfalls of this domain: sign conventions, ν vs μ, reference scales, half-cell boundaries, cancellation, stiffness, degrees vs radians (`math-to-python`) |

## C. What transfers unchanged (and is worth keeping)

* **The seven-phase pipeline with gates**, and the rule that nothing downstream runs on a failed phase.
* **Subagents with written briefs**, so the orchestrator's context stays small and any session can resume from disk.
* **One chapter per session, `/clear` in between.**
* **Knowledge that compounds**: `knowledge/chNN.md` → `CUMULATIVE.md` → read first by the next chapter's analyst; and
  the habit of promoting a general lesson into the skill file itself, so the kit improves while it is used.
* **Publishing as part of "done"**: the chapter is not finished when the tests pass, it is finished when the notebook
  runs headlessly, the page is live and the public-repo check is clean.
* **The independent reviewer**, given fresh context and told to assume there is a bug.

---

## D. Lessons from setting up the fluidpy kit itself (12 September 2026)

| # | What happened | Fix now in the kit |
|---|---|---|
| D1 | A half-finished upgrade left `CLAUDE.md` describing agents, skills and tools that did not exist | every referenced file now exists; `tests/test_machinery.py` + `tools/pipeline_selftest.py` prove the machinery end to end |
| D2 | The explainer "passed" the no-overflow audit while clipped on phones: the `#app` mount had no height, so nothing could overflow | `#app { height: 100% }`; the audit also fails `app-does-not-fill-window` and `app-exceeds-window`; screenshots are always *looked at* |
| D3 | The canvas was sized with the border box → 2 px "overflow" pushed every view to maximum density | stage canvas uses `clientWidth/clientHeight` |
| D4 | A floating Contents button covered the explainer's Next button on the page | hidden while an explainer block is on screen (IntersectionObserver) |
| D5 | `import site` inside `tools/` silently imported the standard-library module | the site builder is `tools/build_site.py` |
| D6 | The venv was created from Anaconda's Python; its `python3` kernelspec could run notebooks in the wrong interpreter | publishing registers and uses the `fluidpy-venv` kernel |
| D7 | Shell heredocs on this Windows machine mangled `\n` inside Python/JS code | CLAUDE.md rule 10: write code files with Write/Edit |
| D8 | The chapter-opener check failed because the book sets "C H A P T E R" letter-spaced; the console crashed on the "ﬂ" ligature | whitespace-insensitive check; ligatures expanded in `chapters/*.txt`; UTF-8 console env in `.claude/settings.json` |
| D9 | `requirements-colab.txt` re-pinned numpy/scipy, which makes Colab upgrade them and demand a runtime restart | it lists only what Colab lacks (`pint`) |
| D10 | ipywidgets freeze on a static page; `to_jshtml` animations are large | plotly `slider_figure`/`animate_figure` for the page, `live` only as an extra; MP4 via imageio-ffmpeg |
