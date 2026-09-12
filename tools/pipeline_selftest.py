"""End-to-end proof that the integration machinery works, without touching any chapter.

Builds a throw-away notebook in ``outputs/_selftest/`` that uses every kind of cell a chapter notebook uses — the
setup cell, a matplotlib figure, an animation (video and frames players), a plotly slider figure, a live-only
ipywidgets cell and an embedded interactive explainer (``templates/viz_example.html``) — then:

  1. executes it headlessly with nbconvert (the same ExecutePreprocessor the publish tool uses);
  2. checks the explainer output is an ``<iframe srcdoc>`` block with the fluidpy-viz markers;
  3. renders the published page with the real ``tools/publish_notebook.py`` transforms and chrome;
  4. audits that page in a headless browser (``tools/shot.py --page``): the explainer block spans the full window
     width and height, the explainer inside passes its own no-overflow audit, the page never scrolls sideways;
  5. checks the plotly figure and both animation players rendered in the page, and the live cell became a note;
  6. writes the Colab twin and checks it has no outputs and the setup cell.

Usage::  .venv/Scripts/python.exe tools/pipeline_selftest.py            (exit 0 = everything works)
Screenshots: reports/viz/_pages/selftest_page/  ·  page: outputs/_selftest/notebooks/selftest_page.html
"""
from __future__ import annotations

import shutil
import sys
import time
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

ST = ROOT / "outputs" / "_selftest"


def build_notebook() -> Path:
    import nbkit

    (ST / "viz" / "ch99").mkdir(parents=True, exist_ok=True)
    (ST / "notebooks").mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "templates" / "viz_example.html", ST / "viz" / "ch99" / "viz_example.html")
    shutil.copytree(ROOT / "assets", ST / "assets", dirs_exist_ok=True)
    shutil.copy(ROOT / "book.yaml", ST / "book.yaml")

    nb = nbf.v4.new_notebook()
    C = nb.cells
    C.append(nbf.v4.new_markdown_cell("# Pipeline self-test\n\nA notebook with every kind of cell.\n\n## 1 Setup"))
    setup = nbkit.SETUP_CODE.format(repo="shammun/fluidpy")
    C.append(nbf.v4.new_code_cell(setup, metadata={"tags": ["setup"]}))
    C.append(nbf.v4.new_markdown_cell("## 2 A matplotlib figure"))
    C.append(nbf.v4.new_code_cell(
        "import numpy as np, matplotlib.pyplot as plt\n"
        "x = np.linspace(0, 2*np.pi, 200)\nfig, ax = plt.subplots()\nax.plot(x, np.sin(x), label='sin')\n"
        "ax.set_xlabel('x [rad]'); ax.set_ylabel('y [-]'); ax.legend(); plt.show()"))
    C.append(nbf.v4.new_markdown_cell("## 3 Animations"))
    C.append(nbf.v4.new_code_cell(
        "from fluidpy.core.anim import animate, show_animation\n"
        "fig, ax = plt.subplots(figsize=(5, 2.5)); (ln,) = ax.plot(x, np.sin(x))\n"
        "def upd(i):\n    ln.set_ydata(np.sin(x - 0.2*i)); return (ln,)\n"
        "show_animation(animate(upd, frames=24, fig=fig), player='video')", metadata={"tags": ["animation"]}))
    C.append(nbf.v4.new_code_cell(
        "fig, ax = plt.subplots(figsize=(5, 2.5)); (ln,) = ax.plot(x, np.cos(x))\n"
        "def upd(i):\n    ln.set_ydata(np.cos(x - 0.2*i)); return (ln,)\n"
        "show_animation(animate(upd, frames=12, fig=fig), player='frames')", metadata={"tags": ["animation"]}))
    C.append(nbf.v4.new_markdown_cell("## 4 A plotly slider figure"))
    C.append(nbf.v4.new_code_cell(
        "from fluidpy.core.interact import slider_figure\n"
        "k = np.linspace(0.05, 3, 200)\n"
        "fig = slider_figure(lambda H: {'c(k)': (k, np.sqrt(9.81/k*np.tanh(k*H)))}, 'H', np.linspace(0.5, 20, 12),\n"
        "                    unit='m', xlabel='k [1/m]', ylabel='c [m/s]', title='slider test')\nfig.show()",
        metadata={"tags": ["plotly"]}))
    C.append(nbf.v4.new_markdown_cell("## 5 A live-only widget"))
    C.append(nbf.v4.new_code_cell(
        "from fluidpy.core.interact import live\n"
        "live(lambda a: print('a =', a), a=(0.0, 1.0, 0.1))", metadata={"tags": ["live-only"]}))
    C.append(nbf.v4.new_markdown_cell("## 6 An interactive explainer\n\nText before the explainer."))
    C.append(nbf.v4.new_code_cell(
        f'show_viz("ch99", "viz_example", root=r"{ST.as_posix()}")', metadata={"tags": ["explainer"]}))
    C.append(nbf.v4.new_markdown_cell("Text after the explainer.\n\n## 7 The end"))
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    path = ST / "notebooks" / "selftest_page.ipynb"
    nbf.write(nb, path)
    return path


def main() -> int:
    import publish_notebook as pub
    import shot

    fails: list[str] = []
    t0 = time.perf_counter()
    path = build_notebook()
    print("1. built", path.relative_to(ROOT).as_posix())

    executed, secs = pub.execute(path)
    probs = pub.execution_problems(executed)
    print(f"2. executed in {secs:.0f} s; problems: {probs or 'none'}")
    fails += probs

    htmls = [pub._html_of(o) or "" for c in executed.cells if c.cell_type == "code" for o in c.outputs]
    viz_out = [h for h in htmls if "fluidpy-viz:BEGIN ch99/viz_example" in h]
    if not viz_out or "srcdoc=" not in viz_out[0]:
        fails.append("explainer output is not an iframe srcdoc block with markers")
    else:
        print("3. explainer output: srcdoc iframe with markers, %.0f kB" % (len(viz_out[0]) / 1e3))
    if not any("<video" in h for h in htmls):
        fails.append("no <video> animation output (is imageio-ffmpeg installed?)")
    if not any("animation" in h and "anim-controls" in h or "_anim_" in h for h in htmls):
        fails.append("no frames (jshtml) animation player output")
    if not any("plotly" in h.lower() for h in htmls):
        fails.append("no plotly html output")

    page_nb, keys = pub.transform(executed, "page")
    row = {"id": "ch99", "number": 99, "title": "Self-test", "slug": "selftest_page", "blurb": "pipeline self-test"}
    page = pub.render_page(page_nb, row=row, name="selftest_page", prev_row=None, next_row=None, title="Self-test")
    page = page.replace('href="../viz/ch99/', 'href="../viz/ch99/')
    out = ST / "notebooks" / "selftest_page.html"
    out.write_text(page, encoding="utf-8")
    print(f"4. page written ({out.stat().st_size / 1e6:.2f} MB), explainer keys: {keys}")
    if "fp-live-note" not in page:
        fails.append("live-only output was not replaced by the note")

    colab = pub.colab_variant(executed, "selftest_page")
    if any(c.get("outputs") for c in colab.cells if c.cell_type == "code"):
        fails.append("colab variant still has outputs")
    if not any("setup" in c.metadata.get("tags", []) for c in colab.cells):
        fails.append("colab variant lost the setup cell")

    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        browser = shot.launch(pw, "auto")
        rep = shot.audit_page(browser, out, ROOT / "reports" / "viz" / "_pages" / "selftest_page")
        ctx = browser.new_context(viewport={"width": 1366, "height": 768})
        page_ = ctx.new_page()
        page_.goto(out.resolve().as_uri(), wait_until="load")
        page_.wait_for_timeout(2500)
        checks = page_.evaluate("""() => ({
            plotly: document.querySelectorAll('.js-plotly-plot .main-svg').length,
            video: document.querySelectorAll('video').length,
            frames: document.querySelectorAll('.animation img, img[id^="_anim_img"]').length,
            toc: !!document.querySelector('.fp-toc-btn'),
            moved: !!document.querySelector('.ce-container > section.fp-viz, main > section.fp-viz, .jp-Notebook > section.fp-viz')
        })""")
        ctx.close()
        browser.close()
    print("5. page audit:", rep["verdict"], rep["failures"] or "", "| in-page checks:", checks)
    fails += rep["failures"]
    if not checks["plotly"]:
        fails.append("plotly figure did not render on the page (renderer / require.js problem)")
    if not checks["video"]:
        fails.append("video animation missing on the page")
    if not checks["frames"]:
        fails.append("frames animation player missing on the page")
    if not checks["moved"]:
        fails.append("explainer block was not moved out of its notebook cell")

    print(f"total {time.perf_counter() - t0:.0f} s")
    if fails:
        print("PIPELINE SELF-TEST FAILED:")
        for f in fails:
            print("  -", f)
        return 1
    print("PIPELINE SELF-TEST PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
