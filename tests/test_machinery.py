"""Fast tests of the project machinery (no browser, no notebook execution). The slow end-to-end checks are
``tools/shot.py templates/viz_example.html`` and ``tools/pipeline_selftest.py``."""
from __future__ import annotations

import json
import re
from pathlib import Path

import nbformat
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_book_yaml_chapter_map_is_consistent():
    book = yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))
    chs = [c for c in book["chapters"] if c.get("do_chapter", True)]
    assert len(chs) == 16
    prev_end = 0
    for c in chs:
        assert re.fullmatch(r"ch\d\d", c["id"]) and re.fullmatch(r"[a-z0-9_]+", c["slug"])
        assert c["pdf_start"] == prev_end + 1 or prev_end == 0, f"gap before {c['id']}"
        assert c["pdf_start"] - c["printed_start"] == book["book"]["printed_offset"]
        assert c["sections"] and c["pdf_end"] >= c["pdf_start"]
        prev_end = c["pdf_end"]


def test_progress_json_has_every_chapter_and_phase():
    prog = json.loads((ROOT / "progress.json").read_text(encoding="utf-8"))
    book = yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))
    assert prog["phases"] == ["analyze", "curate", "design", "implement", "verify", "review", "viz", "notebook",
                              "knowledge", "publish"]
    for c in book["chapters"]:
        row = prog["chapters"][c["id"]]
        assert all(p in row for p in prog["phases"])


def test_project_urls():
    from fluidpy.core import project

    assert project.pages_url("viz/ch07/x.html", ROOT) == "https://shammun.github.io/fluidpy/viz/ch07/x.html"
    assert project.colab_url("ch07_gravity_waves_colab", ROOT).endswith(
        "/github/shammun/fluidpy/blob/main/notebooks/ch07_gravity_waves_colab.ipynb")


def test_show_viz_srcdoc_and_placeholder(tmp_path):
    from fluidpy.core import embed

    (tmp_path / "viz" / "ch99").mkdir(parents=True)
    (tmp_path / "book.yaml").write_text((ROOT / "book.yaml").read_text(encoding="utf-8"), encoding="utf-8")
    src = (ROOT / "templates" / "viz_example.html").read_text(encoding="utf-8")
    (tmp_path / "viz" / "ch99" / "demo.html").write_text(src, encoding="utf-8")
    html = embed.viz_html("ch99", "demo", root=tmp_path)
    assert "srcdoc=" in html and embed.embedded_keys(html) == ["ch99/demo"]
    assert "Drift plus spin" in html                     # title from the viz:title meta tag
    missing = embed.viz_html("ch99", "nope", root=tmp_path)
    assert "fluidpy-viz-missing" in missing and embed.embedded_keys(missing) == ["ch99/nope"]
    with pytest.raises(ValueError):
        embed.viz_html("ch99", "../etc", root=tmp_path)


def test_viz_library_inlined_and_template_lints():
    import sys

    sys.path.insert(0, str(ROOT / "tools"))
    import viz_inline
    import viz_lint

    assert viz_inline.main(["--all", "--check"]) == 0
    msgs = viz_lint.lint_text(ROOT / "templates" / "viz_example.html",
                              (ROOT / "templates" / "viz_example.html").read_text(encoding="utf-8"))
    assert [m for m in msgs if not m.startswith("WARN")] == []


def test_nbkit_builds_a_valid_notebook(tmp_path, monkeypatch):
    import sys

    sys.path.insert(0, str(ROOT / "tools"))
    import nbkit

    monkeypatch.setattr(nbkit, "ROOT", tmp_path)
    (tmp_path / "book.yaml").write_text((ROOT / "book.yaml").read_text(encoding="utf-8"), encoding="utf-8")
    monkeypatch.setattr(nbkit, "_book", lambda: yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8")))
    nb = nbkit.ChapterNotebook("ch07")
    nb.title(big_idea="Waves.", roadmap=["one", "two"]).setup().section("7.2", "Linear waves")
    nb.code("x = 1  # one", explain="Sets x.").figure_notes("a", "b", "c")
    nb.explainer("dispersion_relation", heading="h", why="w", tries=["t"])
    with pytest.raises(ValueError):
        nb.explainer("dispersion_relation", heading="h", why="w", tries=["t"])
    nb.note("Short waves ride on long ones.", equation=r"c = \sqrt{gH}", ref="7.x").pointer("Covered in Ch. 13.")
    with pytest.raises(ValueError, match="7.3"):              # every book section must appear
        nb.save()
    for sec in nb.missing_sections():
        nb.section(sec, "…")
    path = nb.save()
    book = nbformat.read(path, as_version=4)
    assert any("\\sqrt{gH}" in c.source and "(7.x)" in c.source for c in book.cells)
    assert book.cells[0].source.startswith("# Chapter 7")
    assert not book.cells[0].source.splitlines()[2].startswith("    ")          # no accidental code-block indent
    assert nbkit.explainer_calls(path) == ["dispersion_relation"]
    assert any("setup" in c.metadata.get("tags", []) for c in book.cells)


def test_publish_transform_page_and_ipynb():
    import sys

    sys.path.insert(0, str(ROOT / "tools"))
    import publish_notebook as pub
    from fluidpy.core.embed import BEGIN, END

    nb = nbformat.v4.new_notebook()
    c = nbformat.v4.new_code_cell('show_viz("ch07", "demo")')
    c.outputs = [nbformat.v4.new_output("display_data", data={"text/html": BEGIN.format(key="ch07/demo") + "<iframe srcdoc='x'></iframe>" + END})]
    live = nbformat.v4.new_code_cell("live(...)", metadata={"tags": ["live-only"]})
    live.outputs = [nbformat.v4.new_output("display_data", data={"text/plain": "widget"})]
    nb.cells = [c, live]
    page, keys = pub.transform(nb, "page")
    assert keys == ["ch07/demo"]
    assert 'section class="fluidpy-viz fp-viz" id="viz-demo"' in page.cells[0].outputs[0]["data"]["text/html"]
    assert "fp-live-note" in page.cells[1].outputs[0]["data"]["text/html"]
    ipy, _ = pub.transform(nb, "ipynb")
    assert "Open it full-window" in ipy.cells[0].outputs[0]["data"]["text/html"]
    colab = pub.colab_variant(nb, "x")
    assert all(not cc.outputs for cc in colab.cells)


def test_check_public_rules(tmp_path):
    import sys

    sys.path.insert(0, str(ROOT / "tools"))
    import check_public as cp

    assert cp.FORBIDDEN.search("chapters/ch07.txt") and cp.FORBIDDEN.search("book.PDF")
    assert cp.FORBIDDEN.search("reports/viz/ch07/x/a.png") and cp.FORBIDDEN.search("tests/book_values_ch07.json")
    assert not cp.FORBIDDEN.search("chapters/_page_map.json") and not cp.FORBIDDEN.search("viz/ch07/a.html")
