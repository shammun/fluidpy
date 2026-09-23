"""Mechanical audit of interactive explainers and published chapter pages in a real (headless) browser.

EXPLAINER mode — for each ``viz/chNN/<slug>.html`` (or any explainer path):
  * opens it at every viewport in ``SIZES`` (phone portrait/landscape, tablet, laptop, desktop, notebook iframe);
  * waits for ``window.VIZ.ready`` (KaTeX loaded or fallback), then visits EVERY tab and EVERY walkthrough step;
  * fails on: page or panel overflow (the no-scroll rule), console ``VIZ-OVERFLOW`` / ``VIZ-ERROR`` / JS errors,
    visible text below 12 px, tap targets below 24 px on phones, a walkthrough with < 4 or > 8 steps,
    no equations tab, missing ``viz:*`` meta tags; visits every page of every derivation too;
  * pages through EVERY page of every pager in the active panel (walkthrough card, Explore, Explain, Derivation,
    Equations, Code, Check) of every view and fails on any item whose content is clipped by its page (below the
    page, or a slice cut through a line of text / a formula / a code line) or wider than it, or on a page that
    squeezes a visible view below 60 px — reported as
    ``[size] view: pager page k/n: item 'text' clipped N px``; a screenshot ``CLIP__<size>__<view>__pgK.png`` is
    saved for the first few;
  * quality floor: Explain tab with >= 3 numbered sections + an interpretation, Code tab, >= 2 depth features,
    >= 3 questions; every derivation listed in ``viz:derivations`` exists and each step has tex, did, why, plain;
  * runs ``selftest()``: rows with ``expect`` are checked in JS; rows with ``py`` are evaluated here in Python with the
    chapter's fluidpy module(s) importable as ``chNN`` (e.g. ``py: "ch07.phase_speed(0.5, 10.0)"``) plus ``np``,
    ``math``, ``fluidpy`` — this is the JS ↔ Python parity evidence;
  * writes screenshots + ``audit.json`` to ``reports/viz/<chapter>/<slug>/`` (git-ignored) and prints a verdict.

PAGE mode (``--page notebooks/chNN_<slug>.html``) — the published chapter page:
  * every ``.fluidpy-viz`` block is full-bleed (its width = the viewport width) and its iframe fills the window height
    (minus the thin toolbar); the explainer inside passes its own ``VIZ.audit()``; the page never scrolls sideways;
  * screenshots each block scrolled into view at 1366×768 and 390×844.

Usage (repo root)::

    .venv/Scripts/python.exe tools/shot.py viz/ch07/dispersion_relation.html
    .venv/Scripts/python.exe tools/shot.py --chapter ch07                 # all explainers of a chapter
    .venv/Scripts/python.exe tools/shot.py --chapter ch07 --quick         # 4 sizes (360x640, 390x844, notebook, desktop)
    .venv/Scripts/python.exe tools/shot.py --page notebooks/ch07_gravity_waves.html
    .venv/Scripts/python.exe tools/shot.py templates/viz_example.html

Browser: the installed Microsoft Edge or Google Chrome is used through Playwright (no browser download needed);
``--browser chromium`` uses Playwright's own build (``python -m playwright install chromium``).
Exit code 0 = PASS, 1 = FAIL.
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# name: (width, height, is_phone)
SIZES: dict[str, tuple[int, int, bool]] = {
    "phone": (360, 640, True),
    "phone-tall": (390, 844, True),
    "phone-land": (844, 390, True),
    "tablet": (768, 1024, False),
    "laptop": (1280, 720, False),
    "notebook": (1000, 700, False),
    "desktop": (1366, 768, False),
    "full-hd": (1920, 1080, False),
}
QUICK = ("phone", "phone-tall", "notebook", "desktop")   # 360x640 is where pages overflow first
SHOT_STEPS_AT = ("phone-tall", "desktop")   # every step is screenshotted at these sizes
REQUIRED_META = ("chapter", "slug", "title", "summary", "concept", "sections", "equations", "fluidpy", "derivations")

JS_TEXT_AUDIT = r"""
() => {
  const small = [], tiny = [];
  const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  let n;
  while ((n = walker.nextNode())) {
    const t = n.textContent.trim(); if (!t) continue;
    const el = n.parentElement; if (!el || el.closest('.katex, .viz-help:not(.open), script, style')) continue;
    const r = el.getBoundingClientRect(); if (r.width === 0 || r.height === 0) continue;
    const cs = getComputedStyle(el); if (cs.visibility === 'hidden' || cs.display === 'none') continue;
    let scale = 1, a = el; while (a && a !== document.body) { const tr = getComputedStyle(a).transform; if (tr && tr.startsWith('matrix(')) scale *= parseFloat(tr.slice(7)); a = a.parentElement; }
    const px = parseFloat(cs.fontSize) * scale;
    if (px < 11.5) small.push({ text: t.slice(0, 40), px: Math.round(px * 10) / 10, cls: el.className || el.tagName });
  }
  const phoneTargets = [];
  document.querySelectorAll('button, input[type=range], select, .viz-chip-toggle').forEach(b => {
    const r = b.getBoundingClientRect(); if (r.width === 0) return;
    if (r.height < 23.5 || r.width < 23.5) phoneTargets.push({ cls: b.className || b.tagName, w: Math.round(r.width), h: Math.round(r.height) });
  });
  return { small: small.slice(0, 12), nSmall: small.length, targets: phoneTargets.slice(0, 12), nTargets: phoneTargets.length };
}
"""


JS_PAGER_AUDIT = r"""
(shoot) => {
  // Page through EVERY page of every pager in the active panel and measure clipping.
  // shoot = [cardIndex, page] → just show that page (for a screenshot) and return.
  const ATOMS = '.katex-display, .katex, tr, button, input, select, canvas, svg, img';
  const cards = [...document.querySelectorAll('.viz-panel.is-active .viz-paged')].filter(c => c.offsetParent);
  const nav = (card) => {
    const P = card._vizPager;
    if (P && P.pageCount) return { n: P.pageCount(), cur: P.page(), show: k => P.show(k) };
    const prev = card.querySelector('.viz-pager-btn[aria-label="Previous page"]'), next = card.querySelector('.viz-pager-btn[aria-label="Next page"]');
    const lab = card.querySelector('.viz-pager-label'), m = lab && lab.textContent.match(/(\d+)\s*\/\s*(\d+)/);
    const n = m && !card.querySelector('.viz-pager').hidden ? +m[2] : 1, cur = m ? +m[1] - 1 : 0;
    return { n, cur, show: k => { let g = 0; while (prev && !prev.disabled && g++ < 99) prev.click(); for (let i = 0; i < k; i++) next.click(); } };
  };
  if (shoot) { const c = cards[shoot[0]]; if (c) nav(c).show(shoot[1]); return null; }
  const vis = el => { const cs = getComputedStyle(el); return cs.display !== 'none' && cs.visibility !== 'hidden'; };
  // a multi-line formula (aligned, cases) may continue on the next page between its rows: its rows are the atoms
  const katexBoxes = (el, boxes) => [...el.children].forEach(c => {
    if (c.classList.contains('strut') || c.classList.contains('pstrut')) return;
    if (c.classList.contains('mtable'))
      c.querySelectorAll(':scope > span > .vlist-t > .vlist-r > .vlist > span').forEach(row => {
        let t = Infinity, b = -Infinity, l = Infinity, r = -Infinity;
        [...row.children].forEach(x => { if (x.classList.contains('pstrut')) return; const q = x.getBoundingClientRect(); if (q.height) { t = Math.min(t, q.top); b = Math.max(b, q.bottom); l = Math.min(l, q.left); r = Math.max(r, q.right); } });
        if (b > t) boxes.push({ top: t, bottom: b, left: l, right: r });
      }), c.querySelectorAll('.vertical-separator').forEach(x => { const q = x.getBoundingClientRect(); if (q.height) boxes.push(q); });
    else if (c.querySelector('.mtable')) katexBoxes(c, boxes);
    else { const q = c.getBoundingClientRect(); if (q.height) boxes.push(q); }
  });
  const out = [];
  cards.forEach((card, ci) => {
    const N = nav(card), body = card.querySelector('.viz-pages');
    if (!body) return;
    const pages = [];
    for (let k = 0; k < N.n; k++) {
      N.show(k);
      const br = body.getBoundingClientRect(), issues = [];
      [...body.children].forEach(it => {
        if (it.classList.contains('viz-pager-cue') || !vis(it)) return;
        const r = it.getBoundingClientRect(); if (!r.height && !r.width) return;
        const own = getComputedStyle(it).overflow === 'hidden';          // a slice: the item clips itself on purpose
        const clip = own ? Math.min(br.bottom, r.bottom) : br.bottom;
        const txt = (it.textContent || '').trim().replace(/\s+/g, ' ').slice(0, 60);
        const boxes = [];
        const walk = el => {
          for (let c = el.firstChild; c; c = c.nextSibling) {
            if (c.nodeType === 3) {
              if (!c.textContent.trim()) continue;
              const rg = document.createRange(), pre = /pre/.test(getComputedStyle(el).whiteSpace); rg.selectNodeContents(c);
              [...rg.getClientRects()].forEach(q => q.height && boxes.push(pre ? { top: q.top, bottom: q.bottom, left: q.left, right: Math.min(q.right, el.getBoundingClientRect().right) } : q));
            }
            else if (c.nodeType === 1 && vis(c)) { const q = c.getBoundingClientRect(); if (!q.height && !q.width) continue;
              if (c.matches('.katex-display, .katex') && c.querySelector('.mtable')) katexBoxes(c.querySelector('.katex-html') || c, boxes);
              else if (c.matches(ATOMS)) boxes.push(q); else walk(c); }
          }
        };
        walk(it);
        let below = 0, cut = 0, wide = 0;
        boxes.forEach(q => {
          if (own && q.top >= r.bottom - 1) return;                           // the next slice of a sliced item
          if (q.bottom > clip + 1) { if (q.top < clip - 1) cut = Math.max(cut, q.bottom - clip); else below = Math.max(below, q.bottom - clip); }
          if (q.right > br.right + 1) wide = Math.max(wide, q.right - br.right);
        });
        if (r.bottom > br.bottom + 1) below = Math.max(below, r.bottom - br.bottom);
        it.querySelectorAll('.viz-eq-too-wide').forEach(e => { if (e.offsetParent && e.firstChild) wide = Math.max(wide, e.firstChild.getBoundingClientRect().right - e.getBoundingClientRect().right, 1.5); });
        if (below > 1 || cut > 1 || wide > 1)
          issues.push({ item: (it.className || it.tagName).toString().split(' ')[0], text: txt, clipped: Math.round(Math.max(below, cut)), cut: cut > 1, wide: Math.round(wide) });
      });
      // the page shown may move the picture/text split: every visible view must stay >= 60 px on every page
      document.querySelectorAll('.viz-view').forEach(v => {
        const hh = v.querySelector('.viz-view-holder'); if (!v.offsetParent || !hh) return;
        if (hh.clientHeight < 59.5) issues.push({ item: 'view-too-small', text: v.getAttribute('data-view') || v.className, clipped: 0, cut: false, wide: 0, small: Math.round(hh.clientHeight) });
      });
      pages.push(issues);
    }
    N.show(N.cur);
    out.push({ card: ci, cls: card.className.replace('viz-card grow viz-paged', '').trim() || 'pager', n: N.n, pages });
  });
  return out;
}
"""


def launch(pw, browser: str):
    tries = [browser] if browser != "auto" else ["msedge", "chrome", "chromium"]
    last = None
    for b in tries:
        try:
            if b == "chromium":
                return pw.chromium.launch(headless=True)
            return pw.chromium.launch(headless=True, channel=b)
        except Exception as exc:  # noqa: BLE001
            last = exc
    raise SystemExit(f"could not start a browser ({tries}): {last}\n"
                     f"Install one with: .venv/Scripts/python.exe -m playwright install chromium")


def chapter_namespace(chapter: str) -> dict:
    """Names available to selftest ``py`` expressions: np, math, fluidpy, and chNN → the chapter module(s)."""
    import numpy as np

    ns: dict = {"np": np, "math": math}
    try:
        import fluidpy  # noqa: F401
        ns["fluidpy"] = importlib.import_module("fluidpy")
    except Exception:  # noqa: BLE001
        pass
    mods = sorted((ROOT / "fluidpy").glob(f"{chapter}_*.py")) if re.fullmatch(r"ch\d\d", chapter or "") else []
    if len(mods) == 1:
        ns[chapter] = importlib.import_module(f"fluidpy.{mods[0].stem}")
    elif len(mods) > 1:
        class _Multi:  # ch07.<function> searches every chapter module
            def __getattr__(self, name):
                for m in mods:
                    mod = importlib.import_module(f"fluidpy.{m.stem}")
                    if hasattr(mod, name):
                        return getattr(mod, name)
                raise AttributeError(name)
        ns[chapter] = _Multi()
    core = ROOT / "fluidpy" / "core"
    if core.is_dir():
        ns["core"] = importlib.import_module("fluidpy.core")
    return ns


def parity(rows: list[dict], chapter: str) -> tuple[list[dict], list[str]]:
    ns = chapter_namespace(chapter)
    out, fails = [], []
    for r in rows:
        row = dict(r)
        if r.get("py"):
            try:
                val = eval(r["py"], {"__builtins__": {}}, ns)  # noqa: S307 - our own repository's expressions
                val = float(val)
                row["py_value"] = val
                tol = r.get("atol", 0) + r.get("rtol", 1e-6) * abs(val)
                row["ok"] = r.get("js") is not None and abs(float(r["js"]) - val) <= tol
            except Exception as exc:  # noqa: BLE001
                row["ok"] = False
                row["error"] = f"{type(exc).__name__}: {exc}"
        if row.get("ok") is False:
            fails.append(f"selftest '{r.get('name')}': js={r.get('js')} py={row.get('py_value', row.get('expect'))} "
                         f"{row.get('error', '')}".strip())
        out.append(row)
    return out, fails


def audit_explainer(browser, path: Path, sizes: list[str], out_dir: Path, shots: bool) -> dict:
    url = path.resolve().as_uri()
    rep: dict = {"file": path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path), "sizes": {}, "failures": []}
    fails: list[str] = rep["failures"]
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = {}
    for name in sizes:
        w, h, phone = SIZES[name]
        ctx = browser.new_context(viewport={"width": w, "height": h}, device_scale_factor=2 if phone else 1,
                                  is_mobile=False, has_touch=phone)
        page = ctx.new_page()
        logs: list[str] = []
        page.on("console", lambda m, logs=logs: logs.append(f"{m.type}: {m.text}"))
        page.on("pageerror", lambda e, logs=logs: logs.append(f"pageerror: {e}"))
        t0 = time.perf_counter()
        page.goto(url, wait_until="load")
        try:
            page.wait_for_function("() => window.VIZ && window.VIZ.ready", timeout=15000)
            page.evaluate("() => window.VIZ.ready")
        except Exception as exc:  # noqa: BLE001
            fails.append(f"[{name}] window.VIZ never became ready: {exc}")
            ctx.close()
            continue
        info = page.evaluate("() => ({tabs: VIZ.tabs, steps: VIZ.steps, meta: VIZ.meta(), katex: !!window.katex, features: VIZ.features || {}, nCheck: (VIZ.app.cfg.check || []).length, derivations: VIZ.derivations || [], explain: VIZ.explainStats ? VIZ.explainStats() : null, derSteps: (VIZ.app.cfg.derivations || []).map(d => ({id: d.id, goal: !!d.goal, result: !!(d.result && (d.result.tex || d.result.plain)), steps: d.steps.map(x => ({tex: !!x.tex, did: !!x.did, why: (x.why || '').split(' ').filter(Boolean).length, plain: !!x.plain}))}))})")
        rep["features"], rep["n_check"] = info["features"], info["nCheck"]
        rep["derivations"], rep["explain_stats"], rep["der_steps"] = info["derivations"], info["explain"], info["derSteps"]
        meta = info["meta"]
        srep = {"viewport": [w, h], "load_s": round(time.perf_counter() - t0, 2), "katex": info["katex"], "views": []}
        n_clip_shots = [0]
        for old in out_dir.glob(f"CLIP__{name}__*.png"):
            old.unlink()
        for tab in info["tabs"]:
            views = [(tab, None)]
            if tab == "tour":
                views = [(tab, i) for i in range(info["steps"])]
            if tab == "derive":                    # every page of every derivation (goal, steps, result)
                views = [(tab, (di, si)) for di, d in enumerate(info["derivations"]) for si in range(d["pages"])]
            for tab_, step in views:
                if step is None:
                    a = page.evaluate(f"() => VIZ.setTab({json.dumps(tab_)})")
                elif tab_ == "derive":
                    a = page.evaluate(f"() => VIZ.goDerive({step[0]}, {step[1]})")
                else:
                    a = page.evaluate(f"() => VIZ.goStep({step})")
                page.wait_for_timeout(120)
                a = page.evaluate("() => VIZ.audit()")
                label = f"{tab_}" + ("" if step is None else f"-d{step[0] + 1}p{step[1]}" if tab_ == "derive" else f"-step{step + 1}")
                text = page.evaluate(JS_TEXT_AUDIT)
                view = {"view": label, "layout": a["layout"], "dense": a["dense"], "overflow": a["overflow"],
                        "small_text": text["nSmall"]}
                if a["overflow"]:
                    fails.append(f"[{name}] {label}: overflow {a['overflow'][:3]}")
                if text["nSmall"]:
                    fails.append(f"[{name}] {label}: {text['nSmall']} text item(s) below 12px, e.g. {text['small'][:3]}")
                if phone and text["nTargets"] and step in (None, 0, (0, 1)):
                    view["small_targets"] = text["targets"]
                    fails.append(f"[{name}] {label}: {text['nTargets']} tap target(s) below 24px, e.g. {text['targets'][:3]}")
                if shots and (step is None or step == 0 or name in SHOT_STEPS_AT):
                    page.screenshot(path=str(out_dir / f"{name}__{label}.png"))
                # every page of every pager in this panel: nothing may be clipped (the first page is not enough)
                pg = page.evaluate(JS_PAGER_AUDIT, None)
                view["pages"] = sum(p["n"] for p in pg)
                clipped = []
                for p in pg:
                    for k, issues in enumerate(p["pages"]):
                        for it in issues:
                            clipped.append({"pager": p["cls"], "page": k + 1, "of": p["n"], **it})
                            what = (f"clipped {it['clipped']}px" + (" (cut through a line)" if it["cut"] else "")) if it["clipped"] else ""
                            what += (", " if what else "") + (f"{it['wide']}px too wide" if it["wide"] else "")
                            if it["item"] == "view-too-small":
                                what = f"view only {it.get('small')}px high on this page"
                                if a["overflow"] and any(str(o.get("el", "")).startswith("view-too-small") for o in a["overflow"]):
                                    continue            # already reported for page 1 by VIZ.audit()
                            fails.append(f"[{name}] {label}: {p['cls']} page {k + 1}/{p['n']}: {it['item']} "
                                         f"'{it['text'][:40]}' {what}")
                            if shots and n_clip_shots[0] < 6:
                                n_clip_shots[0] += 1
                                page.evaluate(JS_PAGER_AUDIT, [p["card"], k])
                                page.screenshot(path=str(out_dir / f"CLIP__{name}__{label}__pg{k + 1}.png"))
                                page.evaluate(JS_PAGER_AUDIT, [p["card"], 0])
                if clipped:
                    view["clipped"] = clipped
                srep["views"].append(view)
        errs = [l for l in logs if ("VIZ-ERROR" in l or l.startswith("pageerror") or l.startswith("error"))
                and "katex" not in l.lower() and "net::" not in l]
        if errs:
            fails.append(f"[{name}] console errors: {errs[:3]}")
        srep["console"] = logs[-20:]
        rep["sizes"][name] = srep
        if name == sizes[-1] or name == "desktop":
            rep["selftest_js"] = page.evaluate("() => VIZ.selftest()")
            rep["tabs"], rep["steps"] = info["tabs"], info["steps"]
        ctx.close()
    # structural rules
    missing = [k for k in REQUIRED_META if not meta.get(k) or meta.get(k, "").startswith("__")]
    if missing:
        fails.append(f"missing/placeholder viz:* meta tags: {missing}")
    steps = rep.get("steps", 0)
    if not (4 <= steps <= 8):
        fails.append(f"walkthrough has {steps} steps (need 4-8)")
    if "equations" not in rep.get("tabs", []):
        fails.append("no Equations tab (every explainer shows and explains its equations)")
    # quality floor modelled on the reference explainers (skill interactive-viz §4)
    feats = rep.get("features") or {}
    if not (feats.get("explain") or feats.get("calc")):
        fails.append("no Explain tab (explain): compute every displayed number with the reader's settings, then interpret")
    else:
        ex = rep.get("explain_stats") or {}
        if ex.get("sections", 0) < 3:
            fails.append(f"Explain tab has {ex.get('sections', 0)} numbered sections (need >= 3: Viz.work.step(n, title))")
        if ex.get("interpret", 0) < 1:
            fails.append("Explain tab has no interpretation (Viz.work.interpret: what the current setting means)")
    # derivations: meta viz:derivations lists the ids the storyboard asks for ("none" = no derivation-heavy idea here)
    wanted = [x for x in re.split(r"[\s,;]+", meta.get("derivations", "")) if x and x.lower() != "none"]
    have = {d.get("id") for d in rep.get("derivations") or []}
    for w in wanted:
        if w not in have:
            fails.append(f"derivation '{w}' is listed in viz:derivations but not in derivations: [...]")
    for d in rep.get("der_steps") or []:
        if not d["goal"] or not d["result"]:
            fails.append(f"derivation '{d['id']}': needs a goal and a result")
        if len(d["steps"]) < 2:
            fails.append(f"derivation '{d['id']}': {len(d['steps'])} step(s); show every move (>= 2 steps)")
        for k, st in enumerate(d["steps"], 1):
            miss = [f for f in ("tex", "did", "plain") if not st[f]] + (["why (>= 6 words)"] if st["why"] < 6 else [])
            if miss:
                fails.append(f"derivation '{d['id']}' step {k}: missing {', '.join(miss)}")
    if not feats.get("code"):
        fails.append("no Code tab: show the Python behind the picture, synced to the walkthrough")
    depth = [k for k in ("transport", "presets", "status", "terms", "inspect", "notes", "modes") if feats.get(k)]
    if (feats.get("views") or 0) >= 2:
        depth.append("views")
    if len(depth) < 2:
        fails.append(f"only {len(depth)} depth feature(s) {depth}: use at least 2 of linked views, transport, presets, "
                     "status, terms, inspector, notes, modes")
    rep["depth_features"] = depth
    if rep.get("n_check", 0) < 3:
        fails.append(f"{rep.get('n_check', 0)} check-yourself questions (need at least 3)")
    rows = rep.get("selftest_js") or []
    if not rows:
        fails.append("selftest() returned no rows (parity with fluidpy is required)")
    chapter = meta.get("chapter", path.parent.name)
    rows, pfails = parity(rows, chapter)
    rep["selftest"] = rows
    fails.extend(pfails)
    if rows and chapter.startswith("ch") and not any(r.get("py") for r in rows):
        fails.append("selftest() has no py parity row against a fluidpy function")
    rep["verdict"] = "PASS" if not fails else "FAIL"
    (out_dir / "audit.json").write_text(json.dumps(rep, indent=2, default=str), encoding="utf-8")
    return rep


JS_PAGE_AUDIT = r"""
() => {
  const W = window.innerWidth, H = window.innerHeight, out = [];
  document.querySelectorAll('.fluidpy-viz').forEach(b => {
    const r = b.getBoundingClientRect(), f = b.querySelector('iframe'), bar = b.querySelector('.fluidpy-viz-bar');
    const fr = f ? f.getBoundingClientRect() : {width: 0, height: 0};
    out.push({ key: b.getAttribute('data-viz'), blockW: Math.round(r.width), iframeW: Math.round(fr.width), iframeH: Math.round(fr.height),
               barH: bar ? Math.round(bar.getBoundingClientRect().height) : 0, src: f ? (f.getAttribute('src') || 'srcdoc') : null });
  });
  return { W, H, hscroll: document.documentElement.scrollWidth > W + 1, blocks: out };
}
"""


def audit_page(browser, path: Path, out_dir: Path) -> dict:
    url = path.resolve().as_uri()
    rep: dict = {"page": path.relative_to(ROOT).as_posix(), "sizes": {}, "failures": []}
    fails = rep["failures"]
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in ("desktop", "phone-tall"):
        w, h, phone = SIZES[name]
        ctx = browser.new_context(viewport={"width": w, "height": h}, device_scale_factor=1)
        page = ctx.new_page()
        logs: list[str] = []
        page.on("console", lambda m, logs=logs: logs.append(f"{m.type}: {m.text}"))
        page.goto(url, wait_until="load")
        page.wait_for_timeout(800)
        a = page.evaluate(JS_PAGE_AUDIT)
        if a["hscroll"]:
            fails.append(f"[{name}] the page scrolls horizontally")
        if not a["blocks"]:
            fails.append(f"[{name}] no .fluidpy-viz explainer blocks found on the page")
        for i, b in enumerate(a["blocks"]):
            if abs(b["blockW"] - w) > 2:
                fails.append(f"[{name}] {b['key']}: block width {b['blockW']} != viewport {w} (not full-bleed)")
            if abs(b["iframeH"] + b["barH"] - h) > 4:
                fails.append(f"[{name}] {b['key']}: iframe {b['iframeH']}px + bar {b['barH']}px != window {h}px")
            loc = page.locator(".fluidpy-viz").nth(i)
            loc.scroll_into_view_if_needed()
            page.evaluate(f"() => document.querySelectorAll('.fluidpy-viz')[{i}].scrollIntoView({{block: 'start'}})")
            page.wait_for_timeout(900)
            try:
                inner = page.locator(".fluidpy-viz iframe").nth(i).element_handle().content_frame()
                inner.wait_for_function("() => window.VIZ && window.VIZ.ready", timeout=15000)
                inner.evaluate("() => window.VIZ.ready")
                ia = inner.evaluate("() => VIZ.audit()")
                b["inner"] = ia
                if ia["overflow"]:
                    fails.append(f"[{name}] {b['key']}: explainer overflows inside the page: {ia['overflow'][:2]}")
            except Exception as exc:  # noqa: BLE001
                b["inner_error"] = str(exc)
                fails.append(f"[{name}] {b['key']}: explainer did not load inside the page: {exc}")
            page.screenshot(path=str(out_dir / f"page_{name}__{i + 1}_{b['key'].replace('/', '_')}.png"))
        rep["sizes"][name] = a
        errs = [l for l in logs if l.startswith("error") and "net::" not in l]
        if errs:
            rep["sizes"][name]["console_errors"] = errs[:5]
        ctx.close()
    rep["verdict"] = "PASS" if not fails else "FAIL"
    (out_dir / "page_audit.json").write_text(json.dumps(rep, indent=2, default=str), encoding="utf-8")
    return rep


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("paths", nargs="*", help="explainer HTML files")
    ap.add_argument("--chapter", action="append", help="audit every viz/<chapter>/*.html (repeatable)")
    ap.add_argument("--page", help="audit a published chapter page (notebooks/chNN_<slug>.html)")
    ap.add_argument("--quick", action="store_true", help=f"only {QUICK}")
    ap.add_argument("--sizes", help="comma-separated subset of " + ",".join(SIZES))
    ap.add_argument("--no-shots", action="store_true")
    ap.add_argument("--browser", default="auto", choices=["auto", "msedge", "chrome", "chromium"])
    a = ap.parse_args(argv)

    from playwright.sync_api import sync_playwright

    files = [Path(p) if Path(p).is_absolute() else ROOT / p for p in a.paths]
    for ch in a.chapter or []:
        files += sorted(p for p in (ROOT / "viz" / ch).glob("*.html") if p.name != "index.html")
    if not files and not a.page:
        ap.error("give explainer paths, --chapter chNN or --page notebooks/<page>.html")
    sizes = a.sizes.split(",") if a.sizes else (list(QUICK) if a.quick else list(SIZES))
    bad = [s for s in sizes if s not in SIZES]
    if bad:
        ap.error(f"unknown sizes {bad}")

    ok = True
    with sync_playwright() as pw:
        browser = launch(pw, a.browser)
        for f in files:
            chapter = f.parent.name
            out_dir = ROOT / "reports" / "viz" / chapter / f.stem
            rep = audit_explainer(browser, f, sizes, out_dir, shots=not a.no_shots)
            ok &= rep["verdict"] == "PASS"
            n_views = sum(len(s["views"]) for s in rep["sizes"].values())
            print(f"{rep['verdict']}  {rep['file']}  ({len(rep['sizes'])} sizes, {n_views} views, "
                  f"{rep.get('steps', 0)} steps, {len(rep.get('selftest', []))} selftest rows)  -> {out_dir.relative_to(ROOT).as_posix()}")
            for msg in rep["failures"][:25]:
                print("   -", msg)
            if len(rep["failures"]) > 25:
                print(f"   … {len(rep['failures']) - 25} more in audit.json")
        if a.page:
            p = Path(a.page) if Path(a.page).is_absolute() else ROOT / a.page
            out_dir = ROOT / "reports" / "viz" / "_pages" / p.stem
            rep = audit_page(browser, p, out_dir)
            ok &= rep["verdict"] == "PASS"
            print(f"{rep['verdict']}  page {rep['page']}  -> {out_dir.relative_to(ROOT).as_posix()}")
            for msg in rep["failures"]:
                print("   -", msg)
        browser.close()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
