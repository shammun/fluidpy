"""Execute a chapter notebook headlessly (on the fluidpy-venv kernel) and report errors — without publishing anything.

Usage (repo root)::

    .venv/Scripts/python.exe tools/run_notebook.py ch07            # notebooks/ch07_<slug>.ipynb
    .venv/Scripts/python.exe tools/run_notebook.py ch07 --fast     # FLUIDPY_FAST=1
    .venv/Scripts/python.exe tools/run_notebook.py ch07 --save     # also write outputs/ch07/executed.ipynb to inspect

Prints the runtime, the slowest cells, every error (cell index, exception, first traceback line) and any output larger
than 3 MB. Exit 0 = clean.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))


def main(argv: list[str] | None = None) -> int:
    import nbformat
    import yaml

    import publish_notebook as pub

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("chapter")
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--save", action="store_true")
    a = ap.parse_args(argv)
    row = next(c for c in yaml.safe_load((ROOT / "book.yaml").read_text(encoding="utf-8"))["chapters"] if c["id"] == a.chapter)
    path = ROOT / "notebooks" / f"{a.chapter}_{row['slug']}.ipynb"
    nb, secs = pub.execute(path, fast=a.fast)
    probs = pub.execution_problems(nb)
    print(f"{path.relative_to(ROOT).as_posix()}: executed in {secs:.0f} s ({len(nb.cells)} cells, FAST={a.fast})")
    timings = []
    for i, c in enumerate(nb.cells):
        ex = c.get("metadata", {}).get("execution", {})
        if "iopub.execute_input" in ex and "shell.execute_reply" in ex:
            from datetime import datetime
            t0 = datetime.fromisoformat(ex["iopub.execute_input"].replace("Z", "+00:00"))
            t1 = datetime.fromisoformat(ex["shell.execute_reply"].replace("Z", "+00:00"))
            timings.append(((t1 - t0).total_seconds(), i))
        for o in c.get("outputs", []):
            size = len(json.dumps(o))
            if size > 3e6:
                print(f"   big output: cell {i} {size / 1e6:.1f} MB")
            if o.get("output_type") == "error":
                tb = [line for line in o.get("traceback", []) if line.strip()]
                print(f"   ERROR cell {i}: {o.get('ename')}: {o.get('evalue')}")
                if tb:
                    import re
                    print("      ", re.sub(r"\x1b\[[0-9;]*m", "", tb[-1])[:300])
    for t, i in sorted(timings, reverse=True)[:5]:
        print(f"   slow cell {i}: {t:.1f} s")
    if a.save:
        out = ROOT / "outputs" / a.chapter / "executed.ipynb"
        out.parent.mkdir(parents=True, exist_ok=True)
        nbformat.write(nb, out)
        print("   saved", out.relative_to(ROOT).as_posix())
    if probs:
        print("PROBLEMS:")
        for p in probs:
            print("  -", p)
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
