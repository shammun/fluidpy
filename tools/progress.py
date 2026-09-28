"""Update one chapter's phase status in progress.json (keeps LF line endings and formatting).

Usage:
    .venv/Scripts/python.exe tools/progress.py ch08 analyze pass --note "analyze: ..."
    .venv/Scripts/python.exe tools/progress.py ch08 explainers a,b,c
    .venv/Scripts/python.exe tools/progress.py ch08 show
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROGRESS = ROOT / "progress.json"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("chapter")
    ap.add_argument("phase")
    ap.add_argument("status", nargs="?")
    ap.add_argument("--note", default=None, help="text appended to notes ('; ' separated)")
    ap.add_argument("--current", action="store_true", help="also set current_chapter")
    a = ap.parse_args()

    data = json.loads(PROGRESS.read_text(encoding="utf-8"))
    chapters = data.get("chapters", data)
    ch = chapters[a.chapter]
    if a.phase == "show":
        print(json.dumps(ch, indent=1, ensure_ascii=False))
        return
    if a.phase == "explainers":
        ch["explainers"] = [s for s in (a.status or "").split(",") if s]
    elif a.status is not None:
        ch[a.phase] = a.status
    if a.note:
        ch["notes"] = f"{ch['notes']}; {a.note}" if ch.get("notes") else a.note
    ch["updated"] = dt.date.today().isoformat()
    if a.current:
        data["current_chapter"] = a.chapter
    with PROGRESS.open("w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    print(f"{a.chapter}.{a.phase} = {ch.get(a.phase)}")


if __name__ == "__main__":
    main()
