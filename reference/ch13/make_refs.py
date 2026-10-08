"""Write ``reference/ch13/benchmarks.json`` — the cited public numbers the chapter-13 tests compare with.

Only values read first-hand in a public source on the date given are stored (see SOURCES.md for the full citation and
what was read).  Nothing here comes from the textbook.  The cached model runs in this folder (``*.npz``,
``explainer_constants.json``) are our own output and are written by ``scripts/ch13_make_caches.py``, not by this file.

Run: ``.venv/Scripts/python.exe reference/ch13/make_refs.py``  (idempotent: rewrites the same bytes).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_note": "Cited public values for the Eady problem; our own numbers are computed, never typed. See SOURCES.md.",
    "eady_max_growth_coefficient": {
        "value": 0.3098,
        "digits": 4,
        "what": "largest growth rate of the Eady wave in units of f (dU/dz) / N (l = 0)",
        "source": "K. A. Emanuel, MIT OpenCourseWare 12.803 Quasi-Balanced Circulations in Oceans and Atmospheres "
                  "(Fall 2009), Lecture 19 'Baroclinic Instability', text after Eq. (19.17): r c_i = 0.3098",
        "url": "https://ocw.mit.edu/courses/12-803-quasi-balanced-circulations-in-oceans-and-atmospheres-fall-2009/",
        "second_source": "NCAR Command Language, function eady_growth_rate (coefficient 0.3098; cites Eady 1949, "
                         "Tellus 1, 33-52, doi:10.1111/j.2153-3490.1949.tb01265.x, and Lindzen & Farrell 1980, "
                         "J. Atmos. Sci. 37, 1648-1654)",
        "second_url": "https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml",
        "how_obtained": "read in the lecture-note PDF (page 7) and on the NCL documentation page",
        "verified": "2026-10-07",
    },
    "eady_fastest_wavenumber": {
        "value": 1.606,
        "digits": 4,
        "what": "non-dimensional wavenumber (N H K / f) of the fastest-growing Eady wave",
        "source": "Emanuel, MIT OCW 12.803 Lecture 19, text after Eq. (19.17): extremum at r = 1.606",
        "url": "https://ocw.mit.edu/courses/12-803-quasi-balanced-circulations-in-oceans-and-atmospheres-fall-2009/",
        "how_obtained": "read in the lecture-note PDF (page 7)",
        "verified": "2026-10-07",
    },
    "eady_cutoff_wavenumber_rounded": {
        "value": 2.4,
        "digits": 2,
        "what": "non-dimensional wavenumber beyond which Eady waves are neutral, as rounded by the source",
        "source": "Emanuel, MIT OCW 12.803 Lecture 19, text after Eq. (19.17) and under Figure 19.2",
        "url": "https://ocw.mit.edu/courses/12-803-quasi-balanced-circulations-in-oceans-and-atmospheres-fall-2009/",
        "how_obtained": "read in the lecture-note PDF (pages 7 and 8); two digits only, so the tolerance is half a unit "
                        "of the last digit",
        "verified": "2026-10-07",
    },
}


def main() -> int:
    path = HERE / "benchmarks.json"
    path.write_text(json.dumps(BENCHMARKS, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {path} ({path.stat().st_size} bytes, {len(BENCHMARKS) - 1} benchmark values)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
