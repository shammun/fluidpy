"""Write the cited public benchmark tables used by chapter 10 into ``reference/ch10/``.

* ``ghia1982_table1.csv`` — u-velocity along the vertical line through the geometric centre of the lid-driven cavity
  (Ghia, Ghia & Shin 1982, Table I), all seven Reynolds numbers.  Transcription: the public gist by ivan-pi (raw text
  downloaded 2026-09-30); cross-checked against an independent listing of Ghia's centreline minima (Tech Science Press CMC 36,
  1-21 (2013), Table 1: Re = 100 -0.21090 at y = 0.4531; Re = 400 -0.32726 at 0.2813; Re = 1000 -0.38289 at 0.1719; Re = 3200
  -0.41933 at 0.1016 — all four agree digit for digit).  The Re = 3200 entry at y = 0.4531 (-0.86636) is a known typo of the
  original table; it is kept verbatim and flagged in ``SOURCES.md`` (never used).
* ``cavity_vortex_centres.csv`` — primary-vortex centres: Ghia et al. (1982) via Hajabdollahi & Premnath, arXiv:1202.6351,
  table of primary-vortex locations; Hou, Zou, Chen, Doolen & Cogley (1995), J. Comput. Phys. 118, 329, preprint
  arXiv:comp-gas/9401003, Fig. 1 captions (lattice Boltzmann, 256^2).  Both PDFs downloaded and text-extracted 2026-09-30.

Nothing here comes from the textbook (book-quoted values live in the git-ignored ``tests/book_values_ch10.json``).
Run: ``.venv/Scripts/python.exe reference/ch10/make_refs.py`` (idempotent).
"""
from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent

# y, Re = 100, 400, 1000, 3200, 5000, 7500, 10000   (Ghia, Ghia & Shin, J. Comput. Phys. 48, 387-411 (1982), Table I)
GHIA_TABLE1 = """\
1.0000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000
0.9766  0.84123  0.75837  0.65928  0.53236  0.48223  0.47244  0.47221
0.9688  0.78871  0.68439  0.57492  0.48296  0.46120  0.47048  0.47783
0.9609  0.73722  0.61756  0.51117  0.46547  0.45992  0.47323  0.48070
0.9531  0.68717  0.55892  0.46604  0.46101  0.46036  0.47167  0.47804
0.8516  0.23151  0.29093  0.33304  0.34682  0.33556  0.34228  0.34635
0.7344  0.00332  0.16256  0.18719  0.19791  0.20087  0.20591  0.20673
0.6172 -0.13641  0.02135  0.05702  0.07156  0.08183  0.08342  0.08344
0.5000 -0.20581 -0.11477 -0.06080 -0.04272 -0.03039 -0.03800  0.03111
0.4531 -0.21090 -0.17119 -0.10648 -0.86636 -0.07404 -0.07503 -0.07540
0.2813 -0.15662 -0.32726 -0.27805 -0.24427 -0.22855 -0.23176 -0.23186
0.1719 -0.10150 -0.24299 -0.38289 -0.34323 -0.33050 -0.32393 -0.32709
0.1016 -0.06434 -0.14612 -0.29730 -0.41933 -0.40435 -0.38324 -0.38000
0.0703 -0.04775 -0.10338 -0.22220 -0.37827 -0.43643 -0.43025 -0.41657
0.0625 -0.04192 -0.09266 -0.20196 -0.35344 -0.42901 -0.43590 -0.42537
0.0547 -0.03717 -0.08186 -0.18109 -0.32407 -0.41165 -0.43154 -0.42735
0.0000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
"""

CENTRES = [
    # source, method, Re, x, y
    ("Ghia1982", "multigrid stream function-vorticity", 100, 0.6172, 0.7344),
    ("Ghia1982", "multigrid stream function-vorticity", 400, 0.5547, 0.6055),
    ("Hou1995", "lattice Boltzmann 256^2; lid U = 0.01", 100, 0.6196, 0.7373),
    ("Hou1995", "lattice Boltzmann 256^2; lid U = 0.1", 400, 0.5608, 0.6078),
]


# Verifier cross-check (2026-09-30), each value read today from a source independent of the gist above:
# (Re, y) -> u.  CMC 36 (2013) Table 1 lists Ghia's centreline minima (ref. 2); the cfdgasman/lid-driven-cavity repository
# (cavity/ghia.py, raw file downloaded 2026-09-30) is an independent transcription of Table I at Re = 100 (all 17 values agree).
CROSS_CHECK = {
    (100, 0.4531): -0.21090,   # CMC 36 (2013) Table 1, ref. 2; cfdgasman ghia.py U_RE100
    (400, 0.2813): -0.32726,   # CMC 36 (2013) Table 1, ref. 2
    (1000, 0.1719): -0.38289,  # CMC 36 (2013) Table 1, ref. 2
    (100, 0.9766): 0.84123,    # cfdgasman ghia.py U_RE100
    (100, 0.5000): -0.20581,   # cfdgasman ghia.py U_RE100
}


def cross_check() -> None:
    """Assert the transcribed table against the independently read values of CROSS_CHECK."""
    rows = [line.split() for line in GHIA_TABLE1.strip().splitlines()]
    col = {100: 1, 400: 2, 1000: 3}
    for (re, y), val in CROSS_CHECK.items():
        row = [r for r in rows if abs(float(r[0]) - y) < 1e-9][0]
        assert float(row[col[re]]) == val, (re, y, row[col[re]], val)


def main() -> None:
    cross_check()
    rows = [line.split() for line in GHIA_TABLE1.strip().splitlines()]
    out = ["# Ghia, Ghia & Shin (1982) J. Comput. Phys. 48, 387-411, Table I: u on x = 0.5 (lid speed 1, side 1); "
           "Re=3200 at y=0.4531 is a known typo of the original (kept verbatim, do not use)",
           "y,Re100,Re400,Re1000,Re3200,Re5000,Re7500,Re10000"]
    out += [",".join(r) for r in sorted(rows, key=lambda r: float(r[0]))]
    (HERE / "ghia1982_table1.csv").write_text("\n".join(out) + "\n", encoding="utf-8")
    c = ["# primary-vortex centres of the lid-driven cavity (see SOURCES.md)", "source,method,Re,x,y"]
    c += [f"{s},{m},{re},{x},{y}" for s, m, re, x, y in CENTRES]
    (HERE / "cavity_vortex_centres.csv").write_text("\n".join(c) + "\n", encoding="utf-8")
    print("wrote", HERE / "ghia1982_table1.csv", "and", HERE / "cavity_vortex_centres.csv")


if __name__ == "__main__":
    main()
