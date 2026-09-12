---
name: data-and-benchmarks
description: Data and reference-value policy for fluidpy — analytic/synthetic first, then cited published benchmarks, then free online datasets, then manual download into Google Drive; how to fetch, cache, cite and access them on Colab. Load whenever a script, test or notebook needs data or a reference number.
---

# data-and-benchmarks — where numbers come from

Most of a fluid-dynamics textbook needs **no data at all**: the inputs are parameters, and the outputs are fields and
numbers you compute. Data enters in three places — validation benchmarks, fluid properties, and the occasional real
dataset for an application section.

## Tier order (stop at the first tier that satisfies the need)
**Tier 1 — Analytic / synthetic.** Generate it: an exact solution sampled on a grid, a manufactured solution, a
synthetic velocity field with known statistics. Deterministic and seeded, in `fluidpy/core/synth.py` or in the test
itself. This is the default for everything the book derives.

**Tier 2 — Published benchmarks (cited).** Digitised tables from a paper, a standard or a database. They live in
`reference/chNN/` and are produced by `reference/chNN/make_refs.py`, with one row per value in
`reference/chNN/SOURCES.md`:

```
| value/table | what it is | source (authors, journal, year) | DOI/URL | how obtained | verified |
| ghia_re100_u.csv | u on the vertical centreline, lid-driven cavity | Ghia, Ghia & Shin, J. Comput. Phys. 48, 387 (1982) | doi:10.1016/0021-9991(82)90058-4 | table II, typed by hand, checked twice | 2026-09-12 |
```
Rules: **never quote a benchmark from memory** — fetch or cite the source and record the date; say "digitised from a
plot" when it was, and use the looser tolerance; keep the values in a plain `.csv`/`.json` so a reader can check them;
`tools/benchmarks.py` holds the small scalar constants (Blasius f''(0), κ, B, Ra_c, Re_c …) with their citations in the
docstring.

**Tier 3 — Free online datasets (no login).** For application sections: NOAA/NASA/ECMWF open data, the Johns Hopkins
Turbulence Database, NIST thermophysical properties (or CoolProp), USGS river discharge, standard-atmosphere tables.
Add a cached fetcher in `fluidpy/core/refdata.py`:
```python
def fetch(url: str, dest: Path, timeout: int = 60) -> Path:  # cached download, clear error message, records the date
```
Verify every URL with a real request before relying on it and record it in `data/online/SOURCES.md` with its licence.
Prefer public-domain (US federal) sources. Nothing under `data/online/` is committed except `SOURCES.md`.

**Tier 4 — Manual download (login required) → Google Drive.** Emit a precise instruction block and make the code look
for the file with a helpful error:
```
MANUAL DATA NEEDED for chapter N:
  1. Go to <url>; log in (free account).
  2. Search: <exact parameters>; download <file pattern>.
  3. Save to Google Drive: MyDrive/fluidpy/data/manual/chNN/<filename>
  4. In Colab, mount Drive and the notebook will find it at /content/drive/MyDrive/fluidpy/data/manual/chNN/.
  Locally: put it in data/manual/chNN/ (git-ignored).
```

## Values taken from the book itself
These are **not** a data tier — they are private evidence. Put them in `tests/book_values_chNN.json` (git-ignored),
guard the tests that use them with a `skipif`, and report only derived statements ("reproduces Example 5.3 to 0.2 %").
Never put a book table into `reference/`, a notebook, or the published HTML.

## Fluid properties
Use one source, consistently, and name it: a small table in `fluidpy/core/thermo.py` (water and air at a few
temperatures, with the source), or CoolProp if the book needs real equations of state. A property mismatch is a common
reason a "correct" implementation misses a book number by 1–2 %, so the report must say which property set was used.

## Locations and Colab access
- Local: `reference/` (committed, cited), `data/online/` (git-ignored except `SOURCES.md`), `data/synthetic/`,
  `data/manual/` and `data/book/` (git-ignored, private).
- Everything loads through `fluidpy.core.refdata.load(chapter, name)`: `reference/chNN/` → cached download →
  `data/manual/chNN/` → Drive path on Colab → a clear error that says exactly what to do.
- Colab: cell 1 optionally mounts Drive (only for Tier 4), cell 2 clones the repo, cell 3 loads what the chapter needs.

## Never
Fabricate a benchmark value or a citation, silently swap one dataset for another, embed base64 data in notebooks, or
commit anything book-derived.
