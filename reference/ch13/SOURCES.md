# Chapter 13 — reference data and where every number comes from

Two kinds of file live here, and neither comes from the textbook:

1. **Our own model output** (`*.npz`, `explainer_constants.json`), written by `scripts/ch13_make_caches.py` (which calls
   `fluidpy.ch13_geophysical_fluid_dynamics.write_reference_runs`), read back by the same script and compared with the run
   that produced it.
2. **Cited public benchmark values** (`benchmarks.json`), written by `make_refs.py`; each value was read first-hand in the
   source on the date in the last column.

Book-quoted numbers used as test evidence live only in the git-ignored `tests/book_values_ch13.json`.

## 1. Our own model output

| file | what it is | produced by | how to regenerate |
|---|---|---|---|
| `kelvin_basin.npz` | a Kelvin pulse travelling round a closed square basin 10 Rossby radii wide (C-grid shallow-water model, 64 × 64, equivalent depth 1.3 m, 35° N); 24 frames of η as float16 | `ch13.reference_run("kelvin_basin")` | `scripts/ch13_make_caches.py --which kelvin_basin` |
| `turbulence_f.npz` | decaying two-dimensional turbulence, β = 0 (pseudo-spectral barotropic model, 64 × 64, seed 3, non-dimensional); 10 frames of vorticity as float16 and the energy, enstrophy and centroid series | `ch13.reference_run("turbulence_f")` | `--which turbulence_f` |
| `turbulence_beta.npz` | the same initial field with β = 8 (zonal bands) | `ch13.reference_run("turbulence_beta")` | `--which turbulence_beta` |
| `pv_particles.npz` | a nonlinear geostrophic vortex (48 × 48) with 24 marked particles and the potential vorticity they carry; 13 frames | `ch13.reference_run("pv_particles")` | `--which pv_particles` |
| `explainer_constants.json` | constants the explainers compute themselves (Eady cut-off and fastest wave, adjustment energy ratio, the illustrative inputs), for parity rows | `ch13.explainer_constants()` | any run of the script |

All runs are deterministic (seeded). Labels: illustrations of ours ("qualitative"); the conservation and convergence
numbers of the models are in the docstrings of `fluidpy/core/shallow_water.py` and are re-measured in
`tests/test_ch13.py` (the files regenerate to the float16 storage error, 2 × 10⁻³ of each array's range).

## 2. Cited public benchmark values (`benchmarks.json`, written by `make_refs.py`)

| value/table | what it is | source (authors, journal, year) | DOI/URL | how obtained | verified |
|---|---|---|---|---|---|
| `eady_max_growth_coefficient` = 0.3098 | largest growth rate of the Eady wave in units of f (dU/dz)/N | K. A. Emanuel, MIT OpenCourseWare 12.803 *Quasi-Balanced Circulations in Oceans and Atmospheres* (Fall 2009), Lecture 19 "Baroclinic Instability", text after Eq. (19.17) | https://ocw.mit.edu/courses/12-803-quasi-balanced-circulations-in-oceans-and-atmospheres-fall-2009/ | read in the lecture-note PDF, page 7 ("a value of r c_i of 0.3098") | 2026-10-07 |
| the same coefficient, second source | the coefficient in the NCL function `eady_growth_rate` | NCAR Command Language documentation; it cites Eady (1949), Tellus 1, 33–52, doi:10.1111/j.2153-3490.1949.tb01265.x, and Lindzen & Farrell (1980), J. Atmos. Sci. 37, 1648–1654 | https://www.ncl.ucar.edu/Document/Functions/Contributed/eady_growth_rate.shtml | read on the documentation page (the two papers themselves were **not** read) | 2026-10-07 |
| `eady_fastest_wavenumber` = 1.606 | non-dimensional wavenumber N H K/f of the fastest-growing Eady wave | Emanuel, MIT OCW 12.803 Lecture 19, text after Eq. (19.17) | as above | read in the lecture-note PDF, page 7 ("an extremum when r = 1.606") | 2026-10-07 |
| `eady_cutoff_wavenumber_rounded` = 2.4 | non-dimensional wavenumber beyond which Eady waves are neutral, as the source rounds it | Emanuel, MIT OCW 12.803 Lecture 19, text after Eq. (19.17) and under its Figure 19.2 | as above | read in the lecture-note PDF, pages 7–8; two digits only, so the test tolerance is half a unit of the last digit | 2026-10-07 |

Our values (computed, not typed): 0.30982, 1.60612 and 2.39936. They agree with the source to every digit it prints.

## 3. Not used as benchmarks (and why)

- The rotation rate of the Earth and its mean radius are constants of the code (`ROT.OMEGA_EARTH`, `GFD.EARTH_RADIUS_MEAN`);
  comparing a constant with a published value exercises no code and is not counted as evidence.
- The end state of geostrophic adjustment (one third of the released potential energy stays in the jet) is derived and
  tested as "analytic (ours)". The classical treatment (Gill, *Atmosphere–Ocean Dynamics*, 1982, §7.2–7.3) was **not**
  read first-hand, so it is not cited as a benchmark.
- Growth rates of the sech² jet on a β-plane are ours, computed twice (Chebyshev collocation on a complex path and an
  independent shooting integration in the test file). Kuo (1949) and later tabulations were not read.
- Kraichnan (1967) and Rhines (1975) (the −3 range and the Rhines length) were not read; the exponents are checked by
  dimensional analysis only.
