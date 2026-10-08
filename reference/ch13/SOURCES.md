# Chapter 13 — reference data and where every number comes from

Nothing in this folder comes from the textbook or from any external source: every file is **our own model output**,
written by `scripts/ch13_make_caches.py` (which calls `fluidpy.ch13_geophysical_fluid_dynamics.write_reference_runs`),
read back by the same script and compared with the run that produced it. Book-quoted numbers used as test evidence live
only in the git-ignored `tests/book_values_ch13.json`.

| file | what it is | produced by | how to regenerate |
|---|---|---|---|
| `kelvin_basin.npz` | a Kelvin pulse travelling round a closed square basin 10 Rossby radii wide (C-grid shallow-water model, 64 × 64, equivalent depth 1.3 m, 35° N); 24 frames of η as float16 | `ch13.reference_run("kelvin_basin")` | `scripts/ch13_make_caches.py --which kelvin_basin` |
| `turbulence_f.npz` | decaying two-dimensional turbulence, β = 0 (pseudo-spectral barotropic model, 64 × 64, seed 3, non-dimensional); 10 frames of vorticity as float16 and the energy, enstrophy and centroid series | `ch13.reference_run("turbulence_f")` | `--which turbulence_f` |
| `turbulence_beta.npz` | the same initial field with β = 8 (zonal bands) | `ch13.reference_run("turbulence_beta")` | `--which turbulence_beta` |
| `pv_particles.npz` | a nonlinear geostrophic vortex (48 × 48) with 24 marked particles and the potential vorticity they carry; 13 frames | `ch13.reference_run("pv_particles")` | `--which pv_particles` |
| `explainer_constants.json` | constants the explainers compute themselves (Eady cut-off and fastest wave, adjustment energy ratio, the illustrative inputs), for parity rows | `ch13.explainer_constants()` | any run of the script |

All runs are deterministic (seeded). Labels: illustrations of ours ("qualitative"); the conservation and convergence
numbers of the models are in the docstrings of `fluidpy/core/shallow_water.py`.
