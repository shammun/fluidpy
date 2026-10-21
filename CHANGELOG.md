# Changelog

All notable changes to the `fluidpy` library. Versions are git tags (`vX.Y.Z`); Colab notebooks of the AOFD book pin a
tag. Before 1.0, a minor version may add but never remove or rename anything that a Kundu notebook or a published AOFD
chapter imports.

## v0.3.0 — 2026-10-20

- Packaging: `pyproject.toml` (setuptools), so sibling book repos can `pip install -e "../Fluid Dynamics[geo,data,diag]"`
  and Colab can `pip install "fluidpy[geo,data] @ git+https://github.com/shammun/fluidpy@v0.3.0"`.
- Extras: `[geo]` xarray/dask, `[data]` pooch, netCDF4, h5netcdf, zarr, fsspec, gcsfs, cdsapi, earthaccess, gsw, cartopy,
  cmocean, metpy; `[diag]` scikit-image; `[dev]` pytest.
- New empty subpackages `fluidpy.geo`, `fluidpy.data`, `fluidpy.diag` (Vallis-era physics, data access, diagnostics) with
  test folders `tests/geo`, `tests/data`, `tests/diag`.
- Safety: `tools/check_secrets.py` and git hooks block credentials from commits and pushes.
- No change to `fluidpy.core` or the chapter modules.

## v0.2.0

- The Kundu (Fluid Mechanics 5e) library as used by the finished fluidpy site, chapters 1–16.
