# Changelog

All notable changes to the `fluidpy` library. Versions are git tags (`vX.Y.Z`); Colab notebooks of the AOFD book pin a
tag. Before 1.0, a minor version may add but never remove or rename anything that a Kundu notebook or a published AOFD
chapter imports.

## v0.3.1 — 2026-10-22

For aofd ch00 PoC. Additive only: nothing removed or renamed; no change to `fluidpy.core` or the Kundu chapter modules.

### Added
- `fluidpy.geo.vectors`: `dot`, `cross`, `norm`, `angle_between`, `scalar_projection`, `projection`,
  `parallelogram_area`, `triple_scalar`, `k_cross`.
- `fluidpy.geo.frames`: `coriolis_parameter`, `beta_parameter`, `coriolis_acceleration_horizontal`,
  `coriolis_horizontal`, `turntable_puck`, `inertial_circle`, `inertial_oscillation`, `inertial_oscillation_sphere`,
  `to_inertial_lon`, `great_circle`, `great_circle_earth_frame`.
- `fluidpy.geo.atmos`: `geostrophic_wind_from_height`, `scale_height`, `hypsometric_thickness`.
- `fluidpy.geo.constants`: `OMEGA`, `EARTH_RADIUS`, `R_D`, `C_P`, `KAPPA`, `RHO0_OCEAN` (each with its source).
- `fluidpy.data` fetch door: `fetch` (every download passes the disk guard), `DiskGuardError`, `data_root`,
  `aofd_repo`, `disk_report`, `sha256`, `open_remote` (OPeNDAP via `engine="pydap"` on Windows);
  `data.registry` (`registry_path`, `load`, `get`, `evictable`); `data.reduce` (`provenance`, `write_nc`, `write_npz`).
- `fluidpy.data` fetchers (source, licence, citation in each module): `psl` (NCEP/NCAR Reanalysis 1 daily, `z500_day`),
  `naturalearth` (`land_50m`, `bathymetry_band`, `shallow_mask`), `fronts` (Orsi ACC fronts, `orsi_fronts`),
  `gdp` (Global Drifter Program mean velocity via ERDDAP, `southern_ocean_mean_velocity`).

### Tests
- `tests/geo/test_vectors.py`, `tests/geo/test_frames.py`, `tests/geo/test_atmos.py`, `tests/data/test_fetch.py`
  (no downloads in tests). Full suite: 2273 passed (all Kundu chapters included).

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
