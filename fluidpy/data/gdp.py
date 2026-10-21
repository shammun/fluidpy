"""NOAA Global Drifter Program: drifter-derived mean near-surface velocity climatology (key ``CURR-gdp-annual``).

Source
    NOAA AOML ERDDAP, griddap dataset ``drifter_annualmeans`` (0.25°, 72.875°S–84.875°N; U, V, standard errors and
    drifter-days N), ``https://erddap.aoml.noaa.gov/gdp/erddap/griddap/drifter_annualmeans``. Subset server-side.
    Product page: https://www.aoml.noaa.gov/phod/gdp/mean_velocity.php
Licence
    CC-BY 4.0 per the dataset's ``license`` attribute ("This study used data collected and made freely available by the
    NOAA Global Drifter Program (https://www.aoml.noaa.gov/phod/gdp/) accessed from
    https://erddap.aoml.noaa.gov/gdp/erddap/index.html").
Citation
    Laurindo, L., A. Mariano and R. Lumpkin, 2017: An improved near-surface velocity climatology for the global ocean
    from drifter observations. Deep-Sea Res. I, 124, 73-92, doi:10.1016/j.dsr.2017.04.009.
Units
    U, V [m s⁻¹] (eastward, northward); N [drifter-days per square degree]; lon/lat [degrees].
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import numpy as np

from . import fetch, data_root
from .reduce import provenance, write_nc

__all__ = ["ERDDAP", "query", "fetch_mean_velocity", "southern_ocean_mean_velocity", "ATTRIBUTION"]

ERDDAP = "https://erddap.aoml.noaa.gov/gdp/erddap/griddap/drifter_annualmeans.nc"
KEY = "CURR-gdp-annual"
LICENCE = "CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/), per the ERDDAP dataset licence attribute"
CITATION = ("Laurindo, L., Mariano, A., Lumpkin, R., 2017: An improved near-surface velocity climatology for the "
            "global ocean from drifter observations. Deep-Sea Res. I, 124, 73-92, doi:10.1016/j.dsr.2017.04.009")
ATTRIBUTION = ("Drifter mean velocity: NOAA Global Drifter Program (AOML), Laurindo et al. (2017); CC BY 4.0")


def query(vars=("U", "V", "eU", "eV", "N"), lat=(-72.875, -30.125), lon=(-179.875, 179.875)) -> str:
    """ERDDAP griddap constraint string (axis order longitude, latitude) for a lon/lat box [degrees].

    Book: data for Vallis Fig. 21.23 (no equation). Units: [degrees]. Validation: smoke-tested in tests/data/.
    """
    c = f"[({lon[0]}):1:({lon[1]})][({lat[0]}):1:({lat[1]})]"
    return ",".join(v + c for v in vars)


def fetch_mean_velocity(lat=(-72.875, -30.125), vars=("U", "V", "eU", "eV", "N")) -> Path:
    """Fetch the annual-mean drifter velocity for a latitude band, subset on the ERDDAP server, via the disk guard.

    Size estimate = number of cells × variables × 8 bytes (ERDDAP serves float64).

    Book: Vallis §21 Fig. 21.23 (data-driven twin panel); no equation. Units: U, V [m s⁻¹].
    Validation: smoke-tested in tests/data/; real fetch logged in reports/data/CURR-gdp-annual.md.
    """
    nlat = int(round((lat[1] - lat[0]) / 0.25)) + 1
    est = nlat * 1440 * len(vars) * 8 + 50_000
    return fetch(ERDDAP, key=KEY, subset=query(vars, lat), expected_bytes=est,
                 fname=f"drifter_annualmeans_{lat[0]}_{lat[1]}.nc")


def southern_ocean_mean_velocity(out: str | Path | None = None, min_n: float = 0.0) -> Path:
    """Derived product: GDP annual-mean U, V (and errors, N) for 72.875°S–30.125°S at 0.25° (float32, zlib).

    Cells with no drifter data (N ≤ ``min_n`` or NaN) are NaN in U and V. Dimensions are (lat S→N, lon −180…180).

    Book: Vallis §21 Fig. 21.23; no equation. Units: U, V, eU, eV [m s⁻¹]; N [drifter-days per square degree].
    Validation: tests/data/ check ACC eastward mean flow in Drake Passage when the product is present.
    """
    import xarray as xr

    src = fetch_mean_velocity()
    ds = xr.open_dataset(src).load()
    ds = ds.transpose("latitude", "longitude").rename({"latitude": "lat", "longitude": "lon"}).sortby("lat")
    bad = ~(ds["N"] > min_n)
    for v in ("U", "V", "eU", "eV"):
        if v in ds:
            ds[v] = ds[v].where(~bad)
    for v in ds.data_vars:
        ds[v].attrs.pop("_ChunkSizes", None)
    keep = {k: ds.attrs[k] for k in ("title", "license", "references") if k in ds.attrs}
    ds.attrs = {f"source_{k}": str(v) for k, v in keep.items()}
    attrs = provenance(KEY, ERDDAP + "?" + query(), LICENCE, CITATION,
                       "server-side ERDDAP subset 72.875S-30.125S, all longitudes; transposed to (lat, lon), lat "
                       f"sorted S->N; U,V,eU,eV set NaN where N <= {min_n}; float32", attribution=ATTRIBUTION,
                       fetched=_dt.date.today().isoformat())
    out = Path(out) if out else data_root() / "derived" / "ch21" / "fig_21_23_gdp_mean_velocity.nc"
    p = write_nc(ds, out, attrs)
    if not np.isfinite(xr.open_dataset(p)["U"].values).any():
        raise ValueError("GDP subset has no finite velocities")
    return p
