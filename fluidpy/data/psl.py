"""NOAA PSL: NCEP/NCAR Reanalysis 1 daily means on pressure levels (dataset key ``NCEP1-hgt-daily``).

Source
    NOAA Physical Sciences Laboratory (PSL), NCEP/NCAR Reanalysis 1, daily averages on 17 pressure levels, 2.5° grid.
    OPeNDAP ``https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/{var}.{YYYY}.nc``
    (server-side subset; preferred), HTTPS ``https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/
    pressure/{var}.{YYYY}.nc`` (whole year, 40–165 MB; fallback when THREDDS returns 503).
Licence
    US Government work, no restrictions (public domain). PSL asks for the acknowledgement "NCEP Reanalysis data provided
    by the NOAA PSL, Boulder, Colorado, USA, from their website at https://psl.noaa.gov".
Citation
    Kalnay, E., et al., 1996: The NCEP/NCAR 40-Year Reanalysis Project. Bull. Amer. Meteor. Soc., 77, 437–471.
Units
    ``hgt`` geopotential height [m]; ``level`` [hPa]; ``lat``/``lon`` [degrees]. PSL latitudes run N→S: every reader
    here sorts them S→N.
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import numpy as np

from . import fetch, data_root
from .reduce import provenance, write_nc

__all__ = ["OPENDAP", "HTTPS", "daily_url", "fetch_daily", "open_daily", "z500_day", "ATTRIBUTION"]

OPENDAP = "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/"
HTTPS = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/"
KEY = "NCEP1-hgt-daily"
LICENCE = "US Government work, public domain (no restrictions); acknowledgement requested by NOAA PSL"
CITATION = ("Kalnay, E., et al., 1996: The NCEP/NCAR 40-Year Reanalysis Project. Bull. Amer. Meteor. Soc., 77, "
            "437-471, doi:10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2")
ATTRIBUTION = ("NCEP/NCAR Reanalysis 1 data provided by the NOAA PSL, Boulder, Colorado, USA, from their website at "
               "https://psl.noaa.gov")


def daily_url(var: str, year: int, route: str = "opendap") -> str:
    """URL of one yearly daily-mean file (``route`` = ``"opendap"`` or ``"https"``).

    Book: data for Vallis §2 / Fig. 18.22-type 500 hPa maps (no equation). Units: none [dimensionless].
    Validation: smoke-tested in tests/data/ (no network).
    """
    return (OPENDAP if route == "opendap" else HTTPS) + f"{var}.{int(year)}.nc"


def fetch_daily(date: str, level: float = 500.0, var: str = "hgt", route: str = "auto") -> Path:
    """Fetch one day of one NCEP R1 daily-mean variable on one pressure level, through the disk guard.

    Parameters
    ----------
    date : str
        ``"YYYY-MM-DD"``.
    level : float
        Pressure level [hPa].
    var : str
        PSL variable name (``hgt``, ``air``, ``uwnd``, ``vwnd``, ``omega``, ``rhum``, ``shum``).
    route : {"auto", "opendap", "https"}
        ``auto`` tries the OPeNDAP subset (≈ 42 kB) and falls back to the whole yearly file over HTTPS.

    Returns
    -------
    pathlib.Path
        The subset file (OPeNDAP) or the yearly file (HTTPS) under ``<data_root>/raw/NCEP1-hgt-daily/``.

    Book: data for Vallis §2 (geostrophic wind from heights) and Fig. 18.22-type maps; no equation.
    Units: height [m], level [hPa]. Validation: smoke-tested in tests/data/; real fetch logged in reports/data/.
    """
    year = int(date[:4])
    err = None
    if route in ("auto", "opendap"):
        try:
            return fetch(daily_url(var, year, "opendap"), key=KEY,
                         subset={"vars": [var], "sel": {"time": date, "level": float(level)}},
                         fname=f"{var}.{date}.{int(level)}hPa.nc")
        except Exception as e:  # THREDDS sometimes answers 503
            err = e
            if route == "opendap":
                raise
            print(f"[fluidpy.data.psl] OPeNDAP failed ({e!r:.200}); falling back to the whole yearly file (HTTPS)")
    p = fetch(daily_url(var, year, "https"), key=KEY)
    if err is not None:
        print(f"[fluidpy.data.psl] used HTTPS fallback {p}")
    return p


def open_daily(path: str | Path, var: str = "hgt", date: str | None = None, level: float | None = None):
    """Open a fetched file as an xarray DataArray with latitude sorted S→N (PSL stores N→S).

    Book: data for Vallis §2; no equation. Units: as stored (``hgt`` [m]). Validation: smoke-tested in tests/data/.
    """
    import xarray as xr

    da = xr.open_dataset(path)[var]
    if date is not None and "time" in da.dims:
        da = da.sel(time=date)
    if level is not None and "level" in da.dims:
        da = da.sel(level=float(level))
    return da.sortby("lat")


def z500_day(date: str, out: str | Path | None = None, route: str = "auto") -> Path:
    """Derived product: global 2.5° NCEP R1 daily-mean 500 hPa geopotential height for one day (float32, zlib).

    Parameters
    ----------
    date : str
        ``"YYYY-MM-DD"``.
    out : path, optional
        Output file (default ``<data_root>/derived/ch00/ncep_r1_z500_<YYYYMMDD>.nc``).
    route : str
        Passed to :func:`fetch_daily`.

    Returns
    -------
    pathlib.Path
        The derived netCDF: variable ``z`` [m] on (lat S→N, lon 0–357.5 E), attrs with provenance.

    Book: Vallis §2 (geostrophic balance, u_g = -(g/f) dZ/dy, v_g = (g/f) dZ/dx) uses it in PoC F; no equation here.
    Units: z [m], lat/lon [degrees]. Validation: values checked 4600-6000 m range in tests/data/ when present.
    """
    src = fetch_daily(date, 500.0, "hgt", route=route)
    da = open_daily(src, "hgt", date=date, level=500.0)
    da = da.drop_vars([c for c in ("time", "level") if c in da.coords]).astype("float32")
    ds = da.rename("z").to_dataset()
    ds["z"].attrs = {"long_name": "geopotential height of the 500 hPa surface (daily mean)", "units": "m",
                     "standard_name": "geopotential_height"}
    ds["lat"].attrs = {"units": "degrees_north", "long_name": "latitude (sorted south to north)"}
    ds["lon"].attrs = {"units": "degrees_east", "long_name": "longitude (0 to 357.5 E)"}
    if not np.all(np.isfinite(ds["z"].values)):
        raise ValueError("z500 has missing values")
    url = daily_url("hgt", int(date[:4]), "opendap" if "hPa" in Path(src).name else "https")
    attrs = provenance(KEY, url, LICENCE, CITATION,
                       f"daily mean {date} (00-24 UTC), level 500 hPa selected server-side (OPeNDAP); latitude sorted "
                       "S->N; cast to float32", date=date, level_hPa="500", attribution=ATTRIBUTION,
                       title=f"NCEP/NCAR R1 500 hPa geopotential height, daily mean {date}",
                       fetched=_dt.date.today().isoformat())
    out = Path(out) if out else data_root() / "derived" / "ch00" / f"ncep_r1_z500_{date.replace('-', '')}.nc"
    return write_nc(ds, out, attrs)
