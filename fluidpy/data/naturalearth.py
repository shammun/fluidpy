"""Natural Earth land and bathymetry polygons, and a coarse "shallower than D metres" mask (keys ``NE-*``).

Source
    Natural Earth (naturalearthdata.com), files on ``https://naturalearth.s3.amazonaws.com/``:
    ``50m_physical/ne_50m_land.zip`` (land polygons) and ``10m_physical/ne_10m_bathymetry_{band}.zip`` (one depth band,
    e.g. ``H_3000`` = sea floor deeper than 3000 m).
Licence
    Public domain ("All versions of Natural Earth raster + vector map data ... are in the public domain"); credit
    optional: "Made with Natural Earth."
Units
    Longitude/latitude [degrees] (WGS84); depth bands [m].
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import numpy as np

from . import fetch, data_root
from .reduce import provenance, write_nc

__all__ = ["BASE", "BANDS", "land_50m", "bathymetry_band", "read_polygons", "shallow_mask", "ATTRIBUTION"]

BASE = "https://naturalearth.s3.amazonaws.com/"
BANDS = {0: "L_0", 200: "K_200", 1000: "J_1000", 2000: "I_2000", 3000: "H_3000", 4000: "G_4000", 5000: "F_5000",
         6000: "E_6000", 7000: "D_7000", 8000: "C_8000", 9000: "B_9000", 10000: "A_10000"}
LICENCE = "Public domain (https://www.naturalearthdata.com/about/terms-of-use/)"
CITATION = "Natural Earth, free vector and raster map data, naturalearthdata.com (v5.1, accessed 2026-10-21)"
ATTRIBUTION = "Made with Natural Earth (public domain)"


def land_50m() -> Path:
    """Fetch ``ne_50m_land.zip`` (≈ 0.46 MB) through the disk guard; returns the local zip path.

    Book: geography for Vallis Fig. 21.23-type maps (no equation). Units: [degrees].
    Validation: smoke-tested in tests/data/; sha256 pinned in data/registry.yaml.
    """
    return fetch("NE-50m-land")


def bathymetry_band(depth_m: int = 3000) -> Path:
    """Fetch one Natural Earth 10 m bathymetry band (polygons of sea floor deeper than ``depth_m``) [m].

    Book: geography for Vallis Fig. 21.23 (shading of depths shallower than 3000 m); no equation.
    Units: depth [m]. Validation: smoke-tested in tests/data/; sha256 pinned for the 3000 m band.
    """
    band = BANDS[int(depth_m)]
    if int(depth_m) == 3000:
        return fetch("NE-10m-bathy-3000")
    return fetch(f"{BASE}10m_physical/ne_10m_bathymetry_{band}.zip", key="NE-10m-bathy")


def read_polygons(zip_path: str | Path):
    """Read all polygons of a Natural Earth shapefile zip with pyshp → a list of shapely geometries.

    Book: geography helper for Vallis Fig. 21.23 (no equation). Units: [degrees]. Validation: smoke-tested in
    tests/data/.
    """
    import shapefile
    from shapely.geometry import shape

    with shapefile.Reader(str(zip_path)) as r:
        return [shape(s.__geo_interface__) for s in r.shapes()]


def _grid_mask(geoms, lon2, lat2):
    import shapely

    m = np.zeros(lon2.shape, bool)
    for g in geoms:
        if g.is_empty:
            continue
        minx, miny, maxx, maxy = g.bounds
        box = (lon2 >= minx) & (lon2 <= maxx) & (lat2 >= miny) & (lat2 <= maxy)
        if not box.any():
            continue
        shapely.prepare(g)
        m[box] |= shapely.contains_xy(g, lon2[box], lat2[box])
    return m


def shallow_mask(depth_m: int = 3000, lat_range=(-90.0, -30.0), res: float = 0.25,
                 out: str | Path | None = None) -> Path:
    """Derived product: a 0.25° class mask south of 30°S — 0 deep ocean (≥ ``depth_m``), 1 ocean shallower, 2 land.

    Grid cell centres are tested against the Natural Earth polygons (50 m land; 10 m bathymetry band), i.e. a
    point-sampled ("coarsened") mask, not an area fraction.

    Returns
    -------
    pathlib.Path
        netCDF with ``depth_class`` (int8) on (lat S→N, lon −180…180) and provenance attrs.

    Book: Vallis Fig. 21.23 (Southern Ocean map, depths < 3000 m shaded); no equation. Units: [degrees], depth [m].
    Validation: smoke-tested in tests/data/ (class counts, Drake Passage deep, Kerguelen shallow) when present.
    """
    import xarray as xr

    land = read_polygons(land_50m())
    deep = read_polygons(bathymetry_band(depth_m))
    lat = np.arange(lat_range[0] + res / 2, lat_range[1], res)
    lon = np.arange(-180 + res / 2, 180, res)
    lon2, lat2 = np.meshgrid(lon, lat)
    is_land = _grid_mask(land, lon2, lat2)
    is_deep = _grid_mask(deep, lon2, lat2)
    cls = np.where(is_land, 2, np.where(is_deep, 0, 1)).astype("int8")
    ds = xr.Dataset({"depth_class": (("lat", "lon"), cls)}, coords={"lat": lat, "lon": lon})
    ds["depth_class"].attrs = {"long_name": f"0 = ocean deeper than {depth_m} m, 1 = ocean shallower than "
                                            f"{depth_m} m, 2 = land", "flag_values": "0 1 2",
                               "flag_meanings": f"deeper_than_{depth_m}m shallower_than_{depth_m}m land"}
    ds["lat"].attrs = {"units": "degrees_north"}
    ds["lon"].attrs = {"units": "degrees_east"}
    attrs = provenance("NE-10m-bathy-3000", f"{BASE}10m_physical/ne_10m_bathymetry_{BANDS[int(depth_m)]}.zip ; "
                       f"{BASE}50m_physical/ne_50m_land.zip", LICENCE, CITATION,
                       f"cell centres of a {res} deg grid ({lat_range[0]}..{lat_range[1]} N) tested for inside the "
                       f"NE 50m land polygons and the NE 10m bathymetry {BANDS[int(depth_m)]} polygons "
                       "(shapely.contains_xy); class = 2 land, 0 deeper, 1 shallower", attribution=ATTRIBUTION,
                       fetched=_dt.date.today().isoformat())
    out = Path(out) if out else data_root() / "derived" / "ch21" / f"fig_21_23_bathy{depth_m}_mask.nc"
    return write_nc(ds, out, attrs)
