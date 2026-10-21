"""Orsi et al. (1995) Antarctic Circumpolar Current fronts (key ``FRONTS-orsi1995``).

Source
    Australian Antarctic Data Centre (AADC) dataset "Fronts of the Antarctic Circumpolar Current - GIS data"
    (https://data.aad.gov.au/metadata/antarctic_circumpolar_current_fronts): original files from A. Orsi (``*.asc``),
    a shapefile and a CSV. The AADC download needs an e-mail address (a one-time link), so the default route reads the
    same original Orsi files from the AADC's own GitHub package ``AustralianAntarcticDivision/orsifronts``
    (``data-raw/{stf,saf,pf,saccf,sbdy}.txt``, pinned to commit 76d4f6d). A shapefile downloaded from AADC by hand can
    be read with :func:`read_shapefile` (pyshp + shapely).
Licence
    CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/), as stated on the AADC record.
Citation
    Orsi, A. H., T. Whitworth III and W. D. Nowlin Jr., 1995: On the meridional extent and fronts of the Antarctic
    Circumpolar Current. Deep-Sea Res. I, 42, 641-673. Data: Orsi, A. H. and Harris, U. (2019) Fronts of the Antarctic
    Circumpolar Current - GIS data, Ver. 1, Australian Antarctic Data Centre.
Units
    Longitude/latitude [degrees]; polylines separated by NaN.
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import numpy as np

from . import fetch, data_root
from .reduce import provenance, write_npz

__all__ = ["FRONTS", "NAMES", "fetch_front", "read_orsi_txt", "read_shapefile", "orsi_fronts", "ATTRIBUTION"]

COMMIT = "76d4f6d901d6d02b12eae4625a080f615a18542a"
RAW = f"https://raw.githubusercontent.com/AustralianAntarcticDivision/orsifronts/{COMMIT}/data-raw/"
FRONTS = ("stf", "saf", "pf", "saccf", "sbdy")
NAMES = {"stf": "Subtropical Front", "saf": "Subantarctic Front", "pf": "Polar Front",
         "saccf": "southern ACC Front", "sbdy": "Southern Boundary of the ACC"}
KEY = "FRONTS-orsi1995"
LICENCE = "CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/)"
CITATION = ("Orsi, A.H., Whitworth III, T., Nowlin Jr., W.D., 1995: On the meridional extent and fronts of the "
            "Antarctic Circumpolar Current. Deep-Sea Res. I, 42(5), 641-673, doi:10.1016/0967-0637(95)00021-W. "
            "Data: Orsi, A.H., Harris, U. (2019) Fronts of the Antarctic Circumpolar Current - GIS data, Ver. 1, "
            "Australian Antarctic Data Centre, https://data.aad.gov.au/metadata/antarctic_circumpolar_current_fronts")
ATTRIBUTION = "ACC fronts: Orsi et al. (1995), via the Australian Antarctic Data Centre (CC BY 4.0)"


def fetch_front(name: str) -> Path:
    """Fetch one Orsi front file (``stf``, ``saf``, ``pf``, ``saccf``, ``sbdy``) through the disk guard.

    Book: Vallis §21 Fig. 21.23 (ACC fronts); no equation. Units: [degrees].
    Validation: sha256 pinned in data/registry.yaml; smoke-tested in tests/data/.
    """
    from . import registry

    ent = registry.get(KEY) or {}
    h = (ent.get("sha256_files") or {}).get(f"{name}.txt")
    return fetch(RAW + f"{name}.txt", key=KEY, known_hash=f"sha256:{h}" if h else None,
                 expected_bytes=30_000)


def read_orsi_txt(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse an Orsi front file: ``lon lat`` rows, segments separated by ``%`` lines → NaN-separated arrays.

    Book: Vallis §21 Fig. 21.23; no equation. Units: lon, lat [degrees]. Validation: smoke-tested in tests/data/.
    """
    lons, lats = [], []
    seg: list[tuple[float, float]] = []

    def flush():
        if seg:
            a = np.array(seg)
            lons.extend(a[:, 0].tolist() + [np.nan])
            lats.extend(a[:, 1].tolist() + [np.nan])
            seg.clear()

    for line in Path(path).read_text(encoding="latin-1").splitlines():
        s = line.strip()
        if not s or s.startswith("%"):
            flush()
            continue
        parts = s.split()
        if len(parts) >= 2:
            seg.append((float(parts[0]), float(parts[1])))
    flush()
    lon, lat = np.array(lons[:-1]), np.array(lats[:-1])
    lon = np.where(np.isfinite(lon), ((lon + 180.0) % 360.0) - 180.0, np.nan)
    return lon, lat


def read_shapefile(path: str | Path) -> dict:
    """Read the AADC fronts shapefile (``.shp`` or a zip) with pyshp → ``{name: shapely geometry}``.

    Book: Vallis §21 Fig. 21.23; no equation. Units: [degrees]. Validation: not exercised by default (needs the
    manual AADC download).
    """
    import shapefile
    from shapely.geometry import shape

    out = {}
    with shapefile.Reader(str(path)) as r:
        fields = [f[0] for f in r.fields[1:]]
        for sr in r.shapeRecords():
            rec = dict(zip(fields, sr.record))
            name = str(rec.get("NAME") or rec.get("name") or rec.get("FRONT") or len(out))
            out[name] = shape(sr.shape.__geo_interface__)
    return out


def orsi_fronts(out: str | Path | None = None) -> Path:
    """Derived product: the five Orsi fronts as NaN-separated lon/lat polylines in one ``.npz``.

    Keys ``<front>_lon``, ``<front>_lat`` (float32, degrees, lon in −180…180) for ``stf saf pf saccf sbdy``, and
    ``attrs`` (JSON provenance, incl. the long names).

    Book: Vallis §21 Fig. 21.23; no equation. Units: [degrees].
    Validation: tests/data/ check PF ≈ 50°S and SAF ≈ 45°S south of Africa (20°E) when the product is present.
    """
    arrays = {}
    for f in FRONTS:
        lon, lat = read_orsi_txt(fetch_front(f))
        arrays[f"{f}_lon"], arrays[f"{f}_lat"] = lon, lat
    attrs = provenance(KEY, RAW + "{stf,saf,pf,saccf,sbdy}.txt", LICENCE, CITATION,
                       "original Orsi files parsed: '%' lines split segments (NaN separators); lon wrapped to "
                       "-180..180; float32", names=NAMES, attribution=ATTRIBUTION,
                       fetched=_dt.date.today().isoformat())
    out = Path(out) if out else data_root() / "derived" / "ch21" / "fig_21_23_fronts_orsi.npz"
    return write_npz(out, arrays, attrs)
