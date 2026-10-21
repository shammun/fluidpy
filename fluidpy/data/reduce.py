"""Write small derived products (float32, zlib) with provenance attributes.

A derived product is what notebooks, figures and scenes load by default (no keys, seconds on Colab). Each one carries
its dataset key, source URL, licence, citation, processing description, fetch date and the fluidpy git commit.
"""
from __future__ import annotations

import datetime as _dt
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

__all__ = ["provenance", "write_nc", "write_npz", "git_commit"]


def git_commit() -> str:
    """The fluidpy repo's HEAD commit (``+dirty`` when the work tree has changes), or ``"unknown"``.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: none [dimensionless].
    Validation: smoke-tested in tests/data/ (no network).
    """
    repo = Path(__file__).resolve().parents[2]
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=repo, capture_output=True, text=True,
                             timeout=20).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "fluidpy"], cwd=repo, capture_output=True,
                               text=True, timeout=20).stdout.strip()
        return (sha or "unknown") + ("+dirty" if dirty else "")
    except Exception:
        return "unknown"


def provenance(key: str, url: str, licence: str, citation: str, processing: str, **extra: Any) -> dict:
    """Standard provenance attributes for a derived product.

    Parameters
    ----------
    key : str
        Dataset key in ``data/registry.yaml``.
    url, licence, citation, processing : str
        Source URL, licence (with URL), citation, and what was done to the raw data.
    **extra
        Further attributes (period, region, units notes...).

    Returns
    -------
    dict
        Attributes (strings only, netCDF-safe).

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: none [dimensionless].
    Validation: smoke-tested in tests/data/ (no network).
    """
    a = {"dataset_key": key, "source_url": url, "licence": licence, "citation": citation, "processing": processing,
         "created": _dt.datetime.now().isoformat(timespec="seconds"), "fluidpy_commit": git_commit(),
         "producer": "fluidpy.data (AOFD data-scout)"}
    a.update({k: (v if isinstance(v, str) else json.dumps(v)) for k, v in extra.items()})
    return a


def write_nc(ds, path: str | Path, attrs: dict | None = None, complevel: int = 6) -> Path:
    """Write an xarray Dataset as float32 + zlib netCDF with provenance ``attrs``; returns the path.

    Floating variables are cast to float32; integer masks are kept as they are.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: as in ``ds`` [SI].
    Validation: smoke-tested in tests/data/ (no network).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ds = ds.copy()
    for v in ds.data_vars:
        if np.issubdtype(ds[v].dtype, np.floating):
            ds[v] = ds[v].astype("float32")
    if attrs:
        ds.attrs.update(attrs)
    enc = {v: {"zlib": True, "complevel": complevel} for v in ds.data_vars}
    for c in ds.coords:
        if np.issubdtype(ds[c].dtype, np.floating):
            enc[c] = {"dtype": "float32"}
    tmp = path.with_name(path.name + ".part")
    ds.to_netcdf(tmp, encoding=enc)
    tmp.replace(path)
    return path


def write_npz(path: str | Path, arrays: dict, attrs: dict) -> Path:
    """Write arrays (float → float32) to a compressed ``.npz`` with the provenance as a JSON string ``attrs``.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: as given [SI / degrees].
    Validation: smoke-tested in tests/data/ (no network).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = {k: (np.asarray(v, dtype=np.float32) if np.asarray(v).dtype.kind == "f" else np.asarray(v))
           for k, v in arrays.items()}
    out["attrs"] = np.array(json.dumps(attrs))
    np.savez_compressed(path, **out)
    return path
