"""fluidpy.data — open-data access for the AOFD book: registry, fetchers, reductions and labelled pseudo-data.

Every download goes through one door, :func:`fetch`, which estimates the size, checks the disk guard, downloads into
``<data_root>/raw/<key>/…`` (pooch, hash pinned for fixed files), logs the transfer and returns the local path. Raw files
live under the data root (env ``AOFD_DATA_ROOT``; on Colab ``/content/aofd-data``; otherwise
``G:/Geophysical Fluid Dynamics/aofd-data``), never inside a repository. Remote OPeNDAP reads on Windows use
``engine="pydap"`` (netCDF4's bundled curl cannot find TLS certificates there), see :func:`open_remote`.

Disk guard
----------
When the AOFD repository is present (env ``AOFD_REPO``, or the sibling folder ``../AOFD Vallis``), its
``tools/disk_guard.py`` is used. Otherwise a minimal built-in guard applies the same rules: the data drive keeps
50 GB free and the project's data stays below 150 GB; C: keeps 40 GB free; a missing data drive is an error — the
guard never falls back to C:.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import importlib.util
import os
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Any

__all__ = ["fetch", "data_root", "open_remote", "disk_report", "DiskGuardError", "sha256", "aofd_repo"]

GB = 1024 ** 3
_DEFAULT_ROOT = "G:/Geophysical Fluid Dynamics/aofd-data"
_COLAB_ROOT = "/content/aofd-data"
_FLOORS = {"data_floor_gb": 50.0, "project_cap_gb": 150.0, "C_floor_gb": 40.0}


class DiskGuardError(RuntimeError):
    """A download would cross a disk floor or cap, or the data drive is missing."""


# --------------------------------------------------------------------------------------------- locations
def _on_colab() -> bool:
    return "google.colab" in sys.modules or bool(os.environ.get("COLAB_RELEASE_TAG"))


def data_root() -> Path:
    """The data root: env ``AOFD_DATA_ROOT`` → Colab ``/content/aofd-data`` → ``G:/Geophysical Fluid Dynamics/aofd-data``.

    Returns
    -------
    pathlib.Path
        The root (not created here; :func:`fetch` creates sub-folders once the guard has passed).

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: sizes in [bytes] / [GB].
    Validation: smoke-tested in tests/data/ (no network by default).
    """
    env = os.environ.get("AOFD_DATA_ROOT")
    if env:
        return Path(env)
    if _on_colab():
        return Path(_COLAB_ROOT)
    return Path(_DEFAULT_ROOT)


def aofd_repo() -> Path | None:
    """The AOFD book repository (env ``AOFD_REPO``, else the sibling ``../AOFD Vallis`` of the fluidpy repo), or None.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: sizes in [bytes] / [GB].
    Validation: smoke-tested in tests/data/ (no network by default).
    """
    cands = []
    if os.environ.get("AOFD_REPO"):
        cands.append(Path(os.environ["AOFD_REPO"]))
    here = Path(__file__).resolve()
    cands.append(here.parents[2].parent / "AOFD Vallis")
    for c in cands:
        if (c / "tools" / "disk_guard.py").is_file() or (c / "data" / "registry.yaml").is_file():
            return c
    return None


# --------------------------------------------------------------------------------------------- disk guard
_GUARD: Any = None


def _guard_module():
    """Import ``tools/disk_guard.py`` from the AOFD repo (cached); None when it is not available."""
    global _GUARD
    if _GUARD is not None:
        return _GUARD or None
    repo = aofd_repo()
    p = repo / "tools" / "disk_guard.py" if repo else None
    if p is not None and p.is_file():
        spec = importlib.util.spec_from_file_location("aofd_disk_guard", p)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)  # type: ignore[union-attr]
        _GUARD = mod
        return mod
    _GUARD = False
    return None


def _dir_size(p: Path) -> int:
    if not p.exists():
        return 0
    tot = 0
    for dp, _d, fs in os.walk(p):
        for f in fs:
            try:
                tot += os.path.getsize(os.path.join(dp, f))
            except OSError:
                pass
    return tot


def _builtin_status() -> dict:
    root = data_root()
    anchor = Path(root.anchor) if root.anchor else Path("/")
    present = anchor.exists()
    c_anchor = Path(os.environ.get("SystemDrive", "C:") + "\\") if os.name == "nt" else Path("/")
    s = {"data_root": str(root), "data_drive": (root.drive or root.anchor or "/").rstrip("\\/") or "/",
         "data_drive_present": present, "C_free_gb": round(shutil.disk_usage(c_anchor).free / GB, 1), **_FLOORS}
    if present:
        s["data_free_gb"] = round(shutil.disk_usage(anchor).free / GB, 1)
        s["data_used_gb"] = round(_dir_size(root) / GB, 2)
        s["raw_gb"] = round(_dir_size(root / "raw") / GB, 2)
        s["cache_gb"] = round(_dir_size(root / "cache") / GB, 2)
    return s


def _builtin_check(need_bytes: float, what: str) -> dict:
    s = _builtin_status()
    need = need_bytes / GB
    if not s["data_drive_present"]:
        raise DiskGuardError(f"data drive {s['data_drive']} not found ({s['data_root']}). Connect the drive or set "
                             "AOFD_DATA_ROOT; fluidpy never falls back to C:.")
    if s["data_free_gb"] - need < s["data_floor_gb"]:
        raise DiskGuardError(f"{what} needs {need:.2f} GB; {s['data_drive']} would drop below its "
                             f"{s['data_floor_gb']} GB floor ({s['data_free_gb']} GB free)")
    if s["data_used_gb"] + need > s["project_cap_gb"]:
        raise DiskGuardError(f"{what} needs {need:.2f} GB; project data would exceed {s['project_cap_gb']} GB")
    if s["data_drive"].upper().startswith("C") and s["C_free_gb"] - need < s["C_floor_gb"]:
        raise DiskGuardError(f"{what} needs {need:.2f} GB; C: would drop below its {s['C_floor_gb']} GB floor")
    return s


def disk_report() -> str:
    """One line: free space on C: and the data drive, project data size and cache size.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: sizes in [bytes] / [GB].
    Validation: smoke-tested in tests/data/ (no network by default).
    """
    g = _guard_module()
    if g is not None:
        return g.report()
    s = _builtin_status()
    if not s["data_drive_present"]:
        return f"C: {s['C_free_gb']} GB free | data drive {s['data_drive']} MISSING ({s['data_root']})"
    return (f"C: {s['C_free_gb']} GB free (floor {s['C_floor_gb']}) | {s['data_drive']} {s['data_free_gb']} GB free "
            f"(floor {s['data_floor_gb']}) | project data {s['data_used_gb']} GB (raw {s['raw_gb']}, "
            f"cache {s['cache_gb']})")


def _check(need_bytes: float, what: str) -> None:
    """Run the guard; on a refusal try LRU eviction of registered raw files once, then re-check."""
    g = _guard_module()
    if g is None:
        _builtin_check(need_bytes, what)
        return
    try:
        g.check(need_bytes=need_bytes, what=what)
    except g.DiskGuardError as e:
        if "not found" in str(e):
            raise DiskGuardError(str(e)) from None
        from . import registry as _reg

        evicted = g.evict(need_bytes, registry=_reg.evictable())
        if evicted:
            print(f"[fluidpy.data] evicted: {evicted}")
        try:
            g.check(need_bytes=need_bytes, what=what)
        except g.DiskGuardError as e2:
            raise DiskGuardError(f"{e2}. Nothing more can be evicted automatically; free space by hand.") from None


# --------------------------------------------------------------------------------------------- helpers
def sha256(path: str | Path) -> str:
    """Hex sha256 of a file.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: sizes in [bytes] / [GB].
    Validation: smoke-tested in tests/data/ (no network by default).
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _log(line: str) -> None:
    root = data_root()
    try:
        (root / "logs").mkdir(parents=True, exist_ok=True)
        with open(root / "logs" / "downloads.log", "a", encoding="utf-8") as f:
            f.write(line.rstrip() + "\n")
    except OSError:
        pass


def _is_opendap(url: str) -> bool:
    return "/dodsC/" in url or url.startswith(("dap2://", "dap4://"))


def _head_size(url: str) -> int | None:
    """Content-Length from an HTTP HEAD (redirects followed); None if the server does not say."""
    try:
        import requests

        r = requests.head(url, allow_redirects=True, timeout=30, headers={"User-Agent": "fluidpy-data/0.3"})
        n = r.headers.get("Content-Length")
        if r.status_code < 400 and n is not None:
            return int(n)
    except Exception:
        pass
    return None


def open_remote(url: str, **kw: Any):
    """Open a remote dataset lazily with xarray.

    OPeNDAP URLs (``…/dodsC/…``) use ``engine="pydap"`` on Windows, because netCDF4's bundled curl cannot verify TLS
    certificates there; elsewhere xarray's default engine is used unless ``engine`` is given.

    Parameters
    ----------
    url : str
        OPeNDAP or file URL.
    **kw
        Passed to :func:`xarray.open_dataset`.

    Returns
    -------
    xarray.Dataset
        Lazy dataset (only metadata is transferred until values are accessed).

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: sizes in [bytes] / [GB].
    Validation: smoke-tested in tests/data/ (no network by default).
    """
    import xarray as xr

    if _is_opendap(url) and os.name == "nt":
        kw.setdefault("engine", "pydap")
    return xr.open_dataset(url, **kw)


def _select(ds, subset: dict):
    sel = ds[subset["vars"]] if subset.get("vars") else ds
    if subset.get("sel"):
        sel = sel.sel(**subset["sel"])
    if subset.get("isel"):
        sel = sel.isel(**subset["isel"])
    return sel


def _subset_bytes(sel) -> int:
    return int(sum(v.size * v.dtype.itemsize for v in sel.variables.values()))


def _slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")[:120] or "file"


# --------------------------------------------------------------------------------------------- the one door
def fetch(key_or_url: str, *, subset: dict | str | None = None, expected_bytes: int | None = None,
          refresh: bool = False, key: str | None = None, fname: str | None = None,
          known_hash: str | None = None, quiet: bool = False) -> Path:
    """Download one file (or one server-side subset) into ``<data_root>/raw/<key>/`` through the disk guard.

    Parameters
    ----------
    key_or_url : str
        A key of ``data/registry.yaml`` (its ``url``/``raw_url``, ``sha256`` and ``size_bytes`` are used) or a URL.
    subset : dict or str, optional
        * ``dict`` for OPeNDAP URLs: ``{"vars": [...], "sel": {...}, "isel": {...}}``; the selection is read
          remotely (pydap on Windows) and written as a small netCDF file. Size = the selection's bytes.
        * ``str``: a query string appended to the URL (ERDDAP ``griddap`` constraints such as
          ``"U[(-72.875):1:(-30)][(-180):1:(180)]"``).
    expected_bytes : int, optional
        Size used by the guard when HEAD gives no Content-Length (registry ``size_bytes`` is used otherwise).
    refresh : bool
        Download again even when the file exists.
    key : str, optional
        Raw sub-folder name when ``key_or_url`` is a URL (default: the URL's host).
    fname : str, optional
        File name in the raw folder (default: the URL's last path segment, or a slug of the subset).
    known_hash : str, optional
        ``"sha256:<hex>"`` to pin; overrides the registry's ``sha256``. Not used for subsets.
    quiet : bool
        Do not print the disk report lines.

    Returns
    -------
    pathlib.Path
        Local path of the downloaded (or already cached) file.

    Raises
    ------
    DiskGuardError
        Data drive missing, a floor/cap would be crossed, or the size cannot be estimated.

    Notes
    -----
    Free space on C: and the data drive is printed before and after every transfer and written to
    ``<data_root>/logs/downloads.log`` with bytes, seconds and sha256.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: sizes in [bytes] / [GB].
    Validation: smoke-tested in tests/data/ (no network by default).
    """
    from . import registry as _reg

    if re.match(r"^[a-z][a-z0-9+.-]*://", key_or_url, re.I):
        url = key_or_url
        key = key or _slug(url.split("//", 1)[1].split("/", 1)[0])
    else:
        ent = _reg.get(key_or_url)
        if ent is None:
            raise KeyError(f"{key_or_url!r} is neither a URL nor a key of data/registry.yaml")
        url = ent.get("raw_url") or ent["url"]
        key = key or key_or_url
        if known_hash is None and ent.get("sha256") and subset is None:
            known_hash = f"sha256:{ent['sha256']}"
        if expected_bytes is None and ent.get("size_bytes"):
            expected_bytes = int(ent["size_bytes"])
        fname = fname or ent.get("raw_file")

    if isinstance(subset, str):
        q = subset.lstrip("?")
        fname = fname or _slug(url.rsplit("/", 1)[-1].split("?")[0] + "_" + q)
        url = f"{url}{'&' if '?' in url else '?'}{q}"
    elif isinstance(subset, dict):
        tag = "_".join(f"{k}{v}" for k, v in (subset.get("sel") or {}).items())
        fname = fname or _slug(url.rsplit("/", 1)[-1].rsplit(".", 1)[0] + "_" + tag) + ".nc"
    fname = fname or _slug(url.rsplit("/", 1)[-1].split("?")[0])
    root = data_root()
    out = root / "raw" / key / fname
    if out.exists() and not refresh:
        if not quiet:
            print(f"[fluidpy.data] cached {out}  | {disk_report()}")
        return out

    remote = remote_sel = None
    if isinstance(subset, dict):
        remote = open_remote(url)
        remote_sel = _select(remote, subset)
        need = _subset_bytes(remote_sel)
    else:
        need = None if _is_opendap(url) else _head_size(url)
        if need is None:
            need = expected_bytes
        if need is None:
            raise DiskGuardError(f"cannot estimate the size of {url}; pass expected_bytes= or register size_bytes")
    what = f"{key}/{fname} ({url[:160]})"
    before = disk_report()
    if not quiet:
        print(f"[fluidpy.data] before: {before}")
        print(f"[fluidpy.data] need ~ {need / 1e6:.3f} MB for {what}")
    _check(float(need), what)
    out.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    if isinstance(subset, dict):
        ds = remote_sel.load()
        # reserved/remote-only attributes (e.g. _NCProperties from pydap) cannot be written by netCDF4
        ds.attrs = {k: v for k, v in ds.attrs.items() if not k.startswith("_")}
        for v in ds.variables.values():
            v.attrs = {k: a for k, a in v.attrs.items() if not k.startswith("_")}
        enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
        ds.attrs.update({"fluidpy_source_url": url, "fluidpy_subset": repr(subset),
                         "fluidpy_fetched": _dt.datetime.now().isoformat(timespec="seconds")})
        tmp = out.with_name(out.name + ".part")
        ds.to_netcdf(tmp, encoding=enc)
        remote.close()
        os.replace(tmp, out)
    else:
        import pooch

        if out.exists():
            out.unlink()
        pooch.retrieve(url, known_hash=known_hash, fname=fname, path=out.parent, progressbar=False,
                       downloader=pooch.HTTPDownloader(headers={"User-Agent": "fluidpy-data/0.3"}, timeout=300))
    dt_s = time.perf_counter() - t0
    nbytes = out.stat().st_size
    digest = sha256(out)
    after = disk_report()
    _log(f"{_dt.datetime.now().isoformat(timespec='seconds')} key={key} file={fname} bytes={nbytes} "
         f"seconds={dt_s:.1f} sha256={digest} url={url} | before: {before} | after: {after}")
    if not quiet:
        print(f"[fluidpy.data] got {out} - {nbytes:,} bytes in {dt_s:.1f} s, sha256 {digest}")
        print(f"[fluidpy.data] after:  {after}")
    return out
