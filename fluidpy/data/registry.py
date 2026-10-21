"""Read the AOFD dataset registry (``data/registry.yaml`` in the AOFD repo) when it is available.

The registry is the single list of datasets the book uses: URL, access route, licence, citation, size, sha256,
redistributable flag, raw files and derived products. On Colab (fluidpy installed without the AOFD repo) the registry
is empty and :func:`fluidpy.data.fetch` works with plain URLs.
"""
from __future__ import annotations

import os
from pathlib import Path

__all__ = ["registry_path", "load", "get", "evictable"]

_CACHE: dict | None = None
_MTIME: float | None = None


def registry_path() -> Path | None:
    """Path of ``data/registry.yaml`` (env ``AOFD_REGISTRY`` overrides), or None when not found.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: none [dimensionless].
    Validation: smoke-tested in tests/data/ (no network).
    """
    env = os.environ.get("AOFD_REGISTRY")
    if env:
        return Path(env) if Path(env).is_file() else None
    from . import aofd_repo

    repo = aofd_repo()
    p = repo / "data" / "registry.yaml" if repo else None
    return p if p is not None and p.is_file() else None


def load() -> dict:
    """All registry entries as ``{key: entry}`` (``{}`` when there is no registry); re-read when the file changes.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: sizes in [bytes].
    Validation: smoke-tested in tests/data/ (no network).
    """
    global _CACHE, _MTIME
    p = registry_path()
    if p is None:
        return {}
    m = p.stat().st_mtime
    if _CACHE is None or m != _MTIME:
        import yaml

        d = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
        _CACHE, _MTIME = (d.get("datasets") or {}), m
    return _CACHE


def get(key: str) -> dict | None:
    """One registry entry, or None.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: sizes in [bytes].
    Validation: smoke-tested in tests/data/ (no network).
    """
    return load().get(key)


def evictable() -> dict:
    """Raw files the disk guard may evict: ``{relative path under raw/: {"redownloadable", "derived"}}``.

    Only entries with ``redownloadable: true`` and a ``derived_on_root`` product (path relative to the data root)
    qualify; ``tools/disk_guard.evict`` additionally checks that the derived product exists.

    Book: data access for Vallis-book figures (no equation; not §-specific). Units: none [dimensionless].
    Validation: smoke-tested in tests/data/ (no network).
    """
    out = {}
    for key, ent in load().items():
        if not ent or not ent.get("redownloadable", True):
            continue
        ders = ent.get("derived_on_root") or []
        if not ders:
            continue
        for rf in ent.get("raw_files") or []:
            rel = rf[4:] if rf.startswith("raw/") else rf
            out[rel] = {"redownloadable": True, "derived": ders[0], "key": key}
    return out
