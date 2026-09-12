"""Project configuration: where the repository root is and what the public URLs are.

Everything URL-shaped (GitHub repo, GitHub Pages site, raw file links, the Colab badge) is derived from the
``project:`` block of ``book.yaml``, so renaming the repository is a one-line change::

    project:
      package: fluidpy
      repo: shammun/fluidpy        # GitHub user/repo
      branch: main
      site_url: https://shammun.github.io/fluidpy
      site_title: "Fluid Mechanics (Kundu, Cohen & Dowling, 5th ed.) - worked in Python"

Works from a local checkout, from a notebook started in ``notebooks/``, from headless ``nbconvert`` and from Colab
after the repository has been cloned into ``/content/fluidpy``.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

DEFAULTS: dict[str, Any] = {
    "project": {
        "package": "fluidpy",
        "repo": "shammun/fluidpy",
        "branch": "main",
        "site_url": "",
        "site_title": "fluidpy",
    },
    "book": {"title": "", "authors": []},
    "chapters": [],
}

_MARKERS = ("book.yaml",)


def _is_root(p: Path) -> bool:
    return (p / "book.yaml").is_file() or ((p / "fluidpy").is_dir() and (p / "tools").is_dir())


def repo_root(start: str | Path | None = None) -> Path:
    """Return the repository root: the nearest folder (walking up) holding ``book.yaml`` (or ``fluidpy/`` + ``tools/``).

    Search order: ``start`` (if given) → the current working directory → the location of this file →
    ``/content/fluidpy`` (Colab clone). Falls back to the parent of the ``fluidpy`` package.
    """
    bases: list[Path] = []
    if start is not None:
        bases.append(Path(start).resolve())
    bases.append(Path.cwd().resolve())
    bases.append(Path(__file__).resolve().parent)
    for base in bases:
        for cand in (base, *base.parents):
            if _is_root(cand):
                return cand
    colab = Path("/content/fluidpy")
    if _is_root(colab):
        return colab
    return Path(__file__).resolve().parents[2]


def _merge(defaults: dict, data: dict) -> dict:
    out = dict(defaults)
    for k, v in (data or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


@lru_cache(maxsize=8)
def _config_cached(root: str) -> dict:
    path = Path(root) / "book.yaml"
    data: dict = {}
    if path.is_file():
        import yaml

        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    cfg = _merge(DEFAULTS, data)
    proj = cfg["project"]
    user, _, name = str(proj["repo"]).partition("/")
    if not proj.get("site_url"):
        proj["site_url"] = f"https://{user}.github.io/{name}"
    proj["site_url"] = str(proj["site_url"]).rstrip("/")
    return cfg


def config(root: str | Path | None = None) -> dict:
    """The parsed ``book.yaml`` merged over sensible defaults (cached per root)."""
    return _config_cached(str(Path(root).resolve() if root else repo_root()))


def _proj(root: str | Path | None = None) -> dict:
    return config(root)["project"]


def repo_url(root: str | Path | None = None) -> str:
    return f"https://github.com/{_proj(root)['repo']}"


def pages_url(path: str = "", root: str | Path | None = None) -> str:
    """Public GitHub Pages URL of a repo-relative ``path`` (e.g. ``viz/ch07/dispersion.html``)."""
    base = _proj(root)["site_url"]
    return f"{base}/{path.lstrip('/')}" if path else f"{base}/"


def raw_url(path: str, root: str | Path | None = None) -> str:
    """raw.githubusercontent.com URL of a repo-relative ``path`` on the configured branch."""
    p = _proj(root)
    return f"https://raw.githubusercontent.com/{p['repo']}/{p['branch']}/{path.lstrip('/')}"


def github_blob_url(path: str, root: str | Path | None = None) -> str:
    p = _proj(root)
    return f"https://github.com/{p['repo']}/blob/{p['branch']}/{path.lstrip('/')}"


def colab_url(notebook_name: str, root: str | Path | None = None) -> str:
    """Open-in-Colab URL for ``notebooks/<notebook_name>`` (``.ipynb`` added if missing)."""
    p = _proj(root)
    name = notebook_name if notebook_name.endswith(".ipynb") else f"{notebook_name}.ipynb"
    if "/" not in name:
        name = f"notebooks/{name}"
    return f"https://colab.research.google.com/github/{p['repo']}/blob/{p['branch']}/{name}"


def chapters(root: str | Path | None = None) -> list[dict]:
    return list(config(root).get("chapters") or [])


def chapter(cid: str, root: str | Path | None = None) -> dict | None:
    return next((c for c in chapters(root) if c.get("id") == cid), None)
