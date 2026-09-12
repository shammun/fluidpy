"""Comparison helpers for the validation suite: fields, profiles, benchmark tables and conservation residuals.

Import in tests::

    from tools.compare_fields import load_ref, assert_close, compare_profiles, assert_conserved, rel_error

These return a metrics dict and raise ``AssertionError`` carrying those metrics, so a failing test tells the report what
the error actually was instead of just "False is not True".
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def load_ref(path: str | Path, key: str | None = None):
    """Load reference data saved by ``reference/chNN/make_refs.py``: .npz (key required), .npy, .csv or .mat."""
    p = Path(path)
    if p.suffix == ".npz":
        d = np.load(p)
        if key is None:
            raise ValueError(f"{p} is an .npz; give a key (available: {list(d.files)})")
        if key not in d:
            raise KeyError(f"{key} not in {p}; keys = {list(d.files)}")
        return d[key]
    if p.suffix == ".npy":
        return np.load(p)
    if p.suffix in (".csv", ".dat", ".txt"):
        import pandas as pd

        df = pd.read_csv(p)
        return df[key].to_numpy() if key else df
    if p.suffix == ".mat":
        from scipy.io import loadmat

        d = loadmat(str(p))
        if key not in d:
            raise KeyError(f"{key} not in {p}; keys = {[k for k in d if not k.startswith('__')]}")
        return np.asarray(d[key])
    raise ValueError(f"unsupported reference format: {p.suffix}")


def rel_error(num, ref, floor: float = 0.0) -> np.ndarray:
    """Element-wise relative error, using ``max(|ref|, floor)`` as the denominator (floor protects near-zero values)."""
    num = np.asarray(num, dtype=float)
    ref = np.asarray(ref, dtype=float)
    den = np.maximum(np.abs(ref), floor if floor else np.max(np.abs(ref)) * 1e-12)
    return np.abs(num - ref) / den


def assert_close(num, ref, rtol: float = 1e-6, atol: float = 0.0, name: str = "", scale: float | None = None) -> dict:
    """Compare two fields of the same shape; report max/mean absolute and relative error.

    ``scale`` (e.g. the free-stream velocity) makes the relative error meaningful for fields that cross zero.
    """
    num = np.asarray(num, dtype=float)
    ref = np.asarray(ref, dtype=float)
    if num.shape != ref.shape:
        if num.squeeze().shape == ref.squeeze().shape:
            num, ref = num.squeeze(), ref.squeeze()
        else:
            raise AssertionError(f"{name}: shape mismatch {num.shape} vs {ref.shape}")
    err = np.abs(num - ref)
    den = scale if scale is not None else max(float(np.max(np.abs(ref))), 1e-300)
    metrics = {
        "max_abs_err": float(np.nanmax(err)),
        "mean_abs_err": float(np.nanmean(err)),
        "max_rel_err": float(np.nanmax(err) / den),
        "l2_err": float(np.sqrt(np.nanmean(err ** 2))),
    }
    if not np.allclose(num, ref, rtol=rtol, atol=atol, equal_nan=True):
        raise AssertionError(f"{name}: {metrics} (rtol={rtol}, atol={atol})")
    return metrics


def compare_profiles(x_num, y_num, x_ref, y_ref, rtol: float = 0.01, name: str = "", scale: float | None = None) -> dict:
    """Compare a computed profile with a benchmark table sampled at different abscissae.

    The computed profile is interpolated onto the reference abscissae (inside their range only), which is what you want
    when checking e.g. Ghia's cavity centreline values against your own grid.
    """
    x_num = np.asarray(x_num, float)
    y_num = np.asarray(y_num, float)
    x_ref = np.asarray(x_ref, float)
    y_ref = np.asarray(y_ref, float)
    order = np.argsort(x_num)
    inside = (x_ref >= x_num[order][0]) & (x_ref <= x_num[order][-1])
    if inside.sum() == 0 or inside.sum() < 0.5 * x_ref.size:
        raise AssertionError(
            f"{name}: the computed profile covers only {int(inside.sum())}/{x_ref.size} reference points "
            f"(computed range {x_num[order][0]:.4g}..{x_num[order][-1]:.4g}, reference range "
            f"{x_ref.min():.4g}..{x_ref.max():.4g}) — extend the domain or check the abscissa definition"
        )
    y_i = np.interp(x_ref[inside], x_num[order], y_num[order])
    metrics = assert_close(y_i, y_ref[inside], rtol=rtol, atol=0.0, name=name or "profile", scale=scale)
    metrics["n_compared"] = int(inside.sum())
    metrics["n_reference"] = int(x_ref.size)
    return metrics


def assert_conserved(series, tol: float = 1e-12, name: str = "invariant", relative: bool = True) -> dict:
    """Assert that a quantity recorded over time steps does not drift (mass, energy, circulation, ...)."""
    s = np.asarray(series, dtype=float)
    ref = s[0]
    drift = np.abs(s - ref)
    if relative and abs(ref) > 0:
        drift = drift / abs(ref)
    metrics = {"initial": float(ref), "final": float(s[-1]), "max_drift": float(np.max(drift)),
               "relative": relative, "n": int(s.size)}
    if metrics["max_drift"] > tol:
        raise AssertionError(f"{name} not conserved: {metrics} (tol={tol})")
    return metrics
