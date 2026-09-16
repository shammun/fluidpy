"""Structured Cartesian grids with ONE layout for the whole project.

Book: Ch. 2 §2.9 (the ∇ operator (2.22) acting on fields sampled at points); reused by every field computation from
Ch. 3 on.

**Layout (never deviate — a deviation silently transposes curls):**

* 3-D arrays are indexed ``[k, j, i]`` = (z, y, x): **x varies along the last axis**, y along the second-to-last,
  z along the first. Built with ``np.meshgrid(z, y, x, indexing="ij")``.
* 2-D arrays are ``[j, i]`` = (y, x), built with ``np.meshgrid(x, y, indexing="xy")`` (``u[j, i] = u(y_j, x_i)``).
* Vector fields stack the component on axis 0: ``u[c, k, j, i]`` (3-D) or ``u[c, j, i]`` (2-D); second-order tensor
  fields use two leading axes ``T[a, b, ...]``.
* "Direction" d ∈ {0, 1, 2} means (x, y, z) = (x₁, x₂, x₃); the array axis that varies with direction d is
  ``ndim − 1 − d`` (:func:`axis_of_direction`).
* Values live at the **nodes** ``x_i = x0 + i h`` (no half-cell offsets); ``periodic=True`` drops the duplicate end
  node so that ``h = L/n`` and the wrap-around stencil is consistent.

Units: coordinates in metres (or the non-dimensional length the chapter uses; the grid is unit-agnostic).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

__all__ = ["Grid2D", "Grid3D", "grid", "grid2d", "axis_of_direction", "spacing", "evaluate_field"]


def _bounds(bounds, dim: int) -> list[tuple[float, float]]:
    b = np.asarray(bounds, dtype=float)
    if b.shape == (2,):
        return [(float(b[0]), float(b[1]))] * dim
    if b.shape == (dim, 2):
        return [(float(lo), float(hi)) for lo, hi in b]
    raise ValueError(f"bounds must be (lo, hi) or {dim} pairs ((x0, x1), (y0, y1)[, (z0, z1)]); got shape {b.shape}")


def _counts(n, dim: int) -> list[int]:
    if np.ndim(n) == 0:
        return [int(n)] * dim
    if len(n) == dim:
        return [int(v) for v in n]
    raise ValueError(f"n must be an int or {dim} ints (nx, ny[, nz])")


def _axis1d(lo: float, hi: float, n: int, periodic: bool) -> tuple[np.ndarray, float]:
    if n < 3:
        raise ValueError("need at least 3 nodes per direction for second-order stencils")
    if periodic:
        h = (hi - lo) / n
        return lo + h * np.arange(n), h  # nodes lo … hi − h; the node at hi duplicates lo
    x, h = np.linspace(lo, hi, n, retstep=True)
    return x, float(h)


def axis_of_direction(ndim: int, direction: int) -> int:
    """Array axis along which coordinate ``direction`` (0 = x, 1 = y, 2 = z) varies: ``ndim − 1 − direction``."""
    if not 0 <= direction < ndim:
        raise ValueError(f"direction {direction} out of range for a {ndim}-D field")
    return ndim - 1 - direction


def spacing(h, direction: int) -> float:
    """Grid spacing for one direction from a scalar ``h`` or a per-direction sequence ``(hx, hy[, hz])`` [m]."""
    if np.ndim(h) == 0:
        return float(h)
    return float(np.asarray(h, dtype=float)[direction])


@dataclass(frozen=True)
class Grid2D:
    """Nodes of a 2-D grid: ``X[j, i] = x[i]``, ``Y[j, i] = y[j]``; ``h = (hx, hy)`` [m]."""

    x: np.ndarray
    y: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    h: tuple[float, float]
    periodic: bool = False

    @property
    def shape(self) -> tuple[int, int]:
        return self.X.shape  # (ny, nx)

    @property
    def coords(self) -> tuple[np.ndarray, np.ndarray]:
        """``(X, Y)`` — pass as ``fn(*grid.coords)`` to a field callable."""
        return self.X, self.Y


@dataclass(frozen=True)
class Grid3D:
    """Nodes of a 3-D grid: ``X[k, j, i] = x[i]``, ``Y[k, j, i] = y[j]``, ``Z[k, j, i] = z[k]``; ``h = (hx, hy, hz)`` [m]."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    Z: np.ndarray
    h: tuple[float, float, float]
    periodic: bool = False

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.X.shape  # (nz, ny, nx)

    @property
    def coords(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(X, Y, Z)`` — pass as ``fn(*grid.coords)`` to a field callable."""
        return self.X, self.Y, self.Z


def grid(bounds=(-1.0, 1.0), n: int | Sequence[int] = 24, periodic: bool = False) -> Grid3D:
    """Uniform 3-D node grid in the project layout ``[k, j, i]`` = (z, y, x).

    Book: §2.9 tool (fields for (2.22)–(2.25)); §2.12 (box volumes for (2.30)).

    Parameters
    ----------
    bounds : (lo, hi) or ((x0, x1), (y0, y1), (z0, z1))
        Extent in metres (the same pair for all three directions if only one is given).
    n : int or (nx, ny, nz)
        Nodes per direction (≥ 3). FAST notebooks use 24, full runs 48.
    periodic : bool
        Drop the duplicate end node so that ``h = L/n`` and ``operators.partial(..., bc="periodic")`` is exact.

    Returns
    -------
    Grid3D with ``x, y, z`` (1-D), ``X, Y, Z`` (3-D, shape (nz, ny, nx)) and ``h = (hx, hy, hz)``.

    Validation: V1 ``X`` varies only along the last axis, ``Z`` only along the first; ``h == (hi − lo)/(n − 1)``.
    Label: analytic.
    """
    b = _bounds(bounds, 3)
    c = _counts(n, 3)
    x, hx = _axis1d(*b[0], c[0], periodic)
    y, hy = _axis1d(*b[1], c[1], periodic)
    z, hz = _axis1d(*b[2], c[2], periodic)
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")  # [k, j, i] = (z, y, x): x on the LAST axis
    return Grid3D(x, y, z, X, Y, Z, (hx, hy, hz), periodic)


def grid2d(bounds=(-1.0, 1.0), n: int | Sequence[int] = 48, periodic: bool = False) -> Grid2D:
    """Uniform 2-D node grid ``[j, i]`` = (y, x) with ``np.meshgrid(x, y, indexing="xy")``.

    Book: §2.9 (plane fields such as Example 2.3's b × x with b = b e₃), §2.13 (Stokes on a plane).

    Parameters
    ----------
    bounds : (lo, hi) or ((x0, x1), (y0, y1))  [m]
    n : int or (nx, ny)  — nodes per direction (≥ 3)
    periodic : bool

    Returns
    -------
    Grid2D with ``x, y`` (1-D), ``X, Y`` (shape (ny, nx)) and ``h = (hx, hy)``.

    Validation: V1 ``X[j, i] == x[i]``, ``Y[j, i] == y[j]`` (``u[j, i] = u(y_j, x_i)``, math-to-python §3). Label: analytic.
    """
    b = _bounds(bounds, 2)
    c = _counts(n, 2)
    x, hx = _axis1d(*b[0], c[0], periodic)
    y, hy = _axis1d(*b[1], c[1], periodic)
    X, Y = np.meshgrid(x, y, indexing="xy")  # [j, i] = (y, x)
    return Grid2D(x, y, X, Y, (hx, hy), periodic)


def evaluate_field(fn: Callable, g: Grid2D | Grid3D) -> np.ndarray:
    """Sample a callable field ``fn(X, Y[, Z]) -> scalar array or (ncomp, ...) array`` on the grid nodes.

    Constant components (e.g. a zero third component) are broadcast to the grid shape so the result is always a full
    array: scalar field → ``g.shape``; vector field → ``(ncomp, *g.shape)``.
    """
    out = fn(*g.coords)
    if isinstance(out, (list, tuple)) or (isinstance(out, np.ndarray) and out.ndim == len(g.shape) + 1):
        return np.stack([np.broadcast_to(np.asarray(c, dtype=float), g.shape) for c in out])
    return np.broadcast_to(np.asarray(out, dtype=float), g.shape).copy()
