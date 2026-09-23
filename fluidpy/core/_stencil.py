"""Private point-wise finite-difference helpers for field callables (not public API).

Convention (``core.kinematics``): a field is a callable ``F(x, t)`` with ``x`` of shape ``(d,)`` (one point) or
``(d, N)`` (N points, coordinates on axis 0). A scalar field returns ``()`` / ``(N,)``, a vector field ``(d,)`` /
``(d, N)``, a tensor field ``(d, d)`` / ``(d, d, N)``: component axes first, the point axis last. Constant fields may
return a plain number or a constant vector; :func:`ev` broadcasts them.

All stencils are explicit second-order central differences (never ``np.gradient``, which is first order at edges);
``h`` is a length [m], ``ht`` a time [s] — never the same default.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

#: Default time step [s] for ∂/∂t (independent of the length step h), as ``core.kinematics.DEFAULT_HT``.
DEFAULT_HT = 1e-4

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731


def ev(F: Callable, x, t, lead: tuple = ()) -> np.ndarray:
    """Evaluate F at x, t and broadcast the result to ``lead + x.shape[1:]`` (lead = () scalar, (d,) vector …)."""
    x_ = _F(x)
    val = _F(F(x_, t))
    target = tuple(lead) + x_.shape[1:]
    if val.shape != target:
        val = np.broadcast_to(val.reshape(val.shape + (1,) * (len(target) - val.ndim)), target).copy()
    return val


def unit(d: int, i: int, ndim: int) -> np.ndarray:
    """Unit offset e_i shaped (d,) or (d, 1, …) to shift a point array of ``ndim`` dimensions."""
    e = np.zeros(d)
    e[i] = 1.0
    return e.reshape((d,) + (1,) * (ndim - 1))


def ddx(F: Callable, x, t, i: int, h: float, lead: tuple = ()) -> np.ndarray:
    """∂F/∂x_i by the central difference (F(x + h e_i) − F(x − h e_i))/(2h); error O(h²)."""
    x_ = _F(x)
    e = unit(x_.shape[0], i, x_.ndim) * h
    return (ev(F, x_ + e, t, lead) - ev(F, x_ - e, t, lead)) / (2.0 * h)


def d2dx2(F: Callable, x, t, i: int, h: float, lead: tuple = ()) -> np.ndarray:
    """∂²F/∂x_i² by (F(x + h e_i) − 2F(x) + F(x − h e_i))/h²; error O(h²)."""
    x_ = _F(x)
    e = unit(x_.shape[0], i, x_.ndim) * h
    return (ev(F, x_ + e, t, lead) - 2.0 * ev(F, x_, t, lead) + ev(F, x_ - e, t, lead)) / h ** 2


def ddt(F: Callable, x, t, ht: float | None = None, lead: tuple = ()) -> np.ndarray:
    """∂F/∂t at fixed x by (F(x, t + ht) − F(x, t − ht))/(2ht); ``ht`` None → ``DEFAULT_HT``."""
    ht = DEFAULT_HT if ht is None else float(ht)
    return (ev(F, x, t + ht, lead) - ev(F, x, t - ht, lead)) / (2.0 * ht)


def grad(F: Callable, x, t, h: float) -> np.ndarray:
    """Gradient of a scalar field, shape (d,) + points."""
    x_ = _F(x)
    return np.stack([ddx(F, x_, t, i, h) for i in range(x_.shape[0])])


def div(F: Callable, x, t, h: float) -> np.ndarray:
    """Divergence ∂F_i/∂x_i of a vector field, shape points."""
    x_ = _F(x)
    d = x_.shape[0]
    return sum(ddx(lambda X, T, i=i: ev(F, X, T, (d,))[i], x_, t, i, h) for i in range(d))


def laplacian(F: Callable, x, t, h: float, lead: tuple = ()) -> np.ndarray:
    """∇²F = Σ_i ∂²F/∂x_i² (component-wise for vector fields), shape lead + points."""
    x_ = _F(x)
    return sum(d2dx2(F, x_, t, i, h, lead) for i in range(x_.shape[0]))


def grad_vector(F: Callable, x, t, h: float) -> np.ndarray:
    """G[i, j] = ∂F_i/∂x_j of a vector field (column j = derivative along x_j), shape (d, d) + points."""
    x_ = _F(x)
    d = x_.shape[0]
    return np.stack([ddx(F, x_, t, j, h, (d,)) for j in range(d)], axis=1)


def pad3(v: np.ndarray) -> np.ndarray:
    """Embed a 2-component vector array (2, …) into 3 components (third = 0); 3-component arrays are returned."""
    v = _F(v)
    if v.shape[0] == 3:
        return v
    return np.concatenate([v, np.zeros((1,) + v.shape[1:])], axis=0)


def gvec(g, d: int, pts_shape: tuple = ()) -> np.ndarray:
    """Body-force vector g [m/s²] reduced to the first d components and broadcast over the points."""
    g_ = _F(g).ravel()
    if g_.size < d:
        g_ = np.concatenate([g_, np.zeros(d - g_.size)])
    g_ = g_[:d]
    return np.broadcast_to(g_.reshape((d,) + (1,) * len(pts_shape)), (d,) + tuple(pts_shape)).copy()


def as_field(value) -> Callable:
    """A callable (x, t) → value for a constant, or the callable itself."""
    if callable(value):
        return value
    v = float(value)
    return lambda x, t: v
