"""Vector and scalar fields defined by sympy expressions, callable on numpy grids, carrying their exact derivatives.

Book: Ch. 2 §2.9 (Example 2.3: div and curl of a x and b × x), §2.12–2.13 (fields Q, u whose exact ∇·Q, ∇×u the
integral theorems are checked against). Reused by every later chapter that needs a test field with a known answer
(Tier-1 analytic evidence, ``data-and-benchmarks``).

A :class:`VectorField` is a callable ``u(X, Y[, Z]) -> array (ncomp, ...)`` (component on axis 0, broadcast to the
grid shape — the project layout of ``core.grids``) with attributes ``.exprs`` (sympy), ``.div_expr``,
``.curl_expr``, ``.grad_exprs`` and their lambdified twins ``.div_fn``, ``.curl_fn``, ``.grad_fn`` (same calling
convention). Optional ``params`` substitute numeric values for free symbols in the numeric callables while the
symbolic expressions keep the symbols (Example 2.3 with symbolic a, b). ``singular_at`` marks a point where the field
is not differentiable (the irrotational vortex core) so Stokes' theorem can report the hypothesis failure.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
import sympy as sp

from .index_notation import coordinates

__all__ = ["VectorField", "ScalarField"]


def _broadcast_stack(vals, shape) -> np.ndarray:
    return np.stack([np.broadcast_to(np.asarray(v, dtype=float), shape) for v in vals])


class ScalarField:
    """φ(x) from a sympy expression: callable on grids, with exact gradient ``.grad_exprs`` / ``.grad_fn``.

    Book: §2.9 (gradient of a scalar, Fig. 2.7); §2.12 (scalar Q in Gauss' theorem (2.30)).
    """

    def __init__(self, expr, coords: Sequence[sp.Symbol] | None = None, params: dict | None = None, name: str = "phi"):
        self.coords = tuple(coords) if coords is not None else coordinates(3)
        self.dim = len(self.coords)
        self.expr = sp.sympify(expr)
        self.params = dict(params or {})
        self.name = name
        self.grad_exprs = [sp.diff(self.expr, c) for c in self.coords]  # (∇φ)_i = ∂φ/∂x_i
        sub = lambda e: e.subs(self.params)  # noqa: E731
        self._f = sp.lambdify(self.coords, sub(self.expr), "numpy")
        self._g = sp.lambdify(self.coords, [sub(e) for e in self.grad_exprs], "numpy")

    def __call__(self, *coords):
        B = np.broadcast_arrays(*[np.asarray(c, dtype=float) for c in coords])
        return np.broadcast_to(np.asarray(self._f(*B), dtype=float), B[0].shape).copy()

    def grad_fn(self, *coords) -> np.ndarray:
        """Exact gradient, shape ``(dim, ...)``."""
        B = np.broadcast_arrays(*[np.asarray(c, dtype=float) for c in coords])
        return _broadcast_stack(self._g(*B), B[0].shape)

    def __repr__(self) -> str:
        return f"ScalarField({self.name} = {self.expr})"


class VectorField:
    """u(x) from sympy component expressions: callable on grids, with exact divergence, curl and gradient.

    Parameters
    ----------
    exprs : sequence of sympy expressions (length 2 or 3) in ``coords``
    coords : coordinate symbols (default ``coordinates(len(exprs))`` = x1, x2[, x3])
    params : dict {symbol: value} substituted in the numeric callables only
    name : str
    singular_at : point (tuple) where the field is singular, or None

    Book: §2.9 Eqs. (2.23)–(2.25) give ``div_expr`` and ``curl_expr``; the velocity gradient ``grad_exprs[i][j]``
    = ∂u_i/∂x_j (index order of ``core.operators.vector_gradient``). For a 2-component plane field ``curl_expr`` is
    the scalar (∇×u)₃ = ∂u₂/∂x₁ − ∂u₁/∂x₂.

    Validation: V2 Example 2.3: div(a x) = 3a, curl(a x) = 0, div(b × x) = 0, curl(b × x) = 2b symbolically; V1 numeric
    twins agree with ``core.operators`` to truncation. Label: symbolic.
    """

    def __init__(self, exprs, coords: Sequence[sp.Symbol] | None = None, params: dict | None = None,
                 name: str = "u", singular_at=None):
        exprs = [sp.sympify(e) for e in exprs]
        self.coords = tuple(coords) if coords is not None else coordinates(len(exprs))
        self.dim = len(self.coords)
        if len(exprs) != self.dim:
            raise ValueError("VectorField needs one component per coordinate (2-D: 2, 3-D: 3)")
        self.exprs = exprs
        self.params = dict(params or {})
        self.name = name
        self.singular_at = None if singular_at is None else tuple(float(v) for v in singular_at)
        X = self.coords
        self.div_expr = sp.simplify(sum(sp.diff(e, x) for e, x in zip(exprs, X)))  # Eq. (2.23): ∂u_i/∂x_i
        self.grad_exprs = [[sp.diff(ui, xj) for xj in X] for ui in exprs]  # G_ij = ∂u_i/∂x_j
        if self.dim == 3:
            u1, u2, u3 = exprs
            x1, x2, x3 = X
            self.curl_expr = [sp.simplify(sp.diff(u3, x2) - sp.diff(u2, x3)),  # Eq. (2.25)
                              sp.simplify(sp.diff(u1, x3) - sp.diff(u3, x1)),
                              sp.simplify(sp.diff(u2, x1) - sp.diff(u1, x2))]
        else:
            u1, u2 = exprs
            x1, x2 = X
            self.curl_expr = sp.simplify(sp.diff(u2, x1) - sp.diff(u1, x2))  # (∇×u)₃ for a plane field
        sub = lambda e: e.subs(self.params)  # noqa: E731
        self._f = sp.lambdify(X, [sub(e) for e in exprs], "numpy")
        self._div = sp.lambdify(X, sub(self.div_expr), "numpy")
        curl_list = self.curl_expr if self.dim == 3 else [self.curl_expr]
        self._curl = sp.lambdify(X, [sub(e) for e in curl_list], "numpy")
        self._grad = sp.lambdify(X, [[sub(e) for e in row] for row in self.grad_exprs], "numpy")

    # numeric twins ---------------------------------------------------------------------------------------------------
    def _bc(self, coords):
        if len(coords) != self.dim:
            raise ValueError(f"{self.name} is a {self.dim}-D field; got {len(coords)} coordinate arrays")
        return np.broadcast_arrays(*[np.asarray(c, dtype=float) for c in coords])

    def __call__(self, *coords) -> np.ndarray:
        B = self._bc(coords)
        return _broadcast_stack(self._f(*B), B[0].shape)

    def div_fn(self, *coords) -> np.ndarray:
        """Exact ∇·u on the given coordinates (same shape as one coordinate array)."""
        B = self._bc(coords)
        return np.broadcast_to(np.asarray(self._div(*B), dtype=float), B[0].shape).copy()

    def curl_fn(self, *coords) -> np.ndarray:
        """Exact ∇×u: shape ``(3, ...)`` in 3-D, the scalar (∇×u)₃ array in 2-D."""
        B = self._bc(coords)
        out = _broadcast_stack(self._curl(*B), B[0].shape)
        return out if self.dim == 3 else out[0]

    def grad_fn(self, *coords) -> np.ndarray:
        """Exact velocity gradient G[i, j, ...] = ∂u_i/∂x_j."""
        B = self._bc(coords)
        rows = self._grad(*B)
        return np.stack([_broadcast_stack(row, B[0].shape) for row in rows])

    def __repr__(self) -> str:
        return f"VectorField({self.name} = {self.exprs}; div = {self.div_expr}, curl = {self.curl_expr})"
