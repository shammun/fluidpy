"""Order-of-accuracy helpers for the V3 evidence level (see the verify-implementation skill).

Import in tests::

    from tools.convergence import observed_order, refinement_study, richardson, l2_error, mms

Why these exist: a discretisation written from a book's equations can be plausible, produce a pretty picture, and still
be first-order because a boundary condition is applied half a cell off. Measuring the observed order catches that; eyeing
the plot does not.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np


def l2_error(num: np.ndarray, exact: np.ndarray, h: float | None = None) -> float:
    """Discrete L2 norm of the error. With ``h`` given, the grid-normalised norm sqrt(h * sum(e^2))."""
    e = np.asarray(num, dtype=float) - np.asarray(exact, dtype=float)
    if h is None:
        return float(np.sqrt(np.mean(e ** 2)))
    return float(np.sqrt(h * np.sum(e ** 2)))


def observed_order(h: Sequence[float], err: Sequence[float]) -> float:
    """Least-squares slope of log(err) vs log(h): the observed order of accuracy p in err ~ C h^p.

    Use at least three grids. Errors at or below round-off (~1e-14) must be dropped first, or the slope is meaningless.
    """
    h = np.asarray(h, dtype=float)
    err = np.asarray(err, dtype=float)
    keep = (err > 1e-14) & np.isfinite(err)
    if keep.sum() < 2:
        raise ValueError(f"need >=2 usable error values, got {err!r} (all at round-off? use a coarser grid range)")
    p, _ = np.polyfit(np.log(h[keep]), np.log(err[keep]), 1)
    return float(p)


def pairwise_orders(h: Sequence[float], err: Sequence[float]) -> list[float]:
    """Order between each consecutive pair — print these; a drifting sequence means you are not in the asymptotic range."""
    h = np.asarray(h, float)
    err = np.asarray(err, float)
    return [float(np.log(err[i] / err[i + 1]) / np.log(h[i] / h[i + 1])) for i in range(len(h) - 1)]


def richardson(coarse: float, fine: float, r: float = 2.0, p: float = 2.0) -> float:
    """Richardson-extrapolated value from two grids refined by factor ``r`` for a scheme of order ``p``."""
    return float(fine + (fine - coarse) / (r ** p - 1.0))


@dataclass
class Study:
    h: list[float]
    err: list[float]
    order: float
    pairwise: list[float]

    def __str__(self) -> str:  # pragma: no cover - reporting only
        rows = "\n".join(f"    h={hh:<12.6g} err={ee:<12.6g}" for hh, ee in zip(self.h, self.err))
        return f"{rows}\n    observed order = {self.order:.3f}  (pairwise: {[round(p, 3) for p in self.pairwise]})"


def refinement_study(solve: Callable[[int], float], ns: Sequence[int], h_of_n: Callable[[int], float] | None = None) -> Study:
    """Run ``solve(n) -> error`` on each resolution and summarise.

    ``solve`` returns the error for resolution ``n`` (it owns the exact/manufactured solution). ``h_of_n`` defaults to
    ``1/n``. Assert on ``study.order`` in the test and print ``study`` so the report can quote the table.
    """
    h_of_n = h_of_n or (lambda n: 1.0 / n)
    hs = [float(h_of_n(n)) for n in ns]
    es = [float(solve(n)) for n in ns]
    return Study(hs, es, observed_order(hs, es), pairwise_orders(hs, es))


def mms(expr, symbols, operator):
    """Method of manufactured solutions: return (u_exact_fn, source_fn) as fast numpy callables.

    ``expr``    : a sympy expression for the manufactured field, e.g. ``sin(pi*x)*cos(pi*y)``
    ``symbols`` : the sympy symbols in evaluation order, e.g. ``(x, y)``
    ``operator``: a callable applying the *continuous* differential operator symbolically, e.g.
                  ``lambda u: sp.diff(u, x, 2) + sp.diff(u, y, 2)``.

    The source term is ``operator(expr)``; feed it to the solver, solve, and compare with ``u_exact``. Any consistent
    discretisation must then converge at its design order — this works even when no exact solution of the real problem
    is known.
    """
    import sympy as sp

    src = sp.simplify(operator(expr))
    return sp.lambdify(symbols, expr, "numpy"), sp.lambdify(symbols, src, "numpy")
