"""Linear (hydrodynamic) stability toolkit: Chebyshev collocation, constrained generalised eigenproblems, the normal-mode
vocabulary, the Orr–Sommerfeld / Rayleigh / Taylor–Goldstein eigen-solvers, parameter sweeps (maximum growth, neutral
curves, critical points), the classical theorems as checks (Howard's semicircle, inflection points, Squire's map) and the
disturbance-energy budget of a computed mode.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, Ch. 11 — §11.2 (11.1); §11.7 (11.61), (11.62), (11.71), (11.72);
§11.8 (11.77)–(11.80); §11.9 (11.81)–(11.84); §11.10 (11.88).  First used in Ch. 11; reused later for barotropic and
baroclinic instability (Ch. 13: Rayleigh–Kuo, Eady/Charney modes), transition and Reynolds stresses (Ch. 12).
Contract: ``analysis/ch11_design.md`` Part C.1 (names and signatures kept).

Conventions (every solver)
--------------------------
* Normal mode q(x, y, t) = q̂(y) exp{ik(x − ct)}, c = c_r + i c_i complex phase speed, growth rate k c_i (§11.7–§11.9);
  σ = −i|K|c converts to the e^{σt} convention of §11.4/§11.6 (:func:`sigma_from_c`).
* Chebyshev–Gauss–Lobatto nodes ξ_j = cos(jπ/N), j = 0 … N (**descending**: ``x[0]`` / ``grid.y[0]`` is the right/top end).
  Physical coordinate y = y(ξ) through a map ("linear"; "tan" for (−∞, ∞) truncated at ±y_max; "algebraic" for [0, y_max]);
  derivatives by the chain rule D_y = diag(1/y′) D_ξ, D_yy = diag(1/y′²) D_ξ² − diag(y″/y′³) D_ξ.
* Boundary conditions: the solvers **eliminate** boundary unknowns (u_e = −C_e⁻¹ C_k u_k, equations at eliminated nodes
  dropped) — equivalent to row replacement (:func:`apply_bc_rows`, kept for the teaching from-scratch version) but without
  infinite or spurious "boundary" eigenvalues (the row-replacement route gave spurious c ≈ 1 + 0.04i on the mapped Blasius
  grid).  DEVIATION from Part C's wording ("BC rows replaced in A and zeroed in B"): same equations, more robust algebra.
* ``filter=True`` = the **N-convergence filter**: keep only eigenvalues that have a partner within tol·max(1, |c|) in a second
  solve with N₂ = ⌈factor·N⌉ (:func:`converged_mask`).  Never a tight |c| cap alone — the TS mode is weak (c_i ~ 1e-3).
* Inputs are non-dimensional (lengths by the flow's L, velocities by U₀; Re = U₀L/ν); each wrapper in
  ``fluidpy.ch11_instability`` states its own scales.

Validation (planned; the verifier fills the final labels in ``tests/test_ch11.py``): V1 polynomial exactness and
−u″ = λu eigenvalues for :func:`cheb`; V5 Orszag (1971) c = 0.23752649 + 0.00373967i (plane Poiseuille, Re = 10⁴, k = 1);
V3 spectral convergence; V4 the energy budget closes for every computed OS mode.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np
import scipy.linalg as sla
from scipy.optimize import brentq, minimize_scalar

from ._util import as_scalar_if_0d

__all__ = [
    "SpectralGrid", "cheb", "cheb_matrices", "clenshaw_curtis_weights", "cheb_grid",
    "apply_bc_rows", "generalized_eigs", "constrained_eig", "converged_mask", "converged_eigs",
    "normal_mode", "sigma_from_c", "c_from_sigma", "stability_class", "stability_verdict", "marginal_type",
    "orr_sommerfeld_eigs", "os_mode", "rayleigh_eigs", "taylor_goldstein_eigs",
    "max_growth", "neutral_curve", "critical_point",
    "howard_semicircle", "in_howard_semicircle", "inflection_points", "squire_transform",
    "disturbance_energy_budget",
]


# ======================================================================================================================
# Chebyshev collocation
# ======================================================================================================================
def cheb(N: int, domain: Sequence[float] = (-1.0, 1.0)):
    """Chebyshev–Gauss–Lobatto differentiation matrix and nodes (Trefethen, *Spectral Methods in MATLAB*, 2000, `cheb.m`).

    Book: not in the book — the numerical method for the §11.4–§11.9 eigenproblems is ours (the book: "solutions can be
    obtained numerically", p. 515).  DEVIATION (method choice, not physics).
    Parameters
    ----------
    N : int ≥ 1 — polynomial degree (N + 1 nodes).
    domain : (a, b) — interval [a, b] (any units; the matrix carries 1/units).
    Returns
    -------
    D : (N+1, N+1) ndarray — d/dx on the nodes (exact for polynomials of degree ≤ N): off-diagonal c_i(−1)^{i+j}/(c_j(x_i − x_j)),
        c₀ = c_N = 2, else 1; diagonal by the negative-sum trick (rows annihilate constants); scaled by 2/(b − a).
    x : (N+1,) ndarray — nodes a + (b − a)(ξ_j + 1)/2, ξ_j = cos(jπ/N), **descending** (x[0] = b).
    Validation: V1 D x³ = 3x² to 1e-12 (N = 4: 4.4e-16); −u″ = λu, u(±1) = 0 → λ_n = (nπ/2)² (spectral, V3).  Label: analytic.
    """
    N = int(N)
    if N < 1:
        raise ValueError("cheb needs N >= 1")
    a, b = float(domain[0]), float(domain[1])
    if not b > a:
        raise ValueError("domain must be (a, b) with b > a")
    xi = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(xi, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))  # off-diagonal entries
    D = D - np.diag(D.sum(axis=1))  # negative-sum diagonal
    return D * (2.0 / (b - a)), a + (b - a) * (xi + 1.0) / 2.0


def cheb_matrices(N: int, orders: Sequence[int] = (1, 2, 4), domain: Sequence[float] = (-1.0, 1.0)) -> dict:
    """Powers of the Chebyshev differentiation matrix: dict(x, D1, D2, D4) (keys ``D{p}`` for each requested order).

    Book: method (ours).  Parameters: N; orders (positive ints); domain (a, b).  Label: analytic.
    """
    D, x = cheb(N, domain)
    out: dict = {"x": x}
    for p in orders:
        out[f"D{int(p)}"] = np.linalg.matrix_power(D, int(p))
    return out


def clenshaw_curtis_weights(N: int, domain: Sequence[float] = (-1.0, 1.0)) -> np.ndarray:
    """Clenshaw–Curtis quadrature weights on the nodes of :func:`cheb` (Trefethen `clencurt.m`): ∫f ≈ Σ w_j f(x_j).

    Exact for polynomials of degree ≤ N.  Parameters: N; domain (a, b).  Returns w [units of x] (Σw = b − a).
    Validation: V1 ∫₋₁¹ x² = 2/3 to 1e-15.  Label: analytic.
    """
    N = int(N)
    theta = np.pi * np.arange(N + 1) / N
    w = np.zeros(N + 1)
    ii = np.arange(1, N)
    v = np.ones(N - 1)
    if N % 2 == 0:
        w[0] = w[N] = 1.0 / (N ** 2 - 1)
        for k in range(1, N // 2):
            v -= 2.0 * np.cos(2 * k * theta[ii]) / (4 * k ** 2 - 1)
        v -= np.cos(N * theta[ii]) / (N ** 2 - 1)
    else:
        w[0] = w[N] = 1.0 / N ** 2
        for k in range(1, (N - 1) // 2 + 1):
            v -= 2.0 * np.cos(2 * k * theta[ii]) / (4 * k ** 2 - 1)
    w[ii] = 2.0 * v / N
    a, b = float(domain[0]), float(domain[1])
    return w * (b - a) / 2.0


@dataclass
class SpectralGrid:
    """A (possibly mapped) Chebyshev grid: nodes ``y`` (descending), d/dy ``D1``, d²/dy² ``D2``, quadrature weights ``w``.

    ``D3`` = D1·D2 and ``D4`` = D2·D2 on demand; ``map`` ∈ {"linear", "tan", "algebraic"}; ``domain`` the physical ends.
    """

    y: np.ndarray
    xi: np.ndarray
    D1: np.ndarray
    D2: np.ndarray
    w: np.ndarray
    map: str
    domain: tuple
    _cache: dict = field(default_factory=dict, repr=False)

    @property
    def N(self) -> int:
        return len(self.y) - 1

    @property
    def D3(self) -> np.ndarray:
        if "D3" not in self._cache:
            self._cache["D3"] = self.D1 @ self.D2
        return self._cache["D3"]

    @property
    def D4(self) -> np.ndarray:
        if "D4" not in self._cache:
            self._cache["D4"] = self.D2 @ self.D2
        return self._cache["D4"]

    def integrate(self, f):
        """∫ f dy by the grid's quadrature (Clenshaw–Curtis × y′(ξ))."""
        return as_scalar_if_0d(np.sum(self.w * np.asarray(f)))


def cheb_grid(N: int, domain: Sequence[float] = (-1.0, 1.0), map: str = "linear", y_max: float | None = None,
              s: float | None = None) -> SpectralGrid:
    """Build a :class:`SpectralGrid` on a bounded, semi-infinite or infinite interval.

    Book: method (ours) for the bounded (§11.7 lids z = 0, d; §11.8 walls y₁, y₂), semi-infinite (Blasius, §11.11) and
    unbounded (tanh, sech², §11.10) problems.  DEVIATION: ∞ truncated at ``y_max`` (ours); results depending on it are checked
    by changing it (each wrapper states its tested values).
    Parameters
    ----------
    N : degree.   domain : (a, b) — "linear": the interval; "tan": only the midpoint y₀ is used; "algebraic": a is the wall.
    map : "linear"; "tan" y = y₀ + s tan(θξ), θ = arctan(y_max/s) (nodes cluster near y₀; ends y₀ ± y_max);
          "algebraic" y = a + s(1 + ξ)/(β − ξ), β = 1 + 2s/y_max (half the nodes below a + s).
    y_max : truncation for "tan"/"algebraic".   s : map scale (default 1.0 "tan", 2.0 "algebraic"); smaller s clusters more
        nodes at the centre (critical layers of near-neutral modes: s = 0.5 converges tanh/Taylor–Goldstein growth rates to 1e-9
        at N = 120; s = 3 does not).
    Returns SpectralGrid.  Label: analytic.
    """
    D, xi = cheb(N)
    a, b = float(domain[0]), float(domain[1])
    if map == "linear":
        y = a + (b - a) * (xi + 1.0) / 2.0
        yp = np.full_like(xi, (b - a) / 2.0)
        ypp = np.zeros_like(xi)
        dom = (a, b)
    elif map == "tan":
        if y_max is None or y_max <= 0:
            raise ValueError("map='tan' needs y_max > 0 (truncation of the infinite interval)")
        s = 1.0 if s is None else float(s)
        th = math.atan(y_max / s)
        y0 = 0.5 * (a + b)
        y = y0 + s * np.tan(th * xi)
        sec2 = 1.0 / np.cos(th * xi) ** 2
        yp = s * th * sec2
        ypp = 2.0 * s * th ** 2 * sec2 * np.tan(th * xi)
        dom = (y0 - y_max, y0 + y_max)
    elif map == "algebraic":
        if y_max is None or y_max <= 0:
            raise ValueError("map='algebraic' needs y_max > 0 (truncation of the semi-infinite interval)")
        s = 2.0 if s is None else float(s)
        beta = 1.0 + 2.0 * s / y_max
        y = a + s * (1.0 + xi) / (beta - xi)
        yp = s * (beta + 1.0) / (beta - xi) ** 2
        ypp = 2.0 * s * (beta + 1.0) / (beta - xi) ** 3
        dom = (a, a + y_max)
    else:
        raise ValueError('map must be "linear", "tan" or "algebraic"')
    D1 = D / yp[:, None]
    D2 = (D @ D) / yp[:, None] ** 2 - (ypp / yp ** 3)[:, None] * D
    return SpectralGrid(y=y, xi=xi, D1=D1, D2=D2, w=clenshaw_curtis_weights(N) * yp, map=map, domain=dom)


# ======================================================================================================================
# Generalised eigenproblems with boundary conditions
# ======================================================================================================================
def apply_bc_rows(A: np.ndarray, B: np.ndarray, rows: Sequence[tuple[int, np.ndarray]]):
    """Row-replacement boundary conditions: A[r] = vec, B[r] = 0 for each (r, vec) — each becomes an infinite eigenvalue.

    Book: method (ours; the teaching version of C03).  Returns copies (A, B) (complex).  Label: analytic.
    """
    A = np.array(A, dtype=complex, copy=True)
    B = np.array(B, dtype=complex, copy=True)
    for r, vec in rows:
        A[r] = vec
        B[r] = 0.0
    return A, B


def _order(w: np.ndarray, sort: str | None) -> np.ndarray:
    if sort is None:
        return np.arange(len(w))
    if sort == "imag":
        return np.argsort(-w.imag, kind="stable")
    if sort == "real":
        return np.argsort(-w.real, kind="stable")
    if sort == "abs":
        return np.argsort(np.abs(w), kind="stable")
    raise ValueError('sort must be "imag", "real", "abs" or None')


def generalized_eigs(A: np.ndarray, B: np.ndarray, sort: str = "imag", c_max: float | None = None,
                     return_vectors: bool = False):
    """Solve A x = λ B x (``scipy.linalg.eig``), drop non-finite λ and (optionally) |λ| ≥ c_max, sort.

    Parameters: A, B; sort "imag" (descending Im λ — growth order for c), "real" (descending Re λ — order for σ), "abs"
    (ascending |λ|) or None; c_max; return_vectors.  Returns λ (and V, same column order).  Label: analytic.
    """
    if return_vectors:
        w, V = sla.eig(A, B)
    else:
        w, V = sla.eig(A, B, right=False), None
    ok = np.isfinite(w)
    if c_max is not None:
        ok &= np.abs(w) < c_max
    w = w[ok]
    order = _order(w, sort)
    w = w[order]
    if V is not None:
        return w, V[:, ok][:, order]
    return w


def constrained_eig(A: np.ndarray, B: np.ndarray, C: np.ndarray, eliminate: Sequence[int], return_vectors: bool = False,
                    sort: str = "imag"):
    """Solve A u = λ B u subject to C u = 0 by eliminating the unknowns ``eliminate`` (and dropping their equations).

    u_e = −C_e⁻¹ C_k u_k; the reduced (n − m)-square problem has no infinite eigenvalues from boundary rows.
    Parameters: A, B (n × n; rows = collocation equations, row index = unknown index); C (m × n); eliminate (m indices,
    C[:, eliminate] invertible); return_vectors (full vectors U = T u_k); sort.  Returns λ (finite, sorted) [, U].
    Book: method (ours).  Label: analytic.
    """
    A = np.asarray(A)
    n = A.shape[0]
    elim = np.asarray(eliminate, dtype=int)
    keep = np.setdiff1d(np.arange(n), elim)
    C = np.asarray(C)
    Ce = C[:, elim]
    if np.linalg.cond(Ce) > 1e14:
        raise np.linalg.LinAlgError("constrained_eig: the constraint block for the eliminated unknowns is singular")
    X = -np.linalg.solve(Ce, C[:, keep])
    T = np.zeros((n, len(keep)), dtype=complex)
    T[keep, np.arange(len(keep))] = 1.0
    T[elim] = X
    Ar = A[keep] @ T
    Br = np.asarray(B)[keep] @ T
    if return_vectors:
        w, V = sla.eig(Ar, Br)
    else:
        w, V = sla.eig(Ar, Br, right=False), None
    ok = np.isfinite(w)
    w = w[ok]
    order = _order(w, sort)
    w = w[order]
    if return_vectors:
        return w, (T @ V[:, ok])[:, order]
    return w


def converged_mask(w: np.ndarray, w_ref: np.ndarray, tol: float = 1e-6, ci_rtol: float = 0.0) -> np.ndarray:
    """True where w_i has a partner in w_ref within tol·max(1, |w_i|) — or, with ``ci_rtol`` > 0, within ci_rtol·|Im w_i|.

    The second test suits inviscid modes with a critical layer: a physical unstable mode converges slowly (|Δc| ~ 1e-5) but
    by much less than its own c_i, while spurious modes of the discretised continuous spectrum drift by ~30 % of their c_i
    between N and 1.5N.  Label: analytic.
    """
    w = np.atleast_1d(w)
    w_ref = np.atleast_1d(w_ref)
    if w_ref.size == 0:
        return np.zeros(w.shape, dtype=bool)
    d = np.min(np.abs(w[:, None] - w_ref[None, :]), axis=1)
    return (d < tol * np.maximum(1.0, np.abs(w))) | (d < ci_rtol * np.abs(w.imag))


def converged_eigs(eig_fn: Callable[[int], np.ndarray], N: int, factor: float = 1.5, tol: float = 1e-6,
                   sort: str = "imag") -> np.ndarray:
    """Eigenvalues of ``eig_fn(N)`` that reappear (within tol·max(1, |λ|)) in ``eig_fn(round(factor·N))``.

    Book: method (ours) — "refine and keep only the modes that converge" (math-to-python §2); the C03 figure of spurious modes.
    Parameters: eig_fn(N) → eigenvalue array; N; factor; tol; sort.  Returns λ sorted.  Label: converged.
    """
    w = np.atleast_1d(eig_fn(int(N)))
    w2 = np.atleast_1d(eig_fn(int(round(factor * N))))
    w = w[converged_mask(w, w2, tol)]
    return w[_order(w, sort)]


# ======================================================================================================================
# §11.2 normal-mode vocabulary
# ======================================================================================================================
def sigma_from_c(K, c):
    """σ = −i|K|c (growth rate σ_r = |K|c_i, frequency −σ_i = |K|c_r) from the complex phase speed.

    Book: §11.2, Eq. (11.1): exp{ikx + imy + σt} = exp{i|K|(e_K·x − ct)}.  Parameters: K [1/m or –] = |K| > 0; c [m/s or –]
    complex.  Returns σ complex.  Example: K = 2, c = 1 + 0.5i → σ = 1 − 2i.  Validation: V1 round trip.  Label: analytic.
    """
    return as_scalar_if_0d(-1j * np.asarray(K) * np.asarray(c, dtype=complex))  # Eq. (11.1): σ = −i|K|c


def c_from_sigma(K, sigma):
    """c = iσ/|K| (inverse of :func:`sigma_from_c`).  Book: §11.2, Eq. (11.1).  Label: analytic."""
    K = np.asarray(K, dtype=float)
    if np.any(K == 0):
        raise ValueError("c_from_sigma: |K| must be > 0 (a phase speed needs a wavenumber)")
    return as_scalar_if_0d(1j * np.asarray(sigma, dtype=complex) / K)


def normal_mode(x, y, t, uhat, k, m=0.0, sigma=None, c=None):
    """Real normal-mode field u = Re{û exp(ikx + imy + σt)} — Eq. (11.1), either form.

    Book: §11.2, Eq. (11.1): u = û(z) exp{ikx + imy + σt} = û(z) exp{i|K|(e_K·x − ct)}, K = (k, m, 0); real part taken (§7.7).
    Parameters: x, y [m or –] (broadcastable); t [s or –]; uhat complex amplitude at the level of interest; k, m real
    wavenumbers; sigma (first form) or c (second form) — exactly one.  Returns real ndarray (or float).
    Validation: V1 the two forms agree when σ = −i|K|c; parity with ``core.waves.real_field``.  Label: analytic.
    """
    if (sigma is None) == (c is None):
        raise ValueError("give exactly one of sigma (first form of (11.1)) or c (second form)")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    if sigma is not None:
        phase = 1j * k * x + 1j * m * y + sigma * t  # Eq. (11.1), first form
    else:
        K = math.hypot(k, m)
        if K == 0:
            raise ValueError("the travelling-wave form of (11.1) needs |K| > 0")
        phase = 1j * K * ((k * x + m * y) / K - c * t)  # Eq. (11.1), second form
    return as_scalar_if_0d(np.real(np.asarray(uhat) * np.exp(phase)))


def stability_class(sigma=None, c=None, tol: float = 1e-10) -> str:
    """"unstable" (σ_r or c_i > tol), "neutral" (|·| ≤ tol), "stable" (< −tol) for one mode (or the max over an array).

    Book: §11.2 — σ_r < 0 or c_i < 0 stable; = 0 neutrally stable; > 0 unstable.  Parameters: sigma or c (complex/array); tol.
    Returns str.  Label: analytic.
    """
    if (sigma is None) == (c is None):
        raise ValueError("give exactly one of sigma or c")
    vals = np.asarray(sigma if sigma is not None else c, dtype=complex).ravel()
    if vals.size == 0:
        raise ValueError("stability_class: no eigenvalues given")
    g = np.max(vals.real) if sigma is not None else np.max(vals.imag)
    if g > tol:
        return "unstable"
    if g >= -tol:
        return "neutral"
    return "stable"


def stability_verdict(growth, tol: float = 1e-10) -> dict:
    """"For every wavenumber" verdict from growth rates over a k grid: unstable if any > tol, neutral if the max is within
    tol, else stable.

    Book: §11.2 ("the system is unstable if σ_r or c_i are positive for *any* value of the wave number").
    Parameters: growth (array over k: σ_r or kc_i); tol.  Returns dict(verdict, k_index_max, max_growth).  Label: analytic.
    """
    g = np.asarray(growth, dtype=float).ravel()
    if g.size == 0:
        raise ValueError("stability_verdict: empty growth array")
    j = int(np.argmax(g))
    gm = float(g[j])
    verdict = "unstable" if gm > tol else ("neutral" if gm >= -tol else "stable")
    return dict(verdict=verdict, k_index_max=j, max_growth=gm)


def marginal_type(sigma, tol: float = 1e-8) -> str:
    """Marginal state (σ_r = 0): "stationary" if σ_i = 0 (cellular convection, secondary flow), else "oscillatory".

    Book: §11.2, p. 476 (overstability = oscillatory mode).  Parameters: sigma (complex, e^{σt} convention); tol.
    Raises ValueError if |σ_r| > tol·max(1, |σ|) (not marginal).  Label: analytic.
    """
    v = complex(sigma)
    if abs(v.real) > tol * max(1.0, abs(v)):
        raise ValueError(f"not a marginal state: σ_r = {v.real:.3e} ≠ 0")
    return "stationary" if abs(v.imag) <= tol * max(1.0, abs(v)) else "oscillatory"


# ======================================================================================================================
# Orr–Sommerfeld (11.79)–(11.80)
# ======================================================================================================================
def _grid_for(bc: str, domain, N: int, y_max, map, s) -> SpectralGrid:
    if bc == "wall":
        return cheb_grid(N, domain, map or "linear")
    if bc == "semi_infinite":
        ym = 20.0 if y_max is None else float(y_max)
        mp = map or "linear"
        a = float(domain[0])
        if mp == "linear":
            return cheb_grid(N, (a, a + ym), "linear")
        return cheb_grid(N, (a, a + ym), mp, y_max=ym, s=s)
    if bc in ("unbounded", "decay"):
        ym = 30.0 if y_max is None else float(y_max)
        mp = map or "tan"
        y0 = 0.5 * (float(domain[0]) + float(domain[1]))
        if mp == "linear":
            return cheb_grid(N, (y0 - ym, y0 + ym), "linear")
        return cheb_grid(N, (y0 - ym, y0 + ym), mp, y_max=ym, s=s)
    raise ValueError('bc must be "wall", "semi_infinite" or "decay" ("unbounded")')


def _parity_mask(V: np.ndarray, parity: str | None) -> np.ndarray:
    if parity is None:
        return np.ones(V.shape[1], dtype=bool)
    R = V[::-1]
    even = np.linalg.norm(V - R, axis=0)
    odd = np.linalg.norm(V + R, axis=0)
    if parity == "even":
        return even < odd
    if parity == "odd":
        return odd < even
    raise ValueError('parity must be None, "even" or "odd"')


def _os_solve(k, Re, U, Upp, g: SpectralGrid, vectors: bool):
    n = len(g.y)
    I = np.eye(n)
    L = g.D2 - k ** 2 * I
    Uy = np.asarray(U(g.y), dtype=float) * np.ones(n)
    Uppy = np.asarray(Upp(g.y), dtype=float) * np.ones(n)
    visc = (g.D4 - 2.0 * k ** 2 * g.D2 + k ** 4 * I) / (1j * k * Re)
    A = Uy[:, None] * L - np.diag(Uppy) - visc  # Eq. (11.79): (U − c)(φ″ − k²φ) − U″φ = (1/ikRe)(φ⁗ − 2k²φ″ + k⁴φ)
    B = L.astype(complex)
    C = np.vstack([I[0], I[n - 1], g.D1[0], g.D1[n - 1]])  # Eq. (11.80): φ = dφ/dy = 0 at both ends
    return constrained_eig(A, B, C, [0, n - 1, 1, n - 2], return_vectors=vectors)


def _finish(w, V, g, keep, return_vectors, key="phi"):
    w = w[keep]
    if return_vectors:
        return {"c": w, "y": g.y, key: V[:, keep], "grid": g}
    return w


def orr_sommerfeld_eigs(k: float, Re: float, U: Callable, Upp: Callable, domain: Sequence[float] = (-1.0, 1.0),
                        N: int = 100, bc: str = "wall", y_max: float | None = None, map_scale: float | None = None,
                        return_vectors: bool = False, filter: bool = True, parity: str | None = None,
                        map: str | None = None, factor: float = 1.25, tol: float = 1e-6, c_max: float | None = None):
    """Eigenvalues c of the Orr–Sommerfeld problem for a parallel flow U(y), sorted by descending c_i.

    Book: §11.8, Eq. (11.79)  (U − c)(φ″ − k²φ) − U″φ = (1/(ikRe))(φ⁗ − 2k²φ″ + k⁴φ), with no-slip Eq. (11.80)
    φ = φ′ = 0 at y₁, y₂ (û = φ′, v̂ = −ikφ, u = ∂ψ/∂y, v = −∂ψ/∂x — the §11.8 stream-function sign, slip S9).
    Non-dimensional: y by L, U by U₀, t by L/U₀, Re = U₀L/ν (each wrapper states L, U₀).
    Parameters
    ----------
    k > 0, Re > 0 [–].   U, Upp : callables y ↦ U, U″ (exact U″, e.g. −½ff″ for Blasius — never a numerical one).
    domain : walls (y₁, y₂) for "wall"; the wall y₁ = domain[0] for "semi_infinite"; the centre for "decay".
    N : degree (default 100; FAST 60; never > 150).
    bc : "wall"; "semi_infinite" (φ = φ′ = 0 also at y₁ + y_max — DEVIATION from Part C's algebraic map: a **linear** map on
         [0, y_max] by default, because the mapped grid produced spurious modes (c ≈ 1 + 0.04i) at N = 80–120); "decay" /
         "unbounded" (φ = φ′ = 0 at y₀ ± y_max, tan map).
    y_max : truncation [–] (defaults 20 semi-infinite, 30 decay).   map_scale : map scale s (see :func:`cheb_grid`).
    return_vectors : return dict(c, y, phi (columns), grid) instead of c.
    filter : N-convergence filter (second solve at ⌈factor·N⌉, keep c within tol·max(1, |c|)).
    parity : None, "even" or "odd" φ about the domain centre (symmetric problems).   map : override the map name.
    c_max : optional |c| cap after filtering.
    Returns c (complex ndarray) or dict.
    Assumptions: parallel base flow (exact for channels, an approximation for boundary layers, §11.11); 2-D disturbances.
    Validation: V5 plane Poiseuille Re = 10⁴, k = 1 → c = 0.23752649 + 0.00373967i (Orszag 1971; N = 60–120 agree to 1e-8);
    V3; V4 :func:`disturbance_energy_budget` closes for every mode; V7 Re → ∞ approaches :func:`rayleigh_eigs`.
    """
    if k <= 0 or Re <= 0:
        raise ValueError("orr_sommerfeld_eigs: k > 0 and Re > 0 are required")
    g = _grid_for(bc, domain, N, y_max, map, map_scale)
    need_vec = return_vectors or parity is not None
    out = _os_solve(k, Re, U, Upp, g, need_vec)
    w, V = out if need_vec else (out, None)
    keep = np.ones(len(w), dtype=bool)
    if parity is not None:
        keep &= _parity_mask(V, parity)
    if filter:
        g2 = _grid_for(bc, domain, int(math.ceil(factor * N)), y_max, map, map_scale)
        keep &= converged_mask(w, _os_solve(k, Re, U, Upp, g2, False), tol)
    if c_max is not None:
        keep &= np.abs(w) < c_max
    return _finish(w, V, g, keep, return_vectors)


def os_mode(k: float, Re: float, U: Callable, Upp: Callable, domain: Sequence[float] = (-1.0, 1.0), N: int = 100,
            **kw) -> dict:
    """The leading (largest c_i) Orr–Sommerfeld mode, normalised so max|û| = 1 with û real and positive at that point.

    Book: §11.8, Eqs. (11.79)–(11.80); û = dφ/dy, v̂ = −ikφ (§11.8, p. 510).
    Parameters: k, Re; U, Upp callables; domain; N; **kw passed to :func:`orr_sommerfeld_eigs` (bc, y_max, map_scale, parity,
    filter (default False here: the leading mode is checked by the energy budget instead), …).
    Returns dict(c, y, phi, u_hat (= φ′), v_hat (= −ikφ), dphi, grid, k, Re).  Label: converged.
    """
    kw.setdefault("filter", False)
    res = orr_sommerfeld_eigs(k, Re, U, Upp, domain=domain, N=N, return_vectors=True, **kw)
    if len(res["c"]) == 0:
        raise RuntimeError("os_mode: no eigenvalue survived the filters")
    g = res["grid"]
    phi = res["phi"][:, 0]
    dphi = g.D1 @ phi
    j = int(np.argmax(np.abs(dphi)))
    scale = dphi[j]
    phi = phi / scale
    dphi = dphi / scale
    return dict(c=complex(res["c"][0]), y=g.y, phi=phi, u_hat=dphi, v_hat=-1j * k * phi, dphi=dphi, grid=g, k=k, Re=Re)


# ======================================================================================================================
# Rayleigh (11.81)–(11.82) and Taylor–Goldstein (11.61)–(11.62)
# ======================================================================================================================
def _rayleigh_solve(k, U, Upp, g: SpectralGrid, vectors: bool):
    n = len(g.y)
    I = np.eye(n)
    L = g.D2 - k ** 2 * I
    Uy = np.asarray(U(g.y), dtype=float) * np.ones(n)
    Uppy = np.asarray(Upp(g.y), dtype=float) * np.ones(n)
    A = Uy[:, None] * L - np.diag(Uppy)  # Eq. (11.81): (U − c)(φ″ − k²φ) − U″φ = 0  ⇒  [U L − U″] φ = c L φ
    C = np.vstack([I[0], I[n - 1]])  # Eq. (11.82): φ = 0 at y₁, y₂ (or at ±y_max)
    return constrained_eig(A, L, C, [0, n - 1], return_vectors=vectors)


def rayleigh_eigs(k: float, U: Callable, Upp: Callable, domain: Sequence[float] = (-1.0, 1.0), N: int = 120,
                  bc: str = "wall", map_scale: float | None = None, parity: str | None = None,
                  return_vectors: bool = False, filter: bool = True, y_max: float | None = None, map: str | None = None,
                  factor: float = 1.5, tol: float = 1e-6, unstable_only: bool = False, ci_min: float = 1e-4,
                  ci_rtol: float = 0.01):
    """Eigenvalues c of Rayleigh's inviscid equation, sorted by descending c_i.

    Book: §11.9, Eq. (11.81) (U − c)(φ″ − k²φ) − U″φ = 0 with φ = 0 at the walls, Eq. (11.82) (bc="decay": φ = 0 at ±y_max,
    the unbounded tanh/sech² layers); the Re → ∞ limit of (11.79); also (11.95) for v̂ = −ikφ (Ex. 11.12).
    Parameters: k > 0; U, Upp callables; domain; N (default 120; FAST 60); bc ("wall", "decay", "semi_infinite");
    map_scale (tan-map scale; DEVIATION from Part C's y = aξ/√(1 − ξ²) map: a truncated tan map, default s = 1, y_max = 30 —
    same eigenvalues, tested); parity (None/"even"/"odd" φ: "even" = sinuous, "odd" = varicose, Ex. 11.12);
    return_vectors (dict(c, y, phi, grid)); filter (N-convergence, factor 1.5: keep c that moves < tol·max(1, |c|) or < ci_rtol·c_i
    — critical-layer modes converge slowly but by ≪ c_i, spurious ones drift by ~30 % of c_i); y_max; map; unstable_only
    (keep c_i > ci_min, default 1e-4: the continuous-spectrum end points c ≈ U_min, U_max carry c_i ~ 1e-6 junk).
    Returns c or dict.  Assumptions: inviscid disturbances on a (possibly viscous) base profile.  The continuous spectrum
    [U_min, U_max] appears as real c that the filter mostly removes; neutral modes (c_i → 0) converge slowly (critical layer).
    Validation: V1 tanh neutral k = 1, c = 0 (φ = sech y); Bickley k = 2 (sinuous), k = 1 (varicose), c = 2/3;
    V5 tanh most-amplified k = 0.4449, kc_i = 0.18970 (Michalke 1964: 0.4446, 0.1897); V4 conjugate pairs; semicircle (11.72).
    """
    if k <= 0:
        raise ValueError("rayleigh_eigs: k > 0 is required")
    g = _grid_for(bc, domain, N, y_max, map, map_scale)
    need_vec = return_vectors or parity is not None
    out = _rayleigh_solve(k, U, Upp, g, need_vec)
    w, V = out if need_vec else (out, None)
    keep = np.ones(len(w), dtype=bool)
    if parity is not None:
        keep &= _parity_mask(V, parity)
    if unstable_only:
        keep &= w.imag > ci_min
    if filter and keep.any():
        g2 = _grid_for(bc, domain, int(math.ceil(factor * N)), y_max, map, map_scale)
        keep &= converged_mask(w, _rayleigh_solve(k, U, Upp, g2, False), tol, ci_rtol)
    return _finish(w, V, g, keep, return_vectors)


def _tg_solve(k, U, Upp, N2, g: SpectralGrid, vectors: bool):
    n = len(g.y)
    inner = slice(1, n - 1)  # Eq. (11.62): ψ̂ = 0 at both ends → interior unknowns only
    I = np.eye(n)
    L = (g.D2 - k ** 2 * I)[inner, inner]
    Uy = (np.asarray(U(g.y), dtype=float) * np.ones(n))[inner]
    Uppy = (np.asarray(Upp(g.y), dtype=float) * np.ones(n))[inner]
    N2y = (np.asarray(N2(g.y), dtype=float) * np.ones(n))[inner]
    # Eq. (11.61) × (U − c):  (U − c)²(ψ̂″ − k²ψ̂) − U″(U − c)ψ̂ + N²ψ̂ = 0  →  c²M₂ + cM₁ + M₀ = 0
    M2 = L
    M1 = -2.0 * Uy[:, None] * L + np.diag(Uppy)
    M0 = (Uy ** 2)[:, None] * L - np.diag(Uy * Uppy) + np.diag(N2y)
    m = L.shape[0]
    Z = np.zeros((m, m))
    Id = np.eye(m)
    A = np.block([[Z, Id], [-M0, -M1]])  # companion linearisation in v = [ψ̂, cψ̂]
    Bm = np.block([[Id, Z], [Z, M2]])
    if vectors:
        w, V = sla.eig(A, Bm)
    else:
        w, V = sla.eig(A, Bm, right=False), None
    ok = np.isfinite(w)
    w = w[ok]
    order = np.argsort(-w.imag, kind="stable")
    w = w[order]
    if not vectors:
        return w
    Vf = np.zeros((n, int(ok.sum())), dtype=complex)
    Vf[inner] = V[:m][:, ok][:, order]
    return w, Vf


def taylor_goldstein_eigs(k: float, U: Callable, Upp: Callable, N2: Callable, domain: Sequence[float] = (0.0, 1.0),
                          N: int = 80, bc: str = "wall", map_scale: float | None = None, return_vectors: bool = False,
                          filter: bool = True, y_max: float | None = None, map: str | None = None, factor: float = 1.5,
                          tol: float = 1e-5, unstable_only: bool = True, ci_min: float = 1e-4, ci_rtol: float = 0.01):
    """Unstable (by default) eigenvalues c of the Taylor–Goldstein equation, sorted by descending c_i.

    Book: §11.7, Eq. (11.61)  (U − c)(ψ̂″ − k²ψ̂) − U″ψ̂ + N²ψ̂/(U − c) = 0 with rigid lids Eq. (11.62) ψ̂(0) = ψ̂(d) = 0
    (bc="decay": ψ̂ = 0 at ±y_max for unbounded layers); u = ∂ψ/∂z, w = −∂ψ/∂x (the §11.7 sign, slip S9); N² = −(g/ρ₀)dρ̄/dz
    (7.128).  Non-dimensional (z by the layer scale, U by its velocity scale, N² by (velocity/length)²).
    Method: × (U − c) → c²M₂ + cM₁ + M₀ = 0 with M₂ = D² − k², M₁ = −2U(D² − k²) + U″, M₀ = U²(D² − k²) − UU″ + N²;
    companion linearisation (P268).  The linearisation doubles the spurious modes (discretised continuous spectrum, small c_i
    that drift with N): the N-convergence filter is ON by default (factor 1.5, tol 1e-5).
    Parameters: k > 0; U, Upp, N2 callables; domain (default (0, 1)); N (default 80; FAST 60); bc ("wall", "decay");
    map_scale (tan-map scale; default None → 0.5 for "decay" — DEVIATION from Part C's 3.0, which converges poorly near the
    neutral curve: c_i at k = 0.9, J = 0 moves 0.043 → 0.065 from N = 80 to 180 with s = 3, but agrees to 1e-9 with s = 0.5);
    return_vectors (dict(c, y, psi, grid)); filter; y_max (default 30 for "decay"); map; factor; tol, ci_rtol (as
    :func:`rayleigh_eigs`); unstable_only
    (keep c_i > ci_min, default 1e-4 — the continuous-spectrum end points c ≈ U_min, U_max carry c_i ~ 1e-6 junk).
    Returns c or dict.  Assumptions: inviscid, Boussinesq, 2-D (Squire assumed in §11.7).
    Near-neutral modes (critical layer, Frobenius exponents ½ ± √(¼ − Ri)) converge slowly and are dropped by the filter within
    ≈ 0.02 of the neutral curve — growth there is labelled approximate.
    Validation: V1 N² = 0 equals :func:`rayleigh_eigs`; tanh/J sech²z exact neutral curve J = k(1 − k) (neutral mode
    ψ̂ = |tanh z|^{1−k} sech^k z, c = 0); V4 unstable c inside Howard's semicircle (11.72) and closed under conjugation; no
    unstable mode when Ri_min > ¼ (11.67); (11.65), (11.69)–(11.70) residuals ~1e-11 for computed modes.
    """
    if k <= 0:
        raise ValueError("taylor_goldstein_eigs: k > 0 is required")
    if map_scale is None and bc in ("decay", "unbounded"):
        map_scale = 0.5
    g = _grid_for(bc, domain, N, y_max, map, map_scale)
    out = _tg_solve(k, U, Upp, N2, g, return_vectors)
    w, V = out if return_vectors else (out, None)
    keep = np.ones(len(w), dtype=bool)
    if unstable_only:
        keep &= w.imag > ci_min
    if filter and keep.any():
        g2 = _grid_for(bc, domain, int(math.ceil(factor * N)), y_max, map, map_scale)
        keep &= converged_mask(w, _tg_solve(k, U, Upp, N2, g2, False), tol, ci_rtol)
    return _finish(w, V, g, keep, return_vectors, key="psi")


# ======================================================================================================================
# Parameter sweeps
# ======================================================================================================================
def _lead(c) -> complex:
    c = np.atleast_1d(c)
    return complex(c[np.argmax(c.imag)]) if c.size else complex(np.nan, -np.inf)


def max_growth(eig_fn: Callable[[float], complex], k_bounds: Sequence[float], n_scan: int = 24, xatol: float = 1e-7,
               measure: str = "kci") -> dict:
    """Most-amplified wavenumber: maximise k·c_i(k) (temporal growth rate; ``measure="ci"`` maximises c_i).

    Book: §11.7 (kc_i is the growth rate, p. 507), §11.10.  Method (ours): scan + bounded Brent around the best.
    Parameters: eig_fn(k) → leading c (or an array: its largest c_i is used); k_bounds; n_scan; xatol; measure.
    Returns dict(k, c, growth).  Example: tanh Rayleigh → k = 0.4449, kc_i = 0.1897.  Label: converged.
    """
    k_lo, k_hi = float(k_bounds[0]), float(k_bounds[1])

    def obj(k):
        c = _lead(eig_fn(k))
        g = c.imag * (k if measure == "kci" else 1.0)
        return -g if np.isfinite(g) else np.inf

    ks = np.linspace(k_lo, k_hi, int(n_scan))
    vals = np.array([obj(k) for k in ks])
    i = int(np.argmin(vals))
    lo, hi = ks[max(i - 1, 0)], ks[min(i + 1, len(ks) - 1)]
    r = minimize_scalar(obj, bounds=(lo, hi), method="bounded", options=dict(xatol=xatol))
    k = float(r.x) if r.fun <= vals[i] else float(ks[i])
    return dict(k=k, c=_lead(eig_fn(k)), growth=float(-obj(k)))


def neutral_curve(eig_fn2: Callable[[float, float], complex], Re_values: Sequence[float], k_bounds: Sequence[float],
                  n_k: int = 40, xtol: float = 1e-8, k_grid: Sequence[float] | None = None) -> dict:
    """Neutral curve c_i(k, Re) = 0: for each Re the k where c_i changes sign (scan + Brent); lower and upper branch.

    Book: §11.10 (Figs. 11.23, 11.24, 11.26 — marginal curves); §11.2 marginal state.
    Parameters: eig_fn2(k, Re) → leading c (or array); Re_values; k_bounds; n_k (or explicit k_grid); xtol.
    Roots are classified by direction: c_i rising through 0 with k = lower branch, falling = upper branch (NaN where the
    branch lies outside k_bounds — e.g. an adverse-gradient layer unstable down to the smallest scanned k).
    Returns dict(Re, k_lower, k_upper (NaN where none), k_roots (list), c_lower, c_upper).  Label: converged.
    """
    ks = np.asarray(k_grid, dtype=float) if k_grid is not None else np.linspace(k_bounds[0], k_bounds[1], int(n_k))
    ci = lambda k, Re: _lead(eig_fn2(k, Re)).imag  # noqa: E731
    roots, lo, hi, clo, chi = [], [], [], [], []
    for Re in Re_values:
        v = np.array([ci(k, Re) for k in ks])
        idx = [j for j in range(len(ks) - 1) if np.isfinite(v[j]) and np.isfinite(v[j + 1]) and v[j] * v[j + 1] < 0]
        rr = np.array([brentq(lambda k: ci(k, Re), ks[j], ks[j + 1], xtol=xtol) for j in idx])  # noqa: B023
        rising = [r for r, j in zip(rr, idx) if v[j] < 0]  # c_i turns positive: lower branch
        falling = [r for r, j in zip(rr, idx) if v[j] > 0]  # c_i turns negative: upper branch
        roots.append(rr)
        kl = rising[0] if rising else np.nan  # NaN: the lower branch lies below k_bounds (unstable at the scan's start)
        ku = falling[-1] if falling else np.nan
        lo.append(kl)
        hi.append(ku)
        clo.append(_lead(eig_fn2(kl, Re)) if np.isfinite(kl) else complex(np.nan))
        chi.append(_lead(eig_fn2(ku, Re)) if np.isfinite(ku) else complex(np.nan))
    return dict(Re=np.asarray(Re_values, dtype=float), k_lower=np.array(lo), k_upper=np.array(hi), k_roots=roots,
                c_lower=np.array(clo), c_upper=np.array(chi))


def critical_point(eig_fn2: Callable[[float, float], complex], k_bounds: Sequence[float], Re_bounds: Sequence[float],
                   k_guess: float | None = None, n_k: int = 9, xatol: float = 1e-7, rtol_Re: float = 1e-11) -> dict:
    """Critical point: min over k of the neutral Reynolds number Re_n(k) (c_i(k, Re_n) = 0).

    Book: §11.10 (Re_cr of Table 11.1), §11.11.  Method (ours, P263): Re_n(k) by Brent in Re on Re_bounds (needs c_i < 0 at the
    lower and > 0 at the upper bound), bounded Brent in k around the best of a scan (or around ``k_guess``).
    Returns dict(Re_c, k_c, c_r, c (complex)).  Raises RuntimeError with a clear message if no sign change is found.
    Example: plane Poiseuille → 5772.22, 1.02056, 0.26400 (Orszag 1971).  Label: converged.
    """
    Re_lo, Re_hi = float(Re_bounds[0]), float(Re_bounds[1])
    ci = lambda k, Re: _lead(eig_fn2(k, Re)).imag  # noqa: E731

    def Re_n(k):
        a, b = ci(k, Re_lo), ci(k, Re_hi)
        if not (a < 0 < b):
            return np.inf
        return brentq(lambda R: ci(k, R), Re_lo, Re_hi, rtol=rtol_Re, xtol=1e-12 * Re_hi)

    if k_guess is not None:
        span = 0.05 * (float(k_bounds[1]) - float(k_bounds[0]))
        lo, hi = max(k_guess - span, k_bounds[0]), min(k_guess + span, k_bounds[1])
        best_k, best = k_guess, Re_n(k_guess)
    else:
        ks = np.linspace(float(k_bounds[0]), float(k_bounds[1]), int(n_k))
        vals = np.array([Re_n(k) for k in ks])
        if not np.isfinite(vals).any():
            raise RuntimeError("critical_point: c_i does not change sign inside Re_bounds for any k in k_bounds — widen "
                               "the bounds or check the profile/spectrum filter")
        i = int(np.argmin(vals))
        lo, hi = ks[max(i - 1, 0)], ks[min(i + 1, len(ks) - 1)]
        best_k, best = float(ks[i]), float(vals[i])
    r = minimize_scalar(Re_n, bounds=(lo, hi), method="bounded", options=dict(xatol=xatol))
    k_c, Re_c = (float(r.x), float(r.fun)) if r.fun <= best else (float(best_k), float(best))
    if not np.isfinite(Re_c):
        raise RuntimeError("critical_point: no neutral point near k_guess")
    c = _lead(eig_fn2(k_c, Re_c))
    return dict(Re_c=Re_c, k_c=k_c, c_r=c.real, c=c)


# ======================================================================================================================
# Theorems as checks
# ======================================================================================================================
def howard_semicircle(Umin: float, Umax: float, n: int = 200):
    """Boundary of Howard's semicircle (upper half c-plane, diameter [U_min, U_max]).

    Book: §11.7, p. 507 (unnumbered, after (11.72)): [c_r − ½(U_max + U_min)]² + c_i² ≤ [½(U_max − U_min)]²; growth bound
    kc_i < (k/2)(U_max − U_min); Fig. 11.20.  Holds with and without stratification **assuming N² ≥ 0**.
    Returns (c_r, c_i) arrays.  Label: analytic.
    """
    if not Umax > Umin:
        raise ValueError("howard_semicircle: Umax > Umin required")
    th = np.linspace(0.0, np.pi, int(n))
    m, R = 0.5 * (Umax + Umin), 0.5 * (Umax - Umin)
    return m + R * np.cos(th), R * np.sin(th)


def in_howard_semicircle(c, Umin: float, Umax: float, tol: float = 1e-9):
    """True where c lies inside/on Howard's semicircle (p. 507) with c_r in the range of U, Eq. (11.71).

    Book: §11.7, Eqs. (11.71)–(11.72).  Parameters: c (complex, array ok); Umin, Umax; tol (relative to the radius).
    Returns bool (array).  Meaningful for unstable modes (c_i > 0).  Label: analytic.
    """
    c = np.asarray(c, dtype=complex)
    m, R = 0.5 * (Umax + Umin), 0.5 * (Umax - Umin)
    inside = (c.real - m) ** 2 + c.imag ** 2 <= R ** 2 * (1.0 + tol) + tol  # semicircle (completed square of (11.72))
    in_range = (c.real >= Umin - tol * R) & (c.real <= Umax + tol * R)  # Eq. (11.71)
    return as_scalar_if_0d(inside & in_range)


def _second_derivative(U: Callable, y, h: float = 1e-4):
    return (U(y + h) - 2.0 * U(y) + U(y - h)) / h ** 2


def inflection_points(y, U=None, Upp=None, tol: float = 1e-12) -> np.ndarray:
    """Points where U″ changes sign inside the interval (end points excluded).

    Book: §11.9, Eqs. (11.83)–(11.84): an inviscid instability needs U″ to change sign in y₁ < y < y₂ (Rayleigh).
    Parameters: y (array, the bracketing grid); U and/or Upp (callables → Brent; arrays → linear interpolation; only U callable →
    central-difference U″, h = 1e-4); tol (|U″| ≤ tol·max|U″| counts as zero).
    Returns sorted ndarray y_I (empty if none).  Validation: V1 sin y → 0; tanh → 0; Poiseuille → none.  Label: analytic.
    """
    y = np.sort(np.asarray(y, dtype=float))
    if Upp is None and U is None:
        raise ValueError("give U or Upp")
    f = None
    if Upp is not None and callable(Upp):
        f = Upp
    elif Upp is None and callable(U):
        f = lambda yy: _second_derivative(U, yy)  # noqa: E731
    if f is not None:
        vals = np.asarray(f(y), dtype=float) * np.ones(len(y))
    elif Upp is not None:
        vals = np.asarray(Upp, dtype=float)
    else:
        Ua = np.asarray(U, dtype=float)
        vals = np.gradient(np.gradient(Ua, y, edge_order=2), y, edge_order=2)
    scale = max(np.max(np.abs(vals)), 1e-300)
    s = np.where(np.abs(vals) <= tol * scale, 0.0, np.sign(vals))
    roots = []
    for j in range(len(y) - 1):
        if s[j] != 0 and s[j + 1] != 0 and s[j] * s[j + 1] < 0:
            roots.append(brentq(f, y[j], y[j + 1], xtol=1e-13) if f is not None else
                         y[j] - vals[j] * (y[j + 1] - y[j]) / (vals[j + 1] - vals[j]))
        elif s[j + 1] == 0 and 0 < j + 1 < len(y) - 1:  # exact zero on a node: count it if the sign changes across it
            nz = np.nonzero(s[j + 2:])[0]
            if s[j] != 0 and nz.size and s[j] * s[j + 2 + nz[0]] < 0:
                roots.append(y[j + 1])
    return np.array(sorted(set(np.round(roots, 12))))


def squire_transform(k: float, m: float, Re: float, uhat=None, what=None, phat=None) -> dict:
    """Squire's transformation of a 3-D normal mode to the equivalent 2-D problem.

    Book: §11.8, Eq. (11.78): k̄ = √(k² + m²), c̄ = c, k̄ū = kû + mŵ, v̄ = v̂, p̄/k̄ = p̂/k, k̄Re̅ = kRe ⇒ Re̅ = kRe/k̄ ≤ Re,
    2-D growth k̄c_i ≥ kc_i.  Parameters: k > 0, m ≥ 0; Re; optional uhat, what, phat (arrays).
    Returns dict(kbar, Rebar, ubar, pbar, angle_deg (= atan(m/k) in degrees), growth_factor (k̄/k)).
    Example: m/k = tan 30° → Rebar = 0.8660 Re.  Validation: V1 round trip; V7 m = 0; numerical with ``os_3d_eigs``.
    Label: analytic.
    """
    if k <= 0 or m < 0:
        raise ValueError("squire_transform: k > 0, m >= 0")
    kb = math.hypot(k, m)  # Eq. (11.78)
    out = dict(kbar=kb, Rebar=k * Re / kb, ubar=None, pbar=None, angle_deg=math.degrees(math.atan2(m, k)),
               growth_factor=kb / k)  # k̄Re̅ = kRe
    if uhat is not None and what is not None:
        out["ubar"] = (k * np.asarray(uhat) + m * np.asarray(what)) / kb
    if phat is not None:
        out["pbar"] = kb * np.asarray(phat) / k
    return out


# ======================================================================================================================
# Disturbance energy (11.88)
# ======================================================================================================================
def _grid_matching(y, grid: SpectralGrid | None) -> SpectralGrid:
    y = np.asarray(y, dtype=float)
    if grid is not None:
        if len(grid.y) != len(y) or not np.allclose(grid.y, y, rtol=0, atol=1e-10 * max(1.0, np.max(np.abs(y)))):
            raise ValueError("grid does not match y")
        return grid
    g = cheb_grid(len(y) - 1, (float(np.min(y)), float(np.max(y))), "linear")
    if not np.allclose(g.y, y, atol=1e-10 * max(1.0, np.max(np.abs(y)))):
        raise ValueError("y is not a linear Chebyshev–Lobatto grid; pass the solver's grid (return_vectors=True → 'grid')")
    return g


def disturbance_energy_budget(k: float, c: complex, phi, y, Up, Re: float, weights=None, grid: SpectralGrid | None = None
                              ) -> dict:
    """Kinetic-energy budget of a 2-D normal mode, averaged over one wavelength: dE/dt = production − dissipation.

    Book: §11.10, Eq. (11.88)  d/dt∫½u_i²dV = −∫u_iu_j ∂U_i/∂x_j dV − Λ, Λ = ν∫(∂u_i/∂x_j)²dV, and its 2-D form (p. 520)
    d/dt∫½(u² + v²)dV = −∫uv ∂U/∂y dV − Λ (control volume Fig. 11.25: walls or u_i → 0, integer wavelengths).
    With u = Re{ûe^{ik(x−ct)}}, û = φ′, v̂ = −ikφ (§11.8), ⟨ab⟩ = ½Re(âb̂*)e^{2kc_i t}; at t = 0:
    E = ¼∫(|û|² + |v̂|²)dy, P = −½∫Re(ûv̂*)U′dy, Λ = (1/(2Re))∫(k²|û|² + |û′|² + k²|v̂|² + |v̂′|²)dy, dE/dt = 2kc_iE.
    Parameters: k, c; phi (complex, on the nodes y); y; Up (callable or array: U′); Re (np.inf → inviscid); weights (quadrature
    weights; default Clenshaw–Curtis of the grid); grid (the solver's SpectralGrid — needed for mapped/semi-infinite grids;
    default: a linear Chebyshev grid matching y).
    Returns dict(E, dEdt, production, dissipation, residual (dEdt − P + Λ, ~1e-12 for an OS mode), ratio (P/Λ), uv (⟨uv⟩
    profile), production_density, dissipation_density, phase_uv [rad], y).  Example: Poiseuille Re = 10⁴, k = 1: P/Λ = 1.616.
    Validation: V4 residual ≈ 0 for every OS mode; production > dissipation iff c_i > 0.  Label: conserved.
    """
    g = _grid_matching(y, grid)
    w = g.w if weights is None else np.asarray(weights, dtype=float)
    phi = np.asarray(phi, dtype=complex)
    dphi = g.D1 @ phi
    d2phi = g.D2 @ phi
    uh, vh = dphi, -1j * k * phi
    duh, dvh = d2phi, -1j * k * dphi
    Upy = np.asarray(Up(g.y) if callable(Up) else Up, dtype=float) * np.ones(len(g.y))
    uv = 0.5 * np.real(uh * np.conj(vh))  # ⟨uv⟩ = ½Re(ûv̂*)
    E = 0.25 * np.sum(w * (np.abs(uh) ** 2 + np.abs(vh) ** 2))
    prod_density = -uv * Upy  # −⟨uv⟩ dU/dy
    P = np.sum(w * prod_density)  # Eq. (11.88), 2-D production
    nu_eff = 0.0 if not np.isfinite(Re) else 1.0 / Re
    diss_density = 0.5 * nu_eff * (np.abs(duh) ** 2 + k ** 2 * np.abs(uh) ** 2 + np.abs(dvh) ** 2 + k ** 2 * np.abs(vh) ** 2)
    Dd = np.sum(w * diss_density)  # Λ
    dEdt = 2.0 * k * complex(c).imag * E
    return dict(E=float(E), dEdt=float(dEdt), production=float(P), dissipation=float(Dd), residual=float(dEdt - (P - Dd)),
                ratio=float(P / Dd) if Dd != 0 else float("inf"), uv=uv, production_density=prod_density,
                dissipation_density=diss_density, phase_uv=np.angle(uh) - np.angle(vh), y=g.y)
