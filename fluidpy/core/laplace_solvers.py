"""Finite-difference Laplace and Poisson solvers on masked grids (§6.7): the five-point stencil (6.70)–(6.72), single
relaxation sweeps (Jacobi, Gauss–Seidel with the latest values, SOR) for step-by-step explainers, and a solver with a
residual history (Jacobi / Gauss–Seidel / SOR / sparse direct).

Book: Kundu, Cohen & Dowling 5e, §6.7, Eqs. (6.70)–(6.73), Example 6.2 (pages p253–p257). The book specifies the
average rule (6.72) and the Gauss–Seidel sweep "always using the latest available value at each point"; Jacobi, SOR,
the sparse direct solve and the residual-based stopping rule are our additions (labelled). Reused by Ch. 10 (pressure
Poisson, ω–ψ), Ch. 13 (PV inversion on a bounded domain), Ch. 8 (duct flow).

Layout (project rule, ``core.grids``): arrays are ``psi[j, i]`` = ψ(x_i, y_j) (y along axis 0, x along axis 1). A
boolean ``mask`` marks the unknown (interior) nodes; every other node holds a fixed (Dirichlet) value in the array.
Sweep order: ``order="x_outer"`` (alias ``"book"``) visits i (x) in the outer loop and j (y) inside — the book's
FORTRAN (DO I … DO J …); ``"y_outer"`` (alias ``"lex"``, the design's "i fastest, j slowest") is row-major. On
Example 6.2's grid both orders give the same iterates to round-off (checked: 50 sweeps differ by 9e-16). Residual: the defect of the average rule
r = (weighted neighbour average + source term) − ψ at each unknown node (= Δ²∇²_h ψ/4 for Δx = Δy), max-norm.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve, spsolve_triangular

from ._util import as_scalar_if_0d

__all__ = ["laplacian_5pt", "node_update", "jacobi_sweep", "gauss_seidel_sweep", "sor_sweep", "residual_field",
           "residual_norm", "sweep_order", "assemble_laplace", "solve_laplace", "solve_poisson", "optimal_sor_omega",
           "jacobi_spectral_radius"]

_S = as_scalar_if_0d


def _F(a):
    return np.asarray(a, dtype=float)


def laplacian_5pt(psi, dx: float = 1.0, dy: float | None = None, mask=None):
    """Five-point Laplacian (6.70)–(6.71): [ψ_{i+1,j} − 2ψ_{i,j} + ψ_{i−1,j}]/Δx² + [ψ_{i,j+1} − 2ψ_{i,j} + ψ_{i,j−1}]/Δy²
    at interior nodes (second-order accurate — the book calls the half-point differences "first-order", they are
    second-order accurate). Array ``psi[j, i]``; NaN on the outer frame and outside ``mask`` (if given).
    Returns an array like psi [unit of ψ / m²]. Book: §6.7, Eqs. (6.70)–(6.71).
    Validation (planned): V3 observed order 2 on sin(πx) sinh(πy); exact (zero) on xy and x² − y². Label: converged.
    """
    P = _F(psi)
    dy = dx if dy is None else float(dy)
    out = np.full(P.shape, np.nan)
    out[1:-1, 1:-1] = ((P[1:-1, 2:] - 2 * P[1:-1, 1:-1] + P[1:-1, :-2]) / dx ** 2  # Eq. (6.70)
                       + (P[2:, 1:-1] - 2 * P[1:-1, 1:-1] + P[:-2, 1:-1]) / dy ** 2)  # Eq. (6.71)
    if mask is not None:
        out = np.where(np.asarray(mask, bool), out, np.nan)
    return out


def _weights(dx, dy):
    dy = dx if dy is None else float(dy)
    cx, cy = 1.0 / dx ** 2, 1.0 / dy ** 2
    return cx, cy, 2.0 * (cx + cy)


def _with_bc(psi, mask, bc):
    P = _F(psi).copy()
    if bc is not None:
        B = _F(bc)
        fix = ~np.asarray(mask, bool) & np.isfinite(B)
        P[fix] = B[fix]
    return P


def node_update(psi, i: int, j: int, bc=None, dx: float = 1.0, dy: float | None = None, f=None) -> float:
    """The average rule (6.72) at node (i, j): ψ_{i,j} = ¼[ψ_{i−1,j} + ψ_{i+1,j} + ψ_{i,j−1} + ψ_{i,j+1}] for Δx = Δy
    (general Δx, Δy and a Poisson source f: [(ψ_E + ψ_W)/Δx² + (ψ_N + ψ_S)/Δy² − f]/(2/Δx² + 2/Δy²)).
    ``psi[j, i]``; ``bc`` (optional) supplies neighbour values where ``psi`` is NaN. Returns the new value (the array is
    not modified). Book: §6.7, Eq. (6.72). Label: analytic."""
    P = _F(psi)
    if bc is not None:
        P = np.where(np.isfinite(P), P, _F(bc))
    cx, cy, cc = _weights(dx, dy)
    src = 0.0 if f is None else float(_F(f)[j, i])
    return float((cx * (P[j, i - 1] + P[j, i + 1]) + cy * (P[j - 1, i] + P[j + 1, i]) - src) / cc)  # Eq. (6.72)


_ORDER_ALIAS = {"book": "x_outer", "lex": "y_outer", "x_outer": "x_outer", "y_outer": "y_outer"}


def sweep_order(mask, order: str = "lex") -> list:
    """List of (j, i) unknown nodes in sweep order: "x_outer"/"book" = for i: for j (the book's FORTRAN), "y_outer"/"lex"
    = for j: for i (i fastest). Label: analytic."""
    if order not in _ORDER_ALIAS:
        raise ValueError('order must be "lex", "book", "x_outer" or "y_outer"')
    order = _ORDER_ALIAS[order]
    M = np.asarray(mask, bool)
    J, I = np.nonzero(M)
    k = np.lexsort((J, I)) if order == "x_outer" else np.lexsort((I, J))
    return list(zip(J[k].tolist(), I[k].tolist()))


def jacobi_sweep(psi, mask, bc=None, dx: float = 1.0, dy: float | None = None, f=None):
    """One Jacobi sweep: every unknown node replaced by the average rule of the *old* values (our addition, the
    baseline the book's method improves on). ``bc`` (optional) is written into the non-mask nodes first.
    Returns (psi_new, max_change). Book: §6.7 (6.72). Label: analytic."""
    P = _with_bc(psi, mask, bc)
    M = np.asarray(mask, bool)
    cx, cy, cc = _weights(dx, dy)
    src = 0.0 if f is None else _F(f)
    new = P.copy()
    avg = np.zeros_like(P)
    avg[1:-1, 1:-1] = cx * (P[1:-1, :-2] + P[1:-1, 2:]) + cy * (P[:-2, 1:-1] + P[2:, 1:-1])
    new[M] = ((avg - src) / cc)[M]
    return new, float(np.max(np.abs(new[M] - P[M]))) if M.any() else 0.0


def gauss_seidel_sweep(psi, mask, bc=None, order: str = "lex", dx: float = 1.0, dy: float | None = None, f=None):
    """One Gauss–Seidel sweep: the average rule (6.72) node by node, "always using the latest available value at each
    point" (§6.7), a transparent loop. Returns (psi_new, max_change). Book: §6.7, (6.72)–(6.73), Example 6.2.
    Label: analytic."""
    return sor_sweep(psi, mask, bc, 1.0, dx, dy, f, order)


def sor_sweep(psi, mask, bc=None, omega: float = 1.5, dx: float = 1.0, dy: float | None = None, f=None,
              order: str = "lex"):
    """One successive over-relaxation sweep ψ ← ψ + ω(ψ_GS − ψ), 0 < ω < 2 (ω = 1: Gauss–Seidel). Our addition to the
    book's method (labelled). Returns (psi_new, max_change). Label: analytic."""
    if not 0.0 < float(omega) < 2.0:
        raise ValueError("SOR needs 0 < omega < 2")
    P = _with_bc(psi, mask, bc)
    cx, cy, cc = _weights(dx, dy)
    Fs = None if f is None else _F(f)
    ch = 0.0
    for j, i in sweep_order(mask, order):
        src = 0.0 if Fs is None else Fs[j, i]
        gs = (cx * (P[j, i - 1] + P[j, i + 1]) + cy * (P[j - 1, i] + P[j + 1, i]) - src) / cc  # Eq. (6.72)
        new = P[j, i] + omega * (gs - P[j, i])
        ch = max(ch, abs(new - P[j, i]))
        P[j, i] = new
    return P, float(ch)


def residual_field(psi, mask, dx: float = 1.0, dy: float | None = None, f=None):
    """Defect of the average rule at each unknown node, r = avg − ψ (0 elsewhere) [unit of ψ]; equals
    Δx²∇²_hψ/4 for Δx = Δy (Laplace). Label: analytic."""
    P = _F(psi)
    M = np.asarray(mask, bool)
    cx, cy, cc = _weights(dx, dy)
    src = 0.0 if f is None else _F(f)
    avg = np.zeros_like(P)
    avg[1:-1, 1:-1] = cx * (P[1:-1, :-2] + P[1:-1, 2:]) + cy * (P[:-2, 1:-1] + P[2:, 1:-1])
    r = (avg - src) / cc - P
    return np.where(M, r, 0.0)


def residual_norm(psi, mask, dx: float = 1.0, dy: float | None = None, f=None) -> float:
    """max |defect of the average rule| over the unknown nodes = max |ψ_{i−1,j} + ψ_{i+1,j} + ψ_{i,j−1} + ψ_{i,j+1} −
    4ψ_{i,j}|/4 for Δx = Δy (the stopping criterion of :func:`solve_laplace`). Label: analytic."""
    return float(np.max(np.abs(residual_field(psi, mask, dx, dy, f))))


def assemble_laplace(mask, bc, dx: float = 1.0, dy: float | None = None, f=None, order: str = "book"):
    """Sparse system of the average rule for the unknowns in sweep order: ψ_p − Σ c_q ψ_q = b_p, i.e. A = I − N with
    N the neighbour weights (c = (1/Δx²)/(2/Δx² + 2/Δy²), …) and b the known neighbours' contribution (and −f/…).

    Returns dict(A (csr), N (csr), b, nodes (list of (j, i)), index (array −1 or unknown number)). Label: analytic."""
    M = np.asarray(mask, bool)
    B = _F(bc)
    cx, cy, cc = _weights(dx, dy)
    nodes = sweep_order(M, order)
    idx = -np.ones(M.shape, dtype=int)
    for k, (j, i) in enumerate(nodes):
        idx[j, i] = k
    rows, cols, vals = [], [], []
    b = np.zeros(len(nodes))
    Fs = None if f is None else _F(f)
    for k, (j, i) in enumerate(nodes):
        for (jj, ii, c) in ((j, i - 1, cx), (j, i + 1, cx), (j - 1, i, cy), (j + 1, i, cy)):
            if idx[jj, ii] >= 0:
                rows.append(k)
                cols.append(idx[jj, ii])
                vals.append(c / cc)
            else:
                b[k] += c / cc * B[jj, ii]
        if Fs is not None:
            b[k] -= Fs[j, i] / cc
    n = len(nodes)
    N = sps.csr_matrix((vals, (rows, cols)), shape=(n, n))
    A = (sps.identity(n, format="csr") - N).tocsr()
    return {"A": A, "N": N, "b": b, "nodes": nodes, "index": idx}


def optimal_sor_omega(nx: int, ny: int) -> float:
    """ω_opt = 2/(1 + √(1 − ρ_J²)) for the model problem on an nx × ny node rectangle (ρ_J from
    :func:`jacobi_spectral_radius`) — our addition (standard SOR theory). Label: analytic."""
    rJ = jacobi_spectral_radius(nx, ny)
    return float(2.0 / (1.0 + np.sqrt(1.0 - rJ ** 2)))


def jacobi_spectral_radius(nx: int, ny: int) -> float:
    """ρ_J = ½[cos(π/(nx − 1)) + cos(π/(ny − 1))] for Δx = Δy on an nx × ny node rectangle (Dirichlet frame); Gauss–Seidel
    converges as ρ_J² per sweep. Our addition. Label: analytic."""
    return float(0.5 * (np.cos(np.pi / (nx - 1)) + np.cos(np.pi / (ny - 1))))


def solve_laplace(mask, bc, method: str = "gauss_seidel", tol: float = 1e-10, max_iter: int = 100_000,
                  omega: float | None = None, dx: float = 1.0, dy: float | None = None, psi0=None,
                  n_iter: int | None = None, order: str = "book", f=None):
    """Solve the discrete Laplace (or Poisson ∇²ψ = f) equation on a masked grid.

    Parameters
    ----------
    mask : bool array [j, i], True at unknown nodes.   bc : array [j, i] with the fixed (Dirichlet) values elsewhere.
    method : "gauss_seidel" (the book's), "jacobi", "sor" (ω, default :func:`optimal_sor_omega` of the bounding box) or
        "direct" (sparse LU, ``scipy.sparse.linalg.spsolve``).
    tol : stop when the max-norm defect of the average rule ≤ tol [unit of ψ] (the residual, not the last change).
    n_iter : run exactly this many sweeps (the book's fixed iteration count) instead of the tolerance test.
    dx, dy : spacings [m];  psi0 : initial guess at the unknowns (default 0, the book's);  order : sweep order.
    f : optional source for Poisson (same layout).

    Returns
    -------
    (psi, history) — psi [j, i] with the boundary values in place; history dict(residual (list, per sweep, starting
    with the initial guess; empty for "direct"), change (max |Δψ| per sweep), sweeps (= iterations; 0 for "direct"),
    iterations, converged, method, omega, final_residual).

    Iterations are exact matrix forms of the sweeps: Jacobi ψ ← b + Nψ; Gauss–Seidel (I − L)ψ⁺ = b + Uψ; SOR
    (I − ωL)ψ⁺ = ω(b + Uψ) + (1 − ω)ψ (L, U: neighbours before/after in the sweep order) — identical iterates to
    :func:`gauss_seidel_sweep` and :func:`sor_sweep`, fast on large grids.
    Book: §6.7, (6.72)–(6.73), Example 6.2. Validation (planned): V1 4-point system; exact on discrete-harmonic xy;
    all methods agree with "direct" to tol; V3 grid refinement. Label: analytic, converged.
    """
    M = np.asarray(mask, bool)
    B = _F(bc).copy()
    sysd = assemble_laplace(M, B, dx, dy, f, order)
    A, N, b, nodes = sysd["A"], sysd["N"], sysd["b"], sysd["nodes"]
    Jn = np.array([p[0] for p in nodes], dtype=int)
    In = np.array([p[1] for p in nodes], dtype=int)
    n = len(nodes)

    def to_grid(x):
        P = B.copy()
        P[Jn, In] = x
        return P

    x = np.zeros(n) if psi0 is None else (_F(psi0)[Jn, In].copy() if np.ndim(psi0) == 2 else np.full(n, float(psi0)))
    hist = {"residual": [], "change": [], "iterations": 0, "sweeps": 0, "converged": False, "method": method,
            "omega": None}

    def res(xv):
        return float(np.max(np.abs(b + N @ xv - xv))) if n else 0.0

    if method == "direct":
        x = spsolve(A.tocsc(), b) if n else x
        hist["final_residual"] = res(x)
        hist["iterations"] = hist["sweeps"] = 0
        hist["converged"] = True
        return to_grid(x), hist
    L = sps.tril(N, k=-1, format="csr")
    Up = sps.triu(N, k=1, format="csr")
    I_ = sps.identity(n, format="csr")
    if method == "jacobi":
        w = None
    elif method in ("gauss_seidel", "sor"):
        w = 1.0 if method == "gauss_seidel" else (float(omega) if omega is not None
                                                  else optimal_sor_omega(M.shape[1], M.shape[0]))
        if not 0.0 < w < 2.0:
            raise ValueError("SOR needs 0 < omega < 2")
        Lw = (I_ - w * L).tocsr()
    else:
        raise ValueError('method must be "gauss_seidel", "jacobi", "sor" or "direct"')
    hist["omega"] = w
    hist["residual"].append(res(x))
    total = int(n_iter) if n_iter is not None else int(max_iter)
    for k in range(total):
        if method == "jacobi":
            xn = b + N @ x
        else:
            rhs = w * (b + Up @ x) + (1.0 - w) * x
            xn = spsolve_triangular(Lw, rhs, lower=True, unit_diagonal=True)
        hist["change"].append(float(np.max(np.abs(xn - x))) if n else 0.0)
        x = xn
        r = res(x)
        hist["residual"].append(r)
        hist["iterations"] = k + 1
        if n_iter is None and r <= tol:
            hist["converged"] = True
            break
    if n_iter is not None:
        hist["converged"] = hist["residual"][-1] <= tol
    hist["sweeps"] = hist["iterations"]
    hist["final_residual"] = hist["residual"][-1]
    return to_grid(x), hist


def solve_poisson(mask, f, bc, method: str = "direct", **kw):
    """Discrete Poisson ∇²ψ = f on a masked grid (e.g. ω = −∇²ψ (6.4) with f = −ω) — :func:`solve_laplace` with a
    source. Returns (psi, history). Book: §6.2 (6.4), (6.11); §6.7. Label: analytic."""
    return solve_laplace(mask, bc, method=method, f=f, **kw)
