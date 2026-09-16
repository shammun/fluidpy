"""Finite-difference field operators ∇, ∇·, ∇× (and the velocity gradient, tensor divergence, Laplacian) on the
project grid layout of :mod:`fluidpy.core.grids`.

Book: Ch. 2 §2.9, Eqs. (2.22)–(2.25); §2.14 comma notation (2.36). Every stencil is written out explicitly
(second-order central in the interior, second-order one-sided at the edges — never ``np.gradient``, whose edges are
first order) so the observed order can be measured (V3).

Index conventions fixed here (see ``knowledge/notation.md``)
------------------------------------------------------------
* Fields follow ``core.grids``: scalar ``phi[k, j, i]`` (or ``[j, i]``), vector ``u[c, ...]`` with the component on
  axis 0, tensor ``T[a, b, ...]``. Direction d ∈ {0, 1, 2} = (x₁, x₂, x₃) is array axis ``ndim − 1 − d``.
* **Velocity gradient** ``vector_gradient(u)[i, j] = ∂u_i/∂x_j`` (Kundu Ch. 3: S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i)).
* **Tensor divergence** contracts the **second** index by default: ``(∇·τ)_i = ∂τ_ij/∂x_j`` (§2.9) — the opposite
  slot from Cauchy's traction f_i = τ_ji n_j; they agree only for symmetric τ. ``index=0`` gives ∂τ_ji/∂x_j.
* ``h`` is a scalar spacing or a per-direction sequence ``(hx, hy[, hz])`` [m]; ``bc="onesided"`` (default) or
  ``"periodic"`` (grid built with ``periodic=True``).

Units: derivatives carry [field unit / m].
"""
from __future__ import annotations

import numpy as np

from .grids import axis_of_direction, spacing
from .tensors import levi_civita

__all__ = ["partial", "second_partial", "gradient", "directional_derivative", "divergence", "vector_gradient",
           "tensor_divergence", "curl", "curl_components", "laplacian", "is_solenoidal", "is_irrotational"]


def _sl(ndim: int, axis: int, a, b):
    """Slice ``[a:b]`` along ``axis`` of an ``ndim``-array (all other axes full)."""
    idx = [slice(None)] * ndim
    idx[axis] = slice(a, b)
    return tuple(idx)


def partial(f, direction: int, h, bc: str = "onesided", ndim_space: int | None = None) -> np.ndarray:
    """∂f/∂x_d by explicit second-order finite differences along coordinate direction ``d``.

    Book: §2.9, Eq. (2.22) (∇ = e_i ∂/∂x_i — the building block of every operator here); §2.14 comma notation
    A_,i = ∂A/∂x_i (2.36).

    Interior: ``(f_{i+1} − f_{i−1})/(2h)``. Edges (``bc="onesided"``): ``(−3f_0 + 4f_1 − f_2)/(2h)`` and its mirror —
    second order, so the whole field converges at order 2. ``bc="periodic"``: the central stencil wraps around
    (grid built with ``periodic=True``).

    Parameters
    ----------
    f : array_like
        Field on the project grid; leading component axes are allowed (``u[c, k, j, i]``).
    direction : int
        0 = x₁, 1 = x₂, 2 = x₃.
    h : float or sequence
        Spacing [m] (scalar or per direction).
    bc : {"onesided", "periodic"}
    ndim_space : int, optional
        Number of spatial axes (default: all axes of ``f``; pass it when ``f`` carries leading component axes and
        ``direction`` must be counted from the last axis — the default already does this, so it is rarely needed).

    Returns
    -------
    ndarray, same shape as ``f`` [field unit / m].

    Validation: V1 exact (1e-12) on quadratics including the edges; V3 observed order 2.00 ± 0.05 on sin/cos fields
    (n = 16 … 128) for both boundary treatments. Label: analytic, converged.
    """
    F = np.asarray(f, dtype=float)
    nd = F.ndim if ndim_space is None else int(ndim_space)
    axis = F.ndim - nd + axis_of_direction(nd, direction)
    hd = spacing(h, direction)
    n = F.shape[axis]
    if n < 3:
        raise ValueError("need at least 3 nodes along the differentiation axis")
    d = np.empty_like(F)
    if bc == "periodic":
        d[...] = (np.roll(F, -1, axis) - np.roll(F, 1, axis)) / (2.0 * hd)  # central, wrapped
        return d
    if bc != "onesided":
        raise ValueError('bc must be "onesided" or "periodic"')
    S = lambda a, b: _sl(F.ndim, axis, a, b)  # noqa: E731
    d[S(1, -1)] = (F[S(2, None)] - F[S(None, -2)]) / (2.0 * hd)  # central difference, order 2
    d[S(0, 1)] = (-3.0 * F[S(0, 1)] + 4.0 * F[S(1, 2)] - F[S(2, 3)]) / (2.0 * hd)  # one-sided, order 2
    d[S(-1, None)] = (3.0 * F[S(-1, None)] - 4.0 * F[S(-2, -1)] + F[S(-3, -2)]) / (2.0 * hd)
    return d


def second_partial(f, direction: int, h, bc: str = "onesided") -> np.ndarray:
    """∂²f/∂x_d² by explicit second-order stencils (interior ``(f_{i+1} − 2f_i + f_{i−1})/h²``, edges
    ``(2f_0 − 5f_1 + 4f_2 − f_3)/h²``). Tool for :func:`laplacian` (Ch. 4+). Label: analytic, converged."""
    F = np.asarray(f, dtype=float)
    axis = axis_of_direction(F.ndim, direction)
    hd = spacing(h, direction)
    if F.shape[axis] < 4:
        raise ValueError("need at least 4 nodes along the axis for the second-order edge stencil")
    d = np.empty_like(F)
    if bc == "periodic":
        d[...] = (np.roll(F, -1, axis) - 2.0 * F + np.roll(F, 1, axis)) / hd ** 2
        return d
    S = lambda a, b: _sl(F.ndim, axis, a, b)  # noqa: E731
    d[S(1, -1)] = (F[S(2, None)] - 2.0 * F[S(1, -1)] + F[S(None, -2)]) / hd ** 2
    d[S(0, 1)] = (2.0 * F[S(0, 1)] - 5.0 * F[S(1, 2)] + 4.0 * F[S(2, 3)] - F[S(3, 4)]) / hd ** 2
    d[S(-1, None)] = (2.0 * F[S(-1, None)] - 5.0 * F[S(-2, -1)] + 4.0 * F[S(-3, -2)] - F[S(-4, -3)]) / hd ** 2
    return d


def gradient(phi, h, bc: str = "onesided") -> np.ndarray:
    """Gradient of a scalar field: (∇φ)_i = ∂φ/∂x_i, stacked on axis 0.

    Book: §2.9 (unnumbered, after (2.22)): ∇φ = e_i ∂φ/∂x_i; ∇φ is perpendicular to the surfaces φ = const and
    points along the fastest increase.

    Parameters
    ----------
    phi : array_like, shape ``(ny, nx)`` or ``(nz, ny, nx)``  [φ unit]
    h : spacing [m]; bc : boundary treatment

    Returns
    -------
    grad : ndarray, shape ``(d, *phi.shape)``  [φ unit / m]

    Validation: V1 exact on quadratics; ∇φ·t = 0 for t tangent to a level set; V3 order 2; ∇×∇φ = 0 to truncation
    (Exercise 2.20); V7 rotational invariance. Label: analytic, converged.
    """
    P = np.asarray(phi, dtype=float)
    return np.stack([partial(P, d, h, bc) for d in range(P.ndim)])  # (∇φ)_i = ∂φ/∂x_i


def directional_derivative(grad_phi, n):
    """∂φ/∂n = ∇φ·n — the rate of change of φ along the unit direction n.

    Book: §2.9 (unnumbered): "The spatial rate of change of φ in any other direction n is given by ∂φ/∂n = ∇φ·n";
    maximal when n ∥ ∇φ (Fig. 2.7).

    Parameters
    ----------
    grad_phi : array_like, shape ``(d, ...)`` — output of :func:`gradient`
    n : array_like, shape ``(d,)`` or ``(d, ...)`` — direction (normalised internally)

    Returns
    -------
    ndarray, shape ``grad_phi.shape[1:]``  [φ unit / m]

    Validation: V1 zero along a level set, |∇φ| along ∇φ; V7 |∂φ/∂n| ≤ |∇φ| for random n. Label: analytic.
    """
    G = np.asarray(grad_phi, dtype=float)
    N = np.asarray(n, dtype=float)
    N = N / np.linalg.norm(N, axis=0, keepdims=True)
    if N.ndim == 1:
        return np.einsum("i...,i->...", G, N)  # ∇φ·n
    return np.einsum("i...,i...->...", G, N)


def divergence(u, h, bc: str = "onesided") -> np.ndarray:
    """Divergence of a vector field ∇·u = ∂u_i/∂x_i (= u_i,i in comma notation).

    Book: §2.9, Eq. (2.23); §2.14, Eq. (2.36).

    Parameters
    ----------
    u : array_like, shape ``(d, *grid)`` — components on axis 0 (d = 2 or 3 matching the grid)  [u unit]
    h : spacing [m]; bc : boundary treatment

    Returns
    -------
    ndarray, shape ``grid``  [u unit / m]

    Validation: V1 Example 2.3 ∇·(a x) = 3a exactly (linear field), ∇·(b × x) = 0; == trace(``vector_gradient``);
    ∇·(∇×u) = 0 to truncation (Exercise 2.19); V3 order 2; V7 invariance under rotation. Label: analytic, converged.
    """
    U = np.asarray(u, dtype=float)
    d = U.shape[0]
    if U.ndim - 1 != d:
        raise ValueError(f"u has {d} components but {U.ndim - 1} spatial axes; use a matching grid (2-D or 3-D)")
    return sum(partial(U[i], i, h, bc) for i in range(d))  # Eq. (2.23): ∂u_i/∂x_i


def vector_gradient(u, h, bc: str = "onesided") -> np.ndarray:
    """Velocity-gradient tensor G[i, j, ...] = ∂u_i/∂x_j (gradient raises the order by one).

    Book: §2.9 (unnumbered): "the gradient operation increases the order of a tensor by one … i.e., ∂u_i/∂x_j";
    §2.4 lists ∂u_i/∂x_j among the second-order tensors. Index order fixed here for Ch. 3 (S_ij, R_ij).

    Parameters
    ----------
    u : array_like, shape ``(d, *grid)``; a 3-component field on a 2-D grid is allowed (∂/∂x₃ = 0 → G is 3 × 3).

    Returns
    -------
    G : ndarray, shape ``(nc, nc, *grid)`` with ``G[i, j] = ∂u_i/∂x_j``  [u unit / m]

    Validation: V1 exact on linear fields (G = b-cross matrix for u = b × x); V7 transforms per (2.12) under rotation;
    ``antisymmetric_part(G)`` ↔ ½∇×u (sign pinned). Label: analytic.
    """
    U = np.asarray(u, dtype=float)
    nc, nd = U.shape[0], U.ndim - 1
    G = np.zeros((nc, nc) + U.shape[1:])
    for i in range(nc):
        for j in range(min(nc, nd)):
            G[i, j] = partial(U[i], j, h, bc)  # G_ij = ∂u_i/∂x_j (a 3rd column stays 0 on a 2-D grid)
    return G


def tensor_divergence(T, h, bc: str = "onesided", index: int = 1) -> np.ndarray:
    """Divergence of a second-order tensor field: (∇·τ)_i = ∂τ_ij/∂x_j (``index=1``, the book) or ∂τ_ji/∂x_j (``index=0``).

    Book: §2.9 (unnumbered): "(∇·τ)_i = Σ_j ∂τ_ij/∂x_j … the divergence operation decreases the order by one".
    The contracted slot is the **second** index — opposite to Cauchy's traction (2.15); identical for symmetric τ.

    Parameters
    ----------
    T : array_like, shape ``(d, d, *grid)``  [τ unit]
    index : {1, 0} — which index the derivative contracts.

    Returns
    -------
    ndarray, shape ``(d, *grid)``  [τ unit / m]

    Validation: V1 ``tensor_divergence(outer-product field)`` vs sympy; ``index=1`` vs ``index=0`` differ for a
    non-symmetric field and agree for a symmetric one (discrimination); V3 order 2. Label: analytic, converged.
    """
    A = np.asarray(T, dtype=float)
    d = A.shape[0]
    out = np.zeros((d,) + A.shape[2:])
    for i in range(d):
        for j in range(min(d, A.ndim - 2)):
            out[i] += partial(A[i, j], j, h, bc) if index == 1 else partial(A[j, i], j, h, bc)  # ∂τ_ij/∂x_j
    return out


def curl(u, h, bc: str = "onesided") -> np.ndarray:
    """Curl of a vector field in index form: (∇×u)_i = ε_ijk ∂u_k/∂x_j.

    Book: §2.9, Eq. (2.24) (from (2.21) and (2.22)); comma form ε_ijk u_k,j (2.36).

    Parameters
    ----------
    u : array_like
        ``(3, nz, ny, nx)`` → 3-component curl; ``(3, ny, nx)`` (z-independent field) → 3-component curl with
        ∂/∂x₃ = 0; ``(2, ny, nx)`` (plane field) → the scalar (∇×u)₃ = ∂u₂/∂x₁ − ∂u₁/∂x₂.

    Returns
    -------
    ndarray  [u unit / m]

    Validation: V1 == ``curl_components`` (2.25) to 1e-14; Example 2.3 ∇×(a x) = 0, ∇×(b × x) = 2b exactly (linear);
    ∇×∇φ = 0 to truncation; V3 order 2; V7 rotational invariance (as a vector, proper rotations).
    Label: analytic, converged.
    """
    U = np.asarray(u, dtype=float)
    if U.shape[0] == 2 and U.ndim == 3:
        return partial(U[1], 0, h, bc) - partial(U[0], 1, h, bc)  # (∇×u)₃ = ∂u₂/∂x₁ − ∂u₁/∂x₂  (2.25)
    G = vector_gradient(U, h, bc)  # G[k, j] = ∂u_k/∂x_j
    return np.einsum("ijk,kj...->i...", levi_civita(), G)  # Eq. (2.24): ε_ijk ∂u_k/∂x_j


def curl_components(u, h, bc: str = "onesided") -> np.ndarray:
    """The three curl components written out: (∇×u)₁ = ∂u₃/∂x₂ − ∂u₂/∂x₃, (∇×u)₂ = ∂u₁/∂x₃ − ∂u₃/∂x₁,
    (∇×u)₃ = ∂u₂/∂x₁ − ∂u₁/∂x₂.

    Book: §2.9, Eq. (2.25) (D12). Accepts the same shapes as :func:`curl` (3 components required).

    Validation: V1 equals :func:`curl` to machine precision. Label: analytic.
    """
    U = np.asarray(u, dtype=float)
    if U.shape[0] != 3:
        raise ValueError("curl_components needs a 3-component field (use curl for the plane scalar case)")
    nd = U.ndim - 1
    P = lambda comp, d: partial(U[comp], d, h, bc) if d < nd else np.zeros(U.shape[1:])  # noqa: E731  ∂/∂x₃ = 0 on a plane
    c1 = P(2, 1) - P(1, 2)  # Eq. (2.25): ∂u₃/∂x₂ − ∂u₂/∂x₃
    c2 = P(0, 2) - P(2, 0)  #             ∂u₁/∂x₃ − ∂u₃/∂x₁
    c3 = P(1, 0) - P(0, 1)  #             ∂u₂/∂x₁ − ∂u₁/∂x₂
    return np.stack([c1, c2, c3])


def laplacian(phi, h, bc: str = "onesided") -> np.ndarray:
    """∇²φ = ∂²φ/∂x_i∂x_i by explicit second-order stencils (tool for Ch. 4+; ∇·∇φ of §2.9).

    Validation: V1 ∇²(r²) = 2d exactly; V3 order 2 on sin/cos. Label: analytic, converged.
    """
    P = np.asarray(phi, dtype=float)
    return sum(second_partial(P, d, h, bc) for d in range(P.ndim))


def is_solenoidal(u, h, tol: float = 1e-8, bc: str = "onesided"):
    """Divergence-free test: ``(max |∇·u| ≤ tol, max |∇·u|)``. Book: §2.9 ("solenoidal or divergence free if ∇·u = 0")."""
    r = float(np.max(np.abs(divergence(u, h, bc))))
    return bool(r <= tol), r


def is_irrotational(u, h, tol: float = 1e-8, bc: str = "onesided"):
    """Curl-free test: ``(max |∇×u| ≤ tol, max |∇×u|)``. Book: §2.9 ("irrotational or curl free if ∇×u = 0")."""
    r = float(np.max(np.abs(curl(u, h, bc))))
    return bool(r <= tol), r
