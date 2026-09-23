"""Gauss' and Stokes' theorems as computations: volume/surface/line integrals on boxes, spheres, discs and loops,
and the integral definitions of ∇, ∇· and ∇× as small-volume / small-area limits.

Book: Ch. 2 §2.12, Eqs. (2.30)–(2.33) and Example 2.5; §2.13, Eqs. (2.34)–(2.35) and Example 2.6. Read from the
rendered pages chapters/pages/ch02/p085–p089.

Field callables follow ``core.grids``/``core.fields``: ``Q_fn(X, Y, Z)`` (or ``Q_fn(X, Y)`` in 2-D) returns a scalar
array or an array with the component on axis 0, broadcast over the coordinate arrays. Exact derivatives may be
supplied (``div_fn``, ``curl_fn``, ``dQ_fn``, e.g. from :class:`fluidpy.core.fields.VectorField`); otherwise they are
formed by a fourth-order central difference of the callable (step ``fd_step``), which is exact for polynomials of
degree ≤ 4 and accurate to ~1e-11 otherwise — this never enters the *surface* sides of the theorems.

Quadrature: composite **midpoint** rule (order 2; exact for linear integrands, which is what Example 2.5 uses) or
**Simpson** (``rule="simpson"``, order 4, needs an odd node count). Orientation: outward n on closed surfaces; on an
open surface the chosen n fixes the boundary tangent t = n_c × n, where n_c is the in-surface normal to C pointing
*into* A (Fig. 2.10) — t then runs counterclockwise seen from the outside.

Units: coordinates in m; integrals carry [Q unit × m³] (volume), [Q unit × m²] (surface), [u unit × m] (line);
the point definitions (2.31)–(2.33), (2.35) return [field unit / m].
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, NamedTuple

import numpy as np

from ._util import as_scalar_if_0d
from .tensors import levi_civita

__all__ = ["midpoint_nodes", "simpson_nodes", "gauss_legendre_nodes", "volume_integral_box", "flux_through_box_faces", "gauss_gradient_box",
           "divergence_theorem_box", "flux_through_sphere", "divergence_theorem_sphere", "integral_gradient",
           "integral_divergence", "integral_curl", "flux_through_faces", "divergence_theorem_rect2d",
           "integral_divergence_2d", "TiledFlux", "divergence_theorem_tiled", "Loop", "Surface", "boundary_tangent", "planar_loop", "rectangle_loop",
           "planar_disc", "planar_rectangle", "circulation", "curl_flux", "StokesCheck", "stokes_theorem_check",
           "integral_curl_component", "fd_partials", "CurlCheck", "curl_theorem_box"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731


# ----------------------------------------------------------------------------------------------------------------------
# 1-D quadrature building blocks
# ----------------------------------------------------------------------------------------------------------------------
def midpoint_nodes(lo: float, hi: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Cell centres and equal weights of the composite midpoint rule on [lo, hi] with n cells (order 2)."""
    h = (hi - lo) / n
    return lo + h * (np.arange(n) + 0.5), np.full(n, h)


def simpson_nodes(lo: float, hi: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Nodes and weights of composite Simpson on [lo, hi] with n nodes (n odd ≥ 3; order 4)."""
    if n < 3 or n % 2 == 0:
        raise ValueError("Simpson needs an odd number of nodes >= 3")
    x, h = np.linspace(lo, hi, n, retstep=True)
    w = np.ones(n)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return x, w * h / 3.0


def gauss_legendre_nodes(lo: float, hi: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss–Legendre nodes and weights on [lo, hi] (exact for polynomials of degree ≤ 2n − 1; used on the sphere)."""
    x, w = np.polynomial.legendre.leggauss(int(n))
    return 0.5 * (hi - lo) * x + 0.5 * (hi + lo), 0.5 * (hi - lo) * w


def _nodes(lo, hi, n, rule):
    if rule == "midpoint":
        return midpoint_nodes(lo, hi, n)
    if rule == "simpson":
        return simpson_nodes(lo, hi, n)
    if rule == "gauss":
        return gauss_legendre_nodes(lo, hi, n)
    raise ValueError('rule must be "midpoint", "simpson" or "gauss"')


def _bounds3(bounds):
    b = _F(bounds)
    if b.shape == (2,):
        b = np.tile(b, (3, 1))
    if b.shape != (3, 2):
        raise ValueError("bounds must be (lo, hi) or ((x0, x1), (y0, y1), (z0, z1))")
    return b


def _bounds2(bounds):
    b = _F(bounds)
    if b.shape == (2,):
        b = np.tile(b, (2, 1))
    if b.shape != (2, 2):
        raise ValueError("bounds must be (lo, hi) or ((x0, x1), (y0, y1))")
    return b


def _eval(fn: Callable, *coords) -> np.ndarray:
    """Evaluate a field callable and broadcast constant components to the coordinate shape.

    The callable may return a scalar array (shape = points), a list/tuple of components (each a scalar array or a
    constant), or an array with any number of **leading component axes** (vector ``(3, *pts)``, tensor ``(3, 3, *pts)``,
    e.g. ``fd_partials`` of a vector field). Components that are constants are broadcast to the point shape.
    """
    B = np.broadcast_arrays(*[_F(c) for c in coords])
    pts = B[0].shape
    out = fn(*B)
    if isinstance(out, (list, tuple)):
        return np.stack([np.broadcast_to(_F(c), pts) for c in out])
    out = _F(out)
    if out.ndim <= B[0].ndim:
        return np.broadcast_to(out, pts).copy()
    lead = out.shape[: out.ndim - B[0].ndim]  # leading component axes (1 for a vector, 2 for a tensor, …)
    return np.broadcast_to(out, lead + pts).copy()


def fd_partials(fn: Callable, coords, fd_step: float = 1e-3) -> np.ndarray:
    """∂Q/∂x_i of a callable field by the fourth-order central difference (−Q(+2h) + 8Q(+h) − 8Q(−h) + Q(−2h))/(12h).

    Returns ``dQ[i, ...]`` = ∂Q/∂x_i with ``...`` the shape of ``Q`` at the points (scalar Q → (d, *pts); vector Q →
    (d, ncomp, *pts)). Tool for the volume side of (2.30) when no exact derivative is supplied. Label: analytic.
    """
    coords = [_F(c) for c in coords]
    out = []
    for i in range(len(coords)):
        def shifted(s):
            cs = list(coords)
            cs[i] = coords[i] + s * fd_step
            return _eval(fn, *cs)
        out.append((-shifted(2.0) + 8.0 * shifted(1.0) - 8.0 * shifted(-1.0) + shifted(-2.0)) / (12.0 * fd_step))
    return np.stack(out)


# ----------------------------------------------------------------------------------------------------------------------
# §2.12 Gauss' theorem on a box and a sphere
# ----------------------------------------------------------------------------------------------------------------------
def volume_integral_box(f_fn: Callable, bounds, n: int = 24, rule: str = "midpoint"):
    """∭_V f dV over a rectangular box by a product quadrature rule (f scalar or with components on axis 0).

    Book: §2.12, the left side of (2.30). ``n`` nodes per direction (FAST: 16, full: 32 midpoint / 33 Simpson).
    Label: converged (order 2 midpoint, 4 Simpson).
    """
    b = _bounds3(bounds)
    (x, wx), (y, wy), (z, wz) = (_nodes(*b[k], n, rule) for k in range(3))
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")  # [k, j, i] project layout
    W = np.einsum("k,j,i->kji", wz, wy, wx)
    F = _eval(f_fn, X, Y, Z)
    return as_scalar_if_0d(np.einsum("...kji,kji->...", F, W))


def flux_through_box_faces(Q_fn: Callable, bounds, n: int = 24, rule: str = "midpoint") -> dict[str, float]:
    """Outward flux ∬ n·Q dA through each of the six faces of a box (Fig. 2.4 lettering: EADH = +x₁ … DCGH = −x₃).

    Book: §2.12 (right side of the divergence theorem) and Example 2.5 (face-by-face bookkeeping).

    Returns
    -------
    dict with keys ``"+x", "-x", "+y", "-y", "+z", "-z"`` → flux [Q unit × m²]; the sum is ∯ n·Q dA.
    """
    b = _bounds3(bounds)
    out = {}
    names = "xyz"
    for d in range(3):
        others = [k for k in range(3) if k != d]
        (s, ws), (t, wt) = (_nodes(*b[k], n, rule) for k in others)
        T, S = np.meshgrid(t, s, indexing="ij")
        W = np.einsum("j,i->ji", wt, ws)
        for sign, key in ((+1, "+"), (-1, "-")):
            coords = [None, None, None]
            coords[d] = np.full_like(S, b[d, 1] if sign > 0 else b[d, 0])
            coords[others[0]], coords[others[1]] = S, T
            Q = _eval(Q_fn, *coords)
            out[f"{key}{names[d]}"] = float(np.sum(sign * Q[d] * W))  # n·Q with n = ±e_d
    return out


def gauss_gradient_box(Q_fn: Callable, bounds, n: int = 24, dQ_fn: Callable | None = None, rule: str = "midpoint",
                       fd_step: float = 1e-3):
    """Both sides of Gauss' theorem for a scalar or vector field Q on a box: ∭ ∂Q/∂x_i dV = ∬ n_i Q dA.

    Book: §2.12, Eq. (2.30) (a field "of any order"; the surface side pairs the *free* index of n with Q).

    Parameters
    ----------
    Q_fn : callable ``Q(X, Y, Z)`` → scalar array or ``(3, ...)`` array
    bounds : (lo, hi) or three pairs  [m]
    n : nodes per direction
    dQ_fn : callable returning ``dQ[i, ...] = ∂Q/∂x_i`` (exact), or None → :func:`fd_partials`
    rule : {"midpoint", "simpson"}

    Returns
    -------
    (lhs, rhs) : arrays, shape ``(3,)`` for scalar Q (index i) or ``(3, 3)`` for vector Q (``[i, j]`` = ∂Q_j/∂x_i ↔
        n_i Q_j — the derivative index first, as in (2.30)/(2.31); this is the **transpose** of
        ``core.operators.vector_gradient``'s G[i, j] = ∂u_i/∂x_j). A tensor Q ``(3, 3, ...)`` gives ``(3, 3, 3)`` likewise.

    Validation: V1 Q = xyz on the unit cube: lhs == rhs (1e-12, midpoint is exact for the linear face integrands after
    cancellation… checked numerically); vector Q = (xy, yz, zx): lhs == rhs, lhs[i, j] = ∫∂Q_j/∂x_i dV (diag ½,
    [1, 0] = ½, [0, 1] = 0); V3 order 2 (midpoint) / 4 (Simpson) for a sin/cos field. Label: converged.
    """
    b = _bounds3(bounds)
    if dQ_fn is None:
        dfn = lambda X, Y, Z: fd_partials(Q_fn, (X, Y, Z), fd_step)  # noqa: E731
    else:
        dfn = dQ_fn
    lhs = _F(volume_integral_box(dfn, b, n, rule))  # ∭ ∂Q/∂x_i dV
    # surface side: Σ_faces n_i Q dA
    rhs = None
    for d in range(3):
        others = [k for k in range(3) if k != d]
        (s, ws), (t, wt) = (_nodes(*b[k], n, rule) for k in others)
        T, S = np.meshgrid(t, s, indexing="ij")
        W = np.einsum("j,i->ji", wt, ws)
        for sign in (+1, -1):
            coords = [None, None, None]
            coords[d] = np.full_like(S, b[d, 1] if sign > 0 else b[d, 0])
            coords[others[0]], coords[others[1]] = S, T
            Q = _eval(Q_fn, *coords)
            face = sign * np.einsum("...ji,ji->...", Q, W)  # n_d Q dA summed over the face
            if rhs is None:
                rhs = np.zeros((3,) + np.shape(face))
            rhs[d] += face  # Eq. (2.30): ∬ n_i Q dA
    return lhs, rhs


def divergence_theorem_box(Q_fn: Callable, bounds, n: int = 24, div_fn: Callable | None = None, rule: str = "midpoint",
                           fd_step: float = 1e-3):
    """Divergence theorem on a box: ∭_V ∇·Q dV (lhs) and ∯_A n·Q dA (rhs).

    Book: §2.12, the vector form of (2.30): "the volume integral of the divergence of Q is equal to the surface
    integral of the outflux of Q".

    Parameters
    ----------
    Q_fn : callable ``Q(X, Y, Z) → (3, ...)``
    div_fn : exact ``∇·Q(X, Y, Z)`` or None (fourth-order FD of ``Q_fn``)

    Returns
    -------
    (lhs, rhs) : floats [Q unit × m²]

    Validation: V1 Q = (x, y, z) on [0,1]³: 3 == 3; Q = (x², 0, 0): 1 == 1; V4 solenoidal b × x: rhs = 0; V3 order 2 /
    4 for a smooth field; V5 sphere version reproduces 8π/3 (Wikipedia). Label: analytic, converged.
    """
    if div_fn is None:
        def div_fn(X, Y, Z):
            d = fd_partials(Q_fn, (X, Y, Z), fd_step)  # d[i, c, ...] = ∂Q_c/∂x_i
            return d[0, 0] + d[1, 1] + d[2, 2]
    lhs = float(volume_integral_box(div_fn, bounds, n, rule))
    rhs = float(sum(flux_through_box_faces(Q_fn, bounds, n, rule).values()))
    return lhs, rhs


def flux_through_sphere(Q_fn: Callable, R: float, n_theta: int = 48, n_phi: int = 96, center=(0.0, 0.0, 0.0)):
    """∯ n·Q dA over the sphere of radius R (outward n = (x − c)/R): Gauss–Legendre in θ, midpoint in the periodic φ.

    Book: §2.12 (divergence theorem on a closed surface). Both rules converge spectrally for smooth Q (exact for
    polynomial fields once n_theta, n_phi are large enough).

    Validation: V5 Q = (2x, y², z²), R = 1 → 8π/3 (Wikipedia "Divergence theorem" example). Label: benchmark.
    """
    c = _F(center)
    th, wth = gauss_legendre_nodes(0.0, np.pi, n_theta)
    ph, wph = midpoint_nodes(0.0, 2.0 * np.pi, n_phi)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    N = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)])  # outward unit normal
    X = c[:, None, None] + R * N
    Q = _eval(Q_fn, X[0], X[1], X[2])
    dA = R ** 2 * np.sin(TH) * np.einsum("i,j->ij", wth, wph)
    return float(np.sum(np.einsum("cij,cij->ij", N, Q) * dA))  # ∯ n·Q dA


def divergence_theorem_sphere(Q_fn: Callable, R: float, n: int = 48, div_fn: Callable | None = None,
                              center=(0.0, 0.0, 0.0), fd_step: float = 1e-3):
    """Divergence theorem on a ball: ∭ ∇·Q dV in spherical coordinates (Gauss–Legendre in r and θ, midpoint in φ)
    vs ∯ n·Q dA.

    Book: §2.12 (vector form of (2.30)). Returns ``(lhs, rhs)``.

    Validation: V5 F = (2x, y², z²), unit sphere: both sides → 8π/3 (rel. 1e-12 at n = 24). Label: benchmark, converged.
    """
    c = _F(center)
    if div_fn is None:
        def div_fn(X, Y, Z):
            d = fd_partials(Q_fn, (X, Y, Z), fd_step)
            return d[0, 0] + d[1, 1] + d[2, 2]
    r, wr = gauss_legendre_nodes(0.0, R, n)
    th, wth = gauss_legendre_nodes(0.0, np.pi, n)
    ph, wph = midpoint_nodes(0.0, 2.0 * np.pi, 2 * n)
    Rr, TH, PH = np.meshgrid(r, th, ph, indexing="ij")
    X = c[0] + Rr * np.sin(TH) * np.cos(PH)
    Y = c[1] + Rr * np.sin(TH) * np.sin(PH)
    Z = c[2] + Rr * np.cos(TH)
    dV = Rr ** 2 * np.sin(TH) * np.einsum("i,j,k->ijk", wr, wth, wph)
    lhs = float(np.sum(_eval(div_fn, X, Y, Z) * dV))
    rhs = flux_through_sphere(Q_fn, R, n, 2 * n, center)
    return lhs, rhs


# ----------------------------------------------------------------------------------------------------------------------
# §2.12 integral definitions (2.31)–(2.33) on a small cube
# ----------------------------------------------------------------------------------------------------------------------
def _cube_face_sum(Q_fn: Callable, x0, h: float, n_face: int, weight: Callable):
    """Σ_faces ∬ weight(n, Q) dA over the cube of side h centred at x0 (midpoint rule with n_face² nodes per face)."""
    x0 = _F(x0)
    s, ws = midpoint_nodes(-0.5 * h, 0.5 * h, n_face)
    T, S = np.meshgrid(s, s, indexing="ij")
    W = np.einsum("j,i->ji", ws, ws)
    total = None
    for d in range(3):
        others = [k for k in range(3) if k != d]
        for sign in (+1, -1):
            n = np.zeros(3)
            n[d] = sign
            coords = [None, None, None]
            coords[d] = np.full_like(S, x0[d] + sign * 0.5 * h)
            coords[others[0]] = x0[others[0]] + S
            coords[others[1]] = x0[others[1]] + T
            Q = _eval(Q_fn, *coords)
            contrib = np.einsum("...ji,ji->...", weight(n, Q), W)
            total = contrib if total is None else total + contrib
    return total


def integral_gradient(Q_fn: Callable, x0, h: float = 0.1, n_face: int = 8):
    """Generalised field derivative 𝒟Q ≈ (1/V) ∬_A n_i Q dA on a cube of side h centred at x0 (V = h³).

    Book: §2.12, Eq. (2.31) — the small-volume limit of (2.30); "(2.31) defines the gradient of a tensor Q of any
    order". Returns shape ``(3,)`` for scalar Q ((∇Q)_i) or ``(3, 3)`` for vector Q (``[i, j]`` = ∂Q_j/∂x_i, derivative
    index first as in n_i Q_j). **Index order:** this is the transpose of ``core.operators.vector_gradient``, whose
    G[i, j] = ∂u_i/∂x_j (the Ch. 3 velocity-gradient convention) — ``integral_gradient(u, …) ≈ vector_gradient(u, …).T``
    at the point. Exact for linear fields; error O(h²) otherwise (D21, D22).

    Validation: V1 linear Q exact; V3 observed order 2.0 in h against the sympy gradient (h = 0.4 … 0.05).
    Label: analytic, converged.
    """
    def weight(n, Q):
        return np.einsum("i,...->i...", n, Q)  # n_i Q  (free index i first)
    return _cube_face_sum(Q_fn, x0, h, n_face, weight) / h ** 3  # Eq. (2.31)


def integral_divergence(Q_fn: Callable, x0, h: float = 0.1, n_face: int = 8):
    """∇·Q ≈ (1/V) ∬_A n·Q dA on a cube of side h centred at x0 — the divergence as outflux per unit volume.

    Book: §2.12, Eq. (2.32); Example 2.5 (Cartesian formula recovered as h → 0).

    Validation: V1 exact for linear fields (Example 2.3: a x → 3a, b × x → 0); V3 order 2 in h to the sympy divergence;
    agrees with ``core.operators.divergence`` to truncation. Label: analytic, converged.
    """
    def weight(n, Q):
        return np.einsum("i,i...->...", n, Q)  # n·Q
    return float(_cube_face_sum(Q_fn, x0, h, n_face, weight) / h ** 3)  # Eq. (2.32)


def integral_curl(Q_fn: Callable, x0, h: float = 0.1, n_face: int = 8):
    """∇×Q ≈ (1/V) ∬_A n × Q dA on a cube of side h centred at x0.

    Book: §2.12, Eq. (2.33).

    Validation: V1 exact for linear fields (b × x → 2b); V3 order 2 in h to the sympy curl. Label: analytic, converged.
    """
    eps = levi_civita()

    def weight(n, Q):
        return np.einsum("ijk,j,k...->i...", eps, n, Q)  # (n × Q)_i = ε_ijk n_j Q_k
    return _cube_face_sum(Q_fn, x0, h, n_face, weight) / h ** 3  # Eq. (2.33)


# ----------------------------------------------------------------------------------------------------------------------
# 2-D versions (plane field, unit depth) for the notebook figures and the gauss_flux_box explainer
# ----------------------------------------------------------------------------------------------------------------------
def flux_through_faces(Q_fn: Callable, bounds, n: int = 64) -> dict[str, float]:
    """Outward flux ∫ n·Q ds through the four sides of a rectangle for a plane field Q(X, Y) → (2, ...) (unit depth).

    Book: §2.12, the divergence theorem with a z-independent field (Example 2.5's bookkeeping in two directions).
    Keys ``"+x", "-x", "+y", "-y"``  [Q unit × m].
    """
    b = _bounds2(bounds)
    out = {}
    for d, name in enumerate("xy"):
        o = 1 - d
        s, ws = midpoint_nodes(*b[o], n)
        for sign, key in ((+1, "+"), (-1, "-")):
            coords = [None, None]
            coords[d] = np.full_like(s, b[d, 1] if sign > 0 else b[d, 0])
            coords[o] = s
            Q = _eval(Q_fn, *coords)
            out[f"{key}{name}"] = float(np.sum(sign * Q[d] * ws))
    return out


def divergence_theorem_rect2d(Q_fn: Callable, bounds, n: int = 64, div_fn: Callable | None = None,
                              fd_step: float = 1e-3):
    """Plane divergence theorem: ∬ ∇·Q dA (lhs) vs ∮ n·Q ds (rhs) on a rectangle, with the four face fluxes.

    Book: §2.12 (the divergence theorem, 2-D reading of Fig. 2.9 / Example 2.5).

    Returns
    -------
    (lhs, rhs, faces) : (float, float, dict)

    Validation: V1 Q = (x, y) on [0,1]²: 2 == 2; Q = (x², 0): 1 == 1; V3 order 2; V4 solenoidal (−y, x): rhs = 0.
    Label: analytic, converged.
    """
    b = _bounds2(bounds)
    if div_fn is None:
        def div_fn(X, Y):
            d = fd_partials(Q_fn, (X, Y), fd_step)
            return d[0, 0] + d[1, 1]
    x, wx = midpoint_nodes(*b[0], n)
    y, wy = midpoint_nodes(*b[1], n)
    X, Y = np.meshgrid(x, y, indexing="xy")
    W = np.einsum("j,i->ji", wy, wx)
    lhs = float(np.sum(_eval(div_fn, X, Y) * W))
    faces = flux_through_faces(Q_fn, b, n)
    return lhs, float(sum(faces.values())), faces


def integral_divergence_2d(Q_fn: Callable, x0, h: float = 0.1, n_side: int = 8) -> float:
    """∇·Q ≈ (1/A) ∮ n·Q ds on a square of side h centred at x0 for a plane field (Eq. (2.32), unit depth).

    Validation: V1 exact for linear fields; V3 order 2 in h. Label: analytic, converged.
    """
    x0 = _F(x0)
    b = np.array([[x0[0] - 0.5 * h, x0[0] + 0.5 * h], [x0[1] - 0.5 * h, x0[1] + 0.5 * h]])
    return float(sum(flux_through_faces(Q_fn, b, n_side).values()) / h ** 2)


class TiledFlux(NamedTuple):
    """Result of :func:`divergence_theorem_tiled`; a plain tuple ``(sum_tiles, outer, interior)`` for parity rows."""

    sum_tiles: float  # Σ over every tile of its four signed face fluxes
    outer: float  # Σ over the tile faces that lie on the outer boundary  (= ∮_A n·Q ds)
    interior: float  # Σ over the tile faces shared by two tiles (expected 0 to round-off)


def divergence_theorem_tiled(Q_fn: Callable, bounds, tiles: int = 4, n: int = 16) -> TiledFlux:
    """Tile a rectangle into ``tiles × tiles`` sub-boxes and add up their boundary fluxes: the interior faces cancel.

    Book: §2.12, Eq. (2.30) — the step of Gauss' theorem the book leaves to the reader (D25 step 8): if (2.30) holds on
    every small box V_k, summing Σ_k ∬_{A_k} n·Q dA over the tiles counts each interior face twice with opposite
    outward normals (n and −n) and the same Q, so those terms cancel and only the outer boundary A survives:

        Σ_k ∮_{A_k} n·Q ds  =  ∮_A n·Q ds  +  Σ_{interior faces} (n·Q + (−n)·Q) ds  =  ∮_A n·Q ds.      # Eq. (2.30)

    Plane version (z-independent field, unit depth) to match :func:`flux_through_faces` and the E4 explainer.

    Parameters
    ----------
    Q_fn : callable ``Q(X, Y) → (2, ...)``  [Q unit]
    bounds : (lo, hi) or ((x0, x1), (y0, y1))  [m]
    tiles : number of sub-boxes per direction (≥ 1; ``tiles = 1`` reproduces :func:`flux_through_faces`)
    n : midpoint nodes per tile side; the outer boundary therefore carries ``tiles × n`` nodes per side, exactly the
        nodes :func:`divergence_theorem_rect2d` uses with ``n_total = tiles × n``

    Returns
    -------
    TiledFlux(sum_tiles, outer, interior) : floats [Q unit × m]
        ``sum_tiles = outer + interior`` identically; ``interior`` is the explicit sum over the shared faces
        (each visited from both sides), not a difference of two totals, so a non-zero value would be a real bug.

    Assumptions: Q continuous on the box (a shared face sees the *same* Q from both tiles — that is what cancels).
    Validation: V1 ``interior`` = 0 to ~1e-15 and ``outer`` = ``divergence_theorem_rect2d(..., tiles*n)[1]`` to 1e-13
    for Q = (x, y), (x², 0) and the smooth test field; V4 solenoidal (−y, x): every number 0. Label: analytic.
    """
    b = _bounds2(bounds)
    tiles = int(tiles)
    if tiles < 1:
        raise ValueError("tiles must be >= 1")
    xs = np.linspace(b[0, 0], b[0, 1], tiles + 1)
    ys = np.linspace(b[1, 0], b[1, 1], tiles + 1)
    outer = 0.0
    interior = 0.0
    for i in range(tiles):
        for j in range(tiles):
            faces = flux_through_faces(Q_fn, ((xs[i], xs[i + 1]), (ys[j], ys[j + 1])), n)  # (2.30) on tile V_k
            on_outer = {"-x": i == 0, "+x": i == tiles - 1, "-y": j == 0, "+y": j == tiles - 1}
            for key, val in faces.items():
                if on_outer[key]:
                    outer += val  # this tile face is part of A
                else:
                    interior += val  # shared face: its twin (opposite n, same Q) is added by the neighbour
    return TiledFlux(float(outer + interior), float(outer), float(interior))


# ----------------------------------------------------------------------------------------------------------------------
# §2.13 Stokes' theorem: loops, surfaces, circulation
# ----------------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Loop:
    """A closed oriented curve C sampled at midpoints: ``points`` (N, d), unit ``tangents`` t (N, d), arc elements
    ``ds`` (N,), the ``normal`` n of the surface it bounds (d,), ``center`` (d,) and enclosed ``area`` [m²]."""

    points: np.ndarray
    tangents: np.ndarray
    ds: np.ndarray
    normal: np.ndarray
    center: np.ndarray
    area: float

    @property
    def dim(self) -> int:
        return self.points.shape[1]

    @property
    def length(self) -> float:
        return float(np.sum(self.ds))


@dataclass(frozen=True)
class Surface:
    """An oriented surface patch sampled at midpoints: ``points`` (M, d), unit ``normals`` (M, d), areas ``dA`` (M,),
    plus the geometry needed to decide whether a point lies on the patch (``kind``: "disc" or "rectangle")."""

    points: np.ndarray
    normals: np.ndarray
    dA: np.ndarray
    normal: np.ndarray
    center: np.ndarray
    kind: str
    size: tuple

    @property
    def dim(self) -> int:
        return self.points.shape[1]

    @property
    def area(self) -> float:
        return float(np.sum(self.dA))

    def contains(self, point) -> bool:
        """True if ``point`` (projected onto the plane) lies inside the patch."""
        pad = lambda v: np.concatenate([_F(v), np.zeros(3 - len(v))]) if len(v) < 3 else _F(v)  # noqa: E731
        e1, e2, _ = _plane_basis(pad(self.normal))
        rel = pad(point) - pad(self.center)
        a, b = float(rel @ e1), float(rel @ e2)
        if self.kind == "disc":
            return bool(np.hypot(a, b) <= self.size[0])
        return bool(abs(a) <= 0.5 * self.size[0] and abs(b) <= 0.5 * self.size[1])


def boundary_tangent(n_c, n) -> np.ndarray:
    """Unit tangent to the boundary curve: t = n_c × n, with n_c the unit normal to C that is tangent to A and points
    **into** the surface A (Fig. 2.10), and n the surface normal on the chosen outside.

    Book: §2.13 ("the unit tangent vector to C, t, points in the counterclockwise direction when looking at the outside
    of A; it is defined as t = n_c × n"); Fig. 2.10 draws n_c from the rim *into* A — with that choice t = n_c × n runs
    counterclockwise about n (an outward-pointing n_c would make t clockwise). (n_c, n, t) is right-handed.

    Parameters
    ----------
    n_c : array_like, shape (3,) — in-surface normal to C pointing into A (need not be unit)
    n : array_like, shape (3,) — surface normal on the outside of A (need not be unit)

    Returns
    -------
    t : ndarray, shape (3,) — unit tangent to C, counterclockwise about n

    Validation: V1 unit length; det[n_c, n, t] = +1; on a circle (n_c = −radial, n = e₃) t agrees with the derivative
    of the counterclockwise parametrisation. Label: analytic.
    """
    t = np.cross(_F(n_c), _F(n))  # t = n_c × n  (Fig. 2.10)
    return t / np.linalg.norm(t)


def _plane_basis(normal):
    """Orthonormal (e1, e2, n) with e1 × e2 = n (right-handed) for a plane with unit normal n."""
    n = _F(normal)
    n = n / np.linalg.norm(n)
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = a - (a @ n) * n
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)  # so that e1 × e2 = n
    return e1, e2, n


def _embed(center, normal):
    """Handle the 2-D convenience: a 2-vector centre means the plane z = 0 with n = e3."""
    c = _F(center)
    if c.shape[0] == 2:
        return np.array([c[0], c[1], 0.0]), np.array([0.0, 0.0, 1.0]), 2
    if normal is None:
        raise ValueError("a 3-D centre needs a normal")
    return c, _F(normal) / np.linalg.norm(_F(normal)), 3


def _strip(arr: np.ndarray, dim: int) -> np.ndarray:
    return arr[..., :dim] if dim == 2 else arr


def planar_loop(center, normal=None, radius: float = 1.0, n: int = 256) -> Loop:
    """Circle of radius R about ``center`` in the plane ⊥ ``normal``, oriented counterclockwise about ``normal``
    (Fig. 2.10: t = n_c × n with n_c = −(x − c)/R the in-plane normal to C pointing *into* the disc A, towards the
    centre). A 2-vector centre → the plane z = 0, n = e₃.

    Book: §2.13, Fig. 2.10 (orientation of C); the disc/circle pair of the C16 worked number (b × x, Γ = 2|b|πR²).

    Validation: V1 ∮ (x − c) × t ds = 2A n (A = πR²); tangents equal the derivative of the points. Label: analytic.
    """
    c3, n3, dim = _embed(center, normal)
    e1, e2, _ = _plane_basis(n3)
    phi, w = midpoint_nodes(0.0, 2.0 * np.pi, n)
    pts = c3[None, :] + radius * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)
    tan = -np.sin(phi)[:, None] * e1 + np.cos(phi)[:, None] * e2  # d x/dφ normalised: counterclockwise about n
    return Loop(_strip(pts, dim), _strip(tan, dim), radius * w, _strip(n3, dim) if dim == 3 else n3,
                _strip(c3, dim), float(np.pi * radius ** 2))


def rectangle_loop(center, normal=None, a: float = 1.0, b: float = 1.0, n: int = 64) -> Loop:
    """Rectangle with sides a (along e₁) and b (along e₂) in the plane ⊥ ``normal``, counterclockwise about ``normal``,
    ``n`` midpoint samples per side (Example 2.6's contour in the x = const plane when ``normal = e_x``).

    Validation: V1 ∮ (x − c) × t ds = 2ab n; total length 2(a + b). Label: analytic.
    """
    c3, n3, dim = _embed(center, normal)
    e1, e2, _ = _plane_basis(n3)
    s, w = midpoint_nodes(-0.5, 0.5, n)
    pts, tan, ds = [], [], []
    # side 1: e1 direction at −b/2 … side 4 closing the loop counterclockwise about n
    sides = [(0.5 * b * -e2, a * e1, e1, a), (0.5 * a * e1, b * e2, e2, b),
             (0.5 * b * e2, -a * e1, -e1, a), (-0.5 * a * e1, -b * e2, -e2, b)]
    for offset, span, t_dir, length in sides:
        pts.append(c3 + offset + s[:, None] * span)
        tan.append(np.tile(t_dir, (n, 1)))
        ds.append(length * w)
    return Loop(_strip(np.vstack(pts), dim), _strip(np.vstack(tan), dim), np.concatenate(ds),
                _strip(n3, dim) if dim == 3 else n3, _strip(c3, dim), float(a * b))


def planar_disc(center, normal=None, radius: float = 1.0, nr: int = 32, ntheta: int = 64) -> Surface:
    """Disc of radius R ⊥ ``normal`` sampled on a polar midpoint grid; normals all equal ``normal`` (Fig. 2.10's A).

    Validation: V1 Σ dA = πR² (1e-12). Label: analytic.
    """
    c3, n3, dim = _embed(center, normal)
    e1, e2, _ = _plane_basis(n3)
    r, wr = midpoint_nodes(0.0, radius, nr)
    th, wth = midpoint_nodes(0.0, 2.0 * np.pi, ntheta)
    Rr, TH = np.meshgrid(r, th, indexing="ij")
    pts = c3 + (Rr * np.cos(TH))[..., None] * e1 + (Rr * np.sin(TH))[..., None] * e2
    dA = (Rr * np.einsum("i,j->ij", wr, wth)).ravel()
    P = pts.reshape(-1, 3)
    return Surface(_strip(P, dim), _strip(np.tile(n3, (P.shape[0], 1)), dim), dA, _strip(n3, dim) if dim == 3 else n3,
                   _strip(c3, dim), "disc", (float(radius),))


def planar_rectangle(center, normal=None, a: float = 1.0, b: float = 1.0, n: int = 32) -> Surface:
    """Rectangle a × b ⊥ ``normal`` sampled with n × n midpoints (the surface bounded by :func:`rectangle_loop`)."""
    c3, n3, dim = _embed(center, normal)
    e1, e2, _ = _plane_basis(n3)
    s, ws = midpoint_nodes(-0.5 * a, 0.5 * a, n)
    t, wt = midpoint_nodes(-0.5 * b, 0.5 * b, n)
    S, T = np.meshgrid(s, t, indexing="ij")
    pts = c3 + S[..., None] * e1 + T[..., None] * e2
    dA = np.einsum("i,j->ij", ws, wt).ravel()
    P = pts.reshape(-1, 3)
    return Surface(_strip(P, dim), _strip(np.tile(n3, (P.shape[0], 1)), dim), dA, _strip(n3, dim) if dim == 3 else n3,
                   _strip(c3, dim), "rectangle", (float(a), float(b)))


def _field_at_points(fn: Callable, pts: np.ndarray) -> np.ndarray:
    """Evaluate ``fn`` at points (N, d) → (N, ncomp), padding a 2-component plane field to 3 components if needed."""
    U = _eval(fn, *pts.T)  # (ncomp, N)
    U = U.T
    if U.shape[1] < pts.shape[1]:
        U = np.hstack([U, np.zeros((U.shape[0], pts.shape[1] - U.shape[1]))])
    return U


def circulation(u_fn: Callable, loop: Loop) -> float:
    """Circulation ∮_C u·t ds around the loop (midpoint rule on the sampled curve).

    Book: §2.13, the right side of Eq. (2.34) ("In fluid mechanics, the right side of (2.34) is called the circulation
    of u about C").

    Validation: V1 solid-body rotation b × x around a circle of radius R ⊥ b: 2|b|πR² (rel 1e-12); ∇φ: 0
    (Exercise 2.20); V3 order 2 in the loop resolution for a smooth field (circle, midpoint is spectral → exact for
    polynomial fields). Label: analytic, converged.
    """
    U = _field_at_points(u_fn, loop.points)
    return float(np.sum(np.einsum("ni,ni->n", U, loop.tangents) * loop.ds))  # ∮ u·t ds  (2.34 right side)


def curl_flux(u_fn: Callable, surface: Surface, curl_fn: Callable | None = None, fd_step: float = 1e-3) -> float:
    """Flux of the curl through the surface: ∬_A (∇×u)·n dA.

    Book: §2.13, the left side of Eq. (2.34). ``curl_fn(X, Y[, Z])`` may supply the exact curl (3 components in 3-D,
    the scalar (∇×u)₃ in 2-D); otherwise the curl is formed from fourth-order differences of ``u_fn``.

    Validation: V1 b × x through a disc ⊥ b: 2|b|πR²; equals ``circulation`` for smooth fields (V3 order 2). Label:
    analytic, converged.
    """
    P = surface.points
    if curl_fn is not None:
        C = _eval(curl_fn, *P.T)
        if surface.dim == 2:
            return float(np.sum(np.asarray(C) * surface.dA))  # n = e₃ → (∇×u)₃ dA
        return float(np.sum(np.einsum("in,ni->n", C, surface.normals) * surface.dA))
    d = fd_partials(u_fn, tuple(P.T), fd_step)  # d[j, k, n] = ∂u_k/∂x_j
    if surface.dim == 2:
        c3 = d[0, 1] - d[1, 0]  # (∇×u)₃ = ∂u₂/∂x₁ − ∂u₁/∂x₂
        return float(np.sum(c3 * surface.dA))
    C = np.einsum("ijk,jkn->in", levi_civita(), d)  # (∇×u)_i = ε_ijk ∂u_k/∂x_j  (2.24)
    return float(np.sum(np.einsum("in,ni->n", C, surface.normals) * surface.dA))


@dataclass(frozen=True)
class StokesCheck:
    """Both sides of Stokes' theorem and whether its hypothesis (u differentiable on A) holds. Iterates as (lhs, rhs)."""

    lhs: float  # ∬ (∇×u)·n dA
    rhs: float  # ∮ u·t ds
    hypothesis_ok: bool
    note: str

    def __iter__(self):
        yield self.lhs
        yield self.rhs

    @property
    def mismatch(self) -> float:
        return abs(self.lhs - self.rhs)


def stokes_theorem_check(u_fn: Callable, loop: Loop, surface: Surface, curl_fn: Callable | None = None,
                         fd_step: float = 1e-3) -> StokesCheck:
    """Evaluate both sides of Stokes' theorem ∬_A (∇×u)·n dA = ∮_C u·t ds and report whether the hypothesis holds.

    Book: §2.13, Eq. (2.34). If ``u_fn`` carries a ``singular_at`` attribute (e.g. the irrotational vortex u_θ = K/r,
    ``core.fields.VectorField``) and that point lies on the surface, ``hypothesis_ok`` is False and the note says the
    two sides need not agree (the circulation is then 2πK regardless of the loop size).

    Returns
    -------
    StokesCheck (``lhs, rhs = stokes_theorem_check(...)`` also works).

    Validation: V1 b × x on a disc: lhs = rhs = 2|b|πR²; reversing the normal flips both signs; V3 rectangle loop vs
    curl flux for a smooth field, order 2; irrotational vortex with the core inside → hypothesis_ok False.
    Label: analytic, converged.
    """
    lhs = curl_flux(u_fn, surface, curl_fn, fd_step)
    rhs = circulation(u_fn, loop)
    sing = getattr(u_fn, "singular_at", None)
    ok, note = True, "u is smooth on A: the two sides agree to quadrature error"
    if sing is not None and surface.contains(sing):
        ok = False
        note = ("the field is singular at a point inside A: Stokes' theorem does not apply — the curl flux misses the "
                "core and the circulation is fixed by it")
    return StokesCheck(lhs, rhs, ok, note)


def integral_curl_component(u_fn: Callable, x0, n=None, h: float = 0.1, n_side: int = 8) -> float:
    """n·(∇×u) ≈ (1/A) ∮_C u·t ds around a square of side h centred at x0 in the plane ⊥ n (A = h²).

    Book: §2.13, Eq. (2.35) and Example 2.6 (rectangles in the coordinate planes, midpoint values on each side; the
    second bracket of the book's x-component contains u_y, not the printed u_z). A 2-vector x0 means n = e₃.

    Validation: V1 exact for linear fields (b × x → 2 n·b); V3 order 2 in h to n·(∇×u) from sympy / (2.25).
    Label: analytic, converged.
    """
    loop = rectangle_loop(x0, n, h, h, n_side)
    return circulation(u_fn, loop) / (h * h)  # Eq. (2.35)


# ----------------------------------------------------------------------------------------------------------------------
# Ch. 5: Gauss' theorem in curl form (5.15)
# ----------------------------------------------------------------------------------------------------------------------
class CurlCheck(NamedTuple):
    """Both sides of ∭_V ∇×F dV = ∯_A n×F dA (vectors, shape (3,)) and their difference ``diff`` = volume − surface."""

    volume: np.ndarray
    surface: np.ndarray
    diff: np.ndarray


def curl_theorem_box(F_fn: Callable, bounds, n: int = 24, rule: str = "gauss", curl_fn: Callable | None = None,
                     fd_step: float = 1e-3) -> CurlCheck:
    """Gauss' theorem in curl form on a box: ∭_V ∇×F dV = ∯_A n × F dA.

    Book: §5.5, Eq. (5.15): ∫_V′ ∇′×(ω/|x − x′|) d³x′ = ∫_V′ ε_kij ∂/∂x′_i(ω_j/|x − x′|) d³x′ = ∫_A′ ε_kij
    (ω_j/|x − x′|) n_i d²x′ = ∫_A′ n × ω/|x − x′| d²x′ — Gauss (2.30) applied to each Cartesian component with the
    Levi-Civita symbol; here for a general smooth vector field F. The companion of ``divergence_theorem_box``.

    Parameters
    ----------
    F_fn : callable ``F(X, Y, Z) → (3, …)`` (ch02 coordinate convention)
    bounds : ((x0, x1), (y0, y1), (z0, z1)) [m] or (lo, hi)
    n : quadrature nodes per direction;  rule : "gauss" (default; exact for polynomials of degree ≤ 2n − 1),
        "simpson" or "midpoint"
    curl_fn : exact ∇ × F(X, Y, Z) → (3, …), or None (fourth-order differences of ``F_fn`` with step ``fd_step``)

    Returns
    -------
    CurlCheck(volume, surface, diff) — vectors [F unit × m²].

    Validation: V1 polynomial F (Gauss nodes → equal to round-off); F = b × x on the unit cube → 2b·V; V3 order of the
    midpoint rule on a smooth non-polynomial F. Label: analytic, converged.
    """
    b = _bounds3(bounds)
    if curl_fn is None:
        def curl_fn(X, Y, Z):
            d = fd_partials(F_fn, (X, Y, Z), fd_step)  # d[i, j, ...] = ∂F_j/∂x_i
            return np.einsum("kij,ij...->k...", levi_civita(), d)  # (∇×F)_k = ε_kij ∂F_j/∂x_i
    vol = _F(volume_integral_box(curl_fn, b, n, rule))
    eps = levi_civita()
    surf = np.zeros(3)
    for d in range(3):
        others = [k for k in range(3) if k != d]
        (s, ws), (t, wt) = (_nodes(*b[k], n, rule) for k in others)
        T, S = np.meshgrid(t, s, indexing="ij")
        W = np.einsum("j,i->ji", wt, ws)
        for sign in (+1.0, -1.0):
            coords = [None, None, None]
            coords[d] = np.full_like(S, b[d, 1] if sign > 0 else b[d, 0])
            coords[others[0]], coords[others[1]] = S, T
            Fv = _eval(F_fn, *coords)  # (3, …)
            nvec = np.zeros(3)
            nvec[d] = sign
            nxF = np.einsum("kij,i,j...->k...", eps, nvec, Fv)  # (n × F)_k = ε_kij n_i F_j   (5.15)
            surf += (nxF * W).reshape(3, -1).sum(axis=1)
    return CurlCheck(vol, surf, vol - surf)
