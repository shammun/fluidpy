"""Singularity-distribution ("panel") methods: the book's axial singularity method for bodies of revolution (§6.8,
Fig. 6.29) and — our labelled extension for the explainer seed — 2-D constant-strength source panels on a body surface.

Book: Kundu, Cohen & Dowling 5e, §6.8 (6.93)–(6.95) and the text after Fig. 6.29 (page p263):
ψ_m = −Σ_n (k_n/4π)(r^m_{n−1} − r^m_n) + ½U R_m² = 0 at N body points → an N × N linear system.
Source panels: the constant-strength method of J. L. Hess & A. M. O. Smith, "Calculation of potential flow about
arbitrary bodies", Prog. Aerosp. Sci. 8, 1–138 (1967), doi:10.1016/0376-0421(67)90003-6 (method citation only); the
panel influence integrals below are derived in closed form in local panel coordinates (not copied).
Reused by Ch. 14 (vortex panels, thin-airfoil numerics) and Ch. 10 (boundary-element idea).
"""
from __future__ import annotations

import warnings

import numpy as np

from ._util import as_scalar_if_0d
from .potential import LineSource3D

__all__ = ["axial_influence_matrix", "axial_singularity_solve", "axial_singularity_psi", "axial_singularity_velocity",
           "panel_geometry", "panel_induced_velocity", "source_panels", "panel_velocity"]

_S = as_scalar_if_0d
_FOUR_PI = 4.0 * np.pi


def _F(a):
    return np.asarray(a, dtype=float)


# ======================================================================================================================
# axial singularity method (§6.8)
# ======================================================================================================================
def axial_influence_matrix(z_body, R_body, xi_nodes) -> np.ndarray:
    """Influence matrix of unit-density axial line sources on the Stokes stream function at body points:
    A[m, n] = −(1/4π)(r^m_{n−1} − r^m_n), r^m_n = √(R_m² + (z_m − ξ_n)²), so ψ_m = Σ_n A[m, n] k_n + ½UR_m².

    z_body, R_body (M,) [m]; xi_nodes (N + 1,) segment ends on the axis [m]. Returns (M, N) [m].
    Book: §6.8, the ψ_mn formula after Fig. 6.29 (from (6.94), sign flipped for a source). Label: analytic."""
    z = _F(z_body)[:, None]
    R = _F(R_body)[:, None]
    xi = _F(xi_nodes)[None, :]
    r = np.sqrt(R ** 2 + (z - xi) ** 2)  # r^m_n for every node
    return -(r[:, :-1] - r[:, 1:]) / _FOUR_PI  # ψ_mn per unit k_n


def _nodes_from_body(zb):
    """Segment ends that put each collocation point at its segment's midpoint: interior ends halfway between
    consecutive points, the two outer ends mirrored (for body points at the midpoints of uniform segments from nose to
    tail, this returns exactly nose … tail)."""
    zb = np.sort(_F(zb).ravel())
    if zb.size < 2:
        raise ValueError("need at least two body points")
    mid = 0.5 * (zb[:-1] + zb[1:])
    return np.concatenate([[zb[0] - (mid[0] - zb[0])], mid, [zb[-1] + (zb[-1] - mid[-1])]])


def axial_singularity_solve(z_body, R_body, U: float, N=None, z_nodes=None) -> dict:
    """The axial singularity (inverse) method (§6.8): find segment strengths k_n [m²/s] of N uniform line sources on
    the axis so that ψ = 0 at the body points: −Σ_n (k_n/4π)(r^m_{n−1} − r^m_n) + ½UR_m² = 0.

    Parameters: z_body, R_body (M,) body points [m]; U free stream along +z [m/s]; ``N`` number of segments (default
    M) — or, for backward compatibility, an array of segment ends; ``z_nodes`` (N + 1,) segment ends [m]. Default
    segments: from the body's nose to its tail with each body point above its segment's midpoint (Fig. 6.29); if N ≠ M
    the same span is divided uniformly. M = N → square solve; M > N → least squares.
    Returns dict(k (N,), z_nodes (= xi_nodes), xi_mid, residual (ψ at the body points), max_residual, cond (the
    **1-norm** condition number ``np.linalg.cond(A, 1)``, as the explainer computes it), net_strength Σ k_n Δξ_n (→ 0
    for a closed body), A, rhs, U, odd_symmetric). Book: §6.8 (Fig. 6.29, "set ψ_m = 0 … N linear algebraic
    equations"), solved by ``np.linalg.solve`` (or ``lstsq``) — the book allows iteration or matrix inversion.
    ⚠️ Odd N on a fore–aft symmetric body (symmetric body points and segments: Rankine oval, ellipsoid, sphere): the
    reflection z → −z gives A = −PAP (P the flip), so A maps symmetric strength patterns onto antisymmetric ψ patterns,
    and for odd N the symmetric subspace (⌈N/2⌉) is one dimension larger than the antisymmetric one (⌊N/2⌋) — A is
    exactly singular (cond ~ 1e16–1e17). This case is detected (``odd_symmetric`` = True), a RuntimeWarning is issued
    and the minimum-norm least-squares solution (the antisymmetric k) is returned; ψ_m = 0 then holds only in the
    least-squares sense. Use even N. Even N and asymmetric bodies (the airship) are solved by ``np.linalg.solve``.
    Validation: V1 — tests/test_ch06.py: test_axial_method_V1_rankine_oval_moments_and_conditioning,
    test_axial_method_V1_solver_and_field_consistency. Label: analytic.
    """
    zb, Rb = _F(z_body).ravel(), _F(R_body).ravel()
    if N is not None and np.ndim(N) > 0:
        z_nodes, N = N, None
    if z_nodes is None:
        base = _nodes_from_body(zb)
        if N is None or int(N) == zb.size:
            xi = base
        else:
            xi = np.linspace(base[0], base[-1], int(N) + 1)
    else:
        xi = _F(z_nodes).ravel()
    A = axial_influence_matrix(zb, Rb, xi)
    rhs = -0.5 * float(U) * Rb ** 2
    odd_sym = False
    if A.shape[0] == A.shape[1]:
        n = A.shape[0]
        if n % 2 == 1:
            P = np.eye(n)[::-1]
            odd_sym = bool(np.allclose(A, -P @ A @ P, rtol=0.0, atol=1e-12 * np.abs(A).max()))
        if odd_sym:
            # A = −PAP (fore–aft symmetric body and segments): A maps symmetric k onto antisymmetric ψ, whose
            # subspace has one dimension fewer for odd N ⇒ A is exactly singular. The minimum-norm least-squares
            # solution is the antisymmetric (source-fore, sink-aft) k; ψ_m = 0 then holds only in the least-squares
            # sense (⌊N/2⌋ unknown pairs for ⌈N/2⌉ symmetric equations).
            warnings.warn(f"axial_singularity_solve: N = {n} is odd and the target is fore–aft symmetric — the "
                          "square system is exactly singular (one symmetric mode too many); using the minimum-norm "
                          "least-squares solution. Use an even N.", RuntimeWarning, stacklevel=2)
            k = np.linalg.lstsq(A, rhs, rcond=None)[0]
        else:
            k = np.linalg.solve(A, rhs)
        cond = float(np.linalg.cond(A, 1))
    else:
        k = np.linalg.lstsq(A, rhs, rcond=None)[0]
        cond = float(np.linalg.cond(A))
    resid = A @ k - rhs  # = ψ_m
    return {"k": k, "z_nodes": xi, "xi_nodes": xi, "xi_mid": 0.5 * (xi[:-1] + xi[1:]), "residual": resid,
            "max_residual": float(np.max(np.abs(resid))), "cond": cond,
            "net_strength": float(np.sum(k * np.diff(xi))), "A": A, "rhs": rhs, "U": float(U),
            "odd_symmetric": odd_sym}


def _sol_args(args):
    """(sol, R, z) or (R, z, k, xi_nodes[, U]) → (R, z, k, xi, U)."""
    if isinstance(args[0], dict):
        sol, R, z = args[:3]
        return R, z, sol["k"], sol["z_nodes"], sol.get("U", 0.0)
    R, z, k, xi = args[:4]
    U = args[4] if len(args) > 4 else 0.0
    return R, z, k, xi, U


def axial_singularity_psi(*args, U: float | None = None):
    """Stokes stream function [m³/s] of the axial segments (strengths k_n on [ξ_{n−1}, ξ_n]) plus a stream U along z:
    ψ = −Σ (k_n/4π)(r_{n−1} − r_n) + ½UR². Call as ``(sol, R, z)`` with a :func:`axial_singularity_solve` result (its U)
    or ``(R, z, k, xi_nodes, U=0)``. Book: §6.8. Label: analytic."""
    R, z, k, xi, U0 = _sol_args(args)
    U = U0 if U is None else U
    R_, z_ = np.broadcast_arrays(_F(R), _F(z))
    xi = _F(xi).ravel()
    out = 0.5 * float(U) * R_ ** 2
    for kn, a, b in zip(_F(k).ravel(), xi[:-1], xi[1:]):
        out = out + LineSource3D(kn, a, b).psi(R_, z_)
    return _S(out)


def axial_singularity_velocity(*args, U: float | None = None):
    """(u_R, u_z) [m/s] of the axial segments plus the stream U e_z (sum of :class:`LineSource3D` velocities).
    Call as ``(sol, R, z)`` or ``(R, z, k, xi_nodes, U=0)``. Book: §6.8. Label: analytic."""
    R, z, k, xi, U0 = _sol_args(args)
    U = U0 if U is None else U
    R_, z_ = np.broadcast_arrays(_F(R), _F(z))
    xi = _F(xi).ravel()
    uR = np.zeros(R_.shape)
    uz = np.full(R_.shape, float(U))
    for kn, a, b in zip(_F(k).ravel(), xi[:-1], xi[1:]):
        a_, b_ = LineSource3D(kn, a, b).velocity_cyl(R_, z_)
        uR = uR + a_
        uz = uz + b_
    return _S(uR), _S(uz)


# ======================================================================================================================
# 2-D constant-strength source panels (our extension; Hess & Smith 1967)
# ======================================================================================================================
def panel_geometry(xb, yb) -> dict:
    """Panels of a closed polygon given by vertices (xb, yb) (last vertex ≠ first; either orientation).

    Returns dict(x0, y0, x1, y1 (panel ends), xm, ym (midpoints = control points), S (lengths), tx, ty (unit
    tangents in the counterclockwise sense), nx, ny (outward unit normals), ccw (bool: input orientation)).
    Label: analytic."""
    x, y = _F(xb).ravel(), _F(yb).ravel()
    area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
    ccw = bool(area > 0)
    if not ccw:
        x, y = x[::-1], y[::-1]
    x0, y0, x1, y1 = x, y, np.roll(x, -1), np.roll(y, -1)
    S = np.hypot(x1 - x0, y1 - y0)
    tx, ty = (x1 - x0) / S, (y1 - y0) / S
    return {"x0": x0, "y0": y0, "x1": x1, "y1": y1, "xm": 0.5 * (x0 + x1), "ym": 0.5 * (y0 + y1), "S": S,
            "tx": tx, "ty": ty, "nx": ty, "ny": -tx, "ccw": ccw}


def panel_induced_velocity(x, y, geo: dict, self_index=None) -> tuple[np.ndarray, np.ndarray]:
    """Velocity at points (x, y) (P,) induced by unit-density sources on each panel: returns (Ux, Uy) of shape (P, N).

    In panel j's frame (x′ along the panel from its start, y′ along the outward normal), a unit density on
    0 ≤ s ≤ S gives u′ = (1/4π) ln[(x′² + y′²)/((x′ − S)² + y′²)], v′ = (θ₂ − θ₁)/2π with θ₁ = atan2(y′, x′),
    θ₂ = atan2(y′, x′ − S) (our derivation of ∫(1/2π)(x − s)/|x − s|² ds). On the panel itself (``self_index[p] = j``)
    the outward-side limit v′ = ½, u′ = (1/4π) ln[x′²/(S − x′)²] is used. Label: analytic."""
    X, Y = _F(x).ravel()[:, None], _F(y).ravel()[:, None]
    dx, dy = X - geo["x0"][None, :], Y - geo["y0"][None, :]
    tx, ty, nx, ny, S = (geo[k][None, :] for k in ("tx", "ty", "nx", "ny", "S"))
    xp = dx * tx + dy * ty
    yp = dx * nx + dy * ny
    with np.errstate(divide="ignore", invalid="ignore"):
        up = np.log((xp ** 2 + yp ** 2) / ((xp - S) ** 2 + yp ** 2)) / _FOUR_PI
        vp = (np.arctan2(yp, xp - S) - np.arctan2(yp, xp)) / (2.0 * np.pi)
    if self_index is not None:
        p = np.arange(X.shape[0])
        jj = np.asarray(self_index)
        xs = xp[p, jj]
        up[p, jj] = np.log(xs ** 2 / (S[0, jj] - xs) ** 2) / _FOUR_PI
        vp[p, jj] = 0.5
    Ux = up * tx + vp * nx
    Uy = up * ty + vp * ny
    return Ux, Uy


def source_panels(xb, yb, U: float = 1.0, alpha: float = 0.0) -> dict:
    """2-D constant-strength source panels (our extension; method of Hess & Smith 1967): strengths λ_j [m/s] on the N
    panels of a closed polygon such that u·n = 0 at the panel midpoints in a stream U at angle α:
    Σ_j λ_j (U_j·n_i) + U∞·n_i = 0 (self term ½λ_i).

    Parameters: body vertices xb, yb [m] (either orientation, last ≠ first); U [m/s]; α [rad].
    Returns dict(lam, xm, ym (control points), vt (tangential velocity, counterclockwise sense), cp = 1 − (v_t/U)²,
    theta (polar angle of the control points), net_source Σλ_jS_j (→ 0 for a closed body), normals (complex outward
    unit normals), lengths, geo, U, alpha, cond).
    Book: extends §6.7–6.8's singularity-distribution idea to the body surface (analysis row #134).
    Validation: circle: C_p at the control points equals 1 − 4 sin²θ to round-off for every N (a symmetry of the regular
    polygon, not a convergence); ellipse: observed order ≈ 2 against the Zhukhovsky-mapped exact surface speed; Σλ_jS_j
    = 0. V1, V3 — tests/test_ch06.py: test_source_panels_V1_circle_is_exact_and_closed,
    test_source_panels_V3_ellipse_second_order_and_off_body_first_order. Label: analytic, converged.
    """
    geo = panel_geometry(xb, yb)
    N = geo["S"].size
    Ux, Uy = panel_induced_velocity(geo["xm"], geo["ym"], geo, self_index=np.arange(N))
    An = Ux * geo["nx"][:, None] + Uy * geo["ny"][:, None]
    At = Ux * geo["tx"][:, None] + Uy * geo["ty"][:, None]
    Uinf = float(U) * np.array([np.cos(alpha), np.sin(alpha)])
    rhs = -(Uinf[0] * geo["nx"] + Uinf[1] * geo["ny"])
    lam = np.linalg.solve(An, rhs)
    vt = At @ lam + Uinf[0] * geo["tx"] + Uinf[1] * geo["ty"]
    return {"lam": lam, "xm": geo["xm"], "ym": geo["ym"], "vt": vt, "cp": 1.0 - (vt / float(U)) ** 2,
            "theta": np.arctan2(geo["ym"], geo["xm"]), "net_source": float(np.sum(lam * geo["S"])), "geo": geo,
            "normals": geo["nx"] + 1j * geo["ny"], "lengths": geo["S"],
            "U": float(U), "alpha": float(alpha), "cond": float(np.linalg.cond(An))}


def panel_velocity(*args):
    """Velocity (u, v) [m/s] anywhere from a :func:`source_panels` solution (free stream + panel sources): call as
    ``(sol, x, y)`` or ``(x, y, sol)``. Points on a panel use the open-panel formula (evaluate off the body).
    Accuracy: off the body the velocity converges only at about first order in N (the constant-strength panels'
    end effects), unlike the second-order C_p at the control points on a smooth body. Validation: V3 —
    tests/test_ch06.py: test_source_panels_V3_ellipse_second_order_and_off_body_first_order. Label: converged."""
    if isinstance(args[0], dict):
        result, x, y = args[:3]
    else:
        x, y, result = args[:3]
    geo = result["geo"]
    Ux, Uy = panel_induced_velocity(x, y, geo)
    sh = np.shape(_F(x))
    u = (Ux @ result["lam"]).reshape(sh) + result["U"] * np.cos(result["alpha"])
    v = (Uy @ result["lam"]).reshape(sh) + result["U"] * np.sin(result["alpha"])
    return _S(u), _S(v)
