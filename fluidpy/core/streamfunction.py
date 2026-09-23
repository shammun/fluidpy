"""Stream functions: steady mass conservation ∇·(ρu) = 0 built in, in 3-D (two stream functions), plane and
axisymmetric flow.

Book: Kundu, Cohen & Dowling 5e, Ch. 4 §4.3, Eqs. (4.11)–(4.12), Fig. 4.1 and the plane/axisymmetric forms that follow
(rendered pages chapters/pages/ch04/p126–p127).

Conventions (``knowledge/notation.md``)
---------------------------------------
* 3-D: ρu = ∇χ × ∇ψ (from ρu = ∇×Ψ (4.12) with Ψ = χ∇ψ); streamlines are the intersections of χ = const and
  ψ = const; the mass flux through a patch bounded by χ = a, b and ψ = c, d is (b − a)(d − c).
* Plane flow: χ = −z ⇒ ρu = ∂ψ/∂y, ρv = −∂ψ/∂x (the "usual convention"; some texts, e.g. several GFD books, use the
  opposite sign). With ρ = 1 (constant density) ψ is the volume flux per unit depth [m²/s].
* Axisymmetric flow in cylindrical (R, φ, z): χ = −φ ⇒ ρu_R = −(1/R)∂ψ/∂z, ρu_z = (1/R)∂ψ/∂R; ψ [m³/s].
* ``psi`` arguments are a callable ψ(x, y) (or ψ(R, z)) **or a preset name** of :data:`STREAMFUNCTION_PRESETS` with
  its keywords (so a parity row can write ``ch04.velocity_from_streamfunction_2d("cylinder", 2.0, 0.5, U=1.0)``).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import sympy as sp

from ._util import as_scalar_if_0d
from .integral_theorems import gauss_legendre_nodes

__all__ = ["velocity_from_streamfunction_2d", "velocity_from_streamfunction_2d_sym",
           "velocity_from_streamfunction_axisym", "velocity_from_streamfunction_axisym_sym",
           "flux_between_streamlines", "flux_along_path", "mass_flux_from_vector_potential",
           "mass_flux_from_stream_functions", "stream_function_pair", "stream_surface_check",
           "stream_tube_mass_flux", "STREAMFUNCTION_PRESETS", "streamfunction_preset", "velocity_preset"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


# ======================================================================================================================
# presets (E2 and the C03 figures) — our examples
# ======================================================================================================================
#: E2's flows: defaults (SI), singular points and stagnation points (for the default parameters).
STREAMFUNCTION_PRESETS: dict[str, dict] = {
    "uniform": dict(defaults=dict(U=1.0, alpha=0.0), singular=[], stagnation=[]),
    "stagnation": dict(defaults=dict(k=1.0), singular=[], stagnation=[(0.0, 0.0)]),
    "source_stream": dict(defaults=dict(U=1.0, m=1.0), singular=[(0.0, 0.0)],
                          stagnation=["(−m/(2πU), 0)"]),
    "cylinder": dict(defaults=dict(U=1.0, a=1.0), singular=[(0.0, 0.0)], stagnation=["(±a, 0)"]),
    "vortex": dict(defaults=dict(Gamma=1.0), singular=[(0.0, 0.0)], stagnation=[]),
    "shear": dict(defaults=dict(gamma_dot=1.0), singular=[], stagnation=["the line y = 0"]),
    "axisym_uniform": dict(defaults=dict(U=1.0), singular=[], stagnation=[]),
}


def _pp(name: str, p: dict) -> dict:
    if name not in STREAMFUNCTION_PRESETS:
        raise ValueError(f"unknown preset {name!r}; choose from {tuple(STREAMFUNCTION_PRESETS)}")
    q = dict(STREAMFUNCTION_PRESETS[name]["defaults"])
    q.update({k: v for k, v in p.items() if k in q})
    return q


def streamfunction_preset(name: str, x, y, **p):
    """Stream function ψ(x, y) of E2's preset flows (§4.3 convention ρu = ∂ψ/∂y, ρv = −∂ψ/∂x, ρ = 1) [m²/s].

    * "uniform" (U, alpha): ψ = U(y cos α − x sin α) (speed U at angle α);  * "stagnation" (k): ψ = kxy;
    * "source_stream" (U, m): ψ = Uy + (m/2π)θ (Rankine half-body, source strength m [m²/s]);
    * "cylinder" (U, a): ψ = Uy(1 − a²/r²) (= ``ch03.cylinder_streamfunction``);
    * "vortex" (Gamma): ψ = −(Γ/2π) ln r (counter-clockwise for Γ > 0);  * "shear" (gamma_dot): ψ = ½γ̇y² (u = γ̇y);
    * "axisym_uniform" (U): Stokes stream function ψ = ½UR² with (x, y) read as (z, R).

    Book: §4.3 (the stream-function definitions); the flows are ours (Ch. 6 derives them). Scalar-callable.
    Validation: V1 velocities from :func:`velocity_from_streamfunction_2d` equal :func:`velocity_preset` at 200 points.
    Label: analytic.
    """
    q = _pp(name, p)
    x_, y_ = _F(x), _F(y)
    with np.errstate(divide="ignore", invalid="ignore"):  # singular points (listed in STREAMFUNCTION_PRESETS) → inf/nan
        return _S(_psi_value(name, q, x_, y_))


def _psi_value(name, q, x_, y_):
    if name == "uniform":
        out = q["U"] * (y_ * np.cos(q["alpha"]) - x_ * np.sin(q["alpha"]))
    elif name == "stagnation":
        out = q["k"] * x_ * y_
    elif name == "source_stream":
        out = q["U"] * y_ + q["m"] / (2 * np.pi) * np.arctan2(y_, x_)
    elif name == "cylinder":
        out = q["U"] * y_ * (1.0 - q["a"] ** 2 / (x_ ** 2 + y_ ** 2))
    elif name == "vortex":
        out = -q["Gamma"] / (2 * np.pi) * np.log(np.hypot(x_, y_))
    elif name == "shear":
        out = 0.5 * q["gamma_dot"] * y_ ** 2
    else:  # axisym_uniform, (x, y) ≡ (z, R)
        out = 0.5 * q["U"] * y_ ** 2
    return out


def velocity_preset(name: str, x, y, **p):
    """Closed-form velocity (u, v) [m/s] of :func:`streamfunction_preset` (the independent route for parity rows); for
    "axisym_uniform" returns (u_z, u_R) = (U, 0) in the (z, R) plane. Book: §4.3. Label: analytic."""
    q = _pp(name, p)
    x_, y_ = _F(x), _F(y)
    with np.errstate(divide="ignore", invalid="ignore"):
        u, v = _vel_value(name, q, x_, y_)
    return _S(u), _S(v)


def _vel_value(name, q, x_, y_):
    r2 = x_ ** 2 + y_ ** 2
    if name == "uniform":
        u, v = q["U"] * np.cos(q["alpha"]) + 0 * x_, q["U"] * np.sin(q["alpha"]) + 0 * x_
    elif name == "stagnation":
        u, v = q["k"] * x_, -q["k"] * y_
    elif name == "source_stream":
        u, v = q["U"] + q["m"] / (2 * np.pi) * x_ / r2, q["m"] / (2 * np.pi) * y_ / r2
    elif name == "cylinder":
        a, U = q["a"], q["U"]
        u = U * (1.0 - a ** 2 * (x_ ** 2 - y_ ** 2) / r2 ** 2)
        v = -2.0 * U * a ** 2 * x_ * y_ / r2 ** 2
    elif name == "vortex":
        u, v = -q["Gamma"] / (2 * np.pi) * y_ / r2, q["Gamma"] / (2 * np.pi) * x_ / r2
    elif name == "shear":
        u, v = q["gamma_dot"] * y_, 0 * x_
    else:
        u, v = q["U"] + 0 * x_, 0 * x_
    return u, v


def _psi_fn(psi, p: dict) -> Callable:
    if isinstance(psi, str):
        name = psi
        return lambda X, Y: streamfunction_preset(name, X, Y, **p)
    return psi


# ======================================================================================================================
# plane and axisymmetric flow
# ======================================================================================================================
def velocity_from_streamfunction_2d(psi, x, y, rho: float = 1.0, h: float = 1e-5, **p):
    """Velocity of a plane flow from its stream function: ρu = ∂ψ/∂y, ρv = −∂ψ/∂x (§4.3, χ = −z).

    Book: §4.3, text after Fig. 4.1: ρu = −e_z × ∇ψ, i.e. ρu = ∂ψ/∂y and ρv = −∂ψ/∂x (Exercise 4.7).

    Parameters
    ----------
    psi : callable ψ(x, y) [kg/(m s)] (or [m²/s] with rho = 1), or a preset name (keywords **p)
    x, y : position(s) [m] (floats or arrays);  rho : density [kg/m³];  h : central-difference step [m]

    Returns
    -------
    (u, v) [m/s] — floats for scalar input (explainer parity).

    Validation: V1 the cylinder stream function (ch03 (3.2)) gives ``ch03.cylinder_flow``; ψ = kxy gives (kx, −ky);
    V2 :func:`velocity_from_streamfunction_2d_sym` and ∇·u ≡ 0; V3 order 2. Label: analytic, symbolic, converged.
    """
    f = _psi_fn(psi, p)
    x_, y_ = _F(x), _F(y)
    u = (_F(f(x_, y_ + h)) - _F(f(x_, y_ - h))) / (2.0 * h) / rho  # ρu = ∂ψ/∂y
    v = -(_F(f(x_ + h, y_)) - _F(f(x_ - h, y_))) / (2.0 * h) / rho  # ρv = −∂ψ/∂x
    return _S(u), _S(v)


def velocity_from_streamfunction_2d_sym(psi_expr, x: sp.Symbol, y: sp.Symbol, rho=1):
    """Symbolic (u, v) = (∂ψ/∂y, −∂ψ/∂x)/ρ (§4.3). Book: §4.3. Label: symbolic."""
    psi_e = sp.sympify(psi_expr)
    return sp.diff(psi_e, y) / rho, -sp.diff(psi_e, x) / rho


def velocity_from_streamfunction_axisym(psi, R, z, rho: float = 1.0, h: float = 1e-5, **p):
    """Velocity of an axisymmetric flow from its Stokes stream function: ρu_R = −(1/R)∂ψ/∂z, ρu_z = (1/R)∂ψ/∂R.

    Book: §4.3, text after Fig. 4.1 (χ = −φ, ∇χ = −R⁻¹e_φ, ρu = −R⁻¹e_φ × ∇ψ).

    Parameters
    ----------
    psi : callable ψ(R, z) [kg/s] (or [m³/s] with rho = 1), or "axisym_uniform" (U);  R : distance from the axis [m],
    > h;  z [m];  rho [kg/m³];  h : step [m]

    Returns (u_R, u_z) [m/s].
    Assumptions: R > 0 (on the axis use the limit, e.g. the sympy version and L'Hôpital).
    Validation: V1 uniform stream ψ = ½UR² → (0, U); V2 sympy ∇·u = (1/R)∂(Ru_R)/∂R + ∂u_z/∂z ≡ 0. Label: analytic.
    """
    if isinstance(psi, str):
        name = psi
        f = lambda RR, zz: streamfunction_preset(name, zz, RR, **p)  # noqa: E731  presets take (z, R)
    else:
        f = psi
    R_, z_ = _F(R), _F(z)
    if np.any(R_ <= h):
        raise ValueError("velocity_from_streamfunction_axisym needs R > h > 0 (the axis R = 0 is a limit)")
    uR = -(_F(f(R_, z_ + h)) - _F(f(R_, z_ - h))) / (2.0 * h) / (R_ * rho)  # ρu_R = −(1/R)∂ψ/∂z
    uz = (_F(f(R_ + h, z_)) - _F(f(R_ - h, z_))) / (2.0 * h) / (R_ * rho)  # ρu_z = (1/R)∂ψ/∂R
    return _S(uR), _S(uz)


def velocity_from_streamfunction_axisym_sym(psi_expr, R: sp.Symbol, z: sp.Symbol, rho=1):
    """Symbolic (u_R, u_z) = (−(1/R)∂ψ/∂z, (1/R)∂ψ/∂R)/ρ (§4.3). Book: §4.3. Label: symbolic."""
    psi_e = sp.sympify(psi_expr)
    return -sp.diff(psi_e, z) / (R * rho), sp.diff(psi_e, R) / (R * rho)


def flux_along_path(psi, pts, rho: float = 1.0, n: int = 64, h: float = 1e-6, **p) -> float:
    """Flux ∫ρu·n ds across a polyline gate through the points ``pts`` (2, m), n = right-hand normal of the direction of
    travel on each segment — equals ψ(last) − ψ(first) for any path (§4.3).

    Book: §4.3 ("the flux between two streamlines is the difference of their ψ values"; ṁ = ∫χdψ). The velocity comes
    from ψ by central differences and is integrated with n Gauss–Legendre nodes per segment — the independent route.
    Parameters: psi (callable or preset name + **p); pts [m]; rho [kg/m³]; n; h [m]. Returns [kg/(m s)] (m²/s for ρ = 1).
    Validation: V1 = ψ(end) − ψ(start) to 1e-9 for bent gates. Label: analytic.
    """
    P = _F(pts)
    total = 0.0
    for k in range(P.shape[1] - 1):
        a, b = P[:, k], P[:, k + 1]
        L = float(np.hypot(*(b - a)))
        if L == 0.0:
            continue
        s, w = gauss_legendre_nodes(0.0, L, n)
        tv = (b - a) / L
        nv = np.array([tv[1], -tv[0]])  # right-hand normal: flux positive ⇔ ψ increases along the path
        U, V = velocity_from_streamfunction_2d(psi, a[0] + s * tv[0], a[1] + s * tv[1], 1.0, h, **p)
        total += float(np.sum((_F(U) * nv[0] + _F(V) * nv[1]) * w)) * rho  # ∫ ρu·n ds
    return total


def flux_between_streamlines(psi, p1, p2, n: int = 400, rho: float = 1.0, **p) -> float:
    """Flux ∫ρu·n ds across the straight gate p1 → p2 (n = (t_y, −t_x): from the left to the right of the direction of
    travel), computed by integrating the velocity — equal to ψ(p2) − ψ(p1) whatever the gate (§4.3 with χ = −z).

    Book: §4.3: ṁ = ∫ρu·n dA = ∮Ψ·ds = ∫χ dψ; in the plane the flux between two streamlines is the difference of their
    ψ values. Parameters: psi (callable ψ(x, y) or preset name + **p); p1, p2 (x, y) [m]; n Gauss–Legendre nodes;
    rho [kg/m³]. Returns flux per unit depth [kg/(m s)] ([m²/s] for ρ = 1). Scalar-callable (E2 parity).
    Validation: V1 = ψ(p2) − ψ(p1) to 1e-9 for the presets and any gate angle. Label: analytic.
    """
    return flux_along_path(psi, np.stack([_F(p1), _F(p2)], axis=1), rho, n, **p)


# ======================================================================================================================
# three-dimensional flow: two stream functions (Fig. 4.1)
# ======================================================================================================================
def mass_flux_from_vector_potential(Psi_exprs, coords) -> list:
    """Symbolic ρu = ∇ × Ψ, Eq. (4.12), in Cartesian coordinates (x, y, z). Returns [ρu, ρv, ρw]. Book: §4.3.
    Label: symbolic."""
    x, y, z = coords
    P = [sp.sympify(e) for e in Psi_exprs]
    return [sp.diff(P[2], y) - sp.diff(P[1], z), sp.diff(P[0], z) - sp.diff(P[2], x),
            sp.diff(P[1], x) - sp.diff(P[0], y)]  # Eq. (4.12)


def mass_flux_from_stream_functions(chi_expr, psi_expr, coords) -> list:
    """Symbolic ρu = ∇χ × ∇ψ (Eq. (4.12) with Ψ = χ∇ψ), Cartesian. Returns [ρu, ρv, ρw]. Book: §4.3. Label: symbolic."""
    gc = sp.Matrix([sp.diff(chi_expr, c) for c in coords])
    gp = sp.Matrix([sp.diff(psi_expr, c) for c in coords])
    return list(gc.cross(gp))  # ρu = ∇χ × ∇ψ


def stream_function_pair(name: str = "parabolic"):
    """A concrete 3-D pair of stream functions for Fig. 4.1 (our example): "parabolic" χ = y, ψ = z − x², so
    ρu = ∇χ × ∇ψ = (1, 0, 2x): streamlines z = x² + c in the planes y = const, ∇·(ρu) = 0.
    Returns (chi_expr, psi_expr, (x, y, z)). Book: §4.3, Fig. 4.1. Label: symbolic."""
    x, y, z = sp.symbols("x y z", real=True)
    if name != "parabolic":
        raise ValueError("only the 'parabolic' pair is defined")
    return y, z - x ** 2, (x, y, z)


def _pair_callables(chi, psi):
    if isinstance(chi, str):
        ce, pe, cs = stream_function_pair(chi)
        return sp.lambdify([cs], ce, "numpy"), sp.lambdify([cs], pe, "numpy")
    return chi, psi


def _grad_num(f, X, h):
    out = []
    for i in range(3):
        e = np.zeros((3,) + (1,) * (X.ndim - 1))
        e[i] = h
        out.append((_F(f(X + e)) + 0.0 * X[0] - (_F(f(X - e)) + 0.0 * X[0])) / (2 * h))
    return np.stack(out)


def stream_surface_check(chi, psi=None, x=None, h: float = 1e-6) -> dict:
    """Checks that ρu = ∇χ × ∇ψ conserves mass and lies in both stream surfaces: ∇·(ρu), ρu·∇χ, ρu·∇ψ (all 0).

    Book: §4.3, Eqs. (4.11)–(4.12) (div of a curl, Exercise 2.19; ∇χ ⟂ χ-surfaces). Two modes:
    * symbolic — chi, psi sympy expressions and x = the coordinate symbols → simplified expressions (≡ 0);
    * numeric — chi, psi callables of points X (3, N) (or chi = "parabolic") and x points (3, N) → arrays (≈ 0; the
      divergence by nested central differences with step h [m]).
    Returns dict(u_dot_grad_chi, u_dot_grad_psi, div). Label: symbolic / analytic.
    """
    if not isinstance(chi, str) and isinstance(chi, sp.Basic):
        coords = tuple(x)
        m = sp.Matrix(mass_flux_from_stream_functions(chi, psi, coords))
        gc = sp.Matrix([sp.diff(chi, c) for c in coords])
        gp = sp.Matrix([sp.diff(psi, c) for c in coords])
        return {"u_dot_grad_chi": sp.simplify(m.dot(gc)), "u_dot_grad_psi": sp.simplify(m.dot(gp)),
                "div": sp.simplify(sum(sp.diff(m[i], c) for i, c in enumerate(coords)))}
    cf, pf = _pair_callables(chi, psi)
    X = _F(x)
    if X.ndim == 1:
        X = X[:, None]
    m = lambda Y: np.cross(_grad_num(cf, Y, h), _grad_num(pf, Y, h), axis=0)  # noqa: E731
    M = m(X)
    div = 0.0
    for i in range(3):
        e = np.zeros((3, 1))
        e[i] = 1e3 * h
        div = div + (m(X + e)[i] - m(X - e)[i]) / (2e3 * h)
    return {"u_dot_grad_chi": _S(np.sum(M * _grad_num(cf, X, h), axis=0)),
            "u_dot_grad_psi": _S(np.sum(M * _grad_num(pf, X, h), axis=0)), "div": _S(div)}


def stream_tube_mass_flux(chi="parabolic", psi=None, a: float = 0.0, b: float = 1.0, c: float = 0.0, d: float = 1.0,
                          n: int = 64, patch: Callable | None = None, x0: float = 0.5, h: float = 1e-6):
    """Mass flux through the patch bounded by χ = a, b and ψ = c, d, computed by quadrature, vs (b − a)(d − c) (Fig. 4.1).

    Book: §4.3, the Stokes-theorem argument after (4.12): ṁ = ∫_A ρu·n dA = ∫_C χ dψ = b(d − c) + a(c − d)
    = (b − a)(d − c). The numeric value integrates (∇χ × ∇ψ)·(∂X/∂s × ∂X/∂r) over a parametric patch X(s, r),
    (s, r) ∈ [0, 1]² (gradients by central differences) — independent of the closed form.

    Parameters
    ----------
    chi, psi : "parabolic" (χ = y, ψ = z − x²; the patch is built on the plane x = x0: y ∈ [a, b], z = ψ + x0²) or
        callables of X (3, N) together with ``patch``
    a, b, c, d : the bounding values;  n : Gauss nodes per direction;  patch : callable (s, r) → X (3, …)
    x0 : plane of the default patch [m];  h : stencil step [m]

    Returns
    -------
    (numeric, closed_form) [kg/s] — the numeric sign follows ∂X/∂s × ∂X/∂r (the default patch is oriented along +x, the
    stream direction).
    Validation: V1 agreement to 1e-9 (default pair); ``ch04.stream_surface_example`` for a second pair. Label: analytic.
    """
    cf, pf = _pair_callables(chi, psi)
    if patch is None:
        if not (isinstance(chi, str) and chi == "parabolic"):
            raise ValueError("give a patch for user-supplied stream functions")

        def patch(s, r):
            y = a + (b - a) * _F(s)  # χ = y
            z = c + (d - c) * _F(r) + x0 ** 2  # ψ = z − x0²
            return np.stack([x0 + 0.0 * y, y, z])
    s, ws = gauss_legendre_nodes(0.0, 1.0, n)
    S, Rr = np.meshgrid(s, s, indexing="ij")
    S, Rr = S.ravel(), Rr.ravel()
    W = np.einsum("i,j->ij", ws, ws).ravel()
    X = _F(patch(S, Rr)).reshape(3, -1)
    dXs = (_F(patch(S + h, Rr)).reshape(3, -1) - _F(patch(S - h, Rr)).reshape(3, -1)) / (2 * h)
    dXr = (_F(patch(S, Rr + h)).reshape(3, -1) - _F(patch(S, Rr - h)).reshape(3, -1)) / (2 * h)
    dA = np.cross(dXs, dXr, axis=0)
    m = np.cross(_grad_num(cf, X, h), _grad_num(pf, X, h), axis=0)  # ρu = ∇χ × ∇ψ
    return float(np.sum(np.einsum("ik,ik->k", m, dA) * W)), float((b - a) * (d - c))
