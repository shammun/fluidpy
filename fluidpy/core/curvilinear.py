"""Vector-calculus operators in orthogonal curvilinear coordinates (cylindrical, spherical) from scale factors, in
sympy — the Appendix-B forms the book uses for Example 4.5 and later chapters.

Book: Kundu, Cohen & Dowling 5e, Ch. 4 §4.7, Example 4.5 ("Here we have used the results found in the Appendix B for
cylindrical coordinates"; rendered page chapters/pages/ch04/p148); Fig. 3.3 for the coordinate definitions.

Coordinates (q₁, q₂, q₃) and scale factors (h₁, h₂, h₃):
* "cartesian" (x, y, z): (1, 1, 1)
* "cylindrical" (R, φ, z): (1, R, 1)
* "spherical" (r, θ, φ), θ measured from +z (polar angle), φ azimuth: (1, r, r sin θ)

Vectors are lists of three **physical** components along the orthonormal unit vectors (e.g. (u_R, u_φ, u_z)).
General orthogonal-coordinate formulas (h = h₁h₂h₃):
grad f = Σ (1/h_i)∂f/∂q_i e_i; div A = (1/h)Σ ∂(A_i h/h_i)/∂q_i; (curl A)₁ = (1/h₂h₃)[∂(h₃A₃)/∂q₂ − ∂(h₂A₂)/∂q₃]
(cyclic); ∇²f = (1/h)Σ ∂((h/h_i²)∂f/∂q_i)/∂q_i; ∇²A = ∇(∇·A) − ∇×(∇×A); (u·∇)u = ∇(½|u|²) − u × (∇×u).
Unit vectors do not depend on t, so Du/Dt = ∂u/∂t + (u·∇)u component-wise in physical components.
"""
from __future__ import annotations

import sympy as sp

__all__ = ["SYSTEMS", "coordinates", "scale_factors", "gradient", "divergence", "curl", "laplacian",
           "vector_laplacian", "advective_acceleration", "material_acceleration", "strain_rate", "cross"]

SYSTEMS = ("cartesian", "cylindrical", "spherical")


def coordinates(system: str) -> tuple:
    """The coordinate symbols used by default: cartesian (x, y, z); cylindrical (R, phi, z), R > 0; spherical
    (r, theta, phi), r > 0, 0 < θ < π. Build fields with these symbols (sympy assumptions make them distinct from
    plain Symbol("R")). Label: analytic."""
    if system == "cartesian":
        return sp.symbols("x y z", real=True)
    if system == "cylindrical":
        return sp.Symbol("R", positive=True), sp.Symbol("phi", real=True), sp.Symbol("z", real=True)
    if system == "spherical":
        return sp.Symbol("r", positive=True), sp.Symbol("theta", positive=True), sp.Symbol("phi", real=True)
    raise ValueError(f"system must be one of {SYSTEMS}")


def scale_factors(system: str, coords=None) -> tuple:
    """Scale factors (h₁, h₂, h₃) of the system for the given coordinate symbols. Label: analytic."""
    q = coordinates(system) if coords is None else tuple(coords)
    if system == "cartesian":
        return (sp.Integer(1),) * 3
    if system == "cylindrical":
        return sp.Integer(1), q[0], sp.Integer(1)
    if system == "spherical":
        return sp.Integer(1), q[0], q[0] * sp.sin(q[1])
    raise ValueError(f"system must be one of {SYSTEMS}")


def _setup(system, coords):
    q = coordinates(system) if coords is None else tuple(coords)
    return q, scale_factors(system, q)


def _simp(e, simplify):
    return sp.simplify(e) if simplify else e


def gradient(f, system: str = "cylindrical", coords=None, simplify: bool = True) -> list:
    """∇f in physical components. Book: Appendix B via §4.7 Example 4.5. Label: symbolic."""
    q, h = _setup(system, coords)
    return [_simp(sp.diff(f, q[i]) / h[i], simplify) for i in range(3)]


def divergence(A, system: str = "cylindrical", coords=None, simplify: bool = True):
    """∇·A in physical components. Label: symbolic."""
    q, h = _setup(system, coords)
    H = h[0] * h[1] * h[2]
    return _simp(sum(sp.diff(A[i] * H / h[i], q[i]) for i in range(3)) / H, simplify)


def curl(A, system: str = "cylindrical", coords=None, simplify: bool = True) -> list:
    """∇×A in physical components. Label: symbolic."""
    q, h = _setup(system, coords)
    out = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        out.append(_simp((sp.diff(h[k] * A[k], q[j]) - sp.diff(h[j] * A[j], q[k])) / (h[j] * h[k]), simplify))
    return out


def laplacian(f, system: str = "cylindrical", coords=None, simplify: bool = True):
    """Scalar Laplacian ∇²f. Label: symbolic."""
    q, h = _setup(system, coords)
    H = h[0] * h[1] * h[2]
    return _simp(sum(sp.diff(H / h[i] ** 2 * sp.diff(f, q[i]), q[i]) for i in range(3)) / H, simplify)


def cross(a, b) -> list:
    """a × b of physical components in a right-handed orthonormal basis. Label: symbolic."""
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def vector_laplacian(A, system: str = "cylindrical", coords=None, simplify: bool = True) -> list:
    """Vector Laplacian ∇²A = ∇(∇·A) − ∇×(∇×A) in physical components (e.g. cylindrical: ∇²u_R − u_R/R² −
    (2/R²)∂u_φ/∂φ …). Label: symbolic."""
    gd = gradient(divergence(A, system, coords, False), system, coords, False)
    cc = curl(curl(A, system, coords, False), system, coords, False)
    return [_simp(gd[i] - cc[i], simplify) for i in range(3)]


def advective_acceleration(u, system: str = "cylindrical", coords=None, simplify: bool = True) -> list:
    """(u·∇)u = ∇(½|u|²) − u × (∇×u) in physical components (the centripetal −u_φ²/R and u_Ru_φ/R terms appear
    automatically). Book: §4.7 Example 4.5 left sides; Lamb's identity (4.68). Label: symbolic."""
    ke = sp.Rational(1, 2) * sum(c ** 2 for c in u)
    g = gradient(ke, system, coords, False)
    lam = cross(u, curl(u, system, coords, False))
    return [_simp(g[i] - lam[i], simplify) for i in range(3)]


def material_acceleration(u, t: sp.Symbol, system: str = "cylindrical", coords=None, simplify: bool = True) -> list:
    """Du/Dt = ∂u/∂t + (u·∇)u in physical components ((3.9) in curvilinear form). Label: symbolic."""
    adv = advective_acceleration(u, system, coords, False)
    return [_simp(sp.diff(u[i], t) + adv[i], simplify) for i in range(3)]


def strain_rate(u, system: str = "cylindrical", coords=None, simplify: bool = True) -> sp.Matrix:
    """Strain-rate tensor S in physical components of an orthogonal system:
    S_ii = (1/h_i)∂u_i/∂q_i + Σ_{k≠i} u_k/(h_ih_k) ∂h_i/∂q_k;  S_ij = ½[(h_j/h_i)∂(u_j/h_j)/∂q_i +
    (h_i/h_j)∂(u_i/h_i)/∂q_j] (i ≠ j). Cylindrical: S_RR = ∂u_R/∂R, S_φφ = (1/R)∂u_φ/∂φ + u_R/R,
    S_Rφ = ½[R∂(u_φ/R)/∂R + (1/R)∂u_R/∂φ]. Book: (3.12) in curvilinear form (Appendix B). Label: symbolic."""
    q, h = _setup(system, coords)
    S = sp.zeros(3, 3)
    for i in range(3):
        S[i, i] = sp.diff(u[i], q[i]) / h[i] + sum(u[k] / (h[i] * h[k]) * sp.diff(h[i], q[k]) for k in range(3) if k != i)
        for j in range(3):
            if j != i:
                S[i, j] = sp.Rational(1, 2) * (h[j] / h[i] * sp.diff(u[j] / h[j], q[i])
                                               + h[i] / h[j] * sp.diff(u[i] / h[i], q[j]))
    return S.applyfunc(sp.simplify) if simplify else S
