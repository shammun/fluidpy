"""Ideal (potential) flow elements and their superposition — plane flows through the complex potential, axisymmetric and
3-D flows through the Stokes stream function and the velocity potential, and the arbitrarily moving sphere.

Book: Kundu, Cohen & Dowling 5e, Ch. 6 — §6.2 (ψ, φ, (6.3)–(6.23)), §6.3 (elementary flows and superposition
(6.24)–(6.41)), §6.4 (complex potential (6.42)–(6.53)), §6.5 (Blasius and Kutta–Zhukhovsky (6.54)–(6.62)), §6.8
(axisymmetric flow (6.74)–(6.95)), §6.9 (moving sphere and added mass (6.96)–(6.109)). Every equation was transcribed from
the rendered page images (chapters/pages/ch06/p227–p267). Reused by Ch. 7 (wave potentials), Ch. 9 (outer flows),
Ch. 13 (ψ tools), Ch. 14 (airfoils, added mass), Ch. 16 (bubbles).

Conventions (analysis/ch06.md §9, curation §8)
----------------------------------------------
* Plane fields: ``z = x + 1j*y``; complex potential w = φ + iψ (6.42); complex velocity dw/dz = u − iv (6.45);
  u = ∂ψ/∂y, v = −∂ψ/∂x (6.3); u = ∂φ/∂x, v = ∂φ/∂y (6.10). Every plane element has ``w``, ``dwdz``, ``d2wdz2``.
* **Γ is counterclockwise-positive in every project function** ((6.6), (6.8), (6.47)): ``Vortex(Gamma)`` has
  u_θ = +Γ/2πr. The book's (6.36)–(6.40), (6.52), (6.61)–(6.62), (6.68) and Example 6.1 use a *clockwise* Γ (the flow's
  circulation is −Γ, so L = +ρUΓ). Body helpers therefore take keyword-only ``Gamma_cw=`` (the book's) **or**
  ``Gamma_ccw=`` (the project's), ``Gamma_ccw = −Gamma_cw``.
* 2-D doublet: dipole vector **d = Σ x_i m_i points from the sink to the source** (6.28); φ = −d·(x − x′)/2π|x − x′|²
  (6.29); the cylinder needs d = −2πUa² e_x; (6.49)'s scalar d is a dipole −d e_x (:meth:`Doublet.from_book_scalar`).
* Axisymmetric (§6.8): coordinates (R, z) with **z along the stream (horizontal)**, Stokes ψ [m³/s] with
  u_R = −(1/R)∂ψ/∂z, u_z = (1/R)∂ψ/∂R (6.75) — the sign of ``core.streamfunction``; spherical (r, θ) with θ from +z.
  3-D doublet "−d e_z" (6.88); the sphere has d = 2πa³U (6.89). §6.9: ξ = x − x_s, dipole d(t) = +2πa³u_s (6.97).
* Units SI; NaN inside bodies (each body flow carries an ``inside`` test); an element's own centre is singular by design.
* Book slips handled here: (6.104) prints −u_s/a³ inside the bracket — the gradient (6.103) at |ξ| = a needs +u_s/a³
  (``printed_bracket=True`` keeps the printed variant for wrong-variant tests); (6.108) keeps a stray dφ after the
  φ-integration (the result −(2/3)πρa³ du_s/dt is right and is what :func:`added_mass_sphere` returns).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Callable, NamedTuple, Sequence

import numpy as np
import sympy as sp

from ._util import as_scalar_if_0d

__all__ = [
    "log_branch", "power_branch", "as_complex_point",
    "Uniform", "Source", "Vortex", "Doublet", "Corner", "Constant", "element_from_spec",
    "ComplexFlow", "Flow", "FunctionFlow", "stream_function", "velocity_potential",
    "half_body", "cylinder", "gamma_ccw_from", "mirror", "circle_theorem",
    "BlasiusForce", "blasius_force", "blasius_force_report", "contour_crosses_body", "laurent_coefficients",
    "laplacian_residual", "normal_velocity_on", "normal_velocity_values", "far_field_check", "polar_velocity", "polar_velocity_sym",
    "pressure_coefficient", "polygon_signed_area",
    "AxisymUniform", "PointSource3D", "Doublet3D", "LineSource3D", "AxisymFlow", "sphere",
    "axisym_velocity_spherical", "axisym_velocity_spherical_sym", "sphere_potential_vector",
    "moving_sphere_potential", "moving_sphere_velocity", "moving_sphere_surface_velocity", "moving_sphere_dphidt",
    "moving_sphere_surface_pressure", "sphere_force_quadrature", "added_mass_sphere",
]

_S = as_scalar_if_0d
_TWO_PI = 2.0 * np.pi
_FOUR_PI = 4.0 * np.pi


def _F(a):
    return np.asarray(a, dtype=float)


def _C(a):
    return np.asarray(a, dtype=complex)


def as_complex_point(p) -> complex:
    """A point given as a complex number, a real number (on the x-axis) or an (x, y) pair → complex x + iy [m]."""
    if p is None:
        return 0j
    if isinstance(p, (tuple, list, np.ndarray)) and np.size(p) == 2:
        q = np.asarray(p, dtype=float).ravel()
        return complex(q[0], q[1])
    return complex(p)


def log_branch(zeta, cut: float = np.pi, theta_range: str = "(]"):
    """Complex logarithm ln|ζ| + iθ with the branch θ ∈ (cut − 2π, cut] (the cut is the ray at angle ``cut``);
    ``theta_range="[)"`` gives θ ∈ [cut − 2π, cut) instead (the half-body's θ ∈ [0, 2π) with cut = 2π).

    ``cut = π`` is numpy's principal branch. A potential with a logarithm (source, vortex) is multivalued; placing the
    cut inside a body (half-body: ``cut = 2π``) or away from the region of interest keeps ψ (or φ) continuous there.
    Book: §6.4 after (6.47)–(6.48) (θ′ = tan⁻¹((y − y′)/(x − x′)) is multivalued). Label: analytic.
    """
    z = _C(zeta)
    shift = float(cut) - np.pi
    th = np.angle(z * np.exp(-1j * shift)) + shift
    if theta_range == "[)":
        lo = float(cut) - _TWO_PI
        th = lo + np.mod(th - lo, _TWO_PI)
    elif theta_range != "(]":
        raise ValueError('theta_range must be "(]" or "[)"')
    with np.errstate(divide="ignore"):
        return np.log(np.abs(z)) + 1j * th


def power_branch(zeta, p: float, cut: float = np.pi):
    """ζ^p = exp(p·log ζ) on the branch of :func:`log_branch` (cut at angle ``cut``); 0^p = 0 for p > 0.
    Book: §6.4, (6.46). Label: analytic."""
    z = _C(zeta)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.exp(p * log_branch(z, cut))
    if p > 0:
        out = np.where(z == 0, 0.0, out)
    return out


# ======================================================================================================================
# plane elements (§6.2–6.4)
# ======================================================================================================================
class _Element:
    """Mixin: velocity, φ and ψ of one plane element from its complex potential."""
    kind = "element"

    @property
    def singular_points(self) -> tuple:
        """Points where w or dw/dz is singular (complex tuple): the centre of a source, vortex or doublet, the tip of a
        corner flow zⁿ with n < 2 (infinite speed or a non-analytic w there); empty otherwise."""
        if self.kind in ("source", "vortex", "doublet"):
            return (complex(self.z0),)
        if self.kind == "corner" and float(self.n) < 2.0:
            return (0j,)
        return ()

    def w(self, z):  # pragma: no cover - overridden
        raise NotImplementedError

    def dwdz(self, z):  # pragma: no cover - overridden
        raise NotImplementedError

    def d2wdz2(self, z):  # pragma: no cover - overridden
        raise NotImplementedError

    def velocity(self, x, y):
        """(u, v) [m/s] = (Re, −Im) of dw/dz, Eq. (6.45)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            q = self.dwdz(_F(x) + 1j * _F(y))
        return _S(q.real), _S(-q.imag)

    def phi(self, x, y):
        """Velocity potential φ = Re w [m²/s], Eq. (6.42)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return _S(self.w(_F(x) + 1j * _F(y)).real)

    def psi(self, x, y):
        """Stream function ψ = Im w [m²/s], Eq. (6.42)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return _S(self.w(_F(x) + 1j * _F(y)).imag)


@dataclass
class Uniform(_Element):
    """Uniform stream (U, V): ψ = −Vx + Uy (6.7), φ = Ux + Vy (6.14), w = (U − iV)z.

    Book: §6.2, Eqs. (6.7), (6.14). Parameters: U, V [m/s] (x and y components). Validation (planned): V1 Re/Im of w
    equal (6.14)/(6.7). Label: analytic.
    """
    U: float = 1.0
    V: float = 0.0
    kind = "uniform"

    @property
    def strength(self) -> float:
        return float(np.hypot(self.U, self.V))

    def w(self, z):
        return (self.U - 1j * self.V) * _C(z)  # Eq. (6.7), (6.14): φ + iψ = (Ux + Vy) + i(Uy − Vx)

    def dwdz(self, z):
        return np.full(np.shape(z), self.U - 1j * self.V, dtype=complex)

    def d2wdz2(self, z):
        return np.zeros(np.shape(z), dtype=complex)


@dataclass
class Source(_Element):
    """Plane point source (m > 0) or sink (m < 0) at z0: φ = (m/2π) ln r′ (6.15), w = (m/2π) ln(z − z0) (6.48).

    m is the volume flow rate per unit depth [m²/s]; u_r = m/2πr′. ``cut_angle`` places the branch cut of θ′ (ψ = mθ′/2π
    jumps by m across it; ``theta_range`` "(]" or "[)" closes the branch at either end). Book: §6.2 (6.13), (6.15); §6.3 source velocities; §6.4 (6.48).
    Validation (planned): V1 flux of ∇φ through any circle = m; V2 Laplace away from z0. Label: analytic.
    """
    m: float = 1.0
    z0: complex = 0j
    cut_angle: float = np.pi
    theta_range: str = "(]"
    kind = "source"

    def __post_init__(self):
        self.z0 = as_complex_point(self.z0)

    def w(self, z):
        return self.m / _TWO_PI * log_branch(_C(z) - self.z0, self.cut_angle, self.theta_range)  # Eq. (6.48)

    def dwdz(self, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.m / (_TWO_PI * (_C(z) - self.z0))

    def d2wdz2(self, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -self.m / (_TWO_PI * (_C(z) - self.z0) ** 2)


@dataclass
class Vortex(_Element):
    """Plane point (ideal) vortex of **counterclockwise** circulation Γ at z0: ψ = −(Γ/2π) ln r′ (6.8),
    w = −(iΓ/2π) ln(z − z0) (6.47), u_θ = +Γ/2πr′.

    Γ [m²/s] counterclockwise positive (the book's (6.36)–(6.40) Γ is clockwise: pass −Γ here). ``cut_angle`` places the jump
    of φ = Γθ′/2π. Book: §6.2 (6.6), (6.8); §6.3 vortex velocities; §6.4 (6.47).
    Validation (planned): V1 circulation on any circle = +Γ (wrong variant ψ = +(Γ/2π) ln r fails); parity with
    ``core.biot_savart.point_vortex_velocity``. Label: analytic.
    """
    Gamma: float = 1.0
    z0: complex = 0j
    cut_angle: float = np.pi
    kind = "vortex"

    def __post_init__(self):
        self.z0 = as_complex_point(self.z0)

    def w(self, z):
        return -1j * self.Gamma / _TWO_PI * log_branch(_C(z) - self.z0, self.cut_angle)  # Eq. (6.47)

    def dwdz(self, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -1j * self.Gamma / (_TWO_PI * (_C(z) - self.z0))

    def d2wdz2(self, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return 1j * self.Gamma / (_TWO_PI * (_C(z) - self.z0) ** 2)


@dataclass
class Doublet(_Element):
    """Plane doublet of dipole vector d = (d_x, d_y) at z0 (d points from the sink to the source, (6.28)):
    φ = −d·(x − x′)/2π|x − x′|² (6.29), w = −D/2π(z − z0) with D = d_x + i d_y.

    ``d_vec`` [m³/s]. The book's (6.49) w = d/2π(z − z′) is a dipole −d e_x: use :meth:`from_book_scalar`. The cylinder
    of radius a in a stream U needs d = −2πUa² e_x (6.33). Book: §6.3 (6.28)–(6.29), §6.4 (6.49).
    Validation (planned): V1 limit of the source–sink pair (error ∝ ε²); direction (wrong variant d → −d fails).
    Label: analytic.
    """
    d_vec: tuple = (-1.0, 0.0)
    z0: complex = 0j
    kind = "doublet"

    def __post_init__(self):
        dv = np.asarray(self.d_vec, dtype=float).ravel()
        if dv.size != 2:
            raise ValueError("d_vec must be a 2-vector (d_x, d_y) [m³/s]")
        self.d_vec = (float(dv[0]), float(dv[1]))
        self.z0 = as_complex_point(self.z0)

    @classmethod
    def from_book_scalar(cls, d: float, z0=0j) -> "Doublet":
        """(6.49): w = d/2π(z − z′) is the doublet with dipole strength −d e_x."""
        return cls((-float(d), 0.0), z0)

    @property
    def D(self) -> complex:
        return complex(self.d_vec[0], self.d_vec[1])

    def w(self, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -self.D / (_TWO_PI * (_C(z) - self.z0))  # Eq. (6.29) as Re w; (6.49) for d_vec = (−d, 0)

    def dwdz(self, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.D / (_TWO_PI * (_C(z) - self.z0) ** 2)

    def d2wdz2(self, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -self.D / (np.pi * (_C(z) - self.z0) ** 3)


@dataclass
class Corner(_Element):
    """Flow in a corner of angle α = π/n: w = A (z e^{−i·rotate})ⁿ (6.46) (n ≥ ½); for rotate = 0,
    dw/dz = nA z^{n−1} = (Aπ/α) z^{(π−α)/α}.

    A [m^{2−n}/s] may be complex (A → −iA turns (6.24) into (6.25)); ``rotate`` [rad] turns the whole wedge (walls at
    θ = rotate and θ = rotate + α; n = 2 with rotate = π/4 is the 45° turn relating (6.24) and (6.25)). n = 2 is the
    stagnation flow ψ = 2Axy (6.24), φ = A(x² − y²) (6.27); n = ½ the flow round a semi-infinite plate. The branch cut of
    zⁿ is at ``cut_angle`` in the wedge's own frame (default: the middle of the non-fluid sector, π/n + (2π − π/n)/2, so
    the wedge is continuous). Book: §6.3 (6.24)–(6.27), §6.4 (6.46). Validation (planned): V1 ψ = 0 on both walls;
    |dw/dz| ∝ r^{n−1}. Label: analytic.
    """
    A: complex = 1.0
    n: float = 2.0
    cut_angle: float | None = None
    rotate: float = 0.0
    kind = "corner"

    def __post_init__(self):
        if float(self.n) < 0.5:
            raise ValueError("corner exponent n must be ≥ 1/2 (wedge angle α = π/n ≤ 2π)")
        if self.cut_angle is None:
            self.cut_angle = np.pi + 0.5 * self.alpha

    @property
    def alpha(self) -> float:
        """Corner (wedge) angle α = π/n [rad]."""
        return float(np.pi / self.n)

    def _zr(self, z):
        return _C(z) * np.exp(-1j * float(self.rotate))

    def w(self, z):
        return self.A * power_branch(self._zr(z), self.n, self.cut_angle)  # Eq. (6.46)

    def dwdz(self, z):
        e = np.exp(-1j * float(self.rotate))
        return self.n * self.A * power_branch(self._zr(z), self.n - 1.0, self.cut_angle) * e

    def d2wdz2(self, z):
        e = np.exp(-1j * float(self.rotate))
        return self.n * (self.n - 1.0) * self.A * power_branch(self._zr(z), self.n - 2.0, self.cut_angle) * e ** 2


@dataclass
class Constant(_Element):
    """A constant added to w (shifts φ and ψ, no velocity) — e.g. the ln a of (6.36)/(6.52). Label: analytic."""
    c: complex = 0j
    kind = "constant"

    def w(self, z):
        return np.full(np.shape(z), complex(self.c), dtype=complex)

    def dwdz(self, z):
        return np.zeros(np.shape(z), dtype=complex)

    def d2wdz2(self, z):
        return np.zeros(np.shape(z), dtype=complex)


def element_from_spec(spec):
    """Build an element from a plain spec (for explainer parity and the notebook):
    ("uniform", U[, V]) · ("source", m[, z0]) · ("sink", m[, z0]) (a source of −m) · ("vortex", Γ_ccw[, z0]) ·
    ("doublet", d_x, d_y[, z0]) · ("corner", A, n) · ("constant", c); a dict with key "kind" works too; an element is
    returned unchanged. Label: analytic."""
    if isinstance(spec, _Element):
        return spec
    if isinstance(spec, dict):
        s = dict(spec)
        k = s.pop("kind")
        return {"uniform": Uniform, "source": Source, "vortex": Vortex, "doublet": Doublet, "corner": Corner,
                "constant": Constant}[k](**s)
    k, *a = spec
    if k == "uniform":
        return Uniform(*a)
    if k == "source":
        return Source(*a)
    if k == "sink":
        return Source(-float(a[0]), *a[1:])
    if k == "vortex":
        return Vortex(*a)
    if k == "doublet":
        return Doublet((a[0], a[1]), *a[2:])
    if k == "corner":
        return Corner(*a)
    if k == "constant":
        return Constant(*a)
    raise ValueError(f"unknown element kind {k!r}")


# ======================================================================================================================
# flows (sums of elements, or any analytic w)
# ======================================================================================================================
class ComplexFlow:
    """Common API of a plane ideal flow given by an analytic complex potential w(z).

    Subclasses provide ``w(z)``, ``dwdz(z)`` and optionally ``d2wdz2(z)``, an ``inside(x, y)`` body test (NaN there) and
    the free-stream velocity ``u_inf = (U, V)`` [m/s] used by :meth:`cp` and :meth:`pressure`.
    """
    inside: Callable | None = None
    u_inf: tuple = (0.0, 0.0)

    # --- basic fields ---------------------------------------------------------------------------------------------
    def _mask(self, x, y):
        if self.inside is None:
            return np.zeros(np.broadcast(_F(x), _F(y)).shape, dtype=bool)
        return np.asarray(self.inside(_F(x), _F(y)), dtype=bool)

    def complex_velocity(self, x, y):
        """dw/dz = u − iv (6.45) at (x, y); NaN inside the body."""
        with np.errstate(divide="ignore", invalid="ignore"):
            q = _C(self.dwdz(_F(x) + 1j * _F(y)))
        return np.where(self._mask(x, y), np.nan + 0j, q)

    def velocity(self, x, y):
        """(u, v) [m/s] from dw/dz = u − iv, Eq. (6.45); NaN inside the body."""
        q = self.complex_velocity(x, y)
        return _S(q.real), _S(-q.imag)

    def phi(self, x, y):
        """φ = Re w [m²/s] (6.42); NaN inside the body. φ of a vortex jumps across its branch cut."""
        with np.errstate(divide="ignore", invalid="ignore"):
            w = _C(self.w(_F(x) + 1j * _F(y)))
        return _S(np.where(self._mask(x, y), np.nan, w.real))

    def psi(self, x, y):
        """ψ = Im w [m²/s] (6.42); NaN inside the body. ψ of a source jumps by m across its branch cut."""
        with np.errstate(divide="ignore", invalid="ignore"):
            w = _C(self.w(_F(x) + 1j * _F(y)))
        return _S(np.where(self._mask(x, y), np.nan, w.imag))

    def velocity_polar(self, r, theta):
        """(u_r, u_θ) [m/s] at polar (r, θ) about the origin (θ from +x): u_r = u cos θ + v sin θ,
        u_θ = −u sin θ + v cos θ — e.g. (6.34), (6.37) for the cylinder."""
        th = _F(theta)
        u, v = self.velocity(_F(r) * np.cos(th), _F(r) * np.sin(th))
        u, v = _F(u), _F(v)
        return _S(u * np.cos(th) + v * np.sin(th)), _S(-u * np.sin(th) + v * np.cos(th))

    def speed(self, x, y):
        """|u| [m/s]."""
        q = self.complex_velocity(x, y)
        return _S(np.abs(q))

    @property
    def U_inf(self) -> float:
        """Free-stream speed |u∞| [m/s]."""
        return float(np.hypot(*self.u_inf))

    def cp(self, x, y, U: float | None = None):
        """Pressure coefficient C_p = 1 − |u|²/U² (6.32) (U defaults to the free-stream speed). Book: §6.3."""
        U_ = self.U_inf if U is None else float(U)
        if U_ == 0.0:
            raise ValueError("C_p needs a reference speed U > 0 (quiescent far field: pass U)")
        s = _F(self.speed(x, y))
        return _S(1.0 - s ** 2 / U_ ** 2)  # Eq. (6.32)

    def pressure(self, x, y, rho: float = 1.0, p_inf: float = 0.0, U: float | None = None):
        """Steady Bernoulli everywhere (6.18) with const = p∞ + ρU²/2: p = p∞ + ½ρ(U² − |u|²) [Pa]. Book: §6.2–6.3."""
        U_ = self.U_inf if U is None else float(U)
        s = _F(self.speed(x, y))
        return _S(p_inf + 0.5 * rho * (U_ ** 2 - s ** 2))  # Eq. (6.18)

    # --- stagnation points ------------------------------------------------------------------------------------------
    def _d2(self, z):
        d2 = getattr(self, "d2wdz2", None)
        if d2 is not None:
            out = d2(z)
            if out is not None:
                return _C(out)
        h = 1e-6 * max(1.0, float(np.max(np.abs(z))))
        return (_C(self.dwdz(z + h)) - _C(self.dwdz(z - h))) / (2.0 * h)

    def stagnation_points(self, guesses=None, box=None, tol: float = 1e-12, n: int = 9,
                          max_iter: int = 80, keep_inside: bool = False) -> np.ndarray:
        """Stagnation points (dw/dz = 0) by complex Newton z ← z − f/f′ from ``guesses`` (or an n × n grid over
        ``box`` = (x0, x1, y0, y1), default (−3, 3, −3, 3) m); converged roots are de-duplicated and sorted by real then imaginary part.

        Returns a complex array [m]. Points strictly inside the body are dropped unless ``keep_inside``.
        Book: §6.3 (stagnation points of the half-body, (6.38)); analytic d²w/dz² when available.
        Validation (planned): V1 half-body −m/2πU, cylinder sin θ = −Γ/4πaU. Label: analytic.
        """
        box = (-3.0, 3.0, -3.0, 3.0) if box is None else box
        if guesses is None:
            gx = np.linspace(box[0], box[1], n)
            gy = np.linspace(box[2], box[3], n)
            G = (gx[None, :] + 1j * gy[:, None]).ravel()
        else:
            G = np.atleast_1d(_C(guesses)).ravel()
        scale = max(self.U_inf, 1e-300)
        L = max(abs(box[1] - box[0]), abs(box[3] - box[2]), 1.0)
        roots = []
        with np.errstate(all="ignore"):
            for z in G:
                zz = complex(z)
                for _ in range(max_iter):
                    f = complex(np.asarray(self.dwdz(np.array([zz])))[0])
                    fp = complex(self._d2(np.array([zz]))[0])
                    if not np.isfinite(f) or not np.isfinite(fp) or fp == 0:
                        zz = np.nan
                        break
                    step = f / fp
                    zz = zz - step
                    if abs(step) < 1e-14 * max(1.0, abs(zz)):
                        break
                if zz is np.nan or not np.isfinite(zz):
                    continue
                f = complex(np.asarray(self.dwdz(np.array([zz])))[0])
                if abs(f) > tol * scale * 1e3 and abs(f) > 1e-10 * scale:
                    continue
                if abs(zz) > 1e3 * L:
                    continue
                if not keep_inside and self.inside is not None and bool(np.asarray(self.inside(zz.real, zz.imag))):
                    continue
                if all(abs(zz - r) > 1e-8 * L for r in roots):
                    roots.append(zz)
        roots.sort(key=lambda c: (round(c.real, 9), round(c.imag, 9)))
        return np.array(roots, dtype=complex)


class Flow(ComplexFlow):
    """Superposition of plane elements (C03): w = Σ w_k, dw/dz = Σ dw_k/dz — Laplace is linear, so the sum is a flow;
    its boundary conditions belong to the sum (the body is the streamline the sum makes impermeable).

    Parameters: ``elements`` (list of :class:`Uniform`, :class:`Source`, :class:`Vortex`, :class:`Doublet`,
    :class:`Corner`, :class:`Constant` or specs for :func:`element_from_spec`), ``inside`` (optional body test
    (x, y) → bool), ``label``. The free stream is the sum of the uniform elements.
    Book: §6.2 (superposition after (6.17)), §6.3, §6.4. Validation (planned): V1 linearity to round-off. Label: analytic.
    """

    def __init__(self, elements: Sequence, inside: Callable | None = None, label: str = ""):
        self.elements = [element_from_spec(e) for e in elements]
        self.inside = inside
        self.label = label

    @property
    def u_inf(self) -> tuple:
        U = sum(e.U for e in self.elements if isinstance(e, Uniform))
        V = sum(e.V for e in self.elements if isinstance(e, Uniform))
        return (float(U), float(V))

    def w(self, z):
        z = _C(z)
        out = np.zeros(z.shape, dtype=complex)
        for e in self.elements:
            out = out + e.w(z)
        return out

    def dwdz(self, z):
        z = _C(z)
        out = np.zeros(z.shape, dtype=complex)
        for e in self.elements:
            out = out + e.dwdz(z)
        return out

    def d2wdz2(self, z):
        z = _C(z)
        out = np.zeros(z.shape, dtype=complex)
        for e in self.elements:
            out = out + e.d2wdz2(z)
        return out

    def __add__(self, other):
        o = other.elements if isinstance(other, Flow) else [other]
        return Flow(self.elements + list(o), self.inside, self.label)

    @property
    def net_source(self) -> float:
        """Σ m of the sources and sinks [m²/s] (zero ⇒ a closed body, §6.3)."""
        return float(sum(e.m for e in self.elements if isinstance(e, Source)))

    @property
    def circulation_ccw(self) -> float:
        """Σ Γ of the vortices [m²/s], counterclockwise positive."""
        return float(sum(e.Gamma for e in self.elements if isinstance(e, Vortex)))

    @property
    def singular_points(self) -> np.ndarray:
        return np.array([e.z0 for e in self.elements if hasattr(e, "z0")], dtype=complex)

    def contributions(self, x, y) -> list:
        """Per-element velocity at (x, y): list of dict(kind, u, v) — the term bars of E1 (their sum is
        :meth:`velocity`). Label: analytic."""
        out = []
        for e in self.elements:
            u, v = e.velocity(x, y)
            out.append({"kind": e.kind, "u": u, "v": v})
        return out


class FunctionFlow(ComplexFlow):
    """A plane flow given directly by callables w(z), dw/dz(z) (and optionally d²w/dz²) — mapped flows (§6.6),
    circle-theorem flows, any analytic function. ``u_inf`` = free-stream (U, V) [m/s]. Label: analytic."""

    def __init__(self, w: Callable, dwdz: Callable, d2wdz2: Callable | None = None, inside: Callable | None = None,
                 u_inf=(0.0, 0.0), label: str = ""):
        self._w, self._dw, self._d2fn = w, dwdz, d2wdz2
        self.inside = inside
        self.u_inf = (float(u_inf[0]), float(u_inf[1]))
        self.label = label

    def w(self, z):
        return _C(self._w(_C(z)))

    def dwdz(self, z):
        return _C(self._dw(_C(z)))

    def d2wdz2(self, z):
        if self._d2fn is None:
            return None
        return _C(self._d2fn(_C(z)))


def stream_function(flow, x, y):
    """ψ(x, y) [m²/s] of a plane flow (Im w, (6.42)); alias of ``flow.psi``. Book: §6.2 (6.3). Label: analytic."""
    return flow.psi(x, y)


def velocity_potential(flow, x, y):
    """φ(x, y) [m²/s] of a plane flow (Re w, (6.42)); alias of ``flow.phi``. Book: §6.2 (6.10). Label: analytic."""
    return flow.phi(x, y)


# ======================================================================================================================
# the chapter's bodies (§6.3–6.4)
# ======================================================================================================================
def half_body(U: float = 1.0, m: float = 2.0 * np.pi) -> Flow:
    """Rankine half-body: uniform stream U + source m at the origin, φ = Ux + (m/2π) ln r (6.30),
    ψ = Ur sin θ + mθ/2π (6.31), w = Uz + (m/2π) ln z (6.50).

    The source's cut lies along the positive x-axis inside the body, θ ∈ [0, 2π) (numpy's principal θ ∈ (−π, π] is the
    wrong variant: its cut is the upstream axis, where ψ would jump from m/2 to −m/2), so ψ = m/2 on the whole body and
    on the upstream axis. Stagnation point x = −a = −m/2πU; half-width h = m(π − θ)/2πU → h_max = m/2U.
    Parameters: U [m/s] > 0, m [m²/s] > 0. Returns a :class:`Flow` with ``inside`` = the body.
    Book: §6.3, Eqs. (6.30)–(6.31), (6.50), Fig. 6.7. Validation (planned): V1 ψ = m/2 on the body both halves;
    stagnation −m/2πU; V4 mass balance 2Uh_max = m. Label: analytic.
    """
    U, m = float(U), float(m)

    def inside(x, y):
        yy = np.abs(_F(y))
        th = np.arctan2(yy, _F(x))
        psi_u = U * yy + m * th / _TWO_PI
        return psi_u < 0.5 * m * (1.0 - 1e-12)

    return Flow([Uniform(U), Source(m, 0j, cut_angle=_TWO_PI, theta_range="[)")], inside=inside, label="half-body")


def gamma_ccw_from(Gamma_cw: float | None = None, Gamma_ccw: float | None = None) -> float:
    """Resolve the two circulation conventions: returns the project's counterclockwise Γ.

    ``Gamma_cw`` is the book's clockwise Γ of (6.36)–(6.40), (6.52), (6.68); ``Gamma_ccw`` the project's. Exactly one
    may be given (none → 0); Γ_ccw = −Γ_cw. Label: analytic."""
    if Gamma_cw is not None and Gamma_ccw is not None:
        raise ValueError("give Gamma_cw (book, clockwise) or Gamma_ccw (project, counterclockwise), not both")
    if Gamma_ccw is not None:
        return float(Gamma_ccw)
    return -float(Gamma_cw) if Gamma_cw is not None else 0.0


def cylinder(U: float = 1.0, a: float = 1.0, *, Gamma_cw: float | None = None, Gamma_ccw: float | None = None,
             center=0j) -> Flow:
    """Circular cylinder of radius a in a stream U with optional circulation:
    w = U(z + a²/z) + (iΓ_cw/2π) ln(z/a) (6.51)–(6.52), ψ = U(r − a²/r) sin θ + (Γ_cw/2π) ln(r/a) (6.33), (6.36).

    Built as uniform stream + doublet d = −2πUa² e_x + vortex (Γ_ccw = −Γ_cw) + the constant that keeps ψ = 0 on r = a.
    Parameters: U [m/s]; a [m] > 0; keyword-only ``Gamma_cw`` (book, clockwise) or ``Gamma_ccw`` (project) [m²/s];
    ``center`` of the cylinder. Returns a :class:`Flow` with ``inside`` = the disc.
    Book: §6.3 (6.33)–(6.37), §6.4 (6.51)–(6.52). Validation (planned): V1 parity with ``ch03.cylinder_flow``; u·n = 0 on
    r = a for every Γ; sign pinned (Γ_cw > 0 moves the stagnation points down, L = +ρUΓ_cw). Label: analytic.
    """
    U, a = float(U), float(a)
    if a <= 0:
        raise ValueError("cylinder radius a must be > 0")
    c = as_complex_point(center)
    G = gamma_ccw_from(Gamma_cw, Gamma_ccw)
    els = [Uniform(U), Doublet((-_TWO_PI * U * a ** 2, 0.0), c)]
    if G != 0.0:
        els += [Vortex(G, c), Constant(1j * G / _TWO_PI * np.log(a))]  # −(iΓ_ccw/2π)ln(z/a) = (iΓ_cw/2π)ln(z/a)

    def inside(x, y):
        return (_F(x) - c.real) ** 2 + (_F(y) - c.imag) ** 2 < a ** 2 * (1.0 - 1e-12)

    return Flow(els, inside=inside, label="cylinder")


def mirror(elements, wall: str = "y=0") -> list:
    """Method of images: the elements plus their mirror images in a straight wall (Figs. 6.14–6.15).

    Rules (§6.3): vorticity is mirrored with the **opposite** sign (ψ₂ = ψ₁(x, y) − ψ₁(x, −y), ψ₂ = 0 on the wall);
    sources with the **same** sign (φ₂ = φ₁(x, y) + φ₁(x, −y), no normal velocity on the wall); a doublet's vector is
    reflected (it is a source–sink pair). A uniform stream parallel to the wall is kept once; one crossing the wall,
    a corner flow, or an element on the wall raises. ``wall`` = "y=0" or "x=0".
    Returns the list (originals first). Book: §6.3, Figs. 6.14–6.16, (6.41). Validation (planned): V1 zero normal
    velocity on the wall (wrong variant: same-sign vortex image fails). Label: analytic.
    """
    els = [element_from_spec(e) for e in (elements.elements if isinstance(elements, Flow) else elements)]
    if wall not in ("y=0", "x=0"):
        raise ValueError('wall must be "y=0" or "x=0"')

    def refl(z0):
        return complex(z0.real, -z0.imag) if wall == "y=0" else complex(-z0.real, z0.imag)

    out, imgs = list(els), []
    for e in els:
        if isinstance(e, Uniform):
            if (wall == "y=0" and e.V != 0.0) or (wall == "x=0" and e.U != 0.0):
                raise ValueError("a uniform stream crossing the wall has no image")
            continue
        if isinstance(e, Constant):
            continue
        if isinstance(e, Corner):
            raise ValueError("corner flows are not mirrored")
        if (wall == "y=0" and e.z0.imag == 0.0) or (wall == "x=0" and e.z0.real == 0.0):
            raise ValueError("an element on the wall is its own image")
        if isinstance(e, Source):
            imgs.append(Source(e.m, refl(e.z0)))
        elif isinstance(e, Vortex):
            imgs.append(Vortex(-e.Gamma, refl(e.z0)))
        elif isinstance(e, Doublet):
            dx, dy = e.d_vec
            imgs.append(Doublet((dx, -dy) if wall == "y=0" else (-dx, dy), refl(e.z0)))
    return out + imgs


def circle_theorem(flow, a: float, center=0j) -> FunctionFlow:
    """Milne-Thomson circle theorem (our addition, generalising the book's images to circles, §6.3 end):
    W(z) = w(z) + conj(w(c + a²/conj(z − c))) makes |z − c| = a a streamline when w has no singularity inside it.

    ``flow`` is any object with ``w``/``dwdz`` (e.g. a :class:`Flow` of sources and vortices outside the circle) — then
    a :class:`FunctionFlow` is returned — or a bare callable w(z), for which the callable W(z) is returned.
    dW/dz = w′(z) − (a²/(z − c)²)·conj(w′(c + a²/conj(z − c))). Returns a :class:`FunctionFlow` (NaN inside).
    Book: §6.3 (images extend to circular boundaries; Exercises 5.14, 6.26). Validation (planned): V1 u·n = 0 on the
    circle; parity with ``core.biot_savart.circle_image_system`` for vortices. Label: analytic.
    """
    a = float(a)
    c = as_complex_point(center)
    if callable(flow) and not hasattr(flow, "w"):  # a bare w(z): return the callable W(z) (Part C form)
        wf = flow
        return lambda z: _C(wf(z)) + np.conj(_C(wf(c + a ** 2 / np.conj(_C(z) - c))))

    def img(z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return c + a ** 2 / np.conj(_C(z) - c)

    def W(z):
        return _C(flow.w(z)) + np.conj(_C(flow.w(img(z))))

    def dW(z):
        z = _C(z)
        with np.errstate(divide="ignore", invalid="ignore"):
            return _C(flow.dwdz(z)) - a ** 2 / (z - c) ** 2 * np.conj(_C(flow.dwdz(img(z))))

    def inside(x, y):
        return (_F(x) - c.real) ** 2 + (_F(y) - c.imag) ** 2 < a ** 2 * (1.0 - 1e-12)

    return FunctionFlow(W, dW, None, inside, getattr(flow, "u_inf", (0.0, 0.0)), label="circle theorem")


# ======================================================================================================================
# forces (§6.5)
# ======================================================================================================================
class BlasiusForce(NamedTuple):
    """Force per unit depth on the body [N/m]: D along x (drag), L along y (lift)."""
    D: float
    L: float


def polygon_signed_area(pts) -> float:
    """Signed area of a closed polygon [m²] (> 0 when traversed counterclockwise, the orientation of (6.56))."""
    z = _contour(pts)
    return float(0.5 * np.sum(z.real * np.roll(z.imag, -1) - np.roll(z.real, -1) * z.imag))


def _contour(pts) -> np.ndarray:
    p = np.asarray(pts)
    if np.iscomplexobj(p):
        return p.ravel().astype(complex)
    p = _F(p)
    if p.ndim == 2 and p.shape[0] == 2:
        return p[0] + 1j * p[1]
    if p.ndim == 2 and p.shape[1] == 2:
        return p[:, 0] + 1j * p[:, 1]
    raise ValueError("contour: complex array or (2, N) / (N, 2) real array of vertices")


def _dwdz_of(obj) -> Callable:
    return obj.dwdz if hasattr(obj, "dwdz") else obj


def blasius_force(flow_or_dwdz, R: float | None = None, contour=None, rho: float = 1.0, n: int = 256,
                  center=0j) -> BlasiusForce:
    """Blasius theorem (6.60): D − iL = (iρ/2)∮_C (dw/dz)² dz on any counterclockwise contour C enclosing the body
    (and no other singularity of (dw/dz)²).

    Either a circle (``R`` [m], ``center``; periodic trapezoid with ``n`` nodes — exponentially convergent for analytic
    integrands) or a closed polygon ``contour`` (complex array or (2, N)/(N, 2) vertices, last ≠ first; trapezoid on the
    segments — second order; must be counterclockwise, signed area > 0, else ValueError).
    Parameters: flow (object with ``dwdz``) or a callable dw/dz(z) [m/s]; ρ [kg/m³].
    Returns :class:`BlasiusForce` (D, L) [N/m] on the body. Book: §6.5, Eqs. (6.57)–(6.60).
    Validation (planned): V1 cylinder (0, ρUΓ_cw) for every R > a; V3 trapezoid rate; Kutta–Zhukhovsky on the ellipse.
    Label: analytic.
    """
    f = _dwdz_of(flow_or_dwdz)
    if contour is None:
        if R is None:
            raise ValueError("give a circle radius R or a closed contour")
        c = as_complex_point(center)
        th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
        e = np.exp(1j * th)
        z = c + float(R) * e
        dz = 1j * float(R) * e * (_TWO_PI / int(n))
        I = np.sum(_C(f(z)) ** 2 * dz)
    else:
        z = _contour(contour)
        if polygon_signed_area(z) <= 0:
            raise ValueError("the contour must be counterclockwise (signed area > 0), the orientation of (6.56)")
        z2 = np.roll(z, -1)
        g = _C(f(z)) ** 2
        I = np.sum(0.5 * (g + np.roll(g, -1)) * (z2 - z))
    DmiL = 0.5j * rho * I  # Eq. (6.60)
    return BlasiusForce(float(DmiL.real), float(-DmiL.imag))


def contour_crosses_body(inside: Callable | None, contour) -> bool:
    """True if any vertex of the contour lies inside the body (Blasius then integrates through the body's
    singularities and is invalid — E5's warning). Label: analytic."""
    if inside is None:
        return False
    z = _contour(contour)
    return bool(np.any(np.asarray(inside(z.real, z.imag))))


def blasius_force_report(flow, R: float | None = None, contour=None, rho: float = 1.0, n: int = 256,
                         center=0j) -> dict:
    """:func:`blasius_force` plus a validity flag: dict(D, L, valid, reason). ``valid`` is False when the contour
    enters the body (``flow.inside``). Book: §6.5 (6.60). Label: analytic."""
    if contour is None:
        th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
        cz = as_complex_point(center) + float(R) * np.exp(1j * th)
    else:
        cz = _contour(contour)
    bad = contour_crosses_body(getattr(flow, "inside", None), cz)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        F = blasius_force(flow, R=R, contour=contour, rho=rho, n=n, center=center)
    return {"D": F.D, "L": F.L, "valid": not bad,
            "reason": "contour crosses the body" if bad else "contour encloses the body"}


def laurent_coefficients(flow_or_dwdz, R: float, n: int = 128, kmin: int = -6, kmax: int = 2, center=0j) -> dict:
    """Laurent coefficients c_k of f = dw/dz on the circle |z − c| = R, f = Σ c_k (z − c)^k, by FFT of n samples:
    c_k = F_{k mod n}/(n R^k). For a body with clockwise circulation Γ_cw and doublet d (§6.5):
    c₀ = U, c₋₁ = iΓ_cw/2π (= −iΓ_ccw/2π), c₋₂ = −d/2π (book scalar d).

    Returns dict {k: complex c_k} for kmin ≤ k ≤ kmax. Book: §6.5, far-field expansion before (6.61).
    Validation (planned): V1 cylinder coefficients. Label: analytic.
    """
    f = _dwdz_of(flow_or_dwdz)
    c = as_complex_point(center)
    th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
    vals = _C(f(c + float(R) * np.exp(1j * th)))
    Fk = np.fft.fft(vals)
    return {k: complex(Fk[k % int(n)] / (int(n) * float(R) ** k)) for k in range(int(kmin), int(kmax) + 1)}


# ======================================================================================================================
# checks on fields (§6.2)
# ======================================================================================================================
def laplacian_residual(fn: Callable, x, y, h: float = 1e-3):
    """∇²f at (x, y) by the fourth-order 9-point stencil (5 points in x and 5 in y sharing the centre):
    f_xx ≈ (−f₋₂ + 16f₋₁ − 30f₀ + 16f₁ − f₂)/12h², same in y; error O(h⁴).

    ``fn(x, y)`` a scalar field (ψ or φ) [m²/s]; h [m]. Returns ∇²f [1/s]. Zero for harmonic f ((6.5), (6.12)); NaN
    within 5h of a singular point when ``fn`` is a bound method of an element or Flow (e.g. ``flow.psi``).
    Book: §6.2, Eqs. (6.5), (6.12). Validation (planned): V3 observed order 4. Label: converged.
    """
    x_, y_ = _F(x), _F(y)

    def f(a, b):
        return _F(fn(a, b))

    owner = getattr(fn, "__self__", None)
    sing = ()
    if owner is not None:
        sing = tuple(getattr(owner, "singular_points", ())) if not isinstance(owner, Flow) else tuple(
            q for e in owner.elements for q in e.singular_points)
    c = [-1.0, 16.0, -30.0, 16.0, -1.0]
    s = [-2.0, -1.0, 0.0, 1.0, 2.0]
    fxx = sum(ci * f(x_ + si * h, y_) for ci, si in zip(c, s)) / (12.0 * h ** 2)
    fyy = sum(ci * f(x_, y_ + si * h) for ci, si in zip(c, s)) / (12.0 * h ** 2)
    out = fxx + fyy
    for q in sing:  # the stencil must not straddle a singular point: NaN within 5h of one
        out = np.where(np.hypot(x_ - q.real, y_ - q.imag) < 5.0 * h, np.nan, out)
    return _S(out)


def normal_velocity_values(flow, curve, normals=None, closed: bool = True) -> dict:
    """Normal velocity u·n on a sampled curve (the no-through-flow condition (6.16)).

    ``curve``: complex array or (2, N) vertices, counterclockwise about the body; ``normals`` optional unit normals
    (complex array or (2, N); default: from central differences, n = (t_y, −t_x)/|t|; one-sided at the two ends when
    ``closed=False``, e.g. the half-body outline). Returns dict(max, values (N,), rel = max/|u|max).
    Book: §6.2, Eq. (6.16). Label: analytic.
    """
    z = _contour(curve)
    if normals is None:
        t = np.roll(z, -1) - np.roll(z, 1)
        if not closed:
            t[0], t[-1] = z[1] - z[0], z[-1] - z[-2]
        nn = -1j * t / np.abs(t)
    else:
        nv = np.asarray(normals)
        nn = nv.astype(complex) if np.iscomplexobj(nv) else _F(nv)[0] + 1j * _F(nv)[1]
    u, v = flow.velocity(z.real, z.imag)
    un = _F(u) * nn.real + _F(v) * nn.imag
    spd = np.hypot(_F(u), _F(v))
    mx = float(np.nanmax(np.abs(un)))
    return {"max": mx, "values": un, "rel": mx / max(float(np.nanmax(spd)), 1e-300)}


def normal_velocity_on(flow, curve_pts, normals=None, closed: bool = True) -> float:
    """max |u·n| [m/s] on a sampled curve — the no-through-flow condition (6.16) (≈ 0 on a body; see
    :func:`normal_velocity_values` for the values). ``curve_pts`` complex (or (2, N)), counterclockwise; ``normals``
    complex unit normals or None (from the curve); ``closed=False`` for an open outline (half-body).
    Book: §6.2, Eq. (6.16). Validation (planned): V1 cylinder/half-body/ellipse ≈ 0. Label: analytic."""
    return normal_velocity_values(flow, curve_pts, normals, closed)["max"]


def far_field_check(flow, R: float, n: int = 256, U=None, center=0j) -> float:
    """max |(u − iv) − U| [m/s] on the circle |z − c| = R: the far-field condition (6.17) (→ 0 as R → ∞; ∝ 1/R² for a
    closed body, ∝ 1/R for a net source or a circulation). ``U`` = the target complex velocity u∞ − iv∞ (default: the
    flow's free stream). Points inside an infinite body (half-body) are skipped. Book: §6.2, Eq. (6.17).
    Label: analytic."""
    th = np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
    z = as_complex_point(center) + float(R) * np.exp(1j * th)
    u, v = flow.velocity(z.real, z.imag)
    if U is None:
        U, V = flow.u_inf
    else:
        U, V = complex(U).real, -complex(U).imag
    return float(np.nanmax(np.hypot(_F(u) - U, _F(v) - V)))


def polar_velocity(fn: Callable, r, theta, kind: str = "psi", h: float = 1e-5):
    """Polar velocity components from a polar field f(r, θ) (6.21)–(6.22):
    kind "psi": u_r = (1/r)∂ψ/∂θ, u_θ = −∂ψ/∂r;  kind "phi": u_r = ∂φ/∂r, u_θ = (1/r)∂φ/∂θ (central differences, O(h²)).

    Parameters: fn(r, θ) [m²/s]; r [m]; θ [rad] from +x; h (step in r [m] and in θ [rad]). Returns (u_r, u_θ) [m/s].
    Book: §6.2, Eqs. (6.21)–(6.22). Validation (planned): V1 cylinder (6.34), (6.37). Label: analytic.
    """
    r_, t_ = _F(r), _F(theta)
    dfr = (_F(fn(r_ + h, t_)) - _F(fn(r_ - h, t_))) / (2 * h)
    dft = (_F(fn(r_, t_ + h)) - _F(fn(r_, t_ - h))) / (2 * h)
    if kind == "psi":
        return _S(dft / r_), _S(-dfr)  # Eq. (6.21), (6.22)
    if kind == "phi":
        return _S(dfr), _S(dft / r_)
    raise ValueError('kind must be "psi" or "phi"')


def polar_velocity_sym(expr, r: sp.Symbol, theta: sp.Symbol, kind: str = "psi"):
    """sympy twin of :func:`polar_velocity`: (u_r, u_θ) from ψ(r, θ) or φ(r, θ), Eqs. (6.21)–(6.22). Label: symbolic."""
    if kind == "psi":
        return sp.simplify(sp.diff(expr, theta) / r), sp.simplify(-sp.diff(expr, r))
    if kind == "phi":
        return sp.simplify(sp.diff(expr, r)), sp.simplify(sp.diff(expr, theta) / r)
    raise ValueError('kind must be "psi" or "phi"')


def pressure_coefficient(speed, U: float):
    """C_p = (p − p∞)/(½ρU²) = 1 − |u|²/U² (6.32) for steady ideal flow. speed, U [m/s]. Book: §6.3.
    Label: analytic."""
    return _S(1.0 - _F(speed) ** 2 / float(U) ** 2)  # Eq. (6.32)


# ======================================================================================================================
# axisymmetric and 3-D elements (§6.8)
# ======================================================================================================================
class _AxiElement:
    kind = "axisym"

    def velocity_spherical(self, r, theta):
        R, z = _F(r) * np.sin(theta), _F(r) * np.cos(theta)
        uR, uz = self.velocity_cyl(R, z)
        return _sph_from_cyl(uR, uz, theta)


def _sph_from_cyl(uR, uz, theta):
    th = _F(theta)
    uR, uz = _F(uR), _F(uz)
    return _S(uR * np.sin(th) + uz * np.cos(th)), _S(uR * np.cos(th) - uz * np.sin(th))


@dataclass
class AxisymUniform(_AxiElement):
    """Uniform flow U along z: φ = Uz, ψ = ½UR² (cylindrical) = Ur cos θ, ½Ur² sin²θ (spherical), Eq. (6.86).
    U [m/s]; ψ [m³/s]. Book: §6.8 (6.86). Label: analytic."""
    U: float = 1.0
    kind = "uniform"

    def phi(self, R, z):
        return self.U * _F(z) + 0.0 * _F(R)

    def psi(self, R, z):
        return 0.5 * self.U * _F(R) ** 2 + 0.0 * _F(z)  # Eq. (6.86)

    def velocity_cyl(self, R, z):
        sh = np.broadcast(_F(R), _F(z)).shape
        return np.zeros(sh), np.full(sh, float(self.U))


@dataclass
class PointSource3D(_AxiElement):
    """3-D point source of strength Q [m³/s] on the axis at z0: φ = −Q/4πr, ψ = −Q(z − z0)/4πr = −(Q/4π) cos θ
    (6.87), u = Q e_r/4πr². Book: §6.8 (6.87). Validation (planned): V4 flux Q through any sphere. Label: analytic."""
    Q: float = 1.0
    z0: float = 0.0
    kind = "source"

    def _r(self, R, z):
        return np.sqrt(_F(R) ** 2 + (_F(z) - self.z0) ** 2)

    def phi(self, R, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -self.Q / (_FOUR_PI * self._r(R, z))  # Eq. (6.87)

    def psi(self, R, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -self.Q * (_F(z) - self.z0) / (_FOUR_PI * self._r(R, z))  # Eq. (6.87)

    def velocity_cyl(self, R, z):
        r = self._r(R, z)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = self.Q / (_FOUR_PI * r ** 3)
        return f * _F(R), f * (_F(z) - self.z0)


@dataclass
class Doublet3D(_AxiElement):
    """3-D doublet with dipole strength −d e_z at z0 (6.88): φ = d(z − z0)/4πr³ = (d/4πr²) cos θ,
    ψ = −dR²/4πr³ = −(d/4πr) sin²θ. d [m⁴/s]; the sphere of radius a in a stream U has d = 2πa³U (6.89).
    Book: §6.8 (6.88). Label: analytic."""
    d: float = 1.0
    z0: float = 0.0
    kind = "doublet"

    def _r(self, R, z):
        return np.sqrt(_F(R) ** 2 + (_F(z) - self.z0) ** 2)

    def phi(self, R, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.d * (_F(z) - self.z0) / (_FOUR_PI * self._r(R, z) ** 3)  # Eq. (6.88)

    def psi(self, R, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -self.d * _F(R) ** 2 / (_FOUR_PI * self._r(R, z) ** 3)  # Eq. (6.88)

    def velocity_cyl(self, R, z):
        r = self._r(R, z)
        zz = _F(z) - self.z0
        with np.errstate(divide="ignore", invalid="ignore"):
            uR = -3.0 * self.d * zz * _F(R) / (_FOUR_PI * r ** 5)
            uz = self.d * (_F(R) ** 2 - 2.0 * zz ** 2) / (_FOUR_PI * r ** 5)
        return uR, uz


@dataclass
class LineSource3D(_AxiElement):
    """Uniform line source of density k [m²/s] (k < 0: a line sink) on the axis from z1 to z2:
    ψ = −(k/4π)(r₁ − r₂) with r_i = √(R² + (z − z_i)²) — (6.93)–(6.94) with the sign flipped for a source (the book's
    sink from O to A is k_book = −k: ψ_sink = (k_book/4π)(r − r₁)); φ = −(k/4π)[asinh((z₂ − z)/R) − asinh((z₁ − z)/R)];
    u_z = (k/4π)(1/r₂ − 1/r₁), u_R = (k/4πR)[(z − z₁)/r₁ − (z − z₂)/r₂].
    Book: §6.8 (6.93)–(6.94), the axial singularity method. Validation (planned): V1 closed form = quad of (6.93).
    Label: analytic."""
    k: float = 1.0
    z1: float = 0.0
    z2: float = 1.0
    kind = "line_source"

    def _rr(self, R, z):
        R_, z_ = _F(R), _F(z)
        return np.sqrt(R_ ** 2 + (z_ - self.z1) ** 2), np.sqrt(R_ ** 2 + (z_ - self.z2) ** 2)

    def psi(self, R, z):
        r1, r2 = self._rr(R, z)
        return -self.k / _FOUR_PI * (r1 - r2)  # Eq. (6.94) (source form)

    def phi(self, R, z):
        R_, z_ = np.broadcast_arrays(_F(R), _F(z))
        s1, s2 = self.z1 - z_, self.z2 - z_
        with np.errstate(divide="ignore", invalid="ignore"):
            pos = R_ > 0
            Rs = np.where(pos, R_, 1.0)
            a = np.arcsinh(s2 / Rs) - np.arcsinh(s1 / Rs)
            b = np.sign(s2) * (np.log(np.abs(s2)) - np.log(np.abs(s1)))  # on the axis outside the segment
        return -self.k / _FOUR_PI * np.where(pos, a, b)

    def velocity_cyl(self, R, z):
        R_, z_ = np.broadcast_arrays(_F(R), _F(z))
        r1, r2 = self._rr(R_, z_)
        with np.errstate(divide="ignore", invalid="ignore"):
            uz = self.k / _FOUR_PI * (1.0 / r2 - 1.0 / r1)
            uR = np.where(R_ > 0, self.k / (_FOUR_PI * np.where(R_ > 0, R_, 1.0))
                          * ((z_ - self.z1) / r1 - (z_ - self.z2) / r2), 0.0)
        return uR, uz


class AxisymFlow:
    """Superposition of axisymmetric elements (§6.8): ψ, φ, (u_R, u_z), (u_r, u_θ), speed, C_p.

    Coordinates (R, z) with z along the stream; ψ is the Stokes stream function [m³/s]; ``inside(R, z)`` optional body
    test (NaN there). Book: §6.8, Eqs. (6.74)–(6.95). Label: analytic.
    """

    def __init__(self, elements: Sequence, inside: Callable | None = None, label: str = ""):
        self.elements = list(elements)
        self.inside = inside
        self.label = label

    @property
    def U_inf(self) -> float:
        return float(sum(e.U for e in self.elements if isinstance(e, AxisymUniform)))

    def _mask(self, R, z):
        if self.inside is None:
            return np.zeros(np.broadcast(_F(R), _F(z)).shape, dtype=bool)
        return np.asarray(self.inside(_F(R), _F(z)), dtype=bool)

    def psi(self, R, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            v = sum(_F(e.psi(R, z)) for e in self.elements)
        return _S(np.where(self._mask(R, z), np.nan, v))

    def phi(self, R, z):
        with np.errstate(divide="ignore", invalid="ignore"):
            v = sum(_F(e.phi(R, z)) for e in self.elements)
        return _S(np.where(self._mask(R, z), np.nan, v))

    def velocity_cyl(self, R, z):
        """(u_R, u_z) [m/s], Eq. (6.75)/(6.79); NaN inside the body."""
        uR = 0.0
        uz = 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            for e in self.elements:
                a, b = e.velocity_cyl(R, z)
                uR = uR + _F(a)
                uz = uz + _F(b)
        m = self._mask(R, z)
        sh = m.shape
        return _S(np.where(m, np.nan, np.broadcast_to(uR, sh))), _S(np.where(m, np.nan, np.broadcast_to(uz, sh)))

    def velocity_spherical(self, r, theta):
        """(u_r, u_θ) [m/s] with θ from +z, Eq. (6.83)."""
        R, z = _F(r) * np.sin(theta), _F(r) * np.cos(theta)
        uR, uz = self.velocity_cyl(R, z)
        return _sph_from_cyl(uR, uz, theta)

    def speed(self, R, z):
        uR, uz = self.velocity_cyl(R, z)
        return _S(np.hypot(_F(uR), _F(uz)))

    def cp(self, R, z, U: float | None = None):
        """C_p = 1 − |u|²/U² (6.32)/(6.91)."""
        U_ = self.U_inf if U is None else float(U)
        return _S(1.0 - _F(self.speed(R, z)) ** 2 / U_ ** 2)


def sphere(U: float = 1.0, a: float = 1.0) -> AxisymFlow:
    """Sphere of radius a in a stream U e_z = uniform flow + doublet opposing it, d = 2πa³U (6.89):
    ψ = ½Ur²(1 − a³/r³) sin²θ, φ = Ur(1 + a³/2r³) cos θ; velocity (6.90), surface C_p = 1 − (9/4) sin²θ (6.91).
    Parameters: U [m/s], a [m]. Returns an :class:`AxisymFlow` (NaN inside r < a).
    Book: §6.8 (6.89)–(6.91), Fig. 6.27. Validation (planned): V1 ψ = 0 on r = a and on the axis; parity with
    ``core.vortices.hill_stream_function`` exterior; V1 form cross-check McDonald (2015). Label: analytic.
    """
    U, a = float(U), float(a)

    def inside(R, z):
        return _F(R) ** 2 + _F(z) ** 2 < a ** 2 * (1.0 - 1e-12)

    return AxisymFlow([AxisymUniform(U), Doublet3D(_TWO_PI * a ** 3 * U)], inside=inside, label="sphere")


def axisym_velocity_spherical(fn: Callable, r, theta, kind: str = "psi", h: float = 1e-5):
    """Spherical velocity components from ψ(r, θ) or φ(r, θ), Eq. (6.83):
    u_r = (1/r² sin θ)∂ψ/∂θ = ∂φ/∂r, u_θ = −(1/r sin θ)∂ψ/∂r = (1/r)∂φ/∂θ (central differences).
    fn(r, θ) [m³/s or m²/s]; r [m]; θ [rad] from the axis. Returns (u_r, u_θ) [m/s]. Book: §6.8 (6.83).
    Label: analytic."""
    r_, t_ = _F(r), _F(theta)
    dfr = (_F(fn(r_ + h, t_)) - _F(fn(r_ - h, t_))) / (2 * h)
    dft = (_F(fn(r_, t_ + h)) - _F(fn(r_, t_ - h))) / (2 * h)
    if kind == "psi":
        return _S(dft / (r_ ** 2 * np.sin(t_))), _S(-dfr / (r_ * np.sin(t_)))  # Eq. (6.83)
    if kind == "phi":
        return _S(dfr), _S(dft / r_)
    raise ValueError('kind must be "psi" or "phi"')


def axisym_velocity_spherical_sym(expr, r: sp.Symbol, theta: sp.Symbol, kind: str = "psi"):
    """sympy twin of :func:`axisym_velocity_spherical`, Eq. (6.83). Label: symbolic."""
    if kind == "psi":
        return (sp.simplify(sp.diff(expr, theta) / (r ** 2 * sp.sin(theta))),
                sp.simplify(-sp.diff(expr, r) / (r * sp.sin(theta))))
    if kind == "phi":
        return sp.simplify(sp.diff(expr, r)), sp.simplify(sp.diff(expr, theta) / r)
    raise ValueError('kind must be "psi" or "phi"')


def sphere_potential_vector(x, U_vec, d_vec):
    """Coordinate-free sphere potential (6.92): φ = (U − d/4π|x|³)·x.

    x (3,) or (3, N) [m] from the centre; U_vec (3,) free stream [m/s]; ``d_vec`` the dipole vector (3,) [m⁴/s] (the
    sphere of radius a has d = −2πa³U, i.e. φ = U·x(1 + a³/2|x|³)) — or a plain number, taken as the radius a [m] and
    converted. Returns φ [m²/s]. Book: §6.8 (6.92). Label: analytic."""
    X = _F(x)
    U = _F(U_vec).reshape((3,) + (1,) * (X.ndim - 1))
    r = np.sqrt(np.sum(X ** 2, axis=0))
    if np.ndim(d_vec) == 0:
        d = -_TWO_PI * float(d_vec) ** 3 * U
    else:
        d = _F(d_vec).reshape((3,) + (1,) * (X.ndim - 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(np.sum((U - d / (_FOUR_PI * r ** 3)) * X, axis=0))  # Eq. (6.92)


# ======================================================================================================================
# the arbitrarily moving sphere and added mass (§6.9)
# ======================================================================================================================
def _vec(v, ndim):
    return _F(v).reshape((3,) + (1,) * (ndim - 1))


def moving_sphere_potential(x, xs, us, a: float):
    """Potential of a sphere of radius a centred at x_s moving at u_s through fluid at rest (6.96)–(6.97):
    φ = −(a³/2|ξ|³) u_s·ξ, ξ = x − x_s (dipole d(t) = +2πa³u_s).
    x (3,) or (3, N) [m]; xs, us (3,) [m], [m/s]; a [m]. Returns φ [m²/s].
    Book: §6.9, Eqs. (6.96)–(6.97). Validation (planned): V1 u·n = u_s·n on |ξ| = a. Label: analytic."""
    X = _F(x)
    xi = X - _vec(xs, X.ndim)
    r = np.sqrt(np.sum(xi ** 2, axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(-float(a) ** 3 / (2.0 * r ** 3) * np.sum(_vec(us, X.ndim) * xi, axis=0))  # Eq. (6.97)


def moving_sphere_velocity(x, xs, us, a: float):
    """Fluid velocity ∇φ of the moving sphere (6.103): −(a³/2)[−3ξ(u_s·ξ)/|ξ|⁵ + u_s/|ξ|³].
    x (3,) or (3, N) [m]. Returns (3,) or (3, N) [m/s]. Book: §6.9 (6.103). Label: analytic."""
    X = _F(x)
    xi = X - _vec(xs, X.ndim)
    U = _vec(us, X.ndim)
    r = np.sqrt(np.sum(xi ** 2, axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        return -0.5 * float(a) ** 3 * (-3.0 * xi * np.sum(U * xi, axis=0) / r ** 5 + U / r ** 3)  # Eq. (6.103)


def moving_sphere_surface_velocity(e_xi, us, printed_bracket: bool = False):
    """Fluid velocity on the sphere's surface (6.104): u_a = (3/2)(u_s·e_ξ)e_ξ − ½u_s.

    e_xi (3,) or (3, N) unit outward normals; us (3,) [m/s]. ``printed_bracket=True`` evaluates the book's printed middle
    expression, whose last bracket term −u_s/a³ should be +u_s/a³ (it gives (3/2)(u_s·e)e + ½u_s — a wrong variant kept
    for tests). Book: §6.9 (6.104). Label: analytic."""
    E = _F(e_xi)
    U = _vec(us, E.ndim)
    un = np.sum(U * E, axis=0)
    s = 1.0 if printed_bracket else -1.0
    return 1.5 * un * E + s * 0.5 * U  # Eq. (6.104) (corrected sign unless printed_bracket)


def moving_sphere_dphidt(x, xs, us, dus_dt, a: float):
    """∂φ/∂t at a fixed point (6.101): −u·u_s − (a³/2|ξ|³) ξ·du_s/dt (φ depends on t through x_s and u_s).
    Returns [m²/s²]. Book: §6.9 (6.101). Label: analytic."""
    X = _F(x)
    xi = X - _vec(xs, X.ndim)
    r = np.sqrt(np.sum(xi ** 2, axis=0))
    u = moving_sphere_velocity(X, xs, us, a)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(-np.sum(u * _vec(us, X.ndim), axis=0)
                  - float(a) ** 3 / (2.0 * r ** 3) * np.sum(xi * _vec(dus_dt, X.ndim), axis=0))  # Eq. (6.101)


def moving_sphere_surface_pressure(e_xi, us, dus_dt, a: float, rho: float = 1000.0, p_inf: float = 0.0,
                                   split: bool = False):
    """Surface pressure on an arbitrarily moving sphere (6.105):
    (p_a − p∞)/ρ = ½|u_s|²[(9/4)(u_s·e_ξ)²/|u_s|² − 5/4] + (a/2) e_ξ·du_s/dt
    — a speed (steady) part, fore–aft symmetric (6.106), and an acceleration part (high in front of an accelerating
    sphere, low behind).

    e_xi (3,) or (3, N) unit outward normals; us [m/s], dus_dt [m/s²] (3,); a [m]; ρ [kg/m³]; p∞ [Pa].
    Returns p_a [Pa], or with ``split=True`` dict(steady, acceleration, total) where steady/acceleration are the gauge
    parts p − p∞ [Pa] and total = p_a. Book: §6.9, Eqs. (6.100)–(6.106). Validation (planned): V1 parity with
    ``ch04.accelerating_sphere_pressure``; steady part = (6.91) with θ ↔ π − θ_s. Label: analytic.
    """
    E = _F(e_xi)
    U = _vec(us, E.ndim)
    A = _vec(dus_dt, E.ndim)
    un = np.sum(U * E, axis=0)
    u2 = np.sum(U * U, axis=0)
    steady = rho * (9.0 / 8.0 * un ** 2 - 5.0 / 8.0 * u2)  # ½|u_s|²((9/4)cos² − 5/4), Eq. (6.105)
    accel = rho * 0.5 * float(a) * np.sum(E * A, axis=0)  # (a/2) e_ξ·du_s/dt, Eq. (6.105)
    total = p_inf + steady + accel
    if split:
        return {"steady": _S(steady), "acceleration": _S(accel), "total": _S(total)}
    return _S(total)


def sphere_force_quadrature(p_fn: Callable, a: float, n_theta: int = 32, n_phi: int = 64) -> np.ndarray:
    """Pressure force on a sphere F_s = −∮(p − p∞) n dA (6.98) by Gauss–Legendre in cos θ × periodic trapezoid in φ.

    ``p_fn(e)`` returns the surface pressure (or p − p∞) [Pa] at unit normals e (3, N); a [m]. Returns F (3,) [N].
    Exact (to round-off) for pressures that are polynomials of degree < 2 n_theta in e. Book: §6.9, Eqs. (6.98),
    (6.107)–(6.108). Label: analytic."""
    mu, wmu = np.polynomial.legendre.leggauss(int(n_theta))  # μ = cos θ ∈ (−1, 1)
    ph = np.linspace(0.0, _TWO_PI, int(n_phi), endpoint=False)
    MU, PH = np.meshgrid(mu, ph, indexing="ij")
    S = np.sqrt(1.0 - MU ** 2)
    E = np.stack([S * np.cos(PH), S * np.sin(PH), MU]).reshape(3, -1)
    W = (wmu[:, None] * np.full(PH.shape, _TWO_PI / int(n_phi))).ravel() * float(a) ** 2  # dA = a² dμ dφ
    p = _F(p_fn(E)).ravel()
    return -np.sum(p[None, :] * E * W[None, :], axis=1)  # Eq. (6.98)


def added_mass_sphere(a: float, rho: float = 1000.0) -> float:
    """Added (apparent) mass of a sphere M = 2πρa³/3 = ½ × displaced fluid mass [kg] (6.108); F_s = −M du_s/dt.
    Book: §6.9, Eqs. (6.108)–(6.109). Validation (planned): V1 force quadrature; independent kinetic-energy route;
    V5 Wikipedia "Added mass" coefficient ½. Label: analytic."""
    return float(_TWO_PI * float(rho) * float(a) ** 3 / 3.0)  # Eq. (6.108)
