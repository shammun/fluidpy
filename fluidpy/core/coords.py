"""Curvilinear coordinates used from Ch. 3 on: plane polar (r, θ), cylindrical (R, φ, z) and spherical (r, θ, φ) —
conversions, local unit vectors and velocity components as projections on them.

Book: Ch. 3 §3.1, Fig. 3.3 (a)–(d). Definitions as the book uses them (Exercises 3.1–3.2 state the relations, used here
as definitions — the exercise text is not reproduced): R = √(x² + y²), φ = tan⁻¹(y/x) (cylindrical, capital R);
r = √(x² + y² + z²), θ = tan⁻¹(√(x² + y²)/z) the **polar angle from +z**, φ the azimuth (spherical, physics
convention); plane polar (r, θ) with θ from the x-axis (§3.5). Angles in radians; ``np.arctan2`` everywhere (quadrants).
Velocity component names follow the book: (u, v, w), (u_R, u_φ, u_z), (u_r, u_θ, u_φ).

Unit vectors are returned as **rows**: ``E[k]`` is the k-th unit vector in Cartesian components, so the component of
a Cartesian vector u along it is ``E[k] @ u`` and E is the transpose of the direction-cosine matrix C of (2.5)
(C's columns are the new unit vectors). On the axis (R = 0) φ = arctan2(0, 0) = 0 is used (a valid basis).
"""
from __future__ import annotations

import numpy as np

from ._util import as_scalar_if_0d

__all__ = ["polar_from_cartesian", "cartesian_from_polar_coords", "cylindrical_from_cartesian",
           "cartesian_from_cylindrical", "spherical_from_cartesian", "cartesian_from_spherical",
           "unit_vectors_polar", "unit_vectors_cylindrical", "unit_vectors_spherical", "velocity_components",
           "cartesian_components"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def polar_from_cartesian(x, y):
    """Plane polar coordinates (r, θ) of (x, y): r = √(x² + y²), θ = atan2(y, x) ∈ (−π, π] (Fig. 3.3a). Units m, rad.
    Label: analytic.
    Book: §3.1, Fig. 3.3a (plane polar coordinates).
    """
    x_, y_ = _F(x), _F(y)
    return _S(np.hypot(x_, y_)), _S(np.arctan2(y_, x_))


def cartesian_from_polar_coords(r, theta):
    """(x, y) = (r cos θ, r sin θ) (Fig. 3.3a). Label: analytic.
    Book: §3.1, Fig. 3.3a.
    """
    r_, t_ = _F(r), _F(theta)
    return _S(r_ * np.cos(t_)), _S(r_ * np.sin(t_))


def cylindrical_from_cartesian(x, y, z):
    """Cylindrical coordinates (R, φ, z): R = √(x² + y²), φ = tan⁻¹(y/x) (``arctan2``, (−π, π]), z = z.

    Book: §3.1, Fig. 3.3c. Units: m, rad, m. Validation: V1 round trip on 10⁴ random points (1e-13). Label: analytic.
    """
    x_, y_, z_ = _F(x), _F(y), _F(z)
    return _S(np.hypot(x_, y_)), _S(np.arctan2(y_, x_)), _S(z_ + 0.0)


def cartesian_from_cylindrical(R, phi, z):
    """(x, y, z) = (R cos φ, R sin φ, z) (Fig. 3.3c). Label: analytic.
    Book: §3.1, Fig. 3.3c.
    """
    R_, p_, z_ = _F(R), _F(phi), _F(z)
    return _S(R_ * np.cos(p_)), _S(R_ * np.sin(p_)), _S(z_ + 0.0)


def spherical_from_cartesian(x, y, z):
    """Spherical coordinates (r, θ, φ): r = √(x² + y² + z²), θ = tan⁻¹(√(x² + y²)/z) ∈ [0, π] (polar angle from +z,
    via ``arctan2(√(x² + y²), z)``), φ = tan⁻¹(y/x) ∈ (−π, π] (azimuth).

    Book: §3.1, Fig. 3.3d (physics convention: θ from +z). Units: m, rad, rad.
    Validation: V1 round trip; θ = 0 on +z, π on −z, π/2 in the x–y plane. Label: analytic.
    """
    x_, y_, z_ = _F(x), _F(y), _F(z)
    Rh = np.hypot(x_, y_)
    return _S(np.hypot(Rh, z_)), _S(np.arctan2(Rh, z_)), _S(np.arctan2(y_, x_))


def cartesian_from_spherical(r, theta, phi):
    """(x, y, z) = (r sin θ cos φ, r sin θ sin φ, r cos θ), θ polar from +z (Fig. 3.3d). Label: analytic.
    Book: §3.1, Fig. 3.3d.
    """
    r_, t_, p_ = _F(r), _F(theta), _F(phi)
    return _S(r_ * np.sin(t_) * np.cos(p_)), _S(r_ * np.sin(t_) * np.sin(p_)), _S(r_ * np.cos(t_))


def _stack_rows(rows) -> np.ndarray:
    """Stack row vectors whose components may be arrays: result (k, d) or (k, d, N)."""
    return np.array([np.broadcast_arrays(*[_F(c) for c in row]) for row in rows], dtype=float)


def unit_vectors_polar(theta) -> np.ndarray:
    """Rows e_r = (cos θ, sin θ), e_θ = (−sin θ, cos θ) of plane polar coordinates (Fig. 3.3a; ch02 Ex. 2.1).
    Shape (2, 2) or (2, 2, N). Label: analytic.
    Book: §3.1, Fig. 3.3a; ch02 Ex. 2.1.
    """
    t = _F(theta)
    c, s = np.cos(t), np.sin(t)
    return _stack_rows([[c, s], [-s, c]])


def unit_vectors_cylindrical(phi) -> np.ndarray:
    """Rows e_R = (cos φ, sin φ, 0), e_φ = (−sin φ, cos φ, 0), e_z = (0, 0, 1) (Fig. 3.3c).

    Shape (3, 3) (or (3, 3, N) for array φ). Right-handed: e_R × e_φ = e_z. Validation: V1 orthonormal, det +1.
    Label: analytic.
    Book: §3.1, Fig. 3.3c.
    """
    p = _F(phi)
    c, s, z, o = np.cos(p), np.sin(p), np.zeros_like(p), np.ones_like(p)
    return _stack_rows([[c, s, z], [-s, c, z], [z, z, o]])


def unit_vectors_spherical(theta, phi) -> np.ndarray:
    """Rows e_r, e_θ, e_φ of spherical coordinates (θ polar from +z, φ azimuth; Fig. 3.3d):
    e_r = (sin θ cos φ, sin θ sin φ, cos θ), e_θ = (cos θ cos φ, cos θ sin φ, −sin θ), e_φ = (−sin φ, cos φ, 0).

    Right-handed: e_r × e_θ = e_φ. Shape (3, 3) or (3, 3, N). Validation: V1 orthonormal, det +1. Label: analytic.
    Book: §3.1, Fig. 3.3d.
    """
    t, p = np.broadcast_arrays(_F(theta), _F(phi))
    st, ct, sp_, cp = np.sin(t), np.cos(t), np.sin(p), np.cos(p)
    return _stack_rows([[st * cp, st * sp_, ct], [ct * cp, ct * sp_, -st], [-sp_, cp, np.zeros_like(t)]])


def _basis_at(x, system: str) -> np.ndarray:
    x_ = _F(x)
    if system == "polar":
        return unit_vectors_polar(np.arctan2(x_[1], x_[0]))
    if system == "cylindrical":
        return unit_vectors_cylindrical(np.arctan2(x_[1], x_[0]))
    if system == "spherical":
        _, th, ph = spherical_from_cartesian(x_[0], x_[1], x_[2])
        return unit_vectors_spherical(th, ph)
    raise ValueError('system must be "polar", "cylindrical" or "spherical"')


def velocity_components(u_cart, x, system: str = "cylindrical") -> np.ndarray:
    """Components of a Cartesian vector u at position x along the local unit vectors of a curvilinear system.

    Book: §3.1, Fig. 3.3 — (u_r, u_θ) plane polar, (u_R, u_φ, u_z) cylindrical, (u_r, u_θ, u_φ) spherical;
    component = projection on the unit vector (the passive rule x′ = Cᵀx of (2.5), ch02 C03).

    Parameters
    ----------
    u_cart : (d,) or (d, N) — Cartesian components [m/s]
    x : (d,) or (d, N) — position(s) [m] where the basis is taken
    system : "polar" (d = 2), "cylindrical" or "spherical" (d = 3)

    Returns
    -------
    ndarray like ``u_cart`` — curvilinear components.

    Validation: V1 |u| preserved; solid body u = Ω e_z × x has u_φ = ΩR, u_R = u_z = 0. Label: analytic.
    """
    E = _basis_at(x, system)
    return np.einsum("kj...,j...->k...", E, _F(u_cart))


def cartesian_components(u_curv, x, system: str = "cylindrical") -> np.ndarray:
    """Inverse of :func:`velocity_components`: Cartesian components from curvilinear ones at position x
    (u = Σ_k u_k e_k). Label: analytic.
    Book: §3.1, Fig. 3.3 (inverse projection, (2.7)).
    """
    E = _basis_at(x, system)
    return np.einsum("kj...,k...->j...", E, _F(u_curv))
