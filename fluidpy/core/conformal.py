"""Conformal maps and the flows they carry (§6.6): local rotation and stretching of small elements (6.63)–(6.64), grid
images, and the Zhukhovsky (Joukowski) transformation z = ζ + b²/ζ (6.65) with the branch-safe inverse (6.69).

Book: Kundu, Cohen & Dowling 5e, §6.6, Eqs. (6.63)–(6.69), Figs. 6.19–6.21 (transcribed from the page images p249–p252).
Reused by Ch. 14 (Zhukhovsky / Kármán–Trefftz airfoils, the Kutta condition).

⚠️ Branch trap (analysis §9): numpy's principal ``sqrt(z**2 - 4*b**2)`` in (6.69) returns the root *inside* the circle
for Re z < 0 (at z = −3 + 0.5i with b = 1 it gives |ζ| ≈ 0.37, not ≈ 2.70). The correct branch ζ ~ z at infinity is
ζ = ½[z + √(z − 2b)·√(z + 2b)] (product of principal roots: the cut is the slit −2b ≤ x ≤ 2b). Both are callable
(``branch="outside"`` / ``"principal"``) so the wrong variant can be shown and tested.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from ._util import as_scalar_if_0d
from .potential import FunctionFlow, gamma_ccw_from

__all__ = ["conformal_map", "map_elements", "map_elements_info", "angle_preservation", "grid_image", "map_grid_lines",
           "joukowski", "joukowski_derivative", "joukowski_inverse", "joukowski_inverse_derivative", "mapped_flow",
           "circle_flow_zeta"]

_S = as_scalar_if_0d


def _C(a):
    return np.asarray(a, dtype=complex)


def conformal_map(name: str, b: float = 1.0) -> tuple[Callable, Callable]:
    """Named analytic maps (f, df/dz) for the explainers and notebook: "square" z², "exp" e^z, "sin" sin z, "log" ln z,
    "joukowski" z + b²/z, "identity". Book: §6.6 (w = z², w = ln(sin z), (6.65)). Label: analytic."""
    maps = {
        "identity": (lambda z: _C(z), lambda z: np.ones_like(_C(z))),
        "square": (lambda z: _C(z) ** 2, lambda z: 2.0 * _C(z)),
        "exp": (lambda z: np.exp(_C(z)), lambda z: np.exp(_C(z))),
        "sin": (lambda z: np.sin(_C(z)), lambda z: np.cos(_C(z))),
        "log": (lambda z: np.log(_C(z)), lambda z: 1.0 / _C(z)),
        "joukowski": (lambda z: joukowski(z, b), lambda z: joukowski_derivative(z, b)),
    }
    if name not in maps:
        raise ValueError(f"unknown map {name!r}; choose one of {sorted(maps)}")
    return maps[name]


def _fdf(f, dfdz, b):
    if isinstance(f, str):
        return conformal_map(f, b)
    return f, dfdz


def map_elements(f, dfdz=None, z0: complex = 1.0, dz_list=(1e-3, 1e-3j), b: float = 1.0) -> np.ndarray:
    """Images of small elements δz at z0 under an analytic map: δw = (dw/dz) δz (6.63), (6.64) — returns the complex
    array f′(z0)·δz (see :func:`map_elements_info` for the exact images, scale and turn). ``f`` a callable (with
    ``dfdz``) or a :func:`conformal_map` name. Book: §6.6 (6.63)–(6.64). Label: analytic."""
    return np.array(map_elements_info(f, dfdz, z0, dz_list, b)["dw_linear"], dtype=complex)


def map_elements_info(f, dfdz=None, z0: complex = 1.0, dz_list=(1e-3, 1e-3j), b: float = 1.0) -> dict:
    """dict(dw_linear (list of f′(z0)δz), dw_exact (f(z0 + δz) − f(z0)), scale |f′(z0)|, turn arg f′(z0) [rad],
    fprime) for the elements of :func:`map_elements`. Book: §6.6 (6.63)–(6.64). Label: analytic."""
    f, dfdz = _fdf(f, dfdz, b)
    z0 = complex(z0)
    fp = complex(_C(dfdz(np.array([z0])))[0])
    f0 = complex(_C(f(np.array([z0])))[0])
    lin = [fp * complex(d) for d in dz_list]
    ex = [complex(_C(f(np.array([z0 + complex(d)])))[0]) - f0 for d in dz_list]
    return {"dw_linear": lin, "dw_exact": ex, "scale": abs(fp), "turn": float(np.angle(fp)), "fprime": fp}


def angle_preservation(f, dfdz=None, z0: complex = 1.0, dz1: complex = 1.0, dz2: complex = 1j, eps: float = 1e-5,
                       b: float = 1.0) -> tuple[float, float]:
    """Angle α between two small elements δz, δ′z at z0 and the angle β between their images δw, δ′w (6.63)–(6.64):
    α = arg(δ′z/δz), β = arg(δ′w/δw) with the elements scaled to length ``eps`` and mapped exactly by f — as centred
    images ½[f(z0 + δ) − f(z0 − δ)] (error O(eps²)), or one-sided f(z0 + δ) − f(z0) at a critical point, where the
    centred image vanishes.

    α = β wherever f is analytic with f′(z0) ≠ 0; at a critical point (f′ = 0, e.g. z² at 0) the angle is multiplied.
    Returns (alpha, beta) [rad]. Book: §6.6, Fig. 6.19. Validation (planned): V1 z², e^z, sin z, Zhukhovsky; V7 z² at 0
    doubles the angle. Label: analytic."""
    f, _ = _fdf(f, dfdz, b)
    z0 = complex(z0)
    d1 = eps * complex(dz1) / abs(complex(dz1))
    d2 = eps * complex(dz2) / abs(complex(dz2))

    def F(q):
        return complex(_C(f(np.array([q])))[0])

    f0 = F(z0)
    imgs = []
    for d in (d1, d2):
        one = F(z0 + d) - f0
        cen = 0.5 * (F(z0 + d) - F(z0 - d))
        imgs.append(cen if abs(cen) > 1e-3 * abs(one) else one)
    return float(np.angle(d2 / d1)), float(np.angle(imgs[1] / imgs[0]))


def grid_image(w_fn: Callable, xlim=(-2.0, 2.0), ylim=(-2.0, 2.0), n: int = 200, levels: int = 24) -> dict:
    """φ = Re w and ψ = Im w of a complex potential on an (n × n) grid of the z-plane — contour them to see the images
    of the rectangular (φ, ψ) net of the w-plane (Fig. 6.20). Grids use the project layout [j, i] = (y, x).
    Returns dict(X, Y, phi, psi, phi_levels, psi_levels) (``levels`` equally spaced values between the 3rd and 97th
    percentiles, so a singularity does not swamp them; non-finite values masked as NaN).
    Book: §6.6, Fig. 6.20 (w = z²: ψ = 2xy hyperbolae). Label: analytic."""
    x = np.linspace(xlim[0], xlim[1], int(n))
    y = np.linspace(ylim[0], ylim[1], int(n))
    X, Y = np.meshgrid(x, y, indexing="xy")
    with np.errstate(all="ignore"):
        W = _C(w_fn(X + 1j * Y))
    ph = np.where(np.isfinite(W), W.real, np.nan)
    ps = np.where(np.isfinite(W), W.imag, np.nan)

    def lv(A):
        f = A[np.isfinite(A)]
        lo, hi = np.percentile(f, [3, 97]) if f.size else (0.0, 1.0)
        return np.linspace(lo, hi, int(levels))

    return {"X": X, "Y": Y, "phi": ph, "psi": ps, "phi_levels": lv(ph), "psi_levels": lv(ps)}


def map_grid_lines(g: Callable, u_levels, v_levels, u_range, v_range, n: int = 200) -> dict:
    """Images under z = g(ζ) of the straight lines Re ζ = const (``u_levels``) and Im ζ = const (``v_levels``) —
    e.g. the rectangular (φ, ψ) net pushed through z = f⁻¹(w). Returns dict(u_lines, v_lines) of complex arrays.
    Book: §6.6, Fig. 6.20. Label: analytic."""
    s = np.linspace(u_range[0], u_range[1], int(n))
    t = np.linspace(v_range[0], v_range[1], int(n))
    with np.errstate(all="ignore"):
        ul = [_C(g(u + 1j * t)) for u in np.atleast_1d(u_levels)]
        vl = [_C(g(s + 1j * v)) for v in np.atleast_1d(v_levels)]
    return {"u_lines": ul, "v_lines": vl}


# ======================================================================================================================
# Zhukhovsky transformation (6.65)–(6.69)
# ======================================================================================================================
def joukowski(zeta, b: float = 1.0):
    """Zhukhovsky transformation z = ζ + b²/ζ (6.65): identity far away; circle |ζ| = b → slit −2b ≤ x ≤ 2b;
    circle |ζ| = a > b → ellipse with semi-axes a + b²/a, a − b²/a (6.66)–(6.67). ζ, b [m]. Book: §6.6.
    Label: analytic."""
    z = _C(zeta)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(z + float(b) ** 2 / z)  # Eq. (6.65)


def joukowski_derivative(zeta, b: float = 1.0):
    """dz/dζ = 1 − b²/ζ² (zero at the critical points ζ = ±b, where angles are not preserved). Label: analytic."""
    z = _C(zeta)
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(1.0 - float(b) ** 2 / z ** 2)


def joukowski_inverse(z, b: float = 1.0, branch: str = "outside"):
    """Inverse Zhukhovsky map (6.69) ζ = ½z + ½(z² − 4b²)^{1/2}, the root **outside** the circle |ζ| = b.

    ``branch="outside"`` (correct): ζ = ½[z + √(z − 2b)·√(z + 2b)] — ζ ~ z at infinity, |ζ| ≥ b everywhere, the cut on
    the slit. ``branch="principal"`` (⚠️ wrong variant, for E6's toggle and tests): ½[z + sqrt(z² − 4b²)] with numpy's
    principal root, which lands inside the circle for Re z < 0.
    z [m] (complex, any shape), b [m] > 0. Returns ζ. Book: §6.6 (6.69) ("the negative root, which falls inside the
    cylinder, has been excluded"). Validation (planned): V1 joukowski_inverse(joukowski(ζ)) = ζ for |ζ| > b in all four
    quadrants; wrong-variant test. Label: analytic.
    """
    zz = _C(z)
    b = float(b)
    if branch == "outside":
        return _S(0.5 * (zz + np.sqrt(zz - 2.0 * b) * np.sqrt(zz + 2.0 * b)))  # Eq. (6.69), outside root
    if branch == "principal":
        return _S(0.5 * (zz + np.sqrt(zz ** 2 - 4.0 * b ** 2)))  # printed form with numpy's root (wrong for Re z < 0)
    raise ValueError('branch must be "outside" or "principal"')


def joukowski_inverse_derivative(z, b: float = 1.0, branch: str = "outside"):
    """dζ/dz = 1/(dz/dζ) = 1/(1 − b²/ζ²) at ζ = joukowski_inverse(z). Label: analytic."""
    zeta = _C(joukowski_inverse(z, b, branch))
    with np.errstate(divide="ignore", invalid="ignore"):
        return _S(1.0 / (1.0 - float(b) ** 2 / zeta ** 2))


def circle_flow_zeta(U: float, a: float, *, Gamma_cw: float | None = None, Gamma_ccw: float | None = None,
                     alpha: float = 0.0) -> tuple[Callable, Callable]:
    """Complex potential of the ζ-plane flow round a circle of radius a with clockwise circulation Γ (6.68), optionally
    with the stream at angle α to the real axis: w(ζ) = U(ζe^{−iα} + a²e^{iα}/ζ) + (iΓ_cw/2π) ln(ζ/a).
    Returns (w(ζ), dw/dζ). Book: §6.6 (6.68) (α = 0). Label: analytic."""
    U, a, al = float(U), float(a), float(alpha)
    G = gamma_ccw_from(Gamma_cw, Gamma_ccw)
    Gcw = -G
    em, ep = np.exp(-1j * al), np.exp(1j * al)

    def w(zeta):
        s = _C(zeta)
        with np.errstate(divide="ignore", invalid="ignore"):
            return U * (s * em + a ** 2 * ep / s) + 1j * Gcw / (2 * np.pi) * np.log(s / a)  # Eq. (6.68)

    def dw(zeta):
        s = _C(zeta)
        with np.errstate(divide="ignore", invalid="ignore"):
            return U * (em - a ** 2 * ep / s ** 2) + 1j * Gcw / (2 * np.pi * s)

    return w, dw


def mapped_flow(w_zeta: Callable, dw_dzeta: Callable, z_to_zeta: Callable, dzeta_dz: Callable,
                inside_zeta: Callable | None = None, u_inf=(1.0, 0.0), label: str = "mapped flow") -> FunctionFlow:
    """Carry a ζ-plane flow to the z-plane: w(z) = w_ζ(ζ(z)), u − iv = dw/dz = (dw/dζ)(dζ/dz) (6.69).

    ``inside_zeta(ζ)`` → bool marks the body in the ζ-plane (e.g. |ζ| < a); ``u_inf`` the z-plane free stream (the
    Zhukhovsky map is the identity at infinity). Returns a :class:`~fluidpy.core.potential.FunctionFlow`.
    Book: §6.6, (6.68)–(6.69). Label: analytic."""

    def w(z):
        return _C(w_zeta(z_to_zeta(z)))

    def dw(z):
        return _C(dw_dzeta(z_to_zeta(z))) * _C(dzeta_dz(z))  # Eq. (6.69) chain rule

    inside = None
    if inside_zeta is not None:
        def inside(x, y):
            return np.asarray(inside_zeta(_C(z_to_zeta(np.asarray(x, float) + 1j * np.asarray(y, float)))))

    return FunctionFlow(w, dw, None, inside, u_inf, label)
