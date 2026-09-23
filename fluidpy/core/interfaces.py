"""Boundary and interface conditions: jump conditions from a shrinking pillbox, kinematic conditions on moving and
deforming surfaces, and the force balance on a curved interface with surface tension.

Book: Kundu, Cohen & Dowling 5e, Ch. 4 §4.10, Eqs. (4.90)–(4.98), Figs. 4.18–4.19 (rendered pages
chapters/pages/ch04/p164–p169).

Conventions
-----------
* A moving surface is η(x, t) = 0 with η a callable of points (d,) / (d, N) and time, or a preset name:
  "moving_wall" (V [m/s]: η = x − Vt, fluid velocity (V, 0)) and "linear_wave" (a, k, H, g: η = z − a cos(kx − ωt),
  points (x, z), fluid velocity of the linear wave, ω² = gk tanh kH). n = ∇η/|∇η| points into η > 0.
* Fields follow ``core.kinematics`` (components on axis 0). Stencil steps: h [m], ht [s]; g = 9.81 m/s² by default.
* Book typo handled: §4.10 writes the cap boundary curve C as ζ = x²/2R₁ − y²/2R₂; it is + (the cap z = x²/2R₁ + y²/2R₂
  meets the plane z = ζ); the tangent t and the final integrals in the book use + (analysis §9 item 5).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import dblquad

from ._stencil import ddt, ev, grad
from ._util import as_scalar_if_0d
from .thermo import G_BOOK as G0

__all__ = ["SURFACE_PRESETS", "surface_preset", "pillbox_limit", "surface_normal_speed", "kinematic_bc_residual",
           "relative_normal_velocity", "interface_mass_flux", "cap_pressure_force", "cap_surface_tension_force",
           "laplace_jump_from_balance", "capillary_length"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d

SURFACE_PRESETS = {"moving_wall": dict(V=1.0), "linear_wave": dict(a=0.05, k=1.0, H=np.inf, g=G0)}


def surface_preset(name: str, **p):
    """(eta(x, t), u(x, t)) of a preset moving surface (see the module doc); our test surfaces.
    Book: §4.10, Eqs. (4.90)–(4.91). Label: analytic."""
    if name not in SURFACE_PRESETS:
        raise ValueError(f"unknown surface {name!r}; choose from {tuple(SURFACE_PRESETS)}")
    q = dict(SURFACE_PRESETS[name])
    q.update({k: v for k, v in p.items() if k in q})
    if name == "moving_wall":
        V = float(q["V"])

        def eta(X, t):
            return _F(X)[0] - V * t

        def u(X, t):
            X_ = _F(X)
            out = np.zeros_like(X_)
            out[0] = V
            return out
        return eta, u
    a, k, H, g = float(q["a"]), float(q["k"]), float(q["H"]), float(q["g"])
    om = np.sqrt(g * k * (np.tanh(k * H) if np.isfinite(H) else 1.0))

    def eta_w(X, t):
        X_ = _F(X)
        return X_[1] - a * np.cos(k * X_[0] - om * t)

    def u_w(X, t):
        X_ = _F(X)
        ph = k * X_[0] - om * t
        if np.isfinite(H):
            fu = np.cosh(k * (X_[1] + H)) / np.sinh(k * H)
            fw = np.sinh(k * (X_[1] + H)) / np.sinh(k * H)
        else:
            fu = fw = np.exp(k * X_[1])
        return np.stack([a * om * fu * np.cos(ph), a * om * fw * np.sin(ph)])
    return eta_w, u_w


def _surface(eta, u, p):
    if isinstance(eta, str):
        e, uu = surface_preset(eta, **p)
        return e, (uu if u is None else u)
    return eta, u


def pillbox_limit(flux_top, flux_bottom, side_rate, volume_rate, l) -> dict:
    """The shrinking-pillbox balance of Fig. 4.18 for a cylinder of height l straddling an interface (per unit end area).

    Book: §4.10, Fig. 4.18: (F·n)₁ − (F·n)₂ + (side flux) + (volume storage/sources) = 0; the side and volume terms are
    ∝ l, so as l → 0 the normal flux is continuous: ρ₁u₁·n = ρ₂u₂·n; n_iτ_ij continuous (no surface tension); n_iq_i
    continuous (k₁∂T₁/∂n = k₂∂T₂/∂n).

    Parameters: flux_top, flux_bottom : normal fluxes on the end faces in media 1 and 2 [per unit area];
    side_rate : side outflow per unit end area per unit height [per m]; volume_rate : storage/source per unit volume;
    l : pillbox height(s) [m].
    Returns dict(jump = flux_top − flux_bottom, side = side_rate·l, volume = volume_rate·l, residual = jump + side +
    volume) — for a balanced pillbox residual = 0, so jump = −(side + volume) → 0 linearly as l → 0.
    Validation: V3 side and volume ∝ l (slope 1). Label: analytic.
    """
    l_ = _F(l)
    jump = float(flux_top) - float(flux_bottom)
    side, vol = float(side_rate) * l_, float(volume_rate) * l_
    return {"jump": jump, "side": _S(side), "volume": _S(vol), "residual": _S(jump + side + vol)}


def surface_normal_speed(eta, x, t: float = 0.0, h: float = 1e-4, ht: float | None = None, **p):
    """Normal speed of a moving surface η(x, t) = 0: u_s·n = −(∂η/∂t)/|∇η| [m/s] (from (4.90)).

    Book: §4.10, Eq. (4.90): dη/dt = ∂η/∂t + (u_s·∇)η = 0 on η = 0 ⇒ u_s·∇η = −∂η/∂t; only the normal part of u_s is
    defined. Parameters: eta (callable or preset + **p); x surface point(s) [m]; t [s]; h [m]; ht [s].
    Validation: V1 plane wall η = x − Vt: normal speed V. Label: analytic.
    """
    e, _ = _surface(eta, None, p)
    x_ = _F(x)
    gE = grad(e, x_, t, h)
    return _S(-ddt(e, x_, t, ht) / np.sqrt(np.sum(gE ** 2, axis=0)))


def kinematic_bc_residual(eta, u=None, x=None, t: float = 0.0, h: float = 1e-4, ht: float | None = None, **p):
    """Kinematic condition of a surface with no mass flux, Dη/Dt = ∂η/∂t + (u·∇)η on η = 0, Eq. (4.91) [η units/s].

    Book: §4.10, Eq. (4.91) (fluid particles on the surface stay on it; u·∇η = u_s·∇η = −∂η/∂t).
    Parameters: eta (callable or preset + **p); u(x, t) [m/s] (None → the preset's velocity); x surface point(s) [m];
    t [s]; h [m]; ht [s]. Returns Dη/Dt — zero for a material surface.
    Validation: V1 rigid wall η = x − Vt with u·n = V: 0; linear wave: the full condition's residual is O((ka)²) (slope 2
    in ka) while the linearised one vanishes; V3 stencil order 2. Label: analytic, converged.
    """
    e, uf = _surface(eta, u, p)
    x_ = _F(x)
    d = x_.shape[0]
    return _S(ddt(e, x_, t, ht) + np.sum(ev(uf, x_, t, (d,)) * grad(e, x_, t, h), axis=0))  # Eq. (4.91)


def relative_normal_velocity(eta, u=None, x=None, t: float = 0.0, h: float = 1e-4, ht: float | None = None, **p):
    """Fluid velocity relative to a moving surface, normal component: (u_rel)_n = u·n − u_s·n = (1/|∇η|)Dη/Dt, (4.92).

    Book: §4.10, Eq. (4.92), n = ∇η/|∇η|. Returns [m/s] (0 for no mass flux). Label: analytic.
    """
    e, uf = _surface(eta, u, p)
    x_ = _F(x)
    gE = grad(e, x_, t, h)
    return _S(_F(kinematic_bc_residual(e, uf, x_, t, h, ht)) / np.sqrt(np.sum(gE ** 2, axis=0)))  # Eq. (4.92)


def interface_mass_flux(eta, u=None, rho=1000.0, x=None, t: float = 0.0, h: float = 1e-4, ht: float | None = None,
                        **p):
    """Mass flux across a moving surface per unit area, (ρ/|∇η|)Dη/Dt [kg/(m² s)], Eq. (4.93).

    Book: §4.10, Eq. (4.93) (e.g. a moving shock (Ch. 15) or an evaporating interface). ``rho`` a constant or a field.
    Validation: V1 a plane front η = x − Vt in a uniform flow U: ρ(U − V). Label: analytic.
    """
    e, uf = _surface(eta, u, p)
    x_ = _F(x)
    r = _F(rho(x_, t)) if callable(rho) else _F(rho)
    return _S(r * _F(relative_normal_velocity(e, uf, x_, t, h, ht)))  # Eq. (4.93)


# ======================================================================================================================
# surface tension on a curved cap (Fig. 4.19)
# ======================================================================================================================
def cap_pressure_force(dp: float, R1: float, R2: float, zeta: float, exact: bool = True,
                       components: bool = False):
    """Net pressure force on the cap z = x²/2R₁ + y²/2R₂ below the plane z = ζ, with p higher by Δp above, Eq. (4.97).

    Book: §4.10, Eq. (4.97): F_p = −∬_A Δp n dA = −Δp ∫∫ (−x/R₁, −y/R₂, 1) dy dx over the ellipse x²/2R₁ + y²/2R₂ ≤ ζ;
    (F_p)_z = −πΔp√(2R₁ζ)√(2R₂ζ); the x, y components vanish by symmetry. (n dA = (−x/R₁, −y/R₂, 1) dx dy exactly on the
    paraboloid, so the closed form is exact for this cap.)

    Parameters: dp [Pa]; R1, R2 > 0 principal radii [m]; zeta : cap height [m]; exact : True → ``dblquad`` over the cap,
    False → closed form; components : return (F_x, F_y, F_z) instead of F_z.
    Returns F_z [N] (or the 3-vector). Validation: V1 dblquad = closed form to 1e-10 relative. Label: analytic.
    """
    a, b = np.sqrt(2 * R1 * zeta), np.sqrt(2 * R2 * zeta)
    if not exact:
        F = np.array([0.0, 0.0, -np.pi * dp * a * b])  # (F_p)_z = −πΔp√(2R₁ζ)√(2R₂ζ)
    else:
        lo = lambda x: -np.sqrt(max(2 * R2 * zeta - x ** 2 * R2 / R1, 0.0))  # noqa: E731
        hi = lambda x: np.sqrt(max(2 * R2 * zeta - x ** 2 * R2 / R1, 0.0))  # noqa: E731
        import warnings
        comps = []
        for f in (lambda y, x: -x / R1, lambda y, x: -y / R2, lambda y, x: 1.0):
            with warnings.catch_warnings():  # the x, y integrands are odd: their zero value defeats epsrel
                warnings.simplefilter("ignore")
                val, _ = dblquad(f, -a, a, lo, hi, epsabs=1e-13 * a * b, epsrel=1e-12)
            comps.append(-dp * val)  # Eq. (4.97)
        F = np.array(comps)
    return F if components else float(F[2])


def cap_surface_tension_force(sigma: float, R1: float, R2: float, zeta: float, exact: bool = True, n: int = 4096,
                              components: bool = False):
    """Net surface-tension force σ∮_C t × n ds on the rim of the cap z = x²/2R₁ + y²/2R₂ at height ζ, Eq. (4.98).

    Book: §4.10, Eq. (4.98) and the small-ζ evaluation (F_st)_z = πσ√(2R₁ζ)√(2R₂ζ)(1/R₁ + 1/R₂) (the book's quarter-path
    integral after the substitution sin ξ = x/√(2R₁ζ)).

    Parameters
    ----------
    sigma [N/m]; R1, R2 [m]; zeta [m]
    exact : True → integrate t × n round the exact rim ellipse x = a cos θ, y = b sin θ (a = √(2R₁ζ), b = √(2R₂ζ)) with
        the exact paraboloid normal (midpoint rule in θ, spectrally accurate; no singular R₂/y factor); False → the
        book's small-ζ closed form
    n : θ nodes;  components : return the 3-vector (x, y components ≈ 0 by symmetry)

    Returns F_z [N] (t counter-clockwise seen from +z, so t × n points outward and up).
    Validation: V1 exact → closed form as ζ → 0 with (ratio − 1) = O(ζ/R) (slope 1); V2 sympy of the book's quarter
    integral. Label: analytic, converged.
    """
    a, b = np.sqrt(2 * R1 * zeta), np.sqrt(2 * R2 * zeta)
    if not exact:
        F = np.array([0.0, 0.0, np.pi * sigma * a * b * (1.0 / R1 + 1.0 / R2)])
    else:
        th = (np.arange(n) + 0.5) * 2 * np.pi / n
        x, y = a * np.cos(th), b * np.sin(th)
        dX = np.stack([-a * np.sin(th), b * np.cos(th), np.zeros_like(th)])  # dX/dθ on the rim (z = ζ constant)
        ds = np.sqrt(np.sum(dX ** 2, axis=0))
        tvec = dX / ds
        nn = np.stack([-x / R1, -y / R2, np.ones_like(th)])
        nvec = nn / np.sqrt(np.sum(nn ** 2, axis=0))
        F = sigma * np.sum(np.cross(tvec, nvec, axis=0) * ds, axis=1) * (2 * np.pi / n)  # Eq. (4.98)
    return F if components else float(F[2])


def laplace_jump_from_balance(sigma: float, R1: float, R2: float, zeta: float = 1e-6, exact: bool = False) -> float:
    """Pressure jump Δp that makes the pressure force (4.97) and the surface-tension force (4.98) on the cap cancel —
    → σ(1/R₁ + 1/R₂) (Laplace, (1.5)) as ζ → 0.

    Book: §4.10, the balance F_p + F_st = 0 after (4.98), recovering (1.5) (ch01 ``laplace_pressure_jump``); the higher
    pressure is on the side of the centres of curvature. Δp = (F_st)_z/(π√(2R₁ζ)√(2R₂ζ)) (the exact pressure force per
    unit Δp). ``exact=False`` (default) uses the book's small-ζ evaluation of (4.98) (then Δp = σ(1/R₁ + 1/R₂) for any
    ζ); ``exact=True`` integrates round the exact rim: it differs from (1.5) by O(ζ/R) (−0.1 % at ζ/R = 1e-3) and
    converges to it as ζ → 0. Returns Δp [Pa]. Validation: V1 sphere R₁ = R₂ = 1 mm, σ = 0.0728 → 145.6 Pa; V3 exact →
    (1.5) at O(ζ/R) (slope 1). Label: analytic, converged.
    """
    Fst = cap_surface_tension_force(sigma, R1, R2, zeta, exact)
    return float(Fst / (np.pi * np.sqrt(2 * R1 * zeta) * np.sqrt(2 * R2 * zeta)))


def capillary_length(sigma, rho, g: float = G0, rho_other: float = 0.0):
    """Length over which gravity and surface tension are comparable, √(σ/((ρ − ρ_other)g)) [m] (§4.10, before
    Example 4.7; the δ of Example 4.7).

    Book: §4.10, "(σ/ρg)^{1/2}" (the book's air-bubble-in-water number lives in the private JSON).
    Parameters: sigma [N/m]; rho [kg/m³]; g [m/s²]; rho_other : density of the other fluid [kg/m³].
    Validation: V1 ``capillary_length(0.0728, 998.0)`` = 2.7269 mm; V5 water at 20 °C with IAPWS σ = 72.74 mN/m gives
    2.72 mm (Wikipedia "Capillary length": 2.71 mm). Label: analytic, benchmark.
    """
    return _S(np.sqrt(_F(sigma) / ((_F(rho) - float(rho_other)) * float(g))))
