"""Leibniz's theorem and the Reynolds transport theorem with concrete moving, deforming control volumes.

Book: Ch. 3 §3.6, Eqs. (3.30)–(3.35), Figs. 3.17–3.19, Example 3.2 (rendered pages chapters/pages/ch03/p112–p116).

A control volume V*(t) with closed surface A*(t), outward unit normal n and surface velocity b (Fig. 3.18) is an
object with

* ``volume_nodes(t, n)`` → (X (d, M), w (M,)): quadrature points and weights of ∫_{V*(t)} · dV [m³] (area in 2-D);
* ``surface_nodes(t, n)`` → (X (d, K), nrm (d, K), dA (K,), b (d, K)): points, outward normals, area elements [m²]
  (arc length in 2-D) and surface velocities [m/s];
* ``volume(t)`` and ``volume_rate(t)``: exact V* and dV*/dt for the checks.

Every shape is a parametric map X(ξ, t) from fixed reference coordinates ξ; its surface velocity is b = ∂X/∂t at fixed
ξ, so b is the actual velocity of the surface points (b need not follow the fluid; only b·n matters in (3.35)).
Quadrature: Gauss–Legendre in bounded directions, the midpoint rule in periodic angles (spectrally accurate), all
with ``n`` nodes per direction (2n in azimuth). Field callables follow ``core.kinematics``: F(x, t) with x (d, M).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, NamedTuple

import numpy as np
from scipy.integrate import quad

from ._util import as_scalar_if_0d
from .integral_theorems import gauss_legendre_nodes, midpoint_nodes

__all__ = ["LeibnizTerms", "leibniz_terms", "leibniz_check", "ControlVolume", "MovingBox", "GrowingSphere",
           "GrowingCylinder", "GrowingCone", "MovingEllipse2D", "volume_integral", "surface_flux_term",
           "volume_rate_term", "RTTTerms", "reynolds_transport", "volume_integral_rate_fd", "RTTCheck", "rtt_check",
           "swept_volume_integral", "swept_terms", "swept_terms_sphere", "MaterialVolumeRate",
           "material_volume_rate", "rtt_ellipse_2d"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731


# ======================================================================================================================
# Leibniz's theorem (3.30)
# ======================================================================================================================
class LeibnizTerms(NamedTuple):
    """The three parts of Eq. (3.30): ``interior`` ∫ₐᵇ ∂F/∂t dx, ``upper`` (db/dt) F(b, t), ``lower`` (da/dt) F(a, t)
    (entering with a minus sign) and ``total`` = interior + upper − lower.
    Book: §3.6, Eq. (3.30). Label: analytic.
    """

    interior: float
    upper: float
    lower: float
    total: float


def leibniz_terms(F: Callable, dFdt: Callable, a: float, b: float, dadt: float, dbdt: float, t: float,
                  n: int = 64) -> LeibnizTerms:
    """Right-hand side of Leibniz's theorem, Eq. (3.30), term by term.

    Book: §3.6, Eq. (3.30), Fig. 3.17: d/dt ∫_{a(t)}^{b(t)} F(x, t) dx = ∫ₐᵇ ∂F/∂t dx + (db/dt) F(b, t)
    − (da/dt) F(a, t): the change inside, plus the gain at the moving upper limit, minus the loss at the lower one.
    (The book cites Riley et al. for the proof; ours is D21.)

    Parameters
    ----------
    F, dFdt : callables (x, t) with x a float or an (N,) array of positions [m]; units [F], [F/s]
    a, b : limits at time t [m];  dadt, dbdt : their rates [m/s];  t : time [s]
    n : Gauss–Legendre nodes for the interior integral (exact for polynomials of degree ≤ 2n − 1)

    Returns
    -------
    LeibnizTerms(interior, upper, lower, total) [F·m/s]; ``lower`` = ȧ F(a, t) is the term that is **subtracted**.

    Assumptions: F and ∂F/∂t continuous on [a, b]; a, b differentiable.
    Validation: V1/V2 sympy d/dt of closed-form integrals (three cases); V3 ``leibniz_check`` (central difference of the
    integral) agrees at order 2; the slab RTT reduces to it. Label: analytic, symbolic, converged.
    """
    xg, wg = gauss_legendre_nodes(float(a), float(b), n)
    vals = np.broadcast_to(_F(dFdt(xg, t)), xg.shape)
    interior = float(np.sum(vals * wg))  # ∫ₐᵇ ∂F/∂t dx
    upper = float(dbdt) * float(F(float(b), t))  # gain at the moving upper limit
    lower = float(dadt) * float(F(float(a), t))  # loss at the moving lower limit
    return LeibnizTerms(interior, upper, lower, interior + upper - lower)  # Eq. (3.30)


def leibniz_check(F: Callable, a_fn: Callable, b_fn: Callable, t: float, dt: float = 1e-4) -> float:
    """Left side of (3.30) measured directly: the central difference [I(t + dt) − I(t − dt)]/(2dt) of
    I(τ) = ∫_{a(τ)}^{b(τ)} F(x, τ) dx (each integral by ``quad``). Independent of ``leibniz_terms``.

    Units [F·m/s]; error O(dt²) + quad noise/dt. Label: converged.
    Book: §3.6, Eq. (3.30) (left side).
    """
    def I(tau):
        return quad(lambda x: F(x, tau), a_fn(tau), b_fn(tau), epsabs=1e-14, epsrel=1e-13)[0]
    return float((I(t + dt) - I(t - dt)) / (2.0 * dt))


# ======================================================================================================================
# Control volumes (Fig. 3.18)
# ======================================================================================================================
def _lin(v0, rate, t) -> np.ndarray:
    return _F(v0) + _F(rate) * float(t)


@dataclass
class ControlVolume:
    """Base class of a moving, deforming control volume V*(t) with surface A*(t) (Fig. 3.18); see the module doc."""

    dim: int = field(default=3, init=False)

    def volume_nodes(self, t: float, n: int = 24):  # pragma: no cover - interface
        raise NotImplementedError

    def surface_nodes(self, t: float, n: int = 24):  # pragma: no cover - interface
        raise NotImplementedError

    def volume(self, t: float) -> float:  # pragma: no cover - interface
        raise NotImplementedError

    def volume_rate(self, t: float) -> float:  # pragma: no cover - interface
        raise NotImplementedError


@dataclass
class MovingBox(ControlVolume):
    """Box with side lengths L(t) = lengths + rates·t and lower corner c(t) = origin + velocity·t (each face moves
    rigidly with its own velocity: a point at fractional position s has b = velocity + s ⊙ rates). Units m, m/s.
    Fixed box: rates = velocity = 0 (then b = 0 and d/dt passes inside the integral).
    Book: §3.6, Fig. 3.18 (control volume V*, surface A*, velocity b). Label: analytic (volume, volume_rate exact).
    """

    lengths: tuple = (1.0, 1.0, 1.0)
    rates: tuple = (0.0, 0.0, 0.0)
    velocity: tuple = (0.0, 0.0, 0.0)
    origin: tuple = (0.0, 0.0, 0.0)

    @property
    def Ldot(self):
        return self.rates

    @property
    def U(self):
        return self.velocity

    def _geom(self, t):
        return _lin(self.origin, self.velocity, t), _lin(self.lengths, self.rates, t)

    def volume_nodes(self, t, n=16):
        c, L = self._geom(t)
        s, w = gauss_legendre_nodes(0.0, 1.0, n)
        S = np.array(np.meshgrid(s, s, s, indexing="ij")).reshape(3, -1)
        W = np.einsum("i,j,k->ijk", w, w, w).ravel() * np.prod(L)
        return c[:, None] + S * L[:, None], W

    def surface_nodes(self, t, n=16):
        c, L = self._geom(t)
        s, w = gauss_legendre_nodes(0.0, 1.0, n)
        A, B = np.meshgrid(s, s, indexing="ij")
        WA = np.einsum("i,j->ij", w, w).ravel()
        Xs, Ns, dAs, bs = [], [], [], []
        for k in range(3):
            i, j = [m for m in range(3) if m != k]
            for side in (0.0, 1.0):
                S = np.zeros((3, A.size))
                S[i], S[j], S[k] = A.ravel(), B.ravel(), side
                nrm = np.zeros((3, A.size))
                nrm[k] = 1.0 if side == 1.0 else -1.0
                Xs.append(c[:, None] + S * L[:, None])
                Ns.append(nrm)
                dAs.append(WA * L[i] * L[j])
                bs.append(_F(self.U)[:, None] + S * _F(self.Ldot)[:, None])  # b = ∂X/∂t at fixed s
        return np.hstack(Xs), np.hstack(Ns), np.concatenate(dAs), np.hstack(bs)

    def volume(self, t):
        return float(np.prod(self._geom(t)[1]))

    def volume_rate(self, t):
        L, Ld = self._geom(t)[1], _F(self.Ldot)
        return float(Ld[0] * L[1] * L[2] + L[0] * Ld[1] * L[2] + L[0] * L[1] * Ld[2])


@dataclass
class GrowingSphere(ControlVolume):
    """Sphere of radius R(t) = R0 + Rdot t centred at c(t) = center + U t (a balloon; b = U + Rdot e_r). Units m, m/s.
    Book: §3.6, Fig. 3.18. Label: analytic (volume, volume_rate exact).
    """

    R0: float = 1.0
    Rdot: float = 0.0
    center: tuple = (0.0, 0.0, 0.0)
    U: tuple = (0.0, 0.0, 0.0)

    def _geom(self, t):
        return _lin(self.center, self.U, t), float(self.R0 + self.Rdot * t)

    @staticmethod
    def _dirs(n):
        mu, wmu = gauss_legendre_nodes(-1.0, 1.0, n)  # μ = cos θ
        ph, wph = midpoint_nodes(0.0, 2.0 * np.pi, 2 * n)
        MU, PH = np.meshgrid(mu, ph, indexing="ij")
        st = np.sqrt(1.0 - MU ** 2)
        er = np.stack([st * np.cos(PH), st * np.sin(PH), MU]).reshape(3, -1)
        return er, np.einsum("i,j->ij", wmu, wph).ravel()

    def volume_nodes(self, t, n=16):
        c, R = self._geom(t)
        er, wd = self._dirs(n)
        rho, wr = gauss_legendre_nodes(0.0, 1.0, n)
        X = c[:, None, None] + R * rho[None, :, None] * er[:, None, :]
        W = (R ** 3) * np.einsum("i,j->ij", rho ** 2 * wr, wd)
        return X.reshape(3, -1), W.ravel()

    def surface_nodes(self, t, n=16):
        c, R = self._geom(t)
        er, wd = self._dirs(n)
        b = _F(self.U)[:, None] + self.Rdot * er  # b = ∂X/∂t = U + Ṙ e_r
        return c[:, None] + R * er, er, R ** 2 * wd, b

    def volume(self, t):
        return float(4.0 / 3.0 * np.pi * self._geom(t)[1] ** 3)

    def volume_rate(self, t):
        return float(4.0 * np.pi * self._geom(t)[1] ** 2 * self.Rdot)


@dataclass
class GrowingCylinder(ControlVolume):
    """Circular cylinder, axis along z, radius R(t) = R0 + Rdot t, length L(t) = L + Ldot t, bottom-cap centre at
    base + U t (the whole translating at U). Units m, m/s.
    Book: §3.6, Fig. 3.18. Label: analytic (volume, volume_rate exact).
    """

    R0: float = 1.0
    Rdot: float = 0.0
    L: float = 1.0
    Ldot: float = 0.0
    base: tuple = (0.0, 0.0, 0.0)
    U: tuple = (0.0, 0.0, 0.0)

    def _geom(self, t):
        return _lin(self.base, self.U, t), float(self.R0 + self.Rdot * t), float(self.L + self.Ldot * t)

    def volume_nodes(self, t, n=16):
        c, R, L = self._geom(t)
        rho, wr = gauss_legendre_nodes(0.0, 1.0, n)
        ph, wp = midpoint_nodes(0.0, 2.0 * np.pi, 2 * n)
        ze, wz = gauss_legendre_nodes(0.0, 1.0, n)
        P, PH, Z = np.meshgrid(rho, ph, ze, indexing="ij")
        X = np.stack([P * R * np.cos(PH), P * R * np.sin(PH), Z * L]).reshape(3, -1) + c[:, None]
        W = (np.einsum("i,j,k->ijk", rho * wr, wp, wz) * R ** 2 * L).ravel()
        return X, W

    def surface_nodes(self, t, n=16):
        c, R, L = self._geom(t)
        U = _F(self.U)[:, None]
        rho, wr = gauss_legendre_nodes(0.0, 1.0, n)
        ph, wp = midpoint_nodes(0.0, 2.0 * np.pi, 2 * n)
        ze, wz = gauss_legendre_nodes(0.0, 1.0, n)
        # side ρ = 1
        PH, Z = np.meshgrid(ph, ze, indexing="ij")
        eR = np.stack([np.cos(PH), np.sin(PH), np.zeros_like(PH)]).reshape(3, -1)
        Xs = c[:, None] + R * eR + np.stack([0 * PH, 0 * PH, Z * L]).reshape(3, -1)
        bs = U + self.Rdot * eR + np.stack([0 * PH, 0 * PH, Z * self.Ldot]).reshape(3, -1)
        dAs = (np.einsum("i,j->ij", wp, wz) * R * L).ravel()
        # caps ζ = 0 (n = −e_z) and ζ = 1 (n = +e_z)
        P, PH2 = np.meshgrid(rho, ph, indexing="ij")
        eR2 = np.stack([np.cos(PH2), np.sin(PH2), np.zeros_like(PH2)]).reshape(3, -1)
        rr = P.ravel()
        dAc = (np.einsum("i,j->ij", rho * wr, wp) * R ** 2).ravel()
        Xb = c[:, None] + R * rr * eR2
        Xt = Xb + np.array([0.0, 0.0, L])[:, None]
        bb = U + self.Rdot * rr * eR2
        bt = bb + np.array([0.0, 0.0, self.Ldot])[:, None]
        nb = np.tile(np.array([0.0, 0.0, -1.0])[:, None], (1, rr.size))
        return (np.hstack([Xs, Xb, Xt]), np.hstack([eR, nb, -nb]), np.concatenate([dAs, dAc, dAc]),
                np.hstack([bs, bb, bt]))

    def volume(self, t):
        _, R, L = self._geom(t)
        return float(np.pi * R ** 2 * L)

    def volume_rate(self, t):
        _, R, L = self._geom(t)
        return float(np.pi * (2.0 * R * self.Rdot * L + R ** 2 * self.Ldot))


@dataclass
class GrowingCone(ControlVolume):
    """Right circular cone of fixed height h, apex at the origin, axis +z, base (at z = h) radius r(t) = r0 + rdot t
    (Example 3.2, Fig. 3.19). A point at height z on the side moves at b = (z/h) ṙ e_R; points of the base move
    radially in its plane (b = ρ ṙ e_R, so b·n = 0 there — the book says b = 0, only b·n = 0 is true). Units m, m/s.
    Book: §3.6, Example 3.2, Fig. 3.19. Label: analytic (volume, volume_rate exact).
    """

    r0: float = 0.5
    rdot: float = 0.1
    h: float = 1.0

    def radius(self, t):
        return float(self.r0 + self.rdot * t)

    def volume_nodes(self, t, n=16):
        r, h = self.radius(t), float(self.h)
        rho, wr = gauss_legendre_nodes(0.0, 1.0, n)
        ph, wp = midpoint_nodes(0.0, 2.0 * np.pi, 2 * n)
        ze, wz = gauss_legendre_nodes(0.0, 1.0, n)
        P, PH, Z = np.meshgrid(rho, ph, ze, indexing="ij")
        X = np.stack([P * Z * r * np.cos(PH), P * Z * r * np.sin(PH), Z * h]).reshape(3, -1)
        W = (np.einsum("i,j,k->ijk", rho * wr, wp, ze ** 2 * wz) * r ** 2 * h).ravel()  # dV = ρ ζ² r² h dρ dφ dζ
        return X, W

    def surface_nodes(self, t, n=16):
        r, h, rd = self.radius(t), float(self.h), float(self.rdot)
        th = np.arctan2(r, h)  # half-angle: tan θ = r/h
        ph, wp = midpoint_nodes(0.0, 2.0 * np.pi, 2 * n)
        ze, wz = gauss_legendre_nodes(0.0, 1.0, n)
        rho, wr = gauss_legendre_nodes(0.0, 1.0, n)
        PH, Z = np.meshgrid(ph, ze, indexing="ij")
        eR = np.stack([np.cos(PH), np.sin(PH), np.zeros_like(PH)]).reshape(3, -1)
        zz = Z.ravel()
        Xs = np.vstack([zz * r * eR[0], zz * r * eR[1], zz * h])
        ns = np.vstack([np.cos(th) * eR[0], np.cos(th) * eR[1], -np.sin(th) * np.ones_like(zz)])  # n = e_R cos θ − e_z sin θ
        bs = zz * rd * eR  # b = (z/h) ṙ e_R
        dAs = (np.einsum("i,j->ij", wp, ze * wz) * r * h / np.cos(th)).ravel()  # dA = z tan θ dφ dz/cos θ
        P, PH2 = np.meshgrid(rho, ph, indexing="ij")
        eR2 = np.stack([np.cos(PH2), np.sin(PH2), np.zeros_like(PH2)]).reshape(3, -1)
        rr = P.ravel()
        Xb = np.vstack([rr * r * eR2[0], rr * r * eR2[1], np.full(rr.size, h)])
        nb = np.tile(np.array([0.0, 0.0, 1.0])[:, None], (1, rr.size))
        bb = rr * rd * eR2  # base points move radially in the base plane: b·n = 0
        dAb = (np.einsum("i,j->ij", rho * wr, wp) * r ** 2).ravel()
        return np.hstack([Xs, Xb]), np.hstack([ns, nb]), np.concatenate([dAs, dAb]), np.hstack([bs, bb])

    def volume(self, t):
        return float(np.pi * self.h * self.radius(t) ** 2 / 3.0)

    def volume_rate(self, t):
        return float(2.0 / 3.0 * np.pi * self.h * self.radius(t) * self.rdot)


@dataclass
class MovingEllipse2D(ControlVolume):
    """Plane control "volume": the ellipse with semi-axes a(t) (along x) and b(t) (along y) centred at c(t), given as
    callables. Its "volume" is the area [m²], its "surface" the boundary curve (dA = arc length ds [m]).
    Rates ȧ, ḃ, ċ are taken by central differences with step ``h`` [s] (or exactly if ``rates`` callables
    t → (ȧ, ḃ, ċx, ċy) are given). Quadrature with ``n``: n Gauss–Legendre radial × 2n midpoint azimuthal nodes for
    the area, 8n midpoint nodes on the boundary. Book: (3.35) read in two dimensions (E7's stage).
    Label: analytic (area, area rate exact).
    """

    a_fn: Callable = lambda t: 1.0  # noqa: E731
    b_fn: Callable = lambda t: 0.6  # noqa: E731
    center_fn: Callable = lambda t: (0.0, 0.0)  # noqa: E731
    h: float = 1e-6
    rates: Callable | None = None

    def __post_init__(self):
        self.dim = 2

    def _geom(self, t):
        a, b, c = float(self.a_fn(t)), float(self.b_fn(t)), _F(self.center_fn(t))
        if self.rates is not None:
            ad, bd, cx, cy = (float(v) for v in self.rates(t))
            cd = np.array([cx, cy])
        else:
            hh = self.h
            ad = (float(self.a_fn(t + hh)) - float(self.a_fn(t - hh))) / (2 * hh)
            bd = (float(self.b_fn(t + hh)) - float(self.b_fn(t - hh))) / (2 * hh)
            cd = (_F(self.center_fn(t + hh)) - _F(self.center_fn(t - hh))) / (2 * hh)
        return a, b, c, ad, bd, cd

    def volume_nodes(self, t, n=16):
        a, b, c, *_ = self._geom(t)
        rho, wr = gauss_legendre_nodes(0.0, 1.0, n)
        ph, wp = midpoint_nodes(0.0, 2.0 * np.pi, 2 * n)
        P, PH = np.meshgrid(rho, ph, indexing="ij")
        X = np.stack([c[0] + P * a * np.cos(PH), c[1] + P * b * np.sin(PH)]).reshape(2, -1)
        return X, (np.einsum("i,j->ij", rho * wr, wp) * a * b).ravel()

    def surface_nodes(self, t, n=16):
        a, b, c, ad, bd, cd = self._geom(t)
        ph, wp = midpoint_nodes(0.0, 2.0 * np.pi, 8 * n)
        X = np.stack([c[0] + a * np.cos(ph), c[1] + b * np.sin(ph)])
        T = np.stack([-a * np.sin(ph), b * np.cos(ph)])  # dX/dφ
        ds = np.sqrt(np.sum(T ** 2, axis=0))
        nrm = np.stack([b * np.cos(ph), a * np.sin(ph)]) / ds  # outward normal
        bvel = cd[:, None] + np.stack([ad * np.cos(ph), bd * np.sin(ph)])  # b = ∂X/∂t
        return X, nrm, ds * wp, bvel

    def volume(self, t):
        a, b, *_ = self._geom(t)
        return float(np.pi * a * b)

    def volume_rate(self, t):
        a, b, _, ad, bd, _ = self._geom(t)
        return float(np.pi * (ad * b + a * bd))


# ======================================================================================================================
# Reynolds transport theorem (3.31)–(3.35)
# ======================================================================================================================
def _Fv(F: Callable, X: np.ndarray, t: float) -> np.ndarray:
    val = _F(F(X, t))
    return np.broadcast_to(val, (X.shape[1],)) if val.ndim == 0 else val


def volume_integral(F: Callable, cv: ControlVolume, t: float, n: int = 24, time: float | None = None) -> float:
    """∫_{V*(t)} F(x, τ) dV with τ = ``time`` (default t) — the quadrature behind (3.31). Units [F·m³].
    Book: §3.6, Eq. (3.31) (the integrals inside the braces). Label: analytic (Gauss–Legendre, exact for polynomials).
    """
    X, W = cv.volume_nodes(t, n)
    return float(np.sum(_Fv(F, X, t if time is None else time) * W))


def surface_flux_term(F: Callable, cv: ControlVolume, t: float, n: int = 24) -> float:
    """Boundary term of (3.35): ∮_{A*(t)} F(x, t) b·n dA — positive where A* advances (b·n > 0), negative where it
    retreats. Units [F·m³/s]. Label: analytic.
    Book: §3.6, Eq. (3.35) (second term on the right).
    """
    X, N, dA, b = cv.surface_nodes(t, n)
    return float(np.sum(_Fv(F, X, t) * np.einsum("ik,ik->k", b, N) * dA))  # ∮ F b·n dA


def volume_rate_term(dFdt: Callable, cv: ControlVolume, t: float, n: int = 24) -> float:
    """Volume term of (3.35): ∫_{V*(t)} ∂F/∂t dV. Units [F·m³/s].
    Book: §3.6, Eq. (3.35) (first term on the right). Label: analytic.
    """
    return volume_integral(dFdt, cv, t, n)


class RTTTerms(NamedTuple):
    """Eq. (3.35): ``total`` d/dt ∫_{V*} F dV = ``volume_term`` ∫ ∂F/∂t dV + ``surface_term`` ∮ F b·n dA.
    Book: §3.6, Eq. (3.35). Label: analytic.
    """

    volume_term: float
    surface_term: float
    total: float


def reynolds_transport(F: Callable, dFdt: Callable, cv: ControlVolume, t: float, n: int = 24) -> RTTTerms:
    """Right-hand side of the Reynolds transport theorem, Eq. (3.35).

    Book: §3.6, Eq. (3.35): d/dt ∫_{V*(t)} F(x, t) dV = ∫_{V*(t)} ∂F(x, t)/∂t dV + ∫_{A*(t)} F(x, t) b·n dA — "what
    changes in place + what the moving surface sweeps in". d/dt passes inside the integral only for a fixed volume
    (b = 0). With F = 1 it is dV*/dt = ∮ b·n dA (N53).

    Parameters
    ----------
    F, dFdt : callables (x, t), x (d, M) [F], [F/s];  cv : a ControlVolume;  t [s];  n : nodes per direction

    Returns
    -------
    RTTTerms(volume_term, surface_term, total) [F·m³/s].

    Assumptions: F continuous with continuous ∂F/∂t; A* piecewise smooth; b continuous on each face.
    Validation: V1/V3 total equals ``volume_integral_rate_fd`` (the definition (3.31) by finite differences — an
    independent route) for box, cylinder, sphere, cone and ellipse; V2 closed forms by sympy; fixed volume → surface
    term 0; F = const → F dV*/dt. Label: analytic, converged, symbolic.
    """
    vt = volume_rate_term(dFdt, cv, t, n)
    st = surface_flux_term(F, cv, t, n)
    return RTTTerms(vt, st, vt + st)  # Eq. (3.35)


def volume_integral_rate_fd(F: Callable, cv: ControlVolume, t: float, dt: float = 1e-4, n: int = 24) -> float:
    """The definition (3.31) evaluated with a central difference: [I(t + dt) − I(t − dt)]/(2dt),
    I(τ) = ∫_{V*(τ)} F(x, τ) dV (the nodes move with the shape, so both integrals use the same rule).

    Book: §3.6, Eq. (3.31). Units [F·m³/s]; error O(dt²). Label: converged.
    """
    return (volume_integral(F, cv, t + dt, n) - volume_integral(F, cv, t - dt, n)) / (2.0 * dt)


class RTTCheck(NamedTuple):
    """Both sides of (3.35): ``lhs_fd`` (the definition (3.31) by a central difference) and ``rhs_total`` (the RTT).
    Book: §3.6, Eqs. (3.31), (3.35). Label: converged.
    """

    lhs_fd: float
    rhs_total: float

    @property
    def abs_error(self) -> float:
        return abs(self.lhs_fd - self.rhs_total)

    @property
    def rel_error(self) -> float:
        return abs(self.lhs_fd - self.rhs_total) / max(abs(self.rhs_total), 1e-300)


def rtt_check(F: Callable, dFdt: Callable, cv: ControlVolume, t: float, dt: float = 1e-4, n: int = 24) -> RTTCheck:
    """Both sides of (3.35): the measured rate (3.31) by finite differences and the RTT total.

    Returns RTTCheck(lhs_fd, rhs_total) (unpacks as a pair; ``.abs_error``, ``.rel_error`` properties).
    Validation: V3 agreement O(dt²) + quadrature error. Label: converged.
    Book: §3.6, Eqs. (3.31), (3.35).
    """
    return RTTCheck(volume_integral_rate_fd(F, cv, t, dt, n), reynolds_transport(F, dFdt, cv, t, n).total)


def swept_volume_integral(F: Callable, cv: ControlVolume, t: float, dt: float, n: int = 24) -> float:
    """∫_{ΔV} F(x, t) dV over the (signed) volume increment ΔV = V*(t + dt) − V*(t), Eq. (3.34)'s left side, computed
    exactly (to quadrature) as ∫_{V*(t+dt)} F(x, t) dV − ∫_{V*(t)} F(x, t) dV. Where the surface retreats (b·n < 0)
    the sliver counts negatively. Compare with ``surface_flux_term(F, cv, t)·dt``: the gap is O(dt²) (D22).
    Units [F·m³]. Label: converged.
    Book: §3.6, Eq. (3.34) (left side), Fig. 3.18.
    """
    return volume_integral(F, cv, t + dt, n, time=t) - volume_integral(F, cv, t, n, time=t)


def swept_terms(F: Callable, dFdt: Callable, cv: ControlVolume, t: float, dt: float, n: int = 24) -> dict:
    """The four terms of Eq. (3.32) and the pieces of (3.33)–(3.34) for one Δt — orders of smallness made visible.

    Book: §3.6, Eq. (3.32): ∫_{V*(t+Δt)} F(x, t + Δt) dV ≅ T1 + T2 + T3 + T4 with T1 = ∫_{V*(t)} F dV,
    T2 = ∫_{V*(t)} Δt ∂F/∂t dV, T3 = ∫_{ΔV} F dV, T4 = ∫_{ΔV} Δt ∂F/∂t dV. T4 is second order (dropped in (3.33));
    (3.34) replaces T3 by the swept sliver Σ F (bΔt·n) dA.

    Returns
    -------
    dict(T1, T2, T3, T4, exact = ∫_{V*(t+Δt)} F(x, t + Δt) dV, lhs = (exact − T1)/Δt (the difference quotient of
    (3.31), → RTT total as Δt → 0), taylor_residual = exact − (T1+T2+T3+T4) (O(Δt²)), sliver = Δt·∮F b·n dA,
    sliver_error = T3 − sliver (O(Δt²))).

    Validation: V3 T4 and both residuals fall with slope 2 on log–log axes. Label: converged.
    """
    exact = volume_integral(F, cv, t + dt, n)
    T1 = volume_integral(F, cv, t, n)
    T2 = dt * volume_integral(dFdt, cv, t, n)
    T3 = swept_volume_integral(F, cv, t, dt, n)
    T4 = dt * swept_volume_integral(dFdt, cv, t, dt, n)
    sliver = dt * surface_flux_term(F, cv, t, n)  # Eq. (3.34)
    return {"T1": T1, "T2": T2, "T3": T3, "T4": T4, "exact": exact, "lhs": (exact - T1) / dt,
            "taylor_residual": exact - (T1 + T2 + T3 + T4), "sliver": sliver, "sliver_error": T3 - sliver}


def swept_terms_sphere(R: float, Rdot: float, F: Callable, dFdt: Callable, t: float, dt: float, n: int = 24,
                       center=(0.0, 0.0, 0.0)) -> dict:
    """:func:`swept_terms` for a balloon ``GrowingSphere(R, Rdot)``: radius R at t = 0 growing at Rdot, so the radius
    at time t is R + Rdot·t (C15's slider figure, E7's limit view). Units R [m], Rdot [m/s], t, dt [s].
    Label: converged.
    Book: §3.6, Eqs. (3.32)–(3.34).
    """
    return swept_terms(F, dFdt, GrowingSphere(R0=float(R), Rdot=float(Rdot), center=tuple(center)), t, dt, n)


class MaterialVolumeRate(NamedTuple):
    """F = 1, b = u in (3.35): ``surface_flux`` ∮ u·n dA and ``volume_div`` ∫ ∇·u dV (equal by Gauss).
    Book: §3.6, Eq. (3.35) with F = 1, b = u; Eq. (3.14). Label: analytic.
    """

    surface_flux: float
    volume_div: float


def material_volume_rate(u: Callable, cv: ControlVolume, t: float, n: int = 24, h: float = 1e-5) -> MaterialVolumeRate:
    """Rate of change of a material volume (b = u, F = 1 in (3.35)): dV/dt = ∮ u·n dA = ∫ ∇·u dV.

    Book: §3.6 interpretation (Exercise 3.28; our D23): for a small material volume this is (3.14),
    (1/δV) D(δV)/Dt = ∇·u. The surface integral is taken over the current shape of ``cv`` (only its geometry at t is
    used); ∇·u by 2nd-order central differences (step h [m]).

    Returns MaterialVolumeRate(surface_flux, volume_div) [m³/s].
    Validation: V1/V4 the two agree (Gauss); divided by δV they → ∇·u as δV → 0; incompressible u gives 0.
    Label: analytic, conserved.
    """
    Xs, N, dA, _ = cv.surface_nodes(t, n)
    U = _F(u(Xs, t))
    flux = float(np.sum(np.einsum("ik,ik->k", U, N) * dA))
    Xv, W = cv.volume_nodes(t, n)
    d = Xv.shape[0]
    div = 0.0
    for i in range(d):
        e = np.zeros((d, 1))
        e[i] = h
        div = div + (_F(u(Xv + e, t))[i] - _F(u(Xv - e, t))[i]) / (2.0 * h)
    return MaterialVolumeRate(flux, float(np.sum(div * W)))


def rtt_ellipse_2d(F: Callable, dFdt: Callable, a: float, b: float, adot: float, bdot: float, c=(0.0, 0.0),
                   cdot=(0.0, 0.0), t: float = 0.0, n: int = 256) -> RTTTerms:
    """E7's 2-D control volume: (3.35) for an ellipse with semi-axes a, b (along x, y) centred at c, whose sizes and
    centre move at ȧ, ḃ, ċ, at the instant t.

    Book: §3.6, Eq. (3.35) read in the plane: d/dt ∬ F dA = ∬ ∂F/∂t dA + ∮ F b·n ds. The boundary point
    (c_x + a cos s, c_y + b sin s) moves at b = ċ + (ȧ cos s, ḃ sin s); outward normal ∝ (b cos s, a sin s).

    Parameters
    ----------
    F, dFdt : callables (x (2, M), t) [F], [F/s];  a, b [m];  adot, bdot [m/s];  c [m], cdot [m/s] (2-vectors);
    t : time [s];  n : boundary nodes (the area uses n/8 radial × n/4 azimuthal nodes)

    Returns
    -------
    RTTTerms(volume_term, surface_term, total) [F·m²/s]. Scalar-callable (explainer parity).

    Validation: V3 equals the finite difference of the area integral (``volume_integral_rate_fd`` with a
    ``MovingEllipse2D`` moving at these rates). Label: converged.
    """
    a0, b0, t0 = float(a), float(b), float(t)
    c0, cd = _F(c), _F(cdot)
    cv = MovingEllipse2D(a_fn=lambda tt: a0 + float(adot) * (tt - t0), b_fn=lambda tt: b0 + float(bdot) * (tt - t0),
                         center_fn=lambda tt: c0 + cd * (tt - t0),
                         rates=lambda tt: (float(adot), float(bdot), float(cd[0]), float(cd[1])))
    return reynolds_transport(F, dFdt, cv, t0, max(4, int(n) // 8))
