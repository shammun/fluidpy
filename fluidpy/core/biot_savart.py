"""Velocity from vorticity: the Poisson equation ∇²u = −∇ × ω, its Green's-function solution, the Biot–Savart law
(volume, filament, straight segment, circular ring), point vortices with image systems, and vortex sheets.

Book: Ch. 5 §5.5, Eqs. (5.14)–(5.17), Fig. 5.8; §5.7, Figs. 5.11–5.15; §5.8, Fig. 5.16 (rendered pages
chapters/pages/ch05/p208–p218). Exercises 5.9, 5.13, 5.14, 5.19, 5.20 give test cases (their printed answers stay in
the private ``tests/book_values_ch05.json``).

Book slip handled here (analysis §9 note 1): (5.14) is printed with −1/(4π); with ∇²u = −∇ × ω and the Green's
function of ∇² (Exercise 5.9: φ = −(1/4π)∫q/|x − x′|d³x′ solves ∇²φ = q) the factor is **+1/(4π)**. The book's next
line has a compensating slip (it writes +((x − x′)/|x − x′|³) × ω where the product rule gives −(…) × ω), so (5.16)
and (5.17) as printed are correct. :func:`velocity_from_curl_omega` exposes ``sign`` so that a test can show the
printed sign reverses the swirl.

Conventions: 3-D points (3,) or (3, N), 2-D points (2,) or (2, N) (components first); circulation Γ [m²/s]
counterclockwise positive about +z for plane vortices and along e_ω for filaments (right-hand rule); a vortex sheet's
strength γ [m/s] is the circulation per unit length, counterclockwise positive (γ = u_below − u_above, the book text's
u₂ − u₁). Kernels are singular on the vorticity itself: ``eps`` > 0 smooths them (Rosenhead–Moore/Krasny style).

Reuse: Ch. 6 (images, vortex + stream), Ch. 10 (vortex methods, vorticity–stream-function solvers), Ch. 11 (discrete
vortex Kelvin–Helmholtz roll-up), Ch. 13 (point-vortex / PV inversion), **Ch. 14 (lifting line, horseshoe vortices,
induced drag, trailing sheet)**.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipe, ellipk

from ._util import as_scalar_if_0d
from .integral_theorems import gauss_legendre_nodes

__all__ = ["poisson_green_3d", "cylinder_quadrature", "gaussian_tube_fields", "velocity_from_curl_omega",
           "velocity_from_vorticity_fft", "biot_savart_volume", "biot_savart_2d", "segment_induced_velocity",
           "segment_speed", "filament_velocity", "filament_contributions", "filament_velocity_preset",
           "ring_axis_velocity", "ring_ring_velocity", "ring_self_velocity", "RING_CORES", "point_vortex_velocity",
           "point_vortex_rhs", "point_vortex_evolve", "point_vortex_invariants", "centre_of_vorticity",
           "wall_image_system", "circle_image_system", "image_vortices", "channel_image_velocity",
           "channel_image_velocity_exact", "vortex_sheet_velocity", "continuous_sheet_velocity",
           "FILAMENT_PRESETS", "FILAMENT_CLOSED", "filament_preset", "POINT_VORTEX_PRESETS", "point_vortex_preset"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d
_FOUR_PI = 4.0 * np.pi


def _pts(x, d):
    """Points as (d, N) and a flag telling whether a single point was given."""
    x_ = _F(x)
    single = x_.ndim == 1
    return (x_.reshape(d, 1) if single else x_.reshape(d, -1)), single


# ======================================================================================================================
# §5.5 Poisson equation, Green's function, (5.14)–(5.16)
# ======================================================================================================================
def poisson_green_3d(x, xp):
    """Free-space Green's function of the Laplacian, G(x, x′) = −1/(4π|x − x′|): ∇²G = δ(x − x′) (Exercise 5.9a).

    Book: §5.5, the Poisson equation ∇²u = −∇ × ω before (5.14), and Exercise 5.9 (φ = −(1/4π)∫q(x′)/|x − x′| d³x′
    solves ∇²φ = q). Parameters: x, xp points (3,) or (3, N) [m] (broadcast). Returns G [1/m].
    Validation: V2 sympy ∇²(1/r) = 0 for r > 0; V1 the flux of ∇G through any sphere about x′ is 1. Label: symbolic,
    analytic.
    """
    r = np.linalg.norm(_F(x) - _F(xp), axis=0)
    return _S(-1.0 / (_FOUR_PI * r))  # G = −1/(4π|x − x′|)


def cylinder_quadrature(radius: float, z0: float, z1: float, nr: int = 24, nphi: int = 32, nz: int = 24,
                        center=(0.0, 0.0)):
    """Tensor quadrature over a cylinder {R ≤ radius, z0 ≤ z ≤ z1}: Gauss–Legendre in R (weight R) and z, uniform
    (spectral) in φ. Returns (nodes (3, M) [m], weights (M,) [m³]). Used to evaluate (5.14) and (5.16) over a V′ that
    encloses a vortex tube segment (Fig. 5.8). Label: converged. Book: §5.5 (V′ of Fig. 5.8)."""
    r, wr = gauss_legendre_nodes(0.0, float(radius), int(nr))
    z, wz = gauss_legendre_nodes(float(z0), float(z1), int(nz))
    ph = 2.0 * np.pi * (np.arange(nphi) + 0.5) / nphi
    Rr, PH, Z = np.meshgrid(r, ph, z, indexing="ij")
    W = (wr[:, None, None] * Rr * (2.0 * np.pi / nphi) * wz[None, None, :])
    c = _F(center)
    nodes = np.stack([c[0] + Rr * np.cos(PH), c[1] + Rr * np.sin(PH), Z]).reshape(3, -1)
    return nodes, W.ravel()


def gaussian_tube_fields(Gamma: float = 1.0, sigma: float = 0.1, L: float = 4.0) -> dict:
    """A straight Gaussian vortex tube along z of length L (flat ends at z = ±L/2), for the quadrature tests of (5.14)
    and (5.16): ω = (Γ/πσ²)e^{−R²/σ²}e_z for |z| ≤ L/2, and its curl ∇ × ω = (∂ω_z/∂y, −∂ω_z/∂x, 0).

    This flat-ended tube is one segment V′ of a longer tube (Fig. 5.8), not a valid vorticity field on its own: ∇·ω
    has δ-sheets at z = ±L/2 where the tube is cut. It is used only as the V′ of the quadratures.
    The flat ends are normal to ω and the side lies outside the core, so the surface term of (5.15) vanishes and
    (5.14) (sign corrected) and (5.16) give the same velocity. A finite tube is not an infinite line: at the mid-plane
    the velocity is the segment law (our D13) with the core's enclosed circulation, computed exactly here as
    ``u_theta_reference(r)`` (1-D quadrature over the core after integrating z analytically); the infinite-line value
    Γ(1 − e^{−r²/σ²})/2πr is ``u_theta_infinite(r)`` (0.318310 m/s at r = 0.5 for Γ = 1, σ = 0.1 — reached only as
    L → ∞; with L = 4 the tube gives 0.3088 m/s).

    Returns dict(omega, curl_omega (callables of (3, M) points), bounds dict(radius = 6σ, z0 = −L/2, z1 = L/2) for
    :func:`cylinder_quadrature`, u_theta_reference, u_theta_infinite). Book: §5.5, Eqs. (5.14)–(5.16), Fig. 5.8.
    Label: analytic.
    """
    G, s, Lh = float(Gamma), float(sigma), 0.5 * float(L)

    def omega(P, t=0.0):
        P_ = _F(P)
        r2 = P_[0] ** 2 + P_[1] ** 2
        inside = np.abs(P_[2]) <= Lh
        wz = np.where(inside, G / (np.pi * s ** 2) * np.exp(-r2 / s ** 2), 0.0)
        return np.stack([0.0 * wz, 0.0 * wz, wz])

    def curl_omega(P, t=0.0):
        P_ = _F(P)
        w = omega(P_)[2]
        dw = -2.0 * w / s ** 2  # (1/R) dω_z/dR
        return np.stack([dw * P_[1], -dw * P_[0], 0.0 * w])  # ∇×(ω_z e_z) = (∂ω_z/∂y, −∂ω_z/∂x, 0)

    def u_theta_reference(r):
        from scipy.integrate import quad
        r = float(r)

        def ring(Rp):  # ∫ dφ′ ∫ dz′ of the (x − x′) kernel for a ring of vorticity at R′, times ω(R′) R′
            ph = np.linspace(0.0, 2.0 * np.pi, 257)[:-1]
            d2 = r ** 2 + Rp ** 2 - 2.0 * r * Rp * np.cos(ph)
            kz = 2.0 * Lh / (d2 * np.sqrt(d2 + Lh ** 2))  # ∫_{−L/2}^{L/2} dz′/(d² + z′²)^{3/2}
            val = np.mean((r - Rp * np.cos(ph)) * kz) * 2.0 * np.pi
            return G / (np.pi * s ** 2) * np.exp(-Rp ** 2 / s ** 2) * Rp * val
        v, _ = quad(ring, 0.0, min(8.0 * s, 0.999 * r), limit=200)  # vorticity beyond r neglected: use r ≳ 4σ
        return v / _FOUR_PI

    return dict(omega=omega, curl_omega=curl_omega, bounds=dict(radius=6.0 * s, z0=-Lh, z1=Lh),
                u_theta_reference=u_theta_reference,
                u_theta_infinite=lambda r: G * (-np.expm1(-float(r) ** 2 / s ** 2)) / (2.0 * np.pi * float(r)))


def velocity_from_curl_omega(curl_omega_fn: Callable, x, nodes, weights, sign: float = +1.0, t: float = 0.0,
                             chunk: int = 2_000_000):
    """Vorticity-induced velocity from the Green's-function solution of ∇²u = −∇ × ω, Eq. (5.14) with the sign fixed.

    Book: §5.5, Eq. (5.14), printed u(x, t) = −(1/4π)∫_V′ (∇′ × ω(x′, t))/|x − x′| d³x′. With ∇²u = −∇ × ω (the curl of
    ω for ∇·u = 0, identity (B.3.13)) and Exercise 5.9's solution φ = −(1/4π)∫q/|x − x′| of ∇²φ = q, q = −∇′ × ω gives
    u = **+**(1/4π)∫(∇′ × ω)/|x − x′| d³x′ (our D10). ``sign=+1`` (default) is the correct solution; ``sign=−1`` is the
    printed one, kept only so a test can show that it reverses the swirl.

    Parameters
    ----------
    curl_omega_fn : ∇ × ω at points x′ (3, M) → (3, M) [1/(s m)] (``t`` passed as the second argument)
    x : field point(s) (3,) or (3, N) [m];  nodes, weights : quadrature over V′ (e.g. :func:`cylinder_quadrature`)
    sign : +1 (correct) or −1 (as printed);  t : time [s];  chunk : max N·M kernel entries held at once

    Returns
    -------
    u (3,) or (3, N) [m/s].

    Validation: V1 Gaussian vortex tube: equals :func:`biot_savart_volume` on the same nodes and the exact u_θ; the
    printed sign gives −u_θ (wrong-variant check). Label: analytic.
    """
    X, single = _pts(x, 3)
    Pn, w = _F(nodes), _F(weights)
    C = _F(curl_omega_fn(Pn, t)) * w  # (3, M)
    out = np.zeros_like(X)
    step = max(1, chunk // max(Pn.shape[1], 1))
    for i in range(0, X.shape[1], step):
        Xi = X[:, i:i + step]
        r = np.linalg.norm(Xi[:, :, None] - Pn[:, None, :], axis=0)  # (n, M)
        out[:, i:i + step] = np.einsum("km,nm->kn", C, 1.0 / r)
    # DEVIATION: (5.14) is printed with −1/(4π); ∇²u = −∇×ω with Exercise 5.9's Green's function gives +1/(4π)
    # (default sign=+1); the printed sign is kept only as sign=−1 for the discriminating test.
    out *= float(sign) / _FOUR_PI  # u = +(1/4π) ∫ (∇′×ω)/|x − x′| d³x′   (5.14 with the sign corrected)
    return out[:, 0] if single else out


def biot_savart_volume(omega_fn: Callable, x, nodes, weights, eps: float = 0.0, t: float = 0.0,
                       chunk: int = 2_000_000):
    """Biot–Savart law (volume form), Eq. (5.16): u(x, t) = (1/4π)∫_V′ ω(x′, t) × (x − x′)/|x − x′|³ d³x′.

    Book: §5.5, Eq. (5.16), Fig. 5.8 — obtained from (5.14) by the product rule for the curl, Gauss' theorem in curl
    form (5.15) and a V′ whose ends are normal to ω and whose side lies outside the vortex (the surface term vanishes).
    Vorticity here sets the velocity everywhere, like an electric current sets a magnetic field.

    Parameters
    ----------
    omega_fn : ω(x′, t) at nodes (3, M) → (3, M) [1/s];  x : field point(s) (3,) or (3, N) [m]
    nodes, weights : quadrature of V′ [m], [m³];  eps : smoothing length ε [m] (kernel (|r|² + ε²)^{−3/2}; use ε > 0
    only when x lies inside the vorticity);  t : time [s];  chunk : memory control

    Returns
    -------
    u (3,) or (3, N) [m/s].

    Validation: V1 Gaussian/Rankine tube → u_θ(r) of (3.29)/(3.28); far field of a thin ring → (5.17); ∇·u = 0;
    V3 convergence in the node count. Label: analytic, converged.
    """
    X, single = _pts(x, 3)
    Pn, w = _F(nodes), _F(weights)
    Wm = _F(omega_fn(Pn, t)) * w  # ω dV′ (3, M)
    out = np.zeros_like(X)
    step = max(1, chunk // max(Pn.shape[1], 1))
    for i in range(0, X.shape[1], step):
        Xi = X[:, i:i + step]
        r = Xi[:, :, None] - Pn[:, None, :]  # x − x′ (3, n, M)
        k = (np.sum(r * r, axis=0) + eps ** 2) ** -1.5  # |x − x′|^{−3}
        cr = np.cross(Wm[:, None, :], r, axis=0)  # ω × (x − x′)
        out[:, i:i + step] = np.sum(cr * k, axis=2)
    out /= _FOUR_PI  # Eq. (5.16)
    return out[:, 0] if single else out


def biot_savart_2d(omega_z, X, Y, x, eps: float = 0.0, dA: float | None = None):
    """Plane Biot–Savart: velocity of a gridded vorticity ω_z(x′, y′), u = Σ ω dA e_z × (x − x′)/(2π|x − x′|²).

    Book: §5.5, (5.16) integrated along an infinite straight direction (the kernel of an infinite filament, (5.17) →
    Γ/2πr, our D13); §5.7 (the velocity of plane point vortices). Rectangle rule on a uniform grid (spectrally
    accurate for smooth vorticity that decays inside the grid); grid points that coincide with x are skipped.

    Parameters
    ----------
    omega_z : (ny, nx) [1/s] on ``X, Y`` = ``np.meshgrid(x, y, indexing="xy")`` [m]
    x : field point(s) (2,) or (2, N) [m];  eps : smoothing [m];  dA : cell area [m²] (default from the grid)

    Returns
    -------
    (2,) or (2, N) array (u, v) [m/s].

    Accuracy inside the vorticity: with eps = 0 the recovery is accurate only at dual-grid points (cell corners
    midway between nodes); at arbitrary interior points the singular kernel makes the error 1–7 % and it does not
    converge with the grid — use eps ≈ the grid spacing there. Outside the vorticity eps = 0 is spectrally accurate.

    Validation: V1 Gaussian vortex → u_θ(r) = (Γ/2πr)(1 − e^{−r²/σ²}) outside and inside the core; Rankine outside the
    core; V3 grid convergence. Label: analytic, converged.
    """
    P, single = _pts(x, 2)
    Xg, Yg, W = _F(X).ravel(), _F(Y).ravel(), _F(omega_z).ravel()
    if dA is None:
        dA = abs(float(_F(X)[0, 1] - _F(X)[0, 0]) * float(_F(Y)[1, 0] - _F(Y)[0, 0]))
    G = W * dA / (2.0 * np.pi)
    rx = P[0][:, None] - Xg[None, :]
    ry = P[1][:, None] - Yg[None, :]
    r2 = rx ** 2 + ry ** 2
    far = r2 > (1e-12 * np.sqrt(dA)) ** 2  # a node coinciding with x contributes nothing (self term)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(far, 1.0 / np.where(far, r2 + eps ** 2, 1.0), 0.0)
    u = -np.sum(G * ry * inv, axis=1)  # e_z × r = (−r_y, r_x)
    v = np.sum(G * rx * inv, axis=1)
    out = np.stack([u, v])
    return out[:, 0] if single else out


def velocity_from_vorticity_fft(omega, L):
    """Velocity of a periodic vorticity field by solving ∇²u = −∇ × ω spectrally (the Poisson equation of §5.5).

    Book: §5.5, "∇ × ω = ∇ × (∇ × u) = ∇(∇·u) − ∇²u = −∇²u" (identity (B.3.13) with ∇·u = 0). On a periodic box
    the Fourier transform turns it into −k²û = −i k × ω̂, so û = i k × ω̂/k² (the k = 0 mode — a uniform stream — is
    set to 0). In 2-D, ∇²ψ = −ω, ψ̂ = ω̂/k², u = ∂ψ/∂y, v = −∂ψ/∂x.

    Parameters
    ----------
    omega : 2-D ω_z (ny, nx) or 3-D (3, nz, ny, nx) [1/s] on a uniform periodic grid (x along the last axis)
    L : box size [m] (scalar or per-axis tuple in (x, y[, z]) order)

    Returns
    -------
    2-D: (2, ny, nx) array (u, v); 3-D: u (3, nz, ny, nx) [m/s].

    Validation: V1 Taylor–Green (2-D and 3-D) recovered to round-off; a Gaussian vortex in a large box → (3.29) with an
    error that falls with the box size (periodic images). Label: analytic.
    """
    w = _F(omega)
    if w.ndim == 2:
        ny, nx = w.shape
        Lx, Ly = (float(L), float(L)) if np.ndim(L) == 0 else (float(L[0]), float(L[1]))
        kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=Lx / nx)
        ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=Ly / ny)
        KX, KY = np.meshgrid(kx, ky, indexing="xy")
        k2 = KX ** 2 + KY ** 2
        k2[0, 0] = 1.0
        psi = np.fft.fft2(w) / k2  # ψ̂ = ω̂/k²  (∇²ψ = −ω)
        psi[0, 0] = 0.0
        u = np.real(np.fft.ifft2(1j * KY * psi))  # u = ∂ψ/∂y
        v = np.real(np.fft.ifft2(-1j * KX * psi))  # v = −∂ψ/∂x
        return np.stack([u, v])
    if w.ndim == 4 and w.shape[0] == 3:
        _, nz, ny, nx = w.shape
        Ls = (float(L),) * 3 if np.ndim(L) == 0 else tuple(float(v) for v in L)
        kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=Ls[0] / nx)
        ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=Ls[1] / ny)
        kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=Ls[2] / nz)
        KZ, KY, KX = np.meshgrid(kz, ky, kx, indexing="ij")
        K = np.stack([KX, KY, KZ])
        k2 = np.sum(K ** 2, axis=0)
        k2[0, 0, 0] = 1.0
        wh = np.fft.fftn(w, axes=(1, 2, 3))
        uh = 1j * np.cross(K, wh, axis=0) / k2  # û = i k × ω̂ / k²
        uh[:, 0, 0, 0] = 0.0
        return np.real(np.fft.ifftn(uh, axes=(1, 2, 3)))
    raise ValueError("omega must be (ny, nx) or (3, nz, ny, nx)")


# ======================================================================================================================
# §5.5 filaments, (5.17)
# ======================================================================================================================
def segment_induced_velocity(x, xa, xb, Gamma: float, eps: float = 0.0):
    """Velocity induced by straight vortex segment(s) from xa to xb of circulation Γ: (5.17) integrated in closed form.

    Book: §5.5, Eq. (5.17): du = (Γ dl/4π) e_ω × (x − x′)/|x − x′|³. Along a straight segment (our D13) the integral is
    (Γ/4πd)(cos θ_a − cos θ_b) along e_ω × e_d (d the distance from x to the segment's line, θ_a, θ_b the angles
    between e_ω and x − x_a, x − x_b); in vector form u = (Γ/4π) (r₁ × r₂)/|r₁ × r₂|² [r₀·(r₁/|r₁| − r₂/|r₂|)] with
    r₁ = x − x_a, r₂ = x − x_b, r₀ = x_b − x_a. An infinite line (θ_a → 0, θ_b → π) gives Γ/2πd — (5.2) recovered.

    Parameters
    ----------
    x : field point(s) (3,) or (3, N) [m];  xa, xb : segment ends (3,) or (3, M) [m] (M segments are summed)
    Gamma : circulation Γ [m²/s] along xa → xb (scalar or (M,));  eps : core cut-off ε [m] (|r₁ × r₂|² + ε²|r₀|²)

    Returns
    -------
    u (3,) or (3, N) [m/s]; 0 on the segment's line (the kernel is undefined there).

    Validation: V1 long segment → Γ/2πd; midpoint sum of (5.17) converges to it (order 2); V1 form cross-check with
    Wikipedia "Biot–Savart law" (aerodynamics). Label: analytic.
    """
    X, single = _pts(x, 3)
    A, B = _pts(xa, 3)[0], _pts(xb, 3)[0]
    G = np.broadcast_to(_F(Gamma), (A.shape[1],))
    r1 = X[:, :, None] - A[:, None, :]
    r2 = X[:, :, None] - B[:, None, :]
    r0 = (B - A)[:, None, :]
    c = np.cross(r1, r2, axis=0)
    c2 = np.sum(c * c, axis=0) + eps ** 2 * np.sum(r0 * r0, axis=0)
    n1, n2 = np.linalg.norm(r1, axis=0), np.linalg.norm(r2, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        dot = np.sum(r0 * (r1 / n1 - r2 / n2), axis=0)
        fac = np.where(c2 > 1e-300, G[None, :] / _FOUR_PI * dot / np.where(c2 > 1e-300, c2, 1.0), 0.0)
    u = np.sum(c * fac, axis=2)  # (Γ/4π)(r₁×r₂)/|r₁×r₂|² r₀·(r₁/|r₁| − r₂/|r₂|)
    return u[:, 0] if single else u


def segment_speed(d, theta_a, theta_b, Gamma: float):
    """Speed induced by a straight vortex segment, (Γ/4πd)(cos θ_a − cos θ_b) [m/s] (scalar form of our D13).

    Book: §5.5, Eq. (5.17) integrated along a straight segment; θ_a = 0, θ_b = π → Γ/2πd, the ideal line vortex (5.2).
    Parameters: d distance to the segment's line [m]; theta_a, theta_b [rad]; Gamma [m²/s]. Label: analytic.
    """
    return _S(float(Gamma) / (_FOUR_PI * _F(d)) * (np.cos(_F(theta_a)) - np.cos(_F(theta_b))))


def filament_velocity(x, polyline, Gamma: float, closed: bool = True, eps: float = 0.0):
    """Velocity induced by a polygonal vortex filament (a sum of straight segments, (5.17)).

    Book: §5.5, Eq. (5.17) summed along the filament (horseshoe vortices and vortex rings in Ch. 14 are built so).
    Parameters: x (3,) or (3, N) [m]; polyline (3, M) vertices [m] ordered along e_ω; Gamma [m²/s]; closed (join the
    last vertex to the first); eps [m]. Returns u [m/s].
    Validation: V1 regular M-gon ring → ΓR²/2(R² + z²)^{3/2} on the axis with error O(1/M²); square loop centre
    2√2Γ/(πL). Label: analytic, converged.
    """
    P = _F(polyline)
    A = P
    B = np.roll(P, -1, axis=1)
    if not closed:
        A, B = P[:, :-1], P[:, 1:]
    return segment_induced_velocity(x, A, B, Gamma, eps)


def ring_axis_velocity(z, R: float, Gamma: float):
    """Axial velocity on the axis of a circular vortex ring, u_z = ΓR²/[2(R² + z²)^{3/2}] [m/s] (along +z for Γ > 0
    about e_φ). Book: §5.5, Eq. (5.17) integrated round a circle (our worked result; the magnetic-loop analogue).
    Parameters: z axial distance from the ring plane [m]; R ring radius [m]; Gamma [m²/s]. Label: analytic."""
    return _S(float(Gamma) * float(R) ** 2 / (2.0 * (float(R) ** 2 + _F(z) ** 2) ** 1.5))


def ring_ring_velocity(R, z, R0: float, z0: float, Gamma0: float):
    """Velocity (u_R, u_z) at the point (R, z) induced by a coaxial circular vortex filament of radius R₀ in the plane
    z = z₀ with circulation Γ₀ (vorticity along +e_φ), from its Stokes stream function (complete elliptic integrals).

    Book: §5.7, Fig. 5.15 (a vortex ring and its image ring near a wall; leap-frogging rings): (5.17) integrated round
    the circle. Standard closed form (Lamb, Hydrodynamics §161), with Δz = z − z₀,
    D₊ = √((R + R₀)² + Δz²), D₋² = (R − R₀)² + Δz², k² = 4RR₀/D₊²:
    u_z = (Γ₀/2πD₊)[K(k) + (R₀² − R² − Δz²)/D₋² · E(k)],
    u_R = (Γ₀Δz/2πRD₊)[−K(k) + (R₀² + R² + Δz²)/D₋² · E(k)].
    ``scipy.special.ellipk/ellipe`` take the **parameter m = k²** (not the modulus k).

    Parameters: R, z [m] (arrays broadcast); R0, z0 [m]; Gamma0 [m²/s]. Returns (u_R, u_z) [m/s]; on the axis
    u_R = 0 and u_z = ring_axis_velocity; singular on the filament itself.
    Validation: V1 on the axis → ΓR₀²/2(R₀² + Δz²)^{3/2}; V3 equals the polygon Biot–Savart (:func:`filament_velocity`)
    with error O(1/M²). Label: analytic, converged.
    """
    R_, z_ = np.broadcast_arrays(_F(R), _F(z))
    dz = z_ - float(z0)
    Dp2 = (R_ + R0) ** 2 + dz ** 2
    Dm2 = (R_ - R0) ** 2 + dz ** 2
    m = 4.0 * R_ * R0 / Dp2  # parameter m = k²
    K, E = ellipk(m), ellipe(m)
    Dp = np.sqrt(Dp2)
    with np.errstate(divide="ignore", invalid="ignore"):
        uz = Gamma0 / (2.0 * np.pi * Dp) * (K + (R0 ** 2 - R_ ** 2 - dz ** 2) / Dm2 * E)
        uR = np.where(R_ > 1e-12 * R0,
                      Gamma0 * dz / (2.0 * np.pi * np.where(R_ > 0, R_, 1.0) * Dp)
                      * (-K + (R0 ** 2 + R_ ** 2 + dz ** 2) / Dm2 * E), 0.0)
    return _S(uR), _S(uz)


RING_CORES = {"uniform": 0.25, "hollow": 0.5, "gaussian": 0.558}


def ring_self_velocity(R, a, Gamma: float, core: str = "uniform"):
    """Self-induced speed of a thin vortex ring, U = (Γ/4πR)[ln(8R/a) − C] (Kelvin; C = ¼ for a uniform core).

    Book: §5.7 (rings move by their own induced velocity; the book treats rings qualitatively). The straight-filament
    law has no self-induction, a curved one has a logarithmically singular one — a finite core radius a ≪ R is needed.
    ``core``: "uniform" (C = 1/4, Kelvin; a = the core radius), "hollow" (C = 1/2), "gaussian" (C ≈ 0.558, Saffman;
    a is the e-folding radius of the core vorticity ∝ e^{−r²/a²}, a = √(4νt) for a diffusing core).
    Parameters: R ring radius [m]; a core radius [m]; Gamma [m²/s]. Returns U [m/s] along e_z for Γ > 0.
    Validation: V1 form cross-check (published formula): Kelvin's ring speed, Wikipedia "Vortex ring". Label:
    analytic.
    """
    C = RING_CORES[core]
    R_ = _F(R)
    return _S(float(Gamma) / (_FOUR_PI * R_) * (np.log(8.0 * R_ / _F(a)) - C))


# ======================================================================================================================
# §5.7 point vortices and images
# ======================================================================================================================
def point_vortex_velocity(x, xv, Gamma, eps: float = 0.0):
    """Velocity at x induced by plane point vortices (ideal line vortices) at xv with circulations Γ:
    u = Σ_j (Γ_j/2π) e_z × (x − x_j)/(|x − x_j|² + ε²).

    Book: §5.7 (ideal line vortices superpose; each is u_θ = Γ/2πr, (5.2)); a vortex does not move itself — a point
    closer than 1e-14 m to a vortex receives nothing from it (the self term is excluded).

    Parameters: x (2,) or (2, N) [m]; xv (2, M) [m]; Gamma (M,) [m²/s] (counterclockwise +); eps smoothing [m].
    Returns (2,) or (2, N) [m/s].
    Validation: V1 single vortex → Γ/2πr tangential; C12's double loop. Label: analytic.
    """
    X, single = _pts(x, 2)
    V = _F(xv).reshape(2, -1)
    G = np.broadcast_to(_F(Gamma), (V.shape[1],))
    rx = X[0][:, None] - V[0][None, :]
    ry = X[1][:, None] - V[1][None, :]
    r2 = rx ** 2 + ry ** 2
    far = r2 > 1e-28
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(far, G[None, :] / (2.0 * np.pi) / np.where(far, r2 + eps ** 2, 1.0), 0.0)
    u = np.stack([-np.sum(f * ry, axis=1), np.sum(f * rx, axis=1)])  # (Γ/2π) e_z × r / r²
    return u[:, 0] if single else u


def wall_image_system(xv, Gamma, wall_y: float = 0.0):
    """Vortices plus their mirror images in the plane wall y = wall_y (strength −Γ at the mirror point), Fig. 5.14.

    Book: §5.7, method of images: a vortex of equal and opposite strength at the image point B; the wall becomes a
    streamline (zero normal velocity). Returns (xv_all (2, 2M), Gamma_all (2M,)) — originals first.
    Validation: V1 normal velocity on the wall ≡ 0; drift Γ/4πh; wall speed Γh/π(x² + h²). Label: analytic.
    """
    V = _F(xv).reshape(2, -1)
    G = np.broadcast_to(_F(Gamma), (V.shape[1],)).astype(float)
    img = np.stack([V[0], 2.0 * wall_y - V[1]])
    return np.concatenate([V, img], axis=1), np.concatenate([G, -G])


def circle_image_system(xv, Gamma, a: float, center=(0.0, 0.0), inside: bool = True,
                        cylinder_circulation: float = 0.0):
    """Vortices plus images that make the circle |x − c| = a a streamline (Milne-Thomson circle theorem for vortices).

    Each vortex Γ at x gets an image −Γ at the inverse point c + a²(x − c)/|x − c|² (a vortex exactly at the centre
    has its image at infinity: none is added). ``inside=True`` (the fluid is
    inside the circle — the bucket of Fig. 5.13): the inverse-point images alone (no centre vortex is allowed there).
    ``inside=False`` (fluid outside a cylinder, Exercise 5.14, Ch. 6): a vortex ΣΓ + ``cylinder_circulation`` is added
    at the centre so that the circulation round the cylinder is ``cylinder_circulation`` (zero by default).
    Returns (xv_all, Gamma_all) — originals first.

    Book: §5.7 (Fig. 5.13) and Exercise 5.14. Validation: V1 velocity normal to the circle ≡ 0 at 400 points; zero net
    circulation round the cylinder. Label: analytic.
    """
    V = _F(xv).reshape(2, -1)
    G = np.broadcast_to(_F(Gamma), (V.shape[1],)).astype(float)
    c = _F(center).reshape(2, 1)
    d = V - c
    r2 = np.sum(d * d, axis=0)
    fin = r2 > 0.0  # a vortex exactly at the centre has its image at infinity: drop it
    img = c + a ** 2 * d[:, fin] / r2[fin]
    pos, gam = [V, img], [G, -G[fin]]
    if not inside:
        pos.append(c)
        gam.append(np.array([np.sum(G) + float(cylinder_circulation)]))
    return np.concatenate(pos, axis=1), np.concatenate(gam)


def _boundary_spec(boundary, bp):
    """Normalise boundary = None | "wall" | "circle" | "channel" | tuple forms → (kind, params)."""
    if boundary is None:
        return None, {}
    if isinstance(boundary, str):
        return boundary, dict(bp)
    b = tuple(boundary)
    if b[0] == "wall":
        return "wall", dict(wall_y=float(b[1]) if len(b) > 1 else 0.0, **bp)
    if b[0] == "circle":
        return "circle", dict(a=float(b[1]), center=b[2] if len(b) > 2 else (0.0, 0.0), **bp)
    if b[0] == "channel":
        return "channel", dict(H=float(b[1]), **bp)
    raise ValueError(f"unknown boundary {boundary!r}")


def image_vortices(xv, Gamma, boundary, **bp):
    """The image vortices only (positions (2, K), strengths (K,)) for ``boundary`` None, "wall" (``wall_y``),
    "circle" (``a``, ``center``, ``inside``) (tuples ("wall", y), ("circle", a[, center]) also accepted).
    Book: §5.7. Label: analytic."""
    V = _F(xv).reshape(2, -1)
    M = V.shape[1]
    kind, q = _boundary_spec(boundary, bp)
    if kind is None:
        return np.zeros((2, 0)), np.zeros(0)
    if kind == "wall":
        P, G = wall_image_system(V, Gamma, float(q.get("wall_y", 0.0)))
    elif kind == "circle":
        P, G = circle_image_system(V, Gamma, float(q.get("a", 1.0)), q.get("center", (0.0, 0.0)),
                                   bool(q.get("inside", True)), float(q.get("cylinder_circulation", 0.0)))
    else:
        raise ValueError('image_vortices supports "wall" and "circle" (the channel uses a closed-form sum)')
    return P[:, M:], G[M:]


def _channel_velocity(X, V, G, H, self_index=None):
    """Velocity at points X (2, N) of vortices V (2, M) between walls y = 0 and y = H with all their images, in closed
    form: the column of copies at y_j + 2nH (+Γ) and at −y_j + 2nH (−Γ) gives u − iv = (Γ/2πi)(π/2H)coth(π(z − z_j)/2H)
    − (Γ/2πi)(π/2H)coth(π(z − z̄_j)/2H). For X = V (``self_index``) the singular 1/(z − z_k) part of a vortex's own
    column is removed (a vortex does not move itself)."""
    z = X[0] + 1j * X[1]
    zj = V[0] + 1j * V[1]
    k = np.pi / (2.0 * H)
    dz = z[:, None] - zj[None, :]
    dzb = z[:, None] - np.conj(zj)[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        c1 = k / np.tanh(k * dz)
        if self_index is not None:
            n = X.shape[1]
            idx = np.arange(n)
            c1[idx, self_index] = 0.0  # lim_{ζ→0} [k coth(kζ) − 1/ζ] = 0
        c2 = k / np.tanh(k * dzb)
    w = np.sum(G[None, :] / (2j * np.pi) * (c1 - c2), axis=1)  # u − iv
    return np.stack([w.real, -w.imag])


def point_vortex_rhs(xv, Gamma, boundary=None, eps: float = 0.0, **bp):
    """Velocities (2, M) of the point vortices themselves: each is moved by all the others and by the images of the
    boundary (self term excluded) — the right side of dx_k/dt in :func:`point_vortex_evolve`.
    ``boundary``: None, "wall" (``wall_y``=0.0), "circle" (``a``, ``center``, ``inside``=True), "channel" (``H``).
    Book: §5.7. Label: analytic."""
    V = _F(xv).reshape(2, -1)
    G = np.broadcast_to(_F(Gamma), (V.shape[1],)).astype(float)
    kind, q = _boundary_spec(boundary, bp)
    if kind == "channel":
        return _channel_velocity(V, V, G, float(q.get("H", 1.0)), self_index=np.arange(V.shape[1]))
    u = point_vortex_velocity(V, V, G, eps)  # self term excluded (r = 0 → 0)
    if kind is not None:
        Pi, Gi = image_vortices(V, G, kind, **q)
        u = u + point_vortex_velocity(V, Pi, Gi, eps)
    return u


def point_vortex_evolve(xv0, Gamma, t_eval, boundary=None, rtol: float = 1e-11, atol: float = 1e-13,
                        eps: float = 0.0, method: str = "DOP853", **bp):
    """Evolve plane point vortices, each carried by the velocity of all the others (and of the boundary's images):
    dx_k/dt = Σ_{j≠k} (Γ_j/2π) e_z × (x_k − x_j)/|x_k − x_j|² + images.

    Book: §5.7: vortex lines move with the flow (Helmholtz 1) — each ideal line vortex moves with the velocity induced
    at its position by the others; same-sign pairs orbit their centre of vorticity (Fig. 5.11), opposite pairs
    translate (Fig. 5.12), a vortex near a wall slides along it (Fig. 5.14), the knife-blade pair in a bucket separates
    along the wall (Fig. 5.13, circle images).

    Parameters
    ----------
    xv0 : (2, M) initial positions [m];  Gamma : (M,) [m²/s];  t_eval : output times (T,) [s] (first = start)
    boundary : None, "wall" (``wall_y``=0.0), "circle" (``a``, ``center``=(0, 0), ``inside``=True), "channel" (``H``)
    or the tuple forms ("wall", y), ("circle", a[, center]), ("channel", H);  bp : the boundary keywords
    rtol, atol, method : ``solve_ivp`` (DOP853, tight tolerances);  eps : smoothing [m]

    Returns
    -------
    (T, 2, M) positions [m].

    Validation: V4 ΣΓx, ΣΓ|x|² and the Hamiltonian conserved (no boundary) to 1e-10; V1 pair period
    2π/((Γ₁ + Γ₂)/2πh²); V1 form cross-check (published formula): García & Haziot (2023) pair rates.
    Label: conserved, analytic.
    """
    V0 = _F(xv0).reshape(2, -1)
    M = V0.shape[1]
    G = np.broadcast_to(_F(Gamma), (M,)).astype(float)
    te = np.atleast_1d(_F(t_eval))

    def rhs(t, y):
        return point_vortex_rhs(y.reshape(2, M), G, boundary, eps, **bp).ravel()

    if te.size == 1:
        return V0[None]
    sol = solve_ivp(rhs, (te[0], te[-1]), V0.ravel(), method=method, rtol=rtol, atol=atol, t_eval=te)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T.reshape(te.size, 2, M)


def centre_of_vorticity(xv, Gamma) -> np.ndarray:
    """Centre of vorticity ΣΓ_k x_k/ΣΓ_k (2,) [m] (NaN if ΣΓ = 0) — the point G of Fig. 5.11 about which a same-sign
    pair orbits (Exercise 5.18). Book: §5.7. Label: conserved."""
    V = _F(xv).reshape(2, -1)
    G = _F(Gamma)
    tot = float(np.sum(G))
    return V @ G / tot if tot != 0.0 else np.full(2, np.nan)


def point_vortex_invariants(xv, Gamma) -> dict:
    """Invariants of unbounded plane point-vortex motion: linear impulse P = ΣΓ_k x_k (2,) [m³/s], angular impulse
    I = ΣΓ_k|x_k|² [m⁴/s] and the Kirchhoff energy H = −(1/2π)Σ_{j<k}Γ_jΓ_k ln|x_j − x_k| [m⁴/s²]. Accepts (2, M)
    (floats) or (T, 2, M) (arrays over T).

    Book: §5.7 (the centre of vorticity G of Fig. 5.11 stays fixed: P is conserved). Label: conserved.
    """
    V = _F(xv)
    G = _F(Gamma)
    single = V.ndim == 2
    if single:
        V = V[None]
    P = np.einsum("m,tim->ti", G, V)
    I = np.einsum("m,tm->t", G, np.sum(V ** 2, axis=1))
    M = V.shape[2]
    H = np.zeros(V.shape[0])
    for i in range(M):
        for j in range(i + 1, M):
            H -= G[i] * G[j] * np.log(np.linalg.norm(V[:, :, i] - V[:, :, j], axis=1)) / (2.0 * np.pi)
    if single:
        return dict(P=P[0], I=float(I[0]), H=float(H[0]))
    return dict(P=P, I=I, H=H)


def channel_image_velocity_exact(h, H: float, Gamma: float):
    """Closed form of the channel-image series: u = (Γ/4H) cot(πh/H) [m/s] (our partial-fraction sum of Exercise
    5.13's series, using Σ_{n≥1} 1/(n² − a²) = 1/(2a²) − π cot(πa)/(2a)); 0 at mid-channel.
    Book: Exercise 5.13 (§5.7 method of images). Label: analytic."""
    return _S(float(Gamma) / (4.0 * float(H)) / np.tan(np.pi * _F(h) / float(H)))


def channel_image_velocity(h: float, H: float, Gamma: float, tol: float = 1e-12, n_max: int = 10 ** 6) -> float:
    """Speed of a line vortex between parallel walls y = 0 and y = H at height h, from the infinite image series.

    Book: Exercise 5.13 (method of images, §5.7): u(0, h) = (Γ/4πh)[1 − 2Σ_{n≥1} 1/((nH/h)² − 1)]. Summed in blocks
    until the increment (the last term) is below ``tol`` or n reaches ``n_max`` (never a fixed count), then the tail
    beyond n = N is added from the Euler–Maclaurin integral ∫_N^∞ dn/((cn)² − 1) − ½f(N) (c = H/h). The closed form
    (Γ/4H)cot(πh/H) is :func:`channel_image_velocity_exact`. Parameters: h (0 < h < H) [m]; H [m]; Gamma [m²/s].
    Returns u [m/s] (along +x for Γ > 0 nearer the lower wall).
    Validation: V1 h = H/2 → 0 (the sum is ½); h = H/4 → Γ/4H; equals a direct sum over ±20 000 image pairs.
    Label: analytic, converged.
    """
    c = float(H) / float(h)
    total, n0, block = 0.0, 1, 10_000
    while True:
        n = np.arange(n0, min(n0 + block, n_max + 1), dtype=float)
        f = 1.0 / ((c * n) ** 2 - 1.0)
        total += float(np.sum(f))
        n0 = int(n[-1]) + 1
        if f[-1] < tol or n0 > n_max:
            break
        block = min(block * 2, 1_000_000)
    N = n0 - 1
    fN = 1.0 / ((c * N) ** 2 - 1.0)
    total += np.log((c * N + 1.0) / (c * N - 1.0)) / (2.0 * c) - 0.5 * fN  # tail ∫_N^∞ f dn − ½ f(N)
    return float(Gamma) / (4.0 * np.pi * float(h)) * (1.0 - 2.0 * total)  # Exercise 5.13 series


def vortex_sheet_velocity(x, sheet_points, gamma, ds, eps: float = 0.0):
    """Velocity of a vortex sheet represented by point vortices of strength γ ds at ``sheet_points``.

    Book: §5.8, Fig. 5.16 (infinitely many ideal line vortices side by side; the tangential velocity jumps by the
    sheet strength across it, ±γ/2 either side of a long straight sheet, the normal velocity does not).
    Parameters: x (2,) or (2, N) [m]; sheet_points (2, M) [m]; gamma γ [m/s] (scalar or (M,); counterclockwise +);
    ds spacing [m] (scalar or (M,)); eps [m]. Returns (2,) or (2, N) [m/s].
    Label: analytic, converged.
    """
    S = _F(sheet_points).reshape(2, -1)
    Gm = np.broadcast_to(_F(gamma) * _F(ds), (S.shape[1],))
    return point_vortex_velocity(x, S, Gm, eps)


def continuous_sheet_velocity(x, y, gamma: float, L: float = 1.0):
    """Velocity (u, v) of a flat, finite, continuous vortex sheet of strength γ on −L/2 < x′ < L/2, y′ = 0:
    u = −(γ/2π)[arctan((L/2 − x)/y) + arctan((L/2 + x)/y)], v = (γ/4π) ln(((x + L/2)² + y²)/((x − L/2)² + y²)).

    Book: §5.8, Fig. 5.16 (the continuous limit of the row of filaments): u → ∓γ/2 just above/below a long sheet, so
    the jump u_below − u_above = γ (text: u₂ − u₁). Our integral of the plane kernel along the sheet; the reference for
    :func:`discrete_sheet_convergence`. Parameters: x, y [m] (y ≠ 0); gamma [m/s]; L [m]. Returns (u, v) [m/s].
    Label: analytic.
    """
    x_, y_ = _F(x), _F(y)
    with np.errstate(divide="ignore", invalid="ignore"):
        u = -(gamma / (2.0 * np.pi)) * (np.arctan((0.5 * L - x_) / y_) + np.arctan((0.5 * L + x_) / y_))
        u = np.where(y_ == 0.0, 0.0, u)
        v = (gamma / (4.0 * np.pi)) * np.log(((x_ + 0.5 * L) ** 2 + y_ ** 2) / ((x_ - 0.5 * L) ** 2 + y_ ** 2))
    return _S(u), _S(v)


# ======================================================================================================================
# Explainer presets (E6 filaments, E8 point vortices)
# ======================================================================================================================
FILAMENT_PRESETS = ("straight", "semi_infinite", "square", "ring", "helix", "bent")
FILAMENT_CLOSED = {"straight": False, "semi_infinite": False, "square": True, "ring": True, "helix": False,
                   "bent": False}


def filament_preset(name: str, M: int = 64, **p) -> np.ndarray:
    """Filament shapes for explainer E6 (``biot_savart_filament``): vertices (3, M + 1) of an open polyline or (3, M)
    of a closed one, ordered along e_ω (closedness: ``FILAMENT_CLOSED[name]``).

    * "straight" — from (0, 0, −ℓ) to (0, 0, ℓ), ``half_length`` ℓ = 1 m (open);
    * "semi_infinite" — from the origin up the z-axis to ``length`` = 1000 m (open, geometric spacing);
    * "square" — ``side`` = 1 m in z = 0, centred, counterclockwise about +z, M points on the perimeter (closed);
    * "ring" — radius ``R`` = 0.5 m in z = 0, counterclockwise about +z (closed);
    * "helix" — radius ``R`` = 0.3 m, ``pitch`` = 0.4 m per turn, ``turns`` = 3, centred on z = 0 (open);
    * "bent" — an L-shaped filament of arm ``half_length`` (open).
    Book: §5.5, Fig. 5.8 and Eq. (5.17). Label: analytic.
    """
    if name == "straight":
        ell = float(p.get("half_length", 1.0))
        z = np.linspace(-ell, ell, M + 1)
        return np.stack([0 * z, 0 * z, z])
    if name == "semi_infinite":
        Lf = float(p.get("length", 1000.0))
        z = np.concatenate([[0.0], np.geomspace(1e-4 * Lf, Lf, M)])
        return np.stack([0 * z, 0 * z, z])
    if name == "square":
        h = 0.5 * float(p.get("side", 1.0))
        m = max(M // 4, 1)
        s = -h + 2.0 * h * np.arange(m) / m
        sides = [np.stack([s, -h + 0 * s]), np.stack([h + 0 * s, s]), np.stack([-s, h + 0 * s]),
                 np.stack([-h + 0 * s, -s])]
        xy = np.concatenate(sides, axis=1)
        return np.vstack([xy, np.zeros((1, xy.shape[1]))])
    if name == "ring":
        R = float(p.get("R", 0.5))
        th = 2.0 * np.pi * np.arange(M) / M
        return np.stack([R * np.cos(th), R * np.sin(th), 0 * th])
    if name == "helix":
        R, pitch, turns = float(p.get("R", 0.3)), float(p.get("pitch", 0.4)), float(p.get("turns", 3.0))
        th = np.linspace(0.0, 2.0 * np.pi * turns, M + 1)
        return np.stack([R * np.cos(th), R * np.sin(th), pitch * (th - th[-1] / 2) / (2.0 * np.pi)])
    if name == "bent":
        ell = float(p.get("half_length", 1.0))
        s = np.linspace(0.0, ell, M // 2 + 1)
        return np.concatenate([np.stack([0 * s, 0 * s, -ell + s]), np.stack([s[1:], 0 * s[1:], 0 * s[1:]])], axis=1)
    raise ValueError(f"unknown filament {name!r}; choose from {FILAMENT_PRESETS}")


def filament_contributions(x, polyline, Gamma: float, closed: bool = True, eps: float = 0.0) -> np.ndarray:
    """Velocity induced at one point x (3,) by **each** straight segment of a polyline, (3, M) [m/s] (one column per
    segment; their sum is :func:`filament_velocity`) — for E6's per-segment arrows and inspector.
    Book: §5.5, Eq. (5.17). Label: analytic."""
    P = _F(polyline)
    A = P if closed else P[:, :-1]
    B = np.roll(P, -1, axis=1) if closed else P[:, 1:]
    return np.stack([segment_induced_velocity(x, A[:, j], B[:, j], Gamma, eps) for j in range(A.shape[1])], axis=1)


def filament_velocity_preset(name: str, x: float, y: float, z: float, Gamma: float = 1.0, M: int = 64,
                             component: int | None = None, **p):
    """Velocity at (x, y, z) induced by the preset filament ``name`` of :func:`filament_preset` with circulation Γ
    (a float for ``component`` 0/1/2, else (3,)) [m/s] — the parity-friendly wrapper for E6.
    Book: §5.5, Eq. (5.17). Label: analytic."""
    P = filament_preset(name, M, **p)
    u = filament_velocity(np.array([x, y, z], float), P, Gamma, FILAMENT_CLOSED[name])
    return float(u[component]) if component is not None else u


POINT_VORTEX_PRESETS = ("equal_pair", "unequal_pair", "opposite_pair", "near_wall", "knife_bucket", "three_vortices")
_PV_ALIASES = {"wall": "near_wall", "bucket": "knife_bucket", "three": "three_vortices"}


def point_vortex_preset(name: str) -> dict:
    """Configurations for explainer E8 (``point_vortex_lab``) and the C12/C13 animations (SI units):

    * "equal_pair" — Γ = (1, 1) m²/s at (−½, 0), (½, 0): orbit the origin at (Γ₁ + Γ₂)/2πh² (Fig. 5.11);
    * "unequal_pair" — Γ = (1, 3) at the same places: G at h₁ = Γ₂h/(Γ₁ + Γ₂) = 0.75 m from vortex 1;
    * "opposite_pair" — Γ = (1, −1) at (−½, 0), (½, 0): both move in +y at Γ/2πh = 0.159155 m/s (Fig. 5.12);
    * "near_wall" — Γ = 1 at (0, 0.5) above the wall y = 0: drifts at Γ/4πh = 0.159155 m/s (Fig. 5.14);
    * "knife_bucket" — Γ = (−1, 1) at (0, 0.15), (0, −0.15) inside a circle a = 1 about (0.3, 0): the pair moves in −x
      toward the wall and separates along it (Fig. 5.13, interior circle images);
    * "three_vortices" — Γ = (1, 1, −0.5) at (0, 0), (1, 0), (0.3, 0.8) (invariants test).

    Returns dict(xv (2, M), Gamma (M,), boundary (None, "wall", "circle"), bp (boundary keywords), t_end [s], label).
    Book: §5.7. Label: analytic.
    """
    name = _PV_ALIASES.get(name, name)
    if name == "equal_pair":
        return dict(xv=np.array([[-0.5, 0.5], [0.0, 0.0]]), Gamma=np.array([1.0, 1.0]), boundary=None, bp={},
                    t_end=2.0 * 19.7392088, label="same-sign pair: orbits its centre of vorticity")
    if name == "unequal_pair":
        return dict(xv=np.array([[-0.5, 0.5], [0.0, 0.0]]), Gamma=np.array([1.0, 3.0]), boundary=None, bp={},
                    t_end=2.0 * 2.0 * np.pi / (4.0 / (2.0 * np.pi)), label="Γ₂ = 3Γ₁: G sits nearer the stronger")
    if name == "opposite_pair":
        return dict(xv=np.array([[-0.5, 0.5], [0.0, 0.0]]), Gamma=np.array([1.0, -1.0]), boundary=None, bp={},
                    t_end=20.0, label="opposite pair: translates at Γ/2πh")
    if name == "near_wall":
        return dict(xv=np.array([[0.0], [0.5]]), Gamma=np.array([1.0]), boundary="wall", bp=dict(wall_y=0.0),
                    t_end=20.0, label="vortex near a wall: slides along it at Γ/4πh")
    if name == "knife_bucket":
        return dict(xv=np.array([[0.0, 0.0], [0.15, -0.15]]), Gamma=np.array([-1.0, 1.0]), boundary="circle",
                    bp=dict(a=1.0, center=(0.3, 0.0), inside=True), t_end=12.0,
                    label="knife-blade pair in a bucket: separates along the wall")
    if name == "three_vortices":
        return dict(xv=np.array([[0.0, 1.0, 0.3], [0.0, 0.0, 0.8]]), Gamma=np.array([1.0, 1.0, -0.5]), boundary=None,
                    bp={}, t_end=20.0, label="three vortices: the invariants stay flat")
    raise ValueError(f"unknown preset {name!r}; choose from {POINT_VORTEX_PRESETS}")
