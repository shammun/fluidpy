"""Vertical normal modes of a stratified layer: the Sturm–Liouville problem for N²(z) with a free surface or a rigid lid,
equivalent depths, projection and reconstruction, and the exact modes for uniform N.

Book: Kundu, Cohen & Dowling 5e, Ch. 13 §13.9, Eqs. (13.52)–(13.71) (rendered pages chapters/pages/ch13/p671–p675).

The problem
-----------
With [u, v, p/ρ₀] = Σ [u_n, v_n, p_n](x, y, t) ψ_n(z) (13.52), the vertical structure obeys

    d/dz ( (1/N²) dψ_n/dz ) + ψ_n / c_n² = 0                                   (13.56)
    dψ_n/dz = 0                    at the bottom z = −H                          (13.64)
    dψ_n/dz + (N²/g) ψ_n = 0       at the free surface z = 0                     (13.65)
    (rigid lid: dψ_n/dz = 0 at z = 0 instead)

and each mode then moves horizontally like a shallow layer of equivalent depth H_e = c_n²/g (13.62).

Conventions
-----------
* z is measured upward from the free surface: the nodes run from z[0] = −H (bottom) to z[-1] = 0 (surface).
* :class:`Modes` lists the modes by decreasing speed.  With a free surface index 0 is the barotropic mode (the book's
  n = 0) and index n the n-th baroclinic mode.  **With a rigid lid there is no barotropic mode** (its speed is
  infinite, ψ₀ = 1): index 0 is then the first baroclinic mode, index j the book's n = j + 1.
* ψ_n is dimensionless, normalised to ψ_n(0) = 1 (positive at the surface).  Then u_n, v_n are velocities and p_n is
  p/ρ₀ [m²/s²], but w_n multiplies ∫ψ dz (so w_n is in 1/s) and ρ_n multiplies dψ/dz (so ρ_n is in kg/m²) — trap T12.
* **Orthogonality.**  Multiplying (13.56) for mode n by ψ_m, integrating by parts and using (13.64), (13.65) gives
  (1/c_n²) ∫ψ_nψ_m dz = ∫ ψ_n′ψ_m′/N² dz + ψ_n(0)ψ_m(0)/g.  The right-hand side is symmetric in m and n, so
  (1/c_n² − 1/c_m²) ∫ψ_nψ_m dz = 0: **the ψ_n are orthogonal with weight 1 with a free surface as well as with a
  rigid lid**.  The surface term belongs to the *second* ("energy") relation,
  ∫ ψ_n′ψ_m′/N² dz + ψ_n(0)ψ_m(0)/g = 0 for m ≠ n.
* The numerical methods (finite volumes, Chebyshev collocation, shooting) are ours — the book solves only uniform N.
"""
from __future__ import annotations

from typing import Callable, NamedTuple

import numpy as np
from scipy.linalg import eig, eigh_tridiagonal
from scipy.optimize import brentq

from ._util import as_scalar_if_0d
from .stability import cheb
from .thermo import G0

__all__ = ["Modes", "vertical_modes", "vertical_modes_shooting", "modes_uniform_N", "rigid_lid_error", "w_structure",
           "rho_structure", "project", "reconstruct", "orthogonality_matrix", "modal_amplitudes", "wkb_mode_speed"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


class Modes(NamedTuple):
    """A set of vertical modes (a named tuple, so ``modes[0][1]`` is the second speed).

    Fields
    ------
    c : (n_modes,) long-wave speeds [m/s], decreasing.
    psi : (n_modes, nz) mode shapes at the nodes (dimensionless, ψ(z = 0) = 1).
    He : (n_modes,) equivalent depths c²/g [m].
    z : (nz,) nodes [m], from −H to 0.
    weights : (nz,) trapezoid weights [m] — Σ weights·ψ_m·ψ_n approximates ∫ψ_mψ_n dz.
    """

    c: np.ndarray
    psi: np.ndarray
    He: np.ndarray
    z: np.ndarray
    weights: np.ndarray


def _trapz_weights(z: np.ndarray) -> np.ndarray:
    w = np.zeros_like(z)
    h = np.diff(z)
    w[:-1] += 0.5 * h
    w[1:] += 0.5 * h
    return w


def _d1(z: np.ndarray, y: np.ndarray) -> np.ndarray:
    from .gfd import _d1 as d1

    return d1(z, y)


def _check_grid(z, N2):
    z = _F(z)
    if z.ndim != 1 or z.size < 5 or np.any(np.diff(z) <= 0):
        raise ValueError("z must be a strictly increasing 1-D array (bottom first, surface last) with at least 5 nodes")
    N2a = _F(N2(z)) * np.ones_like(z) if callable(N2) else _F(N2) * np.ones_like(z)
    if N2a.shape != z.shape:
        raise ValueError("N2 must be a scalar, a callable or one value per node")
    if np.any(N2a <= 0) or not np.all(np.isfinite(N2a)):
        raise ValueError("N2 must be positive and finite at every node (a statically stable profile)")
    return z, N2a


def _normalise(psi: np.ndarray) -> np.ndarray:
    """Scale each row to ψ(0) = 1 (surface value, last node)."""
    top = psi[:, -1]
    if np.any(np.abs(top) < 1e-14 * np.max(np.abs(psi), axis=1)):
        raise ValueError("a mode vanishes at the surface and cannot be normalised to psi(0) = 1")
    return psi / top[:, None]


def vertical_modes(z, N2, g: float = G0, n_modes: int = 4, lid: str = "free", method: str = "fd",
                   n_cheb: int = 48) -> Modes:
    """Vertical normal modes ψ_n(z) and speeds c_n of a stratification N²(z).

    Book: §13.9, Eq. (13.56) with the boundary conditions (13.64) and (13.65); H_e from Eq. (13.62).
    Parameters
    ----------
    z : nodes [m], strictly increasing from the bottom z[0] = −H to the surface z[-1] = 0 (non-uniform allowed).
    N2 : buoyancy frequency squared [1/s²] > 0 — one value per node, a scalar, or a callable N2(z) (needed for the
        full accuracy of ``method="cheb"``).
    g : gravity [m/s²].   n_modes : number of modes returned, fastest first.
    lid : "free" (surface condition (13.65); index 0 is the barotropic mode) or "rigid" (dψ/dz = 0 at z = 0; the
        barotropic mode, c = ∞, is left out and index 0 is the first baroclinic mode).
    method : "fd" — linear finite elements with a lumped mass matrix (equivalently finite volumes) on the given nodes:
        the weak form ∫ψ′φ′/N² dz + ψ(0)φ(0)/g = (1/c²)∫ψφ dz gives a symmetric tridiagonal generalised eigenproblem,
        solved with ``scipy.linalg.eigh_tridiagonal``; second order in the node spacing.
        "cheb" — Chebyshev collocation (``core.stability.cheb``) with ``n_cheb + 1`` nodes on [−H, 0], boundary rows
        replacing the first and last equations; spectrally accurate for a smooth callable N2; the modes are
        interpolated back to ``z``.  An independent route for the cross-check.
    n_cheb : polynomial degree of the Chebyshev route.
    Returns
    -------
    :class:`Modes` (c, psi, He, z, weights).
    Assumptions: linear, hydrostatic (frequencies ≪ N), flat bottom, no mean shear; with or without rotation.
    Numerics: the eigenvalue is 1/c²; only the lowest ``n_modes`` are returned (the high modes of a grid are not
    converged).  The barotropic eigenvalue is about N²H/(gπ²) times the first baroclinic one (1e-4 for an ocean): the
    symmetric "fd" solver keeps it to 1e-8, the non-symmetric Chebyshev solve loses digits as n_cheb grows —
    see Validation.
    Validation (measured, uniform N = 2.7e-3 rad/s, H = 4200 m, against :func:`modes_uniform_N`): "fd" relative
    error of c₁ 4.1e-5, 1.0e-5, 2.6e-6, 6.4e-7 at 101, 201, 401, 801 nodes (order 2) and of c₀ below 2e-8; "cheb"
    baroclinic speeds to 2e-10 at n_cheb = 32–64, c₀ to 2e-9 at 32 and 2e-7 at 64 (round-off of the differentiation
    matrices grows as n⁴: use "fd" or shooting when c₀ matters); V4 ∫ψ_mψ_n dz off-diagonal at 5e-16 of the diagonal
    for "fd" modes.  Label: converged.
    """
    z, N2a = _check_grid(z, N2)
    if lid not in ("free", "rigid"):
        raise ValueError('lid must be "free" or "rigid"')
    nm = int(n_modes)
    if nm < 1 or nm > z.size - 3:
        raise ValueError("n_modes must be between 1 and len(z) - 3")
    w = _trapz_weights(z)
    skip = 1 if lid == "rigid" else 0          # the rigid lid has a null mode psi = const (c = infinity): drop it
    if method == "fd":
        h = np.diff(z)
        kel = 0.5 * (1.0 / N2a[1:] + 1.0 / N2a[:-1]) / h          # element stiffness: mean of 1/N^2 over the element
        d = np.zeros(z.size)
        d[:-1] += kel
        d[1:] += kel
        if lid == "free":
            d[-1] += 1.0 / g                                       # surface term of Eq. (13.65)
        s = 1.0 / np.sqrt(w)                                       # symmetrise: M^{-1/2} A M^{-1/2}
        lam, vec = eigh_tridiagonal(d * s * s, -kel * s[1:] * s[:-1], select="i", select_range=(skip, nm - 1 + skip))
        psi = (vec * s[:, None]).T
    elif method == "cheb":
        lam, psi = _cheb_modes(z, N2, N2a, g, nm, lid, int(n_cheb), skip)
    else:
        raise ValueError('method must be "fd" or "cheb"')
    if np.any(lam <= 0):
        raise ValueError("vertical_modes: a non-positive eigenvalue was returned (grid too coarse for this profile)")
    c = 1.0 / np.sqrt(lam)                                         # Eq. (13.56): the eigenvalue is 1/c^2
    return Modes(c=c, psi=_normalise(psi), He=c ** 2 / g, z=z, weights=w)


def _cheb_modes(z, N2, N2a, g, nm, lid, n, skip):
    D, zc = cheb(n, (z[0], z[-1]))                 # nodes descending: zc[0] = 0 (surface), zc[-1] = -H (bottom)
    N2c = _F(N2(zc)) * np.ones_like(zc) if callable(N2) else np.interp(zc[::-1], z, N2a)[::-1]
    A = -D @ np.diag(1.0 / N2c) @ D                # -(psi'/N^2)' = lambda psi
    B = np.eye(n + 1)
    A[-1], B[-1] = D[-1], 0.0                      # bottom: psi' = 0, Eq. (13.64)
    A[0] = D[0] + (N2c[0] / g if lid == "free" else 0.0) * np.eye(n + 1)[0]   # surface, Eq. (13.65) or rigid lid
    B[0] = 0.0
    lam, V = eig(A, B)
    Nc = np.sqrt(N2c)
    ref = (np.pi / abs(float(np.sum(0.5 * (Nc[1:] + Nc[:-1]) * np.diff(zc))))) ** 2   # WKB size of 1/c_1^2
    ok = np.isfinite(lam) & (np.abs(lam.imag) <= 1e-6 * np.maximum(np.abs(lam.real), 1e-6 * ref))
    if skip:                                       # rigid lid: drop the null mode (eigenvalue 0 to round-off)
        ok &= np.abs(lam.real) > 1e-6 * ref
    else:
        ok &= lam.real > 0
    lam, V = lam[ok].real, V[:, ok].real
    order = np.argsort(lam)
    lam, V = lam[order], V[:, order]
    if lam.size < nm:
        raise ValueError("Chebyshev route found fewer modes than requested; raise n_cheb")
    from scipy.interpolate import BarycentricInterpolator

    return lam[:nm], np.array([BarycentricInterpolator(zc, V[:, j])(z) for j in range(nm)])


def _shoot(lam: np.ndarray, z: np.ndarray, N2_fn: Callable) -> tuple:
    """Integrate S′ = N² P, P′ = −λ S (S = ψ, P = ψ′/N²) from the bottom (S = 1, P = 0) to the surface with the
    classical RK4 on the given nodes, for an array of trial eigenvalues λ at once.  Returns S (nλ, nz) and P(0)."""
    lam = np.atleast_1d(_F(lam))
    S = np.empty((lam.size, z.size))
    s = np.ones(lam.size)
    p = np.zeros(lam.size)
    S[:, 0] = s
    for j in range(z.size - 1):
        h = z[j + 1] - z[j]
        n0, nh, n1 = float(N2_fn(z[j])), float(N2_fn(z[j] + 0.5 * h)), float(N2_fn(z[j + 1]))
        k1s, k1p = n0 * p, -lam * s
        k2s, k2p = nh * (p + 0.5 * h * k1p), -lam * (s + 0.5 * h * k1s)
        k3s, k3p = nh * (p + 0.5 * h * k2p), -lam * (s + 0.5 * h * k2s)
        k4s, k4p = n1 * (p + h * k3p), -lam * (s + h * k3s)
        s = s + h / 6.0 * (k1s + 2 * k2s + 2 * k3s + k4s)
        p = p + h / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        S[:, j + 1] = s
    return S, p


def vertical_modes_shooting(z, N2_fn: Callable, g: float = G0, n: int = 1, lid: str = "free", scan: int = 60):
    """One vertical mode by shooting: integrate the mode equation up from the bottom and find the eigenvalue that
    satisfies the surface condition.

    Book: §13.9, Eq. (13.56) written as the pair ψ′ = N² P, P′ = −ψ/c² with P = ψ′/N²; bottom condition (13.64)
    P(−H) = 0; surface condition (13.65) P(0) + ψ(0)/g = 0 (rigid lid: P(0) = 0).
    Parameters
    ----------
    z : nodes [m], increasing from −H to 0 (the RK4 steps).   N2_fn : callable z → N² [1/s²] > 0 (scalar argument).
    g [m/s²].   n : the book's mode number (0 = barotropic, free surface only; n ≥ 1 baroclinic).
    lid : "free" or "rigid".   scan : trial eigenvalues per decade used to bracket the root (default 60; adequate for
    n ≤ 20).
    Returns
    -------
    (c_n, psi_n) : speed [m/s] and the mode at the nodes (ψ(0) = 1).
    Numerics (ours): classical RK4 (fourth order in the node spacing); the surface residual is scanned on a
    logarithmic grid of 1/c² from far below the barotropic value to above the WKB estimate of mode n + 1, the n-th
    sign change is refined with ``brentq`` (relative xtol 1e-14).  Independent of the matrix routes of
    :func:`vertical_modes`, and the method the explainer's JavaScript mirrors.
    Raises ValueError if the requested root is not bracketed (n too large for the scan) or for n = 0 with a rigid lid.
    Assumptions: as :func:`vertical_modes`.
    Validation (measured, uniform N, 401 nodes, against Eq. (13.69)): relative error of c_n 2e-16, 3e-11, 5e-10,
    3e-9 for n = 0 … 3.  Label: converged.
    """
    z = _F(z)
    if z.ndim != 1 or z.size < 5 or np.any(np.diff(z) <= 0):
        raise ValueError("z must be strictly increasing from -H to 0 with at least 5 nodes")
    n = int(n)
    if lid not in ("free", "rigid"):
        raise ValueError('lid must be "free" or "rigid"')
    if n < 0 or (lid == "rigid" and n == 0):
        raise ValueError("mode number must be >= 0 (free surface) or >= 1 (rigid lid: the barotropic speed is infinite)")
    H = z[-1] - z[0]
    Nz = np.sqrt(np.array([float(N2_fn(v)) for v in z]))
    if np.any(~np.isfinite(Nz)) or np.any(Nz <= 0):
        raise ValueError("N2_fn must be positive on the grid")
    I = float(np.sum(0.5 * (Nz[1:] + Nz[:-1]) * np.diff(z)))          # WKB phase integral of N
    lam_hi = ((n + 1.6) * np.pi / I) ** 2
    lam_lo = min(1e-3 / (g * H), 1e-3 * (np.pi / I) ** 2)
    npts = max(int(scan * np.log10(lam_hi / lam_lo)) + 2, 20)
    trial = np.geomspace(lam_lo, lam_hi, npts)

    def resid(lam):
        S, p = _shoot(lam, z, N2_fn)
        return p + (S[:, -1] / g if lid == "free" else 0.0)

    r = resid(trial)
    idx = np.nonzero(np.sign(r[:-1]) * np.sign(r[1:]) < 0)[0]
    want = n if lid == "free" else n - 1
    if idx.size <= want:
        raise ValueError(f"vertical_modes_shooting: mode {n} not bracketed (found {idx.size} roots); raise scan")
    a, b = trial[idx[want]], trial[idx[want] + 1]
    lam = brentq(lambda v: float(resid(np.array([v]))[0]), a, b, xtol=1e-14 * b, rtol=1e-14)
    S, _ = _shoot(np.array([lam]), z, N2_fn)
    return float(1.0 / np.sqrt(lam)), S[0] / S[0, -1]


def _uniform_roots(N: float, H: float, g: float, count: int) -> np.ndarray:
    """Roots X_n = N H/c_n of tan X = (N²H/g)/X, n = 0 … count − 1 (one per branch (nπ, nπ + π/2))."""
    eps = N * N * H / g
    X = np.zeros(count)
    for n in range(count):
        lo, hi = n * np.pi + (1e-300 if n == 0 else 0.0), n * np.pi + 0.5 * np.pi
        X[n] = brentq(lambda x: x * np.sin(x) - eps * np.cos(x), lo, hi, xtol=1e-15, rtol=1e-15)   # Eq. (13.69)
    return X


def modes_uniform_N(N: float, H: float, g: float = G0, n_modes: int = 4, lid: str = "free", nz: int = 201) -> Modes:
    """Exact vertical modes of a layer of uniform buoyancy frequency.

    Book: §13.9, Eqs. (13.66)–(13.71): ψ_n = A_n cos(Nz/c_n) + B_n sin(Nz/c_n), B_n = −(c_n N/g) A_n (13.68), with the
    eigenvalues from tan(NH/c_n) = c_n N/g (13.69); c₀ ≈ sqrt(gH) (13.70); rigid lid: c_n = NH/(nπ) (13.71) and
    ψ_n = cos(nπz/H).
    Parameters
    ----------
    N : buoyancy frequency [rad/s].   H : depth [m].   g [m/s²].   n_modes : number of modes, fastest first.
    lid : "free" (index 0 barotropic) or "rigid" (index 0 = first baroclinic mode).
    nz : number of uniform nodes from −H to 0 on which ψ is returned.
    Returns
    -------
    :class:`Modes`.  The roots themselves are X_n = N·H/c_n; the first, X₀ ≈ N sqrt(H/g), is a few hundredths for an
    ocean — the book's sentence "the first root occurs for NH/c_n = 1" should read "≪ 1" (slip #6).
    Numerics: with ε = N²H/g the roots solve X sin X − ε cos X = 0 on (nπ, nπ + π/2) (``brentq``, xtol 1e-15); this
    form has no singularity there, unlike tan X − ε/X.
    Assumptions: uniform N, hydrostatic, flat bottom.
    Validation: V1 residual of (13.69) below 1e-12; c₀ → sqrt(gH) and c_n → NH/(nπ) as ε → 0 (measured for
    N = 2.7e-3 rad/s, H = 4200 m: c₀/sqrt(gH) − 1 = 5.2e-4, c₁/(NH/π) − 1 = −3.2e-4).  Label: analytic.
    """
    N, H = float(N), float(H)
    if N <= 0 or H <= 0:
        raise ValueError("modes_uniform_N: N and H must be positive")
    zz = np.linspace(-H, 0.0, int(nz))
    nm = int(n_modes)
    if lid == "free":
        c = N * H / _uniform_roots(N, H, g, nm)
        psi = np.array([np.cos(N * zz / cn) - cn * N / g * np.sin(N * zz / cn) for cn in c])   # Eqs. (13.67), (13.68)
    elif lid == "rigid":
        n = np.arange(1, nm + 1)
        c = N * H / (n * np.pi)                                                               # Eq. (13.71)
        psi = np.cos(np.outer(n, np.pi * zz / H))
    else:
        raise ValueError('lid must be "free" or "rigid"')
    return Modes(c=c, psi=psi, He=c ** 2 / g, z=zz, weights=_trapz_weights(zz))


def rigid_lid_error(N: float, H: float, g: float = G0, n: int = 1) -> float:
    """Relative change of a baroclinic mode speed between the free surface and the rigid lid:
    (c_n,free − c_n,rigid)/c_n,rigid.

    Book: §13.9 ("the baroclinic modes are negligibly distorted by the rigid lid approximation"; Eqs. (13.69), (13.71)).
    Parameters: N [rad/s]; H [m]; g [m/s²]; n ≥ 1.
    Returns the relative difference (dimensionless; negative and about −N²H/(g n²π²): the free surface makes the
    baroclinic modes slightly slower).
    Assumptions: uniform N.   Validation: V1 against the leading-order estimate (measured −3.16e-4 against −3.16e-4
    for N = 2.7e-3 rad/s, H = 4200 m, n = 1).  Label: analytic.
    """
    if int(n) < 1:
        raise ValueError("rigid_lid_error: n must be >= 1 (the barotropic mode has no rigid-lid counterpart)")
    X = _uniform_roots(float(N), float(H), g, int(n) + 1)[int(n)]
    return float(int(n) * np.pi / X - 1.0)


def w_structure(modes: Modes, n: int):
    """Vertical structure of the vertical velocity of mode ``n``: ∫_{−H}^{z} ψ_n dz′.

    Book: §13.9, Eq. (13.53).
    Parameters: modes :class:`Modes`; n index into it.
    Returns the cumulative integral [m] at the nodes (zero at the bottom; cumulative trapezoid rule, second order).
    Assumptions: none.   Validation: V1 rigid-lid uniform-N mode: (H/nπ) sin(nπz/H), zero at both ends.
    Label: analytic.
    """
    z, psi = modes.z, modes.psi[int(n)]
    return np.concatenate([[0.0], np.cumsum(0.5 * (psi[1:] + psi[:-1]) * np.diff(z))])


def rho_structure(modes: Modes, n: int):
    """Vertical structure of the density perturbation of mode ``n``: dψ_n/dz.

    Book: §13.9, Eq. (13.54).
    Parameters: modes :class:`Modes`; n index into it.
    Returns dψ/dz [1/m] at the nodes (3-point stencils, second order including the ends).
    Assumptions: none.   Validation: V1 cosine mode.  Label: analytic.
    """
    return _d1(modes.z, modes.psi[int(n)])


def orthogonality_matrix(modes: Modes, kind: str = "psi", N2=None, g: float = G0, lid: str = "free"):
    """Normalised matrix of inner products of the modes — the identity when they are orthogonal.

    Book: §13.9 (the statement that the solutions of the Sturm–Liouville equation (13.56) are orthogonal).
    Parameters
    ----------
    modes : :class:`Modes`.
    kind : "psi" (default) — G_mn = ∫ψ_mψ_n dz / sqrt(∫ψ_m² dz ∫ψ_n² dz), **weight 1 and no surface term, for a free
        surface as well as for a rigid lid** (derivation in the module docstring); "energy" — the second relation,
        E_mn = ∫ψ_m′ψ_n′/N² dz + ψ_m(0)ψ_n(0)/g (the surface term only for ``lid="free"``), normalised the same way;
        it needs ``N2`` (values at the nodes or a scalar), ``g`` and ``lid``.
    Returns
    -------
    (n_modes, n_modes) array with unit diagonal.
    Numerics: trapezoid weights for "psi" (exact discrete orthogonality for ``method="fd"`` modes); element-wise
    differences for "energy".
    Assumptions: as :func:`vertical_modes`.
    Validation (measured): "fd" modes — off-diagonal 5e-16; exact uniform-N free-surface modes sampled on 401 nodes —
    off-diagonal 7e-9 with weight 1 (the trapezoid error), which confirms that no surface term is needed.
    Label: analytic.
    """
    psi, z = modes.psi, modes.z
    if kind == "psi":
        G = (psi * modes.weights) @ psi.T
    elif kind == "energy":
        if N2 is None:
            raise ValueError('orthogonality_matrix(kind="energy") needs N2')
        N2a = _F(N2) * np.ones_like(z)
        h = np.diff(z)
        dpsi = np.diff(psi, axis=1) / h
        G = (dpsi * (0.5 * (1.0 / N2a[1:] + 1.0 / N2a[:-1]) * h)) @ dpsi.T
        if lid == "free":
            G = G + np.outer(psi[:, -1], psi[:, -1]) / g
    else:
        raise ValueError('kind must be "psi" or "energy"')
    d = np.sqrt(np.diag(G))
    return G / np.outer(d, d)


def project(modes: Modes, profile):
    """Coefficients of a vertical profile on the modes: a_n = ∫ q ψ_n dz / ∫ ψ_n² dz.

    Book: §13.9, Eq. (13.52) (the expansion of u, v, p/ρ₀ in the modes), inverted with the orthogonality of the ψ_n.
    Parameters: modes :class:`Modes`; profile values at the nodes (nodes on the last axis; leading axes carried).
    Returns coefficients, modes on the last axis (same unit as the profile).
    Assumptions: the modes are orthogonal with weight 1 (true for both lids).  With a rigid lid the depth mean of
    the profile belongs to the omitted barotropic mode and is not represented.
    Validation: V1 projecting a mode returns a unit vector; project ∘ reconstruct = identity.  Label: analytic.
    """
    q = _F(profile)
    return ((q * modes.weights) @ modes.psi.T) / np.sum(modes.psi ** 2 * modes.weights, axis=1)


def reconstruct(modes: Modes, coeffs):
    """Sum of modes Σ a_n ψ_n(z).

    Book: §13.9, Eq. (13.52).
    Parameters: modes :class:`Modes`; coeffs (…, n_modes).   Returns the profile at the nodes (…, nz).
    Assumptions: none.   Validation: V1 inverse of :func:`project` on the span of the modes.  Label: analytic.
    """
    return _F(coeffs) @ modes.psi


def modal_amplitudes(c_n, p_n_t, p_n=None, rho0: float = 1.0, g: float = G0) -> dict:
    """Density and vertical-velocity amplitudes of a mode from its pressure amplitude.

    Book: §13.9, Eqs. (13.60)–(13.61): p_n = −(g/ρ₀) ρ_n, w_n = (1/c_n²) ∂p_n/∂t.
    Parameters
    ----------
    c_n : modal speed [m/s].   p_n_t : ∂p_n/∂t, the rate of change of the modal pressure amplitude p/ρ₀ [m²/s³].
    p_n : the modal pressure amplitude itself [m²/s²] (optional; needed for ``rho_n``).
    rho0 : reference density [kg/m³] (default 1: ``rho_n`` is then per unit reference density).   g [m/s²].
    Returns
    -------
    dict: ``w_n`` = p_n_t/c_n² [1/s] (it multiplies ∫ψ_n dz) and ``rho_n`` = −(ρ₀/g) p_n [kg/m²] (it multiplies
    dψ_n/dz) — ``None`` when ``p_n`` is not given, because Eq. (13.60) needs the amplitude, not its rate of change.
    Assumptions: as :func:`vertical_modes`.   Validation: V1 units; identities.  Label: analytic.
    """
    out = dict(rho_n=None, w_n=_S(_F(p_n_t) / float(c_n) ** 2))   # Eq. (13.61)
    if p_n is not None:
        out["rho_n"] = _S(-rho0 * _F(p_n) / g)                  # Eq. (13.60)
    return out


def wkb_mode_speed(z, N, n: int = 1):
    """WKB estimate of a baroclinic mode speed for a varying N(z): c_n ≈ (1/nπ) ∫_{−H}^{0} N dz.

    Book: **ours — not in the book**; it reduces to Eq. (13.71), c_n = NH/(nπ), for uniform N.
    Parameters: z nodes [m]; N buoyancy frequency at the nodes [rad/s] (not N²); n ≥ 1.
    Returns c_n [m/s] (trapezoid rule).
    Assumptions: N varies slowly over a vertical wavelength of the mode, rigid lid.  Not a good estimate for the
    gravest modes of a sharp thermocline: measured against :func:`vertical_modes` on the default profile of
    ``ch13.thermocline_N2`` (H = 4200 m) the error is quoted in that function's docstring.
    Validation: V1 uniform N exact; the error falls with n.  Label: approximate (WKB).
    """
    if int(n) < 1:
        raise ValueError("wkb_mode_speed: n must be >= 1")
    z, Nz = _F(z), _F(N) * np.ones_like(_F(z))
    return float(np.sum(0.5 * (Nz[1:] + Nz[:-1]) * np.diff(z)) / (int(n) * np.pi))
