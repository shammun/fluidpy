"""One-dimensional diffusion: the model behind molecular transport of species, heat and momentum.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §1.5. The book states the flux laws (1.1) Fick, (1.2) Fourier and
(1.3) Newton and describes profiles relaxing (Figs. 1.2–1.3); it does not write the diffusion equation in Ch. 1.
We use the model PDE ``∂f/∂t = D ∂²f/∂y²`` with D = κ_m (species), k/(ρ C_p) (heat) or ν (momentum), which follows
from each flux law plus conservation (derived formally in Ch. 4). Reused by Ch. 5 (vorticity diffusion), Ch. 8
(Stokes' problems, Couette start-up), Ch. 10 (explicit schemes), Ch. 12.

Grid convention: nodes y_i = y_0 + i Δy, i = 0 … N−1, values live at nodes; boundary nodes are the walls.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from ._util import as_scalar_if_0d, require_nonnegative, require_positive


def stable_time_step(D: float, dy: float, safety: float = 0.9) -> float:
    """Largest stable FTCS time step times a safety factor: ``dt = safety · dy^2 / (2 D)``.

    Book: tool for §1.5 diffusion demos (FTCS stability limit r = D dt/dy^2 <= 1/2; our scheme choice).

    Parameters
    ----------
    D : float
        Diffusivity [m^2/s], > 0.
    dy : float
        Grid spacing [m].
    safety : float, optional
        Safety factor (< 1).

    Returns
    -------
    dt : float
        [s].

    Validation: V1 r = D dt/dy^2 = safety/2 (0.45 for safety 0.9); ``ftcs_stable_time_step`` is the same object.
    Label: analytic.
    """
    require_positive("D", D)
    return safety * dy ** 2 / (2.0 * D)  # FTCS limit: D dt/dy^2 <= 1/2


ftcs_stable_time_step = stable_time_step  # alias


def ftcs_diffusion_1d(f0, D: float, dy: float, dt: float, nsteps: int,
                      bc: Sequence[str] = ("dirichlet", "dirichlet"), values: Sequence[float] | None = None,
                      save_every: int = 1, check_stability: bool = True):
    """March the 1-D diffusion equation with the explicit FTCS scheme.

    Book: §1.5 (molecular diffusion of species, heat and momentum; the relaxation of the profiles in Figs. 1.2–1.3).
    Model PDE (ours in Ch. 1): ``∂f/∂t = D ∂²f/∂y²``;
    scheme: ``f_i^{n+1} = f_i^n + r (f_{i+1}^n − 2 f_i^n + f_{i−1}^n)`` with ``r = D dt/dy^2``.

    Parameters
    ----------
    f0 : array_like, shape (N,)
        Initial profile at the nodes (units of f: kg/kg for Y, K for T, m/s for u).
    D : float
        Diffusivity [m^2/s] (κ_m, k/(ρ C_p) or ν), >= 0.
    dy : float
        Node spacing [m].
    dt : float
        Time step [s].
    nsteps : int
        Number of steps.
    bc : pair of str, optional
        Boundary condition at (y_0, y_{N−1}), case-insensitive: ``"dirichlet"`` (value held fixed), ``"neumann"``
        (zero flux, mirror ghost node) or ``"periodic"`` (must be given at both ends). Any other name raises
        ValueError. Periodic node convention: node N−1 and node 0 are *neighbours* (one spacing dy apart), so the
        period is N·dy and the grid must **not** repeat the endpoint — use ``y = y0 + dy*np.arange(N)`` (or
        ``np.linspace(a, b, N, endpoint=False)``), not ``np.linspace(a, b, N)``.
    values : pair of float, optional
        Dirichlet values; default the end values of ``f0``.
    save_every : int, optional
        Store every k-th step (the first and last states are always stored).
    check_stability : bool, optional
        Raise ValueError if r > 1/2 (set False only to demonstrate the instability).

    Returns
    -------
    F : ndarray, shape (nsave, N)
        Profiles at the saved steps; ``F[0]`` is f0.

    Notes
    -----
    Method: forward Euler in time, second-order central differences in space (our choice; the book specifies no scheme).
    Accuracy O(dt, dy^2); stable for r <= 1/2. With zero-flux ends the trapezoid integral of f is conserved exactly.
    Assumptions: constant D, uniform grid.

    Validation: V3 converges to :func:`couette_startup_profile` at observed order 2.002 (N = 11…81, r = 0.4) and to the
    exact spreading Gaussian with zero-flux ends at order 2.007 (N = 61…481); V4 zero-flux ends conserve the trapezoid
    integral over 2000 steps (drift 2e-16 relative), periodic ends conserve the node sum (1e-10); V7 r = 0.51 raises,
    grows > 1e3× when forced, r = 0.49 decays; mixed periodic/non-periodic ends raise. Label: converged, conserved.
    """
    f = np.array(f0, dtype=float)
    if f.ndim != 1 or f.size < 3:
        raise ValueError("f0 must be a 1-D array with at least 3 nodes")
    require_nonnegative("D", D)
    r = D * dt / dy ** 2
    if check_stability and r > 0.5:
        raise ValueError(f"FTCS unstable: r = D dt/dy^2 = {r:.4g} > 1/2; reduce dt below {0.5 * dy**2 / D:.4g} s")
    if isinstance(bc, str) or len(bc) != 2:
        raise ValueError(f"bc must be a pair of names, got {bc!r}")
    bc = tuple(str(b).strip().lower() for b in bc)
    unknown = [b for b in bc if b not in ("dirichlet", "neumann", "periodic")]
    if unknown:
        raise ValueError(f"unknown boundary condition(s) {unknown}; use 'dirichlet', 'neumann' or 'periodic'")
    if ("periodic" in bc) and bc != ("periodic", "periodic"):
        raise ValueError("periodic boundary conditions must be applied at both ends")
    if values is None:
        values = (f[0], f[-1])
    saved = [f.copy()]
    for n in range(1, nsteps + 1):
        lap = np.empty_like(f)
        lap[1:-1] = f[2:] - 2.0 * f[1:-1] + f[:-2]  # central second difference
        if bc[0] == "periodic":
            lap[0] = f[1] - 2.0 * f[0] + f[-1]
            lap[-1] = f[0] - 2.0 * f[-1] + f[-2]
        else:
            lap[0] = 2.0 * (f[1] - f[0]) if bc[0] == "neumann" else 0.0  # mirror ghost node f[-1] = f[1]
            lap[-1] = 2.0 * (f[-2] - f[-1]) if bc[1] == "neumann" else 0.0
        f = f + r * lap  # FTCS update
        if bc[0] == "dirichlet":
            f[0] = values[0]
        if bc[1] == "dirichlet":
            f[-1] = values[1]
        if n % save_every == 0 or n == nsteps:
            saved.append(f.copy())
    return np.array(saved)


def gaussian_spreading(y, t, D, M=1.0, y0=0.0, t0=0.0):
    """Exact solution of ``∂f/∂t = D ∂²f/∂y²`` on an infinite line for a point release: a spreading Gaussian.

    Book: §1.5 (illustrates diffusive smoothing; not printed in Ch. 1 — the classical heat-kernel solution).
    ``f = M / sqrt(4π D (t + t0)) · exp(−(y − y0)^2 / (4 D (t + t0)))``.

    Parameters
    ----------
    y : float or array_like
        Position [m].
    t : float or array_like
        Time [s]; t + t0 > 0.
    D : float
        Diffusivity [m^2/s].
    M : float, optional
        Integral of f over y [units of f times m].
    y0 : float, optional
        Centre [m].
    t0 : float, optional
        Virtual time origin [s] (t0 > 0 gives a finite-width initial Gaussian of variance 2 D t0).

    Returns
    -------
    f : float or ndarray

    Validation: indirect only — the FTCS scheme converges to it at order 2.007 (N = 61…481), which a wrong kernel
    would not allow; no direct PDE-residual or integral test yet. Label: converged (indirect).
    """
    tau = np.asarray(t, dtype=float) + t0
    require_positive("t + t0", tau)
    y = np.asarray(y, dtype=float)
    return as_scalar_if_0d(M / np.sqrt(4.0 * np.pi * D * tau) * np.exp(-(y - y0) ** 2 / (4.0 * D * tau)))


def couette_startup_profile(y, t, U, h, nu, nterms: int | None = 200, tol: float = 1e-12, nmax: int = 20000):
    """Velocity between plates after the upper plate (y = h) suddenly starts moving at speed U (series solution).

    Book: §1.5 (momentum diffusion, Fig. 1.3) — analytic check of the FTCS run with D = ν. The steady limit is the
    linear Couette profile with uniform shear stress μ U/h (Eq. 1.3). Derived by separation of variables (full
    treatment in Ch. 8):
    ``u = U y/h − (2U/π) Σ_{n≥1} ((−1)^{n+1}/n) sin(nπ y/h) exp(−n²π² ν t/h²)``.

    Parameters
    ----------
    y : float or array_like
        Distance from the fixed lower plate [m], 0 <= y <= h.
    t : float
        Time since start [s], >= 0.
    U : float
        Upper-plate speed [m/s].
    h : float
        Gap width [m].
    nu : float
        Kinematic viscosity [m^2/s].
    nterms : int or None, optional
        Number of series terms (default 200, the count the E2 explainer mirrors). Truncation error is below
        exp(−(nterms π)^2 ν t/h^2) — negligible unless t << h^2/(ν nterms^2). ``None``: add terms adaptively until the
        amplitude bound exp(−n²π²νt/h²)/n < ``tol`` (capped at ``nmax``).
    tol : float, optional
        Adaptive stopping tolerance (only with ``nterms=None``).
    nmax : int, optional
        Hard cap on the adaptive number of terms.

    Returns
    -------
    u : float or ndarray
        [m/s].

    Notes
    -----
    At t = 0 the exact initial state (u = 0 below the plate, U at y = h) is returned; "at the plate" means
    |y − h| <= 1e-9 h (a tolerance relative to the gap, so it works for any h). Assumptions: incompressible
    Newtonian fluid, no pressure gradient, fluid initially at rest, no-slip.

    Validation: V1 t ≫ h^2/ν gives U y/h (1e-12); t = 0 returns the exact initial state; 200 terms and the adaptive
    sum agree (1e-12) at t = 0.2 s; V3 FTCS converges to it at observed order 2.002. Label: analytic, converged.
    """
    y = np.asarray(y, dtype=float)
    if t <= 0.0:
        at_plate = np.abs(y - h) <= 1e-9 * abs(h)  # relative tolerance (np.isclose's absolute 1e-8 fails for tiny h)
        return as_scalar_if_0d(np.where(at_plate, float(U), 0.0))
    a = np.pi ** 2 * nu * t / h ** 2
    # number of terms: exp(−a n²)/n < tol
    if nterms is None:
        n_needed = int(np.ceil(np.sqrt(max(np.log(1.0 / tol), 0.0) / a))) + 1
        n_all = np.arange(1, min(n_needed, nmax) + 1)
    else:
        n_all = np.arange(1, int(nterms) + 1)
    yy = y.ravel()
    series = np.zeros_like(yy)
    for start in range(0, n_all.size, 1000):  # chunks keep memory small for early times
        n = n_all[start:start + 1000]
        series += (((-1.0) ** (n + 1) / n * np.exp(-a * n ** 2))[:, None]
                   * np.sin(np.outer(n, np.pi * yy / h))).sum(axis=0)
    u = U * yy / h - (2.0 * U / np.pi) * series
    return as_scalar_if_0d(u.reshape(y.shape))


def derivative_2nd_order(f, y):
    """First derivative df/dy of sampled data, second-order accurate at every node including the ends.

    Book: tool for §1.5 flux laws (1.1)–(1.3) evaluated on profiles.

    Parameters
    ----------
    f : array_like, shape (N,)
        Samples.
    y : array_like, shape (N,)
        Strictly monotonic coordinates [m]; N >= 3.

    Returns
    -------
    dfdy : ndarray, shape (N,)

    Notes
    -----
    Uniform grid: interior ``(f_{i+1} − f_{i−1})/(2Δy)``, ends ``(∓3 f_0 ± 4 f_1 ∓ f_2)/(2Δy)`` written out explicitly
    (so the order can be measured). Non-uniform grid: ``np.gradient(f, y, edge_order=2)`` (also second order).

    Validation: V1 exact for a quadratic on a non-uniform grid (1e-12) and, through :func:`shear_stress_profile`, on a
    uniform grid; V3 observed order 1.996 on sin(3y) including the one-sided end stencils (N = 21…161).
    Label: analytic, converged.
    """
    f = np.asarray(f, dtype=float)
    y = np.asarray(y, dtype=float)
    if f.shape != y.shape or f.ndim != 1 or f.size < 3:
        raise ValueError("f and y must be 1-D arrays of equal length >= 3")
    dy = np.diff(y)
    if np.allclose(dy, dy[0], rtol=1e-10, atol=0.0):
        h = dy[0]
        d = np.empty_like(f)
        d[1:-1] = (f[2:] - f[:-2]) / (2.0 * h)  # central difference
        d[0] = (-3.0 * f[0] + 4.0 * f[1] - f[2]) / (2.0 * h)  # one-sided, second order
        d[-1] = (3.0 * f[-1] - 4.0 * f[-2] + f[-3]) / (2.0 * h)
        return d
    return np.gradient(f, y, edge_order=2)


__all__ = ["stable_time_step", "ftcs_stable_time_step", "ftcs_diffusion_1d", "gaussian_spreading", "couette_startup_profile",
           "derivative_2nd_order"]
