"""Small numerical models for rotating shallow water and two-dimensional vorticity dynamics:

* :class:`ShallowWater` with :func:`step`, :func:`run` … — the rotating shallow-water equations on an Arakawa C-grid
  (linear or nonlinear, f- or β-plane, doubly periodic, channel or closed basin), with volume, energy and
  potential-vorticity diagnostics and particle tracking;
* :func:`linear_1d_step`, :func:`linear_1d_run` — a forward–backward march of the 1-D linear set (the cheap
  geostrophic-adjustment run);
* :func:`qg_linear_evolve`, :func:`qg_linear_evolve_1d` — exact spectral evolution of the linear quasi-geostrophic
  vorticity equation (Rossby waves);
* :func:`barotropic_run` — a pseudo-spectral model of the barotropic vorticity equation (two-dimensional turbulence,
  with or without β).

Book: Kundu, Cohen & Dowling 5e, Ch. 13 §13.8 Eqs. (13.44)–(13.45); §13.13 Eqs. (13.88)–(13.94); §13.15 Eq. (13.117);
§13.16 Eq. (13.122); §13.18 Eqs. (13.143)–(13.144) (rendered pages chapters/pages/ch13/p670, p685–p687, p700, p703, p713).

**Every numerical scheme in this module is our choice — the book contains no numerical method** (DEVIATION, method only):
the C-grid with an energy-conserving average of the Coriolis/vorticity terms, SSP-RK3 time stepping, the
forward–backward 1-D scheme, the pseudo-spectral RK4 model with 2/3 de-aliasing.  Each docstring states the order of its
scheme and what is conserved, with measured numbers.

Layout (as ``core.mac``): arrays are ``[j, i]`` = (y, x).  Cell (j, i) has its centre at ((i + ½)dx, (j + ½)dy);
η and h live at centres, u at the west face (i dx, (j + ½)dy), v at the south face ((i + ½)dx, j dy), vorticity and
potential vorticity at the south-west corner (i dx, j dy).  All four arrays have shape ``[ny, nx]``.  Walls coincide
with faces: a channel has v = 0 on the row j = 0 (which is both y = 0 and y = Ly), a closed basin also has u = 0 on the
column i = 0.  A model **state** is a plain dict ``{"eta": …, "u": …, "v": …}``.
"""
from __future__ import annotations

import numpy as np

from ._util import as_scalar_if_0d
from .thermo import G0

__all__ = ["ShallowWater", "make_state", "gaussian_bump", "geostrophic_state", "kelvin_state", "continuity_tendency",
           "momentum_tendencies", "relative_vorticity", "potential_vorticity", "sw_potential_vorticity", "energy",
           "volume", "dt_limit", "step", "run", "interpolate", "advect_particles", "particle_potential_vorticity",
           "linear_1d_step", "linear_1d_run", "qg_linear_evolve", "qg_linear_evolve_1d", "spectral_wavenumbers",
           "barotropic_vorticity_rhs", "barotropic_run", "barotropic_time_step", "barotropic_invariants",
           "barotropic_velocity", "barotropic_spectrum", "spectral_centroids", "zonal_energy_fraction",
           "random_vorticity"]

_F = lambda a: np.asarray(a, dtype=float)  # noqa: E731
_S = as_scalar_if_0d


def _xm(a):  # value at i - 1
    return np.roll(a, 1, axis=1)


def _xp(a):  # value at i + 1
    return np.roll(a, -1, axis=1)


def _ym(a):
    return np.roll(a, 1, axis=0)


def _yp(a):
    return np.roll(a, -1, axis=0)


class ShallowWater:
    """Grid and parameters of the rotating shallow-water model (Arakawa C-grid).

    Book: §13.8, Eq. (13.44) (nonlinear continuity) and Eq. (13.45) (the linear set); §13.13, Eqs. (13.88)–(13.90)
    (the nonlinear set), solved in vector-invariant form:
        ∂u/∂t = +q · (h v) − ∂/∂x (g η + K),   ∂v/∂t = −q · (h u) − ∂/∂y (g η + K),   ∂h/∂t = −∇·(h u),
    with q = (ζ + f)/h the potential vorticity of Eq. (13.94), K = ½(u² + v²) and h = H − b + η.  ``linear=True`` in
    the functions below replaces q by f/H, h by the resting depth and drops K: Eq. (13.45).

    Parameters
    ----------
    nx, ny : number of cells.   Lx, Ly : domain size [m].   H : resting depth [m] (or equivalent depth of a mode).
    f0 : Coriolis parameter at mid-channel [1/s].   beta : df/dy [1/(m s)]; f = f0 + beta (y − Ly/2).
    g : gravity or reduced gravity [m/s²].
    bc : "periodic" (both directions), "channel" (periodic in x, walls at y = 0 and Ly) or "closed" (walls all round).
        A β-plane needs walls in y ("channel" or "closed").
    bottom : bottom elevation b(x, y) at the centres [m] (array ``[ny, nx]``) or None for a flat bottom.

    Attributes: ``x_c``, ``y_c`` (centre coordinates), ``x_u``, ``y_v`` (face coordinates), ``dx``, ``dy``, ``H``,
    ``Hrest`` (H − b at centres), ``f_u``, ``f_v``, ``f_q`` (f at u-points, v-points and corners, ``[ny, nx]``),
    ``mask_u``, ``mask_v``, ``mask_q`` (0 on walls), ``bc``, ``g``, ``f0``, ``beta``.

    Scheme (ours): second-order centred differences on the C-grid; the Coriolis/vorticity terms are averaged as
    avg_y(q · avg_x(h v)) and avg_x(q · avg_y(h u)), which makes the spatial discretisation conserve total energy
    (and volume to round-off, continuity being in flux form); potential enstrophy is not conserved.  Time: three-stage
    strong-stability-preserving Runge–Kutta, third order.  Stability: ω_max dt ≤ sqrt(3), ω_max² = f² + 4 c² (1/dx² +
    1/dy²), c = sqrt(g h_max) + |u|_max; :func:`step` raises ValueError beyond it.  Walls are free-slip.
    Assumptions: one homogeneous hydrostatic layer, inviscid.
    Validation (measured): see :func:`run`.  Label: converged.
    """

    def __init__(self, nx: int, ny: int, Lx: float, Ly: float, H: float, f0: float, beta: float = 0.0, g: float = G0,
                 bc: str = "periodic", bottom=None):
        if bc not in ("periodic", "channel", "closed"):
            raise ValueError('bc must be "periodic", "channel" or "closed"')
        if beta != 0.0 and bc == "periodic":
            raise ValueError("a beta-plane is not periodic in y: use bc='channel' or bc='closed'")
        if nx < 4 or ny < 4:
            raise ValueError("need at least 4 cells in each direction")
        self.nx, self.ny, self.Lx, self.Ly = int(nx), int(ny), float(Lx), float(Ly)
        self.dx, self.dy = self.Lx / self.nx, self.Ly / self.ny
        self.H, self.f0, self.beta, self.g, self.bc = float(H), float(f0), float(beta), float(g), bc
        self.x_c = (np.arange(self.nx) + 0.5) * self.dx
        self.y_c = (np.arange(self.ny) + 0.5) * self.dy
        self.x_u = np.arange(self.nx) * self.dx
        self.y_v = np.arange(self.ny) * self.dy
        self.b = np.zeros((self.ny, self.nx)) if bottom is None else _F(bottom) * np.ones((self.ny, self.nx))
        self.Hrest = self.H - self.b
        if np.any(self.Hrest <= 0):
            raise ValueError("the bottom reaches the resting surface (H - b must be positive)")
        self.mask_u = np.ones((self.ny, self.nx))
        self.mask_v = np.ones((self.ny, self.nx))
        self.mask_q = np.ones((self.ny, self.nx))
        if bc in ("channel", "closed"):
            self.mask_v[0, :] = 0.0
            self.mask_q[0, :] = 0.0
        if bc == "closed":
            self.mask_u[:, 0] = 0.0
            self.mask_q[:, 0] = 0.0
        ones = np.ones((1, self.nx))
        self.f_q = (self.f0 + self.beta * (self.y_v - 0.5 * self.Ly))[:, None] * ones    # corners
        self.f_u = (self.f0 + self.beta * (self.y_c - 0.5 * self.Ly))[:, None] * ones    # u-points (and centres)
        self.f_v = self.f_q.copy()                                                        # v-points

    def grids(self) -> dict:
        """Coordinates of the four families of points: dict of (X, Y) meshgrids for "center", "u", "v", "corner"."""
        return dict(center=np.meshgrid(self.x_c, self.y_c), u=np.meshgrid(self.x_u, self.y_c),
                    v=np.meshgrid(self.x_c, self.y_v), corner=np.meshgrid(self.x_u, self.y_v))

    # ---- internal operators (arrays in, arrays out) ---------------------------------------------------------------
    def _corner_avg(self, a):
        return 0.25 * (a + _xm(a) + _ym(a) + _ym(_xm(a)))

    def _deta(self, eta, u, v, linear):
        h = self.Hrest if linear else self.Hrest + eta
        U = 0.5 * (h + _xm(h)) * u * self.mask_u
        V = 0.5 * (h + _ym(h)) * v * self.mask_v
        return -(_xp(U) - U) / self.dx - (_yp(V) - V) / self.dy  # Eq. (13.44)

    def _zeta(self, u, v):
        return ((v - _xm(v)) / self.dx - (u - _ym(u)) / self.dy) * self.mask_q

    def _duv(self, eta, u, v, linear):
        if linear:
            h = self.Hrest
            q = self.f_q / self._corner_avg(h)
            B = self.g * eta
        else:
            h = self.Hrest + eta
            q = (self._zeta(u, v) + self.f_q) / self._corner_avg(h)
            B = self.g * eta + 0.5 * (0.5 * (u * u + _xp(u * u)) + 0.5 * (v * v + _yp(v * v)))
        U = 0.5 * (h + _xm(h)) * u * self.mask_u
        V = 0.5 * (h + _ym(h)) * v * self.mask_v
        Gq = q * 0.5 * (V + _xm(V))            # q * avg_x(h v) at corners
        Hq = q * 0.5 * (U + _ym(U))            # q * avg_y(h u) at corners
        du = (0.5 * (Gq + _yp(Gq)) - (B - _xm(B)) / self.dx) * self.mask_u     # Eq. (13.45b) / (13.88)
        dv = (-0.5 * (Hq + _xp(Hq)) - (B - _ym(B)) / self.dy) * self.mask_v    # Eq. (13.45c) / (13.89)
        return du, dv

    def _tend(self, eta, u, v, linear):
        du, dv = self._duv(eta, u, v, linear)
        return self._deta(eta, u, v, linear), du, dv

    def _energy(self, eta, u, v, linear):
        h = self.Hrest if linear else self.Hrest + eta
        ke = 0.5 * np.sum(0.5 * (h + _xm(h)) * u ** 2 * self.mask_u + 0.5 * (h + _ym(h)) * v ** 2 * self.mask_v)
        pe = 0.5 * self.g * np.sum(eta ** 2)
        dA = self.dx * self.dy
        return float(ke * dA), float(pe * dA)

    def _dt_max(self, eta=None, u=None, v=None):
        hmax = float(np.max(self.Hrest))
        umax = 0.0
        if eta is not None:
            hmax = max(hmax, float(np.max(self.Hrest + eta)))
            umax = float(max(np.max(np.abs(u)), np.max(np.abs(v))))
        c = np.sqrt(self.g * hmax) + umax
        fmax = float(np.max(np.abs(self.f_q)))
        return float(np.sqrt(3.0) / np.sqrt(fmax ** 2 + 4.0 * c ** 2 * (1.0 / self.dx ** 2 + 1.0 / self.dy ** 2)))


def _arr(model: ShallowWater, state: dict):
    return _F(state["eta"]), _F(state["u"]), _F(state["v"])


def make_state(model: ShallowWater, eta=None, u=None, v=None) -> dict:
    """Build a model state ``{"eta", "u", "v"}`` (zeros by default); wall-normal velocities are set to zero.

    Book: the variables of §13.8, Eq. (13.45).   Parameters: model; eta [m] at centres, u [m/s] at west faces, v
    [m/s] at south faces (arrays ``[ny, nx]`` or scalars).   Returns the state dict.
    Assumptions: none.   Validation: V1 shapes and wall values.  Label: analytic.
    """
    z = np.zeros((model.ny, model.nx))
    return dict(eta=z.copy() if eta is None else _F(eta) * np.ones_like(z),
                u=(z.copy() if u is None else _F(u) * np.ones_like(z)) * model.mask_u,
                v=(z.copy() if v is None else _F(v) * np.ones_like(z)) * model.mask_v)


def gaussian_bump(model: ShallowWater, amplitude: float, x0: float | None = None, y0: float | None = None,
                  radius: float | None = None):
    """Surface elevation of a Gaussian bump η = a exp(−r²/2R²) at the cell centres [m].

    Book: not in the book — an initial condition of ours for the adjustment and wave demos of §13.8–§13.13.
    Parameters: model; amplitude a [m]; x0, y0 centre [m] (default mid-domain); radius R [m] (default Lx/10).
    Returns eta ``[ny, nx]``.   Assumptions: none.   Validation: V1 peak value.  Label: analytic.
    """
    X, Y = model.grids()["center"]
    x0 = 0.5 * model.Lx if x0 is None else x0
    y0 = 0.5 * model.Ly if y0 is None else y0
    R = 0.1 * model.Lx if radius is None else radius
    return amplitude * np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2.0 * R * R))


def geostrophic_state(model: ShallowWater, eta) -> dict:
    """A model state in geostrophic balance with a given surface elevation: u = −(g/f) ∂η/∂y, v = (g/f) ∂η/∂x.

    Book: §13.15, Eq. (13.116) (the shallow-water form of Eqs. (13.11)–(13.12)).
    Parameters: model (f must not vanish anywhere); eta at centres [m].
    Returns the state dict, with the gradients averaged to the u- and v-points (second order).
    Assumptions: Ro ≪ 1.   Validation: V1 the linear tendencies of u and v vanish to O(dx²).  Label: converged.
    """
    eta = _F(eta) * np.ones((model.ny, model.nx))
    if np.any(model.f_u == 0) or np.any((model.f_v == 0) & (model.mask_v > 0)):
        raise ValueError("geostrophic_state: f vanishes inside the domain")
    gy = (eta - _ym(eta)) / model.dy * model.mask_v          # at v-points
    gx = (eta - _xm(eta)) / model.dx * model.mask_u          # at u-points
    u = -model.g / model.f_u * 0.25 * (gy + _xm(gy) + _yp(gy) + _yp(_xm(gy)))
    v = model.g / np.where(model.f_v == 0, 1.0, model.f_v) * 0.25 * (gx + _xp(gx) + _ym(gx) + _ym(_xp(gx)))
    return make_state(model, eta, u, v)


def kelvin_state(model: ShallowWater, eta0: float, mode: int = 1, wall: str = "south") -> dict:
    """Initial state of a Kelvin wave leaning on one wall of a channel (or basin).

    Book: §13.12, Eq. (13.87) — sampled on the C-grid with the travel direction for which the wave is trapped against
    the chosen wall (coast on the right of the direction of travel for f > 0, on the left for f < 0).
    Parameters: model with walls in y, f0 ≠ 0 and beta = 0; eta0 amplitude at the wall [m]; mode — number of
    wavelengths along the channel; wall "south" (y = 0) or "north" (y = Ly).
    Returns the state dict (v = 0).
    Assumptions: as ``core.gfd.kelvin_wave``; the channel should be several Rossby radii wide.
    Validation (measured): after one period on a 64 × 128 grid the wave differs from its start by 0.24 % of its
    amplitude, in both hemispheres.  Label: analytic.
    """
    from .gfd import kelvin_wave

    if model.bc == "periodic" or model.f0 == 0:
        raise ValueError("kelvin_state needs walls in y (bc='channel' or 'closed') and f0 != 0")
    k = 2.0 * np.pi * int(mode) / model.Lx
    sgn = 1 if model.f0 > 0 else -1
    d = sgn if wall == "south" else -sgn          # travel direction along x
    Xc, Yc = model.grids()["center"]
    Xu, Yu = model.grids()["u"]
    dist = (lambda Y: Y) if wall == "south" else (lambda Y: model.Ly - Y)
    kw = dict(eta0=eta0, k=k, H=model.H, f=abs(model.f0), g=model.g, direction=1)
    eta, _ = kelvin_wave(d * Xc, dist(Yc), 0.0, **kw)
    _, u = kelvin_wave(d * Xu, dist(Yu), 0.0, **kw)
    return make_state(model, eta, d * u, None)


def continuity_tendency(model: ShallowWater, state: dict, linear: bool = True):
    """∂η/∂t = −∂(h u)/∂x − ∂(h v)/∂y in flux form at the cell centres [m/s].

    Book: §13.8, Eq. (13.44); ``linear=True``: Eq. (13.45a), with the resting depth in the flux.
    Parameters: model; state dict; linear.   Returns an array ``[ny, nx]``.
    Numerics: h is averaged to the faces; the sum over a closed or periodic domain vanishes to round-off, so volume
    is conserved exactly.
    Assumptions: one layer.   Validation: V4 volume.  Label: conserved.
    """
    return model._deta(*_arr(model, state), linear)


def momentum_tendencies(model: ShallowWater, state: dict, linear: bool = True):
    """(∂u/∂t, ∂v/∂t) at the u- and v-points [m/s²].

    Book: §13.8, Eq. (13.45b, c) (``linear=True``) or §13.13, Eqs. (13.88)–(13.89) in vector-invariant form.
    Parameters: model; state dict; linear.   Returns two arrays ``[ny, nx]`` (zero on walls).
    Numerics: see :class:`ShallowWater` (energy-conserving average of the Coriolis/vorticity term, second order).
    Assumptions: one inviscid layer.   Validation: V1 plane waves.  Label: converged.
    """
    return model._duv(*_arr(model, state), linear)


def relative_vorticity(model: ShallowWater, state: dict):
    """Relative vorticity ζ = ∂v/∂x − ∂u/∂y at the cell corners [1/s].

    Book: §13.13 (unnumbered definition of ζ between Eqs. (13.91) and (13.92)).
    Parameters: model; state dict.   Returns ζ ``[ny, nx]`` (zero on walls: free slip).
    Assumptions: none.   Validation: V1 solid-body rotation.  Label: converged.
    """
    _, u, v = _arr(model, state)
    return model._zeta(u, v)


def potential_vorticity(model: ShallowWater, state: dict, linear: bool = False):
    """Potential vorticity q = (ζ + f)/h of a model state at the cell corners [1/(m s)].

    Book: §13.13, Eq. (13.94).  (The closed form for given numbers is ``core.gfd.potential_vorticity(zeta, f, h)``;
    the chapter module exports this gridded version as ``sw_potential_vorticity``.)
    Parameters: model; state dict; linear — return the linearised form (ζ − f η/H)/H + f/H instead.
    Returns q ``[ny, nx]``; NaN on wall corners (the depth average there would straddle the wall).
    Assumptions: one layer.   Validation: V4 constant on particles (:func:`particle_potential_vorticity`).
    Label: conserved.
    """
    eta, u, v = _arr(model, state)
    zeta = model._zeta(u, v)
    Hq = model._corner_avg(model.Hrest)
    if linear:
        q = (zeta - model.f_q * model._corner_avg(eta) / Hq) / Hq + model.f_q / Hq
    else:
        q = (zeta + model.f_q) / model._corner_avg(model.Hrest + eta)  # Eq. (13.94)
    return np.where(model.mask_q > 0, q, np.nan)


sw_potential_vorticity = potential_vorticity


def energy(model: ShallowWater, state: dict, linear: bool = True) -> dict:
    """Kinetic, potential and total energy of a model state per unit density [m⁵/s² = J per kg/m³].

    Book: the invariant of §13.8, Eq. (13.45) (``linear=True``: resting depth in the kinetic energy) or of
    Eqs. (13.88)–(13.90); the discrete form (depth averaged to the faces) is ours — it is the one the spatial scheme
    conserves.
    Parameters: model; state dict; linear.   Returns dict(kinetic, potential, total).
    Assumptions: one layer; potential energy ½ g η² relative to the resting surface.
    Validation: V4 drift falls as dt³ (see :func:`run`).  Label: conserved.
    """
    ke, pe = model._energy(*_arr(model, state), linear)
    return dict(kinetic=ke, potential=pe, total=ke + pe)


def volume(model: ShallowWater, state: dict) -> float:
    """Volume anomaly ∫η dA [m³] — conserved to round-off by the flux-form continuity equation.

    Book: §13.8, Eq. (13.44) integrated over the domain.   Parameters: model; state dict.   Returns a float.
    Assumptions: closed or periodic domain.   Validation: V4.  Label: conserved.
    """
    return float(np.sum(_F(state["eta"])) * model.dx * model.dy)


def dt_limit(model: ShallowWater, courant: float = 0.5, state: dict | None = None) -> float:
    """Time step = ``courant`` × the stability bound sqrt(3)/ω_max of the SSP-RK3 / C-grid scheme [s].

    Book: not in the book (ours).  ω_max² = f_max² + 4 c² (1/dx² + 1/dy²) covers both the gravity-wave CFL condition
    and the inertial limit; c = sqrt(g H_max), plus the largest depth and speed of ``state`` if one is given (use
    that for nonlinear runs).
    Parameters: model; courant ≤ 1 (1 is the bound itself); state (optional).   Returns dt [s].
    Assumptions: none.   Validation: V3 a linear Poincaré wave is stable at courant 1 and :func:`step` refuses 1.01.
    Label: converged.
    """
    if state is None:
        return float(courant * model._dt_max())
    return float(courant * model._dt_max(*_arr(model, state)))


def step(model: ShallowWater, state: dict, dt: float, linear: bool = True) -> dict:
    """One SSP-RK3 step of the shallow-water model (third order in dt).

    Book: §13.8, Eq. (13.45) / §13.13, Eqs. (13.88)–(13.90); scheme ours.
    Parameters: model; state dict; dt [s]; linear.   Returns the new state dict.
    Raises ValueError if dt exceeds the stability bound (:func:`dt_limit` with courant = 1).
    Assumptions: one inviscid layer.   Validation: V3 third order in dt (error ratio 8 per halving).
    Label: converged.
    """
    e, u, v = _arr(model, state)
    lim = model._dt_max() if linear else model._dt_max(e, u, v)
    if dt > lim * (1.0 + 1e-12):
        raise ValueError(f"dt = {dt:.4g} s exceeds the stability limit {lim:.4g} s of the scheme "
                         "(sqrt(3)/omega_max); reduce dt or coarsen the grid")
    k1 = model._tend(e, u, v, linear)
    e1, u1, v1 = e + dt * k1[0], u + dt * k1[1], v + dt * k1[2]
    k2 = model._tend(e1, u1, v1, linear)
    e2 = 0.75 * e + 0.25 * (e1 + dt * k2[0])
    u2 = 0.75 * u + 0.25 * (u1 + dt * k2[1])
    v2 = 0.75 * v + 0.25 * (v1 + dt * k2[2])
    k3 = model._tend(e2, u2, v2, linear)
    return dict(eta=e / 3.0 + 2.0 / 3.0 * (e2 + dt * k3[0]), u=u / 3.0 + 2.0 / 3.0 * (u2 + dt * k3[1]),
                v=v / 3.0 + 2.0 / 3.0 * (v2 + dt * k3[2]))


def run(model: ShallowWater, state: dict, t_end: float, dt: float | None = None, linear: bool = True,
        save_every: int = 1) -> dict:
    """Integrate the shallow-water model from t = 0 to ``t_end`` and keep the history.

    Book: §13.8, Eq. (13.45) / §13.13, Eqs. (13.88)–(13.90); scheme ours (see :class:`ShallowWater`).
    Parameters
    ----------
    model; state : initial state dict.   t_end [s].
    dt : time step [s]; default ``dt_limit(model, 0.5)`` (with the state's speed for nonlinear runs), shortened so
        that an integer number of steps reaches t_end.
    linear : Eq. (13.45) if True, the nonlinear set (13.88)–(13.90) if False.
    save_every : keep every this many steps (the first and the last state are always kept).
    Returns
    -------
    dict: ``t`` (n,), ``eta``, ``u``, ``v`` (n, ny, nx) and ``energy`` (n,) — the total energy per unit density.
    Assumptions: one inviscid layer.
    Validation (measured; H = 4200 m, f at 35°, dt = dt_limit(0.5)): a linear Poincaré wave after one period has
    error 3.8e-2, 9.4e-3, 2.4e-3 of its amplitude on 32², 64², 128² cells (order 2) and energy drift −1.2e-3,
    −1.6e-4, −2.0e-5; at fixed grid the time error falls 8-fold per halving of dt (order 3); volume changes by
    round-off (1e-17 relative).  A nonlinear geostrophic vortex (64², 3 inertial periods) loses 4.5e-6 of its energy.
    Label: converged.
    """
    e, u, v = (a.copy() for a in _arr(model, state))
    if t_end < 0:
        raise ValueError("t_end must be >= 0")
    dt0 = (model._dt_max() if linear else model._dt_max(e, u, v)) * 0.5 if dt is None else float(dt)
    nsteps = max(int(np.ceil(t_end / dt0 - 1e-9)), 0)
    dtt = t_end / nsteps if nsteps else dt0
    every = max(int(save_every), 1)
    cur = dict(eta=e, u=u, v=v)
    ts, es, us, vs, en = [0.0], [e.copy()], [u.copy()], [v.copy()], [sum(model._energy(e, u, v, linear))]
    for n in range(1, nsteps + 1):
        cur = step(model, cur, dtt, linear)
        if n % every == 0 or n == nsteps:
            ts.append(n * dtt)
            es.append(cur["eta"].copy())
            us.append(cur["u"].copy())
            vs.append(cur["v"].copy())
            en.append(sum(model._energy(cur["eta"], cur["u"], cur["v"], linear)))
    return dict(t=np.array(ts), eta=np.array(es), u=np.array(us), v=np.array(vs), energy=np.array(en))


def interpolate(model: ShallowWater, a, xp, yp, loc: str = "center"):
    """Bilinear interpolation of a gridded field to arbitrary positions.

    Book: not in the book (ours; used to follow particles, §13.13).
    Parameters: model; a field ``[ny, nx]`` living at "center", "u", "v" or "corner" points (``loc``); xp, yp [m].
    Returns the interpolated values (shape of xp).
    Numerics: second order; indices wrap periodically, so keep particles at least a cell away from walls.
    Assumptions: smooth field.   Validation: V1 exact for a bilinear field.  Label: converged.
    """
    ox, oy = dict(center=(0.5, 0.5), u=(0.0, 0.5), v=(0.5, 0.0), corner=(0.0, 0.0))[loc]
    xi = _F(xp) / model.dx - ox
    yj = _F(yp) / model.dy - oy
    i0 = np.floor(xi).astype(int)
    j0 = np.floor(yj).astype(int)
    a_, b_ = xi - i0, yj - j0
    i0, i1 = i0 % model.nx, (i0 + 1) % model.nx
    j0, j1 = j0 % model.ny, (j0 + 1) % model.ny
    a = np.asarray(a)
    return ((1 - a_) * (1 - b_) * a[j0, i0] + a_ * (1 - b_) * a[j0, i1] + (1 - a_) * b_ * a[j1, i0]
            + a_ * b_ * a[j1, i1])


def advect_particles(model: ShallowWater, history: dict, xp0, yp0, substeps: int = 2):
    """Paths of marked particles through a saved model run.

    Book: §13.13 — potential vorticity is conserved *following the motion* (Eq. (13.94)); the particles are how that
    statement is checked.
    Parameters
    ----------
    model; history : the dict returned by :func:`run` (``t``, ``u``, ``v``; save every step, or often enough that the
    flow changes little between frames).   xp0, yp0 : initial positions [m] (arrays of n_particles).
    substeps : RK4 steps per saved interval.
    Returns
    -------
    (xp, yp) : arrays (n_frames, n_particles) [m], wrapped in the periodic directions.
    Numerics (ours): classical RK4 in time with the velocity interpolated bilinearly in space and linearly in time
    between frames — fourth order in the sub-step for a steady flow, second order in the frame spacing otherwise.
    Assumptions: particles stay at least a cell away from walls.
    Validation: V1 a particle in solid-body rotation returns to its start.  Label: converged.
    """
    t, U, V = _F(history["t"]), history["u"], history["v"]
    xp, yp = _F(xp0).astype(float).copy(), _F(yp0).astype(float).copy()
    xs, ys = [xp.copy()], [yp.copy()]

    def vel(i, w, x, y):
        u = (1 - w) * interpolate(model, U[i], x, y, "u") + w * interpolate(model, U[i + 1], x, y, "u")
        v = (1 - w) * interpolate(model, V[i], x, y, "v") + w * interpolate(model, V[i + 1], x, y, "v")
        return u, v

    for i in range(len(t) - 1):
        h = (t[i + 1] - t[i]) / int(substeps)
        for s in range(int(substeps)):
            w0, wh, w1 = s / substeps, (s + 0.5) / substeps, (s + 1.0) / substeps
            k1 = vel(i, w0, xp, yp)
            k2 = vel(i, wh, xp + 0.5 * h * k1[0], yp + 0.5 * h * k1[1])
            k3 = vel(i, wh, xp + 0.5 * h * k2[0], yp + 0.5 * h * k2[1])
            k4 = vel(i, w1, xp + h * k3[0], yp + h * k3[1])
            xp = xp + h / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            yp = yp + h / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if model.bc != "closed":
            xp = np.mod(xp, model.Lx)
        if model.bc == "periodic":
            yp = np.mod(yp, model.Ly)
        xs.append(xp.copy())
        ys.append(yp.copy())
    return np.array(xs), np.array(ys)


def particle_potential_vorticity(model: ShallowWater, history: dict, xp, yp, linear: bool = False):
    """Potential vorticity at the positions of marked particles, frame by frame.

    Book: §13.13, Eq. (13.94): D/Dt[(ζ + f)/h] = 0 — each column of the result should be constant in time.
    Parameters: model; history from :func:`run`; xp, yp (n_frames, n_particles) from :func:`advect_particles`;
    linear — use the linearised potential vorticity.
    Returns q (n_frames, n_particles) [1/(m s)] (bilinear interpolation from the corners).
    Assumptions: particles away from walls.
    Validation (measured; nonlinear geostrophic vortex, 3 inertial periods, particles spanning a 55 % range of q):
    largest drift along a path 1.2 % on 32² cells and 0.5 % on 64² (it falls with resolution; the scheme conserves
    energy, not potential enstrophy, so the drift is truncation error of the model and of the interpolation).
    Label: conserved (approximately; numbers above).
    """
    out = []
    for i in range(len(history["t"])):
        q = potential_vorticity(model, dict(eta=history["eta"][i], u=history["u"][i], v=history["v"][i]), linear)
        out.append(interpolate(model, q, xp[i], yp[i], "corner"))
    return np.array(out)


# =====================================================================================================================
# The 1-D linear set (geostrophic adjustment on a line)
# =====================================================================================================================

def linear_1d_step(eta, u, v, dx: float, dt: float, H: float, f: float, g: float = G0):
    """One forward–backward step of the linear shallow-water equations with no variation in y.

    Book: §13.8, Eq. (13.45) with ∂/∂y = 0:  ∂η/∂t + H ∂u/∂x = 0,  ∂u/∂t − f v = −g ∂η/∂x,  ∂v/∂t + f u = 0.
    Grid: η and v at the n cell centres, u at the n + 1 faces; the two end faces are walls (u = 0).
    Scheme (ours):  η ← η − dt H δ_x u and v ← v − dt f ū (both with the old u);  then
    u ← u + dt (f v̄_new − g δ_x η_new) at interior faces, with v̄, ū two-point averages.  First order in time, second
    order in space; neutrally stable for gravity waves when c dt/dx ≤ 1.  Because η and v are advanced with the same
    u, the discrete linear potential vorticity δ_x v − f η̄/H at interior faces is conserved to round-off.
    Parameters: eta (n,) [m]; u (n + 1,) [m/s]; v (n,) [m/s]; dx [m]; dt [s]; H [m]; f [1/s]; g [m/s²].
    Returns (eta, u, v) after one step (new arrays).
    Raises ValueError if sqrt(gH) dt/dx > 1 or |f| dt > 1.
    Assumptions: linear, inviscid, f-plane.   Validation: V1 f = 0 pulse speed sqrt(gH); V3 converges to
    ``core.gfd.geostrophic_adjustment_1d`` (see :func:`linear_1d_run`).  Label: converged.
    """
    eta, u, v = _F(eta), _F(u), _F(v)
    if u.size != eta.size + 1 or v.size != eta.size:
        raise ValueError("linear_1d_step: need len(u) = len(eta) + 1 and len(v) = len(eta)")
    if np.sqrt(g * H) * dt / dx > 1.0 + 1e-12 or abs(f) * dt > 1.0:
        raise ValueError("linear_1d_step: unstable step — need sqrt(g H) dt/dx <= 1 and |f| dt <= 1")
    eta_n = eta - dt * H * (u[1:] - u[:-1]) / dx                                # Eq. (13.45a)
    v_n = v - dt * f * 0.5 * (u[1:] + u[:-1])                                   # Eq. (13.45c), same (old) u
    u_n = u.copy()
    u_n[1:-1] = u[1:-1] + dt * (f * 0.5 * (v_n[1:] + v_n[:-1]) - g * (eta_n[1:] - eta_n[:-1]) / dx)   # Eq. (13.45b)
    return eta_n, u_n, v_n


def linear_1d_run(eta0, dx: float, dt: float, n_steps: int, H: float, f: float, g: float = G0,
                  save_every: int = 1) -> dict:
    """March the 1-D linear shallow-water set from rest for ``n_steps`` steps.

    Book: §13.8, Eq. (13.45) with ∂/∂y = 0; scheme of :func:`linear_1d_step` (ours).
    Parameters: eta0 (n,) initial elevation at the cell centres [m]; dx [m]; dt [s]; n_steps; H [m]; f [1/s];
    g [m/s²]; save_every — keep every this many steps (first and last always kept).
    Returns dict: ``t`` (n_saved,), ``eta`` (n_saved, n), ``u`` (n_saved, n + 1), ``v`` (n_saved, n) and ``pv``
    (2, n − 1): the linear potential vorticity δ_x v − f η̄/H at the interior faces for the first and the last saved
    step (equal to round-off).
    Assumptions: as :func:`linear_1d_step`; the end walls reflect the gravity waves, so make the line long enough
    that they have not come back (travel time 2 × distance to the wall / sqrt(gH)).
    Validation (measured; step of half-height 0.5 m, equivalent depth 1.3 m, f at 35°, line 60 Λ long, run to
    t = 20 Λ/c): the instantaneous η within 5 Λ of the step still differs from the closed-form end state by 7–9 % of
    η₀ (inertia–gravity waves of near-zero group velocity decay only as t^{−1/2}); the average over the last inertial
    period agrees with it to 0.6 % of η₀ and 0.7 % of the jet speed (800 cells; 0.6 % and 0.9 % on 400); volume and
    potential vorticity are conserved to round-off (1e-18 1/s).  Label: converged.
    """
    eta = _F(eta0).copy()
    n = eta.size
    u, v = np.zeros(n + 1), np.zeros(n)

    def pv(e, vv):
        return (vv[1:] - vv[:-1]) / dx - f * 0.5 * (e[1:] + e[:-1]) / H

    every = max(int(save_every), 1)
    ts, es, us, vs = [0.0], [eta.copy()], [u.copy()], [v.copy()]
    pv0 = pv(eta, v)
    for k in range(1, int(n_steps) + 1):
        eta, u, v = linear_1d_step(eta, u, v, dx, dt, H, f, g)
        if k % every == 0 or k == int(n_steps):
            ts.append(k * dt)
            es.append(eta.copy())
            us.append(u.copy())
            vs.append(v.copy())
    return dict(t=np.array(ts), eta=np.array(es), u=np.array(us), v=np.array(vs), pv=np.array([pv0, pv(eta, v)]))


# =====================================================================================================================
# Linear quasi-geostrophic evolution (Rossby waves), exact in time
# =====================================================================================================================

def qg_linear_evolve(eta0, x, y, t, beta, f0, c):
    """Evolve a doubly periodic height field with the linear quasi-geostrophic vorticity equation, exactly in time.

    Book: §13.15, Eq. (13.117): ∂/∂t (∇²η − (f₀²/c²) η) + β ∂η/∂x = 0; each Fourier mode exp[i(kx + ly)] turns with
    its own frequency ω = −βk/(k² + l² + f₀²/c²), Eq. (13.118).
    Parameters
    ----------
    eta0 : initial field ``[ny, nx]`` on a uniform periodic grid [m].   x, y : the grid coordinates [m] (1-D, uniform;
    the periods are nx·dx and ny·dy).   t : time [s].   beta [1/(m s)].   f0 [1/s].   c [m/s] (``np.inf`` for the
    non-divergent barotropic wave; the mean of the field is then left unchanged).
    Returns
    -------
    eta(x, y, t) ``[ny, nx]`` [m].
    Numerics: ``numpy.fft``; exact for the modes resolved by the grid (no time-stepping error).
    Assumptions: quasi-geostrophic, linear, β-plane treated as periodic in y (valid for a packet well inside the box).
    Validation (measured): a single mode moves by ω t to 3e-15.  Label: analytic.
    """
    eta0 = _F(eta0)
    x, y = _F(x), _F(y)
    ny, nx = eta0.shape
    k = 2.0 * np.pi * np.fft.fftfreq(nx, d=x[1] - x[0])[None, :]
    l = 2.0 * np.pi * np.fft.fftfreq(ny, d=y[1] - y[0])[:, None]
    F = 0.0 if np.isinf(c) else (f0 / c) ** 2
    D = k * k + l * l + F
    om = np.where(D > 0, -beta * k / np.where(D > 0, D, 1.0), 0.0)   # Eq. (13.118)
    return np.real(np.fft.ifft2(np.fft.fft2(eta0) * np.exp(-1j * om * float(t))))


def qg_linear_evolve_1d(amplitudes, k, l, t, beta, f0, c, U: float = 0.0, x=None):
    """A Rossby-wave packet as a sum of plane waves, each advanced with its own frequency.

    Book: §13.15, Eq. (13.118) with the mean flow of Eq. (13.120): ω_j = U k_j − β k_j/(k_j² + l² + f₀²/c²); the field
    is η(x, t) = Re Σ_j a_j(t) exp(i k_j x) with a_j(t) = a_j exp(−i ω_j t).
    Parameters: amplitudes a_j (complex or real) [m]; k (same length) [rad/m]; l northward wavenumber [rad/m]; t [s];
    beta; f0; c; U mean flow [m/s]; x — optional positions [m].
    Returns the advanced complex amplitudes a_j(t) (an array); if ``x`` is given, the real field η at x instead
    (a float for scalar x).
    Assumptions: as :func:`qg_linear_evolve`.   Validation: V1 the envelope moves at the group velocity.
    Label: analytic.
    """
    from .gfd import rossby_omega

    a = np.atleast_1d(np.asarray(amplitudes, dtype=complex))
    kk = np.atleast_1d(_F(k))
    adv = a * np.exp(-1j * _F(rossby_omega(kk, l, beta, f0, c, U)) * float(t))
    if x is None:
        return adv
    return _S(np.real(np.exp(1j * np.multiply.outer(_F(x), kk)) @ adv))


# =====================================================================================================================
# Barotropic vorticity equation, pseudo-spectral (two-dimensional turbulence)
# =====================================================================================================================

def spectral_wavenumbers(n: int, L: float):
    """Wavenumber arrays (KX, KY) of an n × n periodic box of side L, in the layout of ``numpy.fft.rfft2``.

    Book: not in the book (ours).   Parameters: n points per side (even); L [m].
    Returns (KX, KY): KX shape (1, n/2 + 1) ≥ 0, KY shape (n, 1), both [rad/m].
    Assumptions: none.   Validation: V1 derivative of a sine.  Label: analytic.
    """
    return (2.0 * np.pi * np.fft.rfftfreq(int(n), d=L / n)[None, :], 2.0 * np.pi * np.fft.fftfreq(int(n), d=L / n)[:, None])


def _spec(zeta_hat, KX, KY):
    n = zeta_hat.shape[0]
    K2 = KX * KX + KY * KY
    return n, K2, np.where(K2 > 0, 1.0 / np.where(K2 > 0, K2, 1.0), 0.0)


def barotropic_velocity(zeta, L: float):
    """Velocity (u, v) of a doubly periodic non-divergent flow from its vorticity: ∇²ψ = ζ, u = −∂ψ/∂y, v = ∂ψ/∂x.

    Book: §13.16 (the stream-function convention of the chapter, u′ = −∂ψ/∂y, v′ = ∂ψ/∂x — opposite to ch04, ch11).
    Parameters: zeta ``[n, n]`` [1/s] on a square periodic box of side L [m].   Returns (u, v) [m/s].
    Numerics: spectral (exact for resolved modes).   Assumptions: zero mean flow.   Validation: V1 a single mode.
    Label: analytic.
    """
    zeta = _F(zeta)
    n = zeta.shape[0]
    KX, KY = spectral_wavenumbers(n, L)
    zh = np.fft.rfft2(zeta)
    _, _, K2inv = _spec(zh, KX, KY)
    psi = -zh * K2inv
    return np.fft.irfft2(-1j * KY * psi, s=(n, n)), np.fft.irfft2(1j * KX * psi, s=(n, n))


def barotropic_vorticity_rhs(zeta_hat, KX, KY, beta: float = 0.0, nu: float = 0.0, dealias: bool = True):
    """Right-hand side of the barotropic vorticity equation in spectral space: the transform of
    −J(ψ, ζ) − β ∂ψ/∂x + ν ∇²ζ.

    Book: §13.16, Eq. (13.122) (∂/∂t + u·∇)(ζ + f) = 0, i.e. ∂ζ/∂t = −u ∂ζ/∂x − v ∂ζ/∂y − β v; the viscous term is
    ours.
    Parameters: zeta_hat = ``rfft2(ζ)`` (shape (n, n/2 + 1)); KX, KY from :func:`spectral_wavenumbers`; beta
    [1/(m s)]; nu [m²/s]; dealias — apply the 2/3 rule to the quadratic term.
    Returns d(zeta_hat)/dt (same shape).
    Numerics: pseudo-spectral — derivatives in spectral space, the product in physical space.
    Assumptions: two-dimensional, non-divergent, doubly periodic.   Validation: V4 inviscid invariants.
    Label: conserved.
    """
    n, K2, K2inv = _spec(zeta_hat, KX, KY)
    psi = -zeta_hat * K2inv
    u = np.fft.irfft2(-1j * KY * psi, s=(n, n))
    v = np.fft.irfft2(1j * KX * psi, s=(n, n))
    zx = np.fft.irfft2(1j * KX * zeta_hat, s=(n, n))
    zy = np.fft.irfft2(1j * KY * zeta_hat, s=(n, n))
    adv = np.fft.rfft2(u * zx + v * zy)
    if dealias:
        kcut = 2.0 / 3.0 * np.max(np.abs(KY))       # strict '<': a mode at exactly 2/3 of Nyquist would still be aliased
        adv = adv * ((np.abs(KX) < kcut * (1 - 1e-12)) & (np.abs(KY) < kcut * (1 - 1e-12)))
    return -adv - beta * 1j * KX * psi - nu * K2 * zeta_hat  # Eq. (13.122)


def barotropic_invariants(zeta, L: float) -> dict:
    """Energy and enstrophy of a doubly periodic two-dimensional flow (domain means).

    Book: §13.18, Eqs. (13.143)–(13.144) (both are conserved by the inviscid flow).
    Parameters: zeta ``[n, n]`` [1/s]; L box side [m].
    Returns dict: ``energy`` = ½ mean(u² + v²) [m²/s²] and ``enstrophy`` = mean(ζ²) [1/s²] (the book's definition:
    the mean-square vorticity, with no factor ½).
    Assumptions: zero mean flow.   Validation: V1 a single mode.  Label: analytic.
    """
    u, v = barotropic_velocity(zeta, L)
    return dict(energy=0.5 * float(np.mean(u * u + v * v)), enstrophy=float(np.mean(_F(zeta) ** 2)))


def _mode_energy(zeta, L):
    zeta = _F(zeta)
    n = zeta.shape[0]
    KX, KY = spectral_wavenumbers(n, L)
    zh = np.fft.rfft2(zeta) / n ** 2
    _, K2, K2inv = _spec(zh, KX, KY)
    w = np.where((KX > 0) & (KX < np.pi * n / L), 2.0, 1.0) * np.ones_like(K2)   # weights of the rfft half-plane
    Z = w * np.abs(zh) ** 2                 # mean-square vorticity per mode
    return K2, 0.5 * Z * K2inv, Z           # K^2, energy per mode, enstrophy per mode


def spectral_centroids(zeta, L: float) -> dict:
    """Mean wavenumbers of the energy and of the enstrophy: K_E = Σ K E_k / Σ E_k and K_Z = Σ K Z_k / Σ Z_k.

    Book: §13.18 (Fjørtoft's argument: energy moves to small wavenumbers, enstrophy to large ones) — the centroids
    are our diagnostic of that statement.
    Parameters: zeta ``[n, n]`` [1/s]; L [m].   Returns dict(K_E, K_Z) [rad/m].
    Assumptions: doubly periodic.   Validation: V1 a single mode gives its own K twice.  Label: analytic.
    """
    K2, E, Z = _mode_energy(zeta, L)
    K = np.sqrt(K2)
    return dict(K_E=float(np.sum(K * E) / np.sum(E)), K_Z=float(np.sum(K * Z) / np.sum(Z)))


def zonal_energy_fraction(zeta, L: float) -> float:
    """Fraction of the kinetic energy held by the zonal-mean flow (the modes with k_x = 0) — how far a β-plane flow
    has organised into east–west jets (small for an isotropic field, 1 for pure jets).

    Book: §13.18 (the anisotropic, zonally elongated end state on a β-plane); the index is ours.
    Parameters: zeta ``[n, n]`` [1/s]; L [m].   Returns the fraction (dimensionless).   Assumptions: doubly periodic.
    Validation: V1 a purely zonal flow gives 1.  Label: analytic.
    """
    _, E, _ = _mode_energy(zeta, L)
    return float(np.sum(E[:, 0]) / np.sum(E))


def barotropic_spectrum(zeta, L: float):
    """Isotropic kinetic-energy spectrum E(K) of a two-dimensional flow, binned in shells of width 2π/L.

    Book: §13.18 (the energy spectrum S(K)).  Normalisation here: Σ E ΔK = ½ mean(u² + v²), the convention of
    ``core.turbstats.shell_spectrum`` — the book writes mean(u²) = ∫S dK with no factor ½ (trap T16), so say which
    one is plotted.
    Parameters: zeta ``[n, n]`` [1/s]; L [m].   Returns (K, E): shell centres [rad/m] and E(K) [m³/s²].
    Assumptions: statistically isotropic.   Validation: V1 Σ E ΔK = energy of :func:`barotropic_invariants`
    (measured equal to 16 digits).  Label: analytic.
    """
    from .turbstats import shell_spectrum

    return shell_spectrum(barotropic_velocity(zeta, L), L)


def random_vorticity(n: int, L: float, k_peak: float = 8.0, u_rms: float = 1.0, seed: int = 0):
    """A seeded random vorticity field with its energy in a ring of wavenumbers — an initial condition for
    :func:`barotropic_run`.

    Book: not in the book (ours).  The curl of ``core.turbstats.synthetic_solenoidal_field`` (kinematic: random
    phases, no dynamics yet) with energy spectrum ∝ (K/K₀)⁴ exp(−2(K/K₀)²), K₀ = k_peak·2π/L, truncated at 2/3 of the
    Nyquist wavenumber and rescaled so that sqrt(mean(u² + v²)) = u_rms.
    Parameters: n points per side (even); L [m]; k_peak (in units of 2π/L); u_rms [m/s]; seed of
    ``numpy.random.default_rng``.
    Returns zeta ``[n, n]`` [1/s].
    Assumptions: none.   Validation: V1 speed_rms; reproducible for a given seed.  Label: analytic.
    """
    from .turbstats import synthetic_solenoidal_field

    K0 = k_peak * 2.0 * np.pi / L
    u, v = synthetic_solenoidal_field(int(n), L, lambda K: (K / K0) ** 4 * np.exp(-2.0 * (K / K0) ** 2), seed=seed, dim=2)
    KX, KY = spectral_wavenumbers(n, L)
    kcut = 2.0 / 3.0 * np.pi * n / L
    zh = (1j * KX * np.fft.rfft2(v) - 1j * KY * np.fft.rfft2(u)) * ((np.abs(KX) < kcut * (1 - 1e-12)) & (np.abs(KY) < kcut * (1 - 1e-12)))
    zeta = np.fft.irfft2(zh, s=(int(n), int(n)))
    return zeta * u_rms / np.sqrt(2.0 * barotropic_invariants(zeta, L)["energy"])


def barotropic_time_step(zeta0, L: float, beta: float = 0.0, nu: float = 0.0, cfl: float = 0.4) -> float:
    """A safe time step for :func:`barotropic_run`: ``cfl``·dx/max|u|, also below 0.5/(β L/2π) and 0.5/(ν K_max²).

    Book: not in the book (ours).   Parameters: zeta0 ``[n, n]`` [1/s]; L [m]; beta; nu; cfl.   Returns dt [s].
    Assumptions: the flow does not speed up (true for decaying turbulence).   Validation: V3 the run is stable.
    Label: converged.
    """
    zeta0 = _F(zeta0)
    n = zeta0.shape[0]
    u, v = barotropic_velocity(zeta0, L)
    umax = float(np.max(np.hypot(u, v)))
    dt = cfl * (L / n) / umax if umax > 0 else np.inf
    if beta:
        dt = min(dt, 0.5 * 2.0 * np.pi / (abs(beta) * L))
    if nu:
        dt = min(dt, 0.5 / (nu * (np.pi * n / L) ** 2))
    if not np.isfinite(dt):
        raise ValueError("barotropic_time_step: the field is at rest and beta = nu = 0 — choose dt yourself")
    return float(dt)


def barotropic_run(n: int, L: float, zeta0, *, t_end: float, dt: float, beta: float = 0.0, nu: float = 0.0,
                   dealias: bool = True, save_every: int = 1) -> dict:
    """Run the barotropic vorticity equation on a doubly periodic box (freely evolving two-dimensional turbulence).

    Book: §13.16, Eq. (13.122), and §13.18 (inverse energy cascade, forward enstrophy cascade, zonal jets on a
    β-plane).  **The model is ours — not in the book.**
    Parameters
    ----------
    n : grid points per side (even, ≥ 8).   L : box side [m].   zeta0 : initial vorticity ``[n, n]`` [1/s] (e.g.
    :func:`random_vorticity`).
    t_end [s], dt [s] : keyword-only, required (:func:`barotropic_time_step` suggests dt); dt is shortened so that an
    integer number of steps reaches t_end.
    beta [1/(m s)].   nu : viscosity [m²/s] (0 = inviscid).   dealias : 2/3 rule.
    save_every : keep every this many steps (the first and last are always kept).
    Returns
    -------
    dict: ``t`` (n_saved,), ``zeta`` (n_saved, n, n) and, per saved step, ``energy`` (½ mean(u² + v²)), ``enstrophy``
    (mean ζ²), ``K_E`` and ``K_Z`` (the centroids of :func:`spectral_centroids`).
    Numerics (ours): pseudo-spectral in space, classical RK4 in time (fourth order), 2/3 de-aliasing of the quadratic
    term.  Inviscid and de-aliased, the truncated equations conserve energy and enstrophy; RK4 adds a drift ∝ dt⁴
    that acts mostly on the highest wavenumbers kept.
    Raises ValueError if dt exceeds the advective limit dx/max|u| of the initial field.
    Assumptions: two-dimensional, non-divergent, doubly periodic, β-plane without boundaries.
    Validation (measured; 64², :func:`random_vorticity` with seed 1, unit u_rms, box 2π, dt from
    :func:`barotropic_time_step`): inviscid to t = 2: energy drift −1.1e-5, enstrophy drift −5.2e-5 (−3.3e-7 and −1.7e-6 with dt halved); a single Rossby
    mode propagates at −βk/K² to 2e-13 of its amplitude; with ν = 2e-4 to t = 20: energy falls to 0.71, enstrophy to
    0.14 of their initial values, K_E from 8.5 to 2.4 (energy to large scales), K_Z from 10.2 to 13.5.
    Label: conserved (invariants), qualitative (cascade).
    """
    n = int(n)
    if n % 2 or n < 8:
        raise ValueError("barotropic_run: n must be even and >= 8")
    zeta = _F(zeta0).copy()
    if zeta.shape != (n, n):
        raise ValueError("barotropic_run: zeta0 must have shape (n, n)")
    KX, KY = spectral_wavenumbers(n, L)
    u, v = barotropic_velocity(zeta, L)
    umax = float(np.max(np.hypot(u, v)))
    if umax > 0 and dt > (L / n) / umax:
        raise ValueError(f"barotropic_run: dt = {dt:.3g} exceeds the advective limit dx/max|u| = {(L / n) / umax:.3g}")
    nsteps = max(int(np.ceil(t_end / dt - 1e-9)), 1)
    dt = t_end / nsteps
    every = max(int(save_every), 1)
    zh = np.fft.rfft2(zeta)
    rhs = lambda q: barotropic_vorticity_rhs(q, KX, KY, beta, nu, dealias)  # noqa: E731
    ts, zs, es, ens, kes, kzs = [], [], [], [], [], []

    def keep(tt, q):
        z = np.fft.irfft2(q, s=(n, n))
        inv, cen = barotropic_invariants(z, L), spectral_centroids(z, L)
        ts.append(tt), zs.append(z), es.append(inv["energy"]), ens.append(inv["enstrophy"])
        kes.append(cen["K_E"]), kzs.append(cen["K_Z"])

    keep(0.0, zh)
    for k in range(1, nsteps + 1):
        k1 = rhs(zh)
        k2 = rhs(zh + 0.5 * dt * k1)
        k3 = rhs(zh + 0.5 * dt * k2)
        k4 = rhs(zh + dt * k3)
        zh = zh + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if k % every == 0 or k == nsteps:
            keep(k * dt, zh)
    return dict(t=np.array(ts), zeta=np.array(zs), energy=np.array(es), enstrophy=np.array(ens), K_E=np.array(kes),
                K_Z=np.array(kzs))
