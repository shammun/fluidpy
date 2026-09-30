"""The MAC (marker-and-cell) projection method on a staggered grid for the non-dimensional incompressible Navier–Stokes
equations: explicit convection–diffusion predictor, discrete pressure Poisson equation, velocity correction.

Book: Kundu, Cohen & Dowling 5e, Ch. 10 §10.4, Eqs. (10.81), (10.111)–(10.128) (operator splitting, MAC split (10.115),
predictor (10.116), projection (10.117)–(10.118), staggered grid Fig. 10.4, (10.119)–(10.124) discrete equations,
(10.125)–(10.126) collocated vs staggered gradients, (10.127)–(10.128) stability).  Pages chapters/pages/ch10/p469–p473.

Grid (cell-centred, staggered — the Arakawa C-grid of ocean/atmosphere models): cells of size Δx × Δy; the pressure
``p[j, i]`` sits at the cell centre ((i + ½)Δx, (j + ½)Δy); ``u[j, i]`` at the vertical face x = iΔx, y = (j + ½)Δy;
``v[j, i]`` at the horizontal face x = (i + ½)Δx, y = jΔy.  Shapes: p (ny, nx); u (ny, nx + 1) with walls in x or (ny, nx)
periodic; v (ny + 1, nx) with walls in y or (ny, nx) periodic.  The book's half indices map as u_{i+1/2, j} ↔ u[j, i + 1]
(0-based cells).  Walls: the normal velocity lives on the wall faces (u[:, 0], u[:, nx], v[0, :], v[ny, :]) and is a
boundary condition, never corrected — which is why the Poisson equation needs no pressure boundary condition (the book's
argument, printed p. 445; slip R4: "(10.120)" there should read "(10.124)").  Tangential no-slip needs values half a cell
outside the wall: ghost values by linear (u_g = 2U − u₀) or quadratic extrapolation — a choice the book leaves open
(# DEVIATION note at ``_pad_u``).
Everything is non-dimensional (lengths L, velocity U, time L/U, pressure ρU², Re = UL/ν as in (10.81)).
Reused by Ch. 11 (nonlinear saturation, Rayleigh–Bénard), Ch. 12 (2-D decaying turbulence), Ch. 13 (C-grid shallow water).
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

__all__ = [
    "MacGrid", "face_coordinates", "new_state", "divergence", "gradient", "curl", "vorticity", "convective_terms",
    "laplacian_faces", "predictor", "pressure_poisson_matrix", "solve_pressure", "pin_pressure", "correct", "project",
    "projection_stages", "mac_projection_summary", "dt_limit", "step", "run", "cavity", "taylor_green",
    "channel_poiseuille", "streamfunction", "primary_vortex_centre", "cavity_centreline", "cavity_centreline_v",
]


@dataclass
class MacGrid:
    """Staggered MAC grid of nx × ny cells on [0, Lx] × [0, Ly] (non-dimensional).

    Book: §10.4, Fig. 10.4.  ``periodic`` = (in x, in y).  ``walls`` = tangential wall speeds {"top": u, "bottom": u,
    "left": v, "right": v} (a lid-driven cavity has top = 1; :func:`step`/:func:`run` take a ``lid`` argument too); the
    normal wall speed is 0 (no through-flow).  ``ghost`` = "linear" | "quadratic" wall extrapolation for the tangential
    velocity (our choice).  Attributes dx, dy, xc, yc (cell centres), xu (u-face x), yv (v-face y), u_shape, v_shape;
    ``zeros()`` gives a state at rest.  Label: analytic.
    """
    nx: int
    ny: int
    Lx: float = 1.0
    Ly: float = 1.0
    periodic: tuple = (False, False)
    walls: dict = field(default_factory=lambda: dict(top=0.0, bottom=0.0, left=0.0, right=0.0))
    ghost: str = "linear"

    @property
    def dx(self) -> float:
        return self.Lx / self.nx

    @property
    def dy(self) -> float:
        return self.Ly / self.ny

    @property
    def u_shape(self):
        return (self.ny, self.nx if self.periodic[0] else self.nx + 1)

    @property
    def v_shape(self):
        return (self.ny if self.periodic[1] else self.ny + 1, self.nx)

    @property
    def xc(self):
        return (np.arange(self.nx) + 0.5) * self.dx

    @property
    def yc(self):
        return (np.arange(self.ny) + 0.5) * self.dy

    @property
    def xu(self):
        return np.arange(self.u_shape[1]) * self.dx

    @property
    def yv(self):
        return np.arange(self.v_shape[0]) * self.dy

    def zeros(self) -> dict:
        return new_state(self)

    def key(self):
        return (self.nx, self.ny, float(self.Lx), float(self.Ly), tuple(self.periodic))


def face_coordinates(g: MacGrid) -> dict:
    """Coordinates of every staggered unknown: dict(xp, yp (cell centres), xu, yu (u faces), xv, yv (v faces)) as 2-D
    meshgrids in the project layout [j, i].  Book: Fig. 10.4.  Label: analytic."""
    XP, YP = np.meshgrid(g.xc, g.yc, indexing="xy")
    XU, YU = np.meshgrid(g.xu, g.yc, indexing="xy")
    XV, YV = np.meshgrid(g.xc, g.yv, indexing="xy")
    return dict(xp=XP, yp=YP, xu=XU, yu=YU, xv=XV, yv=YV)


def new_state(g: MacGrid) -> dict:
    """Fluid at rest: dict(u, v, p, t = 0).  Label: analytic."""
    return dict(u=np.zeros(g.u_shape), v=np.zeros(g.v_shape), p=np.zeros((g.ny, g.nx)), t=0.0)


def _ghost(inner0, inner1, wall, kind):
    if kind == "linear":
        return 2.0 * wall - inner0
    return (8.0 * wall - 6.0 * inner0 + inner1) / 3.0  # quadratic through (wall, inner0, inner1)


def _pad_u(u, g: MacGrid):
    """u with one ghost layer all round: shape (ny + 2, nxu + 2).

    # DEVIATION: tangential no-slip by ghost values half a cell outside the wall (linear u_g = 2U − u₀ or quadratic
    extrapolation) — the book says only "the functions interpolated at the grid locations" (10.119)–(10.120).
    """
    ny, nxu = u.shape
    ue = np.zeros((ny + 2, nxu + 2))
    ue[1:-1, 1:-1] = u
    if g.periodic[0]:
        ue[1:-1, 0], ue[1:-1, -1] = u[:, -1], u[:, 0]
    if g.periodic[1]:
        ue[0, 1:-1], ue[-1, 1:-1] = u[-1, :], u[0, :]
    else:
        ue[0, 1:-1] = _ghost(u[0], u[1], g.walls.get("bottom", 0.0), g.ghost)
        ue[-1, 1:-1] = _ghost(u[-1], u[-2], g.walls.get("top", 0.0), g.ghost)
    if g.periodic[0]:
        ue[0, 0], ue[0, -1], ue[-1, 0], ue[-1, -1] = ue[0, -2], ue[0, 1], ue[-1, -2], ue[-1, 1]
    return ue


def _pad_v(v, g: MacGrid):
    """v with one ghost layer all round: shape (nyv + 2, nx + 2) (tangential ghosts on the left/right walls)."""
    nyv, nx = v.shape
    ve = np.zeros((nyv + 2, nx + 2))
    ve[1:-1, 1:-1] = v
    if g.periodic[1]:
        ve[0, 1:-1], ve[-1, 1:-1] = v[-1, :], v[0, :]
    if g.periodic[0]:
        ve[1:-1, 0], ve[1:-1, -1] = v[:, -1], v[:, 0]
    else:
        ve[1:-1, 0] = _ghost(v[:, 0], v[:, 1], g.walls.get("left", 0.0), g.ghost)
        ve[1:-1, -1] = _ghost(v[:, -1], v[:, -2], g.walls.get("right", 0.0), g.ghost)
    if g.periodic[1]:
        ve[0, 0], ve[-1, 0], ve[0, -1], ve[-1, -1] = ve[-2, 0], ve[1, 0], ve[-2, -1], ve[1, -1]
    return ve


def _u_unknown_mask(g: MacGrid):
    m = np.ones(g.u_shape, dtype=bool)
    if not g.periodic[0]:
        m[:, 0] = m[:, -1] = False
    return m


def _v_unknown_mask(g: MacGrid):
    m = np.ones(g.v_shape, dtype=bool)
    if not g.periodic[1]:
        m[0, :] = m[-1, :] = False
    return m


def _with_lid(g: MacGrid, lid):
    if lid is None or g.periodic[1] or float(lid) == float(g.walls.get("top", 0.0)):
        return g
    w = dict(g.walls)
    w["top"] = float(lid)
    return MacGrid(g.nx, g.ny, g.Lx, g.Ly, g.periodic, w, g.ghost)


# ----------------------------------------------------------------------------------------------------------------------
# discrete operators
# ----------------------------------------------------------------------------------------------------------------------
def divergence(u, v, g: MacGrid) -> np.ndarray:
    """Discrete divergence of each cell (net outflow per unit area), shape (ny, nx).

    Book: §10.4, Eq. (10.123): (u_{i+1/2,j} − u_{i−1/2,j})/Δx + (v_{i,j+1/2} − v_{i,j−1/2})/Δy (second order with four values).
    Label: analytic, conserved.
    """
    du = (np.roll(u, -1, axis=1) - u) if g.periodic[0] else (u[:, 1:] - u[:, :-1])
    dv = (np.roll(v, -1, axis=0) - v) if g.periodic[1] else (v[1:, :] - v[:-1, :])
    return du / g.dx + dv / g.dy  # Eq. (10.123)


def gradient(p, g: MacGrid):
    """Pressure gradient on the faces: ((p_{i+1,j} − p_{i,j})/Δx on u faces, (p_{i,j+1} − p_{i,j})/Δy on v faces).

    Book: §10.4, Eq. (10.126) (second order at the faces).  Wall faces get 0 (they are never corrected).
    Returns (gx of u-shape, gy of v-shape).  Label: analytic.
    """
    gx = np.zeros(g.u_shape)
    gy = np.zeros(g.v_shape)
    if g.periodic[0]:
        gx[:] = (p - np.roll(p, 1, axis=1)) / g.dx  # Eq. (10.126)
    else:
        gx[:, 1:-1] = (p[:, 1:] - p[:, :-1]) / g.dx  # Eq. (10.126)
    if g.periodic[1]:
        gy[:] = (p - np.roll(p, 1, axis=0)) / g.dy
    else:
        gy[1:-1, :] = (p[1:, :] - p[:-1, :]) / g.dy
    return gx, gy


def curl(u, v, g: MacGrid, interior_only: bool = False) -> np.ndarray:
    """Vorticity ω = ∂v/∂x − ∂u/∂y at the cell corners (x = iΔx, y = jΔy).

    Book: §10.4 (the correction −Δt∇p of (10.117) is irrotational: its discrete curl vanishes identically — the
    "irrotational correction field" of the projection method, p. 445).  Shape (ny + 1, nx + 1) with walls (wall corners from
    the ghost values of the grid), (ny, nx) periodic; ``interior_only`` → only corners with four real neighbours.
    Label: analytic.
    """
    ue, ve = _pad_u(u, g), _pad_v(v, g)
    nyv, nx = v.shape
    ny, nxu = u.shape
    ncx = nx if g.periodic[0] else nx + 1
    ncy = ny if g.periodic[1] else ny + 1
    vx = (ve[1:ncy + 1, 1:ncx + 1] - ve[1:ncy + 1, 0:ncx]) / g.dx
    uy = (ue[1:ncy + 1, 1:ncx + 1] - ue[0:ncy, 1:ncx + 1]) / g.dy
    w = vx - uy
    if interior_only:
        return w[slice(None) if g.periodic[1] else slice(1, -1), slice(None) if g.periodic[0] else slice(1, -1)]
    return w


def vorticity(u, v, g: MacGrid) -> np.ndarray:
    """ω = ∂v/∂x − ∂u/∂y at the cell corners (alias of :func:`curl`).  Book: §10.5 (vorticity pictures).  Label: analytic."""
    return curl(u, v, g)


def laplacian_faces(u, v, g: MacGrid):
    """Five-point Laplacians of u and v at their own faces (with the wall ghosts).  Book: §10.4 (10.119)–(10.120).
    Returns (lap_u, lap_v).  Label: analytic, converged."""
    ue, ve = _pad_u(u, g), _pad_v(v, g)
    lu = (ue[1:-1, 2:] - 2 * u + ue[1:-1, :-2]) / g.dx ** 2 + (ue[2:, 1:-1] - 2 * u + ue[:-2, 1:-1]) / g.dy ** 2
    lv = (ve[1:-1, 2:] - 2 * v + ve[1:-1, :-2]) / g.dx ** 2 + (ve[2:, 1:-1] - 2 * v + ve[:-2, 1:-1]) / g.dy ** 2
    return lu, lv


def convective_terms(u, v, g: MacGrid, form: str = "advective"):
    """Convective acceleration (u·∇)u at the u faces and v faces.

    Book: §10.4, (10.119)–(10.120) (u ∂u/∂x + v ∂u/∂y "interpolated at the grid locations"), (10.82) (the conservative form
    ∇·(uu), equal when ∇·u = 0).  ``form="advective"``: centred differences with v averaged from its four neighbouring faces
    (and u likewise at v faces) — our interpolation choice; ``form="conservative"``: Harlow–Welch flux form with products
    formed from two-point averages.  Returns (cu of u-shape, cv of v-shape); only unknown faces are meaningful.
    Label: analytic, converged.
    """
    ue, ve = _pad_u(u, g), _pad_v(v, g)
    ny, nxu = u.shape
    nyv, nx = v.shape
    dx, dy = g.dx, g.dy
    if form == "advective":
        vbar = 0.25 * (ve[1:ny + 1, 0:nxu] + ve[1:ny + 1, 1:nxu + 1] + ve[2:ny + 2, 0:nxu] + ve[2:ny + 2, 1:nxu + 1])
        ux = (ue[1:-1, 2:] - ue[1:-1, :-2]) / (2 * dx)
        uy = (ue[2:, 1:-1] - ue[:-2, 1:-1]) / (2 * dy)
        cu = u * ux + vbar * uy  # (u ∂u/∂x + v ∂u/∂y) at (i + 1/2, j), Eq. (10.119)
        ubar = 0.25 * (ue[0:nyv, 1:nx + 1] + ue[0:nyv, 2:nx + 2] + ue[1:nyv + 1, 1:nx + 1] + ue[1:nyv + 1, 2:nx + 2])
        vx = (ve[1:-1, 2:] - ve[1:-1, :-2]) / (2 * dx)
        vy = (ve[2:, 1:-1] - ve[:-2, 1:-1]) / (2 * dy)
        cv = ubar * vx + v * vy  # (u ∂v/∂x + v ∂v/∂y) at (i, j + 1/2), Eq. (10.120)
        return cu, cv
    if form == "conservative":
        uc_r = 0.5 * (ue[1:-1, 1:-1] + ue[1:-1, 2:])
        uc_l = 0.5 * (ue[1:-1, :-2] + ue[1:-1, 1:-1])
        u_top = 0.5 * (ue[1:-1, 1:-1] + ue[2:, 1:-1])
        u_bot = 0.5 * (ue[:-2, 1:-1] + ue[1:-1, 1:-1])
        v_top = 0.5 * (ve[2:ny + 2, 0:nxu] + ve[2:ny + 2, 1:nxu + 1])
        v_bot = 0.5 * (ve[1:ny + 1, 0:nxu] + ve[1:ny + 1, 1:nxu + 1])
        cu = (uc_r ** 2 - uc_l ** 2) / dx + (u_top * v_top - u_bot * v_bot) / dy  # ∇·(uu), x component, Eq. (10.82)
        vc_t = 0.5 * (ve[1:-1, 1:-1] + ve[2:, 1:-1])
        vc_b = 0.5 * (ve[:-2, 1:-1] + ve[1:-1, 1:-1])
        v_r = 0.5 * (ve[1:-1, 1:-1] + ve[1:-1, 2:])
        v_l = 0.5 * (ve[1:-1, :-2] + ve[1:-1, 1:-1])
        u_r = 0.5 * (ue[0:nyv, 2:nx + 2] + ue[1:nyv + 1, 2:nx + 2])
        u_l = 0.5 * (ue[0:nyv, 1:nx + 1] + ue[1:nyv + 1, 1:nx + 1])
        cv = (u_r * v_r - u_l * v_l) / dx + (vc_t ** 2 - vc_b ** 2) / dy
        return cu, cv
    raise ValueError("form must be advective or conservative")


def predictor(u, v, g: MacGrid, Re: float, dt: float, lid: float | None = None, body=None, form: str = "advective"):
    """Explicit convection–diffusion step to the intermediate velocity u^{n+1/2}.

    Book: §10.4, Eq. (10.116) (u^{n+1/2} − u^n)/Δt + (u^n·∇)u^n − (1/Re)∇²u^n = g^{n+1}, in staggered form (10.119)–(10.120).
    Parameters: u, v (faces), g (grid), Re, dt, lid (top-wall speed; None → the grid's ``walls["top"]``), body = (f, g) body
    force per unit mass (scalars or face arrays), form.
    Returns (u*, v*); wall faces keep their boundary values.  Stable if (10.127)–(10.128) hold (:func:`dt_limit`).
    Label: converged.
    """
    g = _with_lid(g, lid)
    cu, cv = convective_terms(u, v, g, form)
    lu, lv = laplacian_faces(u, v, g)
    fx, fy = (0.0, 0.0) if body is None else body
    us = u + dt * (-cu + lu / Re + fx)  # Eq. (10.119)
    vs = v + dt * (-cv + lv / Re + fy)  # Eq. (10.120)
    return np.where(_u_unknown_mask(g), us, u), np.where(_v_unknown_mask(g), vs, v)


# ----------------------------------------------------------------------------------------------------------------------
# pressure Poisson equation (10.124)
# ----------------------------------------------------------------------------------------------------------------------
def _lap1d(n: int, h: float, periodic: bool):
    """1-D cell-centre Laplacian built link by link: each real neighbour j of cell i adds (p_j − p_i)/h².

    Wall faces contribute no link (normal velocity known — Neumann built in), so a wall cell has one fewer neighbour
    and a single cell between two walls has none (row of zeros). Periodic links wrap; with n = 2 both links of a cell
    reach the same neighbour (entry 2), with n = 1 they reach the cell itself and cancel. Every row sums to zero.
    """
    import scipy.sparse as sps

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for i in range(n):
        for j in (i - 1, i + 1):
            if periodic:
                j %= n
            elif j < 0 or j >= n:
                continue  # wall face: the outside neighbour drops out
            rows += [i, i]
            cols += [i, j]
            vals += [-1.0, 1.0]  # link i–j: −p_i + p_j (duplicates are summed by the COO → CSR conversion)
    A = sps.coo_matrix((vals, (rows, cols)), shape=(n, n))
    return (A / h ** 2).tocsr()


def pressure_poisson_matrix(g: MacGrid, pin: str | int | None = None):
    """Sparse matrix of the discrete Laplacian ∇²_d p of Eq. (10.124) on the cell centres (row index k = j·nx + i).

    Book: §10.4, Eq. (10.124) (from (10.121)–(10.123) — our D20); next to a wall the neighbour outside the domain is absent
    because the wall-face velocity is known and is not expressed through (10.121) (no pressure boundary condition).
    ``pin``: None → the singular matrix (constants in its null space: pure Neumann or periodic); "corner" (or a row
    index k) → that row replaced by p_k = 0 (the gauge), which makes it invertible.
    Built as kron(I_y, L_x) + kron(L_y, I_x) (scipy.sparse).  Label: analytic.
    """
    import scipy.sparse as sps

    Lx = _lap1d(g.nx, g.dx, g.periodic[0])
    Ly = _lap1d(g.ny, g.dy, g.periodic[1])
    A = (sps.kron(sps.identity(g.ny), Lx) + sps.kron(Ly, sps.identity(g.nx))).tocsr()
    if pin is None:
        return A
    k0 = 0 if pin == "corner" else int(pin)
    A = A.tolil()
    A[k0, :] = 0.0
    A[k0, k0] = 1.0
    return A.tocsr()


@lru_cache(maxsize=16)
def _factor(key, pin_index: int):
    from scipy.sparse.linalg import splu

    nx, ny, Lx, Ly, periodic = key
    return splu(pressure_poisson_matrix(MacGrid(nx, ny, Lx, Ly, periodic), pin=pin_index).tocsc())


def pin_pressure(p, where: str = "mean"):
    """Fix the arbitrary constant in p: ``"mean"`` → zero mean, ``"corner"`` → p[0, 0] = 0.

    Book: §10.4, after (10.83): "specify the value of the pressure at one reference point … determined up to a constant".
    Label: analytic.
    """
    p = np.asarray(p, dtype=float)
    return p - (p.mean() if where == "mean" else p[0, 0])


def _fft_poisson(rhs, g: MacGrid):
    ny, nx = rhs.shape
    lam = (-(4.0 / g.dx ** 2) * np.sin(np.pi * np.arange(nx) / nx) ** 2)[None, :] + \
          (-(4.0 / g.dy ** 2) * np.sin(np.pi * np.arange(ny) / ny) ** 2)[:, None]
    lam[0, 0] = 1.0
    ph = np.fft.fft2(rhs) / lam
    ph[0, 0] = 0.0
    return np.real(np.fft.ifft2(ph))


def _sor_neumann(rhs, g: MacGrid, omega: float = 1.8, tol: float = 1e-10, max_iter: int = 20000, p0=None,
                 n_sweeps: int | None = None):
    """Red-black SOR for the pure-Neumann/periodic 5-point problem (teaching option; the ch06 relaxation idea)."""
    ny, nx = rhs.shape
    p = np.zeros_like(rhs) if p0 is None else np.array(p0, dtype=float)
    ax, ay = 1.0 / g.dx ** 2, 1.0 / g.dy ** 2
    cx = np.full((ny, nx), 2.0)
    cy = np.full((ny, nx), 2.0)
    if not g.periodic[0]:
        cx[:, 0] -= 1
        cx[:, -1] -= 1
    if not g.periodic[1]:
        cy[0, :] -= 1
        cy[-1, :] -= 1
    diag = ax * cx + ay * cy
    J, I = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    red = (I + J) % 2 == 0
    A = pressure_poisson_matrix(g)
    it = 0
    for it in range(1, (n_sweeps or max_iter) + 1):
        for colour in (red, ~red):
            pe = np.pad(p, 1, mode="wrap")
            nb = ax * (pe[1:-1, 2:] + pe[1:-1, :-2]) + ay * (pe[2:, 1:-1] + pe[:-2, 1:-1])
            if not g.periodic[0]:
                nb[:, 0] -= ax * pe[1:-1, 0]
                nb[:, -1] -= ax * pe[1:-1, -1]
            if not g.periodic[1]:
                nb[0, :] -= ay * pe[0, 1:-1]
                nb[-1, :] -= ay * pe[-1, 1:-1]
            pgs = (nb - rhs) / diag
            p = np.where(colour, (1 - omega) * p + omega * pgs, p)
        p -= p.mean()
        if n_sweeps is None and it % 20 == 0:
            res = rhs - (A @ p.ravel()).reshape(ny, nx)
            if np.max(np.abs(res)) < tol * max(1.0, np.max(np.abs(rhs))):
                break
    return p, it


def solve_pressure(rhs, g: MacGrid, method: str = "splu", pin: str = "mean", compat_tol: float = 1e-9,
                   return_info: bool = False, n_sweeps: int | None = None):
    """Solve the discrete pressure Poisson equation ∇²_d p = rhs (Eq. (10.124)) for p at the cell centres.

    Book: §10.4, Eq. (10.124).  The pure-Neumann (or periodic) problem is solvable only if Σ rhs·ΔxΔy = 0 (compatibility:
    the net flux out of the whole domain vanishes); this is checked (relative to Σ|rhs|) and a ValueError raised when it
    fails.  ``method``: "splu" (sparse LU of the matrix with one row replaced by p₀ = 0, factorised once per grid and
    cached), "fft" (fully periodic grids), "sor" (red-black SOR, the ch06 relaxation idea — teaching only; ``n_sweeps``
    fixes the number of sweeps).  ``pin``: "mean" (zero-mean p) or "corner".
    Returns p (ny, nx) [and dict(compat, iterations) if return_info].  Label: analytic, converged.
    """
    rhs = np.asarray(rhs, dtype=float)
    compat = float(rhs.sum() * g.dx * g.dy)
    scale = float(np.abs(rhs).sum() * g.dx * g.dy)
    if abs(compat) > compat_tol * max(scale, 1.0):
        raise ValueError(f"incompatible Poisson right-hand side: sum(rhs) dx dy = {compat:.3e} (must be 0 for Neumann/periodic)")
    it = 0
    if method == "fft":
        if not all(g.periodic):
            raise ValueError("fft needs a fully periodic grid")
        p = _fft_poisson(rhs, g)
    elif method == "sor":
        p, it = _sor_neumann(rhs, g, n_sweeps=n_sweeps)
    elif method == "splu":
        b = rhs.ravel().copy()
        b[0] = 0.0
        p = _factor(g.key(), 0).solve(b).reshape(g.ny, g.nx)
    else:
        raise ValueError("method must be splu, fft or sor")
    p = pin_pressure(p, pin)
    if return_info:
        return p, dict(compat=compat, iterations=it)
    return p


def correct(u_s, v_s, p, g: MacGrid, dt: float):
    """Velocity correction with the adjacent-pressure differences.

    Book: §10.4, Eqs. (10.121)–(10.122): u^{n+1}_{i+1/2,j} = u^{n+1/2}_{i+1/2,j} − (Δt/Δx)(p_{i+1,j} − p_{i,j}), same for v.
    Wall faces are not corrected.  Returns (u, v).  Label: analytic.
    """
    gx, gy = gradient(p, g)
    return u_s - dt * gx, v_s - dt * gy  # Eqs. (10.121)–(10.122)


def _set_normal_bc(u, v, g: MacGrid):
    u = np.array(u, dtype=float)
    v = np.array(v, dtype=float)
    if not g.periodic[0]:
        u[:, 0] = 0.0
        u[:, -1] = 0.0
    if not g.periodic[1]:
        v[0, :] = 0.0
        v[-1, :] = 0.0
    return u, v


def project(u_s, v_s, g: MacGrid, dt: float, method: str = "splu", return_parts: bool = False):
    """Projection step: find p^{n+1} and the divergence-free u^{n+1} from the intermediate velocity.

    Book: §10.4, Eqs. (10.117)–(10.118) (continuous), (10.121)–(10.124) (staggered): solve ∇²_d p = ∇_d·u^{n+1/2}/Δt, then
    u^{n+1} = u^{n+1/2} − Δt∇p.  The wall-face normal velocities are reset to their boundary values first (the book: the
    boundary u^{n+1} is known and not expressed through u^{n+1/2} — so the result is independent of the boundary values of
    u^{n+1/2}, as Peyret & Taylor observed).  ``method="auto"`` → FFT on a fully periodic grid, else sparse LU.
    Returns (u, v, p) or, with ``return_parts``, dict(u, v, p, div_before, div_after, rhs, rhs_sum (compatibility), gx, gy
    (∇p on the faces), corr_u, corr_v (the correction −Δt∇p), curl_correction (interior corners)).
    Validation: V4 max|∇_d·u^{n+1}| ≤ 1e-12 (relative); V1 discrete curl of the correction = 0.  Label: conserved, analytic.
    """
    us, vs = _set_normal_bc(u_s, v_s, g)
    div0 = divergence(us, vs, g)
    rhs = div0 / dt  # Eq. (10.124) right-hand side
    if method == "auto":
        method = "fft" if all(g.periodic) else "splu"
    p = solve_pressure(rhs, g, method=method)
    u, v = correct(us, vs, p, g, dt)
    if not return_parts:
        return u, v, p
    cu, cv = u - us, v - vs
    gx, gy = gradient(p, g)
    gtmp = MacGrid(g.nx, g.ny, g.Lx, g.Ly, g.periodic, dict(top=0.0, bottom=0.0, left=0.0, right=0.0), "linear")
    return dict(u=u, v=v, p=p, div_before=div0, div_after=divergence(u, v, g), rhs=rhs, rhs_sum=float(rhs.sum() * g.dx * g.dy),
                gx=gx, gy=gy, corr_u=cu, corr_v=cv, curl_correction=curl(cu, cv, gtmp, interior_only=True))


def projection_stages(us, vs, g: MacGrid, dt: float, method: str = "splu") -> dict:
    """The four stages of one projection, for teaching: divergence of u^{n+1/2}, the Poisson right-hand side and its
    compatibility sum, the pressure and its face gradient, the corrected velocity and its divergence.

    Book: §10.4, (10.117)–(10.124).  (us, vs) is the intermediate field u^{n+1/2} (e.g. from :func:`predictor`).
    Returns dict(div_before, rhs, rhs_sum, p, gx, gy, u, v, div_after, curl_correction, max_div_before, max_div_after).
    Label: conserved.
    """
    parts = project(us, vs, g, dt, method=method, return_parts=True)
    return dict(div_before=parts["div_before"], rhs=parts["rhs"], rhs_sum=parts["rhs_sum"], p=parts["p"], gx=parts["gx"],
                gy=parts["gy"], u=parts["u"], v=parts["v"], div_after=parts["div_after"],
                curl_correction=parts["curl_correction"], max_div_before=float(np.max(np.abs(parts["div_before"]))),
                max_div_after=float(np.max(np.abs(parts["div_after"]))))


def mac_projection_summary(n: int = 8, kind: str = "divergent", dt: float = 1.0) -> dict:
    """One projection of a deterministic intermediate field on an n × n walled unit square (the explainer's parity numbers).

    Book: §10.4, (10.117)–(10.124).  Fields (sampled on the faces, wall-normal faces then set to 0):
    "divergent": u* = sin πx cos πy + ½ sin 2πx, v* = cos πx sin πy;  "shear": u* = sin²πx sin 2πy, v* = 0.3 sin πx sin πy.
    p pinned by zero mean.  Returns dict(div_before_max, div_after_max, rhs_sum, p_max, p_min, curl_correction_max).
    Label: conserved.
    """
    g = MacGrid(n, n)
    c = face_coordinates(g)
    xu, yu, xv, yv = c["xu"], c["yu"], c["xv"], c["yv"]
    if kind == "divergent":
        us = np.sin(np.pi * xu) * np.cos(np.pi * yu) + 0.5 * np.sin(2 * np.pi * xu)
        vs = np.cos(np.pi * xv) * np.sin(np.pi * yv)
    elif kind == "shear":
        us = np.sin(np.pi * xu) ** 2 * np.sin(2 * np.pi * yu)
        vs = 0.3 * np.sin(np.pi * xv) * np.sin(np.pi * yv)
    else:
        raise ValueError("kind must be divergent or shear")
    us, vs = _set_normal_bc(us, vs, g)
    st = projection_stages(us, vs, g, dt)
    return dict(div_before_max=st["max_div_before"], div_after_max=st["max_div_after"], rhs_sum=st["rhs_sum"],
                p_max=float(st["p"].max()), p_min=float(st["p"].min()),
                curl_correction_max=float(np.max(np.abs(st["curl_correction"]))))


# ----------------------------------------------------------------------------------------------------------------------
# time stepping
# ----------------------------------------------------------------------------------------------------------------------
def dt_limit(umax: float, vmax: float, Re: float, dx: float, dy: float | None = None, safety: float = 0.8) -> float:
    """Largest time step allowed by the MAC stability conditions, times a safety factor.

    Book: §10.4, Eq. (10.127) ½(u² + v²)ΔtRe ≤ 1 and Eq. (10.128) 4Δt/(ReΔx²) ≤ 1 (Δx = Δy; for Δx ≠ Δy we use the 2-D
    FTCS diffusion limit 2Δt(1/Δx² + 1/Δy²)/Re ≤ 1, which reduces to (10.128)).  Cited to Peyret & Taylor (1983); checked by
    :func:`fluidpy.core.fd.ftcs2d_max_amplification`.
    Parameters: umax, vmax (non-dimensional speeds), Re, dx, dy, safety.  Returns Δt.  Examples (safety 1):
    (1, 0, 100, 1/32) → 0.02; (1, 0, 100, 1/64) → 0.0061035.  Label: analytic.
    """
    dy = dx if dy is None else dy
    s2 = umax ** 2 + vmax ** 2
    dt_c = 2.0 / (s2 * Re) if s2 > 0 else np.inf  # Eq. (10.127)
    dt_d = Re / (2.0 * (1.0 / dx ** 2 + 1.0 / dy ** 2))  # Eq. (10.128) (= Re Δx²/4 when Δx = Δy)
    return float(safety * min(dt_c, dt_d))


def step(state: dict, g: MacGrid, Re: float, dt: float, lid: float | None = None, body=None, form: str = "advective",
         method: str = "splu") -> dict:
    """One MAC time step: predictor (10.119)–(10.120) → Poisson (10.124) → correction (10.121)–(10.122).

    Book: §10.4 (the algorithm in the paragraph "In summary …", printed p. 445).  First order in time (splitting).
    ``lid``: top-wall speed override (default None → the grid's ``walls["top"]``; pass lid=1.0 or build the grid with
    walls["top"] = 1 for the cavity; ignored for a grid periodic in y).
    Returns the new state dict(u, v, p, t).  Label: converged.
    """
    us, vs = predictor(state["u"], state["v"], g, Re, dt, lid, body, form)
    u, v, p = project(us, vs, g, dt, method=method)
    return dict(u=u, v=v, p=p, t=state["t"] + dt)


def run(state: dict, g: MacGrid, Re: float, dt: float, nsteps: int, callback=None, lid: float | None = None, body=None,
        form: str = "advective", method: str = "splu", tol_steady: float | None = None, save_every: int | None = None,
        check_stability: bool = True, history_every: int = 50) -> dict:
    """March :func:`step` up to ``nsteps`` times (or until steady: max|u^{n+1} − u^n|/Δt < tol_steady).

    Book: §10.4.  ``check_stability``: raise ValueError if dt exceeds the (10.127)–(10.128) limit (safety 1) for the current
    maximum speeds including the lid.  ``lid``: top-wall speed override (default None → the grid's ``walls["top"]``).
    ``callback(k, state)`` after every step.
    Returns dict(state, steps, converged, snapshots (list of (t, u, v) if save_every), history (list of (t, max|Δu|/Δt)
    every ``history_every`` steps), residual).  Label: converged.
    """
    g = _with_lid(g, lid)
    lidv = max(abs(float(x)) for x in g.walls.values()) if g.walls else 0.0
    if check_stability:
        um = max(float(np.max(np.abs(state["u"]))), lidv)
        vm = float(np.max(np.abs(state["v"])))
        lim = dt_limit(max(um, 1e-300), vm, Re, g.dx, g.dy, safety=1.0)
        if dt > lim * (1 + 1e-12):
            raise ValueError(f"dt = {dt:.4g} exceeds the MAC stability limit {lim:.4g} of (10.127)-(10.128)")
    snaps, hist = [], []
    conv, res, k = False, np.inf, 0
    for k in range(1, nsteps + 1):
        new = step(state, g, Re, dt, None, body, form, method)
        if tol_steady or k == nsteps or k % history_every == 0:
            res = max(np.max(np.abs(new["u"] - state["u"])), np.max(np.abs(new["v"] - state["v"]))) / dt
            if k % history_every == 0:
                hist.append((new["t"], float(res)))
        state = new
        if save_every and k % save_every == 0:
            snaps.append((state["t"], state["u"].copy(), state["v"].copy()))
        if callback is not None:
            callback(k, state)
        if not np.isfinite(state["u"]).all() or np.max(np.abs(state["u"])) > 1e6:
            break
        if tol_steady and res < tol_steady:
            conv = True
            break
    return dict(state=state, steps=k, converged=conv, snapshots=snaps, history=hist, residual=float(res))


def _cavity_output(g, s, steps, converged, residual, history, snapshots, dt, Re):
    psi = streamfunction(s["u"], s["v"], g)
    return dict(u=s["u"], v=s["v"], p=s["p"], psi=psi, omega=curl(s["u"], s["v"], g), t=s["t"], steps=steps,
                converged=converged, residual=residual, g=g, grid=g, history=history, snapshots=snapshots, dt=dt, Re=Re,
                centre=primary_vortex_centre(psi, g))


def cavity(Re: float = 100.0, n: int = 32, t_end: float = 20.0, tol_steady: float = 1e-6, dt: float | None = None,
           save_every: int | None = None, cache: bool = True, lid: float = 1.0, form: str = "advective",
           ghost: str = "linear", fast: bool = False) -> dict:
    """Lid-driven square cavity (unit side, lid speed 1 at y = 1) marched from rest to steady state by the MAC scheme.

    Book: §10.5 (the cavity of Fig. 10.6 — the book solves it with MacCormack (``core.maccormack``); here the incompressible
    MAC projection scheme of §10.4 on the same geometry).  Non-dimensional: side L, lid speed U, Re = UL/ν.
    Conventions (design C.10.4, mirrored by the explainer): Δt = 0.8 × :func:`dt_limit` (1, 0, Re, 1/n, safety = 1) unless
    given; steady when max|Δu|/Δt < ``tol_steady`` (0 → run to t_end); the lid enters through the ghost value
    u_ghost = 2U − u_top (``ghost="linear"``; "quadratic" is the second-order alternative); centred convective differences
    with the face averages of :func:`predictor`.
    Parameters: Re; n (cells per side; ``fast`` halves it); t_end; tol_steady; dt; save_every (snapshots, disables the
    cache); cache (npz in outputs/ch10 keyed by the parameters); lid; form; ghost.
    Returns dict(u, v, p, psi (corner streamfunction), omega (corner vorticity), t, steps, converged, residual, g (grid),
    history (list of (t, max|Δu|/Δt)), snapshots, dt, Re, centre (primary vortex dict), cached (path, if loaded)).
    Cost: 32² ≈ 1 s, 64² ≈ 5 s, 128² ≈ 1–2 min.  Validation: V5 Ghia et al. (1982) centreline; V4 ∇·u ≈ 0 each step;
    V3 grid convergence.  Label: converged, benchmark.
    """
    if fast:
        n = max(8, n // 2)
    g = MacGrid(n, n, 1.0, 1.0, (False, False), dict(top=lid, bottom=0.0, left=0.0, right=0.0), ghost)
    if dt is None:
        dt = 0.8 * dt_limit(abs(lid), 0.0, Re, g.dx, g.dy, safety=1.0)
    path = None
    if cache and not save_every:
        from .project import repo_root

        tag = f"Re{Re}_n{n}_t{t_end}_tol{tol_steady}_dt{dt}_lid{lid}_{form}_{ghost}"
        h = hashlib.md5(tag.encode()).hexdigest()[:8]
        path = repo_root() / "outputs" / "ch10" / f"mac_cavity_Re{int(Re)}_n{n}_{h}.npz"
        if path.is_file():
            d = np.load(path, allow_pickle=True)
            s = dict(u=d["u"], v=d["v"], p=d["p"], t=float(d["t"]))
            out = _cavity_output(g, s, int(d["steps"]), bool(d["converged"]), float(d["residual"]),
                                 [tuple(r) for r in d["history"]], [], dt, Re)
            out["cached"] = str(path)
            return out
    r = run(new_state(g), g, Re, dt, int(np.ceil(t_end / dt)), lid=None, form=form, tol_steady=tol_steady,
            save_every=save_every, check_stability=False)
    out = _cavity_output(g, r["state"], r["steps"], r["converged"], r["residual"], r["history"], r["snapshots"], dt, Re)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        s = r["state"]
        np.savez(path, u=s["u"], v=s["v"], p=s["p"], t=s["t"], steps=r["steps"], converged=r["converged"],
                 residual=r["residual"], history=np.array(r["history"], dtype=float).reshape(-1, 2))
    return out


def _tg_run(n, Re, t_end, dt, safety, form):
    L = 2 * np.pi
    g = MacGrid(n, n, L, L, (True, True))
    c = face_coordinates(g)
    st = dict(u=np.sin(c["xu"]) * np.cos(c["yu"]), v=-np.cos(c["xv"]) * np.sin(c["yv"]), p=np.zeros((n, n)), t=0.0)
    if dt is None:
        dt = dt_limit(1.0, 1.0, Re, g.dx, g.dy, safety)
    nst = max(1, int(round(t_end / dt)))
    dt = t_end / nst
    out = run(st, g, Re, dt, nst, lid=None, form=form, method="fft", check_stability=False)
    s = out["state"]
    f = np.exp(-2 * t_end / Re)
    ue = np.sin(c["xu"]) * np.cos(c["yu"]) * f
    ve = -np.cos(c["xv"]) * np.sin(c["yv"]) * f
    pe = 0.25 * (np.cos(2 * c["xp"]) + np.cos(2 * c["yp"])) * f ** 2
    pe = pe - pe.mean()
    return dict(err_u=float(np.max(np.abs(s["u"] - ue))), err_v=float(np.max(np.abs(s["v"] - ve))),
                err_p=float(np.max(np.abs(s["p"] - pe))), grid=g, state=s, dt=dt, steps=out["steps"])


def taylor_green(n: int = 32, Re: float = 100.0, t_end: float = 1.0, periodic: bool = True, dt: float | None = None,
                 safety: float = 0.5, form: str = "advective") -> dict:
    """Decaying Taylor–Green vortex on the periodic box [0, 2π]² — the exact-solution test of the MAC solver.

    Book: §10.4 (code verification; ch04's ``exact_solution("taylor_green")`` in non-dimensional form):
    u = sin x cos y e^{−2t/Re}, v = −cos x sin y e^{−2t/Re}, p = ¼(cos 2x + cos 2y)e^{−4t/Re}.  FFT pressure solve.
    Also runs n/2 with the same Δt rule to report the observed order (the MAC splitting is first order in time).
    Returns dict(err_u, err_v, err_p (max; p mean-removed), order, dt, steps, grid, state).  Validation: V1/V3.
    Label: converged.
    """
    if not periodic:
        raise ValueError("the Taylor–Green test is on the periodic box")
    r = _tg_run(n, Re, t_end, dt, safety, form)
    r2 = _tg_run(max(4, n // 2), Re, t_end, dt, safety, form)
    r["order"] = float(np.log(r2["err_u"] / r["err_u"]) / np.log(2.0)) if r["err_u"] > 0 else float("nan")
    return r


def channel_poiseuille(nx: int = 4, ny: int = 16, Re: float = 10.0, dpdx: float = -1.0, t_end: float | None = None,
                       dt: float | None = None, ghost: str = "quadratic") -> dict:
    """Pressure-gradient-driven channel between walls y = 0, 1 (periodic in x), marched to the steady plane Poiseuille flow.

    Book: §10.4 (code verification; ch04/ch08 exact field): the constant body force G = −dp/dx stands for the imposed mean
    pressure gradient; the steady solution is u = (Re G/2) y(1 − y) (non-dimensional).  With quadratic wall ghosts
    (default) the discrete steady state equals the parabola to round-off; ``ghost="linear"`` shows an O(Δy²) error.
    Parameters: nx, ny, Re, dpdx (< 0 drives +x flow), t_end (default 3Re), dt, ghost.
    Returns dict(y, u (profile), exact, max_err, steps, converged).  Label: analytic.
    """
    G = -dpdx
    g = MacGrid(nx, ny, 1.0, 1.0, (True, False), dict(top=0.0, bottom=0.0, left=0.0, right=0.0), ghost)
    umax = abs(Re * G / 8.0)
    if dt is None:
        dt = dt_limit(max(umax, 1e-12), 0.0, Re, g.dx, g.dy, 0.4)
    t_end = 3.0 * Re if t_end is None else t_end
    out = run(new_state(g), g, Re, dt, int(np.ceil(t_end / dt)), lid=None, body=(G, 0.0), tol_steady=1e-13,
              check_stability=False)
    y = g.yc
    ex = 0.5 * Re * G * y * (1 - y)
    u = out["state"]["u"][:, 0]
    return dict(y=y, u=u, exact=ex, max_err=float(np.max(np.abs(u - ex))), steps=out["steps"], converged=out["converged"])


# ----------------------------------------------------------------------------------------------------------------------
# post-processing
# ----------------------------------------------------------------------------------------------------------------------
def streamfunction(u, v, g: MacGrid) -> np.ndarray:
    """Discrete streamfunction ψ at the cell corners, u = ∂ψ/∂y, v = −∂ψ/∂x, ψ = 0 on the bottom wall.

    Book: §10.5 (streamline pictures, Fig. 10.7; Ch. 4 ψ).  On the staggered grid ψ follows exactly from the face fluxes:
    ψ(x_i, y_j) = Σ_{k<j} u[k, i]Δy (a discretely divergence-free field makes this path-independent and ψ = 0 on every
    wall of a closed cavity).  Shape (ny + 1, nx + 1) with walls.  Label: analytic, conserved.
    """
    ny, nxu = u.shape
    psi = np.zeros((ny + 1, nxu))
    psi[1:, :] = np.cumsum(u, axis=0) * g.dy
    return psi


def primary_vortex_centre(psi, g: MacGrid) -> dict:
    """Location and value of the minimum of ψ (the primary eddy of the cavity, clockwise for a lid moving in +x).

    Book: §10.5 (primary-eddy centres compared with Hou et al. 1995).  Sub-cell location by a quadratic fit through the
    corners around the discrete minimum (three-point parabola in x and in y through the discrete minimum; the fitted
    ψ_min is therefore ≤ the discrete minimum).  Returns dict(x, y, psi_min).
    Validation: V1 exact on a separable paraboloid; the grid study of ψ_min is reported, not asserted here.
    Label: qualitative (demonstration).
    """
    j, i = np.unravel_index(np.argmin(psi), psi.shape)
    x, y, pm = i * g.dx, j * g.dy, float(psi[j, i])
    if 0 < i < psi.shape[1] - 1 and 0 < j < psi.shape[0] - 1:
        fx = psi[j, i - 1: i + 2]
        fy = psi[j - 1: j + 2, i]
        denx = fx[0] - 2 * fx[1] + fx[2]
        deny = fy[0] - 2 * fy[1] + fy[2]
        sx = 0.5 * (fx[0] - fx[2]) / denx if denx != 0 else 0.0
        sy = 0.5 * (fy[0] - fy[2]) / deny if deny != 0 else 0.0
        x += sx * g.dx
        y += sy * g.dy
        # vertex of f(s) = f1 + ½(f2 − f0)s + ½(f0 − 2f1 + f2)s² at s* = (f0 − f2)/(2 den) is f1 + ¼(f2 − f0)s*
        # (separable fit: the x and y corrections add)
        pm = float(psi[j, i] + 0.25 * ((fx[2] - fx[0]) * sx + (fy[2] - fy[0]) * sy))
    return dict(x=float(x), y=float(y), psi_min=pm)


def cavity_centreline(state, g: MacGrid | None = None, lid: float | None = None) -> dict:
    """u(y) on the vertical centreline x = ½ (Fig. 10.8), with the wall values u(0) = 0 and u(1) = lid.

    Book: §10.5, Fig. 10.8.  ``state``: the dict of :func:`cavity` (with "g"), or the u array plus ``g``.  For even nx the
    line x = ½ is a column of u faces; for odd nx the two neighbouring columns are averaged.
    Returns dict(y (ascending, 0 … 1), u).  Label: analytic.
    """
    if isinstance(state, dict):
        u = state["u"]
        g = g or state.get("g") or state.get("grid")
    else:
        u = state
    if g.nx % 2:
        i = g.nx // 2
        uc = 0.5 * (u[:, i] + u[:, i + 1])
    else:
        uc = u[:, g.nx // 2]
    top = g.walls.get("top", 0.0) if lid is None else lid
    return dict(y=np.concatenate([[0.0], g.yc, [g.Ly]]), u=np.concatenate([[g.walls.get("bottom", 0.0)], uc, [top]]))


def cavity_centreline_v(state, g: MacGrid | None = None) -> dict:
    """v(x) on the horizontal centreline y = ½, with v = 0 on the side walls.  Book: §10.5.  Returns dict(x, v).
    Label: analytic."""
    if isinstance(state, dict):
        v = state["v"]
        g = g or state.get("g") or state.get("grid")
    else:
        v = state
    if g.ny % 2:
        j = g.ny // 2
        vc = 0.5 * (v[j, :] + v[j + 1, :])
    else:
        vc = v[g.ny // 2, :]
    return dict(x=np.concatenate([[0.0], g.xc, [g.Lx]]), v=np.concatenate([[0.0], vc, [0.0]]))
