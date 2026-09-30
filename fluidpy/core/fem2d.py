"""Mixed finite elements for the non-dimensional incompressible Navier–Stokes equations on triangles: quadratic (P2)
velocity and linear (P1) pressure (Taylor–Hood), isoparametric curved elements, 7-point quadrature, Newton linearisation,
essential boundary conditions by row replacement or penalty, meshes for the unit square and a cylinder in a channel, and the
discrete inf–sup (LBB) constant.

Book: Kundu, Cohen & Dowling 5e, Ch. 10 §10.4 Eqs. (10.134)–(10.137) (weak NS, strain rate, saddle-point system) and §10.5
Eqs. (10.156)–(10.198) (the cylinder in a channel: Cartesian weak form, time derivative (10.163), Newton (10.164)–(10.167),
expansions (10.168)–(10.169), algebraic equations (10.170)–(10.172), block form (10.173)–(10.183), isoparametric map
(10.184), shape functions (10.185)–(10.187), element matrices (10.188)–(10.197), quadrature (10.198)).  Pages
chapters/pages/ch10/p473–p475, p486–p493.

Conventions.  Non-dimensional (length d or the square side, velocity U, time d/U, pressure ρU², Re = Ud/ν).  Node array
``mesh.nodes`` (N, 2) lists the triangle vertices first (indices 0 … V−1, which are also the P1 pressure nodes) then the edge
mid-nodes (V … V+E−1).  Element connectivity ``mesh.elems`` (T, 6) follows Fig. 10.17: local nodes 1, 2, 3 = vertices
(ξ, η) = (0, 0), (1, 0), (0, 1); 4 = mid(1, 2), 5 = mid(2, 3), 6 = mid(3, 1); ζ = 1 − ξ − η.  Unknown vector [u (N), v (N),
p (V)].  Printed slips: R6 (10.172) second sum printed with u_{A′} — ``assemble_newton_system(printed=True)`` builds that
wrong continuity row (fails the Poiseuille/Stokes tests); R7 (10.186) "v′" for the pressure expansion (coded p′); R10
(10.166) spurious star on β ∂v/∂t(t_n) (known data at t_n).
Reused by Ch. 14 (low-Re bodies), Ch. 16 (Stokes flow in irregular geometry), Ch. 11 (base flows for wake stability).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "p2_shape", "p2_shape_grad", "p1_shape", "p1_shape_grad", "iso_map", "jacobian", "jacobian_matrix", "tri_quad_7pt",
    "integrate_element", "assemble_saddle", "cylinder_steady",
    "Mesh", "mesh_from_triangles", "structured_square_mesh", "structured_tri_mesh_p2", "cylinder_channel_mesh",
    "p2_node_counts", "euler_check", "time_derivative", "element_newton_matrices", "assemble_newton_system",
    "apply_dirichlet", "newton_solve", "stokes_solve", "stokes_p2p1", "stokes_cavity", "kovasznay_exact", "kovasznay_test",
    "poiseuille_test", "infsup_constant", "infsup_table", "channel_bcs", "cylinder_forces", "march_unsteady",
    "nodal_vorticity", "BoundaryConditions",
]


# ======================================================================================================================
# Parent element: shape functions (10.185)–(10.187), map (10.184), quadrature (10.198)
# ======================================================================================================================
def p2_shape(xi, eta) -> np.ndarray:
    """Quadratic triangle shape functions φ₁ … φ₆ at (ξ, η).

    Book: §10.5, Eq. (10.185): φ₁ = ζ(2ζ − 1), φ₂ = ξ(2ξ − 1), φ₃ = η(2η − 1), φ₄ = 4ξζ, φ₅ = 4ξη, φ₆ = 4ηζ, ζ = 1 − ξ − η
    (node numbering of Fig. 10.17).  Returns array (6, …).  Validation: V1 φ_a(node_b) = δ_ab, Σφ_a = 1.  Label: analytic.
    """
    xi = np.asarray(xi, dtype=float)
    eta = np.asarray(eta, dtype=float)
    ze = 1.0 - xi - eta
    return np.array([ze * (2 * ze - 1), xi * (2 * xi - 1), eta * (2 * eta - 1), 4 * xi * ze, 4 * xi * eta, 4 * eta * ze])  # Eq. (10.185)


def p2_shape_grad(xi, eta) -> np.ndarray:
    """∂φ_a/∂ξ and ∂φ_a/∂η of (10.185): array (6, 2, …).  Label: analytic."""
    xi = np.asarray(xi, dtype=float)
    eta = np.asarray(eta, dtype=float)
    ze = 1.0 - xi - eta
    one = np.ones_like(xi)
    dxi = np.array([-(4 * ze - 1), 4 * xi - 1, 0 * one, 4 * (ze - xi), 4 * eta, -4 * eta])
    deta = np.array([-(4 * ze - 1), 0 * one, 4 * eta - 1, -4 * xi, 4 * xi, 4 * (ze - eta)])
    return np.stack([dxi, deta], axis=1)


def p1_shape(xi, eta) -> np.ndarray:
    """Linear pressure shape functions ψ₁ = ζ, ψ₂ = ξ, ψ₃ = η.  Book: §10.5, Eq. (10.187).  Returns (3, …).  Label: analytic."""
    xi = np.asarray(xi, dtype=float)
    eta = np.asarray(eta, dtype=float)
    return np.array([1.0 - xi - eta, xi, eta])  # Eq. (10.187)


def p1_shape_grad(xi, eta) -> np.ndarray:
    """(3, 2, …) derivatives of (10.187).  Label: analytic."""
    one = np.ones_like(np.asarray(xi, dtype=float))
    return np.stack([np.array([-one, one, 0 * one]), np.array([-one, 0 * one, one])], axis=1)


def iso_map(xe, ye, xi, eta):
    """Isoparametric map of the parent triangle onto a (possibly curved) six-node element.

    Book: §10.5, Eq. (10.184): x(ξ, η) = Σ_{a=1}^6 x_a^e φ_a(ξ, η), y likewise (Fig. 10.17).  xe, ye: the six node
    coordinates.  Returns (x, y).  Label: analytic.
    """
    N = p2_shape(xi, eta)
    return np.tensordot(np.asarray(xe, float), N, axes=(0, 0)), np.tensordot(np.asarray(ye, float), N, axes=(0, 0))  # Eq. (10.184)


def jacobian_matrix(xe, ye, xi, eta):
    """Jacobian matrix of (10.184) and its determinant: (J = [[x_ξ, x_η], [y_ξ, y_η]], det J).  Label: analytic."""
    dN = p2_shape_grad(xi, eta)
    xe = np.asarray(xe, float)
    ye = np.asarray(ye, float)
    x_xi, x_eta = np.tensordot(xe, dN[:, 0], axes=(0, 0)), np.tensordot(xe, dN[:, 1], axes=(0, 0))
    y_xi, y_eta = np.tensordot(ye, dN[:, 0], axes=(0, 0)), np.tensordot(ye, dN[:, 1], axes=(0, 0))
    return np.array([[x_xi, x_eta], [y_xi, y_eta]]), x_xi * y_eta - x_eta * y_xi


def jacobian(xe, ye, xi, eta):
    """Jacobian determinant J = x_ξ y_η − x_η y_ξ of the isoparametric map (10.184) at (ξ, η).

    Book: §10.5, below (10.198).  For a straight-sided element J = 2 × area (constant).  Returns float (or array).
    Label: analytic.
    """
    det = jacobian_matrix(xe, ye, xi, eta)[1]
    return float(det) if np.ndim(det) == 0 else det


def tri_quad_7pt():
    """Seven-point quadrature on the parent triangle, degree of precision 5 (Radon's rule).

    Book: §10.5, Eq. (10.198): ∫_{Ωᵉ} f dΩ = ½ Σ_l f(ξ_l, η_l) J(ξ_l, η_l) W_l with Σ W_l = 1 (the ½ is the parent area).
    Points: the centroid (W = 9/40) and two orbits of three points with barycentric coordinates built from
    a = (6 − √15)/21, b = (6 + √15)/21 and weights (155 ∓ √15)/1200.  Returns (points (7, 2) as (ξ, η), W (7,)).
    Validation: V1 integrates every ξ^p η^q with p + q ≤ 5 exactly (1e-15), not degree 6.  Label: analytic.
    """
    s15 = np.sqrt(15.0)
    a, b = (6.0 - s15) / 21.0, (6.0 + s15) / 21.0
    wa, wb = (155.0 - s15) / 1200.0, (155.0 + s15) / 1200.0
    pts = [(1 / 3, 1 / 3), (a, a), (1 - 2 * a, a), (a, 1 - 2 * a), (b, b), (1 - 2 * b, b), (b, 1 - 2 * b)]
    W = [9 / 40, wa, wa, wa, wb, wb, wb]
    return np.array(pts), np.array(W)


def integrate_element(f, xe, ye) -> float:
    """∫_{Ωᵉ} f(x, y) dΩ over one six-node element by (10.198).  f: callable of (x, y) arrays.  Label: analytic."""
    P, W = tri_quad_7pt()
    x, y = iso_map(xe, ye, P[:, 0], P[:, 1])
    det = jacobian(xe, ye, P[:, 0], P[:, 1])
    return float(0.5 * np.sum(f(x, y) * det * W))  # Eq. (10.198)


# ======================================================================================================================
# Meshes
# ======================================================================================================================
@dataclass
class Mesh:
    """Triangular mesh with P2 (or P1) velocity nodes and P1 pressure nodes (the vertices).

    nodes (N, 2); elems (T, 6) for P2 (T, 3) for P1; V vertices (= pressure nodes 0 … V−1); E edges; pair "P2P1"/"P1P1";
    boundary: dict name → sorted node indices; curved: whether some mid-nodes were moved onto a curved boundary.
    """
    nodes: np.ndarray
    elems: np.ndarray
    V: int
    E: int
    pair: str = "P2P1"
    boundary: dict = field(default_factory=dict)
    curved: bool = False
    info: dict = field(default_factory=dict)

    @property
    def T(self) -> int:
        return int(self.elems.shape[0])

    @property
    def N(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def tris(self) -> np.ndarray:
        return self.elems[:, :3]

    def __getitem__(self, key):
        """Dict-style access (design contract): vertices, tris, edges, vel_nodes, p_nodes, boundary_vel, counts."""
        if key == "vertices":
            return self.nodes[:self.V]
        if key == "tris":
            return self.tris
        if key == "edges":
            return self.info.get("edges")
        if key == "vel_nodes":
            return self.nodes
        if key == "p_nodes":
            return self.nodes[:self.V]
        if key == "boundary_vel":
            return self.boundary.get("all", np.unique(np.concatenate(list(self.boundary.values()))))
        if key == "counts":
            c = p2_node_counts(self.V, self.E, self.T)
            return dict(V=self.V, E=self.E, T=self.T, N_u=self.N if self.pair == "P2P1" else self.V, N_p=self.V,
                        euler=c["euler"])
        raise KeyError(key)


def mesh_from_triangles(verts, tris, pair: str = "P2P1", circle=None, edge_order=None) -> Mesh:
    """Build a :class:`Mesh` from a P1 triangulation, adding one mid-node per edge for P2 (Fig. 10.5a).

    Book: §10.4 Fig. 10.5a (six velocity nodes, three pressure vertices), §10.5 (curved sides: ``circle`` = (xc, yc, r) moves
    the mid-node of every edge with both ends on the circle onto it — the isoparametric element of Fig. 10.17).
    Triangles are re-oriented counter-clockwise.  Label: analytic.
    """
    verts = np.asarray(verts, dtype=float)
    tris = np.asarray(tris, dtype=int).copy()
    a, b, c = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
    area2 = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    flip = area2 < 0
    tris[flip] = tris[flip][:, [0, 2, 1]]
    V = verts.shape[0]
    edges = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    es = np.sort(edges, axis=1)
    uniq, inv = np.unique(es, axis=0, return_inverse=True)
    inv = inv.ravel()
    E = uniq.shape[0]
    if edge_order is not None:  # renumber the edges: key(uniq) sorted
        order = np.lexsort(edge_order(verts, uniq)[::-1])
        rank = np.empty(E, dtype=int)
        rank[order] = np.arange(E)
        uniq = uniq[order]
        inv = rank[inv]
    if pair == "P1P1":
        return Mesh(verts.copy(), tris, V, E, "P1P1", info=dict(edges=uniq))
    T = tris.shape[0]
    mid = 0.5 * (verts[uniq[:, 0]] + verts[uniq[:, 1]])
    curved = False
    if circle is not None:
        xc, yc, r = circle
        on = lambda p: np.abs(np.hypot(p[:, 0] - xc, p[:, 1] - yc) - r) < 1e-9 * max(1.0, r)  # noqa: E731
        both = on(verts[uniq[:, 0]]) & on(verts[uniq[:, 1]])
        d = mid[both] - np.array([xc, yc])
        mid[both] = np.array([xc, yc]) + r * d / np.linalg.norm(d, axis=1)[:, None]
        curved = bool(both.any())
    nodes = np.vstack([verts, mid])
    e01, e12, e20 = inv[:T], inv[T:2 * T], inv[2 * T:]
    elems = np.column_stack([tris, V + e01, V + e12, V + e20])
    return Mesh(nodes, elems, V, E, "P2P1", curved=curved, info=dict(edges=uniq))


def structured_square_mesh(n: int, pair: str = "P2P1", x0: float = 0.0, x1: float = 1.0, y0: float = 0.0,
                           y1: float = 1.0, diagonal: str = "right") -> Mesh:
    """n × n squares on [x0, x1] × [y0, y1], each cut into two triangles; boundary sets bottom/top/left/right/all.

    Book: §10.4 (Fig. 10.5 elements on a simple domain — our test mesh for LBB and code verification).
    ``diagonal``: "right" (all diagonals the same way) or "cross" (alternating, a union-jack-like pattern).
    Label: analytic.
    """
    xs = np.linspace(x0, x1, n + 1)
    ys = np.linspace(y0, y1, n + 1)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    verts = np.column_stack([X.ravel(), Y.ravel()])
    idx = lambda i, j: j * (n + 1) + i  # noqa: E731
    tris = []
    for j in range(n):
        for i in range(n):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            if diagonal == "cross" and (i + j) % 2:
                tris += [[a, b, d], [b, c, d]]
            else:
                tris += [[a, b, c], [a, c, d]]
    def order(vv, ed):  # (type, row, column) — horizontal, vertical, diagonal; then y, then x of the midpoint
        a, b = vv[ed[:, 0]], vv[ed[:, 1]]
        typ = np.where(np.abs(a[:, 1] - b[:, 1]) < 1e-12, 0, np.where(np.abs(a[:, 0] - b[:, 0]) < 1e-12, 1, 2))
        mid = 0.5 * (a + b)
        return [typ, np.round(mid[:, 1], 12), np.round(mid[:, 0], 12)]

    m = mesh_from_triangles(verts, np.array(tris), pair, edge_order=order)
    x, y = m.nodes[:, 0], m.nodes[:, 1]
    tol = 1e-12 * max(1.0, abs(x1 - x0))
    m.boundary = dict(bottom=np.where(np.abs(y - y0) < tol)[0], top=np.where(np.abs(y - y1) < tol)[0],
                      left=np.where(np.abs(x - x0) < tol)[0], right=np.where(np.abs(x - x1) < tol)[0])
    m.boundary["all"] = np.unique(np.concatenate(list(m.boundary.values())))
    m.info.update(kind="square", n=n, box=(x0, x1, y0, y1))
    return m


structured_tri_mesh_p2 = structured_square_mesh


def cylinder_channel_mesh(d: float = 1.0, W: float = 5.0, xmin: float = -8.0, xmax: float = 16.0, n_theta: int = 48,
                          n_r: int = 12, h_far: float = 0.4, r_ring: float | None = None, symmetric: bool = True,
                          pair: str = "P2P1", h_wake: float | None = None, wake_length: float = 8.0) -> Mesh:
    """Triangular mesh of a channel of width W with a cylinder of diameter d at the origin, refined near the cylinder.

    Book: §10.5, Fig. 10.15 (inflow Γ₁ at x_min, outflow Γ₂ at x_max, sliding walls Γ₃ (top) and Γ₄ (bottom) at y = ±W/2,
    cylinder Γ₅) and Fig. 10.16 (a mesh finer near the cylinder).  Our geometry (the book's numbers are private): W = 5d,
    x from −8d to 16d.  Our mesh: polar rings around the cylinder with geometric radial growth (square cells at the wall;
    at most n_r rings, up to r_ring), a uniform background of spacing h_far elsewhere (optionally h_wake in a box
    0 < x < wake_length, |y| < W/4 behind the body), Delaunay triangulation (scipy.spatial) of the upper half, mirrored to
    make the mesh up–down symmetric (so a steady solution can be tested for symmetry; the book's mesh is not symmetric,
    p. 466).  Curved P2 edges on Γ₅ (mid-nodes on the circle, (10.184)).
    Parameters: d, W, xmin, xmax (in units of d), n_theta (points round the cylinder, even), n_r (max rings), h_far,
    r_ring (default min(0.8·W/2, radius where the ring spacing reaches h_far)), symmetric, pair, h_wake, wake_length.
    Returns :class:`Mesh` with boundary sets inflow, outflow, top, bottom, cylinder.  Label: analytic.
    """
    from scipy.spatial import Delaunay

    r0 = 0.5 * d
    dth = 2 * np.pi / n_theta
    rmax = min(0.8 * 0.5 * W, h_far / dth) if r_ring is None else r_ring
    radii = [r0]
    while len(radii) < n_r:
        r_next = radii[-1] * (1.0 + dth)
        if r_next > rmax:
            break
        radii.append(r_next)
    th = np.arange(n_theta // 2 + 1) * dth  # upper half 0 … π
    pts = [np.column_stack([r * np.cos(th), r * np.sin(th)]) for r in radii]
    ring = np.vstack(pts)
    rin = radii[-1] + 0.6 * h_far
    nx = int(np.ceil((xmax - xmin) / h_far))
    ny = int(np.ceil(0.5 * W / h_far))
    xs = np.linspace(xmin, xmax, nx + 1)
    ys = np.linspace(0.0, 0.5 * W, ny + 1)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    bg = np.column_stack([X.ravel(), Y.ravel()])
    # stagger alternate rows slightly (better-shaped Delaunay triangles in the background)
    keep = np.hypot(bg[:, 0], bg[:, 1]) > rin
    on_edge = (np.abs(bg[:, 1]) < 1e-12) | (np.abs(bg[:, 1] - 0.5 * W) < 1e-12) | (np.abs(bg[:, 0] - xmin) < 1e-12) | \
              (np.abs(bg[:, 0] - xmax) < 1e-12)
    bg = bg[keep | (on_edge & (np.hypot(bg[:, 0], bg[:, 1]) > radii[-1] + 0.25 * h_far))]
    if h_wake is not None:
        wx = np.arange(0.0, wake_length + 1e-12, h_wake)
        wy = np.arange(0.5 * h_wake, 0.25 * W, h_wake)
        WX, WY = np.meshgrid(wx, wy, indexing="xy")
        wk = np.column_stack([WX.ravel(), WY.ravel()])
        wk = wk[np.hypot(wk[:, 0], wk[:, 1]) > rin]
        # drop background points inside the wake box (keep its edges clean)
        inbox = (bg[:, 0] > -0.5 * h_far) & (bg[:, 0] < wake_length + 0.5 * h_far) & (bg[:, 1] < 0.25 * W + 0.5 * h_far) & \
                (bg[:, 1] > 1e-12)
        bg = np.vstack([bg[~inbox], wk])
    P = np.vstack([ring, bg])
    P[np.abs(P[:, 1]) < 1e-13, 1] = 0.0
    P = np.unique(np.round(P, 12), axis=0)
    tri = Delaunay(P)
    T = tri.simplices
    cen = P[T].mean(axis=1)
    T = T[np.hypot(cen[:, 0], cen[:, 1]) > r0 * (1 + 1e-9)]
    # remove slivers on the straight boundary (all three points on one boundary line)
    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))
    T = T[area > 1e-10]
    if symmetric:
        up = P
        low = P.copy()
        low[:, 1] *= -1
        allp = np.vstack([up, low])
        Tl = T + P.shape[0]
        allT = np.vstack([T, Tl])
        allp_r = np.round(allp, 12)
        uniq, inv = np.unique(allp_r, axis=0, return_inverse=True)
        verts = uniq
        tris = inv.ravel()[allT]
    else:
        verts, tris = P, T
    m = mesh_from_triangles(verts, tris, pair, circle=(0.0, 0.0, r0))
    x, y = m.nodes[:, 0], m.nodes[:, 1]
    tol = 1e-9
    m.boundary = dict(inflow=np.where(np.abs(x - xmin) < tol)[0], outflow=np.where(np.abs(x - xmax) < tol)[0],
                      top=np.where(np.abs(y - 0.5 * W) < tol)[0], bottom=np.where(np.abs(y + 0.5 * W) < tol)[0],
                      cylinder=np.where(np.abs(np.hypot(x, y) - r0) < 1e-9)[0])
    m.info.update(kind="cylinder_channel", d=d, W=W, xmin=xmin, xmax=xmax, n_theta=n_theta, h_far=h_far,
                  n_rings=len(radii), r_ring=radii[-1], h_wake=h_wake)
    return m


def p2_node_counts(V: int, E: int, T: int) -> dict:
    """Unknown counts of a P2–P1 mesh: velocity nodes N_u = V + E, pressure nodes N_p = V, and V − E + T.

    Book: §10.5 (the node counts quoted with Fig. 10.16 and for the finer Re = 100 mesh — values private).  For a planar
    triangulation of a domain with h holes, V − E + T = 1 − h (Euler), so 0 for the channel with one cylinder.
    Returns dict(N_u, N_p, euler).  Label: analytic.
    """
    return dict(N_u=int(V + E), N_p=int(V), euler=int(V - E + T))


def euler_check(mesh: Mesh) -> dict:
    """:func:`p2_node_counts` for a mesh (0 for one hole, 1 for a simply connected domain).  Label: analytic."""
    return p2_node_counts(mesh.V, mesh.E, mesh.T)


# ======================================================================================================================
# Geometry cache: quadrature, Jacobians, physical gradients
# ======================================================================================================================
def _geometry(mesh: Mesh):
    key = "_geo"
    if key in mesh.info:
        return mesh.info[key]
    P, W = tri_quad_7pt()
    xi, eta = P[:, 0], P[:, 1]
    p2 = mesh.pair == "P2P1"
    phi = p2_shape(xi, eta) if p2 else p1_shape(xi, eta)  # (nv, 7)
    dphi = p2_shape_grad(xi, eta) if p2 else p1_shape_grad(xi, eta)  # (nv, 2, 7)
    psi = p1_shape(xi, eta)  # (3, 7)
    X = mesh.nodes[mesh.elems]  # (T, nv, 2)
    x_xi = np.einsum("ta,al->tl", X[:, :, 0], dphi[:, 0])
    x_eta = np.einsum("ta,al->tl", X[:, :, 0], dphi[:, 1])
    y_xi = np.einsum("ta,al->tl", X[:, :, 1], dphi[:, 0])
    y_eta = np.einsum("ta,al->tl", X[:, :, 1], dphi[:, 1])
    det = x_xi * y_eta - x_eta * y_xi  # J of (10.198)
    if np.any(det <= 0):
        raise ValueError("inverted or degenerate element (det J <= 0)")
    xi_x, xi_y, eta_x, eta_y = y_eta / det, -x_eta / det, -y_xi / det, x_xi / det
    gx = dphi[None, :, 0, :] * xi_x[:, None, :] + dphi[None, :, 1, :] * eta_x[:, None, :]  # (T, nv, 7)
    gy = dphi[None, :, 0, :] * xi_y[:, None, :] + dphi[None, :, 1, :] * eta_y[:, None, :]
    w = 0.5 * W[None, :] * det  # (T, 7): ½ J W_l of (10.198)
    xq = np.einsum("ta,al->tl", X[:, :, 0], phi)
    yq = np.einsum("ta,al->tl", X[:, :, 1], phi)
    geo = dict(phi=phi, psi=psi, gx=gx, gy=gy, w=w, xq=xq, yq=yq, pel=mesh.elems[:, :3])
    mesh.info[key] = geo
    return geo


def _scatter(rows_el, cols_el, vals, shape):
    import scipy.sparse as sps

    R = np.broadcast_to(rows_el[:, :, None], vals.shape)
    C = np.broadcast_to(cols_el[:, None, :], vals.shape)
    return sps.coo_matrix((vals.ravel(), (R.ravel(), C.ravel())), shape=shape).tocsr()


def _vec_scatter(rows_el, vals, n):
    out = np.zeros(n)
    np.add.at(out, rows_el.ravel(), vals.ravel())
    return out


# ======================================================================================================================
# Time derivative (10.163) and the Newton system (10.165)–(10.183), (10.188)–(10.197)
# ======================================================================================================================
def time_derivative(u_new, u_old, dudt_old, dt: float, alpha_t: float = 1.0, beta_t: float = 0.0):
    """∂u/∂t(t_{n+1}) ≈ α(u(t_{n+1}) − u(t_n))/Δt − β ∂u/∂t(t_n).

    Book: §10.5, Eq. (10.163): α = 1, β = 0 is first order (backward Euler); α = 2, β = 1 second order (the trapezoidal rule,
    "a variation of Crank–Nicolson").  Returns the rate (unit of u per unit time).  Label: converged.
    """
    return alpha_t * (np.asarray(u_new) - np.asarray(u_old)) / dt - beta_t * np.asarray(dudt_old)  # Eq. (10.163)


def element_newton_matrices(mesh: Mesh, U, Re: float, dt: float | None = None, alpha_t: float = 1.0, beta_t: float = 0.0,
                            U_old=None, dUdt_old=None, convection: bool = True):
    """Element matrices and vectors of the Newton step for every element at once.

    Book: §10.5, Eqs. (10.188)–(10.197) (the element versions of (10.175)–(10.183)):
    A^{euu} = ∫[(α/Δt φ_{a′} + u*φ_{a′,x} + v*φ_{a′,y} + u*_x φ_{a′})φ_a + (1/Re)(2φ_{a′,x}φ_{a,x} + φ_{a′,y}φ_{a,y})],
    A^{euv} = ∫(u*_y φ_{a′}φ_a + (1/Re)φ_{a′,x}φ_{a,y}), A^{evu} = ∫(v*_x φ_{a′}φ_a + (1/Re)φ_{a′,y}φ_{a,x}),
    A^{evv} likewise, B^{eup} = −∫ψ_{b′}φ_{a,x}, B^{evp} = −∫ψ_{b′}φ_{a,y}, and the right sides f^{eu}, f^{ev}, f^{ep} (minus
    the residuals of the weak equations (10.160)–(10.162) at the current guess).  Integrals by (10.198).
    Parameters: mesh; U = (u*, v*, p*) global arrays; Re; dt (None → steady: no time term); alpha_t, beta_t (10.163);
    U_old = (u(t_n), v(t_n)); dUdt_old = (∂u/∂t(t_n), ∂v/∂t(t_n)); convection (False → Stokes).
    Returns dict of arrays: Auu, Auv, Avu, Avv (T, nv, nv), Bup, Bvp (T, nv, 3), fu, fv (T, nv), fp (T, 3).
    Label: analytic.
    """
    g = _geometry(mesh)
    phi, psi, gx, gy, w = g["phi"], g["psi"], g["gx"], g["gy"], g["w"]
    el = mesh.elems
    pel = g["pel"]
    u, v, p = (U["u"], U["v"], U["p"]) if isinstance(U, dict) else U
    if isinstance(U_old, dict):
        dUdt_old = (U_old.get("dudt", np.zeros_like(u)), U_old.get("dvdt", np.zeros_like(u))) if dUdt_old is None else dUdt_old
        U_old = (U_old["u"], U_old["v"])
    ue, ve, pe = u[el], v[el], p[pel]
    us = ue @ phi  # (T, 7)
    vs = ve @ phi
    ps = pe @ psi
    ux = np.einsum("ta,tal->tl", ue, gx)
    uy = np.einsum("ta,tal->tl", ue, gy)
    vx = np.einsum("ta,tal->tl", ve, gx)
    vy = np.einsum("ta,tal->tl", ve, gy)
    cvec = 1.0 if convection else 0.0
    kt = 0.0 if dt is None else alpha_t / dt
    M = np.einsum("al,bl,tl->tab", phi, phi, w)  # ∫ φ_a φ_b
    Kxx = np.einsum("tal,tbl,tl->tab", gx, gx, w)
    Kyy = np.einsum("tal,tbl,tl->tab", gy, gy, w)
    Kxy = np.einsum("tal,tbl,tl->tab", gx, gy, w)  # ∫ φ_{a,x} φ_{b,y}  (test a, trial b)
    Kyx = np.transpose(Kxy, (0, 2, 1))  # ∫ φ_{a,y} φ_{b,x}
    if convection:
        C = np.einsum("al,tl,tbl->tab", phi, w, us[:, None, :] * gx + vs[:, None, :] * gy)  # ∫ φ_a (u*φ_{b,x} + v*φ_{b,y})
        Rux = np.einsum("al,bl,tl->tab", phi, phi, w * ux)
        Ruy = np.einsum("al,bl,tl->tab", phi, phi, w * uy)
        Rvx = np.einsum("al,bl,tl->tab", phi, phi, w * vx)
        Rvy = np.einsum("al,bl,tl->tab", phi, phi, w * vy)
    else:
        C = Rux = Ruy = Rvx = Rvy = 0.0
    Auu = kt * M + C + Rux + (2 * Kxx + Kyy) / Re  # Eq. (10.189)
    Auv = Ruy + Kyx / Re  # Eq. (10.190)
    Avu = Rvx + Kxy / Re  # Eq. (10.191)
    Avv = kt * M + C + Rvy + (Kxx + 2 * Kyy) / Re  # Eq. (10.192)
    Bup = -np.einsum("tal,bl,tl->tab", gx, psi, w)  # Eq. (10.193)
    Bvp = -np.einsum("tal,bl,tl->tab", gy, psi, w)  # Eq. (10.194)
    # right-hand sides (10.195)–(10.197)
    tu = cvec * (us * ux + vs * uy)
    tv = cvec * (us * vx + vs * vy)
    if dt is not None:
        uo, vo = U_old
        uoe, voe = uo[el] @ phi, vo[el] @ phi
        du_o, dv_o = (np.zeros_like(u), np.zeros_like(v)) if dUdt_old is None else dUdt_old
        dtu = kt * (us - uoe) - beta_t * (du_o[el] @ phi)
        dtv = kt * (vs - voe) - beta_t * (dv_o[el] @ phi)  # slip R10: β ∂v/∂t(t_n) is known data (no star)
    else:
        dtu = dtv = 0.0
    fu = (-np.einsum("al,tl->ta", phi, w * (dtu + tu)) + np.einsum("tal,tl->ta", gx, w * ps)
          - np.einsum("tal,tl->ta", gx, w * 2 * ux) / Re - np.einsum("tal,tl->ta", gy, w * (uy + vx)) / Re)  # Eq. (10.195)
    fv = (-np.einsum("al,tl->ta", phi, w * (dtv + tv)) + np.einsum("tal,tl->ta", gy, w * ps)
          - np.einsum("tal,tl->ta", gx, w * (uy + vx)) / Re - np.einsum("tal,tl->ta", gy, w * 2 * vy) / Re)  # Eq. (10.196)
    fp = np.einsum("bl,tl->tb", psi, w * (ux + vy))  # Eq. (10.197)
    return dict(Auu=Auu, Auv=Auv, Avu=Avu, Avv=Avv, Bup=Bup, Bvp=Bvp, fu=fu, fv=fv, fp=fp)


def assemble_newton_system(mesh: Mesh, state_star, state_n=None, Re: float = 1.0, dt: float | None = None,
                           alpha_t: float = 1.0, beta_t: float = 0.0, printed_10_172: bool = False, convection: bool = True,
                           dUdt_old=None, printed: bool | None = None):
    """Global block system of the Newton step (10.173): [[A_uu, A_uv, B_up], [A_vu, A_vv, B_vp], [B_upᵀ, B_vpᵀ, 0]]·(u′, v′, p′)
    = (f_u, f_v, f_p), assembled from :func:`element_newton_matrices` through the local → global node map.

    Book: §10.5, Eqs. (10.170)–(10.183) and the assembly paragraph after (10.198).
    Parameters: mesh; state_star = current guess (u*, v*, p*) as a tuple or dict; state_n = the level-n state (tuple (u, v)
    or dict(u, v, dudt, dvdt)) — None for a steady problem; Re; dt (None → steady); alpha_t, beta_t (10.163);
    printed_10_172 (the printed continuity row of (10.172) with u_{A′} in both sums — slip R6 — i.e. [B_upᵀ, B_upᵀ, 0]);
    convection (False → Stokes); dUdt_old (∂u/∂t, ∂v/∂t at t_n when state_n is a tuple).
    Returns (A (scipy.sparse CSR, (2N + V)²), f (vector)).  Label: analytic.
    """
    import scipy.sparse as sps

    printed = printed_10_172 if printed is None else printed
    e = element_newton_matrices(mesh, state_star, Re, dt, alpha_t, beta_t, state_n, dUdt_old, convection)
    N, V = mesh.N, mesh.V
    el, pel = mesh.elems, mesh.elems[:, :3]
    Auu = _scatter(el, el, e["Auu"], (N, N))
    Auv = _scatter(el, el, e["Auv"], (N, N))
    Avu = _scatter(el, el, e["Avu"], (N, N))
    Avv = _scatter(el, el, e["Avv"], (N, N))
    Bup = _scatter(el, pel, e["Bup"], (N, V))
    Bvp = _scatter(el, pel, e["Bvp"], (N, V))
    row3 = [Bup.T, Bup.T, None] if printed else [Bup.T, Bvp.T, None]  # Eq. (10.172) (printed: u_{A′} twice — slip R6)
    A = sps.bmat([[Auu, Auv, Bup], [Avu, Avv, Bvp], row3], format="csr")  # Eq. (10.173)
    f = np.concatenate([_vec_scatter(el, e["fu"], N), _vec_scatter(el, e["fv"], N), _vec_scatter(pel, e["fp"], V)])
    return A, f


def assemble_saddle(mesh: Mesh, Re: float = 1.0, state_star=None, **kw):
    """The saddle-point system (10.137) [[A, B], [Bᵀ, 0]] at a state (default zero: the Stokes/Oseen matrix) — an alias of
    :func:`assemble_newton_system`.  Book: §10.4, Eq. (10.137).  Returns (A, f).  Label: analytic."""
    if state_star is None:
        state_star = (np.zeros(mesh.N), np.zeros(mesh.N), np.zeros(mesh.V))
    return assemble_newton_system(mesh, state_star, None, Re, **kw)


@dataclass
class BoundaryConditions:
    """Essential (Dirichlet) velocity conditions: node indices and values for u and v; ``pin_pressure`` (vertex index or None)
    fixes the pressure constant when every boundary is Dirichlet (no natural outflow)."""
    nodes_u: np.ndarray
    values_u: np.ndarray
    nodes_v: np.ndarray
    values_v: np.ndarray
    pin_pressure: int | None = None


def apply_dirichlet(A, b, dofs, values, method: str = "replace", penalty: float = 1e12):
    """Impose essential conditions x_k = values_k on a linear system.

    Book: §10.5 (after (10.198)): "use the equation of the boundary condition to replace the corresponding equation" or
    "multiply a large constant by the equation of the boundary condition and add this equation" (penalty).
    Returns (A′, b′) (CSR).  Validation: V1 both methods give the same Poiseuille solution.  Label: analytic.
    """
    import scipy.sparse as sps

    dofs = np.asarray(dofs, dtype=int)
    values = np.asarray(values, dtype=float)
    A = A.tolil() if method == "replace" else A.tocsr().copy()
    b = np.array(b, dtype=float)
    if method == "replace":
        A[dofs, :] = 0.0
        A[dofs, dofs] = 1.0
        b[dofs] = values
        return A.tocsr(), b
    if method == "penalty":
        D = sps.coo_matrix((np.full(dofs.size, penalty), (dofs, dofs)), shape=A.shape)
        b[dofs] += penalty * values
        return (A + D).tocsr(), b
    raise ValueError("method must be 'replace' or 'penalty'")


def _solve(A, b, solver: str = "direct"):
    from scipy.sparse.linalg import gmres, spsolve, splu

    if solver == "direct":
        return spsolve(A.tocsc(), b)
    if solver == "gmres":
        ilu = splu(A.tocsc())
        x, info = gmres(A, b, M=None, rtol=1e-12, restart=200, maxiter=2000, x0=ilu.solve(b))
        return x
    raise ValueError("solver must be direct or gmres")


def newton_solve(mesh: Mesh, state_n=None, Re: float = 1.0, dt: float | None = None, alpha_t: float = 1.0,
                 beta_t: float = 0.0, tol: float = 1e-10, max_iter: int = 12, bc: BoundaryConditions | None = None,
                 U0=None, convection: bool = True, method: str = "replace", printed_10_172: bool = False,
                 solver: str = "direct") -> dict:
    """Newton iteration for the steady, or one fully implicit time step of the, mixed P2–P1 Navier–Stokes equations.

    Book: §10.5, Eq. (10.164) (fields = guess + correction), (10.165)–(10.167) (equations for the corrections, linearised by
    dropping u′·∇u′; their right sides are the residuals "used to monitor the convergence"), (10.173) (block system),
    essential conditions by row replacement or penalty.  The corrections at Dirichlet nodes are value − current guess.
    Parameters: mesh; state_n (the level-n state: dict(u, v, p, dudt, dvdt) or tuple (u, v, p); None = steady); Re; dt (None =
    steady); alpha_t, beta_t (10.163); tol (max|f| of the free rows); max_iter; bc (default: ``mesh.info["bc"]`` or, for a
    cylinder mesh, :func:`channel_bcs`); U0 (initial guess (u, v, p); default state_n, else zero with the boundary values);
    convection; method ("replace" | "penalty"); printed_10_172 (R6 variant); solver ("direct" | "gmres").
    Returns dict(u, v, p, iterations, residuals (list of max|f|), converged).
    Validation: V3 quadratic convergence (residual ratios); V1 Poiseuille exact; V3 Kovasznay orders 3 (u) and 2 (p).
    Label: converged.
    """
    N, V = mesh.N, mesh.V
    if bc is None:
        bc = mesh.info.get("bc")
        if bc is None and mesh.info.get("kind") == "cylinder_channel":
            bc = channel_bcs(mesh)
        if bc is None:
            raise ValueError("boundary conditions needed (bc=)")
    U_old = dUdt_old = None
    if state_n is not None:
        if isinstance(state_n, dict):
            U_old = (np.asarray(state_n["u"], float), np.asarray(state_n["v"], float))
            dUdt_old = (np.asarray(state_n.get("dudt", np.zeros(N)), float), np.asarray(state_n.get("dvdt", np.zeros(N)), float))
            pn = np.asarray(state_n.get("p", np.zeros(V)), float)
        else:
            U_old = (np.asarray(state_n[0], float), np.asarray(state_n[1], float))
            pn = np.asarray(state_n[2], float) if len(state_n) > 2 else np.zeros(V)
            dUdt_old = (np.zeros(N), np.zeros(N))
    if U0 is not None:
        u, v, p = (np.array(a, dtype=float) for a in U0)
    elif U_old is not None:
        u, v, p = U_old[0].copy(), U_old[1].copy(), pn.copy()
    else:
        u, v, p = np.zeros(N), np.zeros(N), np.zeros(V)
    u[bc.nodes_u] = bc.values_u
    v[bc.nodes_v] = bc.values_v
    pin = [] if bc.pin_pressure is None else [np.array([2 * N + bc.pin_pressure])]
    dofs = np.concatenate([bc.nodes_u, N + bc.nodes_v] + pin)
    free = np.ones(2 * N + V, dtype=bool)
    free[dofs] = False
    res, conv, it = [], False, 0
    for it in range(1, max_iter + 1):
        A, f = assemble_newton_system(mesh, (u, v, p), U_old, Re, dt, alpha_t, beta_t, printed_10_172, convection,
                                      dUdt_old=dUdt_old)
        r = float(np.max(np.abs(f[free])))
        res.append(r)
        if r < tol:
            conv = True
            it -= 1
            break
        vals = np.concatenate([bc.values_u - u[bc.nodes_u], bc.values_v - v[bc.nodes_v]]
                              + ([] if bc.pin_pressure is None else [np.array([0.0])]))
        A2, f2 = apply_dirichlet(A, f, dofs, vals, method)
        dx = _solve(A2, f2, solver)  # Eq. (10.173) for the corrections (10.164)
        u += dx[:N]
        v += dx[N:2 * N]
        p += dx[2 * N:]
        if not convection:
            A, f = assemble_newton_system(mesh, (u, v, p), U_old, Re, dt, alpha_t, beta_t, printed_10_172, convection,
                                          dUdt_old=dUdt_old)
            res.append(float(np.max(np.abs(f[free]))))
            conv = res[-1] < max(tol, 1e-9)
            break
    return dict(u=u, v=v, p=p, iterations=it, residuals=res, converged=conv)


def stokes_solve(mesh: Mesh, bc: BoundaryConditions, Re: float = 1.0, method: str = "replace", printed: bool = False,
                 lstsq: bool = False) -> dict:
    """Steady Stokes flow (convection switched off) with the same assembly — one linear solve.

    Book: §10.4 (the Stokes sub-problems of the Θ-scheme steps 1 and 3; the LBB discussion).  ``lstsq``: dense
    least-squares solve (for equal-order pairs whose saddle-point matrix is singular).  Returns dict(u, v, p).
    Label: analytic.
    """
    N, V = mesh.N, mesh.V
    U = (np.zeros(N), np.zeros(N), np.zeros(V))
    A, f = assemble_newton_system(mesh, U, None, Re, convection=False, printed_10_172=printed)
    dofs = np.concatenate([bc.nodes_u, N + bc.nodes_v] + ([] if bc.pin_pressure is None else [np.array([2 * N + bc.pin_pressure])]))
    vals = np.concatenate([bc.values_u, bc.values_v] + ([] if bc.pin_pressure is None else [np.array([0.0])]))
    A2, f2 = apply_dirichlet(A, f, dofs, vals, method)
    if lstsq:
        x = np.linalg.lstsq(A2.toarray(), f2, rcond=None)[0]
    else:
        x = _solve(A2, f2)
    p = x[2 * N:]
    return dict(u=x[:N], v=x[N:2 * N], p=p - p.mean())


def stokes_p2p1(mesh: Mesh, bc: BoundaryConditions, **kw) -> dict:
    """Alias of :func:`stokes_solve` for a P2–P1 mesh (Taylor–Hood, Fig. 10.5a).  Label: analytic."""
    return stokes_solve(mesh, bc, **kw)


def stokes_cavity(n: int = 4, pair: str = "P2P1", lid: float = 1.0) -> dict:
    """Stokes lid-driven cavity on the unit square with the chosen element pair (the explainer's LBB demonstration).

    Book: §10.4 (spurious pressure with equal-order interpolation, Fig. 10.5; LBB).  Mesh :func:`structured_square_mesh`
    (n, pair); velocity Dirichlet on the whole boundary (lid u = 1 on y = 1, the two top corners u = 0); one pressure pinned,
    then the mean removed.  Solved by dense least squares (small n; the P1–P1 matrix is singular).
    Returns dict(p (at the pressure nodes = vertices), u, v (at the velocity nodes), p_range, p_std_checker (amplitude of the
    (−1)^{i+j} component of p), n_u, n_p (unknown counts after the Dirichlet nodes are removed), residual (‖A x − b‖),
    rank_deficiency, mesh).  Label: analytic.
    """
    m = structured_square_mesh(n, pair)
    x, y = m.nodes[:, 0], m.nodes[:, 1]
    allb = m.boundary["all"]
    lidnodes = allb[(np.abs(y[allb] - 1.0) < 1e-12) & (x[allb] > 1e-12) & (x[allb] < 1 - 1e-12)]
    vals_u = np.where(np.isin(allb, lidnodes), lid, 0.0)
    N, V = m.N, m.V
    A, f = assemble_newton_system(m, (np.zeros(N), np.zeros(N), np.zeros(V)), None, 1.0, convection=False)
    dofs = np.concatenate([allb, N + allb, [2 * N]])
    A2, f2 = apply_dirichlet(A, f, dofs, np.concatenate([vals_u, np.zeros(allb.size), [0.0]]))
    Ad = A2.toarray()
    xsol = np.linalg.lstsq(Ad, f2, rcond=None)[0]
    sv = np.linalg.svd(Ad, compute_uv=False)
    rank_def = int(np.sum(sv < 1e-10 * sv.max()))
    p = xsol[2 * N:]
    p = p - p.mean()
    vx, vy = m.nodes[:V, 0], m.nodes[:V, 1]
    sgn = (-1.0) ** (np.round(vx * n).astype(int) + np.round(vy * n).astype(int))
    n_u = 2 * (N - allb.size)
    return dict(p=p, u=xsol[:N], v=xsol[N:2 * N], p_range=float(p.max() - p.min()),
                p_std_checker=float(abs(np.sum(p * sgn)) / V), n_u=int(n_u), n_p=int(V - 1),
                residual=float(np.linalg.norm(Ad @ xsol - f2)), rank_deficiency=rank_def, mesh=m)


# ======================================================================================================================
# Exact test fields
# ======================================================================================================================
def kovasznay_exact(x, y, Re: float = 40.0):
    """Kovasznay's exact steady Navier–Stokes flow behind a grid (Kovasznay 1948):
    u = 1 − e^{λx} cos 2πy, v = (λ/2π) e^{λx} sin 2πy, p = ½(1 − e^{2λx}), λ = Re/2 − √(Re²/4 + 4π²).

    Book: §10.4–10.5 (a code-verification field of ours for the mixed FE solver; not in the book).  Form cross-checked with
    Wikipedia "Kovasznay flow" (2026-09-30); the pressure is verified by a sympy residual in the tests.  Returns (u, v, p).
    Label: analytic.
    """
    lam = 0.5 * Re - np.sqrt(0.25 * Re ** 2 + 4 * np.pi ** 2)
    ex = np.exp(lam * np.asarray(x, dtype=float))
    y = np.asarray(y, dtype=float)
    return 1 - ex * np.cos(2 * np.pi * y), lam / (2 * np.pi) * ex * np.sin(2 * np.pi * y), 0.5 * (1 - ex ** 2)


def _l2_error(mesh: Mesh, vals, exact_fn, kind: str = "velocity"):
    g = _geometry(mesh)
    el = mesh.elems if kind == "velocity" else mesh.elems[:, :3]
    sh = g["phi"] if kind == "velocity" else g["psi"]
    fq = vals[el] @ sh
    ex = exact_fn(g["xq"], g["yq"])
    return float(np.sqrt(np.sum(g["w"] * (fq - ex) ** 2)))


def _kovasznay_errors(n, Re, box):
    x0, x1, y0, y1 = box
    m = structured_square_mesh(n, "P2P1", x0, x1, y0, y1)
    b = m.boundary["all"]
    ue, ve, _ = kovasznay_exact(m.nodes[b, 0], m.nodes[b, 1], Re)
    bc = BoundaryConditions(b, ue, b, ve, pin_pressure=0)
    r = newton_solve(m, None, Re, bc=bc, U0=(np.ones(m.N), np.zeros(m.N), np.zeros(m.V)), max_iter=20)
    g = _geometry(m)
    eu = _l2_error(m, r["u"], lambda X, Y: kovasznay_exact(X, Y, Re)[0])
    ev = _l2_error(m, r["v"], lambda X, Y: kovasznay_exact(X, Y, Re)[1])
    pq = r["p"][m.elems[:, :3]] @ g["psi"]
    pex = kovasznay_exact(g["xq"], g["yq"], Re)[2]
    shift = np.sum(g["w"] * (pex - pq)) / np.sum(g["w"])
    ep = float(np.sqrt(np.sum(g["w"] * (pq + shift - pex) ** 2)))
    return float(np.hypot(eu, ev)), ep, r, m


def kovasznay_test(n: int = 8, Re: float = 40.0, box=(-0.5, 1.0, -0.5, 1.5)) -> dict:
    """Solve the steady NS equations with Kovasznay's velocity on the whole boundary and compare with the exact flow.

    Book: §10.5 (the Newton–mixed-FE machinery (10.164)–(10.198) on a problem with a known answer).  n × n squares on the
    box; also solved on n/2 squares to report the observed orders.
    Returns dict(err_u (L² of u and v combined), err_p (L² after matching the mean), orders (dict u, p from n/2 → n),
    iterations, residuals, h, mesh, sol).  Validation: V1/V3 orders ≈ 3 (velocity) and ≈ 2 (pressure).  Label: converged.
    """
    eu, ep, r, m = _kovasznay_errors(n, Re, box)
    eu2, ep2, _, _ = _kovasznay_errors(max(2, n // 2), Re, box)
    return dict(err_u=eu, err_p=ep, orders=dict(u=float(np.log2(eu2 / eu)), p=float(np.log2(ep2 / ep))),
                iterations=r["iterations"], residuals=r["residuals"], h=(box[1] - box[0]) / n, mesh=m, sol=r)


def poiseuille_test(n: int = 4, Re: float = 10.0, method: str = "replace", printed: bool = False, G: float = 8.0) -> dict:
    """Plane Poiseuille flow through the unit square (walls y = 0, 1; the exact velocity prescribed on the whole boundary;
    one pressure value pinned) — P2 velocity and P1 pressure reproduce it exactly (u quadratic, p linear).

    Book: §10.5 (code verification of the FE machinery; the ch04/ch08 exact field): u = (Re G/2) y(1 − y), v = 0,
    p = −Gx + const.  (With a traction-free outflow instead, the symmetric-gradient weak form (10.156) would impose
    u_y + v_x = 0 at the outlet, which Poiseuille flow does not satisfy — so the whole boundary is Dirichlet here.)
    ``printed``: the R6 continuity row of (10.172) (the test must fail).  ``method``: "replace" | "penalty".
    Returns dict(err_u, err_v (max nodal), err_p (max after removing the mean), iterations, residuals).  Label: analytic.
    """
    m = structured_square_mesh(n, "P2P1")
    y = m.nodes[:, 1]
    uex = 0.5 * Re * G * y * (1 - y)
    dn = m.boundary["all"]
    bc = BoundaryConditions(dn, uex[dn], dn, np.zeros(dn.size), pin_pressure=0)
    r = newton_solve(m, None, Re, bc=bc, U0=(uex * 0, uex * 0, np.zeros(m.V)), method=method, printed_10_172=printed,
                     max_iter=10)
    pex = -G * m.nodes[:m.V, 0]
    ep = float(np.max(np.abs((r["p"] - r["p"].mean()) - (pex - pex.mean()))))
    return dict(err_u=float(np.max(np.abs(r["u"] - uex))), err_v=float(np.max(np.abs(r["v"]))), err_p=ep,
                iterations=r["iterations"], residuals=r["residuals"], p=r["p"], mesh=m)


# ======================================================================================================================
# LBB: the discrete inf–sup constant
# ======================================================================================================================
def infsup_constant(n, pair: str = "P2P1") -> dict:
    """Discrete inf–sup (LBB) constant β_h of the element pair on a mesh with Dirichlet velocity on the whole boundary.

    Book: §10.4 — the Babuška–Brezzi (LBB, inf–sup) condition that the mixed pair must satisfy (stated in words; Oden & Carey).
    β_h² = smallest nonzero eigenvalue λ of Bᵀ A⁻¹ B q = λ M_p q, with A the vector Laplacian (H¹ seminorm) on the interior
    velocity nodes, B the divergence coupling −∫ψ ∇·φ, M_p the pressure mass matrix; the constant pressure (λ = 0) is
    skipped.  For P2–P1 β_h stays bounded as the mesh is refined; for P1–P1 extra (spurious) zero or tiny eigenvalues appear.
    Dense generalised eigenproblem (scipy.linalg.eigh) — fine for ≤ ~1000 pressure nodes.
    Parameters: n (squares per side of :func:`structured_square_mesh`) or a :class:`Mesh`; pair (with an integer n).
    Returns dict(beta (square root of the smallest eigenvalue beyond the constant; 0 when spurious modes exist),
    beta_nonspurious (the smallest one above the spurious ones), n_spurious (eigenvalues < 1e-8 besides the constant),
    eigenvalues (first 8), n_u, n_p (velocity / pressure unknowns), N_u, N_p).  Label: converged.
    """
    from scipy.linalg import eigh

    mesh = n if isinstance(n, Mesh) else structured_square_mesh(int(n), pair)
    g = _geometry(mesh)
    el = mesh.elems
    pel = el[:, :3]
    N, V = mesh.N, mesh.V
    K = _scatter(el, el, np.einsum("tal,tbl,tl->tab", g["gx"], g["gx"], g["w"]) + np.einsum("tal,tbl,tl->tab", g["gy"], g["gy"], g["w"]), (N, N))
    Bx = _scatter(el, pel, -np.einsum("tal,bl,tl->tab", g["gx"], g["psi"], g["w"]), (N, V))
    By = _scatter(el, pel, -np.einsum("tal,bl,tl->tab", g["gy"], g["psi"], g["w"]), (N, V))
    Mp = _scatter(pel, pel, np.einsum("al,bl,tl->tab", g["psi"], g["psi"], g["w"]), (V, V)).toarray()
    bnd = mesh.boundary.get("all")
    if bnd is None:
        bnd = np.unique(np.concatenate(list(mesh.boundary.values())))
    inner = np.setdiff1d(np.arange(N), bnd)
    Ki = K[inner][:, inner].tocsc()
    from scipy.sparse.linalg import splu

    lu = splu(Ki)
    Bxi = Bx[inner].toarray()
    Byi = By[inner].toarray()
    S = Bxi.T @ lu.solve(Bxi) + Byi.T @ lu.solve(Byi)
    S = 0.5 * (S + S.T)
    lam = eigh(S, Mp, eigvals_only=True)
    lam = np.sort(np.abs(lam))
    nz = lam[1:]  # skip the constant pressure mode
    tiny = int(np.sum(nz < 1e-8 * max(1.0, lam.max())))
    beta = float(np.sqrt(nz[tiny])) if tiny < nz.size else 0.0
    b0 = 0.0 if tiny > 0 else float(np.sqrt(max(nz[0], 0.0)))
    return dict(beta=b0, beta_nonspurious=beta, eigenvalues=lam[:8], n_spurious=tiny, n_u=int(2 * inner.size), n_p=int(V),
                N_u=int(2 * inner.size), N_p=int(V))


def infsup_table(pairs=("P1P1", "P2P1"), n_list=(2, 3, 4, 6, 8, 12, 16), diagonal: str = "right") -> dict:
    """β_h for each pair and mesh size on the unit square (the E8 / F7 table; written to reference/ch10/infsup_table.csv by
    the scripts).  Returns dict(pair → list of dict(n, h, beta, beta_nonspurious, n_spurious, n_u, n_p)).  Book: §10.4 (LBB).
    Label: converged."""
    out = {}
    for pr in pairs:
        rows = []
        for n in n_list:
            r = infsup_constant(structured_square_mesh(int(n), pr, diagonal=diagonal))
            rows.append(dict(n=int(n), h=1.0 / n, beta=r["beta"], beta_nonspurious=r["beta_nonspurious"],
                             n_spurious=r["n_spurious"], n_u=r["n_u"], n_p=r["n_p"]))
        out[pr] = rows
    return out


# ======================================================================================================================
# Cylinder in a channel (§10.5)
# ======================================================================================================================
def channel_bcs(mesh: Mesh, U: float = 1.0) -> BoundaryConditions:
    """Boundary conditions of Fig. 10.15: inflow Γ₁ u = U, v = 0; sliding walls Γ₃, Γ₄ u = U, v = 0; cylinder Γ₅ u = v = 0;
    outflow Γ₂ natural (traction-free in the symmetric-gradient weak form (10.156) — the "zero normal stress" of (10.83)).
    Book: §10.5.  Returns :class:`BoundaryConditions` (no pressure pin: the outflow fixes the constant).  Label: analytic.
    """
    b = mesh.boundary
    moving = np.unique(np.concatenate([b["inflow"], b["top"], b["bottom"]]))
    cyl = b["cylinder"]
    nodes = np.concatenate([moving, cyl])
    vals = np.concatenate([np.full(moving.size, U), np.zeros(cyl.size)])
    return BoundaryConditions(nodes, vals, nodes, np.zeros(nodes.size), None)


def cylinder_forces(mesh: Mesh, sol, Re: float, dt: float | None = None, alpha_t: float = 1.0, beta_t: float = 0.0,
                    U_old=None, dUdt_old=None) -> dict:
    """Drag, lift and torque coefficients of the cylinder from the consistent ("reaction") nodal forces.

    Book: §10.5, Fig. 10.21 (drag, lift and torque histories; C_D, C_L as in (4.107)–(4.108) per unit span, torque normalised by
    ½ρU²d²).  The weak momentum residual of the converged solution, evaluated with the test functions of the cylinder nodes,
    equals the traction the wall exerts on the fluid; its sum is minus the force on the body (the ch04 control-volume idea in
    discrete form — more accurate than differentiating the velocity at the wall).
    Returns dict(CD, CL, CM).  Label: qualitative (demonstration; no mesh-convergence order asserted).
    """
    e = element_newton_matrices(mesh, (sol["u"], sol["v"], sol["p"]), Re, dt, alpha_t, beta_t, U_old, dUdt_old)
    # (U_old may be a dict(u, v, dudt, dvdt) too)
    N = mesh.N
    Rx = -_vec_scatter(mesh.elems, e["fu"], N)  # residual = −f  (weak momentum ∫(…)N_A + …)
    Ry = -_vec_scatter(mesh.elems, e["fv"], N)
    cyl = mesh.boundary["cylinder"]
    Fx = -np.sum(Rx[cyl])
    Fy = -np.sum(Ry[cyl])
    xc, yc = mesh.nodes[cyl, 0], mesh.nodes[cyl, 1]
    Mz = -np.sum(xc * Ry[cyl] - yc * Rx[cyl])
    return dict(CD=float(2 * Fx), CL=float(2 * Fy), CM=float(2 * Mz))


def march_unsteady(mesh: Mesh, bc: BoundaryConditions, Re: float, dt: float, nsteps: int, U0=None, alpha_t: float = 2.0,
                   beta_t: float = 1.0, max_iter: int = 4, tol: float = 1e-8, record=None, kick: float = 0.0,
                   kick_time: float = 0.0, start_be: int = 1, verbose: bool = False) -> dict:
    """Fully implicit time marching of the mixed-FE Navier–Stokes equations with Newton at every step.

    Book: §10.5, (10.163) time derivative (α = 2, β = 1 second order by default; the first ``start_be`` steps use α = 1,
    β = 0 because ∂u/∂t(t_0) is unknown), (10.164)–(10.167) Newton, fully implicit.
    Parameters: mesh, bc, Re, dt, nsteps, U0 (initial state, default the steady Stokes-like start from rest with the boundary
    values), alpha_t, beta_t, max_iter, tol, record(k, t, sol) callback, kick (amplitude of a short rotation of the cylinder
    to break the symmetry, a smooth sin² pulse of the rotation rate — our device to shorten the transient; 0 = none),
    kick_time (duration).
    Returns dict(t, CD, CL, CM (histories), sol (final), newton_its).
    Label: qualitative (demonstration; no time- or mesh-convergence order asserted).
    """
    N, V = mesh.N, mesh.V
    if U0 is None:
        u, v, p = np.zeros(N), np.zeros(N), np.zeros(V)
    else:
        u, v, p = (np.array(a, float) for a in U0)
    du, dv = np.zeros(N), np.zeros(N)
    ts, cd, cl, cm, its = [], [], [], [], []
    cyl = mesh.boundary["cylinder"]
    base_u = bc.values_u.copy()
    base_v = bc.values_v.copy()
    is_cyl_u = np.isin(bc.nodes_u, cyl)
    for k in range(1, nsteps + 1):
        t = k * dt
        a_t, b_t = (1.0, 0.0) if k <= start_be else (alpha_t, beta_t)
        if kick and t <= kick_time:
            xc, yc = mesh.nodes[bc.nodes_u, 0], mesh.nodes[bc.nodes_u, 1]
            om = kick * np.sin(np.pi * t / kick_time) ** 2  # smooth rotation-rate pulse (keeps ∂u/∂t consistent)
            bc.values_u = np.where(is_cyl_u, -om * yc, base_u)
            bc.values_v = np.where(is_cyl_u, om * xc, base_v)
        else:
            bc.values_u, bc.values_v = base_u, base_v
        guess = (u + dt * du, v + dt * dv, p) if k > 1 else (u, v, p)
        r = newton_solve(mesh, dict(u=u, v=v, p=p, dudt=du, dvdt=dv), Re, dt, a_t, b_t, tol=tol, max_iter=max_iter,
                         bc=bc, U0=guess)
        f = cylinder_forces(mesh, r, Re, dt, a_t, b_t, (u, v), (du, dv))
        du = time_derivative(r["u"], u, du, dt, a_t, b_t)
        dv = time_derivative(r["v"], v, dv, dt, a_t, b_t)
        u, v, p = r["u"], r["v"], r["p"]
        ts.append(t)
        cd.append(f["CD"])
        cl.append(f["CL"])
        cm.append(f["CM"])
        its.append(r["iterations"])
        if record is not None:
            record(k, t, dict(u=u, v=v, p=p))
        if verbose and k % 50 == 0:
            print(f"  step {k}  t = {t:.2f}  CD = {f['CD']:.4f}  CL = {f['CL']:+.4f}  Newton its {r['iterations']}", flush=True)
    bc.values_u, bc.values_v = base_u, base_v
    return dict(t=np.array(ts), CD=np.array(cd), CL=np.array(cl), CM=np.array(cm), sol=dict(u=u, v=v, p=p),
                newton_its=np.array(its))


def nodal_vorticity(mesh: Mesh, u, v) -> np.ndarray:
    """ω = v_x − u_y at the vertices (area-weighted average of the element values at the quadrature points) — for plotting
    Fig. 10.19-type contours.  Label: analytic."""
    g = _geometry(mesh)
    el = mesh.elems
    w = np.einsum("ta,tal->tl", v[el], g["gx"]) - np.einsum("ta,tal->tl", u[el], g["gy"])
    wa = np.sum(w * g["w"], axis=1)
    area = np.sum(g["w"], axis=1)
    num = np.zeros(mesh.V)
    den = np.zeros(mesh.V)
    for k in range(3):
        np.add.at(num, el[:, k], wa)
        np.add.at(den, el[:, k], area)
    return num / den


def cylinder_steady(Re: float = 40.0, mesh_level: str = "coarse", cache: bool = True, U0=None) -> dict:
    """Steady flow past the cylinder in the channel of :func:`cylinder_channel_mesh` by Newton (our geometry W = 5,
    x ∈ [−8, 16]; sliding walls and uniform inflow U = 1; traction-free outflow).

    Book: §10.5, Figs. 10.18–10.19 (steady Re = 1, 10, 40; the confined configuration), (10.156)–(10.198).
    ``mesh_level``: "coarse" (n_theta 32, h_far 0.5 — fast), "medium" (48, 0.4), "fine" (64, 0.3, wake spacing 0.2).
    Cached in outputs/ch10 (npz keyed by Re and the level).  Returns dict(mesh, u, v, p, residuals, iterations, converged,
    CD, CL, CM, omega (vertex vorticity)).  C_D and C_L are for the confined sliding-wall case — qualitative only against
    unbounded data (slip R9).  Label: qualitative (demonstration; Newton converges, no mesh-convergence order asserted).
    """
    levels = dict(coarse=dict(n_theta=32, h_far=0.5, n_r=10), medium=dict(n_theta=48, h_far=0.4, n_r=12),
                  fine=dict(n_theta=64, h_far=0.3, n_r=16, h_wake=0.2))
    mesh = cylinder_channel_mesh(**levels[mesh_level])
    bc = channel_bcs(mesh)
    path = None
    if cache:
        from .project import repo_root

        path = repo_root() / "outputs" / "ch10" / f"fe_cylinder_steady_Re{Re:g}_{mesh_level}.npz"
        if path.is_file():
            d = np.load(path)
            sol = dict(u=d["u"], v=d["v"], p=d["p"])
            if sol["u"].size == mesh.N:
                f = cylinder_forces(mesh, sol, Re)
                return dict(mesh=mesh, **sol, residuals=list(d["residuals"]), iterations=int(d["iterations"]),
                            converged=bool(d["converged"]), omega=nodal_vorticity(mesh, sol["u"], sol["v"]), cached=str(path), **f)
    if U0 is None and Re > 20:  # continuation from a lower Re (Newton needs a good start at higher Re)
        U0 = tuple(cylinder_steady(min(Re, 10.0), mesh_level, cache)[k] for k in ("u", "v", "p")) if Re > 10 else None
    r = newton_solve(mesh, None, Re, bc=bc, U0=U0, max_iter=20)
    f = cylinder_forces(mesh, r, Re)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, u=r["u"], v=r["v"], p=r["p"], residuals=np.array(r["residuals"]), iterations=r["iterations"],
                 converged=r["converged"])
    return dict(mesh=mesh, u=r["u"], v=r["v"], p=r["p"], residuals=r["residuals"], iterations=r["iterations"],
                converged=r["converged"], omega=nodal_vorticity(mesh, r["u"], r["v"]), **f)
