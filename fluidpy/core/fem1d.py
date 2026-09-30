"""One-dimensional Galerkin finite elements with piecewise-linear "hat" functions for the transport problem (10.1)–(10.3).

Book: Kundu, Cohen & Dowling 5e, Ch. 10 §10.3, Eqs. (10.32)–(10.78) (weak form, Galerkin, hat functions, M ḋ + K d = F,
element matrices and assembly).  Equations transcribed from the page images chapters/pages/ch10/p456–p463.

Conventions.  Nodes x_0 = 0 < x_1 < … < x_n = L (node-based mesh); element e (1-based, e = 1 … n_el) is Ω^e = [x_{e−1}, x_e]
and its local nodes a = 1, 2 map to the global nodes A = e − 1, e (Eq. (10.78)).  Node 0 carries the Dirichlet value g
(essential condition, eliminated from the unknowns); the unknowns are d_A = T^h(x_A), A = 1 … n (Eq. (10.62)).  Arrays of
unknowns are 0-based, so unknown index k ↔ global node A = k + 1.
Printed slips handled here: R1 (10.67)–(10.68) index shift in the element slopes (``shape_slopes(printed=True)`` keeps the
book's labels); R2 the sentence after (10.73) ("A = e or e + 1") — the nonzero entries have A, B ∈ {e − 1, e}.
Units: x, h, L [m]; u [m/s]; D [m²/s]; T in any unit (K in the heat interpretation); q [unit of T per m].
Reused by Ch. 13 (Ekman-layer ODEs with variable eddy viscosity) and Ch. 16.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from ._util import as_scalar_if_0d, require_nonnegative

__all__ = [
    "hat", "hat_basis", "interpolate", "parent_shapes", "map_to_element", "map_to_parent", "shape_slopes",
    "element_matrices_linear", "element_integrals", "element_force", "connectivity", "assemble_1d", "assembly_trace",
    "solve_transport", "solve_steady", "interior_stencil", "bilinear_form", "weak_residual", "galerkin_equations_sympy",
    "gauss_legendre",
]


def gauss_legendre(n: int):
    """Gauss–Legendre points and weights on [−1, 1] (exact for polynomials of degree ≤ 2n − 1).  Label: analytic."""
    return np.polynomial.legendre.leggauss(int(n))


def hat(x, x_nodes, A: int):
    """Global piecewise-linear shape function N_A(x) (the "hat").

    Book: §10.3, Eq. (10.59) (interior nodes, rising on [x_{A−1}, x_A], falling on [x_A, x_{A+1}], zero elsewhere),
    Eq. (10.60) (the last node N_n, rising only), Eq. (10.61) (the Dirichlet node N_0, falling only; N_0(0) = 1 (10.46)).
    N_A(x_B) = δ_AB and Σ_A N_A ≡ 1 (partition of unity).
    Parameters: x [m] (array), x_nodes [m] (increasing, n + 1 nodes), A ∈ {0 … n}.
    Returns N_A(x) (non-dimensional).  Label: analytic.
    """
    x = np.asarray(x, dtype=float)
    xn = np.asarray(x_nodes, dtype=float)
    n = xn.size - 1
    out = np.zeros_like(x)
    if A > 0:
        m = (x >= xn[A - 1]) & (x <= xn[A])
        out[m] = (x[m] - xn[A - 1]) / (xn[A] - xn[A - 1])  # Eq. (10.59) rising part / (10.60)
    if A < n:
        m = (x >= xn[A]) & (x <= xn[A + 1])
        out[m] = (xn[A + 1] - x[m]) / (xn[A + 1] - xn[A])  # Eq. (10.59) falling part / (10.61)
    return as_scalar_if_0d(out)


def hat_basis(x_nodes, x) -> np.ndarray:
    """All hat functions N_0 … N_n evaluated at x: array of shape (len(x), n + 1).

    Book: §10.3, Eqs. (10.43)–(10.47), (10.59)–(10.61); column 0 is N_0 (the lifting g^h = gN_0 (10.47)).
    Validation: V1 partition of unity Σ_A N_A = 1 and N_A(x_B) = δ_AB.  Label: analytic.
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    return np.array([np.atleast_1d(hat(x, x_nodes, A)) for A in range(len(x_nodes))]).T


def interpolate(d, g: float, x_nodes, x):
    """Finite-element field T^h(x) = Σ_{A=1}^n d_A N_A(x) + g N_0(x).

    Book: §10.3, Eq. (10.49) (with (10.48) v^h = Σ d_A N_A); by (10.62) d_A = T^h(x_A).
    Parameters: d (n unknowns, nodes 1 … n), g (Dirichlet value at x = 0), x_nodes [m], x [m].  Returns T^h(x).
    Validation: V1 T^h(x_A) = d_A.  Label: analytic.
    """
    N = hat_basis(x_nodes, x)
    return as_scalar_if_0d(N[:, 1:] @ np.asarray(d, dtype=float) + g * N[:, 0])  # Eq. (10.49)


def parent_shapes(xi):
    """Standard shape functions on the parent element ξ ∈ [−1, 1]: N₁ = ½(1 − ξ), N₂ = ½(1 + ξ).

    Book: §10.3, Eq. (10.64) (Fig. 10.3).  Returns (N1, N2).  Label: analytic.
    """
    xi = np.asarray(xi, dtype=float)
    return as_scalar_if_0d(0.5 * (1.0 - xi)), as_scalar_if_0d(0.5 * (1.0 + xi))  # Eq. (10.64)


def map_to_element(xi, xa: float, xb: float):
    """Map ξ ∈ [−1, 1] to x ∈ [x_a, x_b]: x(ξ) = N₁x₁ᵉ + N₂x₂ᵉ = ½[(x_b − x_a)ξ + x_b + x_a].

    Book: §10.3, Eq. (10.65) with x₁ᵉ = x_{A−1}, x₂ᵉ = x_A.  Returns x [m].  Label: analytic.
    """
    xi = np.asarray(xi, dtype=float)
    return as_scalar_if_0d(0.5 * ((xb - xa) * xi + xb + xa))  # Eq. (10.65)


def map_to_parent(x, xa: float, xb: float):
    """Inverse map ξ(x) = (2x − x_b − x_a)/(x_b − x_a).  Book: §10.3, Eq. (10.66).  Label: analytic."""
    x = np.asarray(x, dtype=float)
    return as_scalar_if_0d((2.0 * x - xb - xa) / (xb - xa))  # Eq. (10.66)


def shape_slopes(xa: float, xb: float, printed: bool = False, A: int | None = None) -> dict:
    """Slopes of the two global shape functions that live on the element [x_a, x_b] = [x_{A−1}, x_A], by the chain rule.

    Book: §10.3, Eqs. (10.67)–(10.68): dN/dx = (dN/dξ)(dξ/dx) with dξ/dx = 2/(x_A − x_{A−1}) from (10.66), dN₁/dξ = −½,
    dN₂/dξ = +½.  Correct labels (consistent with (10.64): N₁ ↔ N_{A−1}, N₂ ↔ N_A): dN_{A−1}/dx = −1/hᵉ, dN_A/dx = +1/hᵉ.
    ``printed=True`` returns the book's printed labels (slip R1): N_A with −1/hᵉ and N_{A+1} with +1/hᵉ — an index shift
    (N_{A+1} is zero on this element), which a chain-rule test against :func:`hat` must reject.
    Parameters: xa, xb [m]; printed; A (optional global index of the right node → numeric ``nodes``).
    Returns dict(labels (tuple of str), slopes (tuple [1/m]), nodes (global indices, if A given)).  Label: analytic.
    """
    h = xb - xa
    dxi_dx = 2.0 / h  # from Eq. (10.66)
    s1, s2 = -0.5 * dxi_dx, 0.5 * dxi_dx  # Eqs. (10.67)–(10.68): −1/h, +1/h
    labels = ("A", "A+1") if printed else ("A-1", "A")
    out = dict(labels=labels, slopes=(s1, s2))
    if A is not None:
        out["nodes"] = (A, A + 1) if printed else (A - 1, A)
    return out


def element_matrices_linear(h: float, u: float, D: float):
    """Element mass and stiffness matrices of a linear element of length h (closed form, our D14).

    Book: §10.3, Eqs. (10.74)–(10.76): m_ab = ∫N_aN_b dx, k_ab = u∫N_{b,x}N_a dx + D∫N_{b,x}N_{a,x} dx (the book never writes
    the numbers).  m = (h/6)[[2, 1], [1, 2]];  k = (u/2)[[−1, 1], [−1, 1]] + (D/h)[[1, −1], [−1, 1]] (row a, column b; the
    convective block is not symmetric).
    Parameters: h [m] (> 0), u [m/s], D [m²/s].  Returns (m [m], k [m/s]) as 2 × 2 arrays.
    Validation: V1 equals 2-point Gauss integration :func:`element_integrals` to 1e-15.  Label: analytic.
    """
    m = h / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])  # Eq. (10.75)
    k = 0.5 * u * np.array([[-1.0, 1.0], [-1.0, 1.0]]) + D / h * np.array([[1.0, -1.0], [-1.0, 1.0]])  # Eq. (10.76)
    return m, k


def element_integrals(xa: float, xb: float, u: float, D: float, quad: int = 2):
    """Element matrices (10.75)–(10.76) by Gauss–Legendre quadrature on the parent element (dx = (h/2)dξ).

    Book: §10.3, Eqs. (10.71)–(10.76) with the map (10.65) and slopes (10.67)–(10.68).  quad = number of Gauss points (2 is
    exact for the quadratic integrand N_aN_b).  Returns (m, k).  Label: analytic.
    """
    xg, wg = gauss_legendre(quad)
    h = xb - xa
    N = np.array(parent_shapes(xg))  # (2, nq)
    dN = np.array([-1.0 / h, 1.0 / h])  # corrected (10.67)–(10.68)
    jac = 0.5 * h
    m = np.einsum("aq,bq,q->ab", N, N, wg) * jac
    k = u * np.einsum("b,aq,q->ab", dN, N, wg) * jac + D * np.outer(dN, dN) * np.sum(wg) * jac
    return m, k


def element_force(e: int, n_el: int, h: float, u: float, D: float, g: float = 0.0, q: float = 0.0) -> np.ndarray:
    """Element force vector fᵉ of Eq. (10.77): −g k¹_{a1} for e = 1 (the Dirichlet column), 0 for interior elements,
    Dq δ_{a2} for e = n_el (the natural Neumann term); both when there is a single element.

    Book: §10.3, Eq. (10.77).  Parameters: e (1-based), n_el, h [m], u [m/s], D [m²/s], g, q.  Returns ndarray(2).
    Label: analytic.
    """
    f = np.zeros(2)
    _, k = element_matrices_linear(h, u, D)
    if e == 1:
        f += -g * k[:, 0]  # Eq. (10.77), e = 1
    if e == n_el:
        f[1] += D * q  # Eq. (10.77), e = n_el
    return f


def connectivity(n_el: int, printed: bool = False) -> np.ndarray:
    """Local → global node map A(e, a): A = e − 1 for a = 1 and A = e for a = 2 (e = 1 … n_el).

    Book: §10.3, Eq. (10.78).  Returns int array (n_el, 2) (row e − 1 = [e − 1, e]); node 0 is the Dirichlet node.
    ``printed=True``: the text after (10.73) ("A = e or e + 1", slip R2) — rows [e, e + 1], which run past the last node
    n_el and never touch node 0.  Label: analytic.
    """
    e = np.arange(1, int(n_el) + 1)
    if printed:
        return np.stack([e, e + 1], axis=1)  # as printed (R2)
    return np.stack([e - 1, e], axis=1)  # Eq. (10.78)


def assemble_1d(x_nodes, u: float, D: float, g: float = 0.0, q: float = 0.0, lumped: bool = False,
                dirichlet_right: float | None = None, sparse: bool = True, right: str | None = None, T_L: float | None = None):
    """Assemble M ḋ + K d = F for the transport problem by scatter-adding element blocks.

    Book: §10.3, Eqs. (10.54)–(10.58) (global definitions), (10.69)–(10.73) (sums of element contributions), (10.74)–(10.77)
    (element blocks; f¹ = −g k¹_{a1} carries the Dirichlet column, f^{n_el} = Dq δ_{a2}), (10.78) (node map).
    Parameters
    ----------
    x_nodes : node positions [m], x_0 = 0 … x_n = L (need not be uniform).   u [m/s], D [m²/s].
    g : Dirichlet value T(0, t) = g (10.2).   q : Neumann slope ∂T/∂x(L, t) = q (10.2) — natural condition, enters F.
    lumped : row-sum ("lumped") mass matrix instead of the consistent one (ours, for comparison with FD).
    dirichlet_right : None (the book's Neumann end) or a value T_L: T(L) = T_L, used for the steady layer (10.84)–(10.85);
        node n is then eliminated like node 0.  (``right``/``T_L`` are accepted as aliases.)
    sparse : return scipy.sparse CSR (True) or dense arrays.
    Returns (M [m], K [m/s], F [unit of T · m/s]) for the unknowns d_1 … d_n (or d_1 … d_{n−1} with a Dirichlet right end).
    Validation: V1 interior rows equal (10.63) × h (:func:`interior_stencil`); V1 steady FE = centred FD on uniform meshes
    (1e-13).  Label: analytic.
    """
    import scipy.sparse as sps

    if dirichlet_right is not None:
        right, T_L = "dirichlet", dirichlet_right
    right = right or "neumann"
    xn = np.asarray(x_nodes, dtype=float)
    n = xn.size - 1
    conn = connectivity(n)
    rows, cols, mv, kv = [], [], [], []
    Ffull = np.zeros(n + 1)
    Mfull_rows = []
    for e in range(1, n + 1):
        A = conn[e - 1]
        h = xn[A[1]] - xn[A[0]]
        m, k = element_matrices_linear(h, u, D)
        if lumped:
            m = np.diag(m.sum(axis=1))
        for a in range(2):
            for b in range(2):
                rows.append(A[a])
                cols.append(A[b])
                mv.append(m[a, b])
                kv.append(k[a, b])
    del Mfull_rows
    Mf = sps.coo_matrix((mv, (rows, cols)), shape=(n + 1, n + 1)).tocsr()  # duplicates summed: Eq. (10.69)
    Kf = sps.coo_matrix((kv, (rows, cols)), shape=(n + 1, n + 1)).tocsr()
    if right == "neumann":
        Ffull[n] += D * q  # Eq. (10.77): f^{n_el} = D q δ_{a2}  (natural boundary condition)
        keep = np.arange(1, n + 1)
    elif right == "dirichlet":
        keep = np.arange(1, n)
    else:
        raise ValueError("right must be 'neumann' or 'dirichlet'")
    M = Mf[keep][:, keep]
    K = Kf[keep][:, keep]
    F = Ffull[keep] - g * Kf[keep, 0].toarray().ravel()  # Eq. (10.57)/(10.77): −g a(N_A, N_0) = −g k¹_{a1}
    if right == "dirichlet":
        TL = 0.0 if T_L is None else float(T_L)
        F = F - TL * Kf[keep, n].toarray().ravel()
    if not sparse:
        return M.toarray(), K.toarray(), F
    return M, K, F


def assembly_trace(n_el: int, u: float, D: float, L: float = 1.0, g: float = 0.0, q: float = 0.0) -> list:
    """Element-by-element record of the assembly on a uniform mesh (the explainer's transport and inspector).

    Book: §10.3, Eqs. (10.69)–(10.78).  The global matrices here are the full (n + 1) × (n + 1) ones *before* node 0 is
    eliminated, so every scatter-add is visible.  Returns a list of dict(e (1-based), nodes (global [A(e,1), A(e,2)]), h,
    m, k (2 × 2 lists), f (element force, (10.77)), M_after, K_after (dense (n + 1)² lists after adding element e)).
    Label: analytic.
    """
    h = L / n_el
    M = np.zeros((n_el + 1, n_el + 1))
    K = np.zeros_like(M)
    out = []
    for e, A in enumerate(connectivity(n_el), start=1):
        m, k = element_matrices_linear(h, u, D)
        M[np.ix_(A, A)] += m  # scatter-add, Eq. (10.69)
        K[np.ix_(A, A)] += k
        out.append(dict(e=e, nodes=[int(A[0]), int(A[1])], h=h, m=m.tolist(), k=k.tolist(),
                        f=element_force(e, n_el, h, u, D, g, q).tolist(), M_after=M.tolist(), K_after=K.tolist()))
    return out


def solve_steady(x_nodes, u: float, D: float, g: float = 0.0, q: float | None = None, T_L: float | None = None) -> dict:
    """Steady FE solution K d = F (the ḋ = 0 case of (10.58)); nodal values T_0 … T_n (T_0 = g).

    Book: §10.3 (10.58); with ``T_L`` the right end is Dirichlet (the steady layer (10.84)–(10.85)), else Neumann slope q.
    Returns dict(x [m], T).  Example: nodes linspace(0, 1, 5), u = 1, D = 0.25, T_L = 1 → [0, 0.025, 0.1, 0.325, 1] (= the
    centred FD solution exactly — the mass matrix drops out when ∂/∂t = 0).  Label: analytic.
    """
    from scipy.sparse.linalg import spsolve

    xn = np.asarray(x_nodes, dtype=float)
    if T_L is not None:
        _, K, F = assemble_1d(xn, u, D, g, 0.0, dirichlet_right=T_L)
        d = spsolve(K.tocsc(), F)
        return dict(x=xn, T=np.concatenate([[g], np.atleast_1d(d), [T_L]]))
    _, K, F = assemble_1d(xn, u, D, g, 0.0 if q is None else q)
    d = spsolve(K.tocsc(), F)
    return dict(x=xn, T=np.concatenate([[g], np.atleast_1d(d)]))


def solve_transport(x_nodes, T0, u: float, D: float, g: float = 0.0, q: float = 0.0, dt: float = 1e-3, nsteps: int = 100,
                    theta: float = 0.5, method: str = "theta", lumped: bool = False, save_every: int | None = None) -> dict:
    """Integrate the semi-discrete system M ḋ + K d = F in time.

    Book: §10.3, Eq. (10.58) "can be integrated by numerical methods, for example Runge–Kutta, or discretized in time by
    finite-difference schemes" — both routes: ``method="theta"`` (θ-scheme, θ = ½ Crank–Nicolson, θ = 1 backward Euler; our
    choice) or ``method="solve_ivp"`` (method of lines, scipy Radau).
    Parameters: x_nodes [m]; T0 (nodal initial values T_0 … T_n, (10.3)); u, D, g, q as in :func:`assemble_1d`; dt [s];
    nsteps; theta; method; lumped; save_every.
    Returns dict(t [s], T (nodal, final, incl. T_0 = g), history, t_saved).
    Validation: V3 order 2 in space (L² at the nodes; nodal superconvergence noted) and in time for θ = ½.  Label: converged.
    """
    from scipy.sparse.linalg import splu

    M, K, F = assemble_1d(x_nodes, u, D, g, q, lumped=lumped)
    d = np.asarray(T0, dtype=float)[1:].copy()
    hist, ts = [np.concatenate([[g], d])], [0.0]
    if method == "solve_ivp":
        from scipy.integrate import solve_ivp

        Minv = splu(M.tocsc())
        t_end = dt * nsteps
        t_eval = np.linspace(0, t_end, (nsteps // (save_every or nsteps)) + 1)
        sol = solve_ivp(lambda t, y: Minv.solve(F - K @ y), (0, t_end), d, method="Radau", t_eval=t_eval,
                        rtol=1e-10, atol=1e-12, jac=(-np.asarray(splu(M.tocsc()).solve(K.toarray()))))
        hist = [np.concatenate([[g], y]) for y in sol.y.T]
        return dict(t=float(sol.t[-1]), T=hist[-1], history=np.array(hist), t_saved=sol.t)
    A = (M + theta * dt * K).tocsc()
    B = (M - (1 - theta) * dt * K).tocsr()
    lu = splu(A)
    for k in range(1, nsteps + 1):
        d = lu.solve(B @ d + dt * F)  # θ-scheme on (10.58)
        if save_every and k % save_every == 0:
            hist.append(np.concatenate([[g], d]))
            ts.append(k * dt)
    T = np.concatenate([[g], d])
    if not save_every or ts[-1] != nsteps * dt:
        hist.append(T)
        ts.append(nsteps * dt)
    return dict(t=nsteps * dt, T=T, history=np.array(hist), t_saved=np.array(ts))


def interior_stencil(h: float, u: float, D: float) -> dict:
    """The row of M and K belonging to an interior node of a uniform mesh.

    Book: §10.3, Eq. (10.63): d/dt(T_{A−1}/6 + 2T_A/3 + T_{A+1}/6) + (u/2h)(T_{A+1} − T_{A−1}) − (D/h²)(T_{A−1} − 2T_A +
    T_{A+1}) = 0, i.e. the FE row divided by h.  Returns dict(M=(h/6, 2h/3, h/6), K=(−u/2 − D/h, 2D/h, u/2 − D/h),
    M_over_h, K_over_h) for (A − 1, A, A + 1).  Validation: V2 sympy (D12).  Label: analytic.
    """
    M = np.array([h / 6.0, 2.0 * h / 3.0, h / 6.0])
    K = np.array([-0.5 * u - D / h, 2.0 * D / h, 0.5 * u - D / h])
    return dict(mass=M, stiff=K, M=M, K=K, M_over_h=M / h, K_over_h=K / h)  # Eq. (10.63) = row / h


def _composite(fn: Callable, breaks, n_quad: int):
    xg, wg = gauss_legendre(n_quad)
    tot = 0.0
    for a, b in zip(breaks[:-1], breaks[1:]):
        xx = 0.5 * (b - a) * xg + 0.5 * (b + a)
        tot += 0.5 * (b - a) * np.sum(wg * fn(xx))
    return float(tot)


def _samples(fn, x):
    return np.asarray(fn(x) if callable(fn) else fn, dtype=float)


def bilinear_form(w, v, x, u: float, D: float) -> float:
    """a(w, v) = u∫₀ᴸ v_x w dx + D∫₀ᴸ v_x w_x dx for piecewise-linear w, v on the nodes x (exact integration).

    Book: §10.3, Eq. (10.42).  w, v: callables of x [m] or node samples; both are taken as their piecewise-linear
    interpolants on ``x`` (exact for hat functions: a(N_A, N_B) = K_AB of :func:`assemble_1d`).
    Returns a [unit of w·v·m/s].  Label: analytic.
    """
    x = np.asarray(x, dtype=float)
    ws, vs = _samples(w, x), _samples(v, x)
    h = np.diff(x)
    vx = np.diff(vs) / h
    wx = np.diff(ws) / h
    wmean = 0.5 * (ws[1:] + ws[:-1])
    return float(np.sum(u * vx * wmean * h + D * vx * wx * h))  # Eq. (10.42), element by element


def _deriv(fn, x, L):
    hh = 1e-5 * L
    return (-fn(x + 2 * hh) + 8 * fn(x + hh) - 8 * fn(x - hh) + fn(x - 2 * hh)) / (12 * hh)


def weak_residual(T_fn, w_fn, u: float, D: float, q: float, L: float = 1.0, g: float = 0.0, n_pieces: int = 256,
                  n_quad: int = 6) -> float:
    """Residual of the steady weak form (10.36): u∫T_x w dx + D∫T_x w_x dx − D q w(L).

    Book: §10.3, Eq. (10.36) with ∂T/∂t = 0 (zero for every w ∈ V when T solves the steady strong problem with ∂T/∂x(L) = q;
    ``g`` only documents T(0) = g).  T_fn, w_fn: callables of x [m] (w(0) = 0); slopes by a fourth-order central difference,
    integrals by composite Gauss on ``n_pieces`` equal pieces (the kinks of hat functions sit on the pieces' ends when the
    mesh nodes are multiples of L/n_pieces).  Returns the residual [unit of T·m/s].  Label: analytic, converged.
    """
    br = np.linspace(0.0, L, n_pieces + 1)
    Tx = lambda x: _deriv(T_fn, x, L)  # noqa: E731
    wx = lambda x: _deriv(w_fn, x, L)  # noqa: E731
    val = _composite(lambda x: u * Tx(x) * w_fn(x) + D * Tx(x) * wx(x), br, n_quad)
    return float(val - D * q * float(np.atleast_1d(w_fn(np.array([L])))[0]))  # Eq. (10.36)


def galerkin_equations_sympy(n: int = 3) -> dict:
    """Galerkin equations (10.53) for n uniform linear elements, written out symbolically, and their matrix form (10.58).

    Book: §10.3, Eqs. (10.50)–(10.58) and (10.63).  Symbols h, u, D, g, q; unknowns d_1 … d_n, rates dd_1 … dd_n.
    Integrals of the hats are done element by element (on element e the only nonzero hats are N_{e−1} = (x_e − x)/h and
    N_e = (x − x_{e−1})/h).
    Returns dict(M, K, F (sympy matrices for A, B = 1 … n), equations (list: row A of M ḋ + K d − F),
    check_interior (row A = 1 minus h × (10.63) with T_0 = g: 0 when reproduced; needs n ≥ 2),
    check_M (M minus (h/6)·tridiag(1, 4, 1) with the last diagonal 2: zero matrix)).
    Validation: V2.  Label: symbolic.
    """
    import sympy as sp

    x, h, u, D, g, q = sp.symbols("x h u D g q", positive=True)

    def local(A, e):
        """Expression of N_A on element e (1-based) or 0."""
        if A == e - 1:
            return (e * h - x) / h
        if A == e:
            return (x - (e - 1) * h) / h
        return sp.Integer(0)

    def integ(A, B, da=False, db=False):
        tot = sp.Integer(0)
        for e in range(1, n + 1):
            fa, fb = local(A, e), local(B, e)
            if da:
                fa = sp.diff(fa, x)
            if db:
                fb = sp.diff(fb, x)
            tot += sp.integrate(fa * fb, (x, (e - 1) * h, e * h))
        return sp.simplify(tot)

    idx = range(1, n + 1)
    M = sp.Matrix([[integ(A, B) for B in idx] for A in idx])  # Eq. (10.55)
    K = sp.Matrix([[u * integ(A, B, db=True) + D * integ(A, B, da=True, db=True) for B in idx] for A in idx])  # Eq. (10.56)
    F = sp.Matrix([D * q * (1 if A == n else 0) - g * u * integ(A, 0, db=True) - g * D * integ(A, 0, da=True, db=True)
                   for A in idx])  # Eq. (10.57)
    d = sp.Matrix(sp.symbols(f"d1:{n + 1}"))
    dd = sp.Matrix(sp.symbols(f"dd1:{n + 1}"))
    eqs = list(M * dd + K * d - F)  # Eq. (10.58): M ḋ + K d = F
    check = None
    if n >= 2:
        T0, T1, T2 = g, d[0], d[1]
        r63 = (2 * dd[0] / 3 + dd[1] / 6) + u / (2 * h) * (T2 - T0) - D / h ** 2 * (T0 - 2 * T1 + T2)  # Eq. (10.63), dT0/dt = 0
        check = sp.simplify(sp.expand(eqs[0] - h * r63.subs(q, 0)))
    Mc = sp.Matrix(n, n, lambda i, j: 4 if i == j else (1 if abs(i - j) == 1 else 0)) * h / 6
    Mc[n - 1, n - 1] = 2 * h / 6
    return dict(G_A=eqs, M=M, K=K, F=F, equations=eqs, interior_row=eqs[0], check_interior=check,
                check_M=sp.simplify(M - Mc),
                symbols=dict(h=h, u=u, D=D, g=g, q=q))


