"""Airfoil sections in ideal flow: Zhukhovsky shapes from off-centre circles, the flow round them at incidence, the
Kutta condition and its circulation, surface speed and pressure, section forces, and thin-airfoil coefficients.

Book: Kundu, Cohen & Dowling 5e, Ch. 14 §14.3–§14.5, Eqs. (14.3)–(14.12) (read from the rendered pages
chapters/pages/ch14/p729–p734) and Exercises 14.2–14.9 (results only; exercise text and numbers are never reproduced).
Built on ``core.conformal`` (the map (14.3) = (6.65), ``circle_flow_zeta(center=)``) and ``core.potential``
(``blasius_force``).

Conventions (``knowledge/notation.md``)
---------------------------------------
* ζ-plane: the circle; z-plane: the section, z = ζ + b²/ζ (14.3).  **x is the mapped chord axis**; the stream comes in
  at the angle of attack α to it (Fig. 14.16), so **lift is not F_y** — forces are returned as (Fx, Fy) *and* (L, D).
* dw/dz = u − iv (the ch06 convention; ch13's Ekman velocity was u + iv).
* **Γ is the clockwise circulation** (the book's sign in (6.62), (14.10), (14.11)); positive lift at positive α needs
  Γ_cw > 0.  Arguments are keyword-only ``Gamma_cw=``; returned values are named ``Gamma_cw``.
* **The circle** (ours, reduces to the book's two cases): centre μ = b − a e^{−iβ}, radius a = b(1 + e)/cos β, so that
  it always passes through ζ = b (the trailing edge).  e = thickness parameter of (14.7) (β = 0: a = b(1 + e), centre
  at −be, Fig. 14.13); β = camber angle of Fig. 14.12 (e = 0: a = b/cos β, centre at i b tan β).  The contour is
  ζ(t) = μ + a e^{it} with t the polar angle about the centre Q from the ξ direction, t = −β + 2πk/n:
  **index 0 (t = −β) is the trailing edge**, traversal counterclockwise (upper surface first).
* Three different θ (trap T4): the polar angle of ζ from the origin in (14.4); the angle from the stream-aligned
  diameter through the centre Q in (14.10) (called ``theta`` in :func:`rear_stagnation_point`); neither is ``t``.
* Angles in radians.  Scalars in, floats out.  ``e`` here is never the span efficiency (that is ``e_span``).
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.optimize import brentq

from ._util import as_scalar_if_0d
from .conformal import circle_flow_zeta, joukowski
from .potential import blasius_force

__all__ = [
    "THICKNESS_RULE_COEFFICIENT", "zhukhovsky_circle", "zhukhovsky_circle_from_centre", "zhukhovsky_airfoil",
    "zhukhovsky_inverse", "section_geometry", "trailing_edge_angle", "zhukhovsky_shape", "circular_arc",
    "zero_lift_angle", "camber_from_beta", "resolve_lift_drag", "section_coefficient", "circulation_from_cl",
    "cl_from_circulation", "normal_force_coefficient", "lift_from_cp", "kutta_circulation", "kutta_circulation_root",
    "stagnation_points", "rear_stagnation_point", "trailing_edge_speed", "ZhukhovskyFlow", "zhukhovsky_flow",
    "flat_plate_flow", "zhukhovsky_lift_coefficient", "zhukhovsky_surface", "surface_circulation",
    "zhukhovsky_force", "flat_plate_forces", "thin_airfoil",
]

_S = as_scalar_if_0d
_TWO_PI = 2.0 * np.pi

#: first-order coefficient of the thickness rule t_max/c ≃ (3√3/4) e of (14.8); the book prints it rounded, "1.3".
THICKNESS_RULE_COEFFICIENT = 3.0 * np.sqrt(3.0) / 4.0


def _F(a):
    return np.asarray(a, dtype=float)


def _C(a):
    return np.asarray(a, dtype=complex)


# ======================================================================================================================
# shapes (§14.4)
# ======================================================================================================================
def zhukhovsky_circle(b: float = 1.0, e: float = 0.0, beta: float = 0.0, *, a: float | None = None) -> tuple:
    """Centre μ and radius a of the ζ-plane circle that the map (14.3) turns into a Zhukhovsky section.

    Book: §14.4, Eq. (14.7) a ≡ b(1 + e) (symmetric case, Fig. 14.13), "QP = a = b/cos β" (circular arc, Fig. 14.12)
    and Fig. 14.14 (both displacements).  Ours (trap T7): μ = b − a e^{−iβ}, a = b(1 + e)/cos β — the one circle
    through ζ = b that reduces to both printed cases.
    Parameters
    ----------
    b : map constant [m] (> 0; chord ≃ 4b).   e : thickness parameter [–] (≥ 0).   beta : camber angle β [rad]
    (|β| < π/2).   a : optional radius [m]; if given, (a, b, β) define the circle and ``e`` is ignored
    (e = a cos β/b − 1) — the parametrisation of the exercise that plots a cambered section.
    Returns
    -------
    (mu, a) : complex centre [m], radius [m].
    Assumptions: none (geometry).
    Validation: V1 e = 0 → μ = i b tan β; β = 0 → μ = −be; |b − μ| = a; V7 round trip with
    :func:`zhukhovsky_circle_from_centre` — tests/test_ch14.py: test_zhukhovsky_circle_V1_reductions_and_round_trip,
    test_circle_flow_zeta_V1_offcentre_circle_is_a_streamline_with_the_right_far_field_and_circulation,
    test_kutta_condition_V1_root_finder_equals_14_11_and_the_edge_point_stagnates.
    Label: analytic.
    """
    b, beta = float(b), float(beta)
    if a is None:
        a = b * (1.0 + float(e)) / np.cos(beta)  # Eq. (14.7) with β = 0; a = b/cos β with e = 0
    a = float(a)
    mu = b - a * np.exp(-1j * beta)
    return complex(mu), a


def zhukhovsky_circle_from_centre(mu: complex, b: float = 1.0) -> tuple:
    """Inverse of :func:`zhukhovsky_circle`: (e, β, a) of the circle centred at μ that passes through ζ = b.

    Book: §14.4, Figs. 14.12–14.14 (ours: the inverse of our parametrisation, for a draggable centre).
    Parameters: mu complex centre [m]; b map constant [m].   Returns (e [–], beta [rad], a [m]).
    Assumptions: Re μ < b (the centre is not beyond the trailing-edge point).
    Validation: V7 round trip — tests/test_ch14.py: test_zhukhovsky_circle_V1_reductions_and_round_trip.
    Label: analytic.
    """
    d = float(b) - complex(mu)  # = a e^{−iβ}
    a = abs(d)
    beta = -float(np.angle(d))
    return float(a * np.cos(beta) / float(b) - 1.0), beta, float(a)


def zhukhovsky_airfoil(b: float = 1.0, e: float = 0.0, beta: float = 0.0, n: int = 400, *,
                       a: float | None = None) -> dict:
    """The Zhukhovsky section: image of the off-centre circle under z = ζ + b²/ζ (14.3).

    Book: §14.4, Eq. (14.3), the unnumbered image z = R e^{iθ} + (b²/R) e^{−iθ} and Eq. (14.4)
    x = (R + b²/R) cos θ, y = (R − b²/R) sin θ, Figs. 14.11–14.14.  e = β = 0 gives the slit −2b ≤ x ≤ 2b (flat plate),
    e = 0 the circular arc, β = 0 the symmetric section, both the cambered airfoil.
    Parameters
    ----------
    b [m]; e [–]; beta [rad]; n : number of contour points (t uniform on [−β, −β + 2π), endpoint excluded);
    a : optional radius [m] (see :func:`zhukhovsky_circle`).
    Returns
    -------
    dict: ``t`` [rad] (polar angle about the centre Q; ζ = μ + a e^{it}), ``zeta``, ``z`` (complex [m]), ``x``, ``y`` [m], ``R`` = |ζ| [m], ``theta`` = arg ζ [rad]
    (the R and θ of (14.4)), ``mu`` (complex), ``a`` [m], ``b``, ``e``, ``beta``.  Index 0 = trailing edge (z = 2b);
    counterclockwise, upper surface first.
    Assumptions: none (geometry); the map is one-to-one outside the circle when ζ = −b is not outside it.
    Validation: V1 slit for e = β = 0; every point of the e = 0 image on the circle of :func:`circular_arc`; V2 (14.4)
    reproduces x, y — tests/test_ch14.py: test_zhukhovsky_airfoil_V1_slit_arc_symmetric_section_and_14_4,
    test_trailing_edge_V7_cusp_angle_tends_to_zero_and_a_wedge_keeps_its_angle, test_cusp_V2_derivation.
    Label: analytic.
    """
    mu, a = zhukhovsky_circle(b, e, beta, a=a)
    b, beta = float(b), float(beta)
    t = -beta + np.linspace(0.0, _TWO_PI, int(n), endpoint=False)
    zeta = mu + a * np.exp(1j * t)
    R, th = np.abs(zeta), np.angle(zeta)
    x = (R + b ** 2 / R) * np.cos(th)  # Eq. (14.4)
    y = (R - b ** 2 / R) * np.sin(th)  # Eq. (14.4)
    z = x + 1j * y  # = joukowski(zeta, b), Eq. (14.3)
    e_out = a * np.cos(beta) / b - 1.0
    return {"t": t, "zeta": zeta, "z": z, "x": x, "y": y, "R": R, "theta": th, "mu": mu, "a": a, "b": b,
            "e": float(e_out), "beta": beta}


def zhukhovsky_inverse(z, b: float = 1.0, mu: complex = 0j, a: float | None = None):
    """Inverse of (14.3) for an off-centre circle: the root ζ of ζ² − zζ + b² = 0 that lies outside |ζ − μ| = a.

    Book: §14.5 (the map is used in both directions; (6.69) is its centred form).  Trap T8: for an off-centre circle
    "outside" means outside the *shifted* circle, so the root is chosen by the larger |ζ − μ| (the two roots are ζ and
    b²/ζ and exactly one of them is outside when z is outside the section) — not by |ζ| ≥ b.
    Parameters: z complex [m] (any shape); b [m]; mu complex centre [m]; a radius [m] or None (no body test).
    Returns ζ (complex); NaN where both roots are inside the circle (z inside the section), if ``a`` is given.
    Assumptions: ζ = −b is not outside the circle (the map is then one-to-one outside it); the square root is cut on
    the slit −2b ≤ x ≤ 2b, so ζ(z) is continuous everywhere outside the section.
    Validation: V1 joukowski(ζ) = z round trip on a grid outside the section; μ = 0 equals
    ``core.conformal.joukowski_inverse(branch="outside")`` — tests/test_ch14.py:
    test_zhukhovsky_inverse_V1_round_trip_and_the_root_outside_the_shifted_circle.
    Label: analytic.
    """
    zz = _C(z)
    b = float(b)
    root = np.sqrt(zz - 2.0 * b) * np.sqrt(zz + 2.0 * b)  # cut on the slit (ch06 lesson: not sqrt(z² − 4b²))
    z1 = 0.5 * (zz + root)
    z2 = 0.5 * (zz - root)
    mu = complex(mu)
    d1, d2 = np.abs(z1 - mu), np.abs(z2 - mu)
    zeta = np.where(d1 >= d2, z1, z2)
    if a is not None:
        zeta = np.where(np.maximum(d1, d2) < float(a) * (1.0 - 1e-12), np.nan + 0j, zeta)
    return _S(zeta)


def _parabolic_peak(xm, x0, xp, fm, f0, fp):
    """Vertex (x*, f*) of the parabola through three samples (falls back to the middle sample)."""
    d1, d2 = (f0 - fm) / (x0 - xm), (fp - f0) / (xp - x0)
    curv = (d2 - d1) / (xp - xm)
    if not np.isfinite(curv) or abs(curv) < 1e-300:
        return x0, f0
    xs = 0.5 * (xm + x0) - 0.5 * d1 / curv
    if not (min(xm, xp) <= xs <= max(xm, xp)):
        return x0, f0
    return xs, f0 + d1 * (xs - x0) + curv * (xs - x0) * (xs - xm)


def section_geometry(x, y, te_index: int | None = None, m: int = 0) -> dict:
    """Chord, camber line, camber and thickness of a closed section contour.

    Book: §14.3, Fig. 14.6 (definitions: chord line joining the leading and trailing edges, camber line midway between
    the two surfaces, camber = its largest height above the chord, thickness t).
    Method (ours): trailing edge = ``te_index`` (default: the point of largest x); leading edge = the contour point
    farthest from it (refined by a parabola through three samples); in the chord frame (ξ from the leading edge along
    the chord, η normal to it) both surfaces are interpolated linearly on common cosine-spaced stations; camber line =
    mean, thickness = difference; the two maxima are refined by a parabola.  Second-order accurate in the point spacing
    (linear interpolation of the surfaces), with a larger constant for a sharp nose.  Second order only when the
    trailing edge is a contour vertex (true for every Zhukhovsky contour, index 0); for a rounded trailing edge that
    falls between two vertices the camber is first order — pass ``te_index`` or refine.
    Parameters
    ----------
    x, y : contour points [m] (closed, either orientation, last point ≠ first).   te_index : index of the trailing
    edge.   m : number of chord stations (0 → max(len(x), 400)).
    Returns
    -------
    dict: ``chord`` [m]; ``camber`` [m] (signed: + toward the upper surface); ``x_camber`` (ξ/c of the maximum);
    ``t_max`` [m]; ``x_tmax`` (ξ/c of the maximum thickness); ``camber_ratio``, ``thickness_ratio`` [–];
    ``camber_line`` = (ξ/c, η_c [m]); ``xc`` (stations), ``y_upper``, ``y_lower`` [m] in the chord frame; ``le``,
    ``te`` (complex [m]); ``chord_angle`` [rad] (direction leading → trailing edge, from +x).
    Assumptions: a section with one leading and one trailing edge, each surface single-valued in ξ.
    Validation: V1 flat plate, circular arc (camber 2b tan β) and symmetric Zhukhovsky sections; V3 second order in the
    number of points — tests/test_ch14.py: test_section_geometry_V1_plate_arc_and_symmetric_sections,
    test_section_geometry_V3_second_order_in_the_number_of_points,
    test_zhukhovsky_lift_V3_pressure_route_converges_to_the_blasius_force_at_second_order.
    Label: analytic, converged.
    """
    z = _F(x) + 1j * _F(y)
    n = z.size
    it = int(np.argmax(z.real)) if te_index is None else int(te_index) % n
    z = np.roll(z, -it)  # trailing edge first
    te = z[0]
    dist = np.abs(z - te)
    il = int(np.argmax(dist))
    le = z[il]
    f = 0.0
    if 0 < il < n - 1:  # refine the farthest point with a parabola in the index
        s_star, _ = _parabolic_peak(il - 1.0, float(il), il + 1.0, dist[il - 1], dist[il], dist[il + 1])
        f = s_star - il
        # the edge itself: the parabola through the three vertices, evaluated at the refined index (third order)
        le = 0.5 * f * (f - 1.0) * z[il - 1] + (1.0 - f * f) * z[il] + 0.5 * f * (f + 1.0) * z[il + 1]
    chord = abs(te - le)
    if chord == 0.0:
        raise ValueError("degenerate contour: leading and trailing edges coincide")
    rot = np.conj((te - le) / chord)
    q = (z - le) * rot  # chord frame: ξ = Re q from the leading edge, η = Im q
    k1 = il + 1 if f >= 0 else il  # the vertex il lies on the first surface when the leading edge is beyond it
    s1 = np.concatenate([q[:k1], [0j]])  # trailing edge → leading edge (the refined edge itself is ξ = η = 0)
    s2 = np.concatenate([[0j], q[k1:], q[:1]])  # leading edge → trailing edge
    m = int(m) if m else max(n, 400)
    xi = 0.5 * chord * (1.0 - np.cos(np.linspace(0.0, np.pi, m)))

    def on_stations(sq):
        xr = np.where(np.abs(sq.real) < 1e-12 * chord, 0.0, sq.real)  # vertices level with the refined edge: ξ = 0
        o = np.lexsort((np.abs(sq.imag), xr))  # by ξ; ties by |η| (the edge itself first)
        v = np.interp(xi, xr[o], sq.imag[o])
        v[0] = 0.0  # both surfaces meet at the leading edge, ξ = η = 0
        return v

    e1, e2 = on_stations(s1), on_stations(s2)
    if np.mean(e1) >= np.mean(e2):
        yu, yl = e1, e2
    else:
        yu, yl = e2, e1
    yc, th = 0.5 * (yu + yl), yu - yl

    def peak(f):
        k = int(np.argmax(np.abs(f)))
        if 0 < k < m - 1:
            xs, fs = _parabolic_peak(xi[k - 1], xi[k], xi[k + 1], f[k - 1], f[k], f[k + 1])
            return float(xs), float(fs)
        return float(xi[k]), float(f[k])

    xcam, cam = peak(yc)
    xt, tm = peak(th)
    return {"chord": float(chord), "camber": cam, "x_camber": float(xcam / chord), "t_max": tm,
            "x_tmax": float(xt / chord), "camber_ratio": float(cam / chord), "thickness_ratio": float(tm / chord), "camber_line": (xi / chord, yc),
            "xc": xi / chord, "y_upper": yu, "y_lower": yl, "le": complex(le), "te": complex(te),
            "chord_angle": float(np.angle(te - le))}


def trailing_edge_angle(x, y, te_index: int | None = None) -> float:
    """Included angle of a section at its trailing edge, from the two contour points next to it [rad].

    Book: §14.4, Fig. 14.15 (finite angle vs cusp).  For a Zhukhovsky section the angle tends to 0 as the contour is
    refined (a cusp: dz/dζ = 1 − b²/ζ² vanishes at ζ = b, so the angle π of the circle doubles to 2π outside).
    Parameters: x, y contour [m]; te_index (default: the point of largest x).   Returns the angle [rad] in [0, π].
    Assumptions: the two neighbours are close enough to stand for the tangents (first order in the spacing).
    Validation: V7 → 0 for a Zhukhovsky contour as n grows; equals the wedge angle for a polygon — tests/test_ch14.py:
    test_trailing_edge_V7_cusp_angle_tends_to_zero_and_a_wedge_keeps_its_angle.
    Label: analytic.
    """
    z = _F(x) + 1j * _F(y)
    n = z.size
    it = int(np.argmax(z.real)) if te_index is None else int(te_index) % n
    v1, v2 = z[(it + 1) % n] - z[it], z[(it - 1) % n] - z[it]
    return float(abs(np.angle(v1 / v2)))


def zhukhovsky_shape(b: float = 1.0, e: float = 0.0, beta: float = 0.0, exact: bool = True, *,
                     a: float | None = None, n: int = 20001) -> dict:
    """Chord, camber and maximum thickness of a Zhukhovsky section: measured on the contour, or the rules of (14.8).

    Book: §14.4, Eq. (14.8) c ≃ 4b, camber ≃ ½βc, t_max/c ≃ 1.3 e ("can be proved for e ≪ 1"), maximum thickness
    "nearly at the quarter chord position x = −b".  Slip #4 (loose): the first-order coefficient is 3√3/4 = 1.299
    (:data:`THICKNESS_RULE_COEFFICIENT`, which is what ``exact=False`` uses; the book rounds it to 1.3), and the rule is
    first order in e — the measured ratio at e = 0.1 is several per cent lower.
    Parameters
    ----------
    b [m]; e [–]; beta [rad]; exact : True → measured with :func:`section_geometry` on an ``n``-point contour;
    False → the three rules.   a : optional radius [m] (then e is derived).
    Returns
    -------
    dict: ``chord`` [m], ``camber`` [m], ``t_max`` [m], ``x_tmax`` [m] (z-plane x of the maximum thickness, ≈ −b),
    ``camber_ratio``, ``thickness_ratio`` [–], ``x_tmax_over_c``, ``x_camber_over_c`` [–] (from the leading edge),
    ``e``, ``beta``, ``a`` [m], ``exact`` (bool).
    Assumptions (rules only): e ≪ 1, β ≪ 1.
    Validation: V1 e = 0: camber = 2b tan β, chord 4b exactly; V7 thickness_ratio/e → 3√3/4 and x_tmax → −b as e → 0;
    planted variant a = b(1 + e) cos β fails — tests/test_ch14.py:
    test_zhukhovsky_shape_V7_rules_of_14_8_are_first_order_and_slip_4, test_thickness_rule_V2_derivation,
    test_section_geometry_V1_plate_arc_and_symmetric_sections.
    Label: analytic.
    """
    mu, a_ = zhukhovsky_circle(b, e, beta, a=a)
    b, beta = float(b), float(beta)
    e_ = a_ * np.cos(beta) / b - 1.0
    if not exact:
        c = 4.0 * b  # Eq. (14.8)
        cam = 0.5 * beta * c  # Eq. (14.8)
        tm = THICKNESS_RULE_COEFFICIENT * e_ * c  # Eq. (14.8) (printed 1.3)
        return {"chord": c, "camber": cam, "t_max": float(tm), "x_tmax": -b, "camber_ratio": 0.5 * beta,
                "thickness_ratio": float(THICKNESS_RULE_COEFFICIENT * e_), "x_tmax_over_c": 0.25, "x_camber_over_c": 0.5,
                "e": float(e_), "beta": beta, "a": a_, "exact": False}
    A = zhukhovsky_airfoil(b, e_, beta, int(n))
    g = section_geometry(A["x"], A["y"], te_index=0)
    xt = (g["le"] + g["x_tmax"] * g["chord"] * np.exp(1j * g["chord_angle"])).real
    return {"chord": g["chord"], "camber": g["camber"], "t_max": g["t_max"], "x_tmax": float(xt),
            "camber_ratio": g["camber_ratio"], "thickness_ratio": g["thickness_ratio"], "x_tmax_over_c": g["x_tmax"],
            "x_camber_over_c": g["x_camber"], "e": float(e_), "beta": beta, "a": a_, "exact": True}


def circular_arc(b: float = 1.0, beta: float = 0.1) -> dict:
    """The circular arc that the raised circle (e = 0) maps to: its circle, height and end angle.

    Book: §14.4, the unnumbered result x² + (y + 2b cot 2β)² = (2b csc 2β)² with cos β = b/a, from (14.5) and (14.6);
    Fig. 14.12 (end tangent at the angle 2β to the chord; Exercise 14.6).
    Parameters: b [m]; beta [rad] (0 < β < π/4 for an arc above the chord with its centre below it).
    Returns dict: ``centre`` (complex, −2ib cot 2β) [m], ``radius`` = 2b csc 2β [m], ``height`` = 2b tan β [m],
    ``end_angle`` = 2β [rad], ``chord`` = 4b [m], ``a`` = b/cos β [m].
    Assumptions: none (geometry).
    Validation: V1 every mapped point satisfies the circle equation (residual 1e-12); V2 sympy
    (``ch14.circular_arc_sym``) — tests/test_ch14.py: test_zhukhovsky_airfoil_V1_slit_arc_symmetric_section_and_14_4,
    test_circular_arc_V2_derivation.
    Label: analytic, symbolic.
    """
    b, beta = float(b), float(beta)
    return {"centre": complex(0.0, -2.0 * b / np.tan(2.0 * beta)), "radius": float(2.0 * b / np.sin(2.0 * beta)),
            "height": float(2.0 * b * np.tan(beta)), "end_angle": 2.0 * beta, "chord": 4.0 * b,
            "a": float(b / np.cos(beta))}


def zero_lift_angle(camber_ratio):
    """Angle of attack of zero lift of a thin Zhukhovsky section, α_{L=0} = −β = −2 camber/chord [rad].

    Book: §14.5 (text after (14.12): "the lift is zero at a negative angle of attack α = −β"; "camber/chord ratio
    β/2"), with camber ≃ ½βc of (14.8).   Parameters: camber_ratio = camber/chord [–].   Returns α_{L=0} [rad].
    Assumptions: thin, slightly cambered section (e ≪ 1, β ≪ 1).
    Validation: V1 against ``thin_airfoil`` for a parabolic arc (−2h/c) — tests/test_ch14.py:
    test_zhukhovsky_lift_V7_small_angle_rule_and_slip_8_thickness_raises_the_slope,
    test_thin_airfoil_V1_plate_parabolic_arc_and_exact_coefficients_of_a_cubic_camber_line,
    test_zhukhovsky_shape_V7_rules_of_14_8_are_first_order_and_slip_4.
    Label: analytic.
    """
    return _S(-2.0 * _F(camber_ratio))


def camber_from_beta(beta, c):
    """Camber of a thin Zhukhovsky section from its camber angle, camber ≃ ½βc (14.8) [m].

    Book: §14.4, Eq. (14.8).   Parameters: beta [rad]; c chord [m].   Assumptions: β ≪ 1 (exactly: 2b tan β).
    Validation: V7 → 2b tan β as β → 0 — tests/test_ch14.py:
    test_zhukhovsky_shape_V7_rules_of_14_8_are_first_order_and_slip_4.
    Label: analytic.
    """
    return _S(0.5 * _F(beta) * _F(c))  # Eq. (14.8)


# ======================================================================================================================
# forces and coefficients (§14.3)
# ======================================================================================================================
def resolve_lift_drag(Fx, Fy, alpha):
    """Lift (normal to the stream) and drag (along it) from the force components in section axes.

    Book: §14.3, Fig. 14.7 (definitions of L and D); ours for a stream at the angle α to x:
    D = Fx cos α + Fy sin α, L = −Fx sin α + Fy cos α.
    Parameters: Fx, Fy [N or N/m]; alpha [rad] (stream direction measured from +x).   Returns (L, D).
    Assumptions: none.
    Validation: V1 L² + D² = Fx² + Fy²; α = 0 identity — tests/test_ch14.py:
    test_resolve_lift_drag_V1_rotation_and_section_coefficients, test_pressure_loop_V2_derivation.
    Label: analytic.
    """
    Fx, Fy, al = _F(Fx), _F(Fy), _F(alpha)
    D = Fx * np.cos(al) + Fy * np.sin(al)
    L = -Fx * np.sin(al) + Fy * np.cos(al)
    return _S(L), _S(D)


def section_coefficient(F_per_span, rho, U, c):
    """Section force coefficient F′/(½ρU²c): (4.107)–(4.108) with the planform area replaced by chord × unit span.

    Book: §14.3 (restated C_D ≡ D/(½ρU²A), C_L ≡ L/(½ρU²A), A = planform area; a section uses c per unit span).
    Parameters: F_per_span [N/m]; rho [kg/m³]; U [m/s]; c [m].   Returns the coefficient [–].
    Assumptions: none.
    Validation: V2 dimensionless; V1 round trip with :func:`circulation_from_cl` — tests/test_ch14.py:
    test_resolve_lift_drag_V1_rotation_and_section_coefficients.
    Label: analytic.
    """
    return _S(_F(F_per_span) / (0.5 * _F(rho) * _F(U) ** 2 * _F(c)))


def circulation_from_cl(CL, U, c):
    """Clockwise circulation of a section from its lift coefficient, Γ = ½ U c C_L [m²/s].

    Book: §14.6, unnumbered line before (14.17): C_L = L/(½ρU²c) = Γ/(½Uc), from L = ρUΓ (6.62).
    Parameters: CL [–]; U [m/s]; c [m].   Returns Γ_cw [m²/s].
    Assumptions: steady ideal two-dimensional flow.
    Validation: V1 round trip — tests/test_ch14.py: test_resolve_lift_drag_V1_rotation_and_section_coefficients,
    test_lifting_line_V1_elliptic_planform_untwisted_gives_the_elliptic_loading_exactly.
    Label: analytic.
    """
    return _S(0.5 * _F(U) * _F(c) * _F(CL))


def cl_from_circulation(Gamma_cw, U, c):
    """Section lift coefficient from the clockwise circulation, C_L = Γ/(½Uc) [–].

    Book: §14.6, unnumbered line before (14.17).   Parameters: Gamma_cw [m²/s]; U [m/s]; c [m].
    Assumptions: steady ideal two-dimensional flow.
    Validation: V1 round trip — tests/test_ch14.py: test_resolve_lift_drag_V1_rotation_and_section_coefficients.
    Label: analytic.
    """
    return _S(_F(Gamma_cw) / (0.5 * _F(U) * _F(c)))


def normal_force_coefficient(xc_upper, cp_upper, xc_lower, cp_lower) -> float:
    """Coefficient of the force normal to the chord axis as the area of the C_p loop, C_y = ∮ C_p d(x/c).

    Book: §14.3, Fig. 14.8 and Exercises 14.3–14.4: F_y = −∫₀ᶜ p_u dx + ∫₀ᶜ p_l dx, so
    C_y = F_y/(½ρU²c) = ∫₀¹ (C_p,lower − C_p,upper) d(x/c).
    Parameters: xc_* chordwise positions x/c of the samples on each surface (increasing); cp_* the pressure
    coefficients there [–].   Returns C_y [–] (trapezoid rule on each surface).
    Assumptions: shear stress neglected; x along the chord axis.
    Validation: V1 exact for piecewise-linear C_p; equals F_y of :func:`zhukhovsky_force` on a Zhukhovsky section —
    tests/test_ch14.py: test_surface_pressure_V1_stagnation_points_bernoulli_and_the_edge_value,
    test_surface_pressure_V4_loop_area_is_the_normal_force_and_the_speed_integral_is_the_circulation.
    Label: analytic.
    """
    lower = np.trapezoid(_F(cp_lower), _F(xc_lower))
    upper = np.trapezoid(_F(cp_upper), _F(xc_upper))
    return float(lower - upper)


def lift_from_cp(xc_upper, cp_upper, xc_lower, cp_lower, alpha: float, suction: str = "ideal") -> float:
    """Lift coefficient from the two surface-pressure curves and the angle of attack.

    Book: §14.3, Exercise 14.3 (C_y and its relation to C_L through α).  The pressure loop gives only the force normal
    to the chord axis; the lift is normal to the stream.  ``suction="ideal"`` (default): the flow is ideal, so the
    resultant is normal to the stream (no drag) and C_L = C_y/cos α — the chordwise (leading-edge suction) force is
    implied.  ``suction="none"``: only the normal force is projected, C_L = C_y cos α (what a thin plate without
    leading-edge suction would give).  The two differ by O(α²).
    Parameters: as :func:`normal_force_coefficient`; alpha [rad].   Returns C_L [–] (based on the same length as x/c).
    Assumptions: as stated per option.
    Validation: V1 against 2Γ/(Uc) on a Zhukhovsky section — tests/test_ch14.py:
    test_surface_pressure_V4_loop_area_is_the_normal_force_and_the_speed_integral_is_the_circulation.
    Label: analytic.
    """
    cy = normal_force_coefficient(xc_upper, cp_upper, xc_lower, cp_lower)
    if suction == "ideal":
        return float(cy / np.cos(alpha))
    if suction == "none":
        return float(cy * np.cos(alpha))
    raise ValueError('suction must be "ideal" or "none"')


# ======================================================================================================================
# Kutta condition and circulation (§14.3, §14.5)
# ======================================================================================================================
def kutta_circulation(U, a, alpha, beta):
    """Clockwise circulation required by the Kutta condition, Γ = 4πUa sin(α + β) (14.11) [m²/s].

    Book: §14.5, Eq. (14.11), from u_θ = −2U sin θ − Γ/(2πa) = 0 (14.10) at B, where θ = −(α + β) (Fig. 14.16).
    Parameters
    ----------
    U : stream speed [m/s].   a : radius of the ζ-plane circle [m] (not b).   alpha : angle of attack [rad].
    beta : camber angle [rad].
    Returns
    -------
    Gamma_cw [m²/s] (clockwise; positive for α + β > 0).
    Assumptions: steady ideal flow; rear stagnation point of the circle at ζ = b (Kutta condition); no stall.
    Validation: V1 equals :func:`kutta_circulation_root`; dw/dζ(b) = 0 with it; V7 zero at α = −β; planted variants
    sin(α − β), radius b, counterclockwise sign fail — tests/test_ch14.py:
    test_kutta_condition_V1_root_finder_equals_14_11_and_the_edge_point_stagnates,
    test_kutta_condition_V5_published_form_of_the_circulation,
    test_stagnation_points_V1_both_points_for_any_circulation.
    Label: analytic.
    """
    return _S(4.0 * np.pi * _F(U) * _F(a) * np.sin(_F(alpha) + _F(beta)))  # Eq. (14.11)


def kutta_circulation_root(U: float, b: float, e: float, beta: float, alpha: float, xtol: float = 1e-14) -> float:
    """The Kutta circulation found numerically: the Γ_cw for which the circle flow has no slip speed at ζ = b.

    Book: §14.3 (Kutta condition) and §14.5 ("arrange that B is a stagnation point in the ζ-plane at which
    dw/dζ = 0").  Independent of (14.11): ``brentq`` on the tangential speed at ζ = b computed from
    ``core.conformal.circle_flow_zeta(center=μ)``, bracket |Γ| ≤ 4πUa·1.01 (beyond 4πUa the stagnation points leave
    the surface).
    Parameters: U [m/s]; b [m]; e [–]; beta [rad]; alpha [rad]; xtol absolute tolerance on Γ [m²/s].
    Returns Gamma_cw [m²/s].
    Assumptions: as :func:`kutta_circulation`.
    Validation: V1 equals (14.11) to 1e-12 — tests/test_ch14.py:
    test_kutta_condition_V1_root_finder_equals_14_11_and_the_edge_point_stagnates,
    test_kutta_condition_V5_published_form_of_the_circulation.
    Label: analytic.
    """
    mu, a = zhukhovsky_circle(b, e, beta)
    tangent = 1j * np.exp(-1j * float(beta))  # counterclockwise unit tangent of the circle at ζ = b

    def slip(G):
        _, dw = circle_flow_zeta(U, a, Gamma_cw=G, alpha=alpha, center=mu)
        return float((complex(dw(float(b))) * tangent).real)  # u_θ at B

    G0 = 4.0 * np.pi * float(U) * a * 1.01
    return float(brentq(slip, -G0, G0, xtol=xtol, rtol=8.9e-16))


def stagnation_points(U: float, b: float, e: float, beta: float, alpha: float, *,
                      Gamma_cw: float | None = None) -> dict:
    """Both stagnation points of the circle flow (and their images on the section) for any circulation.

    Book: §14.5, Eq. (14.10) u_θ = −2U sin θ − Γ/(2πa) = 0 ⇒ sin θ = −Γ/(4πUa), θ from the stream-aligned diameter
    through the centre Q (as (6.38)); Fig. 14.16 (points S and B and their images S′ and B′).
    Parameters: U [m/s]; b [m]; e [–]; beta, alpha [rad]; Gamma_cw [m²/s] (None → the Kutta value).
    Returns dict: ``theta_rear``, ``theta_front`` [rad] (from the stream-aligned diameter; NaN when off the surface);
    ``zeta_rear``, ``zeta_front`` (complex, ζ-plane); ``z_rear``, ``z_front`` (complex, z-plane); ``on_surface``
    (False when |Γ| > 4πUa: the single stagnation point in the fluid is returned as ``zeta_rear`` and the inside root
    as ``zeta_front``); ``Gamma_cw``, ``Gamma_kutta``; ``rear_at_edge`` (bool).
    Assumptions: steady ideal flow.
    Validation: V1 dw/dζ = 0 at both points; with the Kutta value the rear point is ζ = b — tests/test_ch14.py:
    test_stagnation_points_V1_both_points_for_any_circulation,
    test_surface_pressure_V1_stagnation_points_bernoulli_and_the_edge_value.
    Label: analytic.
    """
    mu, a = zhukhovsky_circle(b, e, beta)
    GK = float(kutta_circulation(U, a, alpha, beta))
    G = GK if Gamma_cw is None else float(Gamma_cw)
    k = G / (4.0 * np.pi * float(U) * a)
    # roots of U(1 − a²/σ²) + iΓ/(2πσ) = 0 with σ = (ζ − μ)e^{−iα}:  σ = a[−ik ± sqrt(1 − k²)]
    root = np.sqrt(complex(1.0 - k * k))
    s_rear, s_front = a * (-1j * k + root), a * (-1j * k - root)
    on = abs(k) <= 1.0
    if not on and abs(s_rear) < abs(s_front):
        s_rear, s_front = s_front, s_rear
    rot = np.exp(1j * float(alpha))
    zr, zf = mu + s_rear * rot, mu + s_front * rot
    th_r = float(-np.arcsin(k)) if on else float("nan")  # Eq. (14.10) = 0
    th_f = float(np.pi + np.arcsin(k)) if on else float("nan")
    return {"theta_rear": th_r, "theta_front": th_f, "zeta_rear": complex(zr), "zeta_front": complex(zf),
            "z_rear": complex(joukowski(zr, b)), "z_front": complex(joukowski(zf, b)), "on_surface": bool(on),
            "Gamma_cw": G, "Gamma_kutta": GK, "rear_at_edge": bool(abs(zr - float(b)) < 1e-9 * a)}


def rear_stagnation_point(U: float, a: float, alpha: float, beta: float, *, Gamma_cw: float) -> dict:
    """Where the rear stagnation point of the circle sits for a given circulation.

    Book: §14.3 (Fig. 14.9: no circulation → rear stagnation point on the suction side; Kutta value → at the edge) and
    §14.5, Eq. (14.10).  The trailing-edge point B is at θ = −(α + β).
    Parameters: U [m/s]; a circle radius [m]; alpha, beta [rad]; Gamma_cw [m²/s] (keyword).
    Returns dict: ``theta`` [rad] (from the stream-aligned diameter; −arcsin(Γ/4πUa)); ``phi`` = θ + α [rad] (polar
    angle about Q from the ξ direction); ``offset`` = θ + α + β [rad] (angle about Q ahead of the edge point B:
    > 0 on the upper surface, 0 at the edge, < 0 on the lower surface); ``arc_offset`` = a·offset [m];
    ``on_surface``; ``at_edge``.
    Assumptions: steady ideal flow.
    Validation: V1 offset = 0 exactly with (14.11) — tests/test_ch14.py:
    test_stagnation_points_V1_both_points_for_any_circulation.
    Label: analytic.
    """
    k = float(Gamma_cw) / (4.0 * np.pi * float(U) * float(a))
    on = abs(k) <= 1.0
    th = float(-np.arcsin(k)) if on else float("nan")  # Eq. (14.10) = 0
    off = th + float(alpha) + float(beta)
    return {"theta": th, "phi": th + float(alpha), "offset": off, "arc_offset": float(a) * off, "on_surface": bool(on),
            "at_edge": bool(on and abs(off) < 1e-12)}


def trailing_edge_speed(U, a, b, alpha, beta):
    """Finite speed at the cusped trailing edge when the Kutta condition holds, U (b/a) cos(α + β) [m/s].

    Book: §14.5 (text after (14.9): with dw/dζ = 0 at B, dw/dz is "zero-over-zero" and finite — the value is not
    given).  Ours, by l'Hôpital on (14.9): dw/dz → (d²w/dζ²)·ζ/2 at ζ = b = U (b/a) cos(α + β) e^{2iβ}; the flow leaves
    the edge at the angle −2β to the chord axis.
    Parameters: U [m/s]; a circle radius [m]; b map constant [m]; alpha, beta [rad].   Returns the speed [m/s].
    Assumptions: Kutta circulation (14.11); cusped (Zhukhovsky) edge.
    Validation: V1 limit of |dw/dz| along the cusp's bisector; the difference falls as √d (2.2e-4 of the limit at d =
    1e-6 m, 2.2e-5 at 1e-8 m for the chapter's section) — tests/test_ch14.py:
    test_trailing_edge_speed_V1_finite_limit_of_14_9_with_kutta_and_blow_up_without,
    test_kutta_circulation_V2_derivation, test_surface_pressure_V1_stagnation_points_bernoulli_and_the_edge_value.
    Label: analytic.
    """
    return _S(_F(U) * _F(b) / _F(a) * np.cos(_F(alpha) + _F(beta)))


# ======================================================================================================================
# the flow (§14.5)
# ======================================================================================================================
class ZhukhovskyFlow:
    """Ideal flow at the angle of attack α round a Zhukhovsky section, carried from the circle by (14.3), (14.9).

    Build with :func:`zhukhovsky_flow`.  Attributes: ``U`` [m/s], ``b`` [m], ``e``, ``beta``, ``alpha`` [rad], ``mu``
    (complex), ``a`` [m], ``Gamma_cw`` [m²/s], ``Gamma_kutta`` [m²/s], ``is_kutta`` (bool), ``u_inf`` = (U cos α,
    U sin α).  Methods take complex z (or x, y) of any shape; points inside the section return NaN.

    Book: §14.5, Fig. 14.16, Eq. (14.9) dw/dz = (dw/dζ) ζ²/(ζ² − b²), with the circle flow (6.68) shifted to the centre
    μ and turned by α, and the Kutta circulation Γ = 4πUa sin(α + β) (14.11) as the default.  Every method below is a
    part of that one flow (its Book line names the equation it evaluates).
    Assumptions: steady, incompressible, irrotational, two-dimensional; attached flow (no stall).
    Validation: as :func:`zhukhovsky_flow` — V1 the section is a streamline, far field U e^{−iα}, dw/dz equals
    differences of w(z); V4 ∮u·dx = −Γ_cw; V7 e = β = 0 reproduces the chapter-6 route — tests/test_ch14.py:
    test_kutta_condition_V1_root_finder_equals_14_11_and_the_edge_point_stagnates,
    test_zhukhovsky_flow_V1_section_is_a_streamline_far_field_and_complex_derivative,
    test_zhukhovsky_flow_V4_circulation_is_the_same_in_both_planes.
    Label: analytic.
    """

    def __init__(self, U, b, e, beta, alpha, Gamma_cw=None, a=None):
        self.U, self.b, self.beta, self.alpha = float(U), float(b), float(beta), float(alpha)
        self.mu, self.a = zhukhovsky_circle(b, e, beta, a=a)
        self.e = float(self.a * np.cos(self.beta) / self.b - 1.0)
        self.Gamma_kutta = float(kutta_circulation(self.U, self.a, self.alpha, self.beta))
        self.Gamma_cw = self.Gamma_kutta if Gamma_cw is None else float(Gamma_cw)
        scale = max(abs(self.Gamma_kutta), self.U * self.a, 1e-300)
        self.is_kutta = bool(abs(self.Gamma_cw - self.Gamma_kutta) <= 1e-12 * scale)
        self._w, self._dw = circle_flow_zeta(self.U, self.a, Gamma_cw=self.Gamma_cw, alpha=self.alpha, center=self.mu)
        self.u_inf = (self.U * np.cos(self.alpha), self.U * np.sin(self.alpha))
        self.label = "Zhukhovsky airfoil"

    # --- ζ-plane ---
    def w_zeta(self, zeta):
        """Complex potential w(ζ) = U[(ζ − μ)e^{−iα} + a²e^{iα}/(ζ − μ)] + (iΓ/2π) ln((ζ − μ)/a) [m²/s].

        Book: §14.5 ((6.68) shifted to μ and turned by α).  Label: analytic.
        """
        return self._w(zeta)

    def dwdzeta(self, zeta):
        """Complex velocity of the circle flow dw/dζ [m/s].

        Book: §14.5, first factor of (14.9); on the circle it is (14.10).  Label: analytic.
        """
        return self._dw(zeta)

    def zeta_of(self, z):
        """ζ(z): the root outside the shifted circle (NaN inside the section).

        Book: §14.4, inverse of z = ζ + b²/ζ (14.3) (:func:`zhukhovsky_inverse`).  Label: analytic.
        """
        return zhukhovsky_inverse(z, self.b, self.mu, self.a)

    # --- z-plane ---
    def w(self, z):
        """Complex potential φ + iψ at z [m²/s] (φ jumps by Γ across the cut of the logarithm; ψ is continuous).

        Book: §14.5 (w(ζ(z)): the potential is carried unchanged by the map (14.3)).  Label: analytic.
        """
        return _S(self._w(_C(self.zeta_of(z))))

    def dwdz_from_zeta(self, zeta):
        """dw/dz = (dw/dζ) ζ²/(ζ² − b²) (14.9) at a ζ-plane point; at ζ = b the analytic limit (Kutta) or ∞.

        Book: §14.5, Eq. (14.9); the limit is :func:`trailing_edge_speed` (ours).  Label: analytic.
        """
        zeta = _C(zeta)
        num = _C(self._dw(zeta)) * zeta ** 2
        den = zeta ** 2 - self.b ** 2
        edge = np.abs(zeta - self.b) < 1e-13 * self.a
        with np.errstate(divide="ignore", invalid="ignore"):
            q = num / np.where(edge, 1.0, den)  # Eq. (14.9)
        if np.any(edge):
            lim = (trailing_edge_speed(self.U, self.a, self.b, self.alpha, self.beta) * np.exp(2j * self.beta)
                   if self.is_kutta else np.inf + 0j)
            q = np.where(edge, lim, q)
        return _S(q)

    def dwdz(self, z):
        """Complex velocity u − iv at z [m/s], Eq. (14.9).

        Book: §14.5, Eq. (14.9) dw/dz = (dw/dζ) ζ²/(ζ² − b²).  Label: analytic.
        """
        return self.dwdz_from_zeta(self.zeta_of(z))

    def velocity(self, x, y):
        """(u, v) [m/s] at (x, y) [m].

        Book: §14.5, Eq. (14.9) with dw/dz = u − iv.  Label: analytic.
        """
        q = _C(self.dwdz(_F(x) + 1j * _F(y)))
        return _S(q.real), _S(-q.imag)

    def speed(self, z):
        """|dw/dz| [m/s] at complex z.

        Book: §14.5, modulus of (14.9).  Label: analytic.
        """
        return _S(np.abs(_C(self.dwdz(z))))

    def stream(self, x, y):
        """Stream function ψ = Im w [m²/s] (constant on the section).

        Book: §14.5 (w = φ + iψ, (6.42)).  Label: analytic.
        """
        return _S(_C(self.w(_F(x) + 1j * _F(y))).imag)

    def phi(self, x, y):
        """Velocity potential φ = Re w [m²/s].

        Book: §14.5 (w = φ + iψ, (6.42)).  Label: analytic.
        """
        return _S(_C(self.w(_F(x) + 1j * _F(y))).real)

    def psi(self, x, y):
        """Alias of :meth:`stream`.

        Book: §14.5 (w = φ + iψ, (6.42)).  Label: analytic.
        """
        return self.stream(x, y)

    def inside(self, x, y):
        """True inside the section.

        Book: §14.4 (the section is the image of the circle under (14.3)).  Label: analytic.
        """
        return np.isnan(_C(self.zeta_of(_F(x) + 1j * _F(y))).real)

    def speed_near_edge(self, d):
        """Speed at the distance d [m] behind the trailing edge along the cusp's bisector (direction e^{−2iβ}).

        Book: §14.5, text after (14.9) (finite with the Kutta circulation, infinite otherwise).  Label: analytic.
        """
        z = 2.0 * self.b + _F(d) * np.exp(-2j * self.beta)
        return self.speed(z)

    def cp(self, z):
        """Pressure coefficient 1 − |dw/dz|²/U² (Bernoulli) at complex z.

        Book: §14.3, Fig. 14.8 (C_p round a section), with (14.9) and Bernoulli's equation.  Label: analytic.
        """
        return _S(1.0 - (np.abs(_C(self.dwdz(z))) / self.U) ** 2)


def zhukhovsky_flow(U: float, b: float, e: float, beta: float, alpha: float, *, Gamma_cw: float | None = None,
                    a: float | None = None) -> ZhukhovskyFlow:
    """Ideal flow round a Zhukhovsky section at the angle of attack α with clockwise circulation Γ.

    Book: §14.5, Fig. 14.16 and Eq. (14.9) dw/dz = (dw/dζ) ζ²/(ζ² − b²); the ζ-plane flow is the circle flow (6.68)
    shifted to the centre μ and turned by α (the book writes only the centred case, Exercise 14.5):
    w = U[(ζ − μ)e^{−iα} + a²e^{iα}/(ζ − μ)] + (iΓ/2π) ln((ζ − μ)/a).  The circulation and the far stream are the same
    in both planes (z ≃ ζ at infinity; (dw/dz)dz = (dw/dζ)dζ).
    Parameters
    ----------
    U [m/s]; b [m]; e [–]; beta [rad]; alpha [rad]; Gamma_cw [m²/s] (None → the Kutta value (14.11));
    a : optional circle radius [m] (then e is derived).
    Returns
    -------
    :class:`ZhukhovskyFlow` with ``w(z)``, ``dwdz(z)``, ``dwdzeta(zeta)``, ``velocity(x, y)``, ``speed(z)``,
    ``stream(x, y)``, ``speed_near_edge(d)``, ``zeta_of(z)``, ``inside(x, y)``, ``Gamma_cw``, ``mu``, ``a``.
    Assumptions: steady, incompressible, irrotational, two-dimensional; attached flow (no stall).
    Validation: V1 section is a streamline; far field U e^{−iα}; V4 ∮u·dx = −Γ_cw on any enclosing contour; V7 e = β = 0
    reproduces ``circle_flow_zeta(alpha=)`` through (6.69) — tests/test_ch14.py:
    test_zhukhovsky_flow_V1_section_is_a_streamline_far_field_and_complex_derivative,
    test_zhukhovsky_flow_V4_circulation_is_the_same_in_both_planes,
    test_zhukhovsky_flow_V7_plate_limit_reproduces_the_chapter_6_route.
    Label: analytic.
    """
    return ZhukhovskyFlow(U, b, e, beta, alpha, Gamma_cw, a)


def flat_plate_flow(U: float, b: float, alpha: float, *, Gamma_cw: float | None = None) -> ZhukhovskyFlow:
    """Flat plate of chord 4b at incidence: :func:`zhukhovsky_flow` with e = β = 0 (centred circle of radius b).

    Book: §14.5, Exercise 14.5: w = U(ζe^{−iα} + (b²/ζ)e^{iα}) + (iΓ/2π) ln(ζe^{−iα}), Γ = 4πUb sin α, C_L = 2π sin α.
    (Our logarithm is ln(ζ/b): it differs from the printed one by an additive constant, which changes no velocity.)
    Parameters: U [m/s]; b [m]; alpha [rad]; Gamma_cw [m²/s] (None → Kutta).   Returns a :class:`ZhukhovskyFlow`.
    Assumptions: as :func:`zhukhovsky_flow`; zero thickness (infinite speed at the leading edge).
    Validation: V1 Γ_K = 4πUb sin α; Blasius lift coefficient 2π sin α — tests/test_ch14.py:
    test_zhukhovsky_flow_V7_plate_limit_reproduces_the_chapter_6_route,
    test_flat_plate_V1_leading_edge_suction_closes_the_force_triangle.
    Label: analytic.
    """
    return ZhukhovskyFlow(U, b, 0.0, 0.0, alpha, Gamma_cw)


def zhukhovsky_lift_coefficient(alpha, beta, e: float = 0.0, exact: bool = False):
    """Lift coefficient of a Zhukhovsky section: C_L ≃ 2π(α + β) (14.12), or the exact 8π(a/c) sin(α + β).

    Book: §14.5, Eq. (14.12) C_L = L/(½ρU²c) ≃ 2π(α + β), "where we have used 4a ≃ c, L = ρUΓ and
    sin(α + β) ≃ (α + β)".  ``exact=True`` (ours, slip #8): L = ρUΓ with (14.11) and the chord measured on the contour
    (:func:`zhukhovsky_shape`), C_L = 2Γ/(Uc) = 8π(a/c) sin(α + β); thickness raises the slope above 2π.
    Parameters: alpha [rad]; beta [rad]; e [–] — **used only with** ``exact=True``: the rule (14.12) has no thickness
    in it, so ``exact=False`` ignores ``e`` without a warning; exact (bool).   Returns C_L [–].
    Assumptions: Kutta condition, no stall; (approximate form) small α + β, e ≪ 1.
    Validation: V1 exact = 2Γ_K/(Uc); flat plate 2π sin α; V5 slope 2π per radian (thin-airfoil theory) —
    tests/test_ch14.py: test_zhukhovsky_lift_V1_blasius_force_is_rho_U_gamma_normal_to_the_stream,
    test_zhukhovsky_lift_V5_slope_two_pi_per_radian_and_the_quarter_chord,
    test_zhukhovsky_lift_V7_small_angle_rule_and_slip_8_thickness_raises_the_slope.
    Label: analytic.
    """
    al, be = _F(alpha), _F(beta)
    if not exact:
        return _S(2.0 * np.pi * (al + be))  # Eq. (14.12)
    if be.ndim:
        raise ValueError("exact=True needs a scalar beta (the chord is measured per section)")
    sh = zhukhovsky_shape(1.0, e, float(be), exact=True)
    return _S(8.0 * np.pi * sh["a"] / sh["chord"] * np.sin(al + be))


def zhukhovsky_surface(U: float, b: float, e: float, beta: float, alpha: float, n: int = 400, *,
                       Gamma_cw: float | None = None) -> dict:
    """Speed and pressure coefficient on the surface of a Zhukhovsky section.

    Book: §14.3, Fig. 14.8 (C_p round a section — drawn, no formula) computed with §14.5, Eqs. (14.9), (14.10):
    on the circle the speed is |u_θ|, u_θ = −2U sin θ − Γ/(2πa) (14.10) with θ = t − α; on the section it is
    |u_θ| ζ²/(ζ² − b²) in modulus (14.9); C_p = 1 − |dw/dz|²/U² (Bernoulli).
    Parameters
    ----------
    U [m/s]; b [m]; e [–]; beta, alpha [rad]; n : contour points (index 0 = trailing edge, counterclockwise);
    Gamma_cw [m²/s] (None → Kutta).
    Returns
    -------
    dict: ``t``, ``zeta``, ``z``, ``x``, ``y``; ``s`` arc length from the trailing edge [m] (counterclockwise);
    ``ds_dt`` = |dz/dt| [m]; ``u_theta`` (14.10) on the circle [m/s]; ``u_s`` signed surface speed along the
    counterclockwise traversal [m/s]; ``speed`` [m/s]; ``cp`` [–]; ``xc`` = (x − x_min)/(x_max − x_min);
    ``chord_x`` = x_max − x_min [m]; ``upper``, ``lower``: dicts (``xc``, ``cp``, ``x``, ``y``, ``speed``, ``index``)
    ordered from the leading to the trailing edge; ``stagnation`` (list of the z-plane stagnation points on the surface:
    front, then rear — the rear one is absent at the Kutta value, where its image is the cusp and the speed is finite); ``Gamma_cw``; ``is_kutta``; ``cp_min``, ``xc_cp_min`` (suction peak).
    At the trailing edge: the analytic limit U(b/a)cos(α + β) when Γ is the Kutta value (u_s there is given the
    upper-surface sign, negative), ``np.inf`` otherwise.
    Assumptions: steady ideal flow, attached.
    Validation: V1 C_p = 1 at both stagnation points; V4 −∮u_s ds = Γ_cw; the C_p loop gives F_y of the Blasius route —
    tests/test_ch14.py: test_surface_pressure_V1_stagnation_points_bernoulli_and_the_edge_value,
    test_surface_pressure_V4_loop_area_is_the_normal_force_and_the_speed_integral_is_the_circulation.
    Label: analytic.
    """
    A = zhukhovsky_airfoil(b, e, beta, int(n))
    F = ZhukhovskyFlow(U, b, e, beta, alpha, Gamma_cw)
    t, zeta, z = A["t"], A["zeta"], A["z"]
    a, b, beta = F.a, F.b, F.beta
    u_th = -2.0 * F.U * np.sin(t - F.alpha) - F.Gamma_cw / (_TWO_PI * a)  # Eq. (14.10)
    with np.errstate(divide="ignore", invalid="ignore"):
        jac = np.abs(1.0 - b ** 2 / zeta ** 2)  # |dz/dζ|
        u_s = u_th / jac  # Eq. (14.9) on the circle
    ds_dt = a * jac
    te_lim = float(trailing_edge_speed(F.U, a, b, F.alpha, beta)) if F.is_kutta else np.inf
    u_s[0] = -te_lim
    speed = np.abs(u_s)
    cp = 1.0 - (speed / F.U) ** 2
    seg = np.abs(np.diff(np.concatenate([z, z[:1]])))
    s = np.concatenate([[0.0], np.cumsum(seg)[:-1]])
    x, y = z.real, z.imag
    xmin, xmax = float(x.min()), float(x.max())
    cx = xmax - xmin
    xc = (x - xmin) / cx
    il = int(np.argmin(x))
    iu = np.arange(il, -1, -1)  # leading → trailing edge over the upper surface
    ilow = np.concatenate([np.arange(il, x.size), [0]])  # leading → trailing edge under the lower surface

    def side(ix):
        return {"xc": xc[ix], "cp": cp[ix], "x": x[ix], "y": y[ix], "speed": speed[ix], "index": ix}

    st = stagnation_points(F.U, b, F.e, beta, F.alpha, Gamma_cw=F.Gamma_cw)
    stag = []
    if st["on_surface"]:  # at a cusped edge the image of the rear stagnation point is an ordinary point (finite speed)
        stag = [st["z_front"]] if st["rear_at_edge"] else [st["z_front"], st["z_rear"]]
    fin = np.isfinite(cp)
    k = int(np.argmin(np.where(fin, cp, np.inf)))
    return {"t": t, "zeta": zeta, "z": z, "x": x, "y": y, "s": s, "ds_dt": ds_dt, "u_theta": u_th, "u_s": u_s,
            "speed": speed, "cp": cp, "xc": xc, "chord_x": cx, "upper": side(iu), "lower": side(ilow),
            "stagnation": stag, "Gamma_cw": F.Gamma_cw, "is_kutta": F.is_kutta, "cp_min": float(cp[k]),
            "xc_cp_min": float(xc[k]), "U": F.U, "a": a, "b": b, "alpha": F.alpha, "beta": beta, "e": F.e}


def surface_circulation(surface: dict, method: str = "polygon") -> float:
    """Clockwise circulation as the integral of the surface speed round the contour, Γ_cw = −∮ u_s ds [m²/s].

    Book: §14.3 (lift development = circulation development; in a real flow the circulation is the net vorticity of
    the two boundary layers) with Γ ≡ ∮u·ds (3.18).  ``method="polygon"``: trapezoid rule on the contour polygon with
    the tangential velocity at the vertices (second order; needs a finite trailing-edge speed, i.e. the Kutta value);
    ``method="circle"``: ∮ u_θ a dt in the ζ-plane (exact to rounding for any Γ, because (dw/dz)dz = (dw/dζ)dζ).
    Parameters: surface — the dict of :func:`zhukhovsky_surface`.   Returns Γ_cw [m²/s].
    Assumptions: the contour is traversed counterclockwise.
    Validation: V4 both methods return the Γ_cw the surface was built with (1e-6 at n = 2000 for the polygon) —
    tests/test_ch14.py:
    test_surface_pressure_V4_loop_area_is_the_normal_force_and_the_speed_integral_is_the_circulation.
    Label: analytic, converged.
    """
    if method == "circle":
        dt = _TWO_PI / surface["t"].size
        return float(-np.sum(surface["u_theta"]) * surface["a"] * dt)
    if method != "polygon":
        raise ValueError('method must be "polygon" or "circle"')
    z, us = surface["z"], np.array(surface["u_s"], dtype=float)
    zc = np.concatenate([z, z[:1]])
    ds = np.abs(np.diff(zc))
    us_end = np.concatenate([us, [-us[0]]])  # the trailing edge closes the loop from the lower side: sign flips
    return float(-np.sum(0.5 * (us_end[:-1] + us_end[1:]) * ds))


def zhukhovsky_force(U: float, b: float, e: float, beta: float, alpha: float, rho: float = 1.2,
                     method: str = "blasius", *, Gamma_cw: float | None = None, n: int = 2000) -> dict:
    """Force and moment per unit span on a Zhukhovsky section, by the Blasius contour integral or by surface pressure.

    Book: §14.3 (L = ρUΓ, (6.62)) and §6.5, Eq. (6.60) D_x − iL_y = (iρ/2)∮(dw/dz)² dz on a contour enclosing the
    section; the moment is the companion theorem M₀ = −(ρ/2) Re ∮ z (dw/dz)² dz (counterclockwise about the origin;
    ours, standard).  ``method="pressure"``: F = i∮ p dz and M₀ = Im ∮ z̄ dF with p − p_∞ = ½ρU²C_p on the surface
    polygon (second order in 1/n; the trailing-edge speed must be finite, so Γ must be the Kutta value).
    Parameters
    ----------
    U [m/s]; b [m]; e [–]; beta, alpha [rad]; rho [kg/m³]; method "blasius" | "pressure"; Gamma_cw [m²/s] (None →
    Kutta); n : surface points (pressure route) — the Blasius circle uses 512 nodes of the periodic trapezoid rule.
    Returns
    -------
    dict: ``Fx``, ``Fy`` [N/m] (section axes: x along the mapped chord); ``L``, ``D`` [N/m] (normal to / along the
    stream); ``M0`` [N] (moment about the origin, counterclockwise +); ``M_le`` [N] (about the leading edge,
    counterclockwise +; nose-up is −M_le for a stream from the left); ``x_cp`` [m] (where the line of action crosses
    the chord axis, M0/Fy); ``xcp_over_c`` = (x_cp − x of the leading edge)/chord — measured along the x-axis (the
    mapped chord axis) from the leading edge's x, not along the chord line of a cambered section; ``Gamma_cw``;
    ``CL``, ``CD`` (based on the measured chord); ``chord``.
    Call by keyword — the sibling function :func:`flat_plate_forces` ``(rho, U, c, alpha)`` takes these arguments in
    a different order (here ρ comes after the geometry and the angles).
    Assumptions: steady ideal flow.
    Validation: V1 Blasius L = ρUΓ, D = 0 (1e-10) for α ≠ 0; V3 pressure route converges to it at second order; flat
    plate: centre of pressure at the quarter chord — tests/test_ch14.py:
    test_zhukhovsky_lift_V1_blasius_force_is_rho_U_gamma_normal_to_the_stream,
    test_zhukhovsky_lift_V3_pressure_route_converges_to_the_blasius_force_at_second_order,
    test_zhukhovsky_lift_V5_slope_two_pi_per_radian_and_the_quarter_chord.
    Label: analytic, converged.
    """
    F = ZhukhovskyFlow(U, b, e, beta, alpha, Gamma_cw)
    sh = zhukhovsky_shape(b, F.e, beta, exact=True, n=4001)
    A = zhukhovsky_airfoil(b, F.e, beta, 4001)
    g = section_geometry(A["x"], A["y"], te_index=0)
    le = g["le"]
    if method == "blasius":
        R = 3.0 * max(abs(F.mu) + F.a, 2.0 * F.b)
        bf = blasius_force(F.dwdz, R=R, rho=rho, n=512)  # Eq. (6.60)
        Fx, Fy = bf.D, bf.L
        th = np.linspace(0.0, _TWO_PI, 512, endpoint=False)
        zc = R * np.exp(1j * th)
        dz = 1j * zc * (_TWO_PI / 512)
        M0 = float((-0.5 * rho * np.sum(zc * _C(F.dwdz(zc)) ** 2 * dz)).real)
    elif method == "pressure":
        S = zhukhovsky_surface(U, b, F.e, beta, alpha, int(n), Gamma_cw=F.Gamma_cw)
        p = 0.5 * rho * F.U ** 2 * S["cp"]
        z = S["z"]
        zc, pc = np.concatenate([z, z[:1]]), np.concatenate([p, p[:1]])
        dF = 1j * 0.5 * (pc[:-1] + pc[1:]) * np.diff(zc)  # dF = −p n ds = i p dz
        Ft = np.sum(dF)
        Fx, Fy = float(Ft.real), float(Ft.imag)
        zm = 0.5 * (zc[:-1] + zc[1:])
        M0 = float(np.sum((np.conj(zm) * dF).imag))
    else:
        raise ValueError('method must be "blasius" or "pressure"')
    L, D = resolve_lift_drag(Fx, Fy, F.alpha)
    M_le = M0 - float((np.conj(le) * (Fx + 1j * Fy)).imag)
    x_cp = M0 / Fy if Fy != 0.0 else float("nan")
    q = 0.5 * rho * F.U ** 2 * sh["chord"]
    return {"Fx": float(Fx), "Fy": float(Fy), "L": float(L), "D": float(D), "M0": M0, "M_le": float(M_le),
            "x_cp": float(x_cp), "xcp_over_c": float((x_cp - le.real) / sh["chord"]), "Gamma_cw": F.Gamma_cw,
            "CL": float(L / q), "CD": float(D / q), "chord": sh["chord"], "method": method}


def flat_plate_forces(rho: float, U: float, c: float, alpha: float) -> dict:
    """Lift, plate-normal force and leading-edge suction on a flat plate at incidence in ideal flow.

    Book: §14.3 (a rounded leading edge "recovers the leading edge suction that occurs on foils of negligible
    thickness"), Exercises 14.5 and 14.7: with the Kutta circulation Γ = πUc sin α (c = 4b in (14.11)) the lift is
    L = ρUΓ = πρU²c sin α, normal to the stream.  Pressure can only push normal to the plate, N = L cos α; the missing
    chordwise part is the suction at the sharp leading edge, S = L sin α = πρU²c sin²α, directed along the plate toward
    the leading edge, so that the drag is zero.
    Parameters: rho [kg/m³]; U [m/s]; c chord [m]; alpha [rad].
    Returns dict: ``L``, ``N``, ``S``, ``D`` (= 0) [N/m]; ``Gamma_cw`` [m²/s]; ``CL`` = 2π sin α.
    Call by keyword — the sibling function :func:`zhukhovsky_force` ``(U, b, e, beta, alpha, rho)`` takes these
    arguments in a different order (here ρ comes first, and the length is the chord c = 4b, not b).
    Assumptions: steady ideal flow, Kutta condition, zero thickness.
    Validation: V1 N cos α + S sin α = L and N sin α − S cos α = 0; V4 equals the control-volume force of
    ``ch14.edge_flow_suction`` with U_o²a = U²c sin²α — tests/test_ch14.py:
    test_flat_plate_V1_leading_edge_suction_closes_the_force_triangle,
    test_starting_vortex_V1_steady_limit_is_the_kutta_value_and_first_step_is_one_third.
    Label: analytic.
    """
    G = np.pi * float(U) * float(c) * np.sin(alpha)  # Eq. (14.11) with a = b = c/4, β = 0
    L = float(rho) * float(U) * G  # Eq. (6.62)
    return {"L": float(L), "N": float(L * np.cos(alpha)), "S": float(L * np.sin(alpha)), "D": 0.0,
            "Gamma_cw": float(G), "CL": float(_TWO_PI * np.sin(alpha))}


def thin_airfoil(dyc_dx: Callable, alpha: float, n: int = 200) -> dict:
    """Thin-airfoil theory for a camber line: Fourier coefficients, lift coefficient, zero-lift angle, moment.

    Book: §14.4 names thin-airfoil theory (a vortex sheet on the camber line, "camber line is a streamline", zero sheet
    strength at the trailing edge) and gives no formula.  Ours — the standard result (Glauert): with
    x/c = ½(1 − cos ϑ), γ(ϑ) = 2U[A₀(1 + cos ϑ)/sin ϑ + Σ A_n sin nϑ], A₀ = α − (1/π)∫₀^π (dy_c/dx) dϑ,
    A_n = (2/π)∫₀^π (dy_c/dx) cos nϑ dϑ, C_L = 2π(A₀ + ½A₁), α_{L=0} = −(1/π)∫₀^π (dy_c/dx)(cos ϑ − 1) dϑ,
    C_M about the quarter chord = (π/4)(A₂ − A₁).
    Parameters: dyc_dx callable of x/c ∈ [0, 1] → slope of the camber line [–]; alpha [rad]; n midpoint nodes in ϑ
    (the Gauss–Chebyshev rule: spectrally accurate for smooth slopes).
    Returns dict: ``A0``, ``A1``, ``A2``, ``CL``, ``alpha_L0`` [rad], ``CM_c4``, ``slope`` (= 2π).
    Assumptions: thin section, small camber and α, Kutta condition; NOT a book formula.
    Validation: V1 flat plate C_L = 2πα; parabolic arc of height h: α_{L=0} = −2h/c (= −β of the Zhukhovsky arc);
    C_M,c/4 independent of α — tests/test_ch14.py:
    test_zhukhovsky_lift_V5_slope_two_pi_per_radian_and_the_quarter_chord,
    test_thin_airfoil_V1_plate_parabolic_arc_and_exact_coefficients_of_a_cubic_camber_line,
    test_thin_airfoil_V3_discrete_vortex_method_gives_the_same_lift_and_moment.
    Label: analytic (ours).
    """
    n = int(n)
    th = (np.arange(n) + 0.5) * np.pi / n
    sl = _F(dyc_dx(0.5 * (1.0 - np.cos(th)))) + 0.0 * th
    h = np.pi / n
    A0 = float(alpha) - np.sum(sl) * h / np.pi
    A1 = 2.0 / np.pi * np.sum(sl * np.cos(th)) * h
    A2 = 2.0 / np.pi * np.sum(sl * np.cos(2.0 * th)) * h
    aL0 = -np.sum(sl * (np.cos(th) - 1.0)) * h / np.pi
    return {"A0": float(A0), "A1": float(A1), "A2": float(A2), "CL": float(_TWO_PI * (A0 + 0.5 * A1)),
            "alpha_L0": float(aL0), "CM_c4": float(np.pi / 4.0 * (A2 - A1)), "slope": _TWO_PI}
