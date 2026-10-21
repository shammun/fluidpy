"""Evidence tests for ``fluidpy.geo.frames`` and the rotation constants of ``fluidpy.geo.constants`` (AOFD PoC B and D:
Coriolis on a turntable and a globe, f-plane inertial circles, Coriolis-only motion on the sphere).

Independent routes used here (none shares code with the library): rotating-frame ODEs written from Vallis (2.12),
(3.126) and (2.50a,b) and integrated with scipy DOP853; sympy re-derivations; exact invariants (speed, the Jacobi
integral on the turntable, axial absolute angular momentum on the sphere); great-circle geometry in 3-D Cartesian
form; the f-plane and Omega -> 0 limits; hemisphere mirror symmetry.

Benchmarks (cited, verified 2026-10-21; see AOFD ``reference/ch00/SOURCES.md``):
* Omega = 7 292 115e-11 rad/s — NGA, WGS 84 (TR8350.2) defining parameter.
* Inertial period at 45 deg ~= 16.9 h (cleonis.nl inertial-oscillation table; ~17 h "middle latitudes", Charles
  University hydrodynamics exercise 10) — given to 0.1 h.
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import solve_ivp

from fluidpy.core import gfd as core_gfd
from fluidpy.core import rotating as core_rot
from fluidpy.geo import constants as C
from fluidpy.geo import frames as F
from tools.convergence import observed_order, pairwise_orders

OM, A_E = C.OMEGA, C.EARTH_RADIUS
LAT45 = np.deg2rad(45.0)
RNG = np.random.default_rng(7)


def _wrap(x):
    return (np.asarray(x) + np.pi) % (2 * np.pi) - np.pi


# ---------------------------------------------------------------------------------------------------------- constants
def test_constants_V5_wgs84_omega_and_sidereal_day():  # V5 NGA WGS 84: 7 292 115e-11 rad/s
    assert C.OMEGA == 7292115e-11
    assert C.OMEGA is core_rot.OMEGA_EARTH or C.OMEGA == core_rot.OMEGA_EARTH  # one home, re-exported
    assert C.SIDEREAL_DAY == pytest.approx(2 * np.pi / C.OMEGA, rel=1e-15)
    assert 86164.0 < C.SIDEREAL_DAY < 86164.2  # 23 h 56 min 4 s (2 pi / Omega with the rounded WGS 84 Omega)
    assert C.EARTH_RADIUS == 6.371e6


# --------------------------------------------------------------------------------------------------------- f and beta
def test_coriolis_parameter_V1_field_and_V7_odd():  # V1 f = 2 Omega sin(lat) (Vallis 2.80); V7 odd, f(0) = 0
    lat = np.linspace(-np.pi / 2, np.pi / 2, 721)
    f = F.coriolis_parameter(lat)
    assert np.max(np.abs(f - 2 * OM * np.sin(lat))) < 1e-20
    assert np.array_equal(F.coriolis_parameter(-lat), -f)
    assert F.coriolis_parameter(0.0) == 0.0  # equator limit f -> 0
    assert F.coriolis_parameter(np.pi / 2) == pytest.approx(2 * OM, rel=1e-15)


def test_coriolis_parameter_V5_value_at_45deg():  # V5: f(45 deg) from the WGS 84 Omega
    f45 = F.coriolis_parameter(LAT45)
    assert f45 == pytest.approx(1.0313e-4, rel=1e-4)
    assert f45 == pytest.approx(np.sqrt(2) * 7292115e-11, rel=1e-15)


def test_beta_parameter_V1_complex_step_and_V7_even():  # V1 complex-step d/dy of 2 Omega sin(lat0 + y/a)
    lat = np.linspace(-1.5, 1.5, 61)
    h = 1e-20
    beta_cs = np.imag(2 * OM * np.sin(lat + 1j * h / A_E)) / h
    assert np.allclose(F.beta_parameter(lat), beta_cs, rtol=1e-14, atol=0)
    assert np.array_equal(F.beta_parameter(-lat), F.beta_parameter(lat))
    assert F.beta_parameter(LAT45) == pytest.approx(1.6187e-11, rel=1e-4)
    assert F.beta_parameter(np.pi / 2) == pytest.approx(0.0, abs=1e-26)


# ------------------------------------------------------------------------------------- horizontal Coriolis acceleration
def test_coriolis_horizontal_V2_from_full_rotation_vector():  # V2: horizontal part of -2 Omega x u, traditional approx.
    Om, th, u, v = sp.symbols("Omega vartheta u v", real=True)
    Omv = sp.Matrix([0, Om * sp.cos(th), Om * sp.sin(th)])  # local (east, north, up) components
    acc = -2 * Omv.cross(sp.Matrix([u, v, 0]))
    assert sp.simplify(acc[0] - 2 * Om * sp.sin(th) * v) == 0
    assert sp.simplify(acc[1] + 2 * Om * sp.sin(th) * u) == 0
    # acc[2] = 2 Omega cos(lat) u is the non-traditional vertical term the function (rightly) omits
    ex = sp.lambdify((th, u, v), acc[0].subs(Om, OM)), sp.lambdify((th, u, v), acc[1].subs(Om, OM))
    lat, uu, vv = RNG.uniform(-1.5, 1.5, 300), RNG.normal(0, 20, 300), RNG.normal(0, 20, 300)
    ax, ay = F.coriolis_horizontal(uu, vv, lat)
    assert np.max(np.abs(ax - ex[0](lat, uu, vv))) < 1e-17 and np.max(np.abs(ay - ex[1](lat, uu, vv))) < 1e-17


def test_coriolis_acceleration_V1_perpendicular_magnitude_and_core():  # V1 a.u = 0, |a| = |f||u|; core 3-D form
    u, v, f = RNG.normal(0, 10, 200), RNG.normal(0, 10, 200), RNG.normal(0, 1e-4, 200)
    ax, ay = F.coriolis_acceleration_horizontal(u, v, f)
    assert np.max(np.abs(ax * u + ay * v) / (np.abs(f) * (u * u + v * v))) < 1e-15  # no work done (round-off)
    assert np.allclose(np.hypot(ax, ay), np.abs(f) * np.hypot(u, v), rtol=1e-14)
    full = core_rot.coriolis_acceleration([0.0, 0.0, 0.5e-4], [3.0, -2.0, 0.0])
    assert np.allclose(F.coriolis_acceleration_horizontal(3.0, -2.0, 1e-4), full[:2], rtol=1e-15)


def test_coriolis_horizontal_V7_right_in_NH_left_in_SH_zero_at_equator():  # V7 hemisphere mirror
    for lat, side in ((LAT45, -1), (-LAT45, +1)):
        ax, ay = F.coriolis_horizontal(10.0, 0.0, lat)  # eastward wind
        assert ax == 0.0 and np.sign(ay) == side  # NH: southward = to the right of eastward motion
    axn, ayn = F.coriolis_horizontal(3.0, 4.0, 0.6)
    axs, ays = F.coriolis_horizontal(3.0, 4.0, -0.6)
    assert (axs, ays) == (-axn, -ayn)
    assert F.coriolis_horizontal(3.0, 4.0, 0.0) == (0.0, 0.0)
    assert 3.0 * ayn - 4.0 * axn < 0  # (u x a)_z < 0: clockwise turn seen from above = to the right


# ------------------------------------------------------------------------------------------------------- turntable puck
def _turntable_ode(Om, x0, y0, u0, v0, t):
    """Vallis (2.12) on a flat table: dv/dt = -2 Omega x v - Omega x (Omega x r) — written out independently."""
    def rhs(_, s):
        x, y, u, v = s
        return [u, v, 2 * Om * v + Om**2 * x, -2 * Om * u + Om**2 * y]
    return solve_ivp(rhs, (t[0], t[-1]), [x0, y0, u0, v0], t_eval=t, method="DOP853", rtol=1e-12, atol=1e-14).y


@pytest.mark.parametrize("Om,x0,y0,u0,v0", [(0.5, 1.0, 0.3, 0.2, -0.4), (-1.2, -0.4, 0.8, 0.0, 0.5),
                                            (2.0, 0.0, 0.0, 1.0, 0.0)])
def test_turntable_puck_V1_matches_rotating_frame_ode(Om, x0, y0, u0, v0):  # V1 independent ODE solve of (2.12)
    t = np.linspace(0.0, 10.0, 401)
    r = F.turntable_puck(t, x0, y0, u0, v0, Om)
    y = _turntable_ode(Om, x0, y0, u0, v0, t)
    scale = 1 + np.max(np.hypot(y[0], y[1]))
    assert np.max(np.abs(r["x_rotating"] - y[0])) / scale < 1e-9
    assert np.max(np.abs(r["y_rotating"] - y[1])) / scale < 1e-9


def test_turntable_puck_V4_jacobi_integral_and_radius():  # V4: E_J = |v_R|^2/2 - Omega^2 |r|^2/2 is conserved
    Om, x0, y0, u0, v0 = 0.7, 1.0, -0.5, 0.3, 0.9
    t = np.linspace(0.0, 12.0, 241)
    h = 1e-3
    P = [F.turntable_puck(t + k * h, x0, y0, u0, v0, Om) for k in (-2, -1, 1, 2)]
    d = lambda key: (P[0][key] - 8 * P[1][key] + 8 * P[2][key] - P[3][key]) / (12 * h)  # 4th-order central FD
    r0 = F.turntable_puck(t, x0, y0, u0, v0, Om)
    EJ = 0.5 * (d("x_rotating") ** 2 + d("y_rotating") ** 2) - 0.5 * Om**2 * (r0["x_rotating"] ** 2 + r0["y_rotating"] ** 2)
    EJ0 = 0.5 * (u0**2 + v0**2) - 0.5 * Om**2 * (x0**2 + y0**2)
    assert np.max(np.abs(EJ - EJ0)) / (0.5 * Om**2 * 25) < 1e-9
    # rotation preserves distance from the axis; the rotating-frame launch velocity is (u0, v0) at t = 0
    assert np.allclose(np.hypot(r0["x_rotating"], r0["y_rotating"]), np.hypot(r0["x_inertial"], r0["y_inertial"]),
                       rtol=1e-14)
    assert d("x_rotating")[0] == pytest.approx(u0, abs=1e-10) and d("y_rotating")[0] == pytest.approx(v0, abs=1e-10)


def test_turntable_puck_V7_limits_and_deflection_sense():  # V7 Omega -> 0 straight line; NH right, SH left (mirror)
    t = np.linspace(0, 3, 31)
    r = F.turntable_puck(t, 0.2, 0.1, 1.0, 0.5, 0.0)
    assert np.allclose(r["x_rotating"], 0.2 + t) and np.allclose(r["y_rotating"], 0.1 + 0.5 * t)
    p = F.turntable_puck(t, 0.0, 0.0, 1.0, 0.0, 0.3)
    m = F.turntable_puck(t, 0.0, 0.0, 1.0, 0.0, -0.3)
    assert np.all(p["y_rotating"][1:] < 0)  # anticlockwise table (NH-like): drifts right of +x
    assert np.allclose(m["y_rotating"], -p["y_rotating"], atol=1e-15) and np.allclose(m["x_rotating"], p["x_rotating"])
    # the same physical launch described in either frame gives the same path
    q = F.turntable_puck(t, 0.4, -0.2, 0.3, 0.6, 0.8)
    qi = F.turntable_puck(t, 0.4, -0.2, q["U"], q["V"], 0.8, launch="inertial")
    assert np.allclose(qi["x_rotating"], q["x_rotating"], rtol=0, atol=1e-14)
    with pytest.raises(ValueError):
        F.turntable_puck(t, 0, 0, 1, 0, 1.0, launch="space")


def test_turntable_puck_matches_kundu_projectile_paths():  # cross-check with the Kundu-tested core function
    t = np.linspace(0, 6, 61)
    xi, xr = core_rot.projectile_paths(2.0, 0.4, t)
    r = F.turntable_puck(t, 0, 0, 2.0, 0, 0.4, launch="inertial")
    assert np.allclose(r["x_rotating"], xr[0], atol=1e-14) and np.allclose(r["y_rotating"], xr[1], atol=1e-14)


# ------------------------------------------------------------------------------------------------- f-plane inertial
def test_inertial_circle_V2_symbolic_solution_radius_and_period():  # V2: dsolve-free residual + centripetal balance
    t, f, u0, v0 = sp.symbols("t f u_0 v_0", real=True, nonzero=True)
    w = (u0 + sp.I * v0) * sp.exp(-sp.I * f * t)  # u + i v
    u, v = sp.re(sp.expand_complex(w)), sp.im(sp.expand_complex(w))
    assert sp.simplify(sp.diff(u, t) - f * v) == 0 and sp.simplify(sp.diff(v, t) + f * u) == 0  # Vallis (3.126)
    z = sp.I * (w - (u0 + sp.I * v0)) / f  # displacement x + i y
    zc = sp.I * (-(u0 + sp.I * v0)) / f  # centre
    assert sp.simplify(sp.Abs(z - zc) ** 2 - (u0**2 + v0**2) / f**2) == 0  # radius U/|f|
    assert sp.simplify(w.subs(t, 2 * sp.pi / f) - (u0 + sp.I * v0)) == 0  # period 2 pi/f
    U, r = sp.symbols("U r", positive=True)  # Coriolis f U balances centripetal U^2/r
    assert sp.solve(sp.Eq(f * U, U**2 / r), r) == [U / f]


def test_inertial_circle_V1_numbers_at_45N():  # V1 closed form
    f = F.coriolis_parameter(LAT45)
    ic = F.inertial_circle(10.0, f)
    assert ic["radius"] == pytest.approx(10.0 / f, rel=1e-15)
    assert ic["period"] == pytest.approx(2 * np.pi / f, rel=1e-15)
    assert ic["period"] == pytest.approx(C.SIDEREAL_DAY / (2 * np.sin(LAT45)), rel=1e-14)  # 1 day / (2 sin lat)
    assert ic["radius"] == pytest.approx(96968.68, abs=0.01)


def test_inertial_circle_V5_period_169h_at_45deg():  # V5 cited: ~16.9 h at 45 deg (0.1 h precision)
    T_h = F.inertial_circle(1.0, F.coriolis_parameter(LAT45))["period"] / 3600.0
    assert abs(T_h - 16.9) <= 0.05


def test_inertial_circle_V7_sense_flips_and_equator_limit():  # V7
    f = F.coriolis_parameter(LAT45)
    n, s = F.inertial_circle(10.0, f), F.inertial_circle(10.0, -f)
    assert n["sense"] == "clockwise seen from above" and n["turns"] == "right"
    assert s["sense"] == "anticlockwise seen from above" and s["turns"] == "left"
    assert n["radius"] == s["radius"] and n["period"] == s["period"]
    with pytest.raises(ValueError):
        F.inertial_circle(10.0, 0.0)  # equator: no inertial circle
    lats = np.deg2rad([60.0, 30.0, 10.0, 1.0, 0.1])
    T = [F.inertial_circle(10.0, F.coriolis_parameter(L))["period"] for L in lats]
    assert np.all(np.diff(T) > 0) and T[-1] > 100 * T[0]  # period grows without bound towards the equator


def test_inertial_oscillation_V1_matches_ode_and_V4_speed():  # V1 independent ODE; V4 speed conserved
    f = -1.1e-4
    t = np.linspace(0, 2 * 86400.0, 289)
    io = F.inertial_oscillation(t, 3.0, -7.0, f, x0=1000.0, y0=-500.0)

    def rhs(_, s):
        return [s[2], s[3], f * s[3], -f * s[2]]  # Vallis (3.126)-(3.127)
    y = solve_ivp(rhs, (0, t[-1]), [1000.0, -500.0, 3.0, -7.0], t_eval=t, method="DOP853", rtol=1e-12, atol=1e-10).y
    assert np.max(np.abs(io["x"] - y[0])) < 1e-5 and np.max(np.abs(io["y"] - y[1])) < 1e-5
    assert np.max(np.abs(io["u"] - y[2])) < 1e-9 and np.max(np.abs(io["v"] - y[3])) < 1e-9
    assert np.allclose(np.hypot(io["u"], io["v"]), np.hypot(3.0, -7.0), rtol=1e-14)
    assert np.allclose(np.hypot(io["x"] - io["xc"], io["y"] - io["yc"]), io["radius"], rtol=1e-12)


def test_inertial_oscillation_V7_closes_after_period_and_f0_straight():  # V7 closure, NH centre right, f = 0 line
    f = F.coriolis_parameter(LAT45)
    io = F.inertial_oscillation(np.array([0.0, 2 * np.pi / f]), 10.0, 0.0, f)
    assert np.allclose(io["x"], 0.0, atol=1e-8) and np.allclose(io["y"], 0.0, atol=1e-8)
    assert io["yc"] < 0 and io["xc"] == 0.0  # eastward start: centre to the south (right) in the NH
    st = F.inertial_oscillation(np.array([0.0, 100.0]), 2.0, 1.0, 0.0)
    assert np.allclose([st["x"][1], st["y"][1]], [200.0, 100.0]) and "radius" not in st


# ------------------------------------------------------------------------------------------ Coriolis-only on the sphere
def _sphere_ref(lat0, lon0, u0, v0, T, Om=OM, metric=True, coriolis=True):
    """Vallis (2.50a,b) with grad p = 0 and (2.28), written out here; DOP853 at rtol 1e-13."""
    def rhs(_, s):
        lam, th, u, v = s
        m = (2 * Om * np.sin(th) if coriolis else 0.0) + (u * np.tan(th) / A_E if metric else 0.0)
        return [u / (A_E * np.cos(th)), v / A_E, m * v, -m * u]
    return solve_ivp(rhs, (0, T), [lon0, lat0, u0, v0], method="DOP853", rtol=1e-13, atol=1e-15).y[:, -1]


def test_inertial_oscillation_sphere_V3_rk4_fourth_order():  # V3: classical RK4, design order 4
    T, dts, errs = 4 * 86400.0, [3600.0, 1800.0, 900.0, 450.0], []
    ref = _sphere_ref(LAT45, 0.0, 30.0, 10.0, T)
    for dt in dts:
        s = np.array(F.inertial_oscillation_sphere(LAT45, 0.0, 30.0, 10.0, T, dt=dt))
        e = s - ref
        e[0] = _wrap(e[0])
        errs.append(np.abs(e[:2]).max() * A_E)  # position error [m]
    p = observed_order(dts, errs)
    print("RK4 position errors [m]", errs, "pairwise", pairwise_orders(dts, errs), "order", p)
    assert abs(p - 4.0) < 0.15


def test_inertial_oscillation_sphere_V1_default_step_matches_reference():  # V1 vs the independent ODE solve
    T = 2 * 86400.0
    s = np.array(F.inertial_oscillation_sphere(np.deg2rad(-30.0), 1.0, -5.0, 12.0, T))
    ref = _sphere_ref(np.deg2rad(-30.0), 1.0, -5.0, 12.0, T)
    assert abs(_wrap(s[0] - ref[0])) * A_E < 1e-3 and abs(s[1] - ref[1]) * A_E < 1e-3  # < 1 mm
    assert np.max(np.abs(s[2:] - ref[2:])) < 1e-9
    # an off-grid output time takes one partial step; still on the reference
    s2 = np.array(F.inertial_oscillation_sphere(np.deg2rad(-30.0), 1.0, -5.0, 12.0, 1000.0 + 17.3))
    ref2 = _sphere_ref(np.deg2rad(-30.0), 1.0, -5.0, 12.0, 1017.3)
    assert np.max(np.abs(s2[:2] - ref2[:2])) * A_E < 1e-6


def test_inertial_oscillation_sphere_V4_speed_and_angular_momentum():  # V4: |u| and M = a cos(lat)(u + Omega a cos lat)
    t = np.arange(0.0, 5 * 86400.0 + 1, 600.0)
    lon, lat, u, v = F.inertial_oscillation_sphere(np.deg2rad(40.0), 0.0, 25.0, -15.0, t)
    sp0 = np.hypot(25.0, 15.0)
    assert np.max(np.abs(np.hypot(u, v) / sp0 - 1)) < 1e-10
    M = A_E * np.cos(lat) * (u + OM * A_E * np.cos(lat))
    M0 = A_E * np.cos(np.deg2rad(40.0)) * (25.0 + OM * A_E * np.cos(np.deg2rad(40.0)))
    assert np.max(np.abs(M / M0 - 1)) < 1e-10
    # the invariant really discriminates: the opposite metric sign breaks it at the 1e-6 level
    lonx, latx, ux, vx = (np.array(z) for z in F.inertial_oscillation_sphere(np.deg2rad(40.0), 0.0, 25.0, -15.0, t,
                                                                              metric=False))
    Mx = A_E * np.cos(latx) * (ux + OM * A_E * np.cos(latx))
    assert np.max(np.abs(Mx / M0 - 1)) > 1e-6


def test_inertial_oscillation_sphere_V7_f0_gives_great_circle():  # V7 f = 0 -> closed-form great circle
    t = np.linspace(0, 3 * 86400.0, 73)
    lg, pg, ug, vg = F.inertial_oscillation_sphere(LAT45, 0.3, 40.0, 25.0, t, coriolis=False)
    lc, pc = F.great_circle(LAT45, 0.3, 40.0, 25.0, t)
    assert np.max(np.abs(pg - pc)) * A_E < 1e-3 and np.max(np.abs(_wrap(lg - lc))) * A_E < 1e-3


def test_inertial_oscillation_sphere_V7_small_speed_f_plane_limit():  # V7: U -> 0 recovers the f-plane circle
    f = F.coriolis_parameter(LAT45)
    T = 2 * np.pi / f
    t = np.linspace(0, T, 97)
    rel, eps = [], []
    for U in (2.0, 1.0, 0.5, 0.25):
        lon, lat, u, v = F.inertial_oscillation_sphere(LAT45, 0.0, U, 0.0, t, dt=30.0)
        x = lon * A_E * np.cos(LAT45)
        y = (lat - LAT45) * A_E
        io = F.inertial_oscillation(t, U, 0.0, f)
        rel.append(np.max(np.hypot(x - io["x"], y - io["y"])) / (U / f))
        eps.append(U / (f * A_E))  # the small parameter r/a of the f-plane approximation
    k = np.array(rel) / np.array(eps)
    print("relative departure from the f-plane circle", rel, "rel/(U/(f a))", k)
    assert rel[-1] < 4e-3  # 0.36 % of the radius at U = 0.25 m/s
    assert observed_order(eps, rel) == pytest.approx(1.0, abs=0.05)  # departure is first order in r/a ...
    assert abs(k[-1] / k[-2] - 1) < 2e-3  # ... with a converged coefficient


def test_inertial_oscillation_sphere_V7_hemisphere_mirror_and_equator():  # V7
    t = np.arange(0, 2 * 86400.0 + 1, 900.0)
    ln, pn, un, vn = F.inertial_oscillation_sphere(LAT45, 0.0, 10.0, 3.0, t)
    ls, ps, us, vs = F.inertial_oscillation_sphere(-LAT45, 0.0, 10.0, -3.0, t)
    assert np.allclose(ps, -pn, atol=1e-15) and np.allclose(ls, ln, atol=1e-15)
    assert np.allclose(us, un, atol=1e-12) and np.allclose(vs, -vn, atol=1e-12)
    # first quarter-turn: NH turns right (clockwise), SH left
    assert vn[4] < vn[0] and vs[4] > vs[0]
    # equator: f = 0 and tan(0) = 0 there, so an eastward parcel stays on the equator
    le, pe, ue, ve = F.inertial_oscillation_sphere(0.0, 0.0, 10.0, 0.0, t)
    assert np.max(np.abs(pe)) == 0.0 and np.max(np.abs(ve)) == 0.0


def test_inertial_oscillation_sphere_V7_time_reversal_and_beta_drift():  # V7 backward/forward; westward beta drift
    lb, pb, ub, vb = F.inertial_oscillation_sphere(LAT45, 0.0, 10.0, 0.0, np.array([-7200.0, 0.0]))
    lf, pf, uf, vf = F.inertial_oscillation_sphere(pb[0], lb[0], ub[0], vb[0], 7200.0)
    assert abs(_wrap(lf)) * A_E < 1e-6 and abs(pf - LAT45) * A_E < 1e-6 and abs(uf - 10.0) < 1e-10
    assert pb[1] == LAT45 and lb[1] == 0.0
    f = F.coriolis_parameter(LAT45)
    lon, lat, u, v = F.inertial_oscillation_sphere(LAT45, 0.0, 10.0, 0.0, 3 * 2 * np.pi / f)
    assert lon < 0  # inertial loops on the sphere drift westward (beta effect)


def test_inertial_oscillation_sphere_guards():
    t = np.linspace(0, 86400.0 * 10, 11)
    lon, lat, u, v = F.inertial_oscillation_sphere(np.deg2rad(88.5), 0.0, 0.0, 200.0, t)
    assert np.isnan(lat[-1]) and not np.isnan(lat[0])  # past lat_max_deg the (lon, lat) form stops
    with pytest.raises(ValueError):
        F.inertial_oscillation_sphere(0.5, 0.0, 1.0, 0.0, 10.0, dt=0.0)
    lon, lat, u, v = F.inertial_oscillation_sphere(0.5, 0.0, 1.0, 0.0, 600.0)
    assert all(isinstance(z, float) for z in (lon, lat, u, v))


# ------------------------------------------------------------------------------------ space view, great circles
def test_to_inertial_lon_V1_offset_and_V7_inverse():  # V1 lambda_I = lambda_E + Omega t; V7 inverse with -Omega
    lon = RNG.uniform(-np.pi, np.pi, 50)
    t = RNG.uniform(0, 5 * 86400, 50)
    li = F.to_inertial_lon(lon, t)
    assert np.max(np.abs(_wrap(li - lon - OM * t))) < 1e-12
    assert np.max(np.abs(_wrap(F.to_inertial_lon(li, t, -OM) - lon))) < 1e-12
    assert np.all((li > -np.pi) & (li <= np.pi))
    assert F.to_inertial_lon(0.3, C.SIDEREAL_DAY) == pytest.approx(0.3, abs=1e-12)  # one sidereal day = one turn


def test_to_inertial_lon_V1_agrees_with_3d_rotation():  # V1: the Cartesian R_z route in great_circle_earth_frame
    t = np.linspace(0, 6 * 3600.0, 25)
    g = F.great_circle_earth_frame(t, 0.6, -0.4, 1.0, 150.0)
    assert np.max(np.abs(_wrap(F.to_inertial_lon(g["lon"], t) - g["lon_inertial"]))) < 1e-12


def test_great_circle_V1_plane_through_centre_and_arc_length():  # V1 3-D geometry, independent of the library
    lat0, lon0, u0, v0 = 0.7, -1.2, 120.0, -60.0
    U = np.hypot(u0, v0)
    t = np.linspace(0, 0.9 * np.pi * A_E / U, 200)
    lon, lat = F.great_circle(lat0, lon0, u0, v0, t)
    r = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)
    r0 = np.array([np.cos(lat0) * np.cos(lon0), np.cos(lat0) * np.sin(lon0), np.sin(lat0)])
    east = np.array([-np.sin(lon0), np.cos(lon0), 0.0])
    north = np.cross(r0, east)
    n = np.cross(r0, (u0 * east + v0 * north) / U)  # normal of the great-circle plane
    assert np.max(np.abs(r @ n)) < 1e-13
    arc = np.arccos(np.clip(r @ r0, -1, 1))
    assert np.max(np.abs(arc - U * t / A_E)) < 1e-7  # distance from start = U t / a
    lp, pp = F.great_circle(lat0, lon0, u0, v0, 2 * np.pi * A_E / U)  # one lap closes
    assert abs(pp - lat0) < 1e-12 and abs(_wrap(lp - lon0)) < 1e-12
    l0, p0 = F.great_circle(lat0, lon0, 0.0, 0.0, t)
    assert np.all(p0 == lat0)


def test_great_circle_earth_frame_V1_rotation_and_launch_speed():  # V1 |r| = a, same lat, lon offset -Omega t
    t = np.linspace(0, 3 * 3600.0, 7)
    g = F.great_circle_earth_frame(t, LAT45, 0.0, np.pi / 2, 300.0)
    assert np.allclose(np.linalg.norm(g["xyz_earth"], axis=-1), A_E, rtol=1e-15)
    assert np.array_equal(g["lat"], g["lat_inertial"])
    assert np.max(np.abs(_wrap(g["lon_inertial"] - g["lon"] - OM * t))) < 1e-12
    assert g["speed_inertial"] == pytest.approx(300.0 + OM * A_E * np.cos(LAT45), rel=1e-14)  # Vallis (2.7)
    g0 = F.great_circle_earth_frame(t, LAT45, 0.0, np.pi / 2, 300.0, Omega=0.0)
    lc, pc = F.great_circle(LAT45, 0.0, 300.0, 0.0, t)  # V7 Omega -> 0
    assert np.allclose(g0["lat"], pc, atol=1e-14) and np.allclose(g0["lon"], lc, atol=1e-14)
    with pytest.raises(ValueError):
        F.great_circle_earth_frame(t, 0.0, 0.0, 0.0, 0.0, Omega=0.0)
