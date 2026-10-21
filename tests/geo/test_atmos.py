"""Evidence tests for ``fluidpy.geo.atmos`` and the gas constants of ``fluidpy.geo.constants`` (AOFD PoC F: geostrophic
wind from 500 hPa heights on the sphere; scale height; hypsometric thickness).

Benchmarks (cited, verified 2026-10-21; see AOFD ``reference/ch00/SOURCES.md``):
* U.S. Standard Atmosphere 1976 (NOAA/NASA/USAF, NASA-TM-X-74335), as tabulated by Public Domain Aeronautical Software
  (pdas.com "bigtables", Tables 1-2): pressure scale height 8435 m at sea level and 6382 m at 20 km geometric;
  p = 5.4048e4 Pa and T = 255.676 K at 5 km geometric (geopotential 4996.07 m with r0 = 6 356 766 m).
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import solve_ivp

from fluidpy.core import statics as core_statics
from fluidpy.core import thermo as core_thermo
from fluidpy.geo import atmos as AT
from fluidpy.geo import constants as C
from tools.convergence import observed_order, pairwise_orders

G0, OM, A_E, R_D = C.G0, C.OMEGA, C.EARTH_RADIUS, C.R_D
R0_USSA = 6356766.0  # USSA-1976 effective Earth radius for geopotential altitude [m]


# ---------------------------------------------------------------------------------------------------------- constants
def test_constants_V1_gas_constants_one_home():  # V1 identity with the core values; kappa = 2/7 for gamma = 1.4
    assert C.G0 == 9.80665 and C.G0 == core_thermo.G0
    assert C.R_D == core_thermo.R_AIR and C.R_D == pytest.approx(287.05, abs=0.01)
    assert C.C_P == core_thermo.CP_AIR and C.C_P == pytest.approx(3.5 * C.R_D, rel=1e-15)
    assert C.KAPPA == pytest.approx(2.0 / 7.0, rel=1e-15)
    assert C.P_REF == 1.0e5 and C.RHO0_OCEAN == 1025.0


# --------------------------------------------------------------------------------------------- geostrophic wind from Z
def _sym_geostrophic(Zexpr, th, lam):
    """Vallis (2.198)/(2.215) on the sphere: f u = -(g/a) dZ/dlat, f v = g/(a cos lat) dZ/dlon, f = 2 Omega sin lat."""
    g, Om, a = sp.Rational(980665, 100000), sp.Float(OM, 30), sp.Float(A_E, 30)
    f = 2 * Om * sp.sin(th)
    return -(g / (f * a)) * sp.diff(Zexpr, th), g / (f * a * sp.cos(th)) * sp.diff(Zexpr, lam)


def test_geostrophic_V2_symbolic_solid_body_and_wave_fields():  # V2: exact winds of two analytic height fields
    th, lam, u0, B, Om, a, g, Z0 = sp.symbols("vartheta lambda u_0 B Omega a g Z_0", positive=True)
    f = 2 * Om * sp.sin(th)
    # (i) solid-body rotation u = u0 cos(lat): Z = Z0 - (Omega a u0/g) sin^2(lat)
    Z = Z0 - Om * a * u0 / g * sp.sin(th) ** 2
    assert sp.simplify(-(g / (f * a)) * sp.diff(Z, th) - u0 * sp.cos(th)) == 0
    assert sp.diff(Z, lam) == 0
    # (ii) wave field Z = B sin(lat) cos(lat) sin(2 lon): v_g = g B 2 cos(2 lon)/(2 Omega a), uniform in latitude
    Zw = B * sp.sin(th) * sp.cos(th) * sp.sin(2 * lam)
    vg = g / (f * a * sp.cos(th)) * sp.diff(Zw, lam)
    assert sp.simplify(vg - g * B * sp.cos(2 * lam) / (Om * a)) == 0
    # (iii) the spherical form is f k x u_g = -grad Phi: check the vector identity componentwise
    ug = -(g / (f * a)) * sp.diff(Zw, th)
    grad = sp.Matrix([g * sp.diff(Zw, lam) / (a * sp.cos(th)), g * sp.diff(Zw, th) / a])
    kxu = sp.Matrix([-vg, ug])
    assert sp.simplify(f * kxu + grad) == sp.zeros(2, 1)


def test_geostrophic_V1_exact_on_quadratic_regional_field():  # V1: centred/one-sided 2nd-order stencils are exact here
    lat = np.arange(20.0, 70.01, 2.0)
    lon = np.arange(100.0, 160.01, 2.0)
    th, lam = np.deg2rad(lat)[:, None], np.deg2rad(lon)[None, :]
    al, be, ga = 300.0, -150.0, 80.0  # [m per rad^2]
    Z = 5500.0 + al * th**2 + be * lam**2 + ga * th * lam
    ug, vg = AT.geostrophic_wind_from_height(Z, lat, lon)
    f = 2 * OM * np.sin(th)
    ue = -(G0 / (f * A_E)) * (2 * al * th + ga * lam)
    ve = G0 / (f * A_E * np.cos(th)) * (2 * be * lam + ga * th)
    assert np.max(np.abs(ug - ue)) < 1e-11 and np.max(np.abs(vg - ve)) < 1e-11  # includes every edge node


@pytest.mark.parametrize("field", ["solid_body", "wave"])
def test_geostrophic_V3_second_order_on_global_grid(field):  # V3: design order 2 (centred differences)
    th, lam = sp.symbols("vartheta lambda", real=True)
    if field == "solid_body":
        Zs = 5500 - sp.Float(OM * A_E * 20.0 / G0, 30) * sp.sin(th) ** 2
    else:
        Zs = 5500 + 150 * sp.sin(th) * sp.cos(th) * sp.sin(2 * lam) + 80 * sp.cos(th) ** 2 * sp.cos(3 * lam)
    use, vse = _sym_geostrophic(Zs, th, lam)
    fZ, fu, fv = (sp.lambdify((th, lam), e, "numpy") for e in (Zs, use, vse))
    hs, eu, ev = [], [], []
    for d in (4.0, 2.0, 1.0, 0.5):
        lat = np.arange(-90.0, 90.0 + d / 2, d)
        lon = np.arange(0.0, 360.0, d)
        TH, LA = np.meshgrid(np.deg2rad(lat), np.deg2rad(lon), indexing="ij")
        Z = fZ(TH, LA) + 0 * TH
        ug, vg = AT.geostrophic_wind_from_height(Z, lat, lon)
        sel = (np.abs(lat) >= 20) & (np.abs(lat) <= 70)
        with np.errstate(divide="ignore", invalid="ignore"):  # the exact winds are singular at the poles/equator
            ue, ve = fu(TH, LA) + 0 * TH, fv(TH, LA) + 0 * TH
        hs.append(np.deg2rad(d))
        eu.append(np.max(np.abs(ug - ue)[sel]))
        ev.append(np.max(np.abs(vg - ve)[sel]))
    print(field, "u errors", eu, pairwise_orders(hs, eu), "v errors", ev)
    assert abs(observed_order(hs, eu) - 2.0) < 0.15
    if field == "wave":
        assert abs(observed_order(hs, ev) - 2.0) < 0.15
    else:
        assert max(ev) == 0.0  # zonally uniform field: no meridional wind at all


def _gauss_low(lat, lon, lat_c, lon_c, depth=100.0, w=0.05):
    th, lam = np.deg2rad(lat)[:, None], np.deg2rad(lon)[None, :]
    return 5500.0 - depth * np.exp(-((th - np.deg2rad(lat_c)) ** 2 + (lam - np.deg2rad(lon_c)) ** 2) / w)


def test_geostrophic_V7_cyclonic_round_lows_both_hemispheres():  # V7: anticlockwise round a NH low, clockwise SH
    lat = np.arange(-90.0, 90.01, 1.0)
    lon = np.arange(0.0, 360.0, 1.0)
    for lat_c, sense in ((45.0, +1), (-45.0, -1)):
        ug, vg = AT.geostrophic_wind_from_height(_gauss_low(lat, lon, lat_c, 180.0), lat, lon)
        j0, i0 = int(np.argmin(np.abs(lat - lat_c))), int(np.argmin(np.abs(lon - 180.0)))
        # circulation round a ring of radius 5 deg: sum of tangential components (anticlockwise positive)
        tang = []
        for ang in np.linspace(0, 2 * np.pi, 16, endpoint=False):
            dj, di = int(round(5 * np.sin(ang))), int(round(5 * np.cos(ang)))
            tang.append(-np.sin(ang) * ug[j0 + dj, i0 + di] + np.cos(ang) * vg[j0 + dj, i0 + di])
        assert np.sign(np.sum(tang)) == sense
        assert np.sign(vg[j0, i0 + 5]) == sense  # NH: southerly (northward) wind east of the low
        assert np.sign(vg[j0, i0 - 5]) == -sense  # ... northerly west of it
        assert np.sign(ug[j0 + 5, i0]) == -sense and np.sign(ug[j0 - 5, i0]) == sense  # easterly poleward side (NH)


def test_geostrophic_V7_mirror_gauge_ordering_and_broadcast():  # V7 symmetries of the operator
    lat = np.arange(-90.0, 90.01, 2.5)
    lon = np.arange(0.0, 360.0, 2.5)
    Z = _gauss_low(lat, lon, 50.0, 200.0) + _gauss_low(lat, lon, 30.0, 60.0, depth=-60.0) - 5500.0
    ug, vg = AT.geostrophic_wind_from_height(Z, lat, lon)
    um, vm = AT.geostrophic_wind_from_height(Z[::-1], lat, lon)  # field mirrored about the equator
    assert np.allclose(um[::-1], ug, equal_nan=True, atol=1e-12) and np.allclose(vm[::-1], -vg, equal_nan=True, atol=1e-12)
    uc, vc = AT.geostrophic_wind_from_height(Z + 1234.5, lat, lon)  # adding a constant height changes nothing
    assert np.allclose(uc, ug, equal_nan=True, atol=1e-10) and np.allclose(vc, vg, equal_nan=True, atol=1e-10)
    un, vn = AT.geostrophic_wind_from_height(Z[::-1], lat[::-1], lon)  # PSL files run N -> S
    assert np.allclose(un[::-1], ug, equal_nan=True, atol=1e-12) and np.allclose(vn[::-1], vg, equal_nan=True, atol=1e-12)
    u3, v3 = AT.geostrophic_wind_from_height(np.stack([Z, 2 * Z]), lat, lon)  # leading (time) axis
    assert np.allclose(u3[1], 2 * ug, equal_nan=True, atol=1e-12) and u3.shape == (2,) + Z.shape
    # regional grid interior = global interior
    ur, vr = AT.geostrophic_wind_from_height(Z[:, 20:60], lat, lon[20:60])
    assert np.allclose(ur[:, 1:-1], ug[:, 21:59], equal_nan=True, atol=1e-12)
    # regional d/dlon uses the coordinate differences, global uses 2 dlon: round-off only (|v| ~ 8 m/s here)
    assert np.allclose(vr[:, 1:-1], vg[:, 21:59], equal_nan=True, rtol=0, atol=1e-12 * np.nanmax(np.abs(vg)) * 10)


def test_geostrophic_V7_equator_and_pole_masks():  # V7 f -> 0 at the equator: geostrophy undefined -> NaN
    lat = np.arange(-90.0, 90.01, 2.5)
    lon = np.arange(0.0, 360.0, 2.5)
    Z = _gauss_low(lat, lon, 45.0, 180.0)
    for lm in (5.0, 15.0):
        ug, vg = AT.geostrophic_wind_from_height(Z, lat, lon, lat_min=lm)
        bad = (np.abs(lat) < lm) | (np.abs(lat) == 90.0)
        assert np.all(np.isnan(ug[bad])) and np.all(np.isnan(vg[bad]))
        assert np.all(np.isfinite(ug[~bad])) and np.all(np.isfinite(vg[~bad]))
    with pytest.raises(ValueError):
        AT.geostrophic_wind_from_height(Z[:, :-1], lat, lon)
    with pytest.raises(ValueError):
        AT.geostrophic_wind_from_height(Z, lat, lon[::-1])


def test_geostrophic_V2_dimensions():  # V2 dimensional homogeneity of (g/(f a)) dZ/dlat
    from fluidpy.core.units import Q_
    u = Q_(G0, "m/s**2") / (Q_(1e-4, "1/s") * Q_(A_E, "m")) * Q_(100.0, "m")
    assert u.check("[velocity]")


# -------------------------------------------------------------------------------------------------------- scale height
def test_scale_height_V1_closed_form_and_isothermal_efolding():  # V1 RT/g; independent hydrostatic ODE e-folding
    T = np.linspace(180.0, 320.0, 15)
    assert np.allclose(AT.scale_height(T), R_D * T / G0, rtol=1e-15)
    assert np.allclose(AT.scale_height(T), core_statics.scale_height(T), rtol=1e-15)
    for Tk in (220.0, 288.15):
        H = AT.scale_height(Tk)
        sol = solve_ivp(lambda z, p: -p * G0 / (R_D * Tk), (0, H), [1.0e5], rtol=1e-12, atol=1e-6)
        assert sol.y[0, -1] / 1.0e5 == pytest.approx(np.exp(-1.0), rel=1e-9)  # dp/dz = -rho g, rho = p/(R T)


def test_scale_height_V5_us_standard_atmosphere_1976():  # V5 USSA-1976 (pdas.com Table 2): 8435 m, 6382 m at 20 km
    H0 = AT.scale_height(288.15)
    z = 20000.0
    g20 = G0 * (R0_USSA / (R0_USSA + z)) ** 2  # USSA uses local gravity in its scale height
    H20 = AT.scale_height(216.65, g=g20)
    # printed to 4 significant figures (+-0.5 m); our R_d = 287.058 vs USSA R*/M0 = 287.053 adds 1.7e-5 relative
    assert abs(H0 - 8435.0) <= 0.5 + 1.7e-5 * 8435
    assert abs(H20 - 6382.0) <= 0.5 + 1.7e-5 * 6382
    assert 7000.0 < AT.scale_height(240.0) < 8500.0 < AT.scale_height(300.0)  # "7-8 km" for tropospheric T


def test_scale_height_V2_dimensions():  # V2: R T / g is a length
    from fluidpy.core.units import Q_
    assert (Q_(R_D, "J/(kg*K)") * Q_(250.0, "K") / Q_(G0, "m/s**2")).check("[length]")


# ----------------------------------------------------------------------------------------------- hypsometric thickness
def _ussa_trop_z(p):
    """USSA-1976 troposphere written out here: T = T0 - L H (geopotential H), p = p0 (T/T0)^(g0/(R L))."""
    T0, L, p0, R = 288.15, 0.0065, 101325.0, 8314.32 / 28.9644
    return T0 / L * (1 - (p / p0) ** (R * L / G0)), T0 * (p / p0) ** (R * L / G0), R


def test_hypsometric_V1_isothermal_closed_form_and_profile_route():  # V1 (2.166) at constant T, both code paths
    Tm = 250.0
    assert AT.hypsometric_thickness(1e5, 5e4, Tm) == pytest.approx(R_D * Tm / G0 * np.log(2.0), rel=1e-15)
    pp = np.geomspace(1.05e5, 9e3, 50)
    assert AT.hypsometric_thickness(1e5, 1e4, Tm + 0 * pp, p=pp) == pytest.approx(R_D * Tm / G0 * np.log(10), rel=1e-13)
    # thickness / scale height = ln(p1/p2): one scale height per factor e
    assert AT.hypsometric_thickness(1e5, 1e5 / np.e, Tm) == pytest.approx(AT.scale_height(Tm), rel=1e-15)


def test_hypsometric_V3_profile_integral_second_order():  # V3: trapezoid in ln p, T nonlinear in ln p (USSA troposph.)
    p1, p2 = 100000.0, 30000.0
    z1, _, R = _ussa_trop_z(p1)
    z2, _, _ = _ussa_trop_z(p2)
    exact = z2 - z1
    hs, errs = [], []
    for n in (9, 17, 33, 65, 129):
        lnp = np.linspace(np.log(p2), np.log(p1), n)
        pp = np.exp(lnp)
        _, T, _ = _ussa_trop_z(pp)
        errs.append(abs(AT.hypsometric_thickness(p1, p2, T, p=pp, R=R) - exact))
        hs.append(lnp[1] - lnp[0])
    print("hypsometric errors [m]", errs, pairwise_orders(hs, errs))
    assert abs(observed_order(hs, errs) - 2.0) < 0.15


def test_hypsometric_V5_us_standard_atmosphere_5km():  # V5 USSA-1976: p(5 km geometric) = 54048 Pa
    pp = np.geomspace(101325.0, 50000.0, 2001)
    _, T, R = _ussa_trop_z(pp)
    dz = AT.hypsometric_thickness(101325.0, 54048.0, T, p=pp, R=R)
    H5 = R0_USSA * 5000.0 / (R0_USSA + 5000.0)  # geopotential altitude of 5 km geometric = 4996.07 m
    assert abs(dz - H5) < 0.1  # p printed to 5 figures: +-0.5 Pa is +-0.07 m
    # and the temperature the closed-form profile gives at that pressure is the tabulated 255.676 K
    assert _ussa_trop_z(54048.0)[1] == pytest.approx(255.676, abs=2e-3)


def test_hypsometric_V1_exact_for_T_linear_in_lnp_off_grid_ends():  # V1: ln-p interpolation + trapezoid exact
    pp = np.geomspace(1.0e5, 1.0e4, 7)  # coarse; the bounds fall inside intervals
    c0, c1 = 250.0, 12.0
    T = c0 + c1 * np.log(pp / 1e5)
    p1, p2 = 83000.0, 21000.0
    x1, x2 = np.log(p1 / 1e5), np.log(p2 / 1e5)
    exact = -(R_D / G0) * (c0 * (x2 - x1) + 0.5 * c1 * (x2**2 - x1**2))  # -(R/g) int T dln p, Vallis (2.166)
    assert AT.hypsometric_thickness(p1, p2, T, p=pp) == pytest.approx(exact, rel=1e-13)
    assert AT.hypsometric_thickness(p1, p2, T[::-1], p=pp[::-1]) == pytest.approx(exact, rel=1e-13)  # any order


def test_hypsometric_V7_sign_additivity_and_guards():  # V7 antisymmetric, additive, positive upward
    pp = np.geomspace(1.0e5, 1.0e4, 80)
    T = 290.0 - 60.0 * np.log(1e5 / pp) / np.log(10)
    up = AT.hypsometric_thickness(1e5, 3e4, T, p=pp)
    assert up > 0 and AT.hypsometric_thickness(3e4, 1e5, T, p=pp) == pytest.approx(-up, rel=1e-15)
    a1 = AT.hypsometric_thickness(1e5, 5e4, T, p=pp) + AT.hypsometric_thickness(5e4, 3e4, T, p=pp)
    assert a1 == pytest.approx(up, rel=1e-3)  # additive up to the trapezoid interpolation at the split level
    assert AT.hypsometric_thickness(1e5, 1e5, 250.0) == 0.0
    assert np.allclose(AT.hypsometric_thickness(1e5, 5e4, np.array([220.0, 260.0, 300.0])),
                       R_D * np.array([220.0, 260.0, 300.0]) / G0 * np.log(2))  # warmer layer, thicker
    with pytest.raises(ValueError):
        AT.hypsometric_thickness(1.2e5, 5e4, T, p=pp)
    with pytest.raises(ValueError):
        AT.hypsometric_thickness(1e5, 5e4, T[:-1], p=pp)
