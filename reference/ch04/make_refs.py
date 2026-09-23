"""Write the cited public benchmark values used by ``tests/test_ch04.py`` (V5 evidence) into ``reference/ch04/``.

Chapter 4 (conservation laws) needs no datasets: almost every check is Tier 1 (exact solutions of the Navier–Stokes
equations, sympy re-derivations, conservation budgets, convergence orders). The public sources below pin a handful of
published numbers and forms the chapter's examples touch: the sphere-drag correlation behind our Fig. 4.21 curve
(Morrison 2016), the WGS-84 ellipsoid and rotation rate (Earth's oblateness, Coriolis and centrifugal numbers), the
sharp-orifice contraction coefficient, Bélanger's jump relation (the bore of Example 4.3 in its own frame), the capillary
length of water, Tsiolkovsky's rocket equation (Example 4.4 without gravity), Prandtl numbers of air, water and
monatomic gases, the added mass of a sphere, and the Taylor–Green vortex. Every entry was read from the cited page on
2026-09-23 (see ``SOURCES.md``). Nothing here comes from the textbook; book-quoted numbers live in the git-ignored
``tests/book_values_ch04.json``. Sea-level sound speed (USSA-1976) and the IAPWS surface tension of water are reused
from ``reference/ch01/``.

Run: ``.venv/Scripts/python.exe reference/ch04/make_refs.py``  (idempotent; overwrites ``benchmarks.json``).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_meta": {
        "chapter": "ch04",
        "written_by": "reference/ch04/make_refs.py",
        "verified": "2026-09-23",
        "note": "Public sources only (pages fetched on the date above, see SOURCES.md). No book values.",
    },
    # F. A. Morrison, "Data Correlation for Drag Coefficient for Sphere", Michigan Tech (10 Nov 2016), Eq. (1); from
    # Morrison, An Introduction to Fluid Mechanics, CUP 2013, Fig. 8.13 p. 625. PDF text extracted 2026-09-23:
    # "We do not recommend using equation 1 for Reynolds numbers larger than 10^6"; "At the highest Reynolds numbers,
    # equation 1 plateaus at C_D ≈ 0.14"; "for Re < 2 equation 1 follows the creeping-flow result (C_D = 24/Re)".
    "morrison_sphere_drag": {
        "formula": "24/Re + 2.6*(Re/5.0)/(1 + (Re/5.0)**1.52) + 0.411*(Re/2.63e5)**(-7.94)/(1 + (Re/2.63e5)**(-8.00))"
                   " + 0.25*(Re/1e6)/(1 + Re/1e6)",
        "Re_max": 1.0e6,
        "plateau_CD": 0.14,
        "creeping_limit": "24/Re (Re < 2)",
        "source": "Morrison (2016), Michigan Tech, Eq. (1)",
        "url": "https://pages.mtu.edu/~fmorriso/DataCorrelationForSphereDrag2016.pdf",
    },
    # Wikipedia, "World Geodetic System" (WGS 84 defining parameters).
    "wgs84": {
        "a_m": 6378137.0,
        "inv_f": 298.257223563,
        "b_m": 6356752.314245,
        "omega_rad_s": 72.92115e-6,
        "source": "https://en.wikipedia.org/wiki/World_Geodetic_System (fetched 2026-09-23)",
    },
    # Wikipedia, "Vena contracta": "The typical value may be taken as 0.611 for a sharp orifice (concentric with the
    # flow channel)."
    "vena_contracta_sharp_orifice": {
        "Cc": 0.611,
        "source": "https://en.wikipedia.org/wiki/Vena_contracta (fetched 2026-09-23)",
    },
    # Wikipedia, "Hydraulic jump": h2/h1 = (sqrt(1 + 8 Fr1^2) - 1)/2, Fr1 = v1/sqrt(g h1); "In a frame of reference
    # moving with a surge, the surge is equivalent to a stationary jump."
    "belanger_jump": {
        "h2_over_h1": "(sqrt(1 + 8*Fr1**2) - 1)/2",
        "Fr1": "v1/sqrt(g*h1)",
        "moving_surge": "equivalent to a stationary jump in the frame moving with the surge",
        "source": "https://en.wikipedia.org/wiki/Hydraulic_jump (fetched 2026-09-23)",
    },
    # Wikipedia, "Capillary length": lambda_c = sqrt(gamma/(Delta rho g)); water-air at 20 C: 2.71 mm (the page does
    # not state the surface tension or g used).
    "capillary_length_water_20C": {
        "formula": "sqrt(gamma/(Delta_rho*g))",
        "length_m": 2.71e-3,
        "digits": 3,
        "T_C": 20.0,
        "source": "https://en.wikipedia.org/wiki/Capillary_length (fetched 2026-09-23)",
    },
    # Wikipedia, "Tsiolkovsky rocket equation": Delta v = v_e ln(m0/m_f).
    "tsiolkovsky": {
        "formula": "v_e*log(m0/m_f)",
        "source": "https://en.wikipedia.org/wiki/Tsiolkovsky_rocket_equation (fetched 2026-09-23)",
    },
    # Wikipedia, "Prandtl number": air (250-1000 K) 0.70-0.73; water 10.6 (280 K), 5.9 (300 K), 2.0 (360 K);
    # monatomic gases 2/3 (kinetic theory); Eucken Pr = 4 gamma/(9 gamma - 5).
    "prandtl": {
        "air_range_250_1000K": [0.70, 0.73],
        "water_300K": 5.9,
        "water_280K": 10.6,
        "monatomic": "2/3",
        "eucken": "4*gamma/(9*gamma - 5)",
        "source": "https://en.wikipedia.org/wiki/Prandtl_number (fetched 2026-09-23)",
    },
    # Wikipedia, "Added mass": sphere added mass (2/3) pi r^3 rho_fluid = half the displaced mass.
    "added_mass_sphere": {
        "formula": "2/3*pi*r**3*rho",
        "fraction_of_displaced": 0.5,
        "source": "https://en.wikipedia.org/wiki/Added_mass (fetched 2026-09-23)",
    },
    # Wikipedia, "Taylor–Green vortex": u = U0 sin kx cos ky F, v = -U0 cos kx sin ky F, F = exp(-2 nu k^2 t),
    # p = (rho U0^2/4)(cos 2kx + cos 2ky) F^2.
    "taylor_green": {
        "u": "U0*sin(k*x)*cos(k*y)*exp(-2*nu*k**2*t)",
        "v": "-U0*cos(k*x)*sin(k*y)*exp(-2*nu*k**2*t)",
        "p": "rho*U0**2/4*(cos(2*k*x) + cos(2*k*y))*exp(-4*nu*k**2*t)",
        "source": "https://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex (fetched 2026-09-23)",
    },
}


def main() -> int:
    out = HERE / "benchmarks.json"
    out.write_text(json.dumps(BENCHMARKS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
