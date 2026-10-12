"""Write ``reference/ch14/benchmarks.json`` — the cited public values and forms the chapter-14 tests compare with.

Only statements read first-hand in a public source on the date given are stored (see SOURCES.md for the full citation
and what was read).  Nothing here comes from the textbook.  Most entries are *forms* (a formula as the source writes
it); the tests evaluate our functions against the form at our own inputs.  Numbers are stored exactly as the source
prints them.

Run: ``.venv/Scripts/python.exe reference/ch14/make_refs.py``  (idempotent: rewrites the same bytes).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_note": "Cited public values and forms for aerodynamics; our own numbers are computed, never typed. See SOURCES.md.",
    "thin_airfoil_lift_slope": {
        "form": "c_l = 2 pi alpha (thin symmetric airfoil, alpha in radians); c_l = c_l0 + 2 pi alpha (cambered)",
        "slope_over_pi": 2,
        "what": "section lift-curve slope per radian of a thin airfoil of infinite span",
        "source": "Wikipedia, 'Airfoil', section 'Thin airfoil theory'",
        "url": "https://en.wikipedia.org/wiki/Airfoil",
        "how_obtained": "read on the page (secondary source; the page cites Abbott & von Doenhoff and Clancy, not read)",
        "verified": "2026-10-12",
    },
    "thin_airfoil_centre_of_pressure_symmetric": {
        "value": 0.25,
        "digits": 2,
        "what": "centre of pressure of a thin symmetric airfoil, as a fraction of the chord behind the leading edge "
                "('exactly one quarter of the chord behind the leading edge')",
        "source": "Wikipedia, 'Airfoil', section 'Thin airfoil theory'",
        "url": "https://en.wikipedia.org/wiki/Airfoil",
        "how_obtained": "read on the page (secondary source)",
        "verified": "2026-10-12",
    },
    "joukowsky_kutta_circulation": {
        "form": "Gamma = 4 pi V_inf R sin(alpha + asin(mu_y / R)), R = sqrt((1 - mu_x)^2 + mu_y^2), for z = zeta + 1/zeta "
                "and a circle centred at mu = mu_x + i mu_y through zeta = 1",
        "what": "circulation fixed by the Kutta condition on a Joukowsky airfoil",
        "source": "Wikipedia, 'Joukowsky transform'",
        "url": "https://en.wikipedia.org/wiki/Joukowsky_transform",
        "how_obtained": "read on the page (secondary source)",
        "verified": "2026-10-12",
    },
    "induced_drag_coefficient": {
        "form": "C_D,i = C_L^2 / (pi AR e)",
        "e_elliptic": 1.0,
        "what": "induced drag coefficient; e = 1.0 for an elliptic distribution; AR = span^2 / wing area",
        "source": "NASA Glenn Research Center, Beginner's Guide to Aeronautics, 'Induced Drag Coefficient' (T. Benson)",
        "url": "https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/induced.html",
        "second_source": "Wikipedia, 'Lift-induced drag' (same form)",
        "second_url": "https://en.wikipedia.org/wiki/Lift-induced_drag",
        "how_obtained": "read on both pages",
        "verified": "2026-10-12",
    },
    "induced_drag_force": {
        "form": "D_i = L^2 / (1/2 rho V^2 pi b^2), b the wingspan",
        "what": "induced drag of a wing in terms of its lift and span (span efficiency 1)",
        "source": "Wikipedia, 'Lift-induced drag'",
        "url": "https://en.wikipedia.org/wiki/Lift-induced_drag",
        "how_obtained": "read on the page (secondary source)",
        "verified": "2026-10-12",
    },
    "induced_share_wide_body_cruise": {
        "value": 0.37,
        "digits": 2,
        "what": "induced drag as a fraction of the total drag of a typical twin-engine wide-body aircraft at cruise "
                "('approximately 37 %'; skin friction 'almost 48 %') - evidence for trap T9, not a test of our code",
        "source": "Wikipedia, 'Lift-induced drag'",
        "url": "https://en.wikipedia.org/wiki/Lift-induced_drag",
        "how_obtained": "read on the page (secondary source)",
        "verified": "2026-10-12",
    },
    "wagner_jones_fit": {
        "form": "phi(s) = 1 - C1 exp(-eps1 s) - C2 exp(-eps2 s), s = U t / (semi-chord)",
        "C1": 0.165,
        "C2": 0.335,
        "eps1": 0.0455,
        "eps2": 0.3,
        "phi_at_zero": 0.5,
        "max_abs_error_stated": 0.01,
        "what": "R. T. Jones' two-exponential approximation of Wagner's function (growth of circulatory lift after a step "
                "in incidence); time made non-dimensional with the semi-chord and the free-stream speed",
        "source": "S. T. M. Dawson & S. L. Brunton, 'Improved approximations to the Wagner function using sparse "
                  "identification of nonlinear dynamics', arXiv:2104.15122v2 (2021), Eq. (1) and Table 1; the absolute "
                  "error of the fit is 'less than 1%' (text under Fig. 2)",
        "url": "https://arxiv.org/abs/2104.15122",
        "second_source": "BYU FLOW Lab, Aeroelasticity.jl documentation, 'Wagner's Function' (C1 = 0.165, C2 = 0.335, "
                         "eps1 = 0.0455, eps2 = 0.3; exponent (u/b) t with b the semi-chord)",
        "second_url": "https://flow.byu.edu/Aeroelasticity.jl/dev/models/aerodynamics/wagner",
        "how_obtained": "both read first-hand (the arXiv PDF text was extracted and searched; the BYU page was fetched). "
                        "The primary report (R. T. Jones, NACA TN 667, 1938, as cited by Dawson & Brunton) was not read",
        "verified": "2026-10-12",
    },
    "wagner_garrick_fit": {
        "form": "phi(s) = 1 - 2 / (4 + s)",
        "max_abs_error_stated": 0.02,
        "what": "I. E. Garrick's algebraic approximation of Wagner's function, same non-dimensional time",
        "source": "Dawson & Brunton, arXiv:2104.15122v2, Eq. (2); absolute error '2%' (text under Fig. 2)",
        "url": "https://arxiv.org/abs/2104.15122",
        "how_obtained": "read first-hand in the extracted PDF text",
        "verified": "2026-10-12",
    },
}


def main() -> int:
    path = HERE / "benchmarks.json"
    path.write_text(json.dumps(BENCHMARKS, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {path} ({path.stat().st_size} bytes, {len(BENCHMARKS) - 1} benchmark entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
