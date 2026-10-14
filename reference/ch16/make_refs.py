"""Write ``reference/ch16/benchmarks.json`` — the cited public values and forms the chapter-16 tests compare with.

Only statements read first-hand in a public source on the date given are stored (see SOURCES.md for the full citation
and what was read).  Nothing here comes from the textbook.  Equations are stored as **forms** (a formula as the source
writes it); the tests evaluate our functions against the form at our own inputs.  Numbers are stored exactly as the
source prints them.

Run: ``.venv/Scripts/python.exe reference/ch16/make_refs.py``  (idempotent: rewrites the same bytes).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TODAY = "2026-10-14"

BENCHMARKS = {
    "_note": "Cited public values and forms for biofluid mechanics; our own numbers are computed, never typed. See SOURCES.md.",
    "pries_secomb_pries_2013": {
        "what": "in vitro relative apparent viscosity of blood in glass tubes, Eqs. (3)-(5): "
                "mu_rel = 1 + (mu_045 - 1)[(1 - H_D)^C - 1]/[(1 - 0.45)^C - 1]; "
                "mu_045 = 220 exp(-1.3 D) + 3.2 - 2.44 exp(-0.06 D^0.645); "
                "C = (0.8 + exp(-0.075 D))(-1 + (1 + 1e-11 D^12)^-1) + (1 + 1e-11 D^12)^-1; D in micrometres",
        "coefficients": {"a1": 220.0, "b1": 1.3, "a2": 3.2, "a3": 2.44, "b3": 0.06, "e3": 0.645, "c1": 0.8, "c2": 0.075,
                         "k": 1e-11, "m": 12, "H_ref": 0.45},
        "statement_min_diameter_um": 7.0,
        "statement": "The trend of apparent viscosity to decrease with decreasing diameter continues down to about 7 um. "
                     "At smaller diameters, relative apparent viscosity rises rapidly as the diameter approaches a critical "
                     "minimum diameter, which is about 2.7 um.",
        "source": "T. W. Secomb & A. R. Pries, 'Blood viscosity in microvessels: experiment and theory', "
                  "C. R. Physique 14, 470-478 (2013), Section 2.1, Eqs. (3)-(5)",
        "url": "https://comptes-rendus.academie-sciences.fr/physique/item/10.1016/j.crhy.2013.04.002.pdf",
        "doi": "10.1016/j.crhy.2013.04.002",
        "how_obtained": "PDF fetched; text extracted with pymupdf; equations and the sentence read and typed",
        "verified": TODAY,
    },
    "moens_korteweg_wikipedia": {
        "what": "PWV = sqrt(E_inc h / (2 r rho)), r the vessel radius (form)",
        "source": "Wikipedia, 'Moens-Korteweg equation' (secondary)",
        "url": "https://en.wikipedia.org/wiki/Moens%E2%80%93Korteweg_equation",
        "how_obtained": "page fetched and read",
        "verified": TODAY,
    },
    "windkessel_wikipedia": {
        "what": "two-element Windkessel I(t) = P(t)/R + C dP/dt; diastolic decay P(t) = P(t_d) exp(-(t - t_d)/(R C)) (forms)",
        "source": "Wikipedia, 'Windkessel effect', section 'Two-element' (secondary)",
        "url": "https://en.wikipedia.org/wiki/Windkessel_effect",
        "how_obtained": "page fetched and read",
        "verified": TODAY,
    },
    "murray_wikipedia": {
        "what": "total power P = (8 mu/pi) l Q^2/r^4 + pi lambda l r^2; minimum at Q proportional to r^3; "
                "r^3 = r1^3 + r2^3 at a junction (forms); the page gives no branching angle",
        "source": "Wikipedia, 'Murray's law' (secondary; cites C. D. Murray, PNAS 12, 207 (1926))",
        "url": "https://en.wikipedia.org/wiki/Murray%27s_law",
        "how_obtained": "page fetched and read",
        "verified": TODAY,
    },
    "pulsatile_flow_wikipedia": {
        "what": "Womersley's solution with dp/dx = sum_n P'_n e^{i n omega t}, Lambda_n = alpha n^(1/2) i^(3/2), "
                "alpha = R sqrt(omega rho/mu): u = Re{sum i P'_n/(rho n omega) [1 - J0(Lambda_n r/R)/J0(Lambda_n)] e^{i n omega t}}, "
                "Q = Re{pi R^2 sum i P'_n/(rho n omega) [1 - 2 J1(Lambda_n)/(Lambda_n J0(Lambda_n))] e^{i n omega t}}, "
                "tau_w = mu du/dr|_R = Re{sum P'_n R J1(Lambda_n)/(Lambda_n J0(Lambda_n)) e^{i n omega t}} (forms)",
        "source": "Wikipedia, 'Pulsatile flow', section 'Womersley solution' (secondary; cites J. R. Womersley, "
                  "J. Physiol. 127, 553 (1955))",
        "url": "https://en.wikipedia.org/wiki/Pulsatile_flow",
        "how_obtained": "page fetched and read",
        "verified": TODAY,
    },
    "entrance_length_wikipedia": {
        "what": "hydrodynamic entrance length L_h = 0.0575 Re_D D (laminar, Bergman & Incropera), "
                "L_h = 1.359 D Re_D^(1/4) (turbulent, Cengel)",
        "laminar_coefficient": 0.0575,
        "turbulent_coefficient": 1.359,
        "source": "Wikipedia, 'Entrance length (fluid dynamics)' (secondary)",
        "url": "https://en.wikipedia.org/wiki/Entrance_length_(fluid_dynamics)",
        "how_obtained": "page fetched and read",
        "verified": TODAY,
    },
    "dean_equations_imperial_m4a33": {
        "what": "Dean equations from the Navier-Stokes equations for a slightly curved pipe, eq. (9): "
                "K(psi_z v_x - psi_x v_z) = 1 + lap v; K(psi_z Om_x - psi_x Om_z) = lap Om - 2 K v v_z; Om = -lap psi; "
                "grad psi = 0, v = 0 on the wall; K = (rho U0 a/mu)(a/b)^(1/2), U0 = G a^2/mu (form; no numbers)",
        "source": "Imperial College London, Department of Mathematics, course M4A33 lecture notes 'Flow in curved pipes: "
                  "the Dean equations' (author's page ~ajm8)",
        "url": "https://www.ma.ic.ac.uk/~ajm8/M4A33/Dean.pdf",
        "how_obtained": "PDF fetched; text extracted with pymupdf; equations (7)-(10) read; the vorticity equation checked "
                        "against the curl of the momentum equations by hand",
        "verified": TODAY,
    },
}


def main() -> None:
    path = HERE / "benchmarks.json"
    path.write_text(json.dumps(BENCHMARKS, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", path)


if __name__ == "__main__":
    main()
