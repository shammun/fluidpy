"""Write the cited public benchmark values used by ``tests/test_ch03.py`` (V5 evidence) into ``reference/ch03/``.

Chapter 3 (kinematics) needs no datasets: almost every check is Tier 1 (analytic, symbolic, convergence). The public
sources below pin the *forms* the chapter's vortex models and transport theorems must have (Rankine vortex, Lamb–Oseen
vortex, Reynolds transport theorem, Leibniz integral rule) and one published number: the Lamb–Oseen radius of maximum
tangential velocity, r_max² = α r_c² with α = 1.256 (Canivete Cuissa & Steiner 2022, Eq. (10)). Every entry was read
from the cited page on 2026-09-23 (see ``SOURCES.md``). Nothing here comes from the textbook; book-quoted numbers live
in the git-ignored ``tests/book_values_ch03.json``.

Run: ``.venv/Scripts/python.exe reference/ch03/make_refs.py``  (idempotent; overwrites ``benchmarks.json``).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_meta": {
        "chapter": "ch03",
        "written_by": "reference/ch03/make_refs.py",
        "verified": "2026-09-23",
        "note": "Public sources only (pages fetched on the date above, see SOURCES.md). No book values.",
    },
    # J. R. Canivete Cuissa & O. Steiner (2022), "An innovative and automated method for vortex identification.
    # I. Description of the SWIRL algorithm", arXiv:2210.05223 (A&A), Sect. 2.4, Eq. (10):
    #   v_theta = v_max (1 + 1/(2 alpha)) (r_max/r) [1 - exp(-alpha r^2/r_max^2)],  v_r = 0,
    # "where alpha = 1.256 establishes that r_max is the radius at which the rotational velocity v_theta is maximal
    # and reaches the value v_max". With sigma^2 = r_max^2/alpha this is the book's Gaussian vortex (3.29), so
    # r_max/sigma = sqrt(alpha); the normalisation (1 + 1/(2 alpha))(1 - e^{-alpha}) = 1 is equivalent to e^alpha = 1 + 2 alpha.
    "swirl_lamb_oseen": {
        "velocity": "v_max*(1 + 1/(2*alpha))*(r_max/r)*(1 - exp(-alpha*r**2/r_max**2))",
        "alpha": 1.256,
        "alpha_digits": 3,
        "alpha_half_ulp": 5e-4,
        "meaning": "r_max**2 = alpha * sigma**2 when exp(-alpha r^2/r_max^2) is written exp(-r^2/sigma^2)",
        "source": "Canivete Cuissa & Steiner (2022), arXiv:2210.05223, Sect. 2.4, Eq. (10) (PDF text read 2026-09-23)",
        "url": "https://arxiv.org/abs/2210.05223",
    },
    # Wikipedia, "Rankine vortex": v_theta = (Gamma/2 pi) r/a^2 (r <= a), (Gamma/2 pi)/r (r > a); omega_z = 2 Omega
    # inside (uniform), 0 outside.
    "rankine_vortex": {
        "v_theta_inside": "Gamma/(2*pi)*r/a**2",
        "v_theta_outside": "Gamma/(2*pi)/r",
        "omega_inside": "2*Omega",
        "omega_outside": "0",
        "Omega_of_Gamma": "Gamma/(2*pi*a**2)",
        "source": "https://en.wikipedia.org/wiki/Rankine_vortex (fetched 2026-09-23)",
    },
    # Wikipedia, "Lamb–Oseen vortex": v_theta = Gamma/(2 pi r) (1 - exp(-r^2/(4 nu t))),
    # omega_z = Gamma/(4 pi nu t) exp(-r^2/(4 nu t)).
    "lamb_oseen_vortex": {
        "v_theta": "Gamma/(2*pi*r)*(1 - exp(-r**2/(4*nu*t)))",
        "omega_z": "Gamma/(4*pi*nu*t)*exp(-r**2/(4*nu*t))",
        "source": "https://en.wikipedia.org/wiki/Lamb%E2%80%93Oseen_vortex (fetched 2026-09-23)",
    },
    # Wikipedia, "Reynolds transport theorem": d/dt int_Omega(t) f dV = int ∂f/∂t dV + int_∂Omega(t) (v_b·n) f dA,
    # n outward, v_b the velocity of the area element (not the flow velocity); fixed region (v_b = 0) → d/dt passes
    # inside the integral.
    "reynolds_transport_theorem": {
        "statement": "d/dt int_{Omega(t)} f dV = int_{Omega(t)} df/dt dV + int_{dOmega(t)} (v_b . n) f dA",
        "normal": "outward-pointing unit normal",
        "v_b": "velocity of the area element (not the flow velocity)",
        "fixed_region": "v_b = 0 -> d/dt int f dV = int df/dt dV",
        "source": "https://en.wikipedia.org/wiki/Reynolds_transport_theorem (fetched 2026-09-23)",
    },
    # Wikipedia, "Leibniz integral rule": d/dx int_{a(x)}^{b(x)} f(x, t) dt
    #   = f(x, b(x)) b'(x) - f(x, a(x)) a'(x) + int_{a(x)}^{b(x)} ∂f/∂x dt.
    "leibniz_integral_rule": {
        "statement": "d/dx int_a(x)^b(x) f(x,t) dt = f(x,b(x)) db/dx - f(x,a(x)) da/dx + int_a^b df/dx dt",
        "upper_sign": 1,
        "lower_sign": -1,
        "source": "https://en.wikipedia.org/wiki/Leibniz_integral_rule (fetched 2026-09-23)",
    },
}


def main() -> int:
    out = HERE / "benchmarks.json"
    out.write_text(json.dumps(BENCHMARKS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {out} ({len(BENCHMARKS) - 1} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
