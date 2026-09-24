"""Write the cited public benchmark values and forms used by ``tests/test_ch06.py`` into ``reference/ch06/``.

Chapter 6 (ideal flow) needs no datasets: almost every check is Tier 1 (exact potential-flow solutions, sympy
re-derivations, convergence orders of our own discretisations). The public sources below pin the two published
*numbers* the chapter's code reproduces (the sphere's added-mass coefficient ½ and Besant/Rayleigh's empty-cavity
collapse constant 0.91468 — true V5 benchmarks) and the standard *forms* of the cylinder flow, the Kutta–Joukowski and
Blasius theorems, the Joukowsky transform, the Rankine half-body and the potential flow past a sphere (V1 form
cross-checks: a published closed form reproduced identically is not an independent number). Every entry was read from
the cited page on 2026-09-23 (see ``SOURCES.md``). Nothing here comes from the textbook; book-quoted numbers and printed
closed forms live in the git-ignored ``tests/book_values_ch06.json``.

Run: ``.venv/Scripts/python.exe reference/ch06/make_refs.py``  (idempotent; overwrites ``benchmarks.json``).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_meta": {
        "chapter": "ch06",
        "written_by": "reference/ch06/make_refs.py",
        "verified": "2026-09-23",
        "note": "Public sources only (pages fetched on the date above, see SOURCES.md). No book values.",
    },
    # Wikipedia, "Added mass": "the added mass for a sphere (of radius r) is 2/3 π r³ ρ_fluid, which is half the volume
    # of the sphere times the density of the fluid"; equation of motion (m_p + ρ_c V_p/2) dv/dt = ΣF + ρ_c V_p/2 Du/Dt.
    "added_mass_sphere": {
        "formula": "2/3*pi*r**3*rho",
        "fraction_of_displaced_mass": 0.5,
        "equation_of_motion": "(m_p + rho*V_p/2) dv/dt = sum F (+ rho*V_p/2 Du/Dt for an accelerating fluid)",
        "source": "https://en.wikipedia.org/wiki/Added_mass (fetched 2026-09-23) — secondary",
        "primary": "H. Lamb, Hydrodynamics, 6th ed., Cambridge University Press (1932), Ch. V, §92 (sphere moving through "
                   "liquid: inertia increased by half the mass of liquid displaced) — section number as cited by the "
                   "derivation review; edition confirmed by search, the text itself not fetched",
        "label": "V5",
    },
    # Wikipedia, "Rayleigh–Plesset equation": Besant's empty-cavity filling time t ≈ 0.91468 R0 sqrt(ρ/P∞).
    "rayleigh_collapse": {
        "constant": 0.91468,
        "formula": "t = 0.91468*R0*sqrt(rho/P_inf)",
        "digits": "as printed (5 significant figures)",
        "source": "https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Plesset_equation (fetched 2026-09-23) — secondary",
        "primary": "Lord Rayleigh, On the pressure developed in a liquid during the collapse of a spherical cavity, "
                   "Phil. Mag. (6) 34(200), 94-98 (1917), doi:10.1080/14786440808635681 (citation confirmed by search)",
        "label": "V5",
    },
    # Wikipedia, "Potential flow around a circular cylinder": ϕ = Ur(1 + R²/r²)cos θ, V_r = U(1 − R²/r²)cos θ,
    # V_θ = −U(1 + R²/r²)sin θ; surface pressure coefficient from +1 (θ = 0, π) to −3 (θ = π/2, 3π/2); zero drag.
    "cylinder_form": {
        "phi": "U*r*(1 + R**2/r**2)*cos(theta)",
        "V_r": "U*(1 - R**2/r**2)*cos(theta)",
        "V_theta": "-U*(1 + R**2/r**2)*sin(theta)",
        "cp_surface_max": 1.0,
        "cp_surface_min": -3.0,
        "drag": 0.0,
        "source": "https://en.wikipedia.org/wiki/Potential_flow_around_a_circular_cylinder (fetched 2026-09-23)",
        "label": "V1 form",
    },
    # Wikipedia, "Kutta–Joukowski theorem": L′ = ρ∞V∞Γ with Γ the line integral "followed in the negative (clockwise)
    # direction"; complex force F̄ = (iρ/2)∮ w′² dz.
    "kutta_joukowski_form": {
        "lift": "rho*V*Gamma",
        "circulation_sense": "clockwise (negative) contour",
        "complex_force": "conj(F) = I*rho/2 * contour_integral(w'**2, z)",
        "source": "https://en.wikipedia.org/wiki/Kutta%E2%80%93Joukowski_theorem (fetched 2026-09-23)",
        "label": "V1 form",
    },
    # Wikipedia, "Blasius theorem": F_x − iF_y = (iρ/2)∮_C (dw/dz)² dz; M = Re{−(ρ/2)∮_C z (dw/dz)² dz}.
    "blasius_form": {
        "force": "F_x - I*F_y = I*rho/2 * contour_integral((dw/dz)**2, z)",
        "moment": "M = re(-rho/2 * contour_integral(z*(dw/dz)**2, z))",
        "source": "https://en.wikipedia.org/wiki/Blasius_theorem (fetched 2026-09-23)",
        "label": "V1 form",
    },
    # Wikipedia, "Joukowsky transform": z = ζ + 1/ζ; "the complex unit circle maps to a flat plate on the real-number
    # line from −2 to +2"; W = W̃/(1 − 1/ζ²).
    "joukowsky_form": {
        "map": "z = zeta + 1/zeta",
        "unit_circle_image": [-2.0, 2.0],
        "velocity": "W = W_circle/(1 - 1/zeta**2)",
        "source": "https://en.wikipedia.org/wiki/Joukowsky_transform (fetched 2026-09-23)",
        "label": "V1 form",
    },
    # Wikipedia, "Rankine half body": ψ = Ur sin θ + mθ/2π, stagnation distance U = m/(2πb) (b = m/2πU),
    # body r = b(π − θ)/sin θ.
    "rankine_half_body_form": {
        "psi": "U*r*sin(theta) + m*theta/(2*pi)",
        "b": "m/(2*pi*U)",
        "body": "b*(pi - theta)/sin(theta)",
        "source": "https://en.wikipedia.org/wiki/Rankine_half_body (fetched 2026-09-23)",
        "label": "V1 form",
    },
    # K. T. McDonald, "Pressure in Fluid Flow Past a Sphere", Princeton (March 19, 2015; updated March 25, 2015),
    # eqs. (9), (10), (12), (14): Φ = v(r + a³/2r²)cos θ; u = v(1 − a³/r³)cos θ r̂ − v(1 + a³/2r³)sin θ θ̂;
    # Stokes stream function with u_z = (1/ϱ)∂Ψ/∂ϱ, u_ϱ = −(1/ϱ)∂Ψ/∂z (θ from the stream direction z).
    "sphere_form": {
        "Phi": "v*(r + a**3/(2*r**2))*cos(theta)",
        "u_r": "v*(1 - a**3/r**3)*cos(theta)",
        "u_theta": "-v*(1 + a**3/(2*r**3))*sin(theta)",
        "u_z_cyl": "v*(1 - a**3*(2*z**2 - rho_c**2)/(2*(z**2 + rho_c**2)**(5/2)))",
        "u_rho_cyl": "-v*3*a**3*rho_c*z/(2*(z**2 + rho_c**2)**(5/2))",
        "stream_function_convention": "u_z = (1/rho_c) dPsi/drho_c, u_rho = -(1/rho_c) dPsi/dz",
        "source": "K. T. McDonald, Pressure in Fluid Flow Past a Sphere, Princeton (2015), "
                  "http://kirkmcd.princeton.edu/examples/sphereflow.pdf (fetched 2026-09-23)",
        "label": "V1 form",
    },
    # Method citation only (no number): constant-strength source panels with zero normal velocity at panel midpoints.
    "source_panel_method": {
        "citation": "J. L. Hess & A. M. O. Smith, Calculation of potential flow about arbitrary bodies, "
                    "Prog. Aerosp. Sci. 8, 1-138 (1967), doi:10.1016/0376-0421(67)90003-6",
        "label": "method pedigree (not a number)",
    },
}


def main() -> int:
    HERE.mkdir(parents=True, exist_ok=True)
    (HERE / "benchmarks.json").write_text(json.dumps(BENCHMARKS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {HERE / 'benchmarks.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
