"""Write the cited public benchmark values and forms used by ``tests/test_ch05.py`` into ``reference/ch05/``.

Chapter 5 (vorticity dynamics) needs no datasets: almost every check is Tier 1 (exact solutions of the Euler and
Navier–Stokes equations, sympy re-derivations, invariants of point-vortex motion, convergence orders). The public
sources below pin the handful of published results the chapter touches: Kelvin's thin-ring speed (the only published
*number-producing* formula our ring model relies on), the corotating/counter-rotating point-vortex pair rates of
García & Haziot (2023, §2.1), and the standard forms of the Burgers vortex, Hill's spherical vortex, the aerodynamic
Biot–Savart segment law and the Kelvin/Poincaré–Bjerknes circulation theorems (form cross-checks). Every entry was read
from the cited page on 2026-09-23 (see ``SOURCES.md``). Nothing here comes from the textbook; book-quoted numbers and
printed closed forms live in the git-ignored ``tests/book_values_ch05.json``. Earth's rotation rate (WGS-84) is reused
from ``reference/ch04/``.

Run: ``.venv/Scripts/python.exe reference/ch05/make_refs.py``  (idempotent; overwrites ``benchmarks.json``).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_meta": {
        "chapter": "ch05",
        "written_by": "reference/ch05/make_refs.py",
        "verified": "2026-09-23",
        "note": "Public sources only (pages fetched on the date above, see SOURCES.md). No book values.",
    },
    # Wikipedia, "Vortex ring" (Kelvin's thin-core ring): "U = Γ/(4πR)(ln(8R/a) − 1/4)", R ring radius, a core radius,
    # thin core with uniform vorticity.
    "kelvin_ring_speed": {
        "formula": "Gamma/(4*pi*R)*(log(8*R/a) - 1/4)",
        "core": "uniform (thin core, a << R)",
        "source": "https://en.wikipedia.org/wiki/Vortex_ring (fetched 2026-09-23)",
    },
    # C. García & S. V. Haziot, "Global bifurcation for corotating and counter-rotating vortex pairs", Commun. Math.
    # Phys. (2023), doi:10.1007/s00220-023-04741-6, arXiv:2204.11327 §2.1: point vortices at z1(0) = l, z2(0) = −l with
    # Γ1 = 1, Γ2 = a. a = 1: z_m(t) = e^{iΩ0 t} z_m(0), Ω0 = 1/(4πl²). a = −1: z_m(t) = z_m(0) + V0 t, V0 = −i/(4πl).
    "point_vortex_pairs": {
        "positions": "z1(0) = l, z2(0) = -l",
        "corotating_Gamma": [1.0, 1.0],
        "corotating_Omega0": "1/(4*pi*l**2)",
        "counter_rotating_Gamma": [1.0, -1.0],
        "counter_rotating_V0": "-i/(4*pi*l)",
        "source": "García & Haziot, Commun. Math. Phys. (2023), doi:10.1007/s00220-023-04741-6, arXiv:2204.11327 §2.1",
        "url": "https://arxiv.org/html/2204.11327",
    },
    # Wikipedia, "Burgers vortex": v_r = −αr, v_z = 2αz, v_θ = (Γ/2πr)[1 − exp(−αr²/2ν)], ω_z = (αΓ/2πν)exp(−αr²/2ν).
    "burgers_vortex_form": {
        "v_r": "-alpha_w*r",
        "v_z": "2*alpha_w*z",
        "v_theta": "Gamma/(2*pi*r)*(1 - exp(-alpha_w*r**2/(2*nu)))",
        "omega_z": "alpha_w*Gamma/(2*pi*nu)*exp(-alpha_w*r**2/(2*nu))",
        "parameter_map": "alpha_w (Wikipedia) = alpha/2 (our/book axial stretching rate du_z/dz = alpha)",
        "source": "https://en.wikipedia.org/wiki/Burgers_vortex (fetched 2026-09-23)",
    },
    # Wikipedia, "Hill's spherical vortex": ψ_in = −(3U/4)(1 − r²/a²) r² sin²θ, ψ_out = (U/2)(1 − a³/r³) r² sin²θ,
    # ω_φ = −(15U/2a²) r sin θ; U the free-stream speed far away (frame moving with the vortex).
    "hill_vortex_form": {
        "psi_inside": "-(3*U/4)*(1 - r**2/a**2)*r**2*sin(theta)**2",
        "psi_outside": "(U/2)*(1 - a**3/r**3)*r**2*sin(theta)**2",
        "omega_phi": "-(15*U/(2*a**2))*r*sin(theta)",
        "source": "https://en.wikipedia.org/wiki/Hill%27s_spherical_vortex (fetched 2026-09-23)",
    },
    # Wikipedia, "Biot–Savart law", Aerodynamics: infinite line v = Γ/(2πr); finite segment v = Γ/(4πr)[cos A − cos B].
    "biot_savart_segment_form": {
        "infinite_line": "Gamma/(2*pi*r)",
        "segment": "Gamma/(4*pi*r)*(cos(A) - cos(B))",
        "source": "https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law (fetched 2026-09-23)",
    },
    # Wikipedia, "Kelvin's circulation theorem": barotropic, ideal fluid, conservative body forces, material curve:
    # DΓ/Dt = 0; Poincaré–Bjerknes: Γ = ∫_A (ω + 2Ω)·n dS conserved in a frame rotating at constant Ω.
    "kelvin_bjerknes_form": {
        "kelvin": "DGamma/Dt = 0 (barotropic, ideal, conservative body forces, material curve)",
        "poincare_bjerknes": "int_A (omega + 2*Omega).n dS conserved (constant Omega)",
        "source": "https://en.wikipedia.org/wiki/Kelvin%27s_circulation_theorem (fetched 2026-09-23)",
    },
}


def main() -> int:
    HERE.mkdir(parents=True, exist_ok=True)
    (HERE / "benchmarks.json").write_text(json.dumps(BENCHMARKS, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {HERE / 'benchmarks.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
