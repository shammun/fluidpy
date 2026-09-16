"""Write the cited public benchmark values used by ``tests/test_ch02.py`` (V5 evidence) into ``reference/ch02/``.

Chapter 2 is pure mathematics, so almost all evidence is Tier 1 (analytic / symbolic). The few citable public
statements below pin *conventions* (index placement of Cauchy's formula, active vs passive rotation matrices, the
right-hand rule of Stokes' theorem) and two closed-form numbers (the divergence-theorem example 8π/3 and the
Levi-Civita contraction constants). Every entry was read from the cited page on 2026-09-16 (see ``SOURCES.md``).
Nothing here comes from the textbook; book-quoted numbers live in the git-ignored ``tests/book_values_ch02.json``.

Run: ``.venv/Scripts/python.exe reference/ch02/make_refs.py``  (idempotent; overwrites ``benchmarks.json``).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

BENCHMARKS = {
    "_meta": {
        "chapter": "ch02",
        "written_by": "reference/ch02/make_refs.py",
        "verified": "2026-09-16",
        "note": "Public sources only (Wikipedia pages fetched on the date above, see SOURCES.md). No book values.",
    },
    # Wikipedia, "Divergence theorem", worked example: F = 2x i + y^2 j + z^2 k through the unit sphere.
    # div F = 2(1 + y + z); the odd terms integrate to zero over the ball; flux = 2 * (4 pi / 3) = 8 pi / 3.
    "divergence_theorem_sphere_example": {
        "field": ["2*x", "y**2", "z**2"],
        "divergence": "2*(1 + y + z)",
        "radius": 1.0,
        "flux": 8.0 * math.pi / 3.0,
        "flux_closed_form": "8*pi/3",
        "source": "https://en.wikipedia.org/wiki/Divergence_theorem (fetched 2026-09-16)",
    },
    # Wikipedia, "Levi-Civita symbol", section "Three dimensions": product identities and the determinant formula.
    "levi_civita_identities": {
        "eps_ijk_eps_imn": "delta_jm delta_kn - delta_jn delta_km",
        "eps_jmn_eps_imn": "2 delta_ij",
        "eps_ijk_eps_ijk": 6,
        "contraction_2delta_factor": 2,
        "full_contraction": 6,
        "delta_ii_3d": 3,
        "cross_product": "(a x b)_i = eps_ijk a_j b_k",
        "determinant": "det A = eps_ijk a_1i a_2j a_3k",
        "pseudotensor": "changes sign under an improper rotation (det = -1)",
        "source": "https://en.wikipedia.org/wiki/Levi-Civita_symbol (fetched 2026-09-16)",
    },
    # Wikipedia, "Cauchy stress tensor": index placement of the traction, transformation rule, invariants, Mohr.
    "cauchy_stress_tensor": {
        "traction": "T_j = sigma_ij n_i  (the FIRST index of sigma is contracted with n)",
        "transformation": "sigma'_ij = a_im a_jn sigma_mn, i.e. sigma' = A sigma A^T with A the active rotation (A = C^T in the book's passive convention)",
        "I1": "tr(sigma) = sigma_1 + sigma_2 + sigma_3",
        "I2": "1/2 [ (tr sigma)^2 - tr(sigma^2) ] = sigma_1 sigma_2 + sigma_2 sigma_3 + sigma_3 sigma_1",
        "I3": "det(sigma) = sigma_1 sigma_2 sigma_3",
        "principal_2d": "sigma_1,2 = (sigma_x + sigma_y)/2 +/- sqrt( ((sigma_x - sigma_y)/2)^2 + tau_xy^2 )",
        "max_shear": "tau_max = |sigma_1 - sigma_3| / 2",
        "source": "https://en.wikipedia.org/wiki/Cauchy_stress_tensor (fetched 2026-09-16)",
    },
    # Wikipedia, "Rotation matrix": R(theta) rotates vectors counterclockwise (active); passive = transpose.
    "rotation_matrix": {
        "R_2d": [["cos", "-sin"], ["sin", "cos"]],
        "meaning": "active: rotates points counterclockwise by theta; a rotation of the axes (passive/alias) uses R^T",
        "properties": "R^T R = I, det R = +1",
        "rodrigues": "R = cos(theta) I + sin(theta) [u]_x + (1 - cos(theta)) u u^T",
        "source": "https://en.wikipedia.org/wiki/Rotation_matrix (fetched 2026-09-16)",
    },
    # Wikipedia, "Stokes' theorem": statement and orientation rule (no numeric example on the page).
    "stokes_theorem": {
        "statement": "oint_{dSigma} F . dr = iint_Sigma (curl F) . n dS",
        "orientation": "the positive direction around the boundary and the positive normal are related by the right-hand rule",
        "numeric_example_on_page": False,
        "our_V1_case": "u = b x x, disc of radius R perpendicular to b: circulation = 2 |b| pi R^2",
        "source": "https://en.wikipedia.org/wiki/Stokes%27_theorem (fetched 2026-09-16)",
    },
}


def main() -> None:
    out = HERE / "benchmarks.json"
    out.write_text(json.dumps(BENCHMARKS, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
