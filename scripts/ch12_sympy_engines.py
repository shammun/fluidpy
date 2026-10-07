"""Run the nine symbolic engines of chapter 12 and print their checks (each derived equation minus the book's form must
simplify to zero): Reynolds averaging (12.27)-(12.34), the Reynolds-stress budget (12.35), the two energy budgets (12.46),
(12.47), the temperature-variance budget (12.112) with its printed slip, the overlap matching (12.85)-(12.89), the plane-jet
similarity equation (12.63)-(12.74), the isotropic gradient moments behind (12.43) and the divergence of (12.41).

Run: ``.venv/Scripts/python.exe scripts/ch12_sympy_engines.py --no-show [--no-cache]``
Results are cached in outputs/ch12/cache/sympy_<name>.json (delete the folder after changing an engine).
"""
from __future__ import annotations

from ch12_common import Timer, finish, parse_args, setup

from fluidpy import ch12_turbulence as ch12


def main() -> int:
    args = parse_args(__doc__, extra=lambda ap: ap.add_argument("--no-cache", action="store_true"))
    setup(args)
    bad = []
    expected_false = {"printed_check"}          # the printed form of (12.112) must NOT pass
    for name in ch12.SYMPY_ENGINES:
        with Timer(name):
            s = ch12.sympy_summary(name, cache=not args.no_cache)
        print(f"  {name}: {s['checks']}  (computed in {s['seconds']} s)")
        for k, v in s["checks"].items():
            if v != (k not in expected_false):
                bad.append(f"{name}.{k}")
    gm = ch12.gradient_moments_isotropic_sympy()
    print(f"gradient moments (units u2/lambda_f^2): {gm['m11']}, {gm['m12']}, {gm['cross']}; factor {gm['eps_factor']} (lambda_f), "
          f"{gm['eps_factor_g']} (lambda_g); lambda_g^2/lambda_f^2 = {gm['ratio_lambda']}")
    pj = ch12.plane_jet_similarity_sympy()
    print(f"plane-jet coefficients, power family: {pj['power_family']}; exponential family: {pj['exponential_family']}")
    tv = ch12.temperature_variance_sympy()
    print(f"(12.112) as printed minus derived: {tv['printed_minus_derived']}")
    if bad:
        print("UNEXPECTED:", bad)
    finish(args)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
