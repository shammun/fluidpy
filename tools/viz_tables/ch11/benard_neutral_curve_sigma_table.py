"""Generate the growth-rate table embedded in viz/ch11/benard_neutral_curve.html (constant SIGTAB).

What is tabulated (ours, computed by fluidpy; not book data)
    c(K, x) = sigma_true / sigma_one_mode,   x = Ra / Ra_m(K),
for rigid-rigid and rigid-free walls and the two Prandtl numbers the page offers (water nu/kappa = 1e-6/1.4e-7,
air 1.5e-5/2.2e-5), where
    sigma_true      = ch11.benard_growth_rate(K, Ra, Pr, bc)            leading Chebyshev eigenvalue [kappa/d^2]
    Ra_m(K)         = ch11.benard_marginal_Ra(K, bc, mode)              marginal Rayleigh number
    sigma_one_mode  = root of (s + a2)(s/Pr + a2) = a2^2 x, a2 = (Ra_m K^2)^(1/3)   (the explainer's sigmaEst)
The explainer shows sigma = c * sigma_one_mode with its LIVE Ra_m(K), so sigma = 0 lies exactly on the live neutral curve
(sigma_one_mode vanishes there) and c only has to carry a smooth factor between about 0.7 and 1.05.
The column x = 1 holds the limit of c (mean of x = 1 -/+ 1e-3).

Robustness: benard_growth_rate raises ValueError when no eigenvalue survives its N-convergence filter -> retry with a
larger N and with filter=False, and a value is kept only when two settings agree (see sigma_true); for Ra > 0 the leading sigma is real (exchange of stabilities) and
return_complex=True is used to assert that (|Im| < 1e-6 max(1, |sigma|)); Ra <= 0 is never requested (x >= X[0] > 0).

Usage (repo root):
    .venv/Scripts/python.exe reports/viz/ch11/benard_neutral_curve/gen_sigma_table.py            # table -> sigma_table.json
    .venv/Scripts/python.exe reports/viz/ch11/benard_neutral_curve/gen_sigma_table.py --check 60 # + interpolation error
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT))
from fluidpy import ch11_instability as ch11  # noqa: E402

PRS = {"water": 1.0e-6 / 1.4e-7, "air": 1.5e-5 / 2.2e-5}
BCS = {"rigid": (("rigid", "rigid"), "even"), "rigidfree": (("rigid", "free"), "any")}
NK = 25
KS = np.exp(np.linspace(math.log(0.25), math.log(10.0), NK))       # log-spaced wavenumbers
XS = np.array([0.002, 0.02, 0.06, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.6, 2.0, 2.6, 3.5, 5.0, 7.0, 10.0, 14.0, 20.0, 30.0])


def sigma_est(K: float, Ra: float, Pr: float, RaK: float) -> float:
    a2 = (RaK * K * K) ** (1.0 / 3.0)
    b = a2 * (1.0 + Pr)
    c = Pr * a2 * a2 * (1.0 - Ra / RaK)
    return (-b + math.sqrt(max(0.0, b * b - 4.0 * c))) / 2.0


def sigma_true(K: float, Ra: float, Pr: float, bc, RaK: float | None = None) -> float:
    """Leading growth rate, cross-checked: benard_growth_rate can silently return the SECOND eigenvalue when the leading
    one misses its 1e-6 N-convergence filter (seen at isolated points, e.g. rigid-free K = 1.84, Ra = 1.6 Ra_m, Pr = 7.14),
    so several resolutions are computed and the largest value confirmed by a second setting is kept."""
    cands = []
    for kw in (dict(N=40), dict(N=56), dict(N=56, filter=False), dict(N=72, filter=False)):
        try:
            s = ch11.benard_growth_rate(K, Ra, Pr, bc, return_complex=True, **kw)
        except ValueError:           # no eigenvalue survived the filter at this N: use the other settings
            continue
        if abs(s.imag) > 1e-4 * max(1.0, abs(s)):   # Ra > 0: the leading sigma is real; round-off Im ~ 1e-6 is tolerated
            continue                                # a genuinely complex candidate is not accepted (never used: Ra > 0 only)
        cands.append(float(s.real))
    ok = [u for i, u in enumerate(cands) if any(j != i and abs(u - w) <= 1e-5 * max(1.0, abs(u)) for j, w in enumerate(cands))]
    if not ok:
        raise RuntimeError(f"benard_growth_rate not confirmed at K={K}, Ra={Ra}, Pr={Pr}, bc={bc}: {cands}")
    s = max(ok)
    if RaK is not None and abs(Ra / RaK - 1.0) > 1e-6 and (s > 0) != (Ra > RaK):
        raise RuntimeError(f"sign of sigma {s} contradicts the marginal curve at K={K}, Ra={Ra}, RaK={RaK}")
    return s


def ratio(K: float, x: float, Pr: float, name: str, RaK: float | None = None) -> float:
    bc, mode = BCS[name]
    RaK = ch11.benard_marginal_Ra(K, bc, mode) if RaK is None else RaK
    if abs(x - 1.0) < 1e-12:
        return 0.5 * (ratio(K, 1.0 - 1e-3, Pr, name, RaK) + ratio(K, 1.0 + 1e-3, Pr, name, RaK))
    return sigma_true(K, x * RaK, Pr, bc, RaK) / sigma_est(K, x * RaK, Pr, RaK)


def lagrange4(xs, ys, x):
    """4-point Lagrange interpolation on the nodes around x (the JS `lag4`)."""
    n = len(xs)
    i = int(np.clip(np.searchsorted(xs, x) - 1, 0, n - 2))
    i0 = int(np.clip(i - 1, 0, n - 4))
    out = 0.0
    for a in range(i0, i0 + 4):
        w = 1.0
        for b in range(i0, i0 + 4):
            if b != a:
                w *= (x - xs[b]) / (xs[a] - xs[b])
        out += w * ys[a]
    return out


def interp(tab, K, x):
    lk, lx = np.log(KS), np.log(XS)
    rows = [lagrange4(lx, tab[i], math.log(x)) for i in range(NK)]
    return lagrange4(lk, rows, math.log(K))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", type=int, default=0, help="number of random off-grid points per table to test")
    a = ap.parse_args()
    out = {"note": "ours: c = benard_growth_rate / one-mode estimate at x = Ra/Ra_m(K); 4 s.f.", "K": [round(float(k), 6) for k in KS],
           "x": [float(x) for x in XS], "Pr": PRS, "c": {}}
    for name in BCS:
        bc, mode = BCS[name]
        RaKs = [ch11.benard_marginal_Ra(float(K), bc, mode) for K in KS]
        for pn, Pr in PRS.items():
            tab = [[float(f"{ratio(float(K), float(x), Pr, name, RaKs[i]):.4g}") for x in XS] for i, K in enumerate(KS)]
            out["c"][name + "_" + pn] = tab
            flat = np.array(tab)
            print(name, pn, "c range", flat.min(), flat.max(), flush=True)
            if a.check:
                rng = np.random.default_rng(11)
                worst = 0.0
                for _ in range(a.check):
                    K = float(np.exp(rng.uniform(math.log(0.25 if name == "rigidfree" else 0.5), math.log(10.0))))
                    RaK = ch11.benard_marginal_Ra(K, bc, mode)
                    lo, hi = max(XS[0], 398.0 / RaK / (16 if name == "rigidfree" and K < 0.5 else 1)), min(30.0, 31623.0 / RaK)
                    if hi <= lo:
                        continue
                    x = float(np.exp(rng.uniform(math.log(lo), math.log(hi))))
                    if abs(x - 1) < 0.02:
                        continue
                    st = sigma_true(K, x * RaK, Pr, bc, RaK)
                    si = interp(tab, K, x) * sigma_est(K, x * RaK, Pr, RaK)
                    err = abs(si / st - 1)
                    if err > worst:
                        worst = err
                        wk = (K, x, st, si)
                print("   worst interpolation error", f"{worst:.2e}", "at (K, x, true, table)", wk, flush=True)
    (HERE / "sigma_table.json").write_text(json.dumps(out), encoding="utf-8")
    print("wrote", HERE / "sigma_table.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
