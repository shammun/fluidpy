"""Generate the table embedded in viz/ch11/orr_sommerfeld_neutral_curve.html (constant OS).

What is tabulated (ours, computed by fluidpy; not book data)
    For plane Poiseuille flow, the Blasius boundary layer, the tanh shear layer, the Bickley jet (sinuous) and plane
    Couette flow, on a grid of (Re, k): the Orr-Sommerfeld eigenvalue c = c_r + i c_i (Eqs. (11.79), (11.80)) and the
    energy budget P, Lambda, E (Eq. (11.88), mode normalised to max|u_hat| = 1).
    Poiseuille, Blasius, tanh: ONE mode family, FOLLOWED BY CONTINUATION from the most unstable node (at each new node
    the eigenvalue nearest to the extrapolated neighbour's is kept), on the nodes of reference/ch11/os_grid_<flow>.csv.
    That reference grid holds the least-damped discrete mode instead, which changes family outside the neutral curve
    (wall mode <-> centre mode, ...) and therefore cannot be interpolated smoothly there.
    Bickley jet: the least-damped DISCRETE mode at every node (ch11.os_leading_mode: localised and unchanged in a box
    1.5 times as large — the rule of the reference grid), on a finer grid, 36 Re x 37 k with k log-spaced from 0.02: the
    jet has two unstable bands carried by two
    mode families (ch11.bickley_neutral_curve), so a single continuation misses one of them; the explainer interpolates
    this table bilinearly ("bilin").
    The script prints how each table compares with the reference grid on every reference node.
    Numerics: exactly fluidpy's — the wavelength-scaled box ch11.os_box_numerics (max(profile box, 12/k)) for the open
    flows, the Chebyshev degrees of the reference grids (N = 80; tanh 100).
    Per node a flag: "1" the followed mode is the least-damped localised eigenvalue of its box, "0" it is localised but
    another localised mode decays more slowly, "x" its eigenfunction fills the box (far-field fraction >= 0.2: not a mode
    of the unbounded flow but the box's version of the continuous spectrum).
    Plane Couette: the least-damped mode with c_r >= 0 (modes come in pairs +c_r, -c_r), Re 400 ... 4e4.
    Also: mode shapes (phi, u_hat on 29 heights) on a 5 x 4 sub-grid and at the preset points (direct solves, N = 100;
    Bickley 80), the neutral curves of reference/ch11/os_neutral_*.csv, the Bickley stable gap and long-wave band
    (ch11.bickley_neutral_curve, + two extra Re near the opening of the gap), the critical points of
    reference/ch11/critical_points.json and the Blasius profile U(y), U'(y) in delta* units.

Usage (repo root):
    .venv/Scripts/python.exe tools/viz_tables/ch11/orr_sommerfeld_neutral_curve_table.py            # all flows (~15 min)
    .venv/Scripts/python.exe tools/viz_tables/ch11/orr_sommerfeld_neutral_curve_table.py --flows bickley tanh
Output: outputs/ch11/viz_tables/orr_sommerfeld_neutral_curve_table.json (per-flow pieces are cached beside it as
orr_sommerfeld_neutral_curve_<flow>.json; --flows recomputes only those).
"""
from __future__ import annotations

import argparse
import heapq
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import BarycentricInterpolator

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from fluidpy import ch11_instability as ch11  # noqa: E402
from fluidpy.core import stability as ST  # noqa: E402

REF = ROOT / "reference" / "ch11"
OUT = ROOT / "outputs" / "ch11" / "viz_tables"
FLOWS = ("poiseuille", "blasius", "tanh", "bickley", "couette")
NTAB = dict(poiseuille=80, blasius=80, tanh=100, bickley=80, couette=80)
NY = 29
FAR_TOL = 0.2
GAP_NOSE = (17.5, 0.068)  # ch11.bickley_neutral_curve docstring: the gap opens between Re = 17.4 and 17.6 at k ~ 0.068


def sig(v, n=4):
    if v is None or not np.isfinite(v):
        return None
    return float(f"{float(v):.{n}g}")


def load_csv(p: Path) -> np.ndarray:
    rows = [ln for ln in p.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    return np.array([[float(x) for x in r.split(",")] for r in rows[1:]])


_PR: dict = {}


def profile(name: str) -> dict:
    if name not in _PR:
        _PR[name] = ch11.parallel_profile(name, parity="even") if name == "bickley" else ch11.parallel_profile(name)
    return _PR[name]


def spectrum(name: str, k: float, Re: float, N: int | None = None):
    """all eigenvalues (|c| < 5, descending c_i) with eigenvectors, in fluidpy's box for this k"""
    pr = profile(name)
    nm = ch11.os_box_numerics(pr, float(k), N or NTAB[name])
    res = ch11.orr_sommerfeld_eigs(float(k), float(Re), pr["U"], pr["Upp"], domain=pr["domain"], N=nm["N"], bc=pr["bc"],
                                   y_max=nm["y_max"], map_scale=nm["map_scale"], parity=pr.get("parity"),
                                   return_vectors=True, filter=False)
    keep = np.abs(res["c"]) < 5
    return res["c"][keep], res["phi"][:, keep], res["grid"], nm


def far_of(name: str, phi, g, nm) -> float:
    pr = profile(name)
    if pr["bc"] == "wall":
        return 0.0
    centre = None if pr["bc"] == "semi_infinite" else 0.5 * (pr["domain"][0] + pr["domain"][1])
    return float(ST.far_field_fraction(phi, g.y, nm["y_max"], centre=centre))


def mode_at(name: str, k: float, Re: float, c_guess=None, N: int | None = None, lead_localised: bool = False) -> dict:
    """mode nearest to c_guess (or the leading eigenvalue; ``lead_localised``: the least-damped eigenvalue whose
    eigenfunction is localised), normalised as os_mode, with its energy budget and flag"""
    c, V, g, nm = spectrum(name, k, Re, N)
    if lead_localised:
        j = next((i for i in range(len(c)) if far_of(name, V[:, i], g, nm) < FAR_TOL), 0)
    else:
        j = 0 if c_guess is None else int(np.argmin(np.abs(c - c_guess)))
    phi = V[:, j]
    dphi = g.D1 @ phi
    sc = dphi[int(np.argmax(np.abs(dphi)))]
    phi, dphi = phi / sc, dphi / sc
    pr = profile(name)
    b = ch11.disturbance_energy_budget(float(k), complex(c[j]), phi, g.y, pr["Up"], float(Re), grid=g)
    far = far_of(name, phi, g, nm)
    if far >= FAR_TOL:
        flag = "x"
    else:
        jl = next((i for i in range(len(c)) if far_of(name, V[:, i], g, nm) < FAR_TOL), j)  # least-damped localised
        flag = "1" if abs(c[jl] - c[j]) < 1e-9 else "0"
    return dict(c=complex(c[j]), lead=complex(c[0]), y=g.y, phi=phi, u_hat=dphi, budget=b, flag=flag, far=far)


def lead_discrete(name: str, k: float, Re: float, N: int | None = None) -> dict:
    """ch11.os_leading_mode in this script's record format (flag "1"; "x" with the localised leading eigenvalue of the
    box if fluidpy finds no discrete mode)"""
    lm = ch11.os_leading_mode(profile(name), float(k), float(Re), N or NTAB[name])
    if not lm["discrete"]:
        m = mode_at(name, k, Re, N=N, lead_localised=True)
        m["flag"] = "x"
        return m
    m = lm["mode"]
    return dict(c=complex(lm["c"]), lead=complex(lm["c"]), y=m["y"], phi=m["phi"], u_hat=m["u_hat"], budget=lm["budget"],
                flag="1", far=lm["far"])


def ysample(name: str) -> np.ndarray:
    s = np.linspace(-1, 1, NY)
    if name in ("poiseuille", "couette"):
        return -np.cos(np.pi * np.arange(NY) / (NY - 1))
    if name == "blasius":
        return 8 * np.linspace(0, 1, NY) ** 2
    return 6 * np.sign(s) * np.abs(s) ** 1.5


def sample_mode(name: str, m: dict) -> dict:
    ys = ysample(name)
    if profile(name)["bc"] == "wall":
        f = lambda a: BarycentricInterpolator(m["y"], a)(ys)  # noqa: E731
    else:
        o = np.argsort(m["y"])
        f = lambda a: np.interp(ys, m["y"][o], a[o])  # noqa: E731
    d = dict(pr=f(m["phi"].real), pi=f(m["phi"].imag), ur=f(m["u_hat"].real), ui=f(m["u_hat"].imag))
    return {a: [sig(x, 3) for x in v] for a, v in d.items()}


def track(name: str, Re, k, seed, c_seed=None):
    """continuation over the grid, most reliable steps first (smallest miss of the extrapolated guess)"""
    nR, nk = len(Re), len(k)
    M = [[None] * nk for _ in range(nR)]
    M[seed[0]][seed[1]] = mode_at(name, k[seed[1]], Re[seed[0]], c_seed)
    heap, cnt = [], 0

    def push(i, j):
        nonlocal cnt
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            a, b = i + di, j + dj
            if 0 <= a < nR and 0 <= b < nk and M[a][b] is None:
                p, q = i - di, j - dj
                cg = M[i][j]["c"]
                if 0 <= p < nR and 0 <= q < nk and M[p][q] is not None:
                    cg = 2 * M[i][j]["c"] - M[p][q]["c"]
                heapq.heappush(heap, (M[i][j].get("step", 0.0), cnt, a, b, cg))
                cnt += 1

    push(*seed)
    while heap:
        _, _, a, b, cg = heapq.heappop(heap)
        if M[a][b] is not None:
            continue
        m = mode_at(name, k[b], Re[a], cg)
        m["step"] = abs(m["c"] - cg)
        M[a][b] = m
        push(a, b)
    return M


def walk(name: str, M, Re, k, R: float, kk: float, N: int | None = None, steps=(0.5, 1.0)) -> dict:
    """the followed mode at an off-grid point: continuation from the nearest node"""
    i = int(np.argmin(np.abs(np.log(Re) - np.log(R))))
    j = int(np.argmin(np.abs(np.log(k) - np.log(kk))))
    cg, m = M[i][j]["c"], None
    for f in steps:
        m = mode_at(name, np.exp(np.log(k[j]) + f * (np.log(kk) - np.log(k[j]))), np.exp(np.log(Re[i]) + f * (np.log(R) - np.log(Re[i]))),
                    cg, N=N if f == steps[-1] else None)
        cg = m["c"]
    return m


def bickley_bands() -> dict:
    """main band, stable gap and long-wave band of the sinuous Bickley jet from ch11.bickley_neutral_curve"""
    nc = ch11.bickley_neutral_curve()
    ex = ch11.bickley_neutral_curve([17.6, 18.2])
    Re = np.concatenate([nc["Re"], ex["Re"]])
    o = np.argsort(Re)
    col = lambda key: np.concatenate([np.real(nc[key]), np.real(ex[key])])[o]  # noqa: E731
    Re, ku, kl, kg = Re[o], col("k_upper"), col("k_lower"), col("k_long_upper")
    kmin = float(nc["k_min"])
    gap = [dict(Re=GAP_NOSE[0], top=GAP_NOSE[1], bot=GAP_NOSE[1])]
    main_kl = []
    for R, l, g in zip(Re, kl, kg):
        if R < GAP_NOSE[0]:
            main_kl.append(sig(l))  # None: the single band reaches below k_min
            continue
        main_kl.append(None)  # from the gap on, the band's floor is described by the gap rows
        if np.isfinite(l):
            gap.append(dict(Re=sig(R), top=sig(l), bot=sig(g)))  # bot None: the gap's lower edge is below k_min
        elif gap and gap[-1]["top"] is not None and gap[-1]["top"] > kmin:
            # the gap's upper edge leaves the resolved range: log-log extrapolation of its last two rows to k_min
            a, b = gap[-2], gap[-1]
            slope = np.log(b["top"] / a["top"]) / np.log(b["Re"] / a["Re"])
            R_end = float(b["Re"] * (kmin / b["top"]) ** (1.0 / slope)) if slope < 0 else float(R)
            gap.append(dict(Re=sig(min(R_end, R)), top=kmin, bot=None))
    return dict(nRe=[sig(x) for x in Re], nkl=main_kl, nku=[sig(x) for x in ku], gap=gap, k_min=kmin)


def build_flow(name: str, crit: dict) -> dict:
    t0 = time.time()
    ref = None
    if name == "couette":
        Re, k = np.geomspace(400.0, 4e4, 14), np.linspace(0.2, 2.0, 19)
        neutral = dict(nRe=[], nkl=[], nku=[])
        M = []
        for R in Re:
            row = []
            for kk in k:
                c0 = spectrum(name, kk, R)[0]
                ok = np.where(c0.real >= -1e-7)[0]
                row.append(mode_at(name, kk, R, c0[ok[np.argmax(c0[ok].imag)]]))
            M.append(row)
    else:
        g = load_csv(REF / f"os_grid_{name}.csv")
        nR, nk = len(np.unique(g[:, 0])), len(np.unique(g[:, 1]))
        ref = dict(k=np.unique(g[:, 1]), cr=g[:, 2].reshape(nR, nk), ci=g[:, 3].reshape(nR, nk), P=g[:, 4].reshape(nR, nk))
        lo, hi = dict(poiseuille=(2000.0, 1e6), blasius=(200.0, 6000.0), tanh=(1.0, 400.0), bickley=(2.0, 1000.0))[name]
        Re = np.geomspace(lo, hi, nR)  # the csv rounds Re to 4 s.f.
        ref["Re"] = Re
        if name == "bickley":
            Re = np.geomspace(lo, hi, 36)
            k = np.geomspace(0.02, 1.8, 37)  # log-spaced: the stable gap and the long-wave band live at k < 0.08
            neutral = bickley_bands()
            M = [[lead_discrete(name, kk, R) for kk in k] for R in Re]
        else:
            k = ref["k"]
            n = load_csv(REF / f"os_neutral_{name}.csv")
            neutral = dict(nRe=n[:, 0].tolist(), nkl=[sig(x) for x in n[:, 1]], nku=[sig(x) for x in n[:, 2]])
            seed = tuple(int(x) for x in np.unravel_index(int(np.argmax(ref["ci"])), ref["ci"].shape))
            M = track(name, Re, k, seed)
        if name == "poiseuille":  # above Re ~ 2e5 the N = 80 budget no longer closes to 0.2 %
            Re, M = Re[:18], M[:18]
            ref = {key: (v[:18] if key != "k" else v) for key, v in ref.items()}
    A = lambda f: np.array([[f(m) for m in row] for row in M])  # noqa: E731
    cr, ci = A(lambda m: m["c"].real), A(lambda m: m["c"].imag)
    P, L, E = A(lambda m: m["budget"]["production"]), A(lambda m: m["budget"]["dissipation"]), A(lambda m: m["budget"]["E"])
    lead_ci = A(lambda m: m["lead"].imag)
    flags = [[m["flag"] for m in row] for row in M]
    fl = np.array(flags)
    resid = A(lambda m: abs(m["budget"]["residual"]) / max(abs(m["budget"]["production"]), m["budget"]["dissipation"]))
    steps = A(lambda m: m.get("step", 0.0))
    miss = (lead_ci > 0) & (ci <= 0)
    print(f"{name}: {cr.size} nodes in {time.time() - t0:.0f} s; flags 1/0/x = {(fl == '1').sum()}/{(fl == '0').sum()}/{(fl == 'x').sum()}; "
          f"max continuation miss {steps.max():.3g}; budget residual max {resid.max():.2e}")
    print(f"   unstable nodes: followed {int((ci > 0).sum())}, any eigenvalue {int((lead_ci > 0).sum())}; "
          f"unstable through ANOTHER mode at {int(miss.sum())} nodes"
          + ("" if not miss.any() else ": " + ", ".join(f"(Re {Re[i]:.4g}, k {k[j]:.4g})" for i, j in np.argwhere(miss)[:12])))
    if ref is not None and name != "bickley":
        compare(name, ref, cr, ci)
    if ref is not None and name == "bickley":
        kk = ref["k"]
        cr2, ci2 = np.zeros_like(ref["cr"]), np.zeros_like(ref["ci"])
        for i, R in enumerate(ref["Re"]):
            for j, kv in enumerate(kk):
                m = lead_discrete(name, kv, R)
                cr2[i, j], ci2[i, j] = m["c"].real, m["c"].imag
        compare(name, ref, cr2, ci2)
    mi = [int(round(x)) for x in np.linspace(0, len(Re) - 1, 5)]
    mj = [int(round(x)) for x in np.linspace(0, len(k) - 1, 4)]
    d = dict(Re=[sig(x, 6) for x in Re], k=[sig(x, 6) for x in k], logk=name == "bickley", bilin=name == "bickley",
             cr=[[sig(x) for x in r] for r in cr], ci=[[sig(x) for x in r] for r in ci], P=[[sig(x) for x in r] for r in P],
             L=[[sig(x) for x in r] for r in L], E=[[sig(x) for x in r] for r in E], lead=["".join(r) for r in flags],
             miss=[[sig(Re[i], 4), sig(k[j], 4)] for i, j in np.argwhere(miss)],
             mRe=[sig(Re[i], 6) for i in mi], mk=[sig(k[j], 6) for j in mj], my=[sig(x, 4) for x in ysample(name)],
             modes=[sample_mode(name, M[i][j]) for i in mi for j in mj], **neutral)
    pts = dict(
        poiseuille=[(1e4, 1.0), (crit["poiseuille"]["Re_c"], crit["poiseuille"]["k_c"]), (5000.0, 1.0), (8000.0, 1.2), (6000.0, 0.9)],
        blasius=[(crit["blasius"]["Re_c"], crit["blasius"]["k_c"]), (1000.0, 0.25)],
        tanh=[(10.0, 0.5), (50.0, 0.45)],
        bickley=[(crit["bickley_sinuous"]["Re_c"], crit["bickley_sinuous"]["k_c"]), (20.0, 0.3)],
        couette=[(1e4, 1.0), (1e3, 1.0)],
    )[name]
    ex = []
    for R, kk in pts:
        m = lead_discrete(name, kk, R, N=80) if name == "bickley" else \
            walk(name, M, Re, k, R, kk, N=100, steps=(0.25, 0.5, 0.75, 1.0))
        b = m["budget"]
        ex.append(dict(Re=sig(R, 8), k=sig(kk, 6), cr=sig(m["c"].real, 6), ci=sig(m["c"].imag, 4), P=sig(b["production"], 4),
                       L=sig(b["dissipation"], 4), E=sig(b["E"], 4), lead=m["flag"], **sample_mode(name, m)))
        print(f"   preset {R:.7g} {kk:.6g}: c = {m['c']:.8f}  P/L = {b['ratio']:.5f}  flag {m['flag']}")
    d["exact"] = ex
    return d


def compare(name: str, ref: dict, cr: np.ndarray, ci: np.ndarray) -> None:
    """followed mode against reference/ch11/os_grid_<name>.csv (least-damped discrete mode) on the shared nodes"""
    disc = np.isfinite(ref["P"])
    d = np.hypot(cr - ref["cr"], ci - ref["ci"])
    same = disc & (d < 6e-4 * np.maximum(1.0, np.hypot(ref["cr"], ref["ci"])))  # the csv holds 4 significant figures
    un = disc & (ref["ci"] > 0)
    other = disc & ~same
    print(f"   vs reference grid ({disc.sum()} nodes with a discrete mode, {(~disc).sum()} without): same mode at {same.sum()} "
          f"({100 * same.sum() / max(1, disc.sum()):.0f} %), max |dc| there {d[same].max() if same.any() else 0:.1e}; "
          f"on the reference's {un.sum()} unstable nodes: same at {(un & same).sum()}, max |dc| {d[un].max() if un.any() else 0:.1e}")
    if other.any():
        more = other & (ci < ref["ci"] - 1e-6)
        print(f"   different mode at {other.sum()} nodes: followed mode MORE damped than the reference's at {more.sum()}, "
              f"less damped at {(other & ~more).sum()} (c_r followed {cr[other].min():.3g}..{cr[other].max():.3g}, "
              f"reference {ref['cr'][other].min():.3g}..{ref['cr'][other].max():.3g})")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--flows", nargs="*", default=None, help="recompute only these flows (others from their cached piece)")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    crit = json.loads((REF / "critical_points.json").read_text(encoding="utf-8"))
    data = {}
    for name in FLOWS:
        piece = OUT / f"orr_sommerfeld_neutral_curve_{name}.json"
        if (a.flows is not None and name not in a.flows and piece.exists()):
            data[name] = json.loads(piece.read_text(encoding="utf-8"))
            continue
        data[name] = build_flow(name, crit)
        piece.write_text(json.dumps(data[name], separators=(",", ":")), encoding="utf-8")
    data["crit"] = dict(poiseuille=crit["poiseuille"], blasius=crit["blasius"], bickley=crit["bickley_sinuous"],
                        tanh_inviscid=crit["tanh_inviscid"])
    pr = profile("blasius")
    y = np.linspace(0.0, 8.0, 41)
    data["blasius_profile"] = dict(y=[sig(v) for v in y], U=[sig(v, 5) for v in pr["U"](y)], Up=[sig(v, 5) for v in pr["Up"](y)])
    p = OUT / "orr_sommerfeld_neutral_curve_table.json"
    p.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    print("wrote", p.relative_to(ROOT), p.stat().st_size, "bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
