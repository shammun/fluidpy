"""Shared drawing helpers for the ch06 scripts (flow nets of plane and axisymmetric flows with bodies and stagnation
points). Not a physics module: every number comes from ``fluidpy.ch06_ideal_flow``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fluidpy.core.style import COLORS  # noqa: E402


def parse_args(doc: str):
    """The common command line of every ch06 script: --out, --no-show, --fast."""
    ap = argparse.ArgumentParser(description=doc.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch06"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true")
    return ap.parse_args()


def setup(args):
    """Select a non-interactive backend when not showing, apply the house style, create the output folder."""
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    return out


def save(fig, out: Path, name: str) -> Path:
    path = out / f"{name}.png"
    fig.savefig(path, bbox_inches="tight")
    print(f"  saved {path}")
    return path


def flow_net(ax, flow, xlim=(-3, 3), ylim=(-2, 2), n: int = 241, n_levels: int = 21, mask_body=None, *,
             phi: bool = False, body=None, stag=None, title: str = "", psi_levels=None, cut_ray: bool = False):
    """ψ contours (ink) — and φ contours (dashed teal) if ``phi`` — of a plane flow, the body filled grey with a black
    outline (:func:`draw_body`), stagnation points as orange dots. ``mask_body(x, y) → bool`` blanks extra points (NaN);
    ``body``: complex outline; ``stag``: complex array. Automatic levels span the 3rd–97th percentiles (a singularity
    would swamp them); ``cut_ray`` blanks a thin strip along the positive x-axis where a logarithm's branch cut would
    otherwise draw a bundle of lines."""
    x = np.linspace(*xlim, n)
    y = np.linspace(*ylim, n)
    X, Y = np.meshgrid(x, y, indexing="xy")
    with np.errstate(all="ignore"):
        P = np.asarray(flow.psi(X, Y), float)
    if mask_body is not None:
        P = np.where(np.asarray(mask_body(X, Y), bool), np.nan, P)
    if cut_ray:
        strip = (np.abs(Y) < 1.5 * (y[1] - y[0])) & (X > 0)
        P = np.where(strip, np.nan, P)

    def levels(A):
        f = A[np.isfinite(A)]
        lo, hi = np.percentile(f, [3, 97])
        return np.linspace(lo, hi, n_levels) if hi > lo else n_levels

    lv = psi_levels if psi_levels is not None else levels(P)
    ax.contour(X, Y, P, levels=lv, colors=COLORS["ink"], linewidths=0.8, negative_linestyles="solid")
    if phi:
        with np.errstate(all="ignore"):
            F = np.asarray(flow.phi(X, Y), float)
        if cut_ray:
            F = np.where(strip, np.nan, F)
        ax.contour(X, Y, F, levels=levels(F), colors=COLORS["teal"], linewidths=0.6, linestyles="--")
    if body is not None:
        draw_body(ax, body)
    if stag is not None and len(stag):
        s = np.asarray(stag, complex)
        ax.plot(s.real, s.imag, "o", color=COLORS["orange"], ms=6, zorder=5, label="stagnation point")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    if title:
        ax.set_title(title, fontsize=10)
    return ax


def circle(a: float = 1.0, n: int = 200, center: complex = 0j) -> np.ndarray:
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return center + a * np.exp(1j * th)


def draw_body(ax, pts, **kw):
    """Fill a body outline (complex array, or (2, N)) light grey and draw its edge black; ``kw`` override
    ``facecolor``, ``edgecolor``, ``lw``, ``zorder``. Returns the fill patch list."""
    b = np.asarray(pts)
    b = b if np.iscomplexobj(b) else np.asarray(b[0], float) + 1j * np.asarray(b[1], float)
    fc = kw.pop("facecolor", "#d9dce4")
    ec = kw.pop("edgecolor", COLORS["ink"])
    lw = kw.pop("lw", 1.6)
    z = kw.pop("zorder", 3)
    patches = ax.fill(b.real, b.imag, color=fc, zorder=z, **kw)
    ax.plot(np.append(b.real, b.real[0]), np.append(b.imag, b.imag[0]), color=ec, lw=lw, zorder=z + 1)
    return patches


def pressure_arrows(ax, pts, p, normals, scale: float = 1.0, **kw):
    """Surface pressure as arrows along −p n (the force per area the fluid exerts on the body): gauge pressure
    ``p`` [Pa] at body points ``pts`` (complex) with outward unit ``normals`` (complex); arrow length |p|·scale [m/Pa].
    Over-pressure (p > 0, pushing in) orange, suction (p < 0, pulling out) blue. Returns the quiver."""
    z = np.asarray(pts, complex)
    nn = np.asarray(normals, complex)
    pv = np.asarray(p, float)
    f = -pv * nn * float(scale)  # arrow vector −p n
    tail = z - np.where(pv > 0, f, 0.0)  # over-pressure arrows end on the surface, suction arrows start there
    cols = np.where(pv > 0, COLORS["orange"], COLORS["blue"])
    return ax.quiver(tail.real, tail.imag, f.real, f.imag, color=list(cols), angles="xy", scale_units="xy",
                     scale=1.0, width=kw.pop("width", 0.004), zorder=kw.pop("zorder", 6), **kw)


def separated_cp_band(theta_front_deg, sep_deg: float = 80.0):
    """**Qualitative sketch band, not measured data** (N31, the idea of Fig. 6.10): where a real high-Re cylinder's
    surface C_p lies, against the angle from the front stagnation point [deg]. Before separation (angle < ``sep_deg``)
    it follows the ideal 1 − 4 sin²β (6.35) (floored at −1.05, so the two parts join) within ±0.15; after separation
    it is a flat, low wake pressure between −1.3 and −0.8. Every figure using it must say "qualitative".
    Returns (lo, hi) arrays (lo ≤ hi)."""
    b = np.radians(np.asarray(theta_front_deg, float))
    ideal = np.maximum(1.0 - 4.0 * np.sin(b) ** 2, -1.05)
    before = b < np.radians(sep_deg)
    lo = np.where(before, ideal - 0.15, -1.3)
    hi = np.where(before, np.minimum(ideal + 0.15, 1.0), -0.8)
    return np.minimum(lo, hi), np.maximum(lo, hi)
