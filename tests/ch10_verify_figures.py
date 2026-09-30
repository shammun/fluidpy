"""Verification figures for Chapter 10 (our code only) → outputs/ch10/verify/ (git-ignored).

Run: ``.venv/Scripts/python.exe tests/ch10_verify_figures.py``. Prints the numbers quoted in reports/ch10_verification.md and
writes one PNG per reproduced figure. No book figure is copied: every curve is computed by ``fluidpy`` with our own
parameters (figure numbers are the chapter's; the pictures are ours).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy import ch10_computational_fluid_dynamics as ch10  # noqa: E402
from fluidpy.core import fd as FD  # noqa: E402
from fluidpy.core import fem1d as FEM1  # noqa: E402
from fluidpy.core import fem2d as FEM2  # noqa: E402
from fluidpy.core import mac as MAC  # noqa: E402
from tools.convergence import grid_convergence_index, observed_order  # noqa: E402

OUT = ROOT / "outputs" / "ch10" / "verify"
OUT.mkdir(parents=True, exist_ok=True)
REF = ROOT / "reference" / "ch10"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=100)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def read_csv(name):
    lines = [ln for ln in (REF / name).read_text(encoding="utf-8").splitlines() if ln.strip() and not ln.startswith("#")]
    hdr = lines[0].split(",")
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return {h: data[:, k] for k, h in enumerate(hdr)}


def fig_stencils():
    hs = np.logspace(-8, -0.5, 60)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for kind in ("forward", "backward", "central", "central2", "onesided2"):
        e = [abs(FD.stencil_error("sin", 1.0, h, kind)) for h in hs]
        ax.loglog(hs, e, label=kind)
        k = (hs > 1e-3) & (hs < 0.1)
        print(f"  {kind:10s} order {observed_order(hs[k], np.array(e)[k]):.3f}")
    ax.loglog(hs, 0.5 * hs, "k:", lw=0.8, label="slope 1")
    ax.loglog(hs, 0.2 * hs ** 2, "k--", lw=0.8, label="slope 2")
    ax.set_xlabel("h"), ax.set_ylabel("|stencil − exact| for sin at x₀ = 1"), ax.set_title("(10.6)–(10.7): orders and round-off floor")
    ax.legend(fontsize=7)
    save(fig, "c01_stencil_orders")


def fig_von_neumann():
    th = np.linspace(0, np.pi, 361)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for a, b in ((0.1, 0.2), (0.0, 0.5), (0.0, 0.51), (0.4, 0.2), (0.3, 0.0)):
        ax[0].plot(th, FD.amplification_modulus(th, a, b), label=f"α={a}, β={b}")
    ax[0].axhline(1, color="k", lw=0.8)
    ax[0].set_xlabel("θ = kπΔx [rad]"), ax[0].set_ylabel("|G| (10.24)"), ax[0].legend(fontsize=7), ax[0].set_title("FTCS amplification")
    A = np.linspace(-0.6, 0.6, 121)
    B = np.linspace(0, 0.7, 141)
    S = np.array([[FD.is_von_neumann_stable("ftcs", a, b, n_theta=361) for a in A] for b in B])
    ax[1].contourf(A, B, S.astype(float), levels=[-0.5, 0.5, 1.5], colors=["#f4c7c3", "#c8e6c9"])
    ax[1].plot(A, 2 * A ** 2, "k--", lw=1, label="4α² = 2β (Noye)")
    ax[1].axhline(0.5, color="k", ls=":", lw=1, label="2β = 1")
    ax[1].set_xlabel("α"), ax[1].set_ylabel("β"), ax[1].set_ylim(0, 0.7), ax[1].legend(fontsize=7)
    ax[1].set_title("brute-force max|G| ≤ 1 (green) vs (10.27)")
    save(fig, "c04_von_neumann")


def fig_upwind_maccormack():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for k, prof in enumerate(("square", "gauss")):
        ex = FD.advect_periodic(prof, 0.8, 100, 1.0, "exact")
        ax[k].plot(ex["x"], ex["T0"], "k-", lw=2, alpha=0.3, label="exact (one revolution)")
        for sch, c in (("upwind", "tab:cyan"), ("maccormack", "tab:purple")):
            r = FD.advect_periodic(prof, 0.8, 100, 1.0, sch)
            ax[k].plot(r["x"], r["T"], color=c, label=f"{sch}: amplitude {r['amplitude_ratio']:.3f}")
        ax[k].set_xlabel("x"), ax[k].set_ylabel("T"), ax[k].legend(fontsize=7), ax[k].set_title(f"C = 0.8, N = 100, {prof}")
    save(fig, "c05_c10_upwind_maccormack")


def fig_cell_peclet():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    n, R = 10, 40.0  # R_cell = 4
    xx = np.linspace(0, 1, 400)
    ax[0].plot(xx, FD.steady_cd_exact(xx, R), "k-", alpha=0.4, lw=2, label="exact (10.86), R = 40")
    c = FD.steady_cd_fd(n, R)
    u = FD.steady_cd_fd(n, R, "upwind")
    ax[0].plot(c["x"], c["T"], "o-", color="tab:orange", label="centred (10.91): wiggles")
    ax[0].plot(u["x"], u["T"], "s-", color="tab:cyan", label="upwind (10.93)")
    Dn = FD.numerical_diffusivity(R, 1 / n, D=1.0)
    ax[0].plot(xx, FD.steady_cd_exact(xx, R / (1 + Dn)), "--", color="tab:cyan", label="exact of (10.94)")
    ax[0].set_xlabel("x/L"), ax[0].set_ylabel("T"), ax[0].set_title("R_cell = 4"), ax[0].legend(fontsize=7)
    Rc = np.linspace(0.05, 6, 300)
    ax[1].plot(Rc, [FD.discrete_root(r) if abs(r - 2) > 0.02 else np.nan for r in Rc], color="tab:orange", label="centred r")
    ax[1].plot(Rc, [FD.discrete_root(r, "upwind") for r in Rc], color="tab:cyan", label="upwind r = 1 + R_cell")
    ax[1].axhline(0, color="k", lw=0.8), ax[1].axvline(2, color="k", ls=":")
    ax[1].set_ylim(-10, 10), ax[1].set_xlabel("R_cell"), ax[1].set_ylabel("r"), ax[1].legend(fontsize=7)
    ax[1].set_title("D16: r < 0 ⇔ R_cell > 2")
    save(fig, "c09_cell_peclet")


def fig_fem():
    fig, ax = plt.subplots(figsize=(6, 4))
    xx = np.linspace(0, 1, 400)
    ax.plot(xx, FD.steady_cd_exact(xx, 8.0), "k-", alpha=0.4, lw=2, label="exact, R = 8")
    for n in (4, 8, 16):
        fe = FEM1.solve_steady(np.linspace(0, 1, n + 1), 8.0, 1.0, T_L=1.0)
        fdv = FD.steady_cd_fd(n, 8.0)
        ax.plot(fe["x"], fe["T"], "o-", ms=4, label=f"FE n={n}")
        ax.plot(fdv["x"], fdv["T"], "x", color="k", ms=6)
        print(f"  FE vs FD n = {n}: {np.max(np.abs(fe['T'] - fdv['T'])):.2e}")
    ax.set_xlabel("x"), ax.set_ylabel("T"), ax.legend(fontsize=7), ax.set_title("steady Galerkin FE (o) = centred FD (x) (N33)")
    save(fig, "c07_fe_equals_fd")


def fig_projection():
    g = MAC.MacGrid(16, 16)
    c = MAC.face_coordinates(g)
    us = np.sin(np.pi * c["xu"]) * np.cos(np.pi * c["yu"]) + 0.5 * np.sin(2 * np.pi * c["xu"])
    vs = np.cos(np.pi * c["xv"]) * np.sin(np.pi * c["yv"])
    st = MAC.projection_stages(us, vs, g, 1.0)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    im = ax[0].imshow(st["div_before"], origin="lower", extent=(0, 1, 0, 1))
    fig.colorbar(im, ax=ax[0]), ax[0].set_title(f"∇·u* (max {st['max_div_before']:.2f})")
    im = ax[1].imshow(st["p"], origin="lower", extent=(0, 1, 0, 1), cmap="RdBu_r")
    fig.colorbar(im, ax=ax[1]), ax[1].set_title("p from (10.124)")
    im = ax[2].imshow(np.log10(np.abs(st["div_after"]) + 1e-300), origin="lower", extent=(0, 1, 0, 1))
    fig.colorbar(im, ax=ax[2]), ax[2].set_title(f"log10|∇·u| after (max {st['max_div_after']:.1e})")
    for a in ax:
        a.set_xlabel("x"), a.set_ylabel("y")
    save(fig, "c12_projection")


def fig_infsup():
    t = read_csv("infsup_table.csv")
    fig, ax = plt.subplots(figsize=(6, 4))
    for pair, lab in ((0, "P1–P1 (smallest nonzero)"), (1, "P2–P1")):
        k = t["pair"] == pair
        ax.loglog(t["n"][k], t["beta_nonspurious"][k], "o-", label=lab)
    ax.set_xlabel("squares per side n"), ax.set_ylabel("inf–sup constant β_h"), ax.legend(fontsize=8)
    ax.set_title("LBB: Taylor–Hood bounded, equal order falls (+7 spurious modes)")
    save(fig, "c13_infsup")


def fig_cavity():
    g100 = ch10.ghia_centreline(100)
    g400 = ch10.ghia_centreline(400)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for n in (16, 32, 64, 128):
        d = read_csv(f"cavity_mac_re100_n{n}.csv")
        ax[0].plot(d["u"], d["y"], label=f"MAC {n}²")
    d = read_csv("cavity_mck_re100_n128.csv")
    ax[0].plot(d["u"], d["y"], "--", label="MacCormack 128²")
    ax[0].plot(g100["u"], g100["y"], "ko", ms=4, label="Ghia et al. 1982")
    ax[0].set_xlabel("u(x = ½, y)"), ax[0].set_ylabel("y"), ax[0].legend(fontsize=7), ax[0].set_title("Re = 100 (Fig. 10.8 analogue)")
    for n in (32, 64, 128):
        d = read_csv(f"cavity_mac_re400_n{n}.csv")
        ax[1].plot(d["u"], d["y"], label=f"MAC {n}²")
    d = read_csv("cavity_mck_re400_n64.csv")
    ax[1].plot(d["u"], d["y"], "--", label="MacCormack 64²")
    ax[1].plot(g400["u"], g400["y"], "ko", ms=4, label="Ghia et al. 1982")
    ax[1].set_xlabel("u(x = ½, y)"), ax[1].set_ylabel("y"), ax[1].legend(fontsize=7), ax[1].set_title("Re = 400")
    st = MAC.cavity(100.0, 64)
    g = st["g"]
    xs = np.arange(g.nx + 1) * g.dx
    ax[2].contour(xs, xs, st["psi"], levels=np.linspace(st["psi"].min(), 0, 14), colors="k", linewidths=0.7)
    ax[2].contour(xs, xs, st["psi"], levels=np.linspace(1e-6, max(st["psi"].max(), 2e-6), 4), colors="tab:red", linewidths=0.7)
    c = st["centre"]
    ax[2].plot(c["x"], c["y"], "b+", ms=12, label=f"ours ({c['x']:.4f}, {c['y']:.4f})")
    gc = ch10.ghia_vortex_centre(100)
    ax[2].plot(gc["x"], gc["y"], "gx", ms=10, label="Ghia 1982")
    hx, hy = ch10.hou_centres()[100]
    ax[2].plot(hx, hy, "m1", ms=12, label="Hou et al. 1995")
    ax[2].set_aspect("equal"), ax[2].legend(fontsize=7), ax[2].set_title("MAC 64², Re = 100 streamlines (Fig. 10.7 analogue)")
    save(fig, "c14_cavity")
    for n in (16, 24, 32, 64, 128):
        print(f"  MAC Re100 n={n}: max dev {ch10.cavity_error_vs_ghia(MAC.cavity(100.0, n), 100)['max_dev']:.5f}")


def fig_grid_convergence():
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    cd = []
    for n in (8, 16, 32):
        d = read_csv(f"block_forces_re20_dx{n}.csv")
        ax[0].plot(d["t"], d["CD"], label=f"Δx = 1/{n}")
        cd.append(d["CD"][d["t"] > 20].mean())
    ax[0].set_ylim(3.5, 5.5), ax[0].set_xlabel("t U/D"), ax[0].set_ylabel("C_D"), ax[0].legend(fontsize=7)
    ax[0].set_title("block, Re = 20 (our geometry)")
    gb = grid_convergence_index(cd[2], cd[1], cd[0], 2.0)
    ax[1].plot([1 / 8, 1 / 16, 1 / 32], cd, "o-", label=f"mean C_D, p = {gb['p']:.2f}")
    ax[1].axhline(gb["f_ext"], color="gray", ls=":", label=f"Richardson {gb['f_ext']:.3f}")
    ax[1].set_xlabel("Δx/D"), ax[1].set_ylabel("C_D"), ax[1].legend(fontsize=7), ax[1].set_title("Fig. 10.13 analogue (not asymptotic)")
    print(f"  block C_D {cd}, GCI {gb}")
    ps = []
    hs = []
    for n in (16, 32, 64, 128):
        txt = (REF / f"cavity_mac_re100_n{n}.csv").read_text(encoding="utf-8")
        ps.append(float(txt.split("psi_min = ")[1].split(";")[0]))
        hs.append(1 / n)
    gp = grid_convergence_index(ps[3], ps[2], ps[1], 2.0)
    ax[2].plot(hs, ps, "o-", label=f"ψ_min, p = {gp['p']:.2f}")
    ax[2].axhline(gp["f_ext"], color="gray", ls=":", label=f"Richardson {gp['f_ext']:.5f}")
    ax[2].set_xlabel("h"), ax[2].set_ylabel("ψ_min"), ax[2].legend(fontsize=7), ax[2].set_title("MAC cavity, Re = 100")
    print(f"  cavity psi_min {ps}, GCI {gp}")
    save(fig, "c15_grid_convergence")


def fig_cylinder():
    d = read_csv("cylinder_fe_re100_forces.csv")
    fig, ax = plt.subplots(figsize=(8, 3.8))
    ax.plot(d["t"], d["CD"], label="C_D")
    ax.plot(d["t"], d["CL"], label="C_L")
    ax.plot(d["t"], 100 * d["CM"], label="100 C_M")
    late = d["t"] > d["t"][-1] / 2
    St = ch10.dominant_frequency(d["t"][late], d["CL"][late])
    ax.set_xlabel("t U/d"), ax.legend(fontsize=7), ax.set_title(f"FE cylinder Re = 100 (our W = 5d): St = {St:.4f} (Fig. 10.21 analogue)")
    save(fig, "c13_cylinder_re100_forces")
    print(f"  cylinder St = {St:.4f}")


if __name__ == "__main__":
    for f in (fig_stencils, fig_von_neumann, fig_upwind_maccormack, fig_cell_peclet, fig_fem, fig_projection, fig_infsup,
              fig_cavity, fig_grid_convergence, fig_cylinder):
        f()
