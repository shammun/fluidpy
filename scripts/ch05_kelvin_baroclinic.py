"""§5.2–5.3 Kelvin's circulation theorem, the baroclinic torque and Helmholtz's first theorem (C03, C04, C05, N10–N17,
N48–N51): material loops in E3's scenarios with Γ(t) against the prediction, the split (5.9)/(5.10) of DΓ/Dt, a
hand-written RK4 loop vs ``material_circulation``, the lock exchange's baroclinic field and rate (D07), the pressure
torque on a shrinking disc converging to (∇ρ × ∇p)/ρ² at order 2 (D06), the Kelvin decision table and the frozen-in
check.

Run: ``.venv/Scripts/python.exe scripts/ch05_kelvin_baroclinic.py --no-show``
Figures → outputs/ch05/c03_material_loops.png, c04_lock_exchange.png, c04_torque_convergence.png.
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fluidpy import ch05_vorticity_dynamics as ch05  # noqa: E402
from fluidpy.core.style import COLORS  # noqa: E402
from tools.convergence import observed_order  # noqa: E402


def rk4_loop_circulation(u, pts0, t_end: float, nsteps: int) -> float:
    """From scratch: advect the loop points with classical RK4 and take ∮u·dx by the periodic trapezoid rule."""
    P = np.array(pts0, float)
    dt = t_end / nsteps
    t = 0.0
    for _ in range(nsteps):
        k1 = u(P, t)
        k2 = u(P + 0.5 * dt * k1, t + 0.5 * dt)
        k3 = u(P + 0.5 * dt * k2, t + 0.5 * dt)
        k4 = u(P + dt * k3, t + dt)
        P = P + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
    N = P.shape[1]
    k = 2 * np.pi * np.fft.fftfreq(N, d=1.0 / N)
    dPds = np.real(np.fft.ifft(1j * k * np.fft.fft(P, axis=1), axis=1))  # spectral x′(s)
    return float(np.sum(u(P, t) * dPds) / N)


def loops_figure(fast: bool):
    """Material loops at three times and Γ(t) for four scenarios."""
    import matplotlib.pyplot as plt

    names = ("cellular", "rankine_straddle", "lamb_oseen", "rotating")
    fig, ax = plt.subplots(2, 4, figsize=(16, 7))
    for j, nm in enumerate(names):
        s = ch05.kelvin_scenario(nm)
        te = np.linspace(0.0, s["t_end"], 31 if fast else 61)
        G, loops = ch05.material_circulation(s["u"], s["pts0"], te, return_loops=True)
        for i, c in zip((0, len(te) // 2, -1), (COLORS["muted"], COLORS["accent"], COLORS["orange"])):
            L = np.hstack([loops[i], loops[i][:, :1]])
            ax[0, j].plot(L[0], L[1], color=c, lw=1.2, label=f"t = {te[i]:.1f} s")
        ax[0, j].set_aspect("equal")
        ax[0, j].set_title(nm, fontsize=10)
        ax[0, j].legend(fontsize=7)
        ax[1, j].plot(te, G, color=COLORS["accent"], label="material loop Γ(t)")
        ax[1, j].plot(te, s["Gamma_exact"](te), "--", color=COLORS["muted"], label="prediction")
        if s["Gamma_a_exact"] is not None:
            Ga = [ch05.absolute_circulation(s["u"], L_, s["Omega"], tt)[1] for L_, tt in zip(loops, te)]
            ax[1, j].plot(te, Ga, color=COLORS["amber"], label="Γ_a (5.33)")
        ax[1, j].set_xlabel("t [s]")
        ax[1, j].set_ylabel("Γ [m²/s]")
        ax[1, j].legend(fontsize=7)
    return fig


def lock_figure():
    """ρ and the baroclinic term (∇ρ × ∇p)_z/ρ² of the lock exchange at t = 0⁺."""
    import matplotlib.pyplot as plt

    le = ch05.lock_exchange_fields(nx=256, ny=64)
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    c0 = ax[0].pcolormesh(le["X"], le["Y"], le["rho_grid"], cmap="Blues", shading="auto")
    fig.colorbar(c0, ax=ax[0], label="ρ [kg/m³]")
    c1 = ax[1].pcolormesh(le["X"], le["Y"], le["baroclinic_z"], cmap="Oranges", shading="auto")
    fig.colorbar(c1, ax=ax[1], label="(∇ρ × ∇p)_z/ρ² [1/s²]")
    for a in ax:
        a.set_xlabel("x [m]")
        a.set_ylabel("y [m]")
    ax[0].set_title("heavy (ρ₂) left, light (ρ₁) right")
    ax[1].set_title(f"counterclockwise spin-up, peak {le['baroclinic_z'].max():.3f} s⁻²")
    return fig


def torque_convergence(fast: bool):
    """ratio − 1 of 2·torque/I_G to (∇ρ × ∇p)/ρ² vs disc radius, for a nonlinear ρ (order 2 expected)."""
    import matplotlib.pyplot as plt

    rho = lambda x, y: 1000.0 + 5.0 * x - 3.0 * y + 40.0 * x * y + 30.0 * x ** 2  # noqa: E731
    p = lambda x, y: 1e5 - 9810.0 * y + 20.0 * x + 50.0 * x ** 2 * y  # noqa: E731
    R = np.array([0.1, 0.05, 0.025, 0.0125]) if not fast else np.array([0.1, 0.05, 0.025])
    err = np.array([abs(ch05.pressure_torque_on_element(p, rho, (0.2, 0.1), r, n=2000, nr=64)["ratio"] - 1) for r in R])
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.loglog(R, err, "o-", color=COLORS["accent"], label=f"observed order {observed_order(R, err):.2f}")
    ax.loglog(R, err[0] * (R / R[0]) ** 2, "--", color=COLORS["muted"], label="slope 2")
    ax.set_xlabel("disc radius R [m]")
    ax.set_ylabel("|2·torque/I_G ÷ baroclinic − 1|")
    ax.legend()
    return fig, observed_order(R, err)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch05"))
    ap.add_argument("--no-show", action="store_true")
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for nm in ch05.KELVIN_SCENARIOS:
        s = ch05.kelvin_scenario(nm)
        te = np.linspace(0.0, s["t_end"], 5)
        G, loops = ch05.material_circulation(s["u"], s["pts0"], te, return_loops=True)
        dev = np.abs(G - s["Gamma_exact"](te)).max() if nm != "baroclinic" else float("nan")
        grow = ch05.loop_length(loops[-1]) / ch05.loop_length(loops[0])
        print(f"{nm:17s} Γ(0) = {G[0]:+.7f} m²/s, max |Γ − prediction| = {dev:.1e}, loop length ×{grow:.2f}"
              f"  [{s['label']}]")
    s = ch05.kelvin_scenario("cellular")
    g_rk4 = rk4_loop_circulation(s["u"], s["pts0"], s["t_end"], 600 if args.fast else 1500)
    g_lib = ch05.kelvin_scenario_circulation("cellular", s["t_end"])
    print(f"from scratch (RK4 loop) Γ(t_end) = {g_rk4:.9f}, library {g_lib:.9f}, difference {abs(g_rk4 - g_lib):.1e}")
    for nm in ("cellular", "lamb_oseen", "baroclinic", "rotating"):
        r = ch05.kelvin_scenario_rate(nm, 0.0 if nm != "rotating" else 2.0)
        print(f"(5.9)/(5.10) {nm:10s}: " + ", ".join(f"{k} {v:+.3e}" for k, v in r.items()))
    G, dG = ch05.lamb_oseen_circulation(0.005, 10.0, 0.01, 1e-6)
    print(f"(5.11) Lamb–Oseen r = 5 mm, t = 10 s: Γ = {G:.8f}, ∂Γ/∂t = {dG:.6e}, ∮ν∇²u·dx = "
          f"{ch05.lamb_oseen_viscous_loop_integral(0.005, 10.0, 0.01, 1e-6):.6e} m²/s²")
    le = ch05.lock_exchange_fields()
    bz = ch05.baroclinic_term(le["rho_fn"], le["p_fn"], np.array([0.0, 0.5]), 0.0, 1e-5)[2]
    print(f"lock exchange 1000/1025 kg/m³, δ = 0.1 m: closed form {le['rate']:.6f} s⁻², field route {bz:.6f} s⁻², "
          f"∫ baroclinic dA = {le['circulation_rate']:.6f} m²/s² (≈ gH ln(ρ₂/ρ₁) = {9.81 * np.log(1.025):.6f}; "
          f"mean-density hydrostatic model p = ρ̄g(H − y))")
    pt = ch05.pressure_torque_on_element()
    print(f"pressure torque (R = 1 cm, ∇ρ = (10, 0), ∇p = (0, −9810)): x_G = {pt['x_G'][0]:.3e} m, torque "
          f"{pt['torque']:.5e} N m/m, I_G {pt['I_G']:.5e}, spin-up {pt['spin_up']:.10f} vs baroclinic "
          f"{pt['baroclinic']:.10f} s⁻² (ratio − 1 = {pt['ratio'] - 1:.3e})")
    n_hold = sum(ch05.kelvin_hypotheses(*c)["holds"] for c in itertools.product((True, False), repeat=4))
    print(f"Kelvin decision table: {n_hold} of 16 combinations keep DΓ/Dt = 0; e.g. "
          f"'{ch05.kelvin_hypotheses_text(True, False, True, False)}'")
    fr = ch05.frozen_in_check(ch05.abc_flow_field(), (0.3, 0.2, 0.1), t_span=(0.0, 3.0))
    fb = ch05.frozen_in_check(ch05.burgers_vortex_field(1.0, 1.0, 0.05), (0.3, 0.2, 0.1), t_span=(0.0, 2.0), nu=0.05)
    print(f"Helmholtz 1 (ABC, inviscid): max angle {fr['angle'].max():.1e} rad, |ω|/|δx| drift "
          f"{np.abs(fr['ratio_field'] - 1).max():.1e}; Burgers (ν > 0): |ω|/|δx| falls to {fb['ratio_field'][-1]:.3f}")
    fig = loops_figure(args.fast)
    fig.savefig(out / "c03_material_loops.png", bbox_inches="tight")
    lock_figure().savefig(out / "c04_lock_exchange.png", bbox_inches="tight")
    fig, q = torque_convergence(args.fast)
    fig.savefig(out / "c04_torque_convergence.png", bbox_inches="tight")
    print(f"pressure-torque ratio → 1 with observed order {q:.2f} in R")
    print(f"figures → {out}")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
