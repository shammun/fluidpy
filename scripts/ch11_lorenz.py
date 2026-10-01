"""§11.14 (C15): deterministic chaos — pendulum phase portrait (11.89) with energy conservation, Hopf attractors and the
bifurcation diagram (Fig. 11.28 analogue), the Lorenz system (11.91): fixed points, Jacobian eigenvalues, r_H, X(t) (Fig. 11.29),
the attractor (Fig. 11.30), sensitivity to initial conditions, the RK4 twin of the explainer, the logistic map's period-doubling
tree (Fig. 11.31) and Feigenbaum's ratio.

Run: ``.venv/Scripts/python.exe scripts/ch11_lorenz.py --no-show [--fast]``   Figures -> outputs/ch11/c15_attractors.png,
c15_lorenz.png.  Lorenz's parameters Pr = 10, b = 8/3, r = 28 (Lorenz 1963).
"""
from __future__ import annotations

import numpy as np

from ch11_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    starts = [(0.5, 0.0), (1.5, 0.0), (2.8, 0.0), (-3.1, 1.0), (3.1, 0.5)]
    pp = ch11.phase_portrait(ch11.pendulum_rhs, starts, t_end=20.0, n=800)
    drift = max(float(np.max(np.abs(ch11.pendulum_energy((d["X"], d["Y"])) - ch11.pendulum_energy((d["X"][0], d["Y"][0])))))
                for d in pp)
    print(f"pendulum (11.89): max energy drift over 5 trajectories = {drift:.1e}")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for d in pp:
        ax[0].plot(d["X"], d["Y"], color=COLORS["accent"], lw=1)
    ax[0].set(xlabel="X [rad]", ylabel="Y = Ẋ [rad/s]", title="pendulum phase portrait (11.89)")
    for mu, col in ((-0.3, COLORS["teal"]), (0.5, COLORS["rose"])):
        for d in ch11.phase_portrait(ch11.hopf_normal_form, [(0.05, 0.0), (1.5, 0.0)], t_end=30.0, n=1500, mu=mu):
            ax[1].plot(d["X"], d["Y"], color=col, lw=0.8)
    ax[1].set(xlabel="x", ylabel="y", title="Hopf: point attractor (μ < 0), limit cycle (μ > 0)")
    ax[1].set_aspect("equal")
    mus = np.linspace(-1, 1, 200)
    amp = ch11.limit_cycle_amplitude(mus)
    ax[2].plot(mus, amp, color=COLORS["rose"], label="stable limit cycle ±√μ")
    ax[2].plot(mus, -amp, color=COLORS["rose"])
    ax[2].plot(mus[mus < 0], 0 * mus[mus < 0], color=COLORS["teal"], label="stable fixed point")
    ax[2].plot(mus[mus >= 0], 0 * mus[mus >= 0], color=COLORS["muted"], ls="--", label="unstable fixed point")
    ax[2].set(xlabel="μ ~ R − R_cr", ylabel="extremum x", title="bifurcation diagram (cf. Fig. 11.28c)")
    ax[2].legend(fontsize=8)
    save(fig, out, "c15_attractors")

    print(f"b(k² = π²/2) = {ch11.lorenz_b(np.pi / np.sqrt(2)):.10f};  r(Ra = 27π⁴/4, k = π/√2) = "
          f"{ch11.lorenz_r(ch11.RA_FREE_FREE, np.pi / np.sqrt(2)):.10f}")
    print(f"fixed points r = 28: {ch11.lorenz_fixed_points(28.0)}")
    print(f"eigenvalues at C±: {np.round(ch11.lorenz_eigs(28.0), 4)}; at the origin: {np.round(ch11.lorenz_eigs(28.0, which='O'), 4)}")
    rH = ch11.lorenz_hopf_r()
    print(f"r_H = {rH:.6f}; eigenvalues at C± there: {np.round(ch11.lorenz_eigs(rH), 6)}; ∇·ṡ = {ch11.lorenz_divergence():.4f}")
    print(f"RK4 one step dt = 0.01 from (1, 1, 1): {np.round(ch11.lorenz_rk4((1, 1, 1), 0.01, 1)[1], 6)}")
    rk = ch11.lorenz_rk4((1, 1, 1), 0.005, 1000)
    ad = ch11.lorenz_integrate(t_end=5.0, n=1001)
    print(f"RK4 (dt = 0.005) vs DOP853 at t = 5: |Δ| = {np.linalg.norm(rk[-1] - np.array([ad['X'][-1], ad['Y'][-1], ad['Z'][-1]])):.2e}")
    with Timer("Lorenz run"):
        L = ch11.lorenz_integrate(t_end=25.0 if args.fast else 50.0, n=5001)
    sep = ch11.lorenz_separation(t_end=40.0)
    print(f"separation of runs 1e-8 apart: slope of ln|δ| over t ∈ {tuple(round(v, 2) for v in sep['window'])} = "
          f"{sep['slope']:.3f} (≈ largest Lyapunov exponent ~0.9, qualitative); |δ| > 1 at t = {ch11.lorenz_predictability_time():.2f}")
    if not args.fast:
        with Timer("Lyapunov (cached)"):
            print(f"largest Lyapunov exponent (renormalisation, t = 200): {ch11.lorenz_largest_lyapunov():.4f}")
    with Timer("period doubling"):
        pdp = ch11.period_doubling_points(6)
    print(f"logistic period doublings A_n = {np.round(pdp['A_n'], 7)}; superstable S_n = {np.round(pdp['S_n'], 7)}")
    print(f"Feigenbaum ratios (S_n): {np.round(ch11.feigenbaum_estimate(8), 5)} → δ = 4.669201609")
    fig = plt.figure(figsize=(16, 9))
    a1 = fig.add_subplot(2, 3, (1, 2))
    a1.plot(L["t"], L["X"], color=COLORS["accent"], lw=0.8)
    for p in ch11.lorenz_fixed_points(28.0)[1:]:
        a1.axhline(p[0], color=COLORS["muted"], ls="--", lw=1)
    a1.set(xlabel="t", ylabel="X(t)", title="Lorenz X(t), r = 28 (cf. Fig. 11.29)")
    a2 = fig.add_subplot(2, 3, 3, projection="3d")
    a2.plot(L["X"], L["Y"], L["Z"], lw=0.4, color=COLORS["accent"])
    a2.set(title="attractor (cf. Fig. 11.30)", xlabel="X", ylabel="Y", zlabel="Z")
    a3 = fig.add_subplot(2, 3, 4)
    a3.semilogy(sep["t"], sep["sep"], color=COLORS["rose"])
    tt = np.linspace(*sep["window"], 50)
    a3.semilogy(tt, np.exp(sep["slope"] * (tt - tt[0])) * np.interp(tt[0], sep["t"], sep["sep"]), color=COLORS["ink"], ls="--",
                label=f"slope {sep['slope']:.2f}")
    a3.set(xlabel="t", ylabel="|δ|", title="two runs 10⁻⁸ apart")
    a3.legend(fontsize=8)
    a4 = fig.add_subplot(2, 3, (5, 6))
    bd = ch11.bifurcation_diagram(np.linspace(2.8, 4.0, 900 if not args.fast else 300), 500, 80)
    a4.plot(bd["A"], bd["x"], ",", color=COLORS["ink"], alpha=0.5)
    for A in pdp["A_n"][:5]:
        a4.axvline(A, color=COLORS["rose"], lw=0.8)
    a4.set(xlabel="A", ylabel="x", title="logistic map: period doubling (cf. Fig. 11.31)")
    save(fig, out, "c15_lorenz")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
