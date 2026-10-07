"""§12.3 (C01): ensemble, time and volume averages — our analogues of Figs. 12.2–12.3, the N^{-1/2} scatter of an N-member
mean, the commutation rules (12.4)–(12.9), the product split, and the two window factors of Example 12.1.

Run: ``.venv/Scripts/python.exe scripts/ch12_averaging.py --no-show [--fast]``
Figures -> outputs/ch12/c01_ensemble_average.png, c01_window_factors.png.
Our signal: mean 1.5 exp(-t/4) [m/s], fluctuation sigma = 0.4 m/s with memory tau_c = 0.15 s, seed 0.
"""
from __future__ import annotations

import numpy as np

from ch12_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch12_turbulence as ch12
from fluidpy.core import turbstats as TS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    t = np.linspace(0.0, 10.0, 1001)
    mean = lambda tt: 1.5 * np.exp(-tt / 4.0)  # noqa: E731
    n_big = 512 if args.fast else 4096
    ens = TS.make_ensemble(n_big, t, mean, sigma=0.4, tau_c=0.15, seed=0)
    stat = TS.make_ensemble(1, t, 1.0, sigma=0.4, tau_c=0.15, seed=1)[0]

    print("ensemble mean against the expected value (rms error over the record) and the standard error:")
    Ns = [2, 4, 8, 64, 512] + ([] if args.fast else [4096])
    errs = []
    for N in Ns:
        err = float(np.sqrt(np.mean((TS.ensemble_average(ens[:N]) - mean(t)) ** 2)))
        errs.append(err)
        print(f"  N = {N:5d}: rms error {err:.4f}   sigma/sqrt(N) = {TS.standard_error_of_mean(0.4, N):.4f}")
    slope = np.polyfit(np.log(Ns), np.log(errs), 1)[0]
    print(f"  slope of log(error) against log(N): {slope:.3f}  (expected -0.5)")

    s = TS.statistics(ens[:, 500], normalized=True)
    print(f"statistics at t = 5 s over {n_big} members: mean {s['mean']:.4f} (expected {mean(5.0):.4f}), variance {s['variance']:.4f} "
          f"(0.16), normalised skewness {s['skewness']:.3f} (0), kurtosis {s['kurtosis']:.3f} (3)")
    rules = TS.check_averaging_rules(ens[:64], ens[64:128], t)
    print("commutation rules (12.4)-(12.9), max residual:", {k: f"{v:.1e}" for k, v in rules.items()})
    mp, pm, cov = TS.product_average_split(ens[:, 500], ens[:, 520])
    print(f"product split at one time pair: mean(uv) = {mp:.5f} = {pm:.5f} + {cov:.5f} (product of means + covariance)")

    win = 1.0
    ta = TS.time_average(t, ens[3], win)
    print(f"sliding time average of member 4, window {win} s: rms difference from the expected value "
          f"{np.sqrt(np.nanmean((ta - mean(t)) ** 2)):.4f} (NaN within {win / 2} s of the ends)")

    # Example 12.1 with our numbers
    A, Bc, tau, om = 1.5, 0.4, 4.0, 2.0 * np.pi / 0.5
    sig = A * np.exp(-t / tau) + Bc * np.cos(om * t)
    print("Example 12.1 form, A = 1.5, B = 0.4, tau = 4 s, period 0.5 s:")
    for w in (0.1, 0.25, 0.5, 2.0):
        cf, mf, ff = ch12.time_average_exp_cos(t, A, Bc, tau, om, w)
        num = TS.time_average(t, sig, w)
        print(f"  window {w:4.2f} s: mean factor {mf:.5f}, fluctuation factor {ff:+.5f}; closed form vs sampled signal "
              f"max diff {np.nanmax(np.abs(cf - num)):.1e}")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for k in range(4):
        ax[0].plot(t, ens[k] + 2.0 * k, lw=0.8, color=COLORS["muted"])
    ax[0].plot(t, stat + 8.0, lw=0.8, color=COLORS["blue"])
    ax[0].set(xlabel="t [s]", ylabel="u [m/s] (members offset by 2)", title="four members of a decaying ensemble; top: a stationary record")
    for N, c in zip((2, 8, 64), (COLORS["muted"], COLORS["orange"], COLORS["accent"])):
        ax[1].plot(t, TS.ensemble_average(ens[:N]), lw=1.0, color=c, label=f"N = {N}")
    ax[1].plot(t, mean(t), "k--", lw=1.2, label="expected value")
    ax[1].plot(t, ta, color=COLORS["teal"], lw=1.4, label="time average of one member, 1 s window")
    ax[1].set(xlabel="t [s]", ylabel="mean of u [m/s]", title="ensemble mean (12.10) and sliding time mean (12.2)")
    ax[1].legend(fontsize=8)
    ax[2].loglog(Ns, errs, "o-", color=COLORS["accent"], label="measured rms error")
    ax[2].loglog(Ns, 0.4 / np.sqrt(Ns), "--", color=COLORS["muted"], label=r"$\sigma/\sqrt{N}$")
    ax[2].set(xlabel="members N", ylabel="error of the N-member mean [m/s]", title=f"slope {slope:.2f}")
    ax[2].legend(fontsize=8)
    save(fig, out, "c01_ensemble_average")

    w = np.linspace(1e-3, 4.0, 400)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))
    mf = [ch12.time_average_exp_cos(0.0, A, Bc, tau, om, x)[1] for x in w]
    ff = [ch12.time_average_exp_cos(0.0, A, Bc, tau, om, x)[2] for x in w]
    ax[0].plot(w, mf, color=COLORS["accent"], label=r"mean factor $\sinh(\Delta t/2\tau)/(\Delta t/2\tau)$")
    ax[0].plot(w, ff, color=COLORS["teal"], label=r"fluctuation factor $\sin(\omega\Delta t/2)/(\omega\Delta t/2)$")
    ax[0].axhline(0, color=COLORS["grid"])
    ax[0].set(xlabel=r"window $\Delta t$ [s]", ylabel="factor [-]", title="Example 12.1: what a window does to each part")
    ax[0].legend(fontsize=8)
    for wv, c in ((0.2, COLORS["orange"]), (0.5, COLORS["accent"]), (3.0, COLORS["rose"])):
        ax[1].plot(t, TS.time_average(t, sig, wv), color=c, label=f"window {wv} s")
    ax[1].plot(t, sig, color=COLORS["grid"], lw=0.8, label="signal")
    ax[1].plot(t, A * np.exp(-t / tau), "k--", lw=1.0, label="decaying mean")
    ax[1].set(xlabel="t [s]", ylabel="u [m/s]", title="too short, one whole period, too long")
    ax[1].legend(fontsize=8)
    save(fig, out, "c01_window_factors")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
