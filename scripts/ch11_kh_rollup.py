"""§11.3 (C02, A2): the vortex-sheet roll-up (Fig. 11.6 analogue) from Ch. 5's point-vortex model against the linear growth
e^{kΔU t/2} of (11.20).

Run: ``.venv/Scripts/python.exe scripts/ch11_kh_rollup.py --no-show [--fast]``   Figure -> outputs/ch11/c02_kh_rollup.png.
Our inputs: period L = 1, γ = ΔU = 1, initial amplitude 0.01, Krasny δ = 0.1, N = 200 (fast: 100) vortices.
"""
from __future__ import annotations

import numpy as np

from ch11_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch05_vorticity_dynamics as ch05
from fluidpy import ch11_instability as ch11


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    N = 100 if args.fast else 200
    A0, L, gam = 0.01, 1.0, 1.0
    t = np.linspace(0.0, 2.0, 9)
    with Timer(f"sheet_rollup N = {N}"):
        r = ch05.sheet_rollup(N=N, gamma=gam, amplitude=A0, delta=0.1, t_eval=t, L=L)
    k = 2 * np.pi / L
    c = ch11.vortex_sheet_c(-gam / 2, gam / 2)[0]  # (11.20): c = (U₁ + U₂)/2 + i|U₂ − U₁|/2
    growth = k * c.imag
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
    for j in range(0, len(t), 2):
        ax[0].plot(r["x"][j] % L, r["y"][j] + 0.25 * j / 2, ".", ms=2, color=COLORS["accent"])
    ax[0].set(xlabel="x/L", ylabel="y/L (offset per frame)", title="nonlinear roll-up of a perturbed vortex sheet (Ch. 5)")
    amp = 0.5 * (r["y"].max(axis=1) - r["y"].min(axis=1))
    ax[1].semilogy(t, amp, "o-", color=COLORS["ink"], label="point-vortex sheet: half peak-to-peak")
    ax[1].semilogy(t, A0 * np.exp(growth * t), color=COLORS["rose"],
                   label=f"linear (11.20), singular sheet: A₀e^{{kc_i t}}, kc_i = {growth:.4f}")
    ax[1].set(xlabel="t γ/L", ylabel="amplitude / L", title="linear growth, then saturation")
    ax[1].legend(fontsize=8)
    save(fig, out, "c02_kh_rollup")
    print(f"linear growth rate kc_i = k|ΔU|/2 = {growth:.6f} (π for L = 1, γ = 1)")
    for tt, a in zip(t, amp):
        print(f"  t = {tt:4.2f}: amplitude {a:.5f}   linear {A0 * np.exp(growth * tt):.5f}")
    fit = np.polyfit(t[1:4], np.log(amp[1:4]), 1)[0]
    print(f"fitted early growth rate (t = 0.25–0.75) = {fit:.3f}: below kc_i = π because Krasny's δ = {0.1} regularises the "
          "sheet (finite thickness slows the growth) and a pure displacement excites the decaying root too")
    print(f"mean y conserved: max drift {np.max(np.abs(r['mean_y'] - r['mean_y'][0])):.2e}")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
