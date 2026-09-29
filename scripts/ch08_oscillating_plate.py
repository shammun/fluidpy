"""§8.5 (C11, Fig. 8.16): Stokes' second problem — profiles (8.38) at ωt = 0, π/2, π, 3π/2 inside the envelope
±U e^{−y/δ_e}, the book's δ ~ 4√(ν/ω) vs δ_e = √(2ν/ω), the crest speed √(2νω) measured by tracking a zero crossing,
the Crank–Nicolson solution from rest converging to (8.38) after the transients, and the derivation check in sympy.

Run: ``.venv/Scripts/python.exe scripts/ch08_oscillating_plate.py --no-show [--fast]``
Figure → outputs/ch08/c11_oscillating_plate.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08
from fluidpy.core import diffusion as DIF


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    U, omega, nu = 1.0, 2 * np.pi, 1e-6  # 1 Hz, water
    s = ch08.stokes_layer(nu, omega)
    print("Stokes layer: " + ", ".join(f"{k} = {v:.5g}" for k, v in s.items()))
    de = s["delta_e"]
    y = np.linspace(0, 8 * de, 400)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for ph, c in zip((0, 0.5 * np.pi, np.pi, 1.5 * np.pi), (COLORS["accent"], COLORS["teal"], COLORS["orange"],
                                                             COLORS["blue"])):
        ax[0].plot(ch08.stokes_second_problem(y, ph / omega, U, omega, nu), y / de, lw=2, color=c,
                   label=f"ωt = {ph / np.pi:.1f}π")
    env = ch08.stokes_layer_envelope(y, U, omega, nu)
    ax[0].plot(env, y / de, "k--", lw=0.8)
    ax[0].plot(-env, y / de, "k--", lw=0.8, label="±U e^{−y/δ_e}")
    ax[0].axhline(s["delta_book"] / de, color=COLORS["muted"], ls=":", label="book δ ~ 4√(ν/ω)")
    ax[0].set_xlabel("u/U")
    ax[0].set_ylabel("y/δ_e,  δ_e = √(2ν/ω)")
    ax[0].legend(fontsize=7)
    ax[0].set_title("(8.38) at four phases (Fig. 8.16)")
    # crest speed from tracking the first zero crossing
    tt = np.linspace(0.3, 0.45, 60)
    yf = np.linspace(0, 6 * de, 20001)
    zc = []
    for t in tt:
        u = ch08.stokes_second_problem(yf, t, U, omega, nu)
        i = np.nonzero(np.sign(u[:-1]) != np.sign(u[1:]))[0][0]
        zc.append(yf[i] - u[i] * (yf[i + 1] - yf[i]) / (u[i + 1] - u[i]))
    c_meas = np.polyfit(tt, zc, 1)[0]
    print(f"crest speed from zero-crossing tracking {c_meas:.6e} m/s vs √(2νω) = {s['phase_speed']:.6e} m/s")
    st = ch08.stokes_layer_state(nu, omega, s["delta_book"])
    print("state at y = 4√(ν/ω): " + ", ".join(f"{k} = {v:.4g}" for k, v in st.items()))
    # CN from rest
    yg = np.linspace(0, 12 * de, 481)
    nper = 6 if args.fast else 10
    nsteps = 200 * nper
    t_arr, UU = DIF.crank_nicolson_1d(np.zeros_like(yg), yg, nper * 2 * np.pi / omega / nsteps, nsteps, nu,
                                      lambda t: U * np.cos(omega * t), 0.0, return_all=True, save_every=200)
    err = [np.max(np.abs(UU[k] - ch08.stokes_second_problem(yg, t_arr[k], U, omega, nu))[yg < 6 * de])
           for k in range(len(t_arr))]
    print("CN max error per period (y < 6δ_e): " + ", ".join(f"{e:.2e}" for e in err))
    ax[1].semilogy(t_arr * omega / (2 * np.pi), np.maximum(err, 1e-16), "o-", color=COLORS["accent"])
    ax[1].set_xlabel("periods since start")
    ax[1].set_ylabel("max |u_CN − (8.38)|/U")
    ax[1].set_title("transients die out: CN from rest → (8.38)")
    ax[2].plot(UU[-1], yg / de, color=COLORS["orange"], lw=2, label="Crank–Nicolson")
    ax[2].plot(ch08.stokes_second_problem(yg, t_arr[-1], U, omega, nu), yg / de, "k--", lw=1, label="(8.38)")
    ax[2].set_ylim(0, 8)
    ax[2].set_xlabel("u/U")
    ax[2].set_ylabel("y/δ_e")
    ax[2].legend(fontsize=8)
    ax[2].set_title(f"after {nper} periods")
    d = ch08.stokes_second_sympy()
    print(f"sympy: roots {d['k_roots']} match ±(1 + i)√(ω/2ν): {d['roots_match']}; residuals {d['residuals']}")
    save(fig, out, "c11_oscillating_plate")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
