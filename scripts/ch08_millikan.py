"""§8.6 (C14, Fig. 8.18): settling by Stokes' law — terminal velocities of cloud droplets, sediment and aerosols with
their Re (validity), and a synthetic Millikan oil-drop experiment (seeded) recovering the elementary charge.

Run: ``.venv/Scripts/python.exe scripts/ch08_millikan.py --no-show``
Figure → outputs/ch08/c14_millikan.png.
"""
from __future__ import annotations

import numpy as np

from ch08_common import COLORS, finish, parse_args, save, setup

from fluidpy import ch08_laminar_flow as ch08


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    cases = [("cloud droplet 10 μm in air", 10e-6, 1000.0, 1.2, 1.8e-5), ("drizzle 100 μm in air", 100e-6, 1000.0, 1.2, 1.8e-5),
             ("clay 2 μm in water", 2e-6, 2650.0, 1000.0, 1e-3), ("silt 30 μm in water", 30e-6, 2650.0, 1000.0, 1e-3)]
    for name, a, rp, rf, mu in cases:
        s = ch08.settling_state(a, rp, rf, mu)
        print(f"{name}: U_t = {s['U_t']:.4e} m/s, Re = {s['Re']:.3g} ({'valid' if s['valid'] else 'outside Stokes range'}), "
              f"D = {s['D']:.3e} N (pressure {s['D_pressure'] / s['D']:.3f}, friction {s['D_friction'] / s['D']:.3f})")
    m = ch08.synthetic_millikan(40, seed=0, noise=0.01)
    print(f"synthetic Millikan (40 drops, 1 % noise, seed 0): e = {m['e_est']:.6e} C, relative error "
          f"{m['rel_error']:+.2e}; integers recovered correctly: {np.mean(m['n_est'] == m['n_true']) * 100:.0f} %")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    ax[0].plot(m["q"] / 1e-19, "o", color=COLORS["accent"])
    for k in range(1, 8):
        ax[0].axhline(k * m["e_est"] / 1e-19, color=COLORS["grid"], lw=0.8)
    ax[0].set_xlabel("drop")
    ax[0].set_ylabel("q [10⁻¹⁹ C]")
    ax[0].set_title("charges cluster on multiples of e")
    ax[1].hist(m["q"] / m["e_est"], bins=60, color=COLORS["teal"])
    ax[1].set_xlabel("q / e_est")
    ax[1].set_ylabel("count")
    ax[1].set_title(f"e = {m['e_est']:.4e} C ({m['rel_error'] * 100:+.2f} %)")
    save(fig, out, "c14_millikan")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
