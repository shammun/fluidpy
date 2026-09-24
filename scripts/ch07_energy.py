"""§7.2 (C06) and §7.5 (C09): wave energy per unit horizontal area — kinetic (7.38)–(7.39) and potential (7.40)–(7.41)
by direct double integrals against the closed forms (equipartition at every depth), E = ½ρga² (7.42), and the energy
flux (7.43)–(7.44) = E c_g (7.71); the swell worked number (a = 1 m, T = 10 s, deep water).

Run: ``.venv/Scripts/python.exe scripts/ch07_energy.py --no-show``
Figures → outputs/ch07/c06_energy.png.
"""
from __future__ import annotations

import numpy as np

from ch07_drawings import parse_args, save, setup

from fluidpy import ch07_gravity_waves as ch07
from fluidpy.core.style import COLORS


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    a, k, rho = 1.0, 1.0, 1000.0
    kHs = [0.1, 0.3, 1.0, 3.0, 10.0, np.inf]
    Ek, Ep, F, Fc = [], [], [], []
    for kH in kHs:
        H = kH / k
        q = ch07.wave_energy(a, k, H, rho=rho, method="quad")
        Ek.append(q["Ek"])
        Ep.append(q["Ep"])
        F.append(ch07.energy_flux(a, k, H, rho=rho, method="quad"))
        Fc.append(ch07.energy_flux(a, k, H, rho=rho))
        print(f"kH = {kH:>5}: E_k = {q['Ek']:.6f}, E_p = {q['Ep']:.6f} J/m² (closed ¼ρga² = {0.25 * rho * 9.81:.6f});"
              f" F quad = {F[-1]:.6f}, (7.44) = {Fc[-1]:.6f}, E c_g = "
              f"{ch07.wave_energy_density(a, rho) * ch07.group_velocity(k, H):.6f} W/m")
        assert abs(q["Ek"] / q["Ep"] - 1) < 1e-8 and abs(F[-1] / Fc[-1] - 1) < 1e-8
    ps = ch07.packet_state(ch07.wavenumber_from_omega(2 * np.pi / 10.0), a=1.0)
    print(f"swell a = 1 m, T = 10 s, deep: E = {ps['E'] / 1e3:.2f} kJ/m², c = {ps['c']:.2f} m/s, c_g = {ps['cg']:.2f} m/s,"
          f" F = E c_g = {ps['F'] / 1e3:.1f} kW per metre of crest; energy arrives from 1000 km in {ps['arrival_h']:.1f} h")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    lbl = [("∞" if np.isinf(v) else f"{v:g}") for v in kHs]
    xx = np.arange(len(kHs))
    ax[0].bar(xx - 0.2, Ek, 0.4, color=COLORS["teal"], label="kinetic E_k (7.38), dblquad")
    ax[0].bar(xx + 0.2, Ep, 0.4, color=COLORS["amber"], label="potential E_p (7.40), dblquad")
    ax[0].axhline(0.25 * rho * 9.81 * a ** 2, color=COLORS["ink"], ls="--", lw=1, label="¼ρga² (7.39), (7.41)")
    ax[0].set_xticks(xx, lbl)
    ax[0].set_xlabel("kH")
    ax[0].set_ylabel("energy per unit area [J/m²] (a = 1 m)")
    ax[0].set_title("equipartition at every depth", fontsize=10)
    ax[0].legend(fontsize=8, loc="lower right")
    kk = np.logspace(-2, 1.5, 300)
    c = ch07.phase_speed(kk, 1.0)
    cg = ch07.group_velocity(kk, 1.0)
    ax[1].semilogx(kk, cg / c, color=COLORS["accent"], label="F/(Ec) = c_g/c = ½[1 + 2kH/sinh 2kH] (7.44), (7.69)")
    ax[1].axhline(0.5, color=COLORS["muted"], ls="--", lw=1, label="deep ½ (7.70)")
    ax[1].axhline(1.0, color=COLORS["muted"], ls=":", lw=1, label="shallow 1 (7.70)")
    ax[1].set_xlabel("kH")
    ax[1].set_ylabel("energy speed / phase speed")
    ax[1].set_title("energy travels at c_g, not c", fontsize=10)
    ax[1].legend(fontsize=8)
    save(fig, out, "c06_energy")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
