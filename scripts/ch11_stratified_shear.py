"""§11.7 (C08–C10): stratified shear flows — the Taylor–Goldstein growth map kc_i(k, J) of U = tanh z, N² = J sech²z with the
exact neutral curve J = k(1 − k), Ri(z) and the Miles–Howard bound (11.67), the identities (11.65), (11.69)–(11.70), Howard's
semicircle with computed eigenvalues, and a random-profile check that no mode grows when Ri_min > ¼.

Run: ``.venv/Scripts/python.exe scripts/ch11_stratified_shear.py --no-show [--fast]``   Figure -> outputs/ch11/c08_stratified.png.
The full growth map is cached (outputs/ch11/cache); the first full run takes several minutes.
"""
from __future__ import annotations

import numpy as np

from ch11_common import COLORS, Timer, finish, parse_args, save, setup

from fluidpy import ch11_instability as ch11


def main() -> int:
    args = parse_args(__doc__)
    out = setup(args)
    import matplotlib.pyplot as plt

    with Timer("TG growth map"):
        gm = ch11.tg_growth_map(fast=args.fast)
    j0 = int(np.argmin(np.abs(gm["J"])))
    print(f"J = 0 row: max kc_i = {np.max(gm['kci'][j0]):.4f} at k = {gm['k'][np.argmax(gm['kci'][j0])]:.3f} (Rayleigh: 0.1897 at 0.4449)")
    unstable_beyond = gm["kci"][gm["J"] > 0.25].max() if np.any(gm["J"] > 0.25) else 0.0
    print(f"largest kc_i for J > 1/4 (must be 0, Miles–Howard): {unstable_beyond:.3e}")
    print(f"tg_growth(0.4, 0.1) = {ch11.tg_growth(0.4, 0.1):.5f};  exact neutral J at k = 0.3, 0.5, 0.7: "
          f"{[float(ch11.tg_tanh_neutral_J(k)) for k in (0.3, 0.5, 0.7)]}")
    for k, J in ((0.3, 0.20), (0.3, 0.21), (0.5, 0.24), (0.5, 0.25)):
        print(f"  k = {k}, J = {J}: kc_i = {ch11.tg_growth(k, J):.5f}   (neutral J = {k * (1 - k):.3f})")
    pr = ch11.richardson_profiles("tanh", 0.1)
    res = ch11.taylor_goldstein_eigs(0.4, pr["U"], pr["Upp"], pr["N2"], domain=(-1, 1), N=100, bc="decay", return_vectors=True)
    c, psi, g = res["c"][0], res["psi"][:, 0], res["grid"]
    ri = ch11.richardson_identity_check(0.4, c, psi, g.y, pr["U"], pr["Up"], pr["Upp"], pr["N2"], grid=g)
    hw = ch11.howard_identity_check(0.4, c, psi, g.y, pr["U"], pr["N2"], grid=g, Up=pr["Up"])
    eb = ch11.stratified_energy_budget(0.4, c, psi, g.y, pr["U"], pr["N2"], grid=g, Up=pr["Up"])
    print(f"mode J = 0.1, k = 0.4: c = {c:.6f}; (11.65) residual {ri['residual']:.1e} (Im: {ri['imag_lhs']:.6f} vs "
          f"{ri['imag_rhs']:.6f}); (11.69) {hw['res_69']:.1e}, (11.70) {hw['res_70']:.1e}, c_r = ∫UQ/∫Q = {hw['cr_mean']:.2e}; "
          f"in semicircle {hw['in_semicircle']}; Ex. 11.10 budget residual {eb['residual']:.1e} (KE {eb['KE']:.4f}, APE {eb['APE']:.4f})")
    mh = ch11.miles_howard_stable(np.linspace(-5, 5, 1001), pr["U"], pr["N2"], dUdz=pr["Up"])
    print("Miles–Howard J = 0.1:", mh["text"])
    rng = np.random.default_rng(0)
    worst = 0.0
    n_cases = 10 if args.fast else 30
    with Timer(f"{n_cases} random profiles with Ri_min > 1/4"):
        for _ in range(n_cases):
            a, w = rng.uniform(0.5, 2.0), rng.uniform(0.3, 1.5)
            U = lambda z, a=a: np.tanh(a * z)  # noqa: E731
            Up = lambda z, a=a: a / np.cosh(a * z) ** 2  # noqa: E731
            Upp = lambda z, a=a: -2 * a ** 2 * np.tanh(a * z) / np.cosh(a * z) ** 2  # noqa: E731
            zz = np.linspace(-6, 6, 2001)
            Rimin_target = rng.uniform(0.26, 0.6)
            N2 = lambda z, a=a, w=w, R=Rimin_target: R * a ** 2 / np.cosh(a * z) ** 4 * (1 + w * np.tanh(z) ** 2)  # noqa: E731
            assert ch11.miles_howard_stable(zz, U, N2, dUdz=Up)["guaranteed_stable"]
            for k in (0.2, 0.5, 0.8):
                cc = ch11.taylor_goldstein_eigs(k, U, Upp, N2, domain=(-1, 1), N=60, bc="decay",
                                                y_max=ch11.decay_box(k))  # box scales with 1/k (e^{−k|z|} decay)
                if len(cc):
                    worst = max(worst, float(cc[0].imag))
    print(f"largest converged c_i over random Ri_min > 1/4 profiles: {worst:.2e} (must be 0)")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    cs = ax[0].contourf(gm["k"], gm["J"], gm["kci"], 20, cmap="Reds")
    kk = np.linspace(0, 1, 200)
    ax[0].plot(kk, kk * (1 - kk), color=COLORS["accent"], lw=2, label="exact neutral J = k(1 − k)")
    ax[0].axhline(0.25, color=COLORS["ink"], ls="--", lw=1, label="J = 1/4 (Miles–Howard)")
    ax[0].set(xlabel="k", ylabel="J = Ri(0)", title="TG growth kc_i: U = tanh z, N² = J sech² z")
    ax[0].legend(fontsize=8)
    fig.colorbar(cs, ax=ax[0])
    z = np.linspace(-3, 3, 400)
    for J, col in ((0.1, COLORS["rose"]), (0.25, COLORS["accent"]), (0.4, COLORS["teal"])):
        ax[1].plot(ch11.gradient_richardson(z, np.tanh, ch11.richardson_profiles("tanh", J)["N2"]), z, color=col, label=f"J = {J}")
    ax[1].axvline(0.25, color=COLORS["ink"], ls="--", lw=1)
    ax[1].set(xlim=(0, 2), xlabel="Ri(z) (11.66)", ylabel="z", title="gradient Richardson number")
    ax[1].legend(fontsize=8)
    cr, ci = ch11.howard_semicircle(-1.0, 1.0)
    ax[2].fill(cr, ci, color=COLORS["accent"], alpha=0.15)
    n_semi = 0
    for J in (0.0, 0.1, 0.2):
        p = ch11.richardson_profiles("tanh", J)
        for k in (0.2, 0.4, 0.6, 0.8):
            # box and node clustering both follow the wave, as tg_growth does (near-neutral modes need the clustering)
            cc = ch11.taylor_goldstein_eigs(k, p["U"], p["Upp"], p["N2"], domain=(-1, 1), N=80, bc="decay",
                                            y_max=ch11.decay_box(k), map_scale=ch11.decay_map_scale(k))
            n_semi += len(cc)
            assert np.all(ch11.in_howard_semicircle(cc, -1.0, 1.0)), (J, k, cc)
            ax[2].plot(np.real(cc), np.imag(cc), "o", color=COLORS["rose"], ms=4)
    print(f"semicircle panel: {n_semi} unstable eigenvalues plotted (J = 0, 0.1, 0.2; k = 0.2 … 0.8), all inside (11.72)")
    for kw, what in ((dict(dUdz=0.0, N2=-1.0), "U' = 0, N² < 0"), (dict(dUdz=0.0, N2=0.0), "U' = 0, N² = 0"),
                     (dict(dUdz=0.0, N2=1.0), "U' = 0, N² > 0")):
        v = ch11.miles_howard_stable(0.0, **kw)
        print(f"  shear-free level, {what}: Ri = {ch11.gradient_richardson(0.0, **kw)}, guaranteed stable = "
              f"{v['guaranteed_stable']}")
    ax[2].set(xlabel="c_r", ylabel="c_i", title="unstable c inside Howard's semicircle (p. 507)")
    ax[2].set_aspect("equal")
    save(fig, out, "c08_stratified")
    finish(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
