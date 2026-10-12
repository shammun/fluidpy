"""Reproduce the key plots of chapter 14 with our own code, for the verifier's visual check.

Writes PNGs to ``outputs/ch14/verify/`` (git-ignored).  Every curve is computed by ``fluidpy``; nothing is traced from
the book.  Run: ``.venv/Scripts/python.exe tests/ch14_verify_figures.py``
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fluidpy import ch14_aerodynamics as ch14  # noqa: E402

OUT = ROOT / "outputs" / "ch14" / "verify"
I = ch14.illustrative_inputs()
PI = math.pi
SEC = (I["U_section"], I["b"], I["e"], I["beta"], I["alpha"])


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=110)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def fig_shapes():
    """Figs. 14.11–14.14: four circles and their images; the rule t_max/c = 1.3 e against the measured contour."""
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    cases = [("centred: flat plate", 0.0, 0.0), ("raised: circular arc", 0.0, 0.25), ("shifted: symmetric", 0.12, 0.0), ("both: cambered", 0.12, 0.2)]
    for (name, e, be), col in zip(cases, ("C0", "C1", "C2", "C3")):
        S = ch14.zhukhovsky_airfoil(1.0, e, be, 600)
        ax[0, 0].plot(np.append(S["zeta"].real, S["zeta"].real[0]), np.append(S["zeta"].imag, S["zeta"].imag[0]), col, label=name)
        ax[0, 1].plot(np.append(S["x"], S["x"][0]), np.append(S["y"], S["y"][0]), col, label=name)
    ax[0, 0].plot([1, -1], [0, 0], "k+", ms=10)
    ax[0, 0].set(aspect="equal", title="ζ-plane: circles through ζ = b", xlabel="ξ / b", ylabel="η / b")
    ax[0, 1].set(aspect="equal", title="z-plane: z = ζ + b²/ζ", xlabel="x / b", ylabel="y / b")
    ax[0, 0].legend(fontsize=8)
    arc = ch14.circular_arc(1.0, 0.25)
    S = ch14.zhukhovsky_airfoil(1.0, 0.0, 0.25, 600)
    tt = np.linspace(0, 2 * PI, 400)
    ax[0, 2].plot(arc["centre"].real + arc["radius"] * np.cos(tt), arc["centre"].imag + arc["radius"] * np.sin(tt), "0.7", label="circle x² + (y + 2b cot 2β)² = (2b csc 2β)²")
    ax[0, 2].plot(S["x"], S["y"], "C1", lw=3, label="image of the raised circle")
    ax[0, 2].axhline(arc["height"], color="k", ls=":", lw=1)
    ax[0, 2].set(aspect="equal", xlim=(-2.6, 2.6), ylim=(-1.0, 1.2), title=f"arc: height 2b tan β = {arc['height']:.3f} b", xlabel="x / b", ylabel="y / b")
    ax[0, 2].legend(fontsize=7, loc="lower center")
    es = np.linspace(0.005, 0.25, 26)
    meas = np.array([ch14.zhukhovsky_shape(1.0, float(e), 0.0, n=4001) for e in es])
    ax[1, 0].plot(es, [m["thickness_ratio"] for m in meas], "C2", label="measured on the contour")
    ax[1, 0].plot(es, ch14.THICKNESS_RULE_COEFFICIENT * es, "k--", label="rule (3√3/4) e ≈ 1.3 e")
    ax[1, 0].set(xlabel="e", ylabel="t_max / c", title="thickness rule of (14.8): first order in e")
    ax[1, 0].legend()
    ax[1, 1].plot(es, [m["x_tmax_over_c"] for m in meas], "C2", label="position of t_max")
    ax[1, 1].axhline(0.25, color="k", ls="--", label="quarter chord")
    ax[1, 1].plot(es, [m["chord"] / 4 for m in meas], "C4", label="chord / 4b")
    ax[1, 1].set(xlabel="e", ylim=(0.2, 1.1), title="where the section is thickest; chord against 4b")
    ax[1, 1].legend()
    bs = np.linspace(0.0, 0.4, 21)
    ax[1, 2].plot(bs, [ch14.zhukhovsky_shape(1.0, 0.0, float(b), n=4001)["camber_ratio"] for b in bs], "C1", label="measured (arc)")
    ax[1, 2].plot(bs, bs / 2, "k--", label="rule β/2")
    ax[1, 2].set(xlabel="β [rad]", ylabel="camber / c", title="camber rule of (14.8)")
    ax[1, 2].legend()
    save(fig, "fig_shapes")


def fig_kutta():
    """Fig. 14.9 / 14.16: three circulations; the edge speed; the lift curve (Fig. 14.17) with the generic stalled model."""
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.4))
    F0 = ch14.zhukhovsky_flow(*SEC)
    x = np.linspace(-0.9, 0.9, 260)
    y = np.linspace(-0.45, 0.45, 180)
    X, Y = np.meshgrid(x, y)
    S = ch14.zhukhovsky_airfoil(I["b"], I["e"], I["beta"], 400)
    for k, frac in enumerate((0.0, 0.5, 1.0)):
        F = ch14.zhukhovsky_flow(*SEC, Gamma_cw=frac * F0.Gamma_kutta)
        psi = F.stream(X, Y)
        ax[0, k].contour(X, Y, psi, levels=np.linspace(np.nanmin(psi), np.nanmax(psi), 41), colors="C0", linewidths=0.6)
        ax[0, k].fill(S["x"], S["y"], "0.85", ec="k")
        st = ch14.stagnation_points(*SEC, Gamma_cw=frac * F0.Gamma_kutta)
        ax[0, k].plot([st["z_rear"].real], [st["z_rear"].imag], "ro", ms=6)
        ax[0, k].plot([st["z_front"].real], [st["z_front"].imag], "go", ms=6)
        ax[0, k].set(aspect="equal", title=f"Γ = {frac:g} Γ_Kutta  (rear stagnation point red)", xlabel="x [m]", ylabel="y [m]")
    d = np.logspace(-7, -1.5, 40)
    for frac, col in ((1.0, "C2"), (0.8, "C3"), (1.2, "C1")):
        ax[1, 0].loglog(d, [ch14.edge_speed(*SEC, float(q), Gamma_cw=frac * F0.Gamma_kutta) for q in d], col, label=f"Γ = {frac:g} Γ_K")
    ax[1, 0].axhline(ch14.trailing_edge_speed(F0.U, F0.a, F0.b, F0.alpha, F0.beta), color="k", ls=":", label="U (b/a) cos(α + β)")
    ax[1, 0].set(xlabel="distance behind the edge [m]", ylabel="speed [m/s]", title="trailing-edge speed: finite only with Γ_K")
    ax[1, 0].legend(fontsize=8)
    al = np.deg2rad(np.linspace(-10, 22, 161))
    ax[1, 1].plot(np.rad2deg(al), ch14.zhukhovsky_lift_coefficient(al, I["beta"]), "k--", label="2π(α + β)  (14.12)")
    ax[1, 1].plot(np.rad2deg(al), [ch14.zhukhovsky_lift_coefficient(float(a), I["beta"], I["e"], exact=True) for a in al], "C0", label="exact 8π(a/c) sin(α + β)")
    ax[1, 1].plot(np.rad2deg(al), ch14.generic_lift_curve(al, K=I["K"], beta=I["beta"], alpha_stall=I["alpha_stall"], CL_post=I["CL_post"], width=I["stall_width"]),
                  "C3", label="generic measured curve (model, ours)")
    ax[1, 1].axhline(0, color="0.6", lw=0.6)
    ax[1, 1].axvline(-np.rad2deg(I["beta"]), color="0.6", lw=0.6)
    ax[1, 1].set(xlabel="α [deg]", ylabel="C_L", title="lift curve: zero at α = −β, slope 2π")
    ax[1, 1].legend(fontsize=8)
    r = ch14.impulsive_start_lumped_vortex(I["U_section"], 1.0, I["alpha"], 20 / I["U_section"], 0.05 / I["U_section"])
    ax[1, 2].plot(r["s_chords"], r["ratio"], "C0", label="bound Γ / Γ_Kutta (lumped vortex)")
    ax[1, 2].plot(r["s_chords"], r["Gamma_wake_cw"] / r["Gamma_steady"], "C3", label="shed Γ / Γ_Kutta")
    ax[1, 2].plot(r["s_chords"], r["total"] / r["Gamma_steady"], "C1", lw=3, label="total (Kelvin: zero)")
    ax[1, 2].plot(r["s_chords"], ch14.wagner_function(2 * r["s_chords"]), "k--", label="Wagner: circulatory LIFT / steady lift (Jones' fit) — not a circulation")
    ax[1, 2].set(xlabel="chords travelled", ylabel="circulation / Γ_Kutta  (dashed: lift ratio)", title="starting vortex")
    ax[1, 2].legend(fontsize=8)
    save(fig, "fig_kutta")


def fig_pressure():
    """Fig. 14.8: C_p round the section against chord (negative up); the loop area is the normal force."""
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for al, col in ((0.0, "C0"), (5.0, "C2"), (10.0, "C3")):
        S = ch14.zhukhovsky_surface(I["U_section"], I["b"], I["e"], I["beta"], np.deg2rad(al), 1600)
        Cy = ch14.normal_force_coefficient(S["upper"]["xc"], S["upper"]["cp"], S["lower"]["xc"], S["lower"]["cp"])
        ax[0].plot(S["upper"]["xc"], S["upper"]["cp"], col, label=f"α = {al:g}°: loop area {Cy:.3f}")
        ax[0].plot(S["lower"]["xc"], S["lower"]["cp"], col, ls="--")
    ax[0].invert_yaxis()
    ax[0].axhline(1.0, color="0.6", lw=0.6)
    ax[0].set(xlabel="x / c", ylabel="C_p (negative up)", title="surface pressure (solid upper, dashed lower)")
    ax[0].legend(fontsize=8)
    S = ch14.zhukhovsky_surface(*SEC, 1600)
    ax[1].plot(S["s"], S["u_s"] / I["U_section"], "C0")
    ax[1].axhline(0, color="0.6", lw=0.6)
    ax[1].set(xlabel="arc length from the trailing edge, counterclockwise [m]", ylabel="u_s / U", title="signed surface speed (zero: front stagnation point)")
    al = np.deg2rad(np.linspace(0, 10, 21))
    xs = [ch14.suction_side_separation(float(a), I["e"], I["beta"], 1e6)["x_sep_over_c"] for a in al]
    xp = [ch14.suction_side_separation(float(a), I["e"], I["beta"], 1e6)["xc_peak"] for a in al]
    ax[2].plot(np.rad2deg(al), xs, "C3o-", label="laminar separation (Thwaites) — qualitative")
    ax[2].plot(np.rad2deg(al), xp, "C0s-", label="suction peak")
    ax[2].set(xlabel="α [deg]", ylabel="x / c", title="separation moves forward with α (trend only)")
    ax[2].legend(fontsize=8)
    save(fig, "fig_pressure")


def fig_downwash():
    """Figs. 14.21–14.22: downwash of a horseshoe and of a staircase; downwash ahead of and behind the wing."""
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    s, G1 = 6.0, 1.0
    y = np.linspace(-1.6 * s / 2, 1.6 * s / 2, 801)
    w = ch14.horseshoe_downwash(0.0, y, G1, s)
    ax[0].plot(y, np.clip(w, -0.6, 0.6), "C0", label="one horseshoe (clipped)")
    ax[0].axhline(G1 / (PI * s), color="k", ls=":", label="Γ/(πs) at the centre")
    ax[0].axhline(0, color="0.6", lw=0.6)
    for yy in (-s / 2, s / 2):
        ax[0].axvline(yy, color="0.6", lw=0.6)
    ax[0].set(xlabel="y [m]", ylabel="downwash w [m/s]", title="horseshoe: downwash inboard, upwash outboard")
    ax[0].legend(fontsize=8)
    fn = lambda yy: ch14.elliptic_circulation(yy, s, G1)
    for N, col in ((8, "C1"), (32, "C2"), (128, "C3")):
        ye, Gs = ch14.staircase_from_loading(fn, s, N)
        ym = 0.5 * (ye[:-1] + ye[1:])
        ax[1].plot(ym, ch14.staircase_downwash(ym, ye, Gs), col, marker=".", ms=3, label=f"{N} horseshoes")
    yq = np.linspace(-0.49 * s, 0.49 * s, 25)
    ax[1].plot(yq, ch14.downwash(yq, fn, s), "ko", ms=4, label="quadrature of (14.13)")
    ax[1].axhline(G1 / (2 * s), color="k", ls=":", label="Γ₁/2s  (14.24)")
    ax[1].set(xlabel="y [m]", ylabel="w [m/s]", ylim=(0, 0.2), title="elliptic loading: uniform downwash")
    ax[1].legend(fontsize=8)
    x = np.linspace(-3 * s, 3 * s, 241)
    for Gn, col, name in (([1.0], "C0", "elliptic"), ([1.0, 0.0, 0.15], "C3", "Γ₃ = 0.15 Γ₁")):
        w0 = ch14.downwash_series(Gn, PI / 2, s)
        ax[2].plot(x / s, ch14.downwash_along_x(x, Gn, s) / w0, col, label=name)
    ax[2].axhline(1, color="0.6", lw=0.6)
    ax[2].axhline(2, color="0.6", lw=0.6)
    ax[2].set(xlabel="x / s (downstream +)", ylabel="w(x) / w(0)", title="0 far ahead, w at the line, 2w far behind")
    ax[2].legend(fontsize=8)
    save(fig, "fig_downwash")


def fig_lifting_line():
    """Fig. 14.24 and (14.20)–(14.26): loadings of three planforms, convergence, C_Di = C_L²/(πΛ)."""
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    s, A, K, al = 6.0, 6.0, 2 * PI, np.deg2rad(5.0)
    plans = {"rectangular": ch14.planform("rectangular", s, A / s), "tapered 0.45": ch14.planform("tapered", s, 2 * A / (s * 1.45), 0.45),
             "elliptic": ch14.planform("elliptic", s, 4 * A / (PI * s))}
    gg = np.linspace(0, PI, 201)
    for (name, P), col in zip(plans.items(), ("C0", "C1", "C2")):
        r = ch14.solve_lifting_line(P, s, 1.0, al, K=K, N=40, symmetric=True)
        ax[0].plot(ch14.span_position(gg, s), ch14.circulation_series(r["Gn"], gg), col, label=f"{name}: e_span = {r['e_span']:.4f}, C_L = {r['CL']:.4f}")
    ax[0].set(xlabel="y [m]", ylabel="Γ(y) [m²/s]", title="solutions of (14.20), Λ = 6, α = 5°")
    ax[0].legend(fontsize=8)
    Ns = [4, 8, 16, 32, 64, 128]
    CL = [ch14.solve_lifting_line(1.0, 6.0, 1.0, al, K=K, N=N)["CL"] for N in Ns]
    ax[1].loglog(Ns[:-1], np.abs(np.diff(CL)), "C0o-", label="|C_L(2N) − C_L(N)|, rectangular")
    ax[1].loglog(Ns[:-1], 0.5 * np.array(Ns[:-1], float) ** -4, "k--", label="N⁻⁴")
    ax[1].set(xlabel="N", title="convergence of the collocation solution")
    ax[1].legend(fontsize=8)
    cl = np.linspace(0, 1.4, 50)
    for AR, col in ((4, "C3"), (9.6, "C0"), (20, "C2")):
        ax[2].plot(ch14.drag_polar(cl, CD0=I["CD0"], AR=AR, e_span=I["e_span"]), cl, col, label=f"Λ = {AR:g}")
        bg = ch14.best_glide(I["CD0"], AR, I["e_span"])
        ax[2].plot([0, bg["CD"]], [0, bg["CL"]], col, ls=":", lw=1)
        ax[2].plot([bg["CD"]], [bg["CL"]], col, marker="o")
    ax[2].set(xlabel="C_D", ylabel="C_L", xlim=(0, 0.2), title="polar C_D = C_D0 + C_L²/(πΛ e_span); tangent = best glide")
    ax[2].legend(fontsize=8)
    save(fig, "fig_lifting_line")


def fig_drag_and_sail():
    """§14.7 and §14.9: the U-shaped drag curve; the sail's drive against the course angle; the generic drag curve."""
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    U = np.linspace(22, 80, 200)
    d = ch14.drag_vs_speed(I["W"], I["rho"], I["A"], I["AR"], I["CD0"], U, I["e_span"])
    ax[0].plot(U, d["profile"], "C0", label="profile ∝ U²")
    ax[0].plot(U, d["induced"], "C3", label="induced ∝ 1/U²")
    ax[0].plot(U, d["total"], "k", label="total")
    Us = ch14.speed_min_drag(I["W"], I["rho"], I["A"], I["AR"], I["CD0"], I["e_span"])
    ax[0].axvline(Us, color="0.5", ls=":")
    ax[0].axvline(ch14.stall_speed(I["W"], I["rho"], I["A"], I["CL_max"]), color="C1", ls="--", label="stall speed")
    ax[0].set(xlabel="U [m/s]", ylabel="drag [N]", title=f"level flight: least drag at {Us:.1f} m/s (parts equal)")
    ax[0].legend(fontsize=8)
    th = np.deg2rad(np.linspace(5, 175, 171))
    kw = dict(CL=I["CL_sail"], CD=I["CD_sail"], rho=I["rho"], A=I["A_sail"])
    ax[1].plot(np.rad2deg(th), ch14.sail_drive_curve(th, U=I["U_wind"], V=I["V_boat"], **kw), "C4", label="drive")
    ax[1].plot(np.rad2deg(th), [ch14.sail_forces(I["U_wind"], I["V_boat"], float(t), **kw)["lateral"] for t in th], "C1", label="lateral")
    ax[1].axhline(0, color="0.6", lw=0.6)
    ax[1].set(xlabel="course angle off the true wind θ [deg]", ylabel="force [N]", title="sail: the drive vanishes close to the wind")
    ax[1].legend(fontsize=8)
    al = np.deg2rad(np.linspace(-10, 25, 141))
    kk = dict(K=I["K"], beta=I["beta"], alpha_stall=I["alpha_stall"], CL_post=I["CL_post"], width=I["stall_width"])
    ax[2].plot(np.rad2deg(al), ch14.generic_lift_curve(al, **kk), "C0", label="C_L (model)")
    ax[2].plot(np.rad2deg(al), 5 * ch14.generic_drag_curve(al, CD0=I["CD0"], AR=I["AR"], e_span=I["e_span"], CD_sep=I["CD_sep"], **kk), "C3", label="5 C_D (model)")
    ax[2].set(xlabel="α [deg]", title="generic measured curves — qualitative, constants ours")
    ax[2].legend(fontsize=8)
    save(fig, "fig_drag_and_sail")


def fig_wake():
    """Figs. 14.19–14.20: roll-up of the trailing sheet (qualitative) and its conserved quantities."""
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    r = ch14.wake_rollup([1.0], 6.0, n=160, n_t=41)
    for k, col in ((0, "0.6"), (10, "C0"), (20, "C2"), (40, "C3")):
        ax[0].plot(r["y"][k], r["z"][k], ".", color=col, ms=3, label=f"t = {r['t'][k]:.1f} s")
    ax[0].set(aspect="equal", xlabel="y [m] (seen from behind)", ylabel="z [m]", title="roll-up of the elliptic sheet (qualitative)")
    ax[0].legend(fontsize=8)
    ax[1].plot(r["t"], r["impulse"] / r["impulse"][0], "C0", label="impulse Σσy")
    ax[1].plot(r["t"], r["energy"] / r["energy"][0], "C3--", label="blob Hamiltonian")
    ax[1].plot(r["t"], 2 * r["centroid_y"] / r["spacing_theory"], "C2:", label="centroid spacing / (πs/4)")
    ax[1].set(xlabel="t [s]", ylim=(0.99, 1.01), title="invariants (ratio to the start / to theory)")
    ax[1].legend(fontsize=8)
    ax[2].plot(r["t"], -r["centroid_z"], "C0", label="descent of the half-sheet centroid")
    ax[2].plot(r["t"], r["descent_speed_theory"] * r["t"], "k--", label="rolled-up pair: Γ₀ t / (2π b′)")
    ax[2].plot(r["t"], (1.0 / 6.0) * r["t"], "k:", label="flat sheet: 2w t = Γ₁ t / s")
    ax[2].set(xlabel="t [s]", ylabel="descent [m]", title="sinking of the wake")
    ax[2].legend(fontsize=8)
    save(fig, "fig_wake")


if __name__ == "__main__":
    fig_shapes()
    fig_kutta()
    fig_pressure()
    fig_downwash()
    fig_lifting_line()
    fig_drag_and_sail()
    fig_wake()
