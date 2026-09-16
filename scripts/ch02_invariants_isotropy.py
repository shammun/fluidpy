"""Chapter 2, §2.5 and §2.7: the three invariants I₁, I₂, I₃ stay flat over 50 random rotations while the tensor entries
scatter (Exercise 2.9), the ε_ijk cube unfolded, the epsilon–delta check and the isotropy residuals of δ, ε and a
random tensor (Exercise 2.11; ε flips sign under a reflection).

Run: ``.venv/Scripts/python.exe scripts/ch02_invariants_isotropy.py --no-show``
Figure → outputs/ch02/invariants_isotropy.png.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "outputs" / "ch02"))
    ap.add_argument("--no-show", action="store_true")
    args = ap.parse_args()
    t0 = time.perf_counter()
    import matplotlib

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS, use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    A = rng.standard_normal((3, 3))
    A = A + A.T  # symmetric, so I₁, I₂, I₃ also match Σλ, Σλλ, Πλ
    inv0 = ch02.invariants(A)
    entries, invs = [], []
    for _ in range(50):
        Ap = ch02.transform_tensor(A, ch02.random_rotation(rng))  # Eq. (2.12)
        entries.append(Ap.ravel())
        invs.append(ch02.invariants(Ap))
    entries, invs = np.array(entries), np.array(invs)
    lam = np.linalg.eigvalsh(A)
    vieta = (lam.sum(), lam[0] * lam[1] + lam[0] * lam[2] + lam[1] * lam[2], lam.prod())
    eps = ch02.levi_civita()
    iso = {"δ_ij": ch02.is_isotropic(np.eye(3), rng)[1], "ε_ijk (rotations)": ch02.is_isotropic(eps, rng)[1],
           "ε_ijk (+reflections)": ch02.is_isotropic(eps, rng, proper=False)[1], "random A_ij": ch02.is_isotropic(A, rng)[1]}

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    ax = axs[0]
    for k in range(9):
        ax.plot(entries[:, k], color=COLORS["grid"], lw=1)
    for k, (name, col) in enumerate(zip(("I₁ = A_ii", "I₂ = ½(I₁² − A_ij A_ji)", "I₃ = det A"), (COLORS["accent"], COLORS["teal"], COLORS["orange"]))):
        ax.plot(invs[:, k], color=col, lw=2.2, label=f"{name} = {inv0[k]:.3f}")
    ax.set_xlabel("random rotation #")
    ax.set_ylabel("value")
    ax.set_title("entries scatter (grey), invariants stay flat (Exercise 2.9)")
    ax.legend(fontsize=8)
    ax = axs[1]
    unfolded = np.hstack([eps[i] for i in range(3)])  # three 3×3 slices i = 1, 2, 3 side by side
    ax.imshow(unfolded, cmap="coolwarm", vmin=-1, vmax=1)
    for r in range(3):
        for c in range(9):
            ax.text(c, r, f"{int(unfolded[r, c]):+d}" if unfolded[r, c] else "0", ha="center", va="center", fontsize=9)
    ax.set_xticks(range(9), [f"{i + 1}{j + 1}" for i in range(3) for j in range(3)], fontsize=7)
    ax.set_yticks(range(3), ["k=1", "k=2", "k=3"])
    ax.set_xlabel("slice i (left→right: i = 1, 2, 3), column j")
    ax.set_title(f"ε_ijk unfolded (Eq. 2.18); ε–δ residual = {ch02.epsilon_delta_residual():.0f}")
    ax = axs[2]
    names = list(iso)
    vals = [max(v, 1e-17) for v in iso.values()]
    ax.bar(names, vals, color=[COLORS["teal"], COLORS["teal"], COLORS["rose"], COLORS["rose"]])
    ax.set_yscale("log")
    ax.set_ylabel("max |T′ − T| over 50 frames")
    ax.set_title("isotropy test (Exercise 2.11): δ and ε pass; ε flips under a reflection")
    ax.tick_params(axis="x", labelsize=8)
    fig.savefig(out / "invariants_isotropy.png", bbox_inches="tight")

    print(f"I₁, I₂, I₃ = {np.round(inv0, 6)}; spread over 50 rotations = {np.round(invs.max(0) - invs.min(0), 12)}")
    print(f"Vieta from eigenvalues {np.round(lam, 4)}: Σλ = {vieta[0]:.6f}, Σλλ = {vieta[1]:.6f}, Πλ = {vieta[2]:.6f}")
    print(f"characteristic polynomial coefficients {np.round(ch02.characteristic_polynomial(A), 6)}; roots {np.round(np.sort(np.roots(ch02.characteristic_polynomial(A))), 6)}")
    print(f"ε contractions: ε_pqi ε_pqj = 2δ → {np.einsum('pqi,pqj->ij', eps, eps).tolist()}; ε_pqr ε_pqr = {np.einsum('pqr,pqr', eps, eps):.0f}; δ_ii = {np.trace(np.eye(3)):.0f}")
    print(f"isotropy residuals: {{{', '.join(f'{k}: {v:.1e}' for k, v in iso.items())}}}")
    print(f"saved invariants_isotropy.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
