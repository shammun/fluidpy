"""Chapter 2, §2.3, Eqs. (2.9)–(2.11): the matrix product P_ij = A_ik B_kj with row i of A and column j of B lit up,
assembled term by term (the "matrix view" pattern), and the summation convention printing all nine sums.

Run: ``.venv/Scripts/python.exe scripts/ch02_matrix_product_highlight.py --no-show``
Figure → outputs/ch02/eq2_11_matrix_product.png.
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


def _draw_matrix(ax, M, x0, name, hl_row=None, hl_col=None, color="#f97316"):
    n = M.shape[0]
    for i in range(n):
        for j in range(n):
            face = color if (i == hl_row or j == hl_col) else "white"
            ax.add_patch(plt_rect((x0 + j, n - 1 - i), face))
            ax.text(x0 + j + 0.5, n - 1 - i + 0.5, f"{M[i, j]:g}", ha="center", va="center", fontsize=12)
    ax.text(x0 + n / 2, n + 0.25, name, ha="center", fontsize=13)


def plt_rect(xy, face):
    from matplotlib.patches import Rectangle

    return Rectangle(xy, 1, 1, facecolor=face, edgecolor="#1a1f36", lw=1, alpha=0.9 if face != "white" else 1.0)


DEFAULT_A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
DEFAULT_B = np.array([[9, 8, 7], [6, 5, 4], [3, 2, 1]], dtype=float)


def matrix_product_figure(i: int = 1, j: int = 2, A=None, B=None, figsize=(11, 4.2)):
    """Eq. (2.11) as a picture: P = A·B with row i of A (teal) and column j of B (orange) lit up, the three products
    and their sum written under the highlighted cell P_ij.

    Book: §2.3, Eqs. (2.9)–(2.11): P_ij = A_ik B_kj — "the ij-element of P is determined by multiplying the elements in
    the i-row of A and the j-column of B, and summing" (the dashed boxes of (2.11)). Drawing only — P comes from
    ``ch02.inner(A, B)``.

    Parameters
    ----------
    i, j : **book (1-based) indices** of the element P_ij to highlight (the book's example is P_12, the storyboard's
        "change" is P_31 = 7·1 + 8·4 + 9·7 = 102 for B = A)
    A, B : (3, 3) matrices (default the running 1…9 and 9…1 matrices)
    figsize : matplotlib figure size [in]

    Returns
    -------
    matplotlib Figure.
    """
    import matplotlib.pyplot as plt

    from fluidpy import ch02_cartesian_tensors as ch02
    from fluidpy.core.style import COLORS

    A = DEFAULT_A if A is None else np.asarray(A, dtype=float)
    B = DEFAULT_B if B is None else np.asarray(B, dtype=float)
    if A.shape != (3, 3) or B.shape != (3, 3):
        raise ValueError("matrix_product_figure draws the book's 3 × 3 case")
    if not (1 <= i <= 3 and 1 <= j <= 3):
        raise ValueError("i and j are book indices 1, 2, 3")
    P = ch02.inner(A, B)  # Eq. (2.9): P_ij = A_ik B_kj
    ii, jj = i - 1, j - 1  # 0-based for the arrays
    terms = [A[ii, k] * B[k, jj] for k in range(3)]  # the three products of the hidden sum over k

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-0.5, 12.5)
    ax.set_ylim(-1.6, 4.2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.grid(False)
    _draw_matrix(ax, P, 0, "P", hl_row=None, hl_col=None)
    ax.add_patch(plt_rect((0 + jj, 3 - 1 - ii), COLORS["accent"]))
    ax.text(jj + 0.5, 3 - 1 - ii + 0.5, f"{P[ii, jj]:g}", ha="center", va="center", fontsize=12, color="white")
    ax.text(3.5, 1.5, "=", fontsize=18, ha="center", va="center")
    _draw_matrix(ax, A, 4.5, "A", hl_row=ii, color=COLORS["teal"])
    _draw_matrix(ax, B, 8.5, "B", hl_col=jj, color=COLORS["orange"])
    ax.text(6, -0.9, f"$P_{{{i}{j}}} = A_{{{i}k}}B_{{k{j}}}$ = " + " + ".join(f"{A[ii, k]:g}·{B[k, jj]:g}" for k in range(3))
            + f" = {' + '.join(f'{t:g}' for t in terms)} = {P[ii, jj]:g}", ha="center", fontsize=12)
    ax.set_title(f"Eq. (2.11): the {i}{j}-element of P multiplies row {i} of A (teal) by column {j} of B (orange) and sums over k")
    return fig


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
    from fluidpy.core.style import use_style

    use_style()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    A, B = DEFAULT_A, DEFAULT_B
    P = ch02.inner(A, B)  # Eq. (2.9): P_ij = A_ik B_kj
    assert np.allclose(P, A @ B)
    i, j = 0, 1  # the book's example P_12 = A_11 B_12 + A_12 B_22 + A_13 B_32
    terms = [A[i, k] * B[k, j] for k in range(3)]

    fig = matrix_product_figure(1, 2, A, B)  # book indices: P_12
    fig.savefig(out / "eq2_11_matrix_product.png", bbox_inches="tight")
    fig_aa = matrix_product_figure(i=3, j=1, A=A, B=A)  # the storyboard's "change": P_31 of A·A = 102
    fig_aa.savefig(out / "eq2_11_matrix_product_A_A_31.png", bbox_inches="tight")
    print(f"matrix_product_figure(3, 1, A, A): P_31 = {ch02.inner(A, A)[2, 0]:g} = 7·1 + 8·4 + 9·7")

    print("expand_indices('A_ik B_kj'), the nine sums the convention hides:")
    print(ch02.expand_indices_str("A_ik B_kj", lhs="P"))
    print(f"P_12 numerically: {P[i, j]:g} = {' + '.join(f'{t:g}' for t in terms)}")
    print(f"classify_indices('A_ik B_kj') = {ch02.classify_indices('A_ik B_kj')}; tensor_order = {ch02.tensor_order('A_ik B_kj')}")
    print(f"the four contractions of (2.14): "
          f"{np.allclose(ch02.contract(A, B, 'ij,ki->kj'), B @ A)}, {np.allclose(ch02.contract(A, B, 'ij,ik->jk'), A.T @ B)}, "
          f"{np.allclose(ch02.contract(A, B, 'ij,kj->ik'), A @ B.T)}, {np.allclose(ch02.contract(A, B, 'ij,jk->ik'), A @ B)}")
    print(f"double dot: book A:B = A_ij B_ji = {ch02.double_dot(A, B, 'book'):g} = trace(A·B) = {np.trace(A @ B):g}; "
          f"Frobenius A_ij B_ij = {ch02.double_dot(A, B, 'frobenius'):g} = trace(A·Bᵀ) = {np.trace(A @ B.T):g}")
    print(f"saved eq2_11_matrix_product.png in {out}  ({time.perf_counter() - t0:.1f} s)")
    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
