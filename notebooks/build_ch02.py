"""Build the Chapter 2 teaching notebook: ``notebooks/ch02_cartesian_tensors.ipynb``.

Source of truth: ``analysis/ch02_design.md`` Part A (storyboard, one nbkit call per row), Part E (prerequisite ledger →
25 primers, ch01 reminders, glosses), Part F (the 15 derivations, copied move by move) and ``analysis/ch02_curation.md``
(IDs, depths, section coverage). Physics lives in ``fluidpy.ch02_cartesian_tensors`` (imported as ``ch02``); cells only
call it. Figures that are pure drawing come from ``scripts/ch02_*.py``.

Labels shown to the reader: CORE blocks carry their id (`C05`), notes their id (`N38`), primers their number (`P62` …,
continuing ch01's P01–P61), so cross-references such as "(P26, C05)" can be followed.

Conventions (design header, binding): passive rotation C_ij = e_i·e'_j (columns = new axes), x' = Cᵀx, τ' = CᵀτC;
traction contracts the first index f_i = τ_ji n_j; tensor divergence contracts the second; G[i, j] = ∂u_i/∂x_j;
Ex. 2.4's Γ ≡ S₁₂; Stokes: n_c points into A and t = n_c × n runs counterclockwise about n.

Run:  .venv/Scripts/python.exe notebooks/build_ch02.py
"""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch02")


# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent (same shape as notebooks/build_ch01.py)
# ---------------------------------------------------------------------------------------------------------------------
def core(cid: str, title: str, question: str) -> None:
    """A CORE block heading with its id."""
    nb.core(cid, f"{title} `{cid}`", question=question)


def P(pid: str, term: str, text: str, code: str | None = None) -> None:
    """A 📎 primer; ``term`` is exactly the Part E concept text (the ledger check matches it)."""
    nb.primer(term, f"`{pid}` · {text}", code=code)


def note(nid: str, title: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    """A B/C note with its curation id."""
    nb.note(f"**{title}** `{nid}` — {text}", equation=equation, ref=ref)


def D(key: str, title: str, *, ref: str, goal: str, assumptions: str, start, plan, uses, steps, result,
      interpret: str, check: str, traps: str, check_src: str | None = None) -> None:
    """A Part F derivation: goal + assumptions, plan, tools, one move per step, result, meaning, check, traps."""
    nb.derivation(key, title, ref=ref, goal=f"{goal}\n\n**Assumptions.** {assumptions}", start=start, plan=plan,
                  uses=uses, steps=[dict(did=a, tex=b, why=c, plain=d) for a, b, c, d in steps], result=result,
                  interpret=interpret, check=check, check_src=check_src)
    nb.md(f"> ⚠️ **Common confusion (traps in this derivation):** {traps}")


def see_read_change(see: str, read: str, change: str) -> None:
    nb.figure_notes(see, read, change)


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=(
        "A vector or a tensor is **one physical thing**; its components are what a particular set of axes sees. This "
        "chapter builds the index machinery (a repeated index is a sum), the rotation rule that *defines* vectors and "
        "tensors, the two tensors fluid mechanics runs on (the stress and the velocity gradient), and the two theorems "
        "that turn 'inside a volume' into 'on its boundary' (Gauss and Stokes). Everything here is mathematics, but it "
        "is the mathematics Ch. 3–5 and Ch. 13 are written in."
    ),
    roadmap=[
        "§2.1–2.3 the summation convention (C01), the direction-cosine matrix (C02) and how components rotate (C03)",
        "§2.4–2.6 the stress tensor and its signs (C04), Cauchy's traction f = n·τ (C05), the rule that makes τ a tensor (C06)",
        "§2.5, 2.7, 2.8 contraction and invariants (C07); δ, ε and the cross product (C08)",
        "§2.9 gradient, divergence, curl — on a grid you will reuse in every later chapter (C09–C11)",
        "§2.10–2.11 symmetric + antisymmetric parts (C12); principal axes (C13)",
        "§2.12–2.14 Gauss (C14, C15) and Stokes (C16); comma notation",
    ],
    prerequisites=[
        "vectors and matrices (Ch. 1 primers P53–P58 for determinants, rank, null space)",
        "partial derivatives (Ch. 1 P25) and the first-order Taylor expansion (P26)",
        "the definite integral (P27) and the trapezoid rule (P37)",
        "stress as force per area (Ch. 1 §1.3, P05) and the free-body diagram (P09)",
    ],
)
nb.explainer_index([
    ("rotation_of_axes", "Which rotates — the arrow or the ruler?", "C02 C03 C06: passive rotation, x' = Cᵀx, τ' = CᵀτC"),
    ("cauchy_traction_principal_axes", "Cut the point any way you like: what pushes on the cut?", "C05 C13 C04: f = n·τ, Mohr's circle, principal axes"),
    ("strain_vs_rotation_split", "Is simple shear a rotation?", "C12: G = S + A, the hidden vector of A"),
    ("gauss_flux_box", "What leaks out of a box?", "C14 C15: Gauss' theorem → divergence as outflux per volume"),
    ("stokes_circulation_loop", "How much does the flow go round a loop?", "C16 C11: Stokes' theorem → curl as circulation per area"),
])
nb.setup()
nb.code("""
from fluidpy import ch02_cartesian_tensors as ch02      # the tested chapter-2 module: every function cites its book section and equation
import numpy as np                                      # arrays (Ch. 1 primer P03)
import sympy as sp                                      # symbolic algebra (Ch. 1 primer P40)
import matplotlib.pyplot as plt                         # static figures (Ch. 1 primer P01)
import plotly.graph_objects as go                       # rotatable 3-D figures (Ch. 1 primer P41)
from fluidpy.core.style import COLORS                   # the house palette: teal = old axes, orange = new axes, blue = normal, rose = shear
import logging                                          # standard library: controls library log messages
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless font-substitution notes (as in Ch. 1)
print(f"{len(ch02.__all__)} public functions in fluidpy.ch02_cartesian_tensors")   # the toolbox this notebook calls
""", explain="""
1. `ch02` is the chapter module; it re-exports the reusable primitives of `fluidpy/core/` (tensors, index notation,
   grids, differential operators, integral theorems) so that one name covers the whole chapter.
2. The other imports are the Ch. 1 toolbox: numpy (P03), sympy (P40), matplotlib (P01) and plotly (P41), plus the
   shared colour palette.
""")
nb.md(r"""
**Notation used in this notebook.** Indices $i, j, k, l, m, n$ run over 1, 2, 3 (Python: 0, 1, 2 — we say which
whenever a printed index appears). $\equiv$ means "is defined as". Boldface $\mathbf u$ is the object; $u_i$ its
components. A superscript on an eigenvalue, $\lambda^k$, is a *label*, not a power. Colours: **teal** = the original
axes, **orange** = rotated axes, **blue** = normal stress, **rose** = shear, **purple** = the quantity being taught.
""")
nb.md(r"""
**Where this chapter is used later**

| Result here | Used in |
|---|---|
| the tensor rule (2.12) | Ch. 3 strain rate, Ch. 4 stress, Ch. 12 Reynolds stress |
| Cauchy's formula (2.15) | Ch. 4 §4.3 wall forces, Ch. 8–9 drag and skin friction |
| $\varepsilon_{ijk}$ (2.18) | Ch. 4 Coriolis $2\boldsymbol\Omega\times\mathbf u$, Ch. 5 vorticity identities |
| Gauss (2.30) | Ch. 4 §4.2 continuity, every conservation law in differential form |
| Stokes (2.34) | Ch. 5 Kelvin's theorem, Ch. 6 and 14 lift |

*Climate hook:* the Coriolis term and the geostrophic balance of Ch. 13 are one $\varepsilon_{ijk}$ and one $\nabla p$
away from what is built here.
""")

# =====================================================================================================================
# A.1 §2.1 Scalars, Vectors, Tensors, Notation — C01
# =====================================================================================================================
nb.section("2.1", "Scalars, Vectors, Tensors, Notation", intro="""
**What is this section about?** Three kinds of quantity (one number, three numbers, nine numbers), the index notation
that writes all of them as symbols with subscripts, and the one convention — a repeated index means "add over 1, 2, 3"
— that lets the whole book drop its Σ signs.
""")
note("N01", "A scalar (zero-order tensor)", r"""
is one number at each point, the same whatever axes you choose: pressure $p$, temperature $T$, density $\rho$ — the
fields of Ch. 1. Zero free indices.
""")
note("N02", "A vector (first-order tensor)", r"""
has a size and a direction: position $\mathbf x$, velocity $\mathbf u$, gravity $\mathbf g$. Its three components
$u_1, u_2, u_3$ depend on the axes; the arrow does not. §2.2 sharpens this into a test (C03).
""")
note("N06", "A second-order tensor", r"""
has one component for each *pair* of directions, $3\times3 = 9$: the stress $\tau_{ij}$ (which face, which force
direction — C04). Two free indices.
""")
note("N24", "Hierarchy 1 → 3 → 9", r"""
Newtonian fluid mechanics needs orders 0–2, with one fourth-order exception: the viscosity tensor of Ch. 4 §4.5 (N28).
""")
core("C01", "The summation convention: a repeated index is a sum", r"""
Why does the book never write a Σ sign — and how do you read $a_i b_i$, $A_{ik}B_{kj}$ or $C_{im}C_{jn}\tau_{ij}$
without one?
""")
nb.md(r"""
#### The problem in plain words

Every equation of fluid mechanics is really three or nine equations at once (one per direction, one per pair). Written
out, the momentum equation fills a page; in index notation it is one line. The price is a reading rule you must run in
your head: **if a letter appears twice in a term, add that term up for the letter = 1, 2, 3.** Once the rule is
automatic, $\tau_{ji}n_j$ *is* a matrix–vector product and $\partial u_i/\partial x_i$ *is* the divergence.
""")
nb.md(r"""
#### The idea

```
a_i b_i        →  i is repeated     →  a_1 b_1 + a_2 b_2 + a_3 b_3            (one number:   no free index)
A_ik B_kj      →  k is repeated     →  A_i1 B_1j + A_i2 B_2j + A_i3 B_3j      (nine numbers: free i, j)
C_im C_jn τ_ij →  i, j repeated     →  9 terms for each (m, n)                (nine numbers: free m, n)
```

**Free index = shape of the answer; repeated (dummy) index = a loop.** A dummy letter can be renamed at will
($a_ib_i = a_kb_k$); a free letter must match on both sides of an equation.
""")
note("N03", "The position vector in a basis", r"""
$\mathbf x = \mathbf e_1 x_1 + \mathbf e_2 x_2 + \mathbf e_3 x_3$ — the position vector is its components times the
unit vectors; the subscript on $\mathbf e$ names an *axis*, not a component. In code `ch02.unit_vectors()` returns the
three $\mathbf e_i$ as rows and `ch02.vector_from_components(x, E)` rebuilds the arrow. **Column form** `N04`: the book
also writes $\mathbf x = [x_1\ x_2\ x_3]^{\rm T}$ (the transpose ᵀ turns a row into a column) or a triplet
$(x_1, x_2, x_3)$ — in code `x = np.array([1., 2., 3.])`, `x.T`.
""", equation=r"\mathbf x = \mathbf e_1 x_1 + \mathbf e_2 x_2 + \mathbf e_3 x_3", ref="2.1")
P("P62", "np.einsum index strings", """
The code form of the summation convention: `np.einsum('i,i', a, b)` sums over the repeated letter i; `'ik,kj->ij'`
sums over k and keeps i, j; the letters after `->` are the free indices (the shape of the answer). Nothing else is
needed to write any formula of this chapter.
""", code="""
a, b = np.array([1., 2., 3.]), np.array([4., 5., 6.])      # two vectors (any unit)
print(np.einsum('i,i', a, b))                              # 32.0 = 4 + 10 + 18: the repeated i is summed
A = np.arange(9.).reshape(3, 3)                            # a 3×3 matrix 0…8
print(np.einsum('ik,kj->ij', A, A)[0, 1])                  # P_12 = A_1k A_k2 (Python row 0, column 1) = 0·1 + 1·4 + 2·7 = 18
""")
P("P63", "matrix multiplication, transpose and identity", r"""
Row-by-column: $(AB)_{ij} = \sum_k A_{ik}B_{kj}$ — walk along row i of A and down column j of B, multiply pairs, add.
In Python `A @ B`; `A.T` swaps rows and columns ($A^{\rm T}_{ij} = A_{ji}$); `np.eye(3)` is the identity (1 on the
diagonal, 0 elsewhere) — it leaves any matrix unchanged.
""", code="""
A = np.array([[1., 2.], [3., 4.]]); B = np.array([[0., 1.], [1., 0.]])   # two 2×2 matrices
print(A @ B)             # [[2, 1], [4, 3]]: row i of A · column j of B
print(A.T)               # [[1, 3], [2, 4]]: rows and columns swapped
print(np.eye(2) @ A)     # A again: the identity changes nothing
""")
nb.md(r"""
#### The maths, step by step (three equations expanded by hand)

1. **Dot product (2.2):** $\mathbf a\cdot\mathbf b = a_1b_1 + a_2b_2 + a_3b_3 = \sum_{i=1}^{3} a_ib_i \equiv a_ib_i$ —
   one repeated letter, one sum, a scalar.
2. **Matrix product (2.9):** $P_{ij} = \sum_{k=1}^{3} A_{ik}B_{kj} \equiv A_{ik}B_{kj}$ — k summed, i and j free, nine
   numbers.
3. **Kronecker substitution (2.17):** $\delta_{ij}u_j = \delta_{i1}u_1 + \delta_{i2}u_2 + \delta_{i3}u_3 = u_i$ — only
   the term with j = i survives, so δ *replaces* its summed index.
""")
note("N07", "The dot product as an implied sum", r"""
the ≡ is the summation convention at work: the last form *means* the middle one. **Also** `N45`: $\mathbf u\cdot\mathbf v
= uv\cos\theta$ (the projection of one vector on the other), and $\mathbf u\cdot\mathbf v$ is the trace of the nine
products $u_iv_j$ — the contraction idea of §2.5 (C07); `ch02.angle_between(u, v)` returns θ.
""", equation=r"\mathbf a\cdot\mathbf b = a_1 b_1 + a_2 b_2 + a_3 b_3 = \sum_{i=1}^{3} a_i b_i \equiv a_i b_i", ref="2.2")
note("N40", "The Kronecker delta", r"""
$\delta_{ij} = 1$ if $i = j$, 0 otherwise — the identity matrix in index clothes (`ch02.kronecker_delta()` is
`np.eye(3)`). We need it now because the rotation matrix's orthogonality (C02) is written with it.
""", equation=r"\delta_{ij} = \begin{cases} 1 & i = j \\ 0 & i \neq j\end{cases}", ref="2.16")
note("N41", "δ substitutes its summed index", r"""
$\delta_{ij}u_j = u_i$, $\delta_{ij}A_{jk} = A_{ik}$, and — the classic slip — $\delta_{ii} = 3$ (a summed pair), not 1.
""", equation=r"\delta_{ij} u_j = u_i", ref="2.17")
note("N21", "Matrix product as an index sum", r"""
$P_{ij} = A_{ik}B_{kj}$ sums the *adjacent* index k; `ch02.inner(A, B)` is `np.einsum('ik,kj->ij')` and equals `A @ B`.
**Single dot** `N22`: (2.10) $\mathbf P = \mathbf A\cdot\mathbf B$ — a single dot means one index summed, the same dot
that later reads $\mathbf n\cdot\boldsymbol\tau$ (C05) and $\nabla\cdot\mathbf u$ (C10). **Explicitly** `N23`: (2.11)
$P_{12} = A_{11}B_{12} + A_{12}B_{22} + A_{13}B_{32}$ — row 1 of A against column 2 of B; the figure below boxes them.
""", equation=r"P_{ij} = \sum_{k=1}^{3} A_{ik} B_{kj} \equiv A_{ik} B_{kj}", ref="2.9")
note("N08", "Boldface for meaning, indices for manipulation", r"""
$\mathbf a\mathbf b$ is ambiguous for tensors ($\mathbf a\mathbf b \neq \mathbf b\mathbf a$) and hides the order;
indices never do: **order = number of free indices**. `ch02.tensor_order('A_ij B_kl')` → 4. **Free vs dummy** `N12`:
$x_iC_{ij} = x_kC_{kj}$ (rename the summed letter freely); the free j must appear on both sides. A letter appearing
three times in one term is an error (`ch02.classify_indices` raises).
""")
nb.worked_example("(1, 2, 3)·(4, 5, 6) and δ_ij u_j for i = 2", r"""
1. $a_ib_i$ with i = 1, 2, 3: $1\times4 + 2\times5 + 3\times6 = 4 + 10 + 18 = 32$.
2. $\delta_{2j}u_j = \delta_{21}u_1 + \delta_{22}u_2 + \delta_{23}u_3 = 0 + u_2 + 0 = u_2$.
3. $P_{12}$ for $A = B = \begin{bmatrix}1&2&3\\4&5&6\\7&8&9\end{bmatrix}$: $1\cdot2 + 2\cdot5 + 3\cdot8 = 2 + 10 + 24 = 36$.
""")
nb.code("""
a, b = np.array([1., 2., 3.]), np.array([4., 5., 6.])           # the two vectors of the tiny example
print(ch02.expand_indices_str("a_i b_i"))                       # the expander writes the hidden sum of (2.2): a_1*b_1 + a_2*b_2 + a_3*b_3
print(ch02.dot(a, b))                                           # (2.2) evaluated: 32.0 (dot = np.einsum('i,i'))
print(ch02.expand_indices_str("delta_ij u_j"))                  # (2.17): one line per free i — δ collapses the sum to u_i
A = np.arange(1., 10.).reshape(3, 3)                            # the 1…9 matrix of the example
print(ch02.inner(A, A)[0, 1], (A @ A)[0, 1])                    # (2.9)/(2.11): P_12 (Python [0, 1]) = 36.0, and `@` gives the same
print(ch02.classify_indices("x_i C_ij"), ch02.tensor_order("A_ij B_kl"))   # (free, dummy) = (['j'], ['i']); order 4
print(ch02.rename_dummy("x_i C_ij", "i", "k"))                  # renaming the summed letter: x_k C_kj — the same nine numbers
""", explain="""
1. `expand_indices_str` parses the subscripts, finds the repeated letter and writes the hidden sum.
2. `dot` is `np.einsum('i,i')` — (2.2) as one call.
3. δ collapses the sum: for each free i only one term survives.
4. `inner` is (2.9) and equals the `@` operator (P63).
5. `classify_indices` splits free from dummy letters; `tensor_order` counts the free ones.
6. Renaming a dummy changes nothing (N12).
""")
nb.check_agree("""
mine = 0.0                                                   # from scratch: the sum the convention hides
for i in range(3):                                           # i = 1, 2, 3 (Python 0, 1, 2)
    mine += a[i] * b[i]                                      # add a_i b_i
assert np.allclose(mine, ch02.dot(a, b))                     # same 32 → dot is exactly this loop (Ch. 1 primer P15)
P = np.zeros((3, 3))                                         # P_ij = A_ik A_kj, built term by term
for i in range(3):                                           # free index i …
    for j in range(3):                                       # … free index j: the shape of the answer
        for k in range(3):                                   # dummy k: the loop
            P[i, j] += A[i, k] * A[k, j]                     # one term of (2.9)
assert np.allclose(P, ch02.inner(A, A))                      # the library does the same
assert np.allclose(P, np.einsum('ik,kj->ij', A, A))          # and so does einsum (P62)
print("the triple loop IS the summation convention; einsum only writes it faster")
""")
nb.figure("""
from scripts.ch02_matrix_product_highlight import matrix_product_figure   # drawing helper (no physics) for Eq. (2.11)
fig = matrix_product_figure(i=1, j=2, A=A, B=A)                          # book indices: row 1 of A (teal) × column 2 of B (orange) → P_12
plt.show()                                                               # the three boxed products and their sum 36
""", see="Three boxed pairs, one sum: the highlighted row of A and column of B, the three products written out and their total under the cell $P_{12}$.",
   read="The highlighted row and column are the two copies of the dummy k walking together: $A_{1k}$ and $B_{k2}$ for k = 1, 2, 3.",
   change="pick i = 3, j = 1: row 3 and column 1 light up and $P_{31} = 7\\cdot1 + 8\\cdot4 + 9\\cdot7 = 102$ (`matrix_product_figure(i=3, j=1, A=A, B=A)`).")
P("P64", "plotly 3-D arrows, lines and meshes", """
Ch. 1 used `go.Surface` and `go.Scatter3d` points. Arrows are `go.Cone(x, y, z, u, v, w)` (position and direction),
edges are `go.Scatter3d(mode='lines')`, faces are `go.Mesh3d` (vertices + triangles). All stay rotatable on the
published page. `scripts.ch02_drawings.plotly_arrow` bundles one shaft line and one cone into an arrow.
""", code="""
fig = go.Figure(go.Cone(x=[0], y=[0], z=[0], u=[1], v=[0], w=[0], sizemode='absolute', sizeref=0.5))   # one cone at the origin pointing along x
fig.update_layout(height=250, scene_aspectmode='cube', margin=dict(l=0, r=0, t=0, b=0))               # a small cubic scene
fig.show()                                                                                            # drag to rotate
""")
nb.plotly("""
from scripts.ch02_drawings import plotly_arrow, plotly_layout_3d   # drawing helpers (P64): arrow = shaft + cone; equal-aspect layout
x = np.array([1., 2., 3.])                                         # the running position vector OP (any unit)
E = ch02.unit_vectors()                                            # rows e_1, e_2, e_3 (2.1)
assert np.allclose(ch02.vector_from_components(x, E), x)           # (2.1): x = e_1 x_1 + e_2 x_2 + e_3 x_3 rebuilds the same array
fig = go.Figure()                                                  # our own Fig. 2.1
for k in range(3):                                                 # the three unit vectors, teal
    fig.add_traces(plotly_arrow((0, 0, 0), E[k], COLORS["teal"], name="unit vectors e_i", showlegend=(k == 0)))
fig.add_traces(plotly_arrow((0, 0, 0), x, COLORS["accent"], name="x = OP", showlegend=True, width=7))   # the arrow itself, purple
foot = np.array([x[0], x[1], 0.])                                  # the shadow of P on the x_1 x_2 plane
for a_, b_ in (((0, 0, 0), foot), (foot, x), ((x[0], 0, 0), foot), ((0, x[1], 0), foot)):   # dashed projections
    fig.add_trace(go.Scatter3d(x=[a_[0], b_[0]], y=[a_[1], b_[1]], z=[a_[2], b_[2]], mode="lines",
                               line=dict(color=COLORS["muted"], dash="dash"), showlegend=False, hoverinfo="skip"))
plotly_layout_3d(fig, ((-0.3, 3.3),) * 3, title="Eq. (2.1): one arrow, three shadows x₁ = 1, x₂ = 2, x₃ = 3", height=480)
fig.show()                                                         # rotate it: the shadows are the components
""", explain="""
1. `unit_vectors` gives the basis $\\mathbf e_i$ as rows; `vector_from_components` performs (2.1) and returns the same
   array we started from — the equation is a reconstruction rule.
2. The rest is drawing (P64): three teal unit arrows, the purple arrow OP and the dashed shadows that mark $x_1$, $x_2$
   and $x_3$. `N05`: this is our own version of the book's Fig. 2.1.
""")
see_read_change("One purple arrow from the origin to P = (1, 2, 3); its dashed shadows drop onto the $x_1x_2$ plane and then onto the axes; three short teal arrows are the unit vectors.",
                "The numbers (1, 2, 3) are the *shadows* of one arrow on three rulers. Rotate the view: the arrow is the same object however you look at it.",
                "…§2.2 changes the rulers: the same arrow will cast different shadows, and the rule connecting the two sets of shadows is the whole point of the chapter.")
nb.md(r"""
**What would change if…** you meet $\varepsilon_{ijk}u_iv_j$ (§2.7)? The same rule: i and j are summed (nine terms), k
is free — a vector. And $\partial u_i/\partial x_i$ (§2.9) is a sum of three derivatives, a scalar. Every later formula
is read this way; try `ch02.expand_indices_str` on it when in doubt.
""")

# =====================================================================================================================
# A.2 §2.2 Rotation of Axes — C02 (D02, D03), C03 (D01, E1)
# =====================================================================================================================
nb.section("2.2", "Rotation of Axes: Formal Definition of a Vector", intro="""
**What is this section about?** Keep the arrow, turn the rulers. The nine cosines between old and new axes form a
matrix C; the components of *any* vector change by the same rule x' = Cᵀx — and that rule becomes the definition of
a vector.
""")
core("C02", "The direction-cosine matrix C_ij = e_i · e'_j", r"""
Two sets of axes share an origin. What is the one table of numbers that relates them — and what does a single entry
$C_{ij}$ measure?
""")
nb.md(r"""
#### The problem in plain words

A weather model uses axes east/north/up; a wind sensor on a tilted mast reports its own three components; a stress
calculation wants axes normal and tangent to a wall. Same wind, same stress, different numbers. We need the dictionary
between two frames — and it turns out to be nine cosines.
""")
nb.md(r"""
#### The idea

```
old axes  e_1, e_2, e_3  (teal)        new axes  e'_1, e'_2, e'_3  (orange), same origin O
C_ij = e_i · e'_j = cos(angle between old axis i and new axis j)
row i    = "how much of old axis i lies along each new axis"
column j = "the new axis e'_j written in old components"     ← the whole matrix, one picture
```

For a turn by θ about $\mathbf e_3$: $\mathbf e'_1 = (\cos\theta, \sin\theta, 0)$, $\mathbf e'_2 = (-\sin\theta,
\cos\theta, 0)$, $\mathbf e'_3 = \mathbf e_3$ — these are the **columns** of C.
""")
note("N09", "The rotated frame", r"""
O1'2'3' shares the origin; its unit vectors are $\mathbf e'_j$; the vector $\mathbf x$ is *the same arrow* with new
components $x'_j$. **Fig. 2.2** `N17` is our own 3-D drawing below (two frames, one arrow) and the animation after it.
""")
note("N10", "The same vector in the primed basis", r"""
$\mathbf x = x'_1\mathbf e'_1 + x'_2\mathbf e'_2 + x'_3\mathbf e'_3$ — the same vector spelled in the primed basis; in
code `ch02.vector_from_components(xp, E_new)` rebuilds the very same array as `ch02.vector_from_components(x, E_old)`.
""", equation=r"\mathbf x = x'_1 \mathbf e'_1 + x'_2 \mathbf e'_2 + x'_3 \mathbf e'_3", ref="2.3")
P("P65", "orthonormal basis, projection and completeness", r"""
A basis is orthonormal when each vector has length 1 and any two are perpendicular: $\mathbf e_i\cdot\mathbf e_j =
\delta_{ij}$. The component of a vector along a unit vector is the dot product (its projection). **Completeness**: a
vector is the sum of its projections, $\mathbf v = \sum_j (\mathbf v\cdot\mathbf e'_j)\,\mathbf e'_j$ — nothing is
lost, because three perpendicular unit vectors span 3-D space.
""", code="""
E = ch02.rotation_matrix_3d([0, 0, 1], np.deg2rad(30)).T    # rows = the new unit vectors of a frame turned 30° about e_3
print(np.round(E @ E.T, 12))                                # the identity: each row has length 1, rows are perpendicular → orthonormal
v = np.array([1., 2., 3.])                                  # any vector
print(sum((v @ E[j]) * E[j] for j in range(3)))            # [1 2 3]: the sum of the three projections is v again (completeness)
""")
P("P66", "cosines of angles between unit vectors", r"""
For unit vectors the dot product *is* the cosine of the angle between them. Two identities we use: $\cos(\pi/2 -
\theta) = \sin\theta$ and $\cos(\theta + \pi/2) = -\sin\theta$. Python trigonometry works in radians: `np.deg2rad(30)`
= 0.5236.
""", code="""
th = np.deg2rad(30)                                  # 30° in radians
print(np.cos(np.pi/2 - th), np.sin(th))              # 0.5 0.5: the angle between e_2 and e'_1 is 90° − θ
print(np.cos(th + np.pi/2), -np.sin(th))             # -0.5 -0.5: the angle between e_1 and e'_2 is 90° + θ
""")
P("P67", "np.linalg.norm and np.linalg.qr", r"""
`np.linalg.norm(v)` is the length $\sqrt{v_iv_i}$. `np.linalg.qr(M)` factors a matrix into an orthogonal Q and a
triangular R — we use Q of a random matrix as a random rotation (after fixing det Q = +1), which is what
`ch02.random_rotation` does.
""", code="""
rng = np.random.default_rng(0)                       # seeded generator (Ch. 1 primer P10): reproducible "random" numbers
Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))         # Q: the orthogonal factor of a random 3×3 matrix
print(np.round(Q.T @ Q, 12))                         # the identity: Q is orthogonal
print(np.linalg.norm(Q @ [3., 4., 0.]))              # 5.0: an orthogonal matrix keeps lengths (3-4-5 triangle)
""")
nb.md(r"""
#### The maths, step by step

1. Define $C_{ij} \equiv \mathbf e_i\cdot\mathbf e'_j$ (nine numbers; dimensionless).
2. Column j of C is $\mathbf e'_j$ in old components: $(\mathbf e'_j)_i = \mathbf e_i\cdot\mathbf e'_j = C_{ij}$.
3. Row i of C is $\mathbf e_i$ in new components.
4. Because both bases are orthonormal, C is **orthogonal**: $C_{ij}C_{kj} = \delta_{ik}$ and $C_{ji}C_{jk} =
   \delta_{ik}$ — matrix form $\mathbf C\mathbf C^{\rm T} = \mathbf C^{\rm T}\mathbf C = \mathbf I$ (D02 below).
5. Hence $\mathbf C^{-1} = \mathbf C^{\rm T}$ and $\det\mathbf C = +1$ for a rotation (−1 would be a mirror).
""")
D("D02", "Why C is orthogonal: C Cᵀ = Cᵀ C = I and det C = +1", ref="Exercise 2.8",
  goal="""Show that the direction-cosine matrix undoes itself when transposed, so that going back to the old components
  is free (D03), lengths are preserved, and contracted index pairs of C become δ (D18). The book leaves this to
  Exercise 2.8.""",
  assumptions="""Both bases orthonormal and complete in 3-D (steps 2, 4, 5) · "rotation" means reachable from the identity
  by turning (step 8).""",
  start=(r"C_{ij} = \mathbf e_i\cdot\mathbf e'_j", "the definition of the direction cosines (C02)."),
  plan=["Expand an old unit vector in the new basis (completeness).",
        "Dot with another old unit vector: δ on one side, a product of two C's on the other.",
        "Repeat with the roles of the bases swapped.",
        "Take determinants for the sign."],
  uses=[r"completeness of an orthonormal basis, $\mathbf v = (\mathbf v\cdot\mathbf e'_j)\mathbf e'_j$ (P65)",
        "orthonormality of both bases (P65)", "bilinearity of the dot product (it distributes over sums; scalars pull out — D01 step 2 says it again)",
        "matrix product as an index sum (2.9) (N21)", "det(AB) = det A det B (Ch. 1 primer P53)",
        "continuity argument for the sign (gloss in step 8)"],
  steps=[
      ("Write the product we want to evaluate",
       r"C_{ij}C_{kj} = (\mathbf e_i\cdot\mathbf e'_j)(\mathbf e_k\cdot\mathbf e'_j)",
       r"Just the definition of C inserted twice; j is summed. We choose this product because it is $(\mathbf C\mathbf C^{\rm T})_{ik}$.",
       "Row i of C dotted with row k of C."),
      (r"Expand $\mathbf e_i$ in the new basis",
       r"\mathbf e_i = (\mathbf e_i\cdot\mathbf e'_j)\,\mathbf e'_j = C_{ij}\,\mathbf e'_j",
       "Completeness: any vector is the sum of its projections on an orthonormal basis. We expand an *old* axis in *new* axes because that produces exactly a row of C.",
       "An old axis is a sum of new axes weighted by its row of C."),
      (r"Dot with $\mathbf e_k$",
       r"\mathbf e_i\cdot\mathbf e_k = C_{ij}\,(\mathbf e'_j\cdot\mathbf e_k) = C_{ij}C_{kj}",
       r"Bilinearity pulls $C_{ij}$ out; $\mathbf e'_j\cdot\mathbf e_k = \mathbf e_k\cdot\mathbf e'_j = C_{kj}$ by the definition of C. The right side is now the product of step 1.",
       "The same product equals a dot product of two old axes."),
      ("Use orthonormality of the old basis",
       r"C_{ij}C_{kj} = \delta_{ik}\qquad(\mathbf C\mathbf C^{\rm T} = \mathbf I)",
       r"$\mathbf e_i\cdot\mathbf e_k = \delta_{ik}$ for perpendicular unit vectors. This is the first orthogonality relation; in matrix form the summed second indices make $\mathbf C\mathbf C^{\rm T}$.",
       "The rows of C are perpendicular unit vectors."),
      ("Swap the roles: expand a new axis in the old basis",
       r"\mathbf e'_i = (\mathbf e'_i\cdot\mathbf e_j)\,\mathbf e_j = C_{ji}\,\mathbf e_j",
       "Completeness again, now of the old basis. The second relation does not follow from the first without this second argument — a common gap.",
       "A new axis is a sum of old axes weighted by its column of C."),
      (r"Dot with $\mathbf e'_k$ and use orthonormality of the new basis",
       r"C_{ji}C_{jk} = \delta_{ik}\qquad(\mathbf C^{\rm T}\mathbf C = \mathbf I)",
       r"Bilinearity and $\mathbf e_j\cdot\mathbf e'_k = C_{jk}$ on the right; $\mathbf e'_i\cdot\mathbf e'_k = \delta_{ik}$ on the left. Now the *first* indices are summed.",
       "The columns of C are perpendicular unit vectors too."),
      ("Take the determinant of step 6",
       r"(\det\mathbf C)^2 = 1\ \Rightarrow\ \det\mathbf C = \pm1",
       "det(AB) = det A · det B and det Cᵀ = det C (P53), while det I = 1. We take determinants because the sign tells a rotation from a mirror.",
       "C cannot stretch volumes; it may or may not flip them."),
      ("Fix the sign for a rotation",
       r"\det\mathbf C = +1",
       "Turning the axes gradually from θ = 0 changes det C continuously; it starts at det I = +1 and can never pass through 0 (it is ±1), so it stays +1. A mirror (−1) cannot be reached by turning.",
       "A genuine rotation keeps the handedness of the axes."),
  ],
  result=(r"C_{ij}C_{kj} = C_{ji}C_{jk} = \delta_{ik}, \qquad \det\mathbf C = +1",
          "Cᵀ is the inverse of C, and C is a proper rotation."),
  interpret="""A rotation is a relabelling of directions that loses nothing: lengths, angles and volumes are kept; the
  inverse costs a transpose (D03). Every time two C's share a summed index in a formula, they collapse to δ — the engine
  of D06 and D18. It fails for skewed or non-unit axes (then $\\mathbf C^{-1} \\neq \\mathbf C^{\\rm T}$).""",
  check="""θ = 30° in 2-D: row 1 · row 1 = 0.75 + 0.25 = 1, row 1 · row 2 = 0.866(0.5) + (−0.5)(0.866) = 0 ✓; det =
  0.866² + 0.5² = 1 ✓ (`orthogonality_residual` ≈ 1e-17, `is_proper_rotation` True below). Limit: C = I trivially ✓. 200
  random rotations from `random_rotation` all pass (code below).""",
  traps="""proving CᵀC = I and assuming CCᵀ = I "obviously" — for square matrices it is true, but the basis argument
  needs steps 5–6 to say why. A det of −1 passes the orthogonality test and is still not a rotation.""")
note("N15", "Orthogonality of C", r"""
the result of D02 — checked in code with `ch02.is_orthogonal(C)`, `ch02.orthogonality_residual(C)` (≈ 1e-16) and
`ch02.is_proper_rotation(C)` (det +1).
""", equation=r"C_{ij} C_{kj} = C_{ji} C_{jk} = \delta_{ik}, \qquad \det \mathbf C = +1", ref="Exercise 2.8")
nb.md(r"""
The next derivation starts from the forward rule (2.5), $x'_j = x_iC_{ij}$. (2.5) is derived in the next block (C03,
D01); here we only need that it holds.
""")
D("D03", "The inverse transformation x_j = x'_i C_ji", ref="2.7",
  goal="""Get the old components back from the new ones — the book says "it can be shown (Exercise 2.2)".""",
  assumptions="Same as D01/D02: orthonormal bases (enters through D02 in step 2).",
  start=(r"x'_j = x_iC_{ij}", "the forward rule (2.5) of D01 (derived in C03; here we only need that it holds)."),
  plan=["Multiply by a C with the free index and sum.", "Two C's with a shared summed index become δ (D02).",
        "Let δ substitute.", "Rename."],
  uses=["dummy renaming (N12)", r"orthogonality $C_{ij}C_{kj} = \delta_{ik}$ (D02, C02)", "Kronecker substitution (2.17) (N41)"],
  steps=[
      (r"Multiply both sides by $C_{kj}$ and sum on j",
       r"x'_j\,C_{kj} = x_i\,C_{ij}C_{kj}",
       r"Multiplying an identity by the same factor on both sides keeps it true; summing on the now-repeated j too. We choose $C_{kj}$ so that the C-pair on the right can collapse.",
       "Contract the new components with a row of C."),
      ("Collapse the C-pair with D02",
       r"x'_j\,C_{kj} = x_i\,\delta_{ik}",
       r"Orthogonality $C_{ij}C_{kj} = \delta_{ik}$ (second indices summed) from D02 step 4. This is where \"C is orthogonal\" does its work.",
       "The two C's cancel each other."),
      ("Let δ substitute its summed index",
       r"x_k = x'_j\,C_{kj}",
       r"Kronecker substitution (2.17): $x_i\delta_{ik} = x_k$. We now have an old component alone.",
       "An old component is a weighted sum of the new ones."),
      ("Rename the letters to the book's (k → j, j → i)",
       r"x_j = x'_i\,C_{ji}",
       "Renaming a free index on both sides and a dummy anywhere changes nothing (N12); this is (2.7). Note the summed index of C is now its *second* one.",
       "Going back uses C with its indices the other way round: x = C x'."),
  ],
  result=(r"x_j = x'_iC_{ji}, \qquad \mathbf x = \mathbf C\,\mathbf x'",
          r"the inverse of $\mathbf x' = \mathbf C^{\rm T}\mathbf x$ is multiplication by C itself, because $\mathbf C^{-1} = \mathbf C^{\rm T}$."),
  interpret="""Rotating back is as cheap as rotating forward — no matrix inversion. In (2.5) the first index of C is
  summed, in (2.7) the second: the two index placements are the whole difference between C and Cᵀ.""",
  check="""Round trip with the D01 numbers: $x_1 = 1.866(0.866) + 1.232(-0.5) = 1.616 - 0.616 = 1.000$ ✓, $x_2 =
  1.866(0.5) + 1.232(0.866) = 0.933 + 1.067 = 2.000$ ✓ (`inverse_transform_vector`). Limit θ = 0: identity ✓.""",
  traps="""using (2.5) with the indices swapped and calling it the inverse without D02 — it happens to be right only
  because C is orthogonal.""")
note("N14", "Back to the old components", r"""
$x_j = x'_iC_{ji}$ — now the *second* index of C is summed (matrix form $\mathbf x = \mathbf C\mathbf x'$). Round trip
in code: `ch02.inverse_transform_vector(ch02.transform_vector(x, C), C)` returns x to 1e-16 (C03).
""", equation=r"x_j = \sum_{i=1}^{3} x'_i C_{ji} \equiv x'_i C_{ji}", ref="2.7")
nb.worked_example("a 30° turn about e₃", r"""
1. $\cos 30° = 0.866$, $\sin 30° = 0.5$.
2. New axes in old components: $\mathbf e'_1 = (0.866, 0.5, 0)$, $\mathbf e'_2 = (-0.5, 0.866, 0)$, $\mathbf e'_3 = (0, 0, 1)$.
3. Columns → $C = \begin{bmatrix}0.866 & -0.5 & 0\\ 0.5 & 0.866 & 0\\ 0 & 0 & 1\end{bmatrix}$.
4. Check one entry of $C^{\rm T}C$: column 1 · column 1 = $0.866^2 + 0.5^2 = 0.75 + 0.25 = 1$; column 1 · column 2 =
   $0.866\cdot(-0.5) + 0.5\cdot0.866 = 0$ ✓.
5. $\det C = \cos^2 30° + \sin^2 30° = 1$ ✓ (a rotation, not a mirror).
""")
nb.code("""
C = ch02.rotation_matrix_3d([0, 0, 1], np.deg2rad(30))      # passive C for a frame turned 30° about e_3 (Rodrigues' formula; columns = new axes)
print(np.round(C, 3))                                       # the 3×3 of the tiny example
E_old = ch02.unit_vectors()                                 # rows e_1, e_2, e_3
E_new = C.T                                                 # rows e'_1, e'_2, e'_3 = the columns of C
C2 = ch02.direction_cosines(E_old, E_new)                   # C_ij = e_i·e'_j from nine dot products
assert np.allclose(C, C2)                                   # the same matrix both ways
print(ch02.orthogonality_residual(C), ch02.is_proper_rotation(C))   # max|CᵀC − I| ≈ 1e-17, det +1 → True  (D02)
rng = np.random.default_rng(0)                              # seeded generator (Ch. 1 P10)
R = ch02.random_rotation(rng)                               # a random proper rotation for later invariance tests (P67)
print(ch02.is_orthogonal(R), np.round(np.linalg.det(R), 12))   # True 1.0
ok = all(ch02.is_proper_rotation(ch02.random_rotation(rng)) for _ in range(200))   # D02's check: 200 random rotations
print("200 random rotations pass CᵀC = I, det = +1:", ok)
""", explain="""
1. `rotation_matrix_3d` builds the passive C (columns = new axes) for a turn by 30° about $\\mathbf e_3$.
2. `direction_cosines` computes the same C from the nine dot products $\\mathbf e_i\\cdot\\mathbf e'_j$ — the definition.
3. The orthogonality residual and the determinant are D02's two results, as numbers.
4. `random_rotation` (QR of a Gaussian matrix, P67) gives rotations we will use to test invariance later; 200 of them pass D02.
""")
nb.check_agree("""
C_mine = np.zeros((3, 3))                                   # from scratch: nine cosines
for i in range(3):                                          # old axis i
    for j in range(3):                                      # new axis j
        C_mine[i, j] = E_old[i] @ E_new[j]                  # C_ij = e_i · e'_j (a dot product of unit vectors = a cosine, P66)
assert np.allclose(C_mine, ch02.direction_cosines(E_old, E_new))          # same matrix
assert np.allclose(C_mine[:2, :2], ch02.rotation_matrix_2d(np.deg2rad(30)))   # the 2×2 block is the plane rotation matrix
print("nine dot products = the direction-cosine matrix; its upper 2×2 block is the 2-D C")
""")
nb.md(r"""
> ⚠️ **Common confusion (passive vs active):** Wikipedia's rotation matrix $R(\theta)$ rotates a *vector*
> counterclockwise; it has exactly the entries of our C for a frame turned by +θ. But the book uses C *passively*: the
> arrow stays, the axes turn, so the new components are $\mathbf x' = \mathbf C^{\rm T}\mathbf x$ — the transpose.
> Turning the axes by +θ looks, from the axes' point of view, like turning the arrow by −θ.
> `scipy.spatial.transform.Rotation.as_matrix()` is active too. Code: every C-returning function's docstring says
> "passive"; the explainer of C03 has a "what rotates?" toggle that shows both.
""")
nb.plotly("""
from scripts.ch02_fig2_2_rotated_axes import rotated_axes_figure   # drawing helper (P64): two frames, one arrow, both component sets
fig = rotated_axes_figure(C, x=[1, 2, 3])                          # teal old frame, orange new frame (columns of C), the arrow x and its shadows in both
fig.show()                                                         # rotate: the arrow is drawn once
""", explain="""
1. The helper draws the old frame (teal), the new frame (orange, the columns of C) and the single arrow x = (1, 2, 3).
2. The two dashed "staircases" are the components $x_i$ (teal) and $x'_j$ (orange) — the numbers the two observers
   would write down; the title lists both and confirms $|\\mathbf x| = |\\mathbf x'|$ and $\\mathbf C^{\\rm T}\\mathbf C = \\mathbf I$.
""")
see_read_change("Two sets of three axes through one origin; one purple arrow; two dashed staircases climbing to its tip, one per frame.",
                "The arrow is drawn once; two sets of shadows fall from it. Each orange shadow is a weighted sum of the teal ones with the cosines of C as weights (C03).",
                "…C were a turn about $\\mathbf e_1$ instead: the orange 1-axis would stay on the teal one and the 2- and 3-shadows would reshuffle (`ch02.rotation_matrix_3d([1, 0, 0], θ)`).")
nb.md(r"""
🔁 *Reminder (Ch. 1 primer P16):* `animate(update, frames, fig)` calls `update(i)` once per frame; `show_animation(…,
player="frames")` embeds it with ◀ ▮▮ ▶ step buttons.
""")
nb.animation("""
from fluidpy.core.anim import animate                                     # animation helper (Ch. 1 primer P16)
from scripts.ch02_fig2_2_rotated_axes import frame_rotation_frames        # drawing helper: fixed arrow, turning orange axes, live C and bars
n_fr = 10 if not FAST else 6                                              # frames at 0°, 10°, … 90° (FAST: 6 frames)
thetas = np.deg2rad(np.linspace(0, 90, n_fr))                             # the frame angles [rad]
fig, update = frame_rotation_frames(thetas, x=(1.0, 2.0))                 # left: the plane with the arrow (1, 2); right: bars x_i (teal) and x'_j (orange)
fig.set_layout_engine("none")                                             # fixed layout: the automatic layout engine would re-run on every frame
show_animation(animate(update, frames=n_fr, fig=fig, interval=500), player="frames")   # step frame by frame
""", explain="""
1. Each frame turns the orange axes to the next θ and recomputes `ch02.rotation_matrix_2d(θ)` and
   `ch02.transform_vector((1, 2), C)` (2.5).
2. The left panel keeps the arrow fixed; the right panel's teal bars ($x_i$) never move, the orange bars ($x'_j$) do;
   the text box prints the live 2×2 C and max|CᵀC − I|.
""")
see_read_change("The arrow never moves; the orange rulers turn under it. The orange bars change; the teal ones do not; the C box updates at every frame.",
                "At θ = 90° the roles swap: $x'_1 = x_2 = 2$, $x'_2 = -x_1 = -1$ — read it off the last frame's bars.",
                "…you turned by −30° instead: C becomes its transpose (columns and rows swap), and the orange bars go the other way.")
nb.md(r"""
**What would change if…** the second frame were *mirrored* ($\mathbf e'_3 = -\mathbf e_3$)? C is still orthogonal but
det C = −1: not a rotation. §2.7 will show that ε changes sign under such a frame — the book's "rotation" always
means det C = +1.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C03", "How components transform: x'_j = x_i C_ij — the formal definition of a vector", r"""
If the arrow does not change but its numbers do, exactly how do the numbers change — and which triples of numbers
deserve the name "vector"?
""")
nb.md(r"""
#### The problem in plain words

Any three numbers can be stacked into a column. Are (temperature, pressure, density) a vector? Are the three squares
$(x_1^2, x_2^2, x_3^2)$? Physics says no: turn the axes and those triples do not turn like an arrow. The book's answer
is operational — **a vector is anything whose components change under a rotation exactly as the position vector's
do.**
""")
nb.md(r"""
#### The idea

The same arrow casts two shadows: onto the old axes and onto the new. Each new shadow is a weighted sum of the old
ones, the weights being cosines: $x'_j = x_1C_{1j} + x_2C_{2j} + x_3C_{3j}$. Projection does the work: dot $\mathbf x$
with $\mathbf e'_j$ in both spellings of $\mathbf x$.
""")
note("N11", "The first projection", r"""
For θ = 30° and x = (1, 2, 0): $1\times0.866 + 2\times0.5 + 0 = 1.866$.
""", equation=r"\mathbf x\cdot\mathbf e'_1 = x_1\,\mathbf e_1\cdot\mathbf e'_1 + x_2\,\mathbf e_2\cdot\mathbf e'_1 + x_3\,\mathbf e_3\cdot\mathbf e'_1 = x'_1", ref="2.4")
D("D01", "The transformation rule for components", ref="2.5",
  goal="""Find how the three numbers that describe a vector change when the axes are rotated — the rule that will
  *define* what a vector is.""",
  assumptions="""Both bases orthonormal (step 3) · same origin, pure rotation — no translation (start line: the same arrow
  has the same tail).""",
  start=(r"\mathbf x = x_i\mathbf e_i = x'_j\mathbf e'_j",
         "the same arrow written once in the old basis (2.1) and once in the new basis (2.3); only the numbers differ."),
  plan=["Dot both spellings with one new unit vector.", "On the primed side orthonormality kills every term but one.",
        "Name the cosines that remain: that is C.", "Do it for a general new axis, not only e'_1."],
  uses=["dot product (2.2) (N07)", "bilinearity of the dot product (gloss in step 2)",
        r"orthonormality $\mathbf e'_i\cdot\mathbf e'_j = \delta_{ij}$ (P65)", "Kronecker substitution (2.17) (N41)",
        r"the direction cosines $C_{ij} = \mathbf e_i\cdot\mathbf e'_j$ (C02)"],
  steps=[
      ("Write the arrow in both bases",
       r"x_i\,\mathbf e_i = x'_j\,\mathbf e'_j",
       "Equations (2.1) and (2.3) describe one and the same arrow, so their right sides are equal. We put them side by side to compare coefficients.",
       "Two spellings, one vector."),
      (r"Dot both sides with $\mathbf e'_k$",
       r"x_i\,(\mathbf e_i\cdot\mathbf e'_k) = x'_j\,(\mathbf e'_j\cdot\mathbf e'_k)",
       r"The dot product is bilinear: it distributes over the sums and lets the scalars $x_i$, $x'_j$ pull out. We dot with a *new* unit vector to isolate one new component.",
       "Project both spellings onto the new axis k."),
      ("Use orthonormality of the new basis",
       r"x_i\,(\mathbf e_i\cdot\mathbf e'_k) = x'_j\,\delta_{jk}",
       r"The new unit vectors are perpendicular and of length 1, so $\mathbf e'_j\cdot\mathbf e'_k$ is 1 for j = k and 0 otherwise — that is δ. This is the move the book leaves silent.",
       "On the primed side only the k-th term can survive."),
      ("Let δ substitute its summed index",
       r"x_i\,(\mathbf e_i\cdot\mathbf e'_k) = x'_k",
       r"Kronecker substitution (2.17): $x'_j\delta_{jk} = x'_k$, because the sum over j has one nonzero term. We now have one new component alone on one side.",
       "The new component k is a weighted sum of the old ones."),
      ("Name the weights as direction cosines",
       r"x'_k = x_i\,C_{ik},\qquad C_{ik} \equiv \mathbf e_i\cdot\mathbf e'_k",
       "Each weight is the cosine of the angle between old axis i and new axis k (unit vectors); giving them a name turns the projection into a matrix. For k = 1 this line is (2.4).",
       "The weights are the nine cosines between old and new axes."),
      ("Rename the free index k → j",
       r"x'_j = x_i\,C_{ij}",
       r"A free index may be renamed on both sides at once (N12); k = 1, 2, 3 are three equations and this one line holds them all — no \"similarly\" needed.",
       "The new components are the old ones multiplied into C, first index summed."),
  ],
  result=(r"x'_j = x_iC_{ij}, \qquad \mathbf x' = \mathbf C^{\rm T}\mathbf x",
          "each new component is the old components weighted by the cosines in column j of C."),
  interpret="""The rule contains no property of x except its components: **anything** with three components that obey
  it is a Cartesian vector (2.8); anything that does not is a list of numbers. It fails for triples like $(x_1^2, x_2^2,
  x_3^2)$ and for components in a non-orthonormal basis (then δ is replaced by a metric — Appendix B).""",
  check="""Units: both sides carry the unit of x (C is dimensionless) ✓. Limit θ = 0: C = I, x' = x ✓. Number: x = (1, 2),
  θ = 30°: $x'_1 = 1(0.866) + 2(0.5) = 1.866$, $x'_2 = 1(-0.5) + 2(0.866) = 1.232$; lengths 2.236 both ✓
  (`transform_vector` below).""",
  traps="""writing $C_{ji}$ for $C_{ij}$ (row = old axis, column = new axis). Forgetting step 3, which is why only one
  primed term survives. Doing j = 1 and saying "similarly" — one indexed line does all three.""")
note("N13", "The same rule with other letters", r"""
(2.6) $x'_i = x_kC_{ki}$ says the same nine things as (2.5) — `ch02.expand_indices_str` of both prints identical sums
(the §2.3 cell below prints the (2.6) form).
""")
note("N16", "Matrix forms", r"""
$\mathbf x' = \mathbf C^{\rm T}\mathbf x$ (2.5: the first index of C summed = rows of Cᵀ) and $\mathbf x = \mathbf C
\mathbf x'$ (2.7). In code `ch02.transform_vector(x, C)` is `C.T @ x` — asserted below.
""", equation=r"\mathbf x' = \mathbf C^{\rm T}\cdot\mathbf x, \qquad \mathbf x = \mathbf C\cdot\mathbf x'", ref="§2.3")
note("N18", "The formal definition of a vector", r"""
$\mathbf u$ is a Cartesian vector if and only if $u'_j = u_iC_{ij}$ for every rotation. Pass: $\mathbf x$, $3\mathbf x$,
$\mathbf b\times\mathbf x$, $\nabla\phi$. Fail: $(x_1^2, x_2^2, x_3^2)$ and $(|\mathbf x|, 0, 0)$ — the residual table
below shows the failures are not small.
""", equation=r"u'_j = \sum_{i=1}^{3} u_i C_{ij} \equiv u_i C_{ij}", ref="2.8")
nb.worked_example("x = (1, 2) seen from axes turned by 30°", r"""
1. $C = \begin{bmatrix}0.866 & -0.5\\ 0.5 & 0.866\end{bmatrix}$ (columns = new axes).
2. $x'_1 = x_1C_{11} + x_2C_{21} = 1\times0.866 + 2\times0.5 = 1.866$.
3. $x'_2 = x_1C_{12} + x_2C_{22} = 1\times(-0.5) + 2\times0.866 = 1.232$.
4. Length check: $1^2 + 2^2 = 5$ and $1.866^2 + 1.232^2 = 3.482 + 1.518 = 5.000$ ✓ — a rotation cannot change a length.
5. Back with (2.7): $x_1 = x'_1C_{11} + x'_2C_{12} = 1.866\times0.866 + 1.232\times(-0.5) = 1.616 - 0.616 = 1.000$ ✓.
""")
note("N19", "Polar components (Ex. 2.1)", r"""
are (2.5) with j ∈ {r, θ}: $u_r = u_1\cos\theta + u_2\sin\theta$, $u_\theta = -u_1\sin\theta + u_2\cos\theta$, i.e.
$C = \begin{bmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{bmatrix}$ — the 2-D matrix of the tiny
example (stated; the book works it out). With u = (1, 2), θ = 30°: $u_r = 1.866$, $u_\theta = 1.232$; the point's
polar frame *is* our rotated frame. **Fig. 2.3** `N20` is the "polar" mode of the explainer below (our drawing).
""")
nb.code("""
C2 = ch02.rotation_matrix_2d(np.deg2rad(30))                # the passive 2-D C (columns = new axes)
x = np.array([1., 2.])                                      # the arrow (any unit)
xp = ch02.transform_vector(x, C2)                           # (2.5): x'_j = x_i C_ij
print(np.round(xp, 3), np.round(C2.T @ x, 3))               # [1.866 1.232] twice: transform_vector IS Cᵀx (N16)
print(np.round(ch02.inverse_transform_vector(xp, C2), 12))  # (2.7) round trip: [1. 2.]
print(np.round(ch02.polar_components(1., 2., np.deg2rad(30)), 3))   # Ex. 2.1 is the same call: (u_r, u_θ) = (1.866, 1.232)
print(np.linalg.norm(x), np.linalg.norm(xp))                # 2.236 2.236: a rotation keeps lengths (P67)
""", explain="""
1. `rotation_matrix_2d` is the passive C of the tiny example.
2. (2.5) as `transform_vector` and as `C.T @ x` — identical numbers (N16).
3. (2.7) brings the components back exactly (D03).
4. Ex. 2.1's polar components are the same transformation (N19).
5. The lengths agree — D02's orthogonality at work.
""")
nb.check_agree("""
xp_mine = np.array([sum(x[i] * C2[i, j] for i in range(2)) for j in range(2)])   # (2.5) by hand: x'_j = Σ_i x_i C_ij (first index of C summed)
assert np.allclose(xp_mine, ch02.transform_vector(x, C2))                       # the library does exactly this
print(np.round(xp_mine, 3), "= the tiny example's numbers")
""")
nb.md(r"""
🔁 *Reminder (Ch. 1 primers P23, P29):* a `dict` maps names to values, and a `lambda x: …` is a one-line function that
can be stored in it — here a dictionary of candidate "vector formulas" to test.
""")
nb.code("""
rng = np.random.default_rng(1)                                     # seeded (Ch. 1 P10)
R = ch02.random_rotation(rng)                                      # one random rotation C
pts = rng.normal(size=(20, 3))                                     # 20 random points x
b = np.array([0., 0., 1.])                                         # a fixed vector b for the b × x candidate
cands = {"x":         (lambda x: x, ()),                            # the position vector itself
         "3x":        (lambda x: 3 * x, ()),                        # a scalar multiple
         "b×x":       (lambda x, bb: np.cross(bb, x), (b,)),        # b × x — b is a vector parameter, so the test rotates it too
         "x_i²":      (lambda x: x**2, ()),                         # the three squares — NOT a vector
         "(|x|,0,0)": (lambda x: np.array([np.linalg.norm(x), 0, 0]), ())}   # a length stacked into a column — NOT a vector
for name, (fn, params) in cands.items():                           # each candidate: evaluate in the old frame, transform with (2.8), compare with the formula in the new frame
    print(f"{name:10s} residual of the vector test (2.8): {ch02.transforms_as_vector(fn, R, pts, vector_params=params):.2e}")
""", explain="""
1. One random rotation and 20 random points.
2. For each candidate, `transforms_as_vector` evaluates the triple in the old frame and transforms it with (2.8), and
   evaluates the same *formula* in the new frame's coordinates (rotating any vector parameter such as b as well).
3. The residual is the largest difference — zero (to round-off) only for true vectors: the first three pass, the
   squares and the stacked length fail by order one.
""")
nb.figure("""
names = list(cands)                                                                     # the five candidates
res = [ch02.transforms_as_vector(fn, R, pts, vector_params=p) for fn, p in cands.values()]   # their (2.8) residuals
fig, ax = plt.subplots(figsize=(6.5, 3.2))                                              # one panel
cols = [COLORS["teal"] if r < 1e-10 else COLORS["rose"] for r in res]                   # teal = passes, rose = fails
ax.bar(names, np.maximum(res, 1e-16), color=cols)                                       # floor at 1e-16 so zero shows on a log axis
ax.set_yscale("log"); ax.set_ylim(1e-17, 10)                                            # residuals span 17 decades
ax.axhline(1e-12, ls="--", color=COLORS["muted"]); ax.text(-0.4, 3e-12, "round-off level", ha="left", color=COLORS["muted"], fontsize=9)   # the noise floor
ax.set_ylabel("max |u' − Cᵀu| [–]"); ax.set_title("Three formulas are vectors; two are just triples of numbers")
plt.show()
""", see="Three bars at the round-off floor (teal), two bars at order one (rose).",
   read="A vector's formula gives the same arrow whichever axes you compute it in, so its residual is zero; a non-vector's components disagree between observers by numbers of order one — not by a rounding error.",
   change="…you replaced `x**2` by `x**3`: it fails too; `x * (x @ x)` passes (a scalar, $|\\mathbf x|^2$, times a vector is a vector).")
nb.md(r"""
🔁 *Reminder (Ch. 1 primer P18):* `show_viz("ch02", slug)` embeds a chapter explainer (the local file in Jupyter, the
GitHub Pages copy in Colab, a full-window block on the published page).
""")
nb.explainer("rotation_of_axes", heading="Which rotates — the arrow or the ruler?",
             why="""One object, two sets of numbers: only motion shows that the arrow stays put while C, $x'_j$ and (later)
             $\\tau'_{mn}$ move together — and the toggle settles the passive/active confusion once.""",
             tries=["Drag θ to 90° and read C: the columns are the new axes (0, 1) and (−1, 0).",
                    "Flip 'what rotates' to the arrow: the same C now moves the vector the other way (x' = Cx).",
                    "Polar mode: move the point round its circle and watch u_r, u_θ — Ex. 2.1 live.",
                    "Tensor mode, pure shear: find the angle where τ'₁₂ = 0 (45°) — that is §2.11's principal axis."])
nb.md(r"""
**What would change if…** you applied the rule *twice*, once per index? You would get the transformation of a 3×3
array — and that is precisely the definition of a second-order tensor, (2.12), which we reach in §2.4 after meeting the
stress tensor and Cauchy's formula.
""")

# =====================================================================================================================
# A.3 §2.3 Multiplication of Matrices — no A item
# =====================================================================================================================
nb.section("2.3", "Multiplication of Matrices", intro="""
**What is this section about?** The matrix product is the summation convention with one adjacent index — we already
used it: (2.9)–(2.11) in C01 and the matrix forms x' = Cᵀx, x = Cx' in C03. This short section only collects them.
""")
nb.md(r"""
📝 (2.9) $P_{ij} = A_{ik}B_{kj}$ and its boxes (2.11) were expanded in **C01** (N21–N23); the rule "a single dot sums
one index" is (2.10) (N22). (2.6) written with the summed index adjacent, $x'_i = C^{\rm T}_{ik}x_k$, gives the matrix
form of **C03** (N16): $\mathbf x' = \mathbf C^{\rm T}\cdot\mathbf x$, and (2.7) is $\mathbf x = \mathbf C\cdot
\mathbf x'$.
""")
nb.code("""
C = ch02.rotation_matrix_3d([0, 0, 1], np.deg2rad(30)); x = np.array([1., 2., 3.])   # the 30° frame and an arrow
assert np.allclose(ch02.transform_vector(x, C), C.T @ x)        # (2.5) as a matrix product: x' = Cᵀ·x  (N16)
assert np.allclose(ch02.inverse_transform_vector(C.T @ x, C), C @ (C.T @ x))   # (2.7): x = C·x'
assert np.allclose(ch02.inner(C.T, C), np.eye(3))               # CᵀC = I as the (2.9) product `inner`  (D02)
print(ch02.expand_indices_str("C_ki x_k"))                      # the (2.6) form: for each free i, C_1i*x_1 + C_2i*x_2 + C_3i*x_3
""", explain="""
1–2. The two matrix forms of C03, asserted against the index-form functions.
3. $\\mathbf C^{\\rm T}\\mathbf C = \\mathbf I$ written as the (2.9) product `inner`.
4. The adjacent-index form (2.6) expanded: the same nine numbers as (2.5) with other letters (N13).
""")

# =====================================================================================================================
# A.4 §2.4, 2.6 Second-Order Tensors and the Force on a Surface — C04, C05 (D05, E2), C06 (D06)
# =====================================================================================================================
nb.section("2.4, 2.6", "Second-Order Tensors and the Force on a Surface", intro="""
**What is this section about?** The stress tensor: nine numbers that tell the force per area on *any* plane through a
point. We meet its sign convention (C04), derive Cauchy's formula f = n·τ from a tiny tetrahedron (C05, the book's
§2.6), and only then prove the rule that makes τ a tensor, (2.12) (C06). The book states (2.12) first and cites a
tetrahedron argument it never shows; deriving (2.15) first lets us *prove* (2.12) in eight moves.
""")
core("C04", "The stress tensor τ_ij and its sign convention", r"""
Nine stress components sit on a cube. Which index is the face, which is the force — and which way do positive arrows
point?
""")
nb.md(r"""
#### The problem in plain words

Ch. 1 split stress into "normal" (push/pull) and "shear" (slide) on one surface (🔁 Ch. 1 primer P05: stress is force
per area). But at a point in a fluid there are infinitely many surfaces. Pick three perpendicular ones — the faces of a
tiny cube — and write down the force per area on each: three faces × three force directions = nine numbers. That array
is the stress tensor. Ch. 4 will balance forces on this cube to get the equation of motion, so every sign here matters.
""")
nb.md(r"""
#### The idea (rule table)

```
τ_ij :  i = the face (its outward normal is ±e_i)      j = the direction of the force per area
face with outward normal +e_i : positive τ_ij points along +e_j   (tensile normal stress is positive)
face with outward normal −e_i : the arrows reverse (traction −τ_i·)   ← Newton III at a point (R01)
i = j : normal stress (τ_11, τ_22, τ_33)      i ≠ j : shear stress (τ_12 = force along 2 on the face ⊥ 1)
```
""")
nb.recap("R01", "Opposite faces carry equal and opposite stresses", """
Ch. 1 §1.7 balanced pressure on the two faces of a box (D05, P09): at a point the two sides of one plane push on each
other with equal and opposite force — Newton's third law. So the −e₂ face of the book's Fig. 2.4 carries exactly
−(τ₂₁, τ₂₂, τ₂₃). Across a *finite* cube the values differ by a Taylor step, and that difference is Ch. 4's net force.
""", where="Ch. 1 §1.7 C20")
nb.current_core = "C04"   # nb.recap() closes the current CORE block; the design places R01 *inside* C04, so re-open it
nb.md(r"""
🔁 *Reminder (P64 above):* 3-D arrows, edges and faces are `go.Cone`, `go.Scatter3d(mode='lines')` and `go.Mesh3d`.
""")
nb.md(r"""
#### The maths, step by step

1. On the face whose outward normal is $+\mathbf e_2$, the force per area is the vector $(\tau_{21}, \tau_{22},
   \tau_{23})$ [Pa]: the first index names the face, the second the component.
2. On the opposite face ($-\mathbf e_2$) it is $-(\tau_{21}, \tau_{22}, \tau_{23})$ (R01).
3. Collect the three "+" faces as rows:
   $\boldsymbol\tau = \begin{bmatrix}\tau_{11}&\tau_{12}&\tau_{13}\\ \tau_{21}&\tau_{22}&\tau_{23}\\ \tau_{31}&\tau_{32}
   &\tau_{33}\end{bmatrix}$ — **row i = traction on the +e_i face**.
4. Pressure will enter Ch. 4 as $\tau_{ij} = -p\,\delta_{ij} + (\text{viscous})$: a pressure pushes *inward*, so its
   normal stress is negative.
""")
note("N26", "Array vs tensor", r"""
The 3×3 array above is what the code carries everywhere (`tau[i, j]`, Python indices 0–2); it is a *tensor* only
because it obeys (2.12) — C06. **Fig. 2.4** `N25` is our rotatable cube below: nine arrows on the three visible
faces, a toggle (`show_hidden=True`) for the hidden faces with their reversed arrows.
""")
nb.worked_example("a pressure p with one shear a", r"""
τ = [[−p, a, 0], [a, −p, 0], [0, 0, −p]] with p = 3 Pa, a = 1 Pa.
1. Face $+\mathbf e_1$: traction = row 1 = (−3, 1, 0) Pa — pushed *into* the cube by 3 Pa (compression, negative normal
   stress) and dragged along $+\mathbf e_2$ by 1 Pa.
2. Face $-\mathbf e_1$: (+3, −1, 0) Pa (reversed).
3. Face $+\mathbf e_3$: (0, 0, −3) Pa: pure pressure, no shear.
4. $\tau_{12} = \tau_{21} = 1$ Pa: the shear on the face ⊥ 1 along 2 equals the shear on the face ⊥ 2 along 1 — Ch. 4
   shows this symmetry is forced by angular momentum; here we simply chose it.
""")
nb.code("""
p, a = 3.0, 1.0                                                   # pressure and shear level [Pa]
tau = np.array([[-p, a, 0.], [a, -p, 0.], [0., 0., -p]])          # the running stress state of C04, C05 (PRESSURE_SHEAR)
faces = ch02.cube_face_tractions(tau)                             # six faces → outward normal and traction by the sign rule (row i, sign of the normal)
for name, f in faces.items():                                     # keys '+1' … '-3' (book axis numbers)
    print(f"face {name}: n = {f['normal']}, traction = {f['traction']} Pa")
total = sum(f["traction"] for f in faces.values())                # add the six tractions
print("sum over the six faces:", total)                           # [0 0 0]: stresses AT A POINT carry no net force (R01)
print(ch02.stress_component_meaning(2, 3))                        # the rule table's sentence, generated for τ_23
""", explain="""
1. The running stress state: a pressure of 3 Pa on the diagonal (negative = compression) and one shear pair a = 1 Pa.
2. `cube_face_tractions` applies the sign rule: the +e_i face carries row i, the −e_i face carries −row i.
3. The six tractions cancel — at a point, stresses exert no net force (R01); Ch. 4's net force comes from how τ
   *changes* across a finite cube.
4. `stress_component_meaning` says in words what τ₂₃ is: the face ⊥ 2, the force along 3 — a shear stress.
""")
nb.check_agree("""
mine = {}                                                         # from scratch: the sign rule by hand
for i in range(3):                                                # face index i (Python 0–2 = book 1–3)
    e = np.eye(3)[i]                                              # the unit vector e_i
    mine[f"+{i + 1}"] = (e, tau[i])                               # +e_i face: traction = row i of τ
    mine[f"-{i + 1}"] = (-e, -tau[i])                             # −e_i face: the reversed row (R01)
for k in mine:                                                    # every face
    assert np.allclose(mine[k][1], faces[k]["traction"])          # the library applies exactly this rule
print("row i of τ is the traction on the +e_i face; the −e_i face reverses it")
""")
nb.plotly("""
from scripts.ch02_fig2_4_stress_cube import stress_cube_figure   # drawing helper (P64): one cone per non-zero component on the + faces
fig = stress_cube_figure(tau, show_hidden=False)                 # normal stresses blue, shears rose; the book's face letters A–H on hover
fig.show()                                                       # rotate the cube
""", explain="""
1. The helper reads `ch02.cube_face_tractions(tau)` and draws one arrow per non-zero component on each visible + face.
2. Blue arrows are normal stresses (here pointing *into* the cube: compression), rose arrows are shears lying in the
   face; `show_hidden=True` adds the three − faces with reversed grey arrows.
""")
see_read_change("On each visible face one arrow points straight in (the normal stress −p) and one lies in the face (the shear a) — on the +3 face only the inward arrow.",
                "The face's index is the first subscript, the arrow's direction the second; a negative normal stress points *into* the cube.",
                "…you set a = 0: only the three inward pressure arrows remain, equal on every face — an isotropic stress, Ch. 1's pressure (`stress_cube_figure(-3 * np.eye(3))`).")
nb.md(r"""
**What would change if…** you cut the point with a *slanted* plane? None of the nine arrows sits on it — yet the
force per area there is fixed by them. That is Cauchy's formula, next.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C05", "Cauchy's traction formula f = n·τ", r"""
Given the nine stresses on the coordinate faces, what is the force per unit area on a plane with an arbitrary unit
normal $\mathbf n$?
""")
nb.md(r"""
#### The problem in plain words

The drag on a wing, the friction on a river bed, the pressure on a dam: each is the force per area on a surface whose
normal points wherever the surface happens to face. The nine $\tau_{ij}$ were measured on three special planes.
Cauchy's formula says the answer on *any* plane is a matrix–vector product, $f_i = \tau_{ji}n_j$, and every wall force
of the book is $\oint \mathbf n\cdot\boldsymbol\tau\,dA$.
""")
nb.md(r"""
#### The idea

Slice a corner off the cube: a tiny tetrahedron with three coordinate faces (outward normals $-\mathbf e_1, -\mathbf
e_2, -\mathbf e_3$) and one slanted face with normal $\mathbf n$. Balance forces on it (🔁 Ch. 1 primer P09: Newton's
second law on a free-body diagram). Face forces scale like (size)², volume forces (weight, inertia) like (size)³ —
shrink it and only the faces matter. The slanted face's traction must then balance the three coordinate faces'.
""")
P("P68", "limits and orders of smallness", r"""
When an element of size h shrinks, a quantity ∝ h² (an area) shrinks more slowly than one ∝ h³ (a volume); their ratio
h³/h² = h → 0. So in a balance of face terms (∝ h²) and volume terms (∝ h³), the volume terms drop out in the limit —
Ch. 1 §1.7 used the same idea once, glossed. Here it is the hinge of the derivation.
""", code="""
for h in [1e-1, 1e-2, 1e-3]:                 # element size [m]
    print(h, h**3 / h**2)                    # the volume/area ratio is h itself → 0
""")
P("P69", "vector area of a closed surface", r"""
Give every patch of a closed surface the vector n dA (outward normal times area). Their sum is zero: the surface has no
net "direction" (think of a closed box: each face's vector is cancelled by the opposite face's). For the tetrahedron
this gives $\mathbf n\,dA = \mathbf e_1 dA_1 + \mathbf e_2 dA_2 + \mathbf e_3 dA_3$, i.e. $dA_i = n_i\,dA$ — the
coordinate face i is the shadow of the slanted face.
""", code="""
n = np.array([0.6, 0.0, 0.8]); dA = 2.0                     # a unit normal and the slanted face's area [m²]
print(ch02.tetrahedron_face_areas(n, dA))                   # [1.2 0. 1.6] = n * dA: the three shadows
print(n * dA - ch02.tetrahedron_face_areas(n, dA))          # closed-surface sum n dA − Σ e_i dA_i: zero
""")
note("N34", "A surface element is a vector", r"""
$d\mathbf A = \mathbf n\,dA$ — size dA, direction the outward normal. Fluxes (C14) and forces (here) are dot products
with it.
""")
note("N35", "The tetrahedron balance", r"""
in the 1-direction: the book's two lines, which are steps 5–7 of the derivation below (`ch02.tetrahedron_face_areas`).
""", equation=r"f_1\,dA = \tau_{11}\,dA_1 + \tau_{21}\,dA_2 + \tau_{31}\,dA_3, \qquad dA_i = n_i\,dA", ref="§2.6")
D("D05", "Cauchy's traction formula from a shrinking tetrahedron", ref="2.15",
  goal="""Find the force per unit area on a plane with *any* unit normal n from the nine stresses on the three
  coordinate planes.""",
  assumptions="""Stresses continuous at the point (step 3: the face tractions are the point's τ) · body force and
  acceleration finite (step 4) · the tetrahedron shrinks to the point (steps 8–9).""",
  start=(r"\sum\mathbf F = m\,\mathbf a \ \text{ on a tiny tetrahedron}",
         "Newton's second law for the fluid inside a tetrahedron with three faces on the coordinate planes and one slanted face with outward normal n."),
  plan=["List the forces: the slanted face, the three coordinate faces, weight and inertia.",
        "Write the coordinate faces' tractions with C04's sign rule.", "Relate the face areas to n.",
        "Shrink: volume terms vanish faster than face terms."],
  uses=["free-body diagram and Newton's second law (Ch. 1 primer P09)", "the stress sign convention on the −e_j faces (C04)",
        r"vector area of a closed surface, $dA_j = n_j\,dA$ (P69)", "orders of smallness h² vs h³ (P68)", "limits"],
  steps=[
      ("Write the force balance on the tetrahedron",
       r"\mathbf F_{\rm slant} + \mathbf F_1 + \mathbf F_2 + \mathbf F_3 + \mathbf F_{\rm body} = \rho\,dV\,\mathbf a",
       "Newton's second law for the mass ρ dV inside; every force on the boundary and in the volume is listed (free-body diagram, P09).",
       "Four face forces plus the body force equal mass times acceleration."),
      ("Write the slanted face's force",
       r"\mathbf F_{\rm slant} = \mathbf f\,dA",
       r"By definition $\mathbf f$ is the force per unit area on the plane with normal n, and dA is that face's area. This is the unknown we want.",
       "The slanted face carries the traction we are after."),
      ("Write the coordinate faces' forces with the sign rule",
       r"(F_j)_i = -\,\tau_{ji}\,dA_j\quad(j = 1, 2, 3\text{, no sum})",
       r"The face perpendicular to axis j has *outward* normal $-\mathbf e_j$; by C04's rule the traction on a − face is the reversed row j, so its i-component is $-\tau_{ji}$.",
       "Each coordinate face pushes with minus its row of τ."),
      ("Write body force and inertia together",
       r"\mathbf F_{\rm body} - \rho\,dV\,\mathbf a = \rho(\mathbf g - \mathbf a)\,dV",
       "Both are proportional to the mass, hence to the volume dV; we group them because they will share the same fate when the element shrinks.",
       "Weight and inertia scale with the volume."),
      ("Assemble the i-component",
       r"f_i\,dA - \tau_{1i}dA_1 - \tau_{2i}dA_2 - \tau_{3i}dA_3 + \rho(g_i - a_i)\,dV = 0",
       "Steps 2–4 substituted into step 1, component i. The book's first line is this one without the volume term (and for i = 1).",
       "Slanted face minus the three coordinate faces plus volume terms is zero."),
      ("Relate the face areas to the normal",
       r"dA_j = n_j\,dA",
       r"The vector area of a closed surface vanishes: $\mathbf n\,dA - \mathbf e_1dA_1 - \mathbf e_2dA_2 - \mathbf e_3dA_3 = 0$ (P69); reading component j gives this. The coordinate face j is the shadow of the slanted face.",
       "Each coordinate face is the slanted face times a cosine."),
      ("Substitute the areas and write the sum with an index",
       r"f_i\,dA - \tau_{ji}\,n_j\,dA + \rho(g_i - a_i)\,dV = 0",
       r"Step 6 into step 5; the three face terms share the pattern $\tau_{ji}n_j dA$ with j summed (summation convention, C01).",
       "The three faces together give τ contracted with n on its first index."),
      ("Divide by dA and compare sizes",
       r"f_i - \tau_{ji}n_j + \rho(g_i - a_i)\,\frac{dV}{dA} = 0,\qquad \frac{dV}{dA}\propto h",
       "dA ≠ 0. For a tetrahedron of size h, dA ∝ h² and dV ∝ h³, so dV/dA ∝ h (P68); ρ, g, a stay finite.",
       "The volume term is smaller than the face terms by the element's size."),
      ("Let the element shrink to the point",
       r"f_i = \tau_{ji}\,n_j",
       "As h → 0 the last term of step 8 vanishes while the others do not depend on h (stresses continuous at the point). This is (2.15).",
       "The traction on any plane is τ contracted with the normal on its first index."),
  ],
  result=(r"f_i = \tau_{ji}n_j, \qquad \mathbf f = \mathbf n\cdot\boldsymbol\tau",
          "the nine numbers on the coordinate planes determine the force per area on every plane through the point."),
  interpret="""Stress is nine numbers, not a vector, precisely because the force per area depends on the cut; yet those
  nine suffice for every cut. Every wall force in the book is $\\oint\\mathbf n\\cdot\\boldsymbol\\tau\\,dA$. The formula
  contracts the *first* index; $\\boldsymbol\\tau\\cdot\\mathbf n$ equals it only when τ is symmetric (Ch. 4 proves that for
  the stress). It fails where stresses are discontinuous (a shock, an interface with surface tension — then the jump
  conditions of Ch. 4 apply).""",
  check="""Units: Pa on both sides ✓. Limit n = e₃: $f_i = \\tau_{3i}$, the traction on the coordinate face — the formula
  returns its inputs ✓. Number (Ex. 2.2, a = 1 Pa, φ = 30°): $f = (0\\cdot0.866 + 1\\cdot0.5,\\ 1\\cdot0.866 + 0\\cdot0.5)
  = (0.5, 0.866)$ Pa ✓ (`traction`, `example_2_2` below). Pure pressure τ = −pδ: f = −p n for every n (Ch. 1's isotropic
  pressure) ✓.""",
  traps="""forgetting the minus signs of step 3 (outward normal −e_j). Dropping the volume terms without saying why (they
  are O(h³) against O(h²)). Writing $\\tau_{ij}n_j$ — that needs symmetry.""")
note("N37", "Contract with the normal", r"""
(2.15) $\mathbf f = \mathbf n\cdot\boldsymbol\tau$ is to τ what $u_n = \mathbf u\cdot\mathbf n$ is to $\mathbf u$ —
contract with the normal — except that $\mathbf f$ is a vector.
""")
nb.md(r"""
> ⚠️ **Common confusion (which index?):** $f_i = \tau_{ji}n_j$ contracts the *first* index (row = face).
> $\boldsymbol\tau\cdot\mathbf n$ contracts the second and equals $\mathbf n\cdot\boldsymbol\tau$ only when τ is symmetric
> — true for the stress (Ch. 4), false for a general tensor; the code never assumes it (`traction` uses
> `einsum('ji,j->i')`).
""")
nb.worked_example("Ex. 2.2 with a = 1 Pa, φ = 30°", r"""
τ = [[0, a], [a, 0]] (pure shear in a channel, $x_1$ along the flow).
1. $\mathbf n = (\cos30°, \sin30°) = (0.866, 0.5)$.
2. $f_1 = \tau_{11}n_1 + \tau_{21}n_2 = 0 + 1\times0.5 = 0.5$ Pa.
3. $f_2 = \tau_{12}n_1 + \tau_{22}n_2 = 1\times0.866 + 0 = 0.866$ Pa.
4. $|\mathbf f| = \sqrt{0.25 + 0.75} = 1 = |a|$.
5. Direction: $\theta = \arctan2(0.866, 0.5) = 60°$; for a = −1 Pa both components flip → 240°.
6. Normal part $\sigma_n = \mathbf f\cdot\mathbf n = 0.5\times0.866 + 0.866\times0.5 = 0.866 = a\sin60°$; shear part
   $\tau_s = |\mathbf f - \sigma_n\mathbf n| = 0.5 = a\cos60°$ — in general $\sigma_n = a\sin2\phi$, $\tau_s = a\cos2\phi$
   (stated; the Mohr-circle preview of C13).
""")
P("P70", "np.arctan2", r"""
`np.arctan2(y, x)` returns the angle of the point (x, y) in the correct quadrant (−π, π]; `np.arctan(y/x)` cannot tell
(−1, −1) from (1, 1). Add 2π (or take `% 360` in degrees) when you want 0–360°.
""", code="""
print(np.rad2deg(np.arctan2(0.866, 0.5)))                   # 60.0: the direction of f for a > 0
print(np.rad2deg(np.arctan2(-0.866, -0.5)) % 360)           # 240.0: the direction of f for a < 0 (opposite quadrant)
""")
note("N38", "Ex. 2.2, stated", r"""
the numbers of the tiny example, plus the (2.12) route the book also takes: $\tau'_{11} = \sqrt3 a/2 = 0.866a$,
$\tau'_{12} = a/2$ — the same two numbers as $\sigma_n$ and $\tau_s$, because the rotated frame's 1'-axis *is* n;
repeated inside C06. **Fig. 2.6** `N39` (channel, element at 30°) is the "channel shear" preset of the explainer and
the slider figure below.
""")
nb.code("""
tau2 = ch02.shear_flow_stress(1.0)                               # Ex. 2.2's 2×2 pure shear [[0, a], [a, 0]], a = 1 Pa
n = np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))])   # the unit normal at 30° from the flow direction
f = ch02.traction(tau2, n)                                       # (2.15): f_i = τ_ji n_j  (einsum 'ji,j->i')
print(np.round(f, 4))                                            # [0.5 0.866] Pa
res = ch02.example_2_2(1.0, np.deg2rad(30))                      # Ex. 2.2 packaged: magnitude, direction, normal and shear parts
print(np.round(res["magnitude"], 4), np.round(np.rad2deg(res["angle_rad"]), 1), np.round(res["sigma_n"], 4), np.round(res["tau_s"], 4))   # 1.0 60.0 0.866 0.5
res_neg = ch02.example_2_2(-1.0, np.deg2rad(30))                 # the other half of the channel: a < 0
print(np.round(np.rad2deg(res_neg["angle_rad"]), 1))             # 240.0: f flips through the origin
f3 = ch02.traction(tau, [0.6, 0, 0.8])                           # the same formula on C04's 3-D state with n = (0.6, 0, 0.8)
print(f3)                                                        # [-1.8 0.6 -2.4] Pa = (−3·0.6, 1·0.6, −3·0.8)
""", explain="""
1. `shear_flow_stress` builds the 2×2 stress of a channel flow at one point.
2. The unit normal at 30°.
3. (2.15) by `einsum('ji,j->i')` — the first index of τ is contracted.
4. `example_2_2` returns the magnitude (= |a|), the direction (60°, P70) and the normal/shear split.
5. a < 0 flips the direction to 240°.
6. The same formula on the 3-D state of C04: each component is a row-wise weighted sum.
""")
nb.check_agree("""
f_mine = np.array([sum(tau2[j, i] * n[j] for j in range(2)) for i in range(2)])   # (2.15) by hand: f_i = Σ_j τ_ji n_j
assert np.allclose(f_mine, ch02.traction(tau2, n))                               # the library contracts the first index
tau_ns = np.array([[0., 2.], [0., 0.]])                                          # a NON-symmetric 2×2 (not a stress; a test)
print(ch02.traction(tau_ns, n), tau_ns @ n)                                      # [0 1.732] vs [1 0]: n·τ ≠ τ·n
print("different: the index order is not a formality; only a symmetric τ lets you forget it")
""")
nb.plotly("""
from scripts.ch02_fig2_5_tetrahedron import tetrahedron_figure   # drawing helper (P64): Fig. 2.5 as a force balance
fig = tetrahedron_figure(tau, n=[0.6, 0, 0.8])                   # coordinate faces with reversed arrows; slanted face with f = (−1.8, 0.6, −2.4) Pa
fig.show()                                                       # rotate
""", explain="""
1. The three coordinate faces carry the reversed rows of τ (C04's − faces); the slanted face carries
   $\\mathbf f = \\mathbf n\\cdot\\boldsymbol\\tau$ as one bold arrow, with n as a thin arrow.
2. The face areas are `ch02.tetrahedron_face_areas(n, dA)` = (0.6, 0, 0.8) dA — D05 step 6.
""")
see_read_change("Three coordinate faces (shadows) and one slanted face; the bold arrow on the slanted face is the balance of the three faces' arrows.",
                "f is not along n: it has a normal part (blue, along n) and a shear part (rose, in the plane) — here $\\sigma_n = \\mathbf f\\cdot\\mathbf n = -3$ Pa (the pressure) and a shear of 0.6 Pa.",
                "…n = $\\mathbf e_3$: f = row 3 = (0, 0, −3) Pa — Cauchy's formula returns the coordinate face (`tetrahedron_figure(tau, n=[0, 0, 1])`).")
nb.md(r"""
🔁 *Reminders (Ch. 1 primers P17, P06):* `slider_figure(fn, name, values)` precomputes `fn(value) → {trace: (x, y)}` for
every slider position so the plotly figure keeps working on the published page; `np.linspace(a, b, n)` gives n evenly
spaced values.
""")
nb.plotly("""
from fluidpy.core.interact import slider_figure                          # precomputed-slider plotly figures (Ch. 1 P17)
sq = np.array([[-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1]]) * 0.5             # the square element's outline (drawing units)
def cut(phi_deg):                                                        # one slider position: the cut with normal at φ
    d = ch02.traction_2d(tau2, np.deg2rad(phi_deg))                      # (2.15) + the split: n, f, σ_n, τ_s for a = +1 Pa
    n_, f_ = d["n"], d["f"]                                              # unit normal and traction [Pa]
    sn = d["sigma_n"] * n_                                               # the normal part σ_n n (blue)
    return {"element": (sq[0], sq[1]),                                   # the square (fixed)
            "cut (plane ⊥ n)": (np.array([-n_[1], n_[1]]) * 0.7, np.array([n_[0], -n_[0]]) * 0.7),   # the cut line through the centre
            "n (unit normal)": ([0, n_[0]], [0, n_[1]]),                 # n
            "f = n·τ  (|f| = a)": ([0, f_[0]], [0, f_[1]]),              # the traction
            "σ_n n = a sin2φ · n": ([0, sn[0]], [0, sn[1]]),             # normal part
            "shear part, |·| = |a cos2φ|": ([sn[0], f_[0]], [sn[1], f_[1]])}   # from the tip of σ_n n to the tip of f
phis = np.linspace(0, 180, 37 if not FAST else 25)                       # 37 slider positions (FAST: 25)
fig = slider_figure(cut, "φ", phis, unit="°", xlabel="x₁ (flow direction)", ylabel="x₂",
                    title="Ex. 2.2: turn the cut — f turns twice as fast as n (a = 1 Pa)",
                    xrange=[-1.3, 1.3], yrange=[-1.3, 1.3], height=520, active=6)   # start at φ = 30°
fig.update_yaxes(scaleanchor="x", scaleratio=1)                          # equal axes so right angles look right
fig.show()                                                               # drag φ
""", explain="""
1. `traction_2d` evaluates (2.15) for the normal at angle φ and splits f into $\\sigma_n\\mathbf n$ and the shear part.
2. Each slider position redraws the cut, n, f, the blue normal part and the rose shear part (the leg from the tip of
   $\\sigma_n\\mathbf n$ to the tip of f). `N36`: this is our version of the book's traction sketch.
""")
see_read_change("A square element, a cut line through its centre, the normal n and the traction f drawn from the centre with its normal and shear legs.",
                "f turns twice as fast as n: at 45° f is along n (pure normal stress a, no shear); at 0° and 90° f is perpendicular to n (pure shear). $\\sigma_n = a\\sin2\\phi$, $\\tau_s = a\\cos2\\phi$.",
                "…a < 0 (the other half of the channel): every f arrow is mirrored through the origin — `ch02.traction_2d(ch02.shear_flow_stress(-1.0), φ)`.")
nb.explainer("cauchy_traction_principal_axes", heading="Cut the point any way you like: what pushes on the cut?",
             why="""The traction is a function of the cut direction; only dragging n and watching f, $\\sigma_n$, $\\tau_s$
             and the Mohr point move together makes (2.15) a picture — and lets you *find* the shear-free planes before
             §2.11 proves they are eigenvectors.""",
             tries=["Channel-shear preset: drag φ to 45° — the shear vanishes and σ_n = a.",
                    "Hydrostatic preset: f stays along n whatever φ; every plane is principal.",
                    "Toggle 'contract second index' on the non-symmetric demo: f changes — the index order matters.",
                    "Read the Explain tab at φ = 30°: f₁ = 0·0.866 + 1·0.5 = 0.5 Pa, term by term."])
nb.md(r"""
**What would change if…** a second observer used rotated axes? Both $\mathbf f$ and $\mathbf n$ are vectors (2.8),
and (2.15) must hold for both observers — that forces the nine $\tau_{ij}$ to transform in one particular way. That
way is the definition of a tensor: C06.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C06", "The transformation rule of a second-order tensor: τ' = Cᵀ τ C", r"""
How do the nine stress components change when the axes turn — and why does *that* rule, not the array itself, define
a tensor?
""")
nb.md(r"""
#### The problem in plain words

The same stress state must give the same force on the same plane whichever axes you compute in. Two observers, two
3×3 arrays. The only way they can agree about every plane is if the arrays are related by one C per index:
$\tau'_{mn} = C_{im}C_{jn}\tau_{ij}$. Ch. 3 will need it for the strain rate, Ch. 4 for the stress, Ch. 12 for the
Reynolds stress — it is the working definition of "tensor".
""")
nb.md(r"""
#### The idea

(2.8) once per index. A vector: $u'_n = u_iC_{in}$. A tensor has two slots, each behaves like a vector: $\tau'_{mn} =
C_{im}\,C_{jn}\,\tau_{ij}$. Matrix form: $\boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau\,\mathbf C$ — Cᵀ on the left
for the first index, C on the right for the second.

The book states (2.12) and points to a tetrahedron argument (Sommerfeld); we have the tetrahedron result (2.15) in
hand, so we derive (2.12) from it.
""")
D("D06", "The tensor transformation rule from Cauchy's formula", ref="2.12",
  goal="""Show how the nine stress components must change under a rotation of axes, given only that force per area and
  the normal are vectors. The book states (2.12) and cites Sommerfeld's tetrahedron; we have the tetrahedron result
  (2.15) already, so we derive (2.12) from it.""",
  assumptions="""τ is *defined* in each frame by (2.15) (start and step 5) · C orthogonal (enters through (2.7) in
  step 3).""",
  start=(r"f_i = \tau_{ji}n_j \ \text{(old frame)}, \qquad f'_n = \tau'_{mn}n'_m \ \text{(new frame)}",
         "Cauchy's formula holds for every observer, each with their own components of τ."),
  plan=["Transform f as a vector.", "Insert Cauchy's formula in the old frame.",
        "Express the old n through the new n'.", "Compare with Cauchy's formula in the new frame, for every n'."],
  uses=["f and n are vectors, (2.8) (C03)", "Cauchy's formula (2.15) (C05)", "the inverse transformation (2.7) (D03)",
        "dummy renaming (N12)", "\"true for every n' ⇒ the coefficients agree\" (gloss in step 7)"],
  steps=[
      ("Transform the traction as a vector",
       r"f'_n = f_i\,C_{in}",
       "f is a physical force per area, hence a vector, so its components obey (2.8) with the first index of C summed. We start from f because both frames must agree on it.",
       "The new observer's traction components come from the old ones by C."),
      ("Insert Cauchy's formula in the old frame",
       r"f'_n = \tau_{ji}\,n_j\,C_{in}",
       r"(2.15) gives $f_i$ in terms of the old τ and the old n; substituting expresses the new traction through old quantities.",
       "The new traction, written with the old stresses and the old normal."),
      ("Write the old normal through the new one",
       r"n_j = n'_m\,C_{jm}",
       "n is a unit vector, so the inverse rule (2.7) applies (second index of C summed). We need n' because the new observer measures n'.",
       "The old components of the normal from the new ones."),
      ("Substitute the normal",
       r"f'_n = C_{jm}\,C_{in}\,\tau_{ji}\,n'_m",
       "Step 3 into step 2; the factors are reordered (scalars commute) so that both C's sit in front. Everything on the right is now old τ, two C's and the new normal.",
       "The new traction equals two C's times the old τ times the new normal."),
      ("Write Cauchy's formula in the new frame",
       r"f'_n = \tau'_{mn}\,n'_m",
       "The new observer sees the same physics, so (2.15) holds with primed components everywhere — this *defines* the new observer's τ'.",
       "The new observer's own stress tensor gives the same traction."),
      ("Equate steps 4 and 5",
       r"\big(\tau'_{mn} - C_{jm}C_{in}\tau_{ji}\big)\,n'_m = 0",
       r"Both are the same $f'_n$; subtracting and factoring the common $n'_m$ (a sum on m) collects the difference. The bracket does not depend on n'.",
       "A certain matrix, applied to any unit normal, gives zero."),
      ("Conclude the bracket vanishes",
       r"\tau'_{mn} = C_{jm}\,C_{in}\,\tau_{ji}",
       "Step 6 holds for every unit n'; choose n' = e'_1, e'_2, e'_3 in turn and each picks out one column of the bracket, which must be zero (compare coefficients).",
       "The new stress components are fixed by the old ones and C."),
      ("Rename the dummies (j → i, i → j)",
       r"\tau'_{mn} = C_{im}\,C_{jn}\,\tau_{ij}\qquad(\boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau\,\mathbf C)",
       r"Swapping the names of two summed letters changes nothing (N12); this is the book's (2.12). In matrix form $C_{im} = (\mathbf C^{\rm T})_{mi}$, so Cᵀ stands on the left and C on the right.",
       "One C per index — the vector rule applied twice."),
  ],
  result=(r"\tau'_{mn} = C_{im}C_{jn}\tau_{ij}, \qquad \boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau\,\mathbf C",
          "a second-order tensor transforms with one direction-cosine matrix per index; whatever obeys this is a tensor."),
  interpret="""The rule contains nothing about stress: any nine numbers that obey it are the components of a physical
  object, and any that do not are not (N27). Ch. 3's strain rate, Ch. 4's stress and Ch. 12's Reynolds stress all pass.
  The book's order (2.12 stated first, 2.15 derived later) is reversed here on purpose.""",
  check="""Units: Pa ✓. Limit C = I: τ' = τ ✓. Number (Ex. 2.2, a = 1, θ = 30°): $\\tau'_{11} = 2(0.866)(0.5) = 0.866$,
  $\\tau'_{12} = 0.75 - 0.25 = 0.5$, $\\tau'_{22} = -0.866$ ✓ (`transform_tensor` = `C.T @ tau @ C` below). Trace 0
  before and after (D18).""",
  traps="""writing $\\mathbf C\\boldsymbol\\tau\\mathbf C^{\\rm T}$ — that is the *active* convention; ours is passive.
  Forgetting the final rename and comparing $C_{jm}C_{in}\\tau_{ji}$ with the book. Assuming symmetry of τ anywhere — not
  needed.""")
note("N27", "Tensor ≠ matrix", r"""
Any 3×3 array is a matrix; its entries are the components of a tensor only if they obey (2.12) in every frame. The code
below declares an array "the same in every frame" and measures the (2.12) residual: not small.
""")
note("N28", "Fourth order", r"""
one C per index, four of them, $3^4 = 81$ components — `ch02.transform_tensor(A4, C)` handles any order; the viscosity
tensor of Ch. 4 §4.5 is the one we will meet.
""", equation=r"A'_{mnpq} = C_{im} C_{jn} C_{kp} C_{lq} A_{ijkl}", ref="2.13")
note("N29", "Examples that pass", r"""
the stress $\tau_{ij}$, the velocity gradient $\partial u_i/\partial x_j$ (N52), and the outer product $u_iv_j$ of two
vectors (Exercise 2.10) — asserted below with `ch02.outer`.
""")
nb.worked_example("Ex. 2.2's rotated stress, a = 1 Pa, θ = 30°", r"""
1. $C = \begin{bmatrix}0.866 & -0.5\\ 0.5 & 0.866\end{bmatrix}$.
2. $\tau'_{11} = C_{i1}C_{j1}\tau_{ij} = C_{11}C_{21}\tau_{12} + C_{21}C_{11}\tau_{21} = 2\times0.866\times0.5\times1 = 0.866$ Pa.
3. $\tau'_{12} = C_{11}C_{22}\tau_{12} + C_{21}C_{12}\tau_{21} = 0.866^2 - 0.5^2 = 0.75 - 0.25 = 0.5$ Pa.
4. $\tau'_{22} = 2C_{12}C_{22}\tau_{12} = -2\times0.5\times0.866 = -0.866$ Pa.
5. So $\tau' = \begin{bmatrix}0.866 & 0.5\\ 0.5 & -0.866\end{bmatrix}$: the 1'-axis is n, so $\tau'_{11} = \sigma_n$
   and $\tau'_{12} = \tau_s$ of C05 ✓; the trace is still 0 (C07).
""")
nb.code("""
C2 = ch02.rotation_matrix_2d(np.deg2rad(30))                        # the passive 2-D C
tp = ch02.transform_tensor(tau2, C2)                                # (2.12): τ'_mn = C_im C_jn τ_ij  (einsum 'im,jn,ij->mn')
print(np.round(tp, 4))                                              # [[0.866 0.5], [0.5 −0.866]] Pa
print(np.round(C2.T @ tau2 @ C2, 4))                                # the matrix form CᵀτC — identical
print(np.round(ch02.example_2_2(1.0, np.deg2rad(30))["tau_rot"], 4))   # Ex. 2.2's own route — identical again (N38)
R = ch02.random_rotation(np.random.default_rng(2))                  # a random 3-D rotation
print(ch02.transforms_as_tensor(lambda x, bb: ch02.outer(x, np.cross(bb, x)), R, pts, vector_params=(b,)))   # outer product of two vector fields x and b × x (b rotated too): ≈ 1e-15 → a tensor (N29)
print(ch02.transforms_as_tensor(lambda x: np.diag([1., 2., 3.]), R, pts))   # a fixed array declared "the same in every frame": order 1 → NOT a tensor (N27)
A4 = np.random.default_rng(3).normal(size=(3, 3, 3, 3))            # 81 random numbers
print(ch02.transform_tensor(A4, R).shape)                           # (3, 3, 3, 3): (2.13), four C's  (N28)
""", explain="""
1. (2.12) as `transform_tensor` — an einsum with one C per index.
2. The matrix form $\\mathbf C^{\\rm T}\\boldsymbol\\tau\\mathbf C$ — identical numbers.
3. Ex. 2.2's route gives the same τ'.
4. The outer product $x_i(\\mathbf b\\times\\mathbf x)_j$ of two vector fields passes the tensor test (residual at
   round-off); as in C03, the fixed vector b is handed over as a `vector_params` entry so that the test rotates it too.
5. A fixed array "the same in every frame" fails by order one — a matrix, not a tensor.
6. (2.13): four C's, 81 numbers.
""")
nb.check_agree("""
tp_mine = np.zeros((2, 2))                                          # (2.12) by hand
for m in range(2):                                                  # free m
    for n_ in range(2):                                             # free n
        for i in range(2):                                          # dummy i
            for j in range(2):                                      # dummy j
                tp_mine[m, n_] += C2[i, m] * C2[j, n_] * tau2[i, j]   # one term of C_im C_jn τ_ij
assert np.allclose(tp_mine, ch02.transform_tensor(tau2, C2))        # the library's einsum is this quadruple loop
assert np.allclose(tp_mine, C2.T @ tau2 @ C2)                       # and so is the matrix form
print(np.round(tp_mine, 4), "= the tiny example")
""")
nb.figure("""
thetas = np.deg2rad(np.linspace(0, 180, 181))                       # frame angles 0…180°
presets = {"pure shear a = 1 (Ex. 2.2)": tau2, "uniaxial [[1, 0], [0, 0]]": np.array([[1., 0.], [0., 0.]])}   # two stress states
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)       # one panel per preset
for ax, (name, T) in zip(axs, presets.items()):                     # for each preset …
    tps = np.array([ch02.transform_tensor(T, ch02.rotation_matrix_2d(t)) for t in thetas])   # … τ'(θ) by (2.12) at every angle
    ax.plot(np.rad2deg(thetas), tps[:, 0, 0], color=COLORS["blue"], label="τ'₁₁ (normal)")     # normal stress on the 1' face
    ax.plot(np.rad2deg(thetas), tps[:, 1, 1], color=COLORS["blue"], ls="--", label="τ'₂₂ (normal)")
    ax.plot(np.rad2deg(thetas), tps[:, 0, 1], color=COLORS["rose"], label="τ'₁₂ (shear)")      # shear stress in the new frame
    k30 = 30; ax.plot([30] * 3, tps[k30, [0, 1, 0], [0, 1, 1]], "o", color=COLORS["ink"], ms=4)   # the θ = 30° values of the tiny example
    ax.axhline(0, color=COLORS["muted"], lw=0.8); ax.set_xlabel("frame rotation θ [°]"); ax.set_title(name, fontsize=10)
axs[0].set_ylabel("stress component [Pa]"); axs[0].legend(fontsize=8)
fig.suptitle("Three cosine curves of period 180°: the components move, the state does not")
plt.show()
""", see="Three cosine curves of period 180° in each panel; the dots at θ = 30° are the tiny example's 0.866, 0.5 and −0.866 Pa.",
   read="τ'₁₂ crosses zero at 45° for pure shear (and at 0° for the uniaxial state) — a principal axis (C13); the shear curve's amplitude is the largest shear any plane can carry.",
   change="…you added a pressure (−p on the diagonal): both normal curves shift down by p, the shear curve does not move.")
nb.plotly("""
def mohr(theta_deg):                                                 # one slider position: the frame turned by θ
    C_ = ch02.rotation_matrix_2d(np.deg2rad(theta_deg))              # passive C
    T = ch02.transform_tensor(tau2, C_)                              # τ' = CᵀτC (2.12) for pure shear a = 1 Pa
    c, r = ch02.mohr_circle_2d(tau2)                                 # Mohr's circle: centre I₁/2 = 0, radius 1 Pa
    ang = np.linspace(0, 2 * np.pi, 181)                             # to draw the circle
    return {"Mohr circle (all planes)": (c + r * np.cos(ang), r * np.sin(ang)),          # fixed circle
            "(τ'₁₁, τ'₁₂): the 1' face": ([c, T[0, 0]], [0, T[0, 1]]),                    # radius to the current point
            "(τ'₂₂, −τ'₁₂): the 2' face": ([c, T[1, 1]], [0, -T[0, 1]])}                   # its antipode: the perpendicular face
fig = slider_figure(mohr, "θ", np.linspace(0, 180, 37 if not FAST else 25), unit="°", xlabel="normal stress σ_n [Pa]",
                    ylabel="shear stress τ_s [Pa]", title="Mohr's circle: as the frame turns by θ, the point runs round at 2θ",
                    xrange=[-1.3, 1.3], yrange=[-1.3, 1.3], height=480, active=6)   # start at 30°
fig.update_yaxes(scaleanchor="x", scaleratio=1)                      # a circle should look like a circle
fig.show()                                                           # drag θ
""", explain="""
1. For each θ the stress is rotated with (2.12) and the pair $(\\tau'_{11}, \\tau'_{12})$ is plotted as a point.
2. `mohr_circle_2d` gives the circle every such point lies on: centre $I_1/2$ (here 0), radius 1 Pa. The second radius
   marks the perpendicular face $(\\tau'_{22}, -\\tau'_{12})$ — always the antipode.
""")
see_read_change("A circle of radius 1 Pa centred at the origin; two radii from the centre to opposite points, turning as you drag θ.",
                "One turn of the frame (180°) is one full turn of the point (360°): the double angle of the cosine curves above. The circle's radius is the largest shear any plane can carry; where the radius hits the σ axis (θ = 45°) the shear is zero — a principal axis (C13).",
                "…you added a pressure −p: the circle slides left by p; its radius (the shear) is unchanged.")
nb.md(r"""
**What would change if…** you contracted the two indices of τ' with each other? $\tau'_{mm} = C_{im}C_{jm}\tau_{ij} =
\delta_{ij}\tau_{ij} = \tau_{ii}$: the trace is the same in every frame. That is the next section — contraction makes
invariants.
""")

# =====================================================================================================================
# A.5 §2.5 Contraction and Multiplication — C07 (D18)
# =====================================================================================================================
nb.section("2.5", "Contraction and Multiplication", intro="""
**What is this section about?** Setting two indices equal and summing (contraction) lowers the order by two;
multiplying raises it. Contracted quantities with no free index are the same in every frame — the invariants of a
tensor — and the einsum patterns here are used in every later notebook.
""")
core("C07", "Contraction: the trace, the three invariants and the double dot", r"""
Which combinations of the nine numbers $\tau_{ij}$ are the same for every observer — and why exactly those?
""")
nb.md(r"""
#### The problem in plain words

The pressure in a flowing fluid is defined in Ch. 4 as $-\tau_{ii}/3$; the rate of volume expansion in Ch. 3 is
$S_{ii}$; turbulence models in Ch. 12 are built from invariants of the Reynolds stress. Each is a number that must not
depend on how you set up your axes. Contraction is the recipe that produces such numbers.
""")
nb.md(r"""
#### The idea

```
A_ij  →  set j = i and sum  →  A_ii = A_11 + A_22 + A_33        (a chain of indices closing on itself)
under a rotation each index brings one C; a closed pair  C_im C_jm  is  δ_ij  (D02)  →  the C's vanish
A_ij A_ji = a closed chain of two  →  invariant, and it is the one that enters I₂
A_ij A_ij = tr(A Aᵀ) is also invariant*  — but it is NOT the λ-coefficient of the characteristic polynomial
```

(*both scalars survive a rotation — $\mathbf C^{\rm T}\mathbf A\mathbf C$ leaves $\mathrm{tr}(\mathbf A\mathbf A^{\rm T})$
unchanged too. The distinction is that $I_2 = \tfrac12(I_1^2 - A_{ij}A_{ji})$ needs the closed chain $A_{ij}A_{ji}$;
the two coincide when A is symmetric.)
""")
P("P71", "Vieta's formulas", r"""
For a cubic with roots $\lambda^1, \lambda^2, \lambda^3$: $(\lambda - \lambda^1)(\lambda - \lambda^2)(\lambda -
\lambda^3) = \lambda^3 - (\sum\lambda^k)\lambda^2 + (\sum_{k<l}\lambda^k\lambda^l)\lambda - \lambda^1\lambda^2\lambda^3$
— the coefficients are the sum, the pair-sums and the product of the roots. (🔁 Ch. 1 primer P40: sympy `expand`.)
""", code="""
l, a_, b_, c_ = sp.symbols('lambda a b c')                    # the variable and three roots
print(sp.expand((l - a_) * (l - b_) * (l - c_)))              # lambda**3 − (a+b+c) lambda**2 + (ab+ac+bc) lambda − abc
""")
nb.md(r"""
#### The maths, step by step

1. **Contraction:** $\sum_j A_{jj} \equiv A_{jj} = A_{11} + A_{22} + A_{33}$, the trace.
2. **Invariance:** apply (2.12) and set n = m — D18 below shows every closed index chain survives the rotation unchanged.
3. **The three invariants:** $I_1 = A_{ii}$, $I_2 = \tfrac12(I_1^2 - A_{ij}A_{ji})$, $I_3 = \det\mathbf A$; they are the
   coefficients of $\det(\mathbf A - \lambda\boldsymbol\delta) = 0$ (Exercise 2.9), so the eigenvalues of C13 are
   invariants too.
4. **Products** raise the order (N30); contractions of products (2.14) are matrix products in disguise (N31).
""")
D("D18", "Why the trace, I₂ and the determinant do not depend on the axes", ref="Exercise 2.9",
  goal="""Show that three combinations of the nine components of a tensor are the same for every observer, that they
  are the coefficients of the characteristic polynomial, and what they become in the principal frame. The book gives
  hints in Exercise 2.9.""",
  assumptions="A any second-order tensor (steps 1–7); for step 8 A symmetric (so that a principal frame exists).",
  start=(r"A'_{mn} = C_{im}\,C_{jn}\,A_{ij}", "how any second-order tensor's components change under a rotation (2.12, D06)."),
  plan=["Contract (2.12) on m = n: a C-pair collapses to δ.", "Do the same for the chain A_ij A_ji.",
        "Take determinants.", "Expand det(A − λδ) and identify the coefficients; evaluate them in the principal frame."],
  uses=["the transformation rule (2.12) (C06)", r"orthogonality as \"a contracted C-pair is δ\", $C_{im}C_{jm} = \delta_{ij}$ (D02)",
        "Kronecker substitution (N41)", "det(AB) = det A det B (Ch. 1 primer P53)", "Vieta's formulas (P71)",
        "the principal frame (D17 fact 3, forward pointer; step 8 only needs \"τ' = diag(λ)\")"],
  steps=[
      ("Contract (2.12) on its free indices",
       r"A'_{mm} = C_{im}\,C_{jm}\,A_{ij}",
       "Set n = m and sum (contraction, C07). The two C's now share their summed second index m. We contract because a fully contracted quantity has no free index — a candidate scalar.",
       "The trace in the new frame, written with the old components."),
      ("Collapse the C-pair",
       r"A'_{mm} = \delta_{ij}\,A_{ij} = A_{ii}",
       r"$C_{im}C_{jm} = \delta_{ij}$ (D02, second indices summed), then δ substitutes. The trace is the same in both frames: $I_1 = A_{ii}$ is invariant.",
       "Every observer measures the same trace."),
      ("Transform the two-link chain",
       r"A'_{mn}A'_{nm} = C_{im}C_{jn}A_{ij}\;C_{kn}C_{lm}A_{kl}",
       r"Apply (2.12) to each factor with fresh dummy letters (i, j) and (k, l) — reusing a letter would wrongly sum it. We choose $A_{ij}A_{ji}$ because its indices close into a loop.",
       "The chain in the new frame, with four C's."),
      ("Pair the C's and collapse",
       r"A'_{mn}A'_{nm} = (C_{im}C_{lm})(C_{jn}C_{kn})\,A_{ij}A_{kl} = \delta_{il}\delta_{jk}A_{ij}A_{kl} = A_{ij}A_{ji}",
       r"Scalars commute, so group the C's sharing m and those sharing n; each pair is a δ (D02); the δ's substitute l → i and k → j. Hence $I_2 = \tfrac12(I_1^2 - A_{ij}A_{ji})$ is invariant.",
       "Any closed chain of indices survives a rotation unchanged."),
      ("Take the determinant of the matrix form",
       r"\det\mathbf A' = \det(\mathbf C^{\rm T}\mathbf A\,\mathbf C) = (\det\mathbf C)^2\det\mathbf A = \det\mathbf A",
       r"det of a product is the product of dets (P53); det Cᵀ = det C and (det C)² = 1 (D02 step 7). So $I_3 = \det\mathbf A$ is invariant.",
       "The determinant is the third frame-independent number."),
      ("Expand the characteristic determinant in λ",
       r"\det(\mathbf A - \lambda\boldsymbol\delta) = -\lambda^3 + A_{ii}\,\lambda^2 - M\,\lambda + \det\mathbf A",
       r"Expanding a 3×3 determinant whose diagonal entries are $A_{ii} - \lambda$ (cofactors, P53): the λ³ and λ² coefficients come from the diagonal product, the constant term is det A, and M is the sum of the three principal 2×2 minors.",
       "A cubic whose coefficients are built from A alone."),
      ("Identify M with I₂",
       r"M = (A_{11}A_{22} - A_{12}A_{21}) + (A_{22}A_{33} - A_{23}A_{32}) + (A_{11}A_{33} - A_{13}A_{31}) = \tfrac12(A_{ii}A_{jj} - A_{ij}A_{ji})",
       r"Expand $A_{ii}A_{jj} = (A_{11} + A_{22} + A_{33})^2$ and subtract $A_{ij}A_{ji}$ (diagonal squares and the pairs $A_{12}A_{21}$ …): the diagonal squares cancel, each cross product appears twice, hence the ½. So the cubic is $\lambda^3 - I_1\lambda^2 + I_2\lambda - I_3 = 0$ (sign flipped).",
       "The middle coefficient is exactly the second invariant."),
      ("Evaluate in the principal frame",
       r"I_1 = \lambda^1 + \lambda^2 + \lambda^3,\quad I_2 = \lambda^1\lambda^2 + \lambda^2\lambda^3 + \lambda^3\lambda^1,\quad I_3 = \lambda^1\lambda^2\lambda^3",
       "For symmetric A the principal frame has A' = diag(λ) (D17 fact 3); the invariants may be computed in any frame, so use this one: trace, sum of 2×2 minors, determinant of a diagonal matrix — Vieta's formulas for the cubic's roots.",
       "The invariants are the sum, the pair-sums and the product of the eigenvalues — so the eigenvalues are frame-independent too."),
  ],
  result=(r"I_1 = A_{ii},\quad I_2 = \tfrac12(I_1^2 - A_{ij}A_{ji}),\quad I_3 = \det\mathbf A \ \text{ unchanged by (2.12)};\qquad \det(\mathbf A - \lambda\boldsymbol\delta) = -(\lambda^3 - I_1\lambda^2 + I_2\lambda - I_3)",
          "three numbers, and hence the eigenvalues, belong to the tensor, not to the observer; in the principal frame they are the symmetric functions of the λ's."),
  interpret="""The pressure $-\\tau_{ii}/3$ (Ch. 4), the volume strain rate $S_{ii} = \\nabla\\cdot\\mathbf u$ (Ch. 3), and
  the turbulence invariants of Ch. 12 are physical because they are contractions. $A_{ij}A_{ij} = \\mathrm{tr}(\\mathbf A
  \\mathbf A^{\\rm T})$ is also frame-independent, but it is not the coefficient of λ in the cubic — $I_2$ needs the closed
  chain $A_{ij}A_{ji}$ (equal to it only for symmetric A).""",
  check="""Units: A, A², A³ respectively ✓. Limit A = δ: I = (3, 3, 1), polynomial $(1 - \\lambda)^3$ ✓. Number: A =
  [[2, 1, 0], [1, 3, 1], [0, 1, 4]]: $I_1 = 9$, $A_{ij}A_{ji} = 29 + 4 = 33$, $I_2 = ½(81 − 33) = 24$, $I_3 = 18$; roots
  1.268, 3.000, 4.732 sum to 9.000 and multiply to 18.00 ✓ (`invariants`, `characteristic_polynomial` below,
  `principal_axes` in C13).""",
  traps="""reusing a dummy letter in step 3. Using $A_{ij}A_{ij}$ for $I_2$. The sign pattern of the cubic. Reading
  $\\lambda^k$ as a power.""")
note("N30", "Multiplication raises the order", r"""
$P_{ijkl} = A_{ij}B_{kl}$ has four free indices (`ch02.tensor_product(A, B)` = `np.multiply.outer`); it transforms by
(2.13).
""")
note("N31", "The four contractions of (2.14) are matrix products", r"""
| index string | matrix form |
|---|---|
| `A_ij B_ki` | $(\mathbf B\cdot\mathbf A)_{kj}$ |
| `A_ij B_ik` | $(\mathbf A^{\rm T}\cdot\mathbf B)_{jk}$ |
| `A_ij B_kj` | $(\mathbf A\cdot\mathbf B^{\rm T})_{ik}$ |
| `A_ij B_jk` | $(\mathbf A\cdot\mathbf B)_{ik}$ |

Rearrange until the summed index is adjacent, then it is a matrix product; each row is asserted below with
`ch02.contract(A, B, pattern)`.
""", equation=r"A_{ij}B_{ki} = (\mathbf B\cdot\mathbf A)_{kj}, \quad A_{ij}B_{ik} = (\mathbf A^{\rm T}\cdot\mathbf B)_{jk}, \quad A_{ij}B_{kj} = (\mathbf A\cdot\mathbf B^{\rm T})_{ik}, \quad A_{ij}B_{jk} = (\mathbf A\cdot\mathbf B)_{ik}", ref="2.14")
note("N32", "Tensor · vector, two ways", r"""
$A_{ij}u_j = (\mathbf A\cdot\mathbf u)_i$ vs $A_{ij}u_i = (\mathbf A^{\rm T}\cdot\mathbf u)_j$ —
`ch02.dot_tensor_vector(A, u, index=1 or 0)`; exactly the distinction of C05 ($\mathbf n\cdot\boldsymbol\tau$ vs
$\boldsymbol\tau\cdot\mathbf n$).
""")
note("N33", "⚠️ Double dot, two conventions", r"""
the book's $\mathbf A:\mathbf B = A_{ij}B_{ji}$ (transpose pairing) vs the common Frobenius $A_{ij}B_{ij}$ (Wikipedia,
most continuum texts). For A = [[1, 2], [3, 4]], B = [[0, 1], [1, 0]]: book 2·1 + 3·1 = 5, Frobenius 2·1 + 3·1 = 5 too —
a symmetric B hides the difference; take B = [[0, 1], [0, 0]]: book $A_{12}B_{21} + A_{21}B_{12} = 3$, Frobenius
$A_{12}B_{12} = 2$. They agree whenever one operand is symmetric (Ch. 4's dissipation $\tau_{ij}\,\partial u_i/\partial
x_j$). `ch02.double_dot(A, B, convention=…)` — always say which.
""")
nb.worked_example("invariants of A = [[2, 1, 0], [1, 3, 1], [0, 1, 4]]", r"""
1. $I_1 = 2 + 3 + 4 = 9$.
2. $A_{ij}A_{ji}$ = sum of squares (A symmetric) = $4 + 9 + 16 + 2(1 + 0 + 1) = 33$, so $I_2 = \tfrac12(81 - 33) = 24$.
3. $I_3 = 2(3\cdot4 - 1\cdot1) - 1(1\cdot4 - 0) + 0 = 22 - 4 = 18$.
4. Characteristic polynomial $\lambda^3 - 9\lambda^2 + 24\lambda - 18 = 0$; its roots multiply to 18 and add to 9
   (Vieta) — C13 finds them numerically: 1.27, 3.00, 4.73 (sum 9.00 ✓).
""")
nb.code("""
A = np.array([[2., 1., 0.], [1., 3., 1.], [0., 1., 4.]])            # the running symmetric example (EXAMPLE_TENSOR; reused in C13)
print(ch02.trace(A), ch02.invariants(A))                            # 9.0 (9.0, 24.0, 18.0): I₁, I₂, I₃
print(ch02.characteristic_polynomial(A))                            # [1 −9 24 −18]: λ³ − I₁λ² + I₂λ − I₃  (D18)
R = ch02.random_rotation(np.random.default_rng(4))                  # a random rotation
Ap = ch02.transform_tensor(A, R)                                    # the same tensor seen by another observer (2.12)
print(np.round(Ap, 3))                                              # nine different numbers …
print(np.round(ch02.invariants(Ap), 10))                            # … the same three invariants (9, 24, 18)
print(ch02.double_dot(A, Ap, convention="book"), ch02.double_dot(A, Ap, convention="frobenius"))   # equal: both operands symmetric (N33)
B = np.array([[0., 1.], [0., 0.]]); A2 = np.array([[1., 2.], [3., 4.]])   # N33's counter-example
print(ch02.double_dot(A2, B, "book"), ch02.double_dot(A2, B, "frobenius"))   # 3.0 2.0: the conventions differ for non-symmetric operands
""", explain="""
1–2. `invariants` and `characteristic_polynomial` compute $I_1, I_2, I_3$ and the cubic's coefficients (D18 steps 6–7).
3. Rotating the tensor changes every entry but leaves the three invariants untouched (D18 steps 2, 4, 5).
4–5. The two double-dot conventions agree for symmetric operands and differ otherwise (N33).
""")
nb.check_agree("""
I1 = sum(A[i, i] for i in range(3))                                                   # I₁ = A_ii
I2 = 0.5 * (I1**2 - sum(A[i, j] * A[j, i] for i in range(3) for j in range(3)))       # I₂ = ½(I₁² − A_ij A_ji): the closed chain
I3 = (A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])                               # I₃ = det A by cofactors along row 1 (Ch. 1 P53)
      - A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
      + A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]))
assert np.allclose((I1, I2, I3), ch02.invariants(A))                                  # (9, 24, 18)
B3 = np.random.default_rng(7).normal(size=(3, 3))                                     # a random second operand for (2.14)
for pat, mat in [("ij,ki->kj", B3 @ A), ("ij,ik->jk", A.T @ B3), ("ij,kj->ik", A @ B3.T), ("ij,jk->ik", A @ B3)]:   # N31's four rows
    assert np.allclose(ch02.contract(A, B3, pat), mat)                                # each index string IS the stated matrix product
u = np.array([1., 0., 0.])                                                            # a test vector
assert np.allclose(ch02.dot_tensor_vector(A, u, index=1), A @ u) and np.allclose(ch02.dot_tensor_vector(A, u, index=0), A.T @ u)   # N32
print("invariants by hand = library; the four (2.14) patterns are the four matrix products")
""")
nb.figure("""
rng = np.random.default_rng(5)                                       # seeded
n_rot = 50                                                           # observers
ents, invs = [], []                                                  # the nine entries and the three invariants per observer
for k in range(n_rot):                                               # one random rotation per observer
    Ak = ch02.transform_tensor(A, ch02.random_rotation(rng))         # (2.12)
    ents.append(Ak.ravel()); invs.append(ch02.invariants(Ak))        # collect
ents, invs = np.array(ents), np.array(invs)                          # shapes (50, 9), (50, 3)
fig, ax = plt.subplots(figsize=(7, 3.6))                             # one panel
ax.plot(np.repeat(np.arange(n_rot), 9), ents.ravel(), ".", color=COLORS["muted"], ms=4, label="the nine entries A'_ij")   # grey dots
for c, (lab, col) in enumerate([("I₁ (trace)", COLORS["accent"]), ("I₂", COLORS["teal"]), ("I₃ (det)", COLORS["orange"])]):
    ax.plot(invs[:, c], color=col, lw=2, label=lab)                  # three flat lines
ax.set_ylim(-4, 32); ax.legend(fontsize=8, ncol=4, loc="upper center")   # room above I₂ = 24 for the legend
ax.set_xlabel("observer (random rotation number)"); ax.set_ylabel("value [–]")
ax.set_title("50 observers: the entries wander, the three invariants do not")
plt.show()
""", see="Grey dots wander over roughly [−2, 5]; three lines stay perfectly flat at 9, 24 and 18.",
   read="Each dot is one observer's opinion about one component; the lines are facts about the tensor itself (D18).",
   change="…you used a non-symmetric A and also plotted $A_{ij}A_{ij}$: it too stays flat (it is $\\mathrm{tr}(\\mathbf A\\mathbf A^{\\rm T})$, invariant), but it is no longer equal to $A_{ij}A_{ji}$ — the one that belongs in $I_2$.")
nb.md(r"""
**What would change if…** you contracted a *third-order* object? $\varepsilon_{ijk}$ contracted with two vectors is
the cross product — the alternating tensor of the next section.
""")

# =====================================================================================================================
# A.6 §2.7, 2.8 Kronecker Delta, Alternating Tensor, Cross Product — C08 (D09)
# =====================================================================================================================
nb.section("2.7, 2.8", "Kronecker Delta and Alternating Tensor; Vector, Dot, and Cross Products", intro="""
**What is this section about?** Two constant tensors that look the same to every observer: δ (met in C01) and the
alternating tensor ε_ijk. ε turns the cross product, the curl and the rotation tensor into index sums, and one identity
— the epsilon–delta relation — closes every vector identity you will ever need.
""")
core("C08", "The alternating tensor ε_ijk", r"""
How can 27 numbers, all 0 or ±1, encode "right-handed perpendicular" — and give the cross product, the curl and the
Coriolis term as index sums?
""")
nb.md(r"""
#### The problem in plain words

The Coriolis acceleration $2\boldsymbol\Omega\times\mathbf u$, the vorticity $\nabla\times\mathbf u$, the torque
$\mathbf r\times\mathbf F$ — fluid mechanics is full of cross products, and the component formula (2.20) is a mess to
manipulate. $\varepsilon_{ijk}$ is the bookkeeping device that turns "cross" into "sum", and (2.19) is the one identity
that turns products of two ε's back into δ's.
""")
nb.md(r"""
#### The idea

```
ε_ijk = +1  for ijk = 123, 231, 312   (cyclic: walk round 1→2→3→1)
      = −1  for ijk = 321, 213, 132   (anticyclic: walk backwards)
      =  0  whenever two indices agree                 → 6 nonzero entries out of 27
swap two indices → sign flips;  move one index two places (a cyclic shift) → sign kept
```
""")
P("P72", "permutations, cyclic order and parity", r"""
A permutation of (1, 2, 3) is a reordering. It is *even* if it takes an even number of pairwise swaps to reach from 123
(123, 231, 312 — the cyclic shifts) and *odd* otherwise (132, 213, 321). $\varepsilon_{ijk}$ is +1 on even, −1 on odd
permutations, 0 if not a permutation. `itertools.permutations` lists them (sibling of Ch. 1's `combinations`, P55).
""", code="""
from itertools import permutations                         # standard library: every ordering of a tuple
for p in permutations((1, 2, 3)):                          # the six orderings of (1, 2, 3)
    print(p, ch02.permutation_sign(*p))                    # (1,2,3) 1, (1,3,2) −1, (2,1,3) −1, (2,3,1) 1, (3,1,2) 1, (3,2,1) −1
""")
P("P73", "numpy arrays with three axes and np.transpose", r"""
`eps[i, j, k]` indexes a 3×3×3 array (shape (3, 3, 3)); `np.transpose(eps, (1, 2, 0))` reorders the axes so that
`new[i, j, k] = eps[j, k, i]` — the code form of "move an index".
""", code="""
eps = ch02.levi_civita()                                                   # ε_ijk as a (3, 3, 3) array, Python 0-based
print(eps.shape, eps[0, 1, 2], eps[0, 2, 1], eps[0, 0, 1])                 # (3,3,3) 1 −1 0 = ε_123, ε_132, ε_112
print(np.array_equal(np.transpose(eps, (1, 2, 0)), eps))                   # True: a two-place move keeps ε  (np.array_equal: exact equality)
""")
P("P74", "right-hand rule and orientation", r"""
Point the fingers of your right hand along the first vector, curl them toward the second: the thumb gives the direction
of the cross product. $\mathbf e_1\times\mathbf e_2 = \mathbf e_3$ fixes a right-handed frame; the same rule orients
the boundary of a surface in §2.13 (thumb along n, fingers along t).
""", code="""
print(np.cross([1, 0, 0], [0, 1, 0]))    # [0 0 1] = e_3: right-handed
""")
nb.md(r"""
#### The maths, step by step

1. **Definition (2.18)** as in the idea box.
2. **Index moves** (N43): swapping two indices flips the sign, $\varepsilon_{ijk} = -\varepsilon_{ikj}$; moving one
   index two places keeps it, $\varepsilon_{ijk} = \varepsilon_{jki} = \varepsilon_{kij}$.
3. **Cross product (2.21):** $(\mathbf u\times\mathbf v)_k = \varepsilon_{ijk}u_iv_j$ — nine terms, two survive for
   each k; k = 1: $\varepsilon_{231}u_2v_3 + \varepsilon_{321}u_3v_2 = u_2v_3 - u_3v_2$ ✓ (2.20).
4. **Epsilon–delta (2.19):** $\varepsilon_{ijk}\varepsilon_{klm} = \delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}$;
   contracted once $\varepsilon_{pqi}\varepsilon_{pqj} = 2\delta_{ij}$, twice $\varepsilon_{pqr}\varepsilon_{pqr} = 6$ —
   D09.
5. **Isotropy** (N42): δ and ε (up to a factor) are the only constant tensors of order 2 and 3 unchanged by every
   proper rotation.
""")
note("N46", "Cross product, the geometric definition", r"""
$\mathbf w = \mathbf u\times\mathbf v$ has length $uv\sin\theta$, is perpendicular to both, and $(\mathbf u, \mathbf v,
\mathbf w)$ is right-handed (P74); $\mathbf u\times\mathbf v = -\mathbf v\times\mathbf u$; $\mathbf e_1\times\mathbf e_2
= \mathbf e_3$.
""")
note("N47", "Cross product in components (stated)", r"""
because $\mathbf e_i\times\mathbf e_j = \varepsilon_{ijk}\mathbf e_k$ and the product distributes; `ch02.cross(u, v)`
is asserted equal to `np.cross` below. **Determinant form** `N48`: $\mathbf u\times\mathbf v = \det[\mathbf e_1\
\mathbf e_2\ \mathbf e_3;\ u_1\ u_2\ u_3;\ v_1\ v_2\ v_3]$ — one sympy line (`sp.Matrix(u).cross(sp.Matrix(v))`)
reproduces (2.20).
""", equation=r"\mathbf u\times\mathbf v = (u_2v_3 - u_3v_2)\mathbf e_1 + (u_3v_1 - u_1v_3)\mathbf e_2 + (u_1v_2 - u_2v_1)\mathbf e_3", ref="2.20")
note("N49", "Cross product in index form", r"""
`ch02.cross_einsum` is `np.einsum('ijk,i,j->k')`; `ch02.expand_indices_str("eps_ijk u_i v_j")` prints (2.20).
**k = 1 check** `N50`: only (i, j) = (2, 3) and (3, 2) survive — visible in the printed expansion.
""", equation=r"(\mathbf u\times\mathbf v)_k = \sum_{i=1}^{3}\sum_{j=1}^{3}\varepsilon_{ijk}u_iv_j \equiv \varepsilon_{ijk}u_iv_j = \varepsilon_{kij}u_iv_j", ref="2.21")
note("N43", "Index moves on ε", r"""
with `np.transpose` assertions below: axis orders `(1, 2, 0)` and `(2, 0, 1)` return ε (two-place moves); `(0, 2, 1)`
returns −ε (one swap).
""")
note("N42", "Isotropic tensors", r"""
$\delta'_{ij} = C_{im}C_{jn}\delta_{mn} = C_{im}C_{jm} = \delta_{ij}$ (D02) in every frame; ε is isotropic under proper
rotations (det C = +1) and flips sign under a reflection — a *pseudo*tensor, which is why the book says "rotation".
`ch02.is_isotropic(delta)` True, `is_isotropic(eps, proper=True)` True, `is_isotropic(eps, proper=False)` False, a
random tensor False. Pointer: the isotropic Newtonian stress model of Ch. 4 §4.5 is built from δ alone.
""")
note("N44", "The epsilon–delta relation", r"""
the book says "verify by choosing values"; D09 below proves it, and `ch02.epsilon_delta_residual()` checks all 81
cases; the vector triple product $\mathbf a\times(\mathbf b\times\mathbf c) = (\mathbf a\cdot\mathbf c)\mathbf b -
(\mathbf a\cdot\mathbf b)\mathbf c$ is its first payoff.
""", equation=r"\sum_{k=1}^{3}\varepsilon_{ijk}\varepsilon_{klm} \equiv \varepsilon_{ijk}\varepsilon_{klm} = \delta_{il}\delta_{jm} - \delta_{im}\delta_{jl}", ref="2.19")
D("D09", "The epsilon–delta relation, its contractions and the triple product", ref="2.19",
  goal="""Prove the identity that turns a product of two alternating tensors into Kronecker deltas — the tool that
  closes every vector identity (curl of a curl, a × (b × c), the vorticity algebra of Ch. 5). The book says "verify by
  choosing some values".""",
  assumptions="Three dimensions (the count \"exactly one k outside {i, j}\" is 3-D; step 4).",
  start=(r"L_{ijlm} \equiv \varepsilon_{ijk}\,\varepsilon_{klm}\quad(\text{sum on } k)",
         "the left side of (2.19) as a four-index object we will evaluate case by case."),
  plan=["Both sides vanish when i = j or l = m.", "For i ≠ j only one k survives, and both ε's are nonzero only if {l, m} = {i, j}.",
        "Evaluate the two surviving cases: +1 and −1.", "Contract the result once, twice; apply it to a × (b × c)."],
  uses=["definition of ε (2.18) and its index moves (N43)", "antisymmetry: swapping two indices flips the sign",
        "Kronecker substitution and δ_ii = 3 (N41)", "the index form of the cross product (2.21) (N49)", "permutation parity (P72)"],
  steps=[
      ("Fix the free indices and expand the sum",
       r"L_{ijlm} = \varepsilon_{ij1}\varepsilon_{1lm} + \varepsilon_{ij2}\varepsilon_{2lm} + \varepsilon_{ij3}\varepsilon_{3lm}",
       "The summation convention on the repeated k (C01). We write the three terms out so that \"which k survives\" can be argued.",
       "Three products, one per value of the summed index."),
      ("Dispose of i = j",
       r"i = j:\quad L_{iilm} = 0 = \delta_{il}\delta_{im} - \delta_{im}\delta_{il}",
       "ε with two equal indices is zero (2.18), so every term of step 1 vanishes; the right side of (2.19) is a difference of two equal products. Both sides are antisymmetric in (i, j).",
       "Equal first indices: zero on both sides."),
      ("Dispose of l = m the same way",
       r"l = m:\quad L_{ijll} = 0 = \delta_{il}\delta_{jl} - \delta_{il}\delta_{jl}",
       r"$\varepsilon_{kll} = 0$ kills the left side; the right side is again a difference of identical products. From now on i ≠ j and l ≠ m.",
       "Equal last indices: zero on both sides."),
      ("Identify the only surviving k",
       r"i \neq j:\quad \varepsilon_{ijk} \neq 0 \iff k = k^*,\ \{i, j, k^*\} = \{1, 2, 3\}",
       "In 3-D, with i and j distinct, exactly one value of k differs from both; for it ε is ±1, for the other two ε is 0 (2.18). The sum in step 1 has one term.",
       "Only the third index completing the set contributes."),
      ("Require the second ε to be nonzero too",
       r"\varepsilon_{k^*lm} \neq 0 \iff \{l, m\} = \{i, j\}",
       r"l and m must be distinct from $k^*$ and from each other, and the only two values left are i and j. If {l, m} ≠ {i, j}, L = 0 — and then the right side is 0 as well, since l or m equals $k^*$ and every δ there has a zero factor.",
       "Both ε's are alive only when the last pair is the first pair, in some order."),
      ("Evaluate the case l = i, m = j",
       r"L_{ijij} = \varepsilon_{ijk^*}\,\varepsilon_{k^*ij} = \varepsilon_{ijk^*}^2 = 1",
       r"Moving $k^*$ two places ($\varepsilon_{k^*ij} = \varepsilon_{ijk^*}$, N43) makes the two factors equal; a nonzero ε squared is 1. Right side: $\delta_{ii}\delta_{jj} - \delta_{ij}\delta_{ji} = 1\cdot1 - 0 = 1$ (here $\delta_{ii}$ with i fixed is 1, not 3).",
       "Same pair, same order: plus one on both sides."),
      ("Evaluate the case l = j, m = i",
       r"L_{ijji} = \varepsilon_{ijk^*}\,\varepsilon_{k^*ji} = -\varepsilon_{ijk^*}^2 = -1",
       r"$\varepsilon_{k^*ji} = -\varepsilon_{k^*ij}$ (one swap flips the sign), then as in step 6. Right side: $\delta_{ij}\delta_{ji} - \delta_{ii}\delta_{jj} = 0 - 1 = -1$. All 81 cases are now covered: (2.19) holds.",
       "Same pair, swapped order: minus one on both sides."),
      ("Contract (2.19) on j = m",
       r"\varepsilon_{ijk}\,\varepsilon_{klj} = \delta_{il}\delta_{jj} - \delta_{ij}\delta_{jl} = 3\delta_{il} - \delta_{il} = 2\delta_{il}",
       r"Set m = j and sum: now $\delta_{jj} = 3$ (a summed pair, N41) and $\delta_{ij}\delta_{jl} = \delta_{il}$ by substitution. Two-place moves rewrite this as the book's $\varepsilon_{pqi}\varepsilon_{pqj} = 2\delta_{ij}$.",
       "One contraction leaves twice a delta."),
      ("Contract again on i = l",
       r"\varepsilon_{pqr}\,\varepsilon_{pqr} = 2\delta_{rr} = 6",
       r"Set the remaining free pair equal and sum; $\delta_{rr} = 3$. This counts the nonzero entries of ε: six of them, each squared to 1.",
       "Fully contracted, ε with itself gives 6."),
      ("Write a × (b × c) with two ε's",
       r"[\mathbf a\times(\mathbf b\times\mathbf c)]_m = \varepsilon_{mpq}\,a_p\,\varepsilon_{qij}\,b_i c_j",
       r"(2.21) twice: $(\mathbf b\times\mathbf c)_q = \varepsilon_{ijq}b_ic_j = \varepsilon_{qij}b_ic_j$ and $(\mathbf a\times\mathbf d)_m = \varepsilon_{pqm}a_pd_q = \varepsilon_{mpq}a_pd_q$ (two-place moves). The shared summed index q is adjacent in both ε's.",
       "A double cross product is two ε's contracted on one index."),
      ("Apply (2.19) and let the deltas substitute",
       r"[\mathbf a\times(\mathbf b\times\mathbf c)]_m = (\mathbf a\cdot\mathbf c)\,b_m - (\mathbf a\cdot\mathbf b)\,c_m",
       r"$\varepsilon_{mpq}\varepsilon_{qij} = \delta_{mi}\delta_{pj} - \delta_{mj}\delta_{pi}$ by (2.19) with (i, j, k, l, m) → (m, p, q, i, j); then $\delta_{mi}b_i = b_m$, $\delta_{pj}a_pc_j = \mathbf a\cdot\mathbf c$, etc.",
       "The \"BAC − CAB\" rule falls out in one line."),
  ],
  result=(r"\varepsilon_{ijk}\varepsilon_{klm} = \delta_{il}\delta_{jm} - \delta_{im}\delta_{jl};\quad \varepsilon_{pqi}\varepsilon_{pqj} = 2\delta_{ij};\quad \varepsilon_{pqr}\varepsilon_{pqr} = 6;\quad \mathbf a\times(\mathbf b\times\mathbf c) = (\mathbf a\cdot\mathbf c)\mathbf b - (\mathbf a\cdot\mathbf b)\mathbf c",
          "a product of two ε's sharing an index is a difference of two δ-pairs, \"first with first, second with second, minus the crossed pair\"."),
  interpret="""Every identity of vector calculus that contains two crosses — ∇ × (∇ × u) = ∇(∇·u) − ∇²u, the Lagrange
  identity |u × v|² = u²v² − (u·v)², the curl of b × x in Ex. 2.3, the ω ↔ R map of D15 — is this one line plus index
  bookkeeping. It is a 3-D statement (in 2-D the alternating symbol has two indices and $\\varepsilon_{ij}\\varepsilon_{kl}
  = \\delta_{ik}\\delta_{jl} - \\delta_{il}\\delta_{jk}$).""",
  check="""Dimensionless ✓. `epsilon_delta_residual()` = 0 over all 81 cases; `einsum('pqi,pqj')` = 2I,
  `einsum('pqr,pqr')` = 6 ✓ (below). Number: a = (1, 0, 0), b = (0, 1, 0), c = (1, 0, 0): b × c = (0, 0, −1),
  a × (b × c) = (0, 1, 0); right side (a·c)b − (a·b)c = 1·b − 0 = (0, 1, 0) ✓ (`triple_product`).""",
  traps="""$\\delta_{ii} = 3$ when summed, 1 when the index is fixed (steps 6 vs 8). The sign when moving the summed k
  next to its partner. Treating the 81-case table as the proof — it is the check; the proof is the antisymmetry
  argument.""")
nb.worked_example("one case of (2.19) and one cross product", r"""
1. Take i = 1, j = 2, l = 1, m = 2: left side $\sum_k\varepsilon_{12k}\varepsilon_{k12} = \varepsilon_{123}\varepsilon_{312}
   = 1\times1 = 1$ (k = 3 is the only term); right side $\delta_{11}\delta_{22} - \delta_{12}\delta_{21} = 1 - 0 = 1$ ✓.
2. Swap l and m (l = 2, m = 1): left $\varepsilon_{123}\varepsilon_{321} = -1$, right $\delta_{12}\delta_{21} -
   \delta_{11}\delta_{22} = -1$ ✓.
3. Cross product (1, 2, 3) × (4, 5, 6): k = 1: $2\cdot6 - 3\cdot5 = -3$; k = 2: $3\cdot4 - 1\cdot6 = 6$; k = 3:
   $1\cdot5 - 2\cdot4 = -3$ → (−3, 6, −3); check ⊥: $(1, 2, 3)\cdot(-3, 6, -3) = -3 + 12 - 9 = 0$ ✓.
""")
nb.code("""
eps = ch02.levi_civita()                                              # ε_ijk (2.18) as a (3, 3, 3) array
print(eps[0, 1, 2], eps[2, 1, 0], np.count_nonzero(eps))              # 1 −1 6: ε_123, ε_321, and six nonzero entries out of 27
print(ch02.epsilon_delta_residual())                                  # 0.0: all 81 cases of (2.19) hold exactly
print(np.einsum('pqi,pqj->ij', eps, eps)); print(np.einsum('pqr,pqr', eps, eps))   # 2δ_ij (D09 step 8) and 6 (step 9)
u, v = np.array([1., 2., 3.]), np.array([4., 5., 6.])                 # the tiny example's vectors
print(ch02.cross(u, v), ch02.cross_einsum(u, v), np.cross(u, v))      # [−3 6 −3] three times: (2.20), (2.21), numpy
print(ch02.expand_indices_str("eps_ijk u_i v_j"))                     # (2.21) expanded: the three lines of (2.20) — two terms survive per k (N50)
a_, b_, c_ = np.random.default_rng(8).normal(size=(3, 3))             # three random vectors
print(np.allclose(ch02.triple_product(a_, b_, c_), (a_ @ c_) * b_ - (a_ @ b_) * c_))   # True: a × (b × c) = (a·c)b − (a·b)c (D09 step 11)
rng = np.random.default_rng(9)                                        # for the isotropy tests
print(ch02.is_isotropic(np.eye(3), rng)[0], ch02.is_isotropic(eps, rng, proper=True)[0],
      ch02.is_isotropic(eps, rng, proper=False)[0], ch02.is_isotropic(A, rng)[0])   # True True False False (N42)
""", explain="""
1. Six nonzero entries; Python's `eps[0, 1, 2]` is the book's $\\varepsilon_{123}$.
2. All 81 cases of (2.19) hold exactly (D09's check).
3. The two contractions: $2\\delta_{ij}$ and 6.
4. Three cross products agree — the component formula, the einsum and numpy's.
5. The expansion shows the k = 1 check of the worked example.
6. The triple-product identity for random vectors.
7. Isotropy: δ always, ε only under proper rotations (a pseudotensor), a generic tensor never.
""")
nb.check_agree("""
from itertools import permutations, product                           # permutations (P72); product = nested loops in one line
eps_mine = np.zeros((3, 3, 3))                                        # from scratch: fill the 6 nonzero entries
for p in permutations(range(3)):                                      # the six orderings of (0, 1, 2)
    eps_mine[p] = ch02.permutation_sign(*[q + 1 for q in p])          # +1 even, −1 odd (book indices 1–3)
assert np.array_equal(eps_mine, eps)                                  # identical to levi_civita()
for i, j, l, m in product(range(3), repeat=4):                        # all 81 (i, j, l, m) cases of (2.19)
    lhs = sum(eps[i, j, k] * eps[k, l, m] for k in range(3))          # the sum over k
    rhs = (i == l) * (j == m) - (i == m) * (j == l)                   # δ_il δ_jm − δ_im δ_jl
    assert lhs == rhs                                                 # exact integers
print("ε from parity = levi_civita(); (2.19) holds in all 81 cases")
""")
nb.figure("""
from matplotlib.colors import ListedColormap                          # a colormap from a list of colours
fig, axs = plt.subplots(1, 3, figsize=(9.5, 3.4), sharey=True)         # one slice per k
labels = {1: "+1", -1: "−1"}                                          # cell labels
for k, ax in enumerate(axs):                                          # slice k (Python) = ε_ij(k+1)
    sl = eps[:, :, k]                                                 # the 3×3 slice
    ax.imshow(sl, cmap=ListedColormap([COLORS["teal"], "white", COLORS["accent"]]), vmin=-1, vmax=1)   # −1 teal, 0 white, +1 purple
    ax.grid(False)                                                    # no grid lines through the cells
    for i in range(3):                                                # label the nonzero cells with their index triple
        for j in range(3):
            if sl[i, j] != 0:
                ax.text(j, i, f"ε{i + 1}{j + 1}{k + 1}\\n{labels[int(sl[i, j])]}", ha="center", va="center", fontsize=9, color="white", fontweight="bold")   # two-line label: index triple, then the value
    ax.set_xticks(range(3), ["j=1", "j=2", "j=3"]); ax.set_yticks(range(3), ["i=1", "i=2", "i=3"]); ax.set_title(f"slice k = {k + 1}")
fig.suptitle("ε_ijk unfolded: six coloured cells in 27 — an antisymmetric 2-D pattern in every slice")
plt.show()
""", see="Three 3×3 slices, six coloured cells in all: purple +1, teal −1, white 0; the coloured pair shifts by one place from slice to slice.",
   read="Each slice k holds the 2-D antisymmetric pattern of the remaining two indices — the matrix (2.26) of C12 in disguise: $-\\varepsilon_{ijk}\\omega_k$ fills exactly these cells.",
   change="…you relabelled 1 ↔ 2 everywhere (a mirror): every colour flips — ε changes sign under a reflection (N42).")
nb.md(r"""
**What would change if…** the second vector were the operator $\partial/\partial x_j$? $\varepsilon_{ijk}\,\partial
u_k/\partial x_j$ — the curl (C11). And $-\varepsilon_{ijk}\omega_k$ packs a vector into an antisymmetric matrix
(C12). ε is the bridge between vectors and antisymmetric tensors.
""")

# =====================================================================================================================
# A.7 §2.9 Gradient, Divergence, and Curl — R02, C09, C10, C11 (D12, E5)
# =====================================================================================================================
nb.section("2.9", "Gradient, Divergence, and Curl", intro="""
**What is this section about?** The del operator acting three ways: on a scalar it gives the direction of steepest
climb (gradient), on a vector either how much the field spreads (divergence) or how much it spins (curl). We compute
all three on a grid you will meet again in every chapter, and check the stencils converge at second order.
""")
nb.recap("R02", "The del operator ∇", r"""
Ch. 1 primed the partial derivative (P25: the slope in one variable with the others held fixed) and ∇ as the vector
that "points uphill", and finite differences as its numerical form (P21: the central difference is second order; P22:
`np.gradient`'s one-sided edges). New here: the index form $\nabla = \mathbf e_i\,\partial/\partial x_i$ (2.22) — an
operator with one free index, so it raises the order of whatever it acts on by one, or lowers it by one when
contracted.
""", where="Ch. 1 §1.5 primers P21, P22, P25")
core("C09", "The gradient ∇φ: perpendicular to the level sets, along the steepest climb", r"""
On a contour map of a scalar φ, which way is steepest, how steep is it, and how steep is it in some *other* direction
n?
""")
nb.md(r"""
#### The problem in plain words

Pressure maps drive the wind: air accelerates from high to low pressure, fastest where the isobars are closest — down
$-\nabla p$. A temperature map drives heat flux down $-\nabla T$ (Ch. 1, Fourier). The gradient turns a scalar map into
a vector field: direction of fastest increase, magnitude the rate.
""")
nb.md(r"""
#### The idea

```
contours of φ (level sets)  ──  ∇φ is perpendicular to them  (moving along a contour changes nothing)
|∇φ| = the slope in that steepest direction;  in any other unit direction n:  ∂φ/∂n = ∇φ · n = |∇φ| cos(angle)
```
""")
P("P75", "level sets and the directional derivative", r"""
A level set (contour line in 2-D, surface in 3-D) is where φ takes one value. The directional derivative ∂φ/∂n is the
rate of change of φ per metre when you walk along the unit vector n; by the chain rule (🔁 Ch. 1 P49) it equals ∇φ·n,
so it is zero along a contour (n ⊥ ∇φ) and largest along ∇φ.
""", code="""
g_ = np.array([2 * 1.0, 2 * 2.0])                                  # ∇φ of φ = x² + y² at (1, 2) = (2, 4)
for n_ in ([1, 0], [0, 1], g_ / np.linalg.norm(g_), np.array([-2, 1]) / np.sqrt(5)):   # four unit directions
    print(np.round(g_ @ np.asarray(n_), 3))                        # 2.0 4.0 4.472 0.0 — the last one walks along the contour
""")
P("P76", "np.meshgrid and the project grid layout", r"""
`np.meshgrid` turns 1-D coordinate arrays into arrays that hold the coordinate of every grid point. **Project convention
(every later chapter):** 3-D arrays are indexed `[k, j, i]` = (z, y, x) — x on the last axis — built with
`np.meshgrid(z, y, x, indexing='ij')`; 2-D arrays are `[j, i]` = (y, x) from `indexing='xy'`. Vector fields stack the
component on axis 0 (`np.stack`): `u[c, k, j, i]`. `ch02.grid` and `ch02.grid2d` do this once so nobody transposes a
curl by accident.
""", code="""
g2 = ch02.grid2d(((-1, 1), (-1, 1)), 5)                    # a 5×5 grid on [−1, 1]²
print(g2.X.shape, g2.X[0], g2.Y[:, 0])                     # (5, 5); x varies along the LAST axis, y along the first
g3 = ch02.grid(((0, 1), (0, 2), (0, 3)), 4)                # a 4×4×4 grid with different extents
print(g3.X.shape, g3.h)                                    # (4, 4, 4) laid out as (z, y, x); spacings (hx, hy, hz) = (0.333, 0.667, 1.0)
""")
P("P77", "numpy broadcasting", r"""
An operation between a full grid array and a scalar (or a lower-dimensional array) is applied at every grid point at
once; `X**2 + Y**2` on meshgrid arrays evaluates φ everywhere without a loop. Slices like `phi[:, 2:] - phi[:, :-2]`
subtract neighbours along one axis for the whole grid.
""", code="""
X, Y = g2.X, g2.Y                                          # the 5×5 coordinate arrays
phi = X**2 + Y**2                                          # φ at all 25 points, no loop
print((phi[:, 2:] - phi[:, :-2]).shape)                    # (5, 3): the differences of x-neighbours two apart (interior only)
""")
P("P78", "plt.contour, plt.quiver and plt.streamplot", r"""
`ax.contour(X, Y, phi, levels)` draws level curves; `ax.quiver(X, Y, U, V)` draws an arrow at each point;
`ax.streamplot(x, y, U, V)` draws curves tangent to a vector field (1-D x, y arrays; the field on an `[j, i]` grid —
our layout).
""", code="""
fig, ax = plt.subplots(figsize=(3, 3))                     # a small square panel
ax.contour(X, Y, phi, levels=6)                            # six level curves of φ = x² + y² (circles, coarse on a 5×5 grid)
ax.quiver(X, Y, 2 * X, 2 * Y, scale=20)                    # the gradient (2x, 2y) as arrows
ax.set_aspect('equal'); plt.show()                         # arrows cross the circles at right angles
""")
nb.md(r"""
#### The maths, step by step

1. (2.22) $\nabla = \mathbf e_i\,\partial/\partial x_i$; on a scalar, $(\nabla\phi)_i = \partial\phi/\partial x_i$ — a
   vector (it passes (2.8) by the chain rule; the residual table of C03 said so).
2. Along a contour φ is constant, so its rate of change along the tangent t is $\nabla\phi\cdot\mathbf t = 0$: ∇φ ⊥
   level sets.
3. In direction n: $\partial\phi/\partial n = \nabla\phi\cdot\mathbf n = |\nabla\phi|\cos\alpha$ — maximal
   ($|\nabla\phi|$) for n ∥ ∇φ (the book's Fig. 2.7).
4. On the grid: the central difference $(\phi_{i+1} - \phi_{i-1})/2h$ with error ∝ h² (P21); at the edges a one-sided
   second-order stencil $(-3\phi_i + 4\phi_{i+1} - \phi_{i+2})/2h$ keeps the order.
""")
nb.worked_example("φ = x² + y² at (1, 2)", r"""
1. $\partial\phi/\partial x = 2x = 2$, $\partial\phi/\partial y = 2y = 4$: ∇φ = (2, 4), $|\nabla\phi| = \sqrt{20} = 4.472$.
2. Along $\mathbf n = (1, 0)$: $\partial\phi/\partial n = 2$.
3. Along the contour direction $(-2, 1)/\sqrt5$: $(-4 + 4)/\sqrt5 = 0$.
4. Along ∇φ itself: 4.472 — the steepest.
5. Central difference with h = 0.1 at x = 1: $(1.1^2 - 0.9^2)/0.2 = (1.21 - 0.81)/0.2 = 2.000$ — exact, because a
   quadratic has no third derivative.
""")
nb.md(r"""
*Gloss — observed order of convergence:* the slope of log(error) vs log(h) (🔁 Ch. 1 P13 read slopes on log–log axes);
`tools.convergence.observed_order(h, err)` fits it by least squares. A second-order stencil has slope 2: halving h
quarters the error.
""")
nb.code("""
from tools.convergence import observed_order                       # least-squares slope of log(err) vs log(h) (gloss above)
g = ch02.grid2d(((-2, 2), (-2, 2)), 41)                            # 41×41 nodes on [−2, 2]²: h = 0.1 (P76)
phi = g.X**2 + g.Y**2                                              # φ at every node by broadcasting (P77)
grad = ch02.gradient(phi, g.h)                                     # (2.22) on a scalar: (∇φ)_i = ∂φ/∂x_i, components stacked on axis 0 as (x, y)
print(grad.shape, grad[:, 30, 30], (g.x[30], g.y[30]))            # (2, 41, 41); ∇φ at node (i, j) = (30, 30) = the point (1, 1): [2. 2.] — exact for a quadratic
print(ch02.directional_derivative(grad, [1, 0])[30, 30],           # ∂φ/∂n along e_1 at (1, 1): 2.0
      ch02.directional_derivative(grad, np.array([-1, 1]) / np.sqrt(2))[30, 30])   # along the contour direction at (1, 1): 0.0
ns = (8, 16, 32) if FAST else (8, 16, 32, 64)                      # grid sizes for the convergence study (FAST: three levels)
study = {op: ch02.operator_convergence(op, ns) for op in ("gradient", "divergence", "curl")}   # max error vs h on a sin/cos field (cached; reused in C15, C16)
print({op: round(s["order"], 2) for op, s in study.items()})       # observed orders ≈ 2.0 for all three
""", explain="""
1. A 41×41 grid on [−2, 2]² (spacing h = 0.1) and φ at every node.
2. `gradient` applies `partial` along each coordinate direction; components are stacked on axis 0 as (x, y).
3. At the node (1, 1) the gradient is (2, 2) exactly (central differences are exact for quadratics).
4. `directional_derivative` is ∇φ·n: 2 along $\\mathbf e_1$, 0 along the contour.
5. `operator_convergence` measures the max error of gradient, divergence and curl on a sin/cos field against the exact
   sympy derivatives for n = 8 … 64 nodes; the log–log slope is 2 — the stencils are second order, edges included.
""")
nb.check_agree("""
h = g.h[0]                                                                          # spacing in x
dphidx = np.zeros_like(phi)                                                         # ∂φ/∂x from scratch
dphidx[:, 1:-1] = (phi[:, 2:] - phi[:, :-2]) / (2 * h)                              # central difference inside (P21), along the LAST axis = x
dphidx[:, 0] = (-3 * phi[:, 0] + 4 * phi[:, 1] - phi[:, 2]) / (2 * h)               # second-order one-sided stencil at the left edge
dphidx[:, -1] = (3 * phi[:, -1] - 4 * phi[:, -2] + phi[:, -3]) / (2 * h)            # … and at the right edge
assert np.allclose(dphidx, ch02.partial(phi, 0, h))                                 # direction 0 = x; the library maps it to array axis 1 of the [j, i] layout
assert np.allclose(dphidx, grad[0])                                                 # and it is the first component of `gradient`
print("x is the LAST array axis in the [j, i] layout — `partial` takes the coordinate direction, not the axis")
""")
nb.figure("""
phi7 = g.X**2 + g.Y**2 / 4                                                   # Fig. 2.7's scalar: ellipses as level sets
grad7 = ch02.gradient(phi7, g.h)                                             # its gradient on the grid
probe = np.array([1.0, 1.0]); gp = np.array([2 * probe[0], probe[1] / 2])    # the probe point and ∇φ there = (2, 0.5)
n30 = np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30))])             # a direction n at 30°
dphidn = float(gp @ n30)                                                     # ∂φ/∂n = ∇φ·n = 2·0.866 + 0.5·0.5 = 1.98
fig, ax = plt.subplots(figsize=(6.2, 5.2))                                   # one square panel
cs = ax.contour(g.X, g.Y, phi7, levels=12, colors=COLORS["muted"], linewidths=1); ax.clabel(cs, fmt="%.1f", fontsize=7)   # level sets (P78)
sk = 4; ax.quiver(g.X[::sk, ::sk], g.Y[::sk, ::sk], grad7[0][::sk, ::sk], grad7[1][::sk, ::sk], color=COLORS["accent"], scale=28, width=0.004)   # ∇φ arrows
ax.annotate("", probe + 0.35 * gp, probe, arrowprops=dict(arrowstyle="->", color=COLORS["accent"], lw=2.5))       # bold ∇φ at the probe
ax.annotate("", probe + 0.7 * n30, probe, arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=2.5))       # the direction n
box = dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9)                                             # white boxes keep the labels readable over the contours
ax.text(1.2, 0.72, "∇φ = (2, 0.5)", color=COLORS["accent"], bbox=box); ax.text(1.15, 1.75, "n at 30°", color=COLORS["orange"], bbox=box)
ax.text(-1.9, -1.85, f"∂φ/∂n = ∇φ·n = 2·0.866 + 0.5·0.5 = {dphidn:.2f}", fontsize=9, bbox=box)                    # the number
ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_title("∇φ crosses every level set at right angles and lengthens where they crowd")
plt.show()
""", see="Elliptical contours of φ = x² + y²/4; purple arrows crossing every contour at right angles, longer where the contours crowd (left and right); at the probe (1, 1) a bold ∇φ and an orange n at 30° with ∂φ/∂n = 1.98.",
   read="|∇φ| is the local contour density; n's projection on ∇φ is the slope you would feel walking along n — 1.98 here, against 2.06 = |∇φ| straight uphill. `N51`: this is our own Fig. 2.7.",
   change="…you rotated n to lie along a contour (perpendicular to ∇φ, about 104° here): ∂φ/∂n → 0.")
nb.plotly("""
ang = np.linspace(0, 2 * np.pi, 121)                                              # to draw the ellipses
def probe_dir(alpha_deg):                                                        # one slider position: n at angle α
    n_ = np.array([np.cos(np.deg2rad(alpha_deg)), np.sin(np.deg2rad(alpha_deg))])   # unit direction
    d = float(gp @ n_)                                                           # ∂φ/∂n = ∇φ·n at the probe
    out = {f"level set φ = {c}": (np.sqrt(c) * np.cos(ang), 2 * np.sqrt(c) * np.sin(ang)) for c in (0.5, 1.25, 2.5)}   # three fixed ellipses of x² + y²/4 = c
    out["∇φ at (1, 1)"] = ([1, 1 + 0.4 * gp[0]], [1, 1 + 0.4 * gp[1]])           # the gradient arrow (fixed)
    out["n (drag α)"] = ([1, 1 + 0.8 * n_[0]], [1, 1 + 0.8 * n_[1]])              # the turning direction
    out["∂φ/∂n · n  (the slope you feel)"] = ([1, 1 + 0.4 * d * n_[0]], [1, 1 + 0.4 * d * n_[1]])   # ∇φ·n drawn along n: positive uphill, negative downhill
    return out
fig = slider_figure(probe_dir, "α", np.linspace(0, 350, 36 if not FAST else 18), unit="°", xlabel="x [m]", ylabel="y [m]",
                    title="∂φ/∂n = ∇φ·n: largest along ∇φ, zero along the contour, negative downhill",
                    xrange=[-2.2, 3.4], yrange=[-3.4, 3.4], height=520, active=1)
fig.update_yaxes(scaleanchor="x", scaleratio=1)                                  # equal axes
fig.show()
""", explain="""
1. Three level sets of φ = x² + y²/4 are drawn as fixed ellipses; the gradient at the probe (1, 1) is (2, 0.5).
2. Dragging α turns the direction n; the third arrow has length ∇φ·n along n — it is the directional derivative made
   visible (positive uphill, shrinking to nothing along the contour, negative downhill).
""")
see_read_change("Three ellipses, a fixed purple gradient arrow at (1, 1), and an orange direction n that turns with the slider; a third arrow along n whose length is ∂φ/∂n.",
                "The third arrow is longest (2.06) when n ∥ ∇φ (α ≈ 14°), vanishes along the contour (α ≈ 104°) and reverses downhill (α ≈ 194°).",
                "…the probe sat at (0, 1): ∇φ = (0, 0.5) — a shorter arrow pointing straight up; the ellipses are flattest there.")
nb.figure("""
fig, ax = plt.subplots(figsize=(6, 3.8))                                        # one panel, log–log
for op, col in zip(("gradient", "divergence", "curl"), (COLORS["accent"], COLORS["teal"], COLORS["orange"])):   # the three cached studies
    s = study[op]                                                                # dict(n, h, err, order)
    ax.loglog(s["h"], s["err"], "o-", color=col, label=f"{op}: observed order {s['order']:.2f}")   # error vs spacing
h_ref = np.array([study["gradient"]["h"].min(), study["gradient"]["h"].max()])   # for the reference line
ax.loglog(h_ref, study["gradient"]["err"].max() * (h_ref / h_ref.max())**2, "--", color=COLORS["muted"], label="slope 2 reference")   # ∝ h²
ax.set_xlabel("grid spacing h [m]"); ax.set_ylabel("max |numerical − exact| [–]"); ax.legend(fontsize=8)
ax.set_title("Halving h quarters the error: second order, edges included")
plt.show()
""", see="Three straight lines of slope 2 on log–log axes, parallel to the dashed reference.",
   read="Halving h quarters the error — the stencil is second order also at the edges, thanks to the one-sided second-order stencil (the from-scratch cell above).",
   change="…you replaced the edge stencil by a first-order one (the `np.gradient` default of Ch. 1 P22): the edge error would fall only like h and the slope would drop towards 1.")
nb.md(r"""
**What would change if…** ∇ acted on a vector instead? Two choices: contract the operator's index with the vector's
(divergence, C10) or cross them with ε (curl, C11).
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C10", "The divergence ∇·u = ∂u_i/∂x_i", r"""
How much does a velocity field *spread out* at a point — and why is that one number the most used quantity in the
book?
""")
nb.md(r"""
#### The problem in plain words

Mass conservation (Ch. 4): if more fluid leaves a tiny box than enters, the density inside must drop. "More leaves than
enters" per unit volume is the divergence. For water and slow air ∇·u = 0 — the incompressibility condition that shapes
every chapter after this one.
""")
nb.md(r"""
#### The idea

```
∂u_1/∂x_1 : does u_1 grow along x_1?  (stretching along 1)        ∇·u = sum of the three stretching rates
a box with more leaving than entering  →  ∇·u > 0  (source);  the same in and out  →  ∇·u = 0  (solenoidal)
```
""")
nb.md(r"""
#### The maths, step by step

1. Contract ∇'s index with u's: $\nabla\cdot\mathbf u = \partial u_i/\partial x_i = \partial u_1/\partial x_1 +
   \partial u_2/\partial x_2 + \partial u_3/\partial x_3$ (2.23) — a scalar.
2. It is the trace of the velocity gradient $G_{ij} = \partial u_i/\partial x_j$ (N52): $\nabla\cdot\mathbf u =
   G_{ii}$, hence an invariant (C07).
3. Ex. 2.3: $\mathbf u = a\mathbf x$: $\partial(ax_i)/\partial x_i = a\,\delta_{ii} = 3a$. *Gloss:* $\partial x_i/
   \partial x_j = \delta_{ij}$ — each coordinate depends only on itself. The book's third term reads $ax_3\mathbf e_3$
   (its "$ax_2\mathbf e_2$" printed twice is a misprint). $\mathbf u = \mathbf b\times\mathbf x$: $\partial(\varepsilon_
   {lmk}b_lx_m)/\partial x_k = \varepsilon_{lmk}b_l\delta_{mk} = \varepsilon_{lkk}b_l = 0$.
""")
note("N52", "⚠️ Two index-order conventions to fix once", r"""
The *gradient of a vector* raises the order: `ch02.vector_gradient(u, h)[i, j] = ∂u_i/∂x_j` (Kundu's Ch. 3 usage,
$S_{ij} = \tfrac12(\partial u_i/\partial x_j + \partial u_j/\partial x_i)$). The *divergence of a tensor* lowers it and
contracts the **second** index, $(\nabla\cdot\boldsymbol\tau)_i = \partial\tau_{ij}/\partial x_j$ — the opposite index
from Cauchy's $f_i = \tau_{ji}n_j$; the two agree only for symmetric τ (Ch. 4's momentum equation uses
$\partial\tau_{ij}/\partial x_j$). `ch02.tensor_divergence(T, h, index=1)` exposes the choice. (Comma notation for
both: §2.14.)
""", equation=r"(\nabla\cdot\boldsymbol\tau)_i = \sum_{j=1}^{3} \frac{\partial\tau_{ij}}{\partial x_j} \equiv \frac{\partial\tau_{ij}}{\partial x_j}", ref="§2.9")
nb.worked_example("Ex. 2.3 with a = 1: u = x", r"""
1. $\mathbf u = (x_1, x_2, x_3)$.
2. $\partial u_1/\partial x_1 = 1$, likewise for 2 and 3.
3. $\nabla\cdot\mathbf u = 1 + 1 + 1 = 3 = 3a$ ✓ — every direction stretches at unit rate.
4. On the grid with h ≈ 0.04: central differences of a linear field are exact, so the code prints 3.000 at every point.
""")
nb.code("""
g3 = ch02.grid(((-1, 1),) * 3, 48 if not FAST else 24)                # a 3-D grid on [−1, 1]³, layout [k, j, i] (P76); 48³ nodes (FAST: 24³)
u_rad = ch02.radial_field(1.0)(g3.X, g3.Y, g3.Z)                      # Ex. 2.3's u = a x with a = 1, evaluated by broadcasting: shape (3, n, n, n)
div = ch02.divergence(u_rad, g3.h)                                    # (2.23): ∂u_i/∂x_i on the grid
print(div.min(), div.max())                                           # 3.0 3.0 (to round-off): 3a everywhere
u_rot = ch02.solid_body_rotation_field([0, 0, 1.0])(g3.X, g3.Y, g3.Z) # Ex. 2.3's u = b × x with b = e_3
print(np.abs(ch02.divergence(u_rot, g3.h)).max(), ch02.is_solenoidal(u_rot, g3.h))   # ~1e-15, (True, residual): solenoidal
G = ch02.vector_gradient(u_rad, g3.h)                                 # G[i, j] = ∂u_i/∂x_j at every node (N52)
print(G[:, :, 5, 5, 5])                                               # the identity matrix: ∂x_i/∂x_j = δ_ij
print(np.allclose(np.einsum('ii...', G), div))                        # True: the divergence is the trace of G (C07)
X3 = sp.Matrix(sp.symbols('x1 x2 x3'))                                # sympy coordinates (Ch. 1 P40, P61)
print(ch02.exact_div_curl(sp.Symbol('a') * X3))                       # the symbolic twin: (3a, [0, 0, 0])
""", explain="""
1. A 3-D grid (kept for C11 and C12).
2. `radial_field(1.0)` is a callable field with sympy behind it; evaluated on the grid arrays it returns the three
   components stacked on axis 0.
3. Its divergence is 3 everywhere — exact, because the field is linear.
4. The rotation field is solenoidal: its divergence is zero to round-off.
5. The velocity gradient of u = x is δ_ij at every node; its trace is the divergence.
6. `exact_div_curl` computes the same symbolically: (3a, 0).
""")
nb.check_agree("""
div_mine = (ch02.partial(u_rad[0], 0, g3.h[0])                        # ∂u_1/∂x_1: component 0 differentiated along coordinate direction 0 (x)
            + ch02.partial(u_rad[1], 1, g3.h[1])                      # ∂u_2/∂x_2
            + ch02.partial(u_rad[2], 2, g3.h[2]))                     # ∂u_3/∂x_3
assert np.allclose(div_mine, div)                                     # (2.23) written out = the library's divergence
print("component c is differentiated along coordinate direction c; inside `partial`, direction 0 (x) is array axis 2 of [k, j, i]")
""")
nb.figure("""
k0 = g3.X.shape[0] // 2                                               # the slice index closest to z = 0
sk = 3 if not FAST else 2                                             # quiver thinning
fig, axs = plt.subplots(1, 2, figsize=(9, 4.2))                       # left: u = a x, right: u = b × x
for ax, U, ttl in zip(axs, (u_rad, u_rot), ("u = a x:  ∇·u = 3a everywhere (a source at every point)", "u = b × x:  ∇·u = 0 (spins, does not spread)")):
    dv = ch02.divergence(U, g3.h)[k0]                                 # the z = 0 slice of ∇·u
    im = ax.imshow(dv, origin="lower", extent=(-1, 1, -1, 1), cmap="Oranges", vmin=0, vmax=3.2)   # heatmap: 0 white … 3 orange
    ax.quiver(g3.X[k0, ::sk, ::sk], g3.Y[k0, ::sk, ::sk], U[0][k0, ::sk, ::sk], U[1][k0, ::sk, ::sk], color=COLORS["ink"], scale=25, width=0.003)   # the field's arrows
    ax.set_title(ttl, fontsize=9); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_aspect("equal")
fig.colorbar(im, ax=axs, shrink=0.8, label="∇·u [1/s]")
plt.show()
""", see="Left: spreading arrows and a uniform orange square (∇·u = 3). Right: circling arrows and a blank white square (∇·u = 0).",
   read="Divergence measures spreading, not speed — the rotating field is fast at the rim yet has zero divergence everywhere.",
   change="…u = (x, −y, 0): stretch in x, squeeze in y, divergence 1 − 1 = 0 — solenoidal without rotating (a pure strain).")
nb.md(r"""
*Forward pointer:* the exact tracer positions in the next animation are $\mathbf x(t) = e^{\mathbf Gt}\mathbf x_0$
from `ch02.linear_flow_map(G, t)` — the matrix exponential is primed in C12 (P79); here it is only a way to move the
dots exactly.
""")
nb.animation("""
from fluidpy.core.anim import animate                                     # (Ch. 1 P16)
n_fr = 40 if not FAST else 24                                             # frames (FAST: 24)
ts = np.linspace(0, 0.6, n_fr)                                            # time [s]
ang = np.linspace(0, 2 * np.pi, 40, endpoint=False)                       # a ring of tracers …
ring0 = 0.5 * np.stack([np.cos(ang), np.sin(ang)])                        # … of radius 0.5 m
G_src = 1.0 * np.eye(2)                                                   # velocity gradient of u = a x (a = 1/s): G = a I
G_rot = np.array([[0., -1.], [1., 0.]])                                   # velocity gradient of u = b × x (b = 1/s): G = [b×]
fig, axs = plt.subplots(1, 2, figsize=(8, 3.9))                           # left: source, right: rotation
fig.set_layout_engine("none")                                             # fixed layout: the automatic layout engine would re-run on every frame
dots = []                                                                 # the scatter artists to move
for ax, ttl in zip(axs, ("u = a x:  ∇·u = 3a, ∇×u = 0", "u = b × x:  ∇·u = 0, ∇×u = 2b")):
    ax.plot(*ring0, ".", color=COLORS["muted"], ms=3)                     # the ring at t = 0 (grey ghost)
    dots.append(ax.scatter(*ring0, s=14, color=COLORS["accent"]))         # the moving ring
    ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1); ax.set_aspect("equal"); ax.set_title(ttl, fontsize=9); ax.set_xlabel("x [m]")
    ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])                  # few tick labels: text is what makes video frames slow to draw
txt = fig.text(0.5, 0.01, "", ha="center")                                # the clock
def update(i):                                                            # frame i
    for d, G_ in zip(dots, (G_src, G_rot)):                               # each panel
        d.set_offsets((ch02.linear_flow_map(G_, ts[i]) @ ring0).T)        # exact tracer positions x(t) = e^{Gt} x_0
    txt.set_text(f"t = {ts[i]:.2f} s")
    return (*dots, txt)
show_animation(animate(update, frames=n_fr, fig=fig, interval=60), player="video")   # smooth MP4
""", explain="""
1. Forty tracer dots on a ring of radius 0.5 m; each frame moves them exactly with the linear flow's map.
2. Left: u = a x — the ring expands (∇·u = 3a; in this plane the ring's area grows at rate 2a per unit area). Right:
   u = b × x — the ring turns rigidly and keeps its size (∇·u = 0); its spin is half the curl 2b (C11).
""")
see_read_change("One ring grows; the other spins without changing size.",
                "Divergence is the rate at which the ring's area grows (2a per unit area in this plane, 3a per unit volume in 3-D); curl is twice the ring's spin rate.",
                "…a < 0: the ring shrinks — a sink; a mixed G = a I + [b×] would do both at once.")
nb.md(r"""
**What would change if…** you crossed ∇ with u instead of dotting? The curl — spin instead of spread.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C11", "The curl ∇×u = ε_ijk ∂u_k/∂x_j", r"""
How much does a velocity field *spin* at a point — and why does a perfectly straight shear flow have curl?
""")
nb.md(r"""
#### The problem in plain words

Vorticity — the curl of the velocity — is the subject of Ch. 5 and half of Ch. 13: hurricanes, ocean gyres and the
planetary vorticity $2\boldsymbol\Omega$ are all $\nabla\times\mathbf u$. A paddle wheel dropped into the flow turns at
half the curl. It turns in a whirlpool — and also in a straight river whose speed varies across the channel.
""")
nb.md(r"""
#### The idea

```
(∇×u)_3 = ∂u_2/∂x_1 − ∂u_1/∂x_2 :  does u_2 grow to the right?  (+)   does u_1 grow upward?  (−)
a paddle wheel spins if the flow on one side is faster than on the other   →  shear flow has curl
ε_ijk ∂_j u_k = "cross ∇ with u" with the operator on the left: ∂ acts on u_k, never the other way round
```
""")
nb.md(r"""
#### The maths, step by step

1. (2.21) with $u_i \to \partial/\partial x_j$: $(\nabla\times\mathbf u)_i = \varepsilon_{ijk}\,\partial u_k/\partial x_j$
   (2.24), operator ordering: ∂ acts on u.
2. D12 below expands it into the three components (2.25).
3. Ex. 2.3: $\mathbf u = a\mathbf x$: $\varepsilon_{ijk}\partial(ax_k)/\partial x_j = a\varepsilon_{ijk}\delta_{jk} =
   a\varepsilon_{ijj} = 0$ — irrotational. $\mathbf u = \mathbf b\times\mathbf x$: $\varepsilon_{ijk}\partial_j(
   \varepsilon_{lmk}b_lx_m) = \varepsilon_{ijk}\varepsilon_{lmk}b_l\delta_{mj} = \varepsilon_{ijk}\varepsilon_{ljk}b_l =
   2\delta_{il}b_l = 2b_i$ (the 2δ contraction of D09): **the curl of a solid-body rotation is twice its angular
   velocity** — Ch. 3's vorticity = 2Ω.
""")
D("D12", "The three components of the curl", ref="2.25",
  goal="""Turn the index form of the curl into the three explicit component formulas — and see why only two of the nine
  terms survive for each component.""",
  assumptions="u differentiable (the partial derivatives exist; every step).",
  start=(r"(\nabla\times\mathbf u)_i = \varepsilon_{ijk}\,\frac{\partial u_k}{\partial x_j}",
         "(2.24): the cross product (2.21) with the first vector replaced by the operator ∇ = e_j ∂/∂x_j, acting to the right."),
  plan=["Fix i = 1.", "Keep only the (j, k) pairs where ε ≠ 0.", "Insert ±1.", "Relabel cyclically for i = 2, 3."],
  uses=["enumeration of ε (2.18) (C08)", "the cross product in index form (2.21) (N49)",
        "operator ordering: ∂ acts on u_k (gloss in the start line)", "partial derivative (Ch. 1 P25)"],
  steps=[
      ("Write (2.24) for i = 1 with both sums",
       r"(\nabla\times\mathbf u)_1 = \sum_{j=1}^{3}\sum_{k=1}^{3}\varepsilon_{1jk}\,\frac{\partial u_k}{\partial x_j}",
       "The summation convention (C01) hides a double sum over j and k — nine terms. We write it out to see which terms are nonzero.",
       "Nine candidate terms for the first component."),
      ("Drop the terms with a repeated index",
       r"(\nabla\times\mathbf u)_1 = \varepsilon_{123}\,\frac{\partial u_3}{\partial x_2} + \varepsilon_{132}\,\frac{\partial u_2}{\partial x_3}",
       "ε vanishes whenever two indices agree (2.18); with i = 1 fixed, only (j, k) = (2, 3) and (3, 2) avoid a repeat. Seven terms die.",
       "Only the two \"other\" directions contribute."),
      ("Insert the values of ε",
       r"(\nabla\times\mathbf u)_1 = \frac{\partial u_3}{\partial x_2} - \frac{\partial u_2}{\partial x_3}",
       r"123 is cyclic, so $\varepsilon_{123} = +1$; 132 is anticyclic, so $\varepsilon_{132} = -1$ (2.18). The minus sign comes from ε, not from the derivative.",
       "The first component is a difference of two cross-derivatives."),
      ("Relabel cyclically for i = 2",
       r"(\nabla\times\mathbf u)_2 = \frac{\partial u_1}{\partial x_3} - \frac{\partial u_3}{\partial x_1}",
       "Replacing every index by its cyclic successor 1 → 2 → 3 → 1 leaves ε unchanged (a two-place move, N43), so step 3 maps to this line. The survivors are now (j, k) = (3, 1) and (1, 3).",
       "Second component: the same pattern one step round the cycle."),
      ("Relabel once more for i = 3 and assemble",
       r"(\nabla\times\mathbf u)_3 = \frac{\partial u_2}{\partial x_1} - \frac{\partial u_1}{\partial x_2}",
       r"Another cyclic step; survivors (1, 2) and (2, 1). The three lines together are (2.25). Throughout, ∂/∂x_j acts on u_k — $u_k\,\partial_j$ would be an operator, not a number.",
       "Third component: does u₂ grow along 1 more than u₁ grows along 2?"),
  ],
  result=(r"(\nabla\times\mathbf u)_1 = \partial_2u_3 - \partial_3u_2,\quad (\nabla\times\mathbf u)_2 = \partial_3u_1 - \partial_1u_3,\quad (\nabla\times\mathbf u)_3 = \partial_1u_2 - \partial_2u_1",
          "each component of the curl is the \"cross-derivative\" of the other two velocity components."),
  interpret="""The curl measures the local spin: a paddle wheel turns at half of it (D15 + Ex. 2.3). It is nonzero for a
  straight shear flow and zero for the irrotational vortex away from its core — spin is about differences of speed
  across a point, not about curved paths. Requires differentiable u; at a vortex core the formula is undefined and the
  circulation (D26) takes over.""",
  check="""Units: (m/s)/m = 1/s ✓. Number (Ex. 2.3, u = b × x = (−x₂, x₁, 0)): $(\\nabla\\times\\mathbf u)_3 = 1 - (-1) = 2
  = 2b_3$, the other two 0 ✓ (`curl_components` = `curl` below). Shear flow $u_1 = \\Gamma x_2$: $(\\nabla\\times\\mathbf u)_3
  = 0 - \\Gamma = -\\Gamma$ ✓.""",
  traps="""writing $u_k\\partial_j$. Losing the minus from $\\varepsilon_{132}$. Doing i = 1 and saying "similarly" — the
  cyclic relabelling is the reason, and it is said.""")
note("N53", "The three components (2.25)", r"""
D12's result; `ch02.curl_components(u, h)` writes them out and is asserted equal to the einsum form `ch02.curl` below.
""", equation=r"(\nabla\times\mathbf u)_1 = \frac{\partial u_3}{\partial x_2} - \frac{\partial u_2}{\partial x_3}, \quad (\nabla\times\mathbf u)_2 = \frac{\partial u_1}{\partial x_3} - \frac{\partial u_3}{\partial x_1}, \quad (\nabla\times\mathbf u)_3 = \frac{\partial u_2}{\partial x_1} - \frac{\partial u_1}{\partial x_2}", ref="2.25")
note("N54", "Solenoidal and irrotational", r"""
**Solenoidal** (∇·u = 0; the word comes from magnetism — a field with no sources) and **irrotational** (∇×u = 0).
Incompressible flow (Ch. 4) is solenoidal; potential flow (Ch. 6) is both; Ch. 3 explains why "irrotational" means no
local spin. `ch02.is_solenoidal`, `ch02.is_irrotational`.
""")
note("N55", "Ex. 2.3, stated", r"""
$a\mathbf x$: divergence 3a, curl 0; $\mathbf b\times\mathbf x$: divergence 0, curl 2b — the worked numbers of C10 and
here, by hand above and by sympy below; pointer: vorticity = 2Ω (Ch. 3 §3.5).
""")
nb.worked_example("Ex. 2.3 with b = e₃: u = (−x₂, x₁, 0)", r"""
1. $(\nabla\times\mathbf u)_3 = \partial u_2/\partial x_1 - \partial u_1/\partial x_2 = 1 - (-1) = 2$.
2. $(\nabla\times\mathbf u)_1 = \partial u_3/\partial x_2 - \partial u_2/\partial x_3 = 0 - 0 = 0$, likewise component 2.
3. $\nabla\times\mathbf u = (0, 0, 2) = 2\mathbf b$ ✓.
4. A straight shear flow $u_1 = \Gamma x_2$: $(\nabla\times\mathbf u)_3 = 0 - \Gamma = -\Gamma$ — curl without a single
   curved streamline (the paddle wheel turns clockwise for Γ > 0).
""")
nb.code("""
w = ch02.curl(u_rot, g3.h)                                                     # (2.24) on the grid: ε_ijk ∂u_k/∂x_j, shape (3, n, n, n)
print(w[:, 5, 5, 5], np.abs(w[2] - 2).max())                                   # [0 0 2] and ~1e-14: ∇×(b × x) = 2b everywhere
print(np.allclose(np.stack(ch02.curl_components(u_rot, g3.h)), w))             # True: (2.25) written out = the einsum form (N53)
print(ch02.is_irrotational(u_rad, g3.h)[0], ch02.is_irrotational(u_rot, g3.h)[0])   # True False: a x has no spin, b × x does
print(ch02.exact_div_curl(sp.Matrix([0, 0, 1]).cross(X3)))                     # the sympy twin: (0, [0, 0, 2])
Xs = ch02.coordinates(3)                                                       # sympy symbols x1, x2, x3 (the field's coordinates)
shear3 = ch02.VectorField([Xs[1], 0, 0], Xs, name="simple shear u1 = Γ x2, Γ = 1")   # the straight shear flow as a 3-D field (Γ = 1/s)
u_sh = shear3(g3.X, g3.Y, g3.Z)                                                # on the grid (kept for C12)
print(ch02.curl(u_sh, g3.h)[2, 5, 5, 5])                                       # −1.0 = −Γ: straight streamlines, nonzero curl
""", explain="""
1. The curl of the rotation field is (0, 0, 2b) at every node — exact for a linear field.
2. The three explicit components (2.25) agree with the einsum form (2.24).
3. The radial field is irrotational, the rotating one is not.
4. The sympy twin gives (0, [0, 0, 2]) symbolically.
5. The straight shear flow $u_1 = \\Gamma x_2$ (built as a `VectorField` from sympy expressions) has curl −Γ.
""")
nb.check_agree("""
dx = lambda f, c: ch02.partial(f, c, g3.h[c])                                  # ∂f/∂x_c along coordinate direction c (Ch. 1 P29: a one-line function)
w_mine = np.stack([dx(u_rot[2], 1) - dx(u_rot[1], 2),                          # (∇×u)_1 = ∂u_3/∂x_2 − ∂u_2/∂x_3   (2.25)
                   dx(u_rot[0], 2) - dx(u_rot[2], 0),                          # (∇×u)_2 = ∂u_1/∂x_3 − ∂u_3/∂x_1
                   dx(u_rot[1], 0) - dx(u_rot[0], 1)])                         # (∇×u)_3 = ∂u_2/∂x_1 − ∂u_1/∂x_2
assert np.allclose(w_mine, w)                                                  # = the library's einsum curl
print("the three cross-derivative differences of D12 = curl()")
""")
nb.figure("""
fields = [(u_rot, "u = b × x: (∇×u)₃ = 2b (spins everywhere)"),                 # the whirlpool
          (u_sh, "u₁ = Γx₂: (∇×u)₃ = −Γ (straight, yet it spins)"),              # the straight shear flow
          (u_rad, "u = a x: (∇×u)₃ = 0 (spreads, no spin)")]                     # the source
fig, axs = plt.subplots(1, 3, figsize=(12, 4.2))                                # three heatmaps
for ax, (U, ttl) in zip(axs, fields):                                           # one panel per field
    c3 = ch02.curl(U, g3.h)[2, k0]                                              # the z = 0 slice of (∇×u)_3
    im = ax.imshow(c3, origin="lower", extent=(-1, 1, -1, 1), cmap="PRGn_r", vmin=-2.2, vmax=2.2)   # purple positive, green negative
    ax.quiver(g3.X[k0, ::sk, ::sk], g3.Y[k0, ::sk, ::sk], U[0][k0, ::sk, ::sk], U[1][k0, ::sk, ::sk], color=COLORS["ink"], scale=25, width=0.003)   # the flow's arrows
    wheel = plt.Circle((0.55, -0.55), 0.18, fill=False, color=COLORS["rose"], lw=2); ax.add_patch(wheel)   # a paddle wheel glyph
    spin = float(c3[k0 // 2, -k0 // 2 - 1]) if c3.shape[0] > 4 else 0.0        # its local (∇×u)₃ sets the arrow's sense
    if abs(spin) > 1e-6:                                                        # curved arrow: counterclockwise for positive curl, clockwise for negative
        ax.annotate("", (0.55 + 0.18 * np.cos(2.2), -0.55 + 0.18 * np.sin(2.2)), (0.55 + 0.18 * np.cos(0.8), -0.55 + 0.18 * np.sin(0.8)),   # from angle 0.8 to 2.2 rad on the wheel
                    arrowprops=dict(arrowstyle="->", color=COLORS["rose"], lw=2, connectionstyle=f"arc3,rad={0.5 if spin > 0 else -0.5}"))   # the arc's sense follows the sign
    ax.set_title(ttl, fontsize=9); ax.set_xlabel("x [m]"); ax.set_aspect("equal")   # labels
fig.colorbar(im, ax=axs, shrink=0.8, label="(∇×u)₃ [1/s]")                    # one shared colour scale
plt.show()                                                                      # display
""", see="Three heatmaps of (∇×u)₃ with the flow's arrows: uniform purple for the whirlpool (+2), uniform green for the straight shear flow (−1), white for the source (0); a paddle-wheel glyph shows the spin sense.",
   read="Curl is the *difference* of speeds across the wheel, not curvature of streamlines: the straight shear flow is as coloured as the whirlpool, and the spreading source has none.",
   change="…you used the irrotational vortex $u_\\theta = K/r$ (the explainer's preset): curl 0 everywhere except the centre — a whirlpool with no local spin, the case Stokes' theorem (C16) cannot handle across the core.")
nb.explainer("stokes_circulation_loop", heading="How much does the flow go round a loop?",
             why="""Curl and circulation are two views of one thing; dragging a loop over a shear flow (curl without curves)
             and over an irrotational vortex (curves without curl) removes the confusion no static figure can. (Its
             theorem, Stokes', is proved in C16 — this explainer is embedded once, here where the curl first bites.)""",
             tries=["Shear preset: the streamlines are straight, yet the wheel turns and Γ ≠ 0.",
                    "Irrotational-vortex preset: shrink the loop away from the centre — Γ/A → 0; enclose the centre — Γ = 2πK whatever the size (the ⚠️ status).",
                    "Flip the orientation: both sides of (2.34) change sign together.",
                    "Rectangle mode: read the four side sums of Ex. 2.6 in the Explain tab."])
nb.md(r"""
**What would change if…** you split the velocity gradient $G_{ij}$ into its symmetric and antisymmetric halves? The
antisymmetric half *is* the curl, packed into a matrix — C12.
""")

# =====================================================================================================================
# A.8 §2.10 Symmetric and Antisymmetric Tensors — C12 (D14, D15, E3)
# =====================================================================================================================
nb.section("2.10", "Symmetric and Antisymmetric Tensors", intro="""
**What is this section about?** Every second-order tensor splits uniquely into a symmetric part (6 numbers) and an
antisymmetric part (3 numbers). The antisymmetric part is a vector in disguise, and a symmetric tensor cannot see it —
the two facts Ch. 3 (strain vs rotation) and Ch. 4 (dissipation) are built on.
""")
core("C12", "Every tensor is a symmetric plus an antisymmetric part — and the antisymmetric part is a vector", r"""
A fluid element both stretches and spins. How does the velocity gradient split those two motions apart — and why does
the spin hide inside three numbers?
""")
nb.md(r"""
#### The problem in plain words

Ch. 3 will take the velocity gradient $G_{ij} = \partial u_i/\partial x_j$ (nine numbers) and ask: how fast is the fluid
element being deformed, and how fast is it turning? The answer is a split: $\mathbf G = \mathbf S + \mathbf A$, S
symmetric (the strain rate — deformation), A antisymmetric (the rotation rate — spin). Ch. 4 then shows viscous friction
costs energy only through S. So the split is not algebra for its own sake; it separates the motion that heats the fluid
from the motion that does not.
""")
nb.md(r"""
#### The idea

```
B = ½(B + Bᵀ) + ½(B − Bᵀ) = S + A       S_ij = S_ji (6 free numbers)     A_ij = −A_ji, A_ii = 0 (3 free numbers)
R =  [  0   −ω₃   ω₂ ]                  three numbers → one vector ω ;  R·x = ω × x  (R rotates x about ω)
     [  ω₃   0   −ω₁ ]
     [ −ω₂   ω₁   0  ]
simple shear u₁ = Γx₂ :  G = [[0, Γ],[0, 0]] = ½[[0, Γ],[Γ, 0]] + ½[[0, Γ],[−Γ, 0]]  = stretch at 45° + spin at Γ/2
```

*Two closely related objects, kept apart on purpose:* the book's **rotation tensor** $R_{ij} = \partial u_i/\partial
x_j - \partial u_j/\partial x_i = G - G^{\rm T}$ (no ½) packs the **vorticity** $\boldsymbol\omega = \nabla\times\mathbf u$
by (2.26)–(2.27); the **antisymmetric part** $\mathbf A = \tfrac12(\mathbf G - \mathbf G^{\rm T}) = \tfrac12\mathbf R$
packs $\tfrac12\boldsymbol\omega = \tfrac12\nabla\times\mathbf u$, the angular velocity of a material line.
""")
P("P79", "scipy.linalg.expm", r"""
For a linear velocity field u = G·x the trajectory is $\mathbf x(t) = e^{\mathbf Gt}\mathbf x_0$, where the matrix
exponential $e^{\mathbf Gt} = \mathbf I + \mathbf Gt + (\mathbf Gt)^2/2 + \dots$. `scipy.linalg.expm(G*t)` computes it;
for a pure rotation matrix it returns the rotation by angle |ω|t. `ch02.linear_flow_map(G, t)` wraps it.
""", code="""
from scipy.linalg import expm                              # the matrix exponential
G_rot = np.array([[0., -1.], [1., 0.]])                    # u = (−y, x): solid-body rotation at 1 rad/s
print(np.round(expm(G_rot * np.pi / 2), 12))               # rotation by 90°: [[0, −1], [1, 0]]
""")
nb.md(r"""
#### The maths, step by step

1. Add and subtract $\tfrac12B_{ji}$: $B_{ij} = \tfrac12(B_{ij} + B_{ji}) + \tfrac12(B_{ij} - B_{ji}) = S_{ij} + A_{ij}$
   — D14 shows the split is unique and frame-independent.
2. Counting: a symmetric 3×3 has 3 diagonal + 3 off-diagonal = 6 independent entries; an antisymmetric one has zero
   diagonal and 3 (N56).
3. Pack a vector into an antisymmetric tensor: $R_{ij} = -\varepsilon_{ijk}\omega_k$ (2.26)–(2.27) and unpack it:
   $\omega_k = -\tfrac12\varepsilon_{ijk}R_{ij}$ — D15, which also shows $\mathbf R\cdot\mathbf x = \boldsymbol\omega
   \times\mathbf x$.
4. A symmetric τ against any B: $\tau_{ij}B_{ij} = \tau_{ij}S_{ij}$, because $\tau_{ij}A_{ij} = 0$ (N59–N61).
""")
note("N56", "Definitions", r"""
symmetric $B_{ij} = B_{ji}$ (6 independent), antisymmetric $B_{ij} = -B_{ji}$ (zero diagonal, 3 independent);
`ch02.independent_components` → 6 / 3 / 9.
""")
D("D14", "The unique split into symmetric and antisymmetric parts", ref="§2.10",
  goal="""Show that any second-order tensor is the sum of a symmetric and an antisymmetric tensor in exactly one way,
  and that both parts are themselves tensors — the basis of Ch. 3's strain rate and rotation rate.""",
  assumptions="None beyond B being a tensor (step 6).",
  start=(r"B_{ij}\ \text{— any second-order tensor}", "nine numbers that transform by (2.12)."),
  plan=["Add and subtract half the transpose.", "Check the symmetry of each half.",
        "Suppose a second split and show it coincides.", "Transform each half with (2.12)."],
  uses=[r"transpose $B_{ji}$ (P63)", "the definitions of symmetric and antisymmetric (N56)", "linearity of (2.12) (C06)",
        "dummy renaming (N12)", "\"both symmetric and antisymmetric ⇒ zero\" (step 5)"],
  steps=[
      (r"Add and subtract $\tfrac12B_{ji}$",
       r"B_{ij} = \tfrac12(B_{ij} + B_{ji}) + \tfrac12(B_{ij} - B_{ji})",
       r"The right side is $\tfrac12B_{ij} + \tfrac12B_{ij}$ plus and minus the same $\tfrac12B_{ji}$ — an identity. We choose the transpose because swapping i and j is the operation whose behaviour we want to separate.",
       "Split B into \"the part that ignores the swap\" and \"the part that flips\"."),
      ("Name the first half and swap its indices",
       r"S_{ij} \equiv \tfrac12(B_{ij} + B_{ji}),\qquad S_{ji} = S_{ij}",
       r"Swapping i ↔ j in the definition gives $\tfrac12(B_{ji} + B_{ij})$, the same sum in the other order. S is symmetric by construction.",
       "The first half does not notice the swap."),
      ("Name the second half and swap its indices",
       r"A_{ij} \equiv \tfrac12(B_{ij} - B_{ji}),\qquad A_{ji} = -A_{ij},\quad A_{ii} = 0",
       r"Swapping gives $\tfrac12(B_{ji} - B_{ij}) = -A_{ij}$; for i = j the difference vanishes. A is antisymmetric with zero diagonal — three independent entries.",
       "The second half changes sign under the swap and has an empty diagonal."),
      ("Suppose another split and subtract",
       r"B = \tilde S + \tilde A\ \Rightarrow\ \tilde S_{ij} - S_{ij} = A_{ij} - \tilde A_{ij} \equiv X_{ij}",
       r"If a second symmetric $\tilde S$ and antisymmetric $\tilde A$ also add to B, subtracting the two splits and moving terms gives one object X written two ways. This is the uniqueness half the book omits.",
       "Any other split differs from ours by one tensor X."),
      ("Show X is zero",
       r"X_{ij} = X_{ji} = -X_{ij}\ \Rightarrow\ X_{ij} = 0",
       r"X is a difference of symmetric tensors (so symmetric) and a difference of antisymmetric ones (so antisymmetric); a number equal to its own negative is zero. Hence $\tilde S = S$, $\tilde A = A$.",
       "Only the zero tensor is both symmetric and antisymmetric — the split is unique."),
      ("Transform S with (2.12)",
       r"S'_{mn} = \tfrac12(B'_{mn} + B'_{nm}) = C_{im}C_{jn}\,\tfrac12(B_{ij} + B_{ji}) = C_{im}C_{jn}S_{ij}",
       r"(2.12) is linear, so it applies to each term; in the second term $B'_{nm} = C_{in}C_{jm}B_{ij}$, and renaming the dummies i ↔ j turns it into $C_{im}C_{jn}B_{ji}$. The same lines with a minus give A.",
       "The symmetric part of the transformed B is the transformed symmetric part — S and A are tensors."),
  ],
  result=(r"B_{ij} = S_{ij} + A_{ij},\quad \mathbf S = \tfrac12(\mathbf B + \mathbf B^{\rm T}),\quad \mathbf A = \tfrac12(\mathbf B - \mathbf B^{\rm T})\ \text{— unique, both tensors; } 6 + 3 = 9",
          "every tensor is a stretch-like part plus a spin-like part, and every observer agrees on the split."),
  interpret="""For the velocity gradient $\\partial u_i/\\partial x_j$, S is the strain rate (how a fluid element deforms)
  and A the rotation rate (how it spins) — Ch. 3 §3.4; Ch. 4's viscous dissipation sees only S (N59–N61).
  Frame-independence (step 6) is what lets us compute the split in any convenient axes.""",
  check="""Units: those of B ✓. Limit: B symmetric ⇒ A = 0; B antisymmetric ⇒ S = 0 ✓. Number: B = [[0, 2], [0, 0]]:
  S = [[0, 1], [1, 0]], A = [[0, 1], [−1, 0]], S + A = B ✓ (`symmetric_part`, `antisymmetric_part` below);
  `independent_components` 6, 3, 9.""",
  traps="""thinking the split depends on the axes. Forgetting the uniqueness half. Counting 3 + 3 = 6 for the
  antisymmetric part (the diagonal is zero: three numbers).""")
note("N57", "The antisymmetric tensor of a vector ω", r"""
`ch02.antisymmetric_from_vector(omega)` builds it. Pointer: with $\boldsymbol\omega = \nabla\times\mathbf u$ this R is
Ch. 3's rotation tensor $R_{ij} = \partial u_i/\partial x_j - \partial u_j/\partial x_i$ (`ch02.rotation_tensor(G)`,
no ½).
""", equation=r"\mathbf R = \begin{bmatrix}0 & -\omega_3 & \omega_2\\ \omega_3 & 0 & -\omega_1\\ -\omega_2 & \omega_1 & 0\end{bmatrix}", ref="2.26")
note("N58", "The two-way map", r"""
the book prints the lower limit of the first sum as "i−1"; it is i = 1. Round trip in code to 1e-16
(`ch02.vector_from_antisymmetric`); and $\mathbf R\cdot\mathbf x = \boldsymbol\omega\times\mathbf x$ (D15).
""", equation=r"R_{ij} = -\varepsilon_{ijk}\omega_k, \qquad \omega_k = -\tfrac12\,\varepsilon_{ijk}R_{ij}", ref="2.27")
D("D15", "The vector hidden in an antisymmetric tensor", ref="2.26, 2.27",
  goal="""Show that the three independent numbers of an antisymmetric tensor form a vector ω, give the two-way map, and
  show that R acting on x rotates it: R·x = ω × x. The book checks one entry; the inverse and the sign are ours.""",
  assumptions="Three dimensions (ε has three indices; everywhere).",
  start=(r"R_{ij} \equiv -\varepsilon_{ijk}\,\omega_k", "(2.27, first part): pack a vector into a 3×3 array with the alternating tensor."),
  plan=["Check R is antisymmetric and read off its entries.", "Contract R with ε to get ω back (D09's 2δ).",
        "Apply R to a vector and recognise a cross product."],
  uses=["definition and index moves of ε (N43)", r"the 2δ contraction $\varepsilon_{ijl}\varepsilon_{ijk} = 2\delta_{lk}$ (D09 step 8)",
        "Kronecker substitution (N41)", "the cross product in index form (2.21) (N49)"],
  steps=[
      ("Swap the indices of R",
       r"R_{ji} = -\varepsilon_{jik}\,\omega_k = +\varepsilon_{ijk}\,\omega_k = -R_{ij}",
       "Swapping two indices of ε flips its sign (one-place move, N43). So R is antisymmetric — the object we want to represent.",
       "The packed array is antisymmetric, as promised."),
      ("Read off the entries",
       r"R_{12} = -\varepsilon_{123}\omega_3 = -\omega_3,\quad R_{13} = -\varepsilon_{132}\omega_2 = +\omega_2,\quad R_{23} = -\varepsilon_{231}\omega_1 = -\omega_1",
       "For each (i, j) only the k different from both survives; ε_123 = ε_231 = +1, ε_132 = −1 (2.18). The diagonal is zero by step 1. With step 1 for the lower triangle this is the matrix (2.26).",
       "Above the diagonal: −ω₃, +ω₂, −ω₁."),
      ("Contract R with ε to invert",
       r"\varepsilon_{ijl}\,R_{ij} = -\varepsilon_{ijl}\,\varepsilon_{ijk}\,\omega_k",
       r"Multiply the definition by $\varepsilon_{ijl}$ and sum on i, j. We do this because a product of two ε's sharing two indices collapses to a δ (D09).",
       "Contract both sides with the alternating tensor."),
      ("Apply the 2δ contraction",
       r"\varepsilon_{ijl}\,R_{ij} = -2\,\delta_{lk}\,\omega_k = -2\,\omega_l",
       r"$\varepsilon_{ijl}\varepsilon_{ijk} = 2\delta_{lk}$ (D09 step 8, the book's $\varepsilon_{pqi}\varepsilon_{pqj} = 2\delta_{ij}$), then δ substitutes.",
       "The contraction returns minus twice the vector."),
      ("Solve for ω and rename",
       r"\omega_k = -\tfrac12\,\varepsilon_{ijk}\,R_{ij}",
       "Divide by −2 and rename the free index l → k. This is the second half of (2.27); the sum runs over i = 1..3 and j = 1..3 (the book's lower limit \"i−1\" is a misprint for i = 1).",
       "Half the contraction of R with ε, sign reversed, gives ω back."),
      ("Apply R to a vector x",
       r"R_{ij}\,x_j = -\varepsilon_{ijk}\,\omega_k\,x_j = +\varepsilon_{ikj}\,\omega_k\,x_j",
       "Insert the definition; then swap the last two indices of ε (one place, sign flips) to bring k before j — the order the cross-product formula needs.",
       "R times x is ε contracted with ω and x in the cross-product order."),
      ("Recognise the cross product",
       r"R_{ij}\,x_j = (\boldsymbol\omega\times\mathbf x)_i",
       r"(2.21) in the form $(\mathbf u\times\mathbf v)_i = \varepsilon_{ikj}u_kv_j$ (free index first, then u's, then v's — a two-place move of $\varepsilon_{kji}$). With u = ω, v = x this is step 6.",
       "R rotates x about ω: R·x = ω × x."),
  ],
  result=(r"R_{ij} = -\varepsilon_{ijk}\omega_k,\qquad \omega_k = -\tfrac12\varepsilon_{ijk}R_{ij},\qquad \mathbf R\cdot\mathbf x = \boldsymbol\omega\times\mathbf x",
          "an antisymmetric tensor *is* a vector, and acting with it means \"cross with that vector\"."),
  interpret="""For the velocity gradient $\\mathbf G$, the book's rotation tensor $\\mathbf R = \\mathbf G - \\mathbf G^{\\rm T}$
  (no ½) packs exactly the vorticity $\\boldsymbol\\omega = \\nabla\\times\\mathbf u$ (Ch. 3, (3.15)/(3.17)); the antisymmetric
  part $\\mathbf A = \\tfrac12\\mathbf R$ therefore carries $\\tfrac12\\boldsymbol\\omega = \\tfrac12\\nabla\\times\\mathbf u$, and
  $\\mathbf A\\cdot\\mathbf x = \\tfrac12\\boldsymbol\\omega\\times\\mathbf x$ is the velocity of a solid-body rotation at angular
  velocity $\\tfrac12\\boldsymbol\\omega$ (Ex. 2.3 read backwards). That is *why* the antisymmetric part means "spin": a
  material line turns at half the vorticity. The sign in (2.27) is chosen so that this works; with the opposite sign
  R·x = −ω × x.""",
  check="""Units: those of ω ✓. Round trip: ω = (1, 2, 3) → R → ω to 1e-16 (`antisymmetric_from_vector`,
  `vector_from_antisymmetric`). Number: ω = (0, 0, −1), x = (1, 0, 0): R = [[0, 1, 0], [−1, 0, 0], [0, 0, 0]], R·x =
  (0, −1, 0); ω × x = (0·0 − (−1)·0, (−1)·1 − 0·0, 0) = (0, −1, 0) ✓ (the sign-discrimination test).""",
  traps="""the sign in step 6 (a one-place move). The lower limit misprint. Forgetting that $\\delta_{lk}\\omega_k =
  \\omega_l$ (not 3ω). Calling the vector of A "ω": it is ½ω — the book's R has no ½, A does.""")
note("N59", "A symmetric tensor against a general one", r"""
`ch02.symmetric_double_contraction(τ, B)` returns (P, P_S, P_A) with the Frobenius pairing $\tau_{kl}B_{kl}$.
""", equation=r"P = \tau_{kl}B_{kl} = \tau_{kl}(S_{kl} + A_{kl}) = \tau_{ij}S_{ij} + \tau_{ij}A_{ij}", ref="2.28")
note("N60", "The antisymmetric part drops out (stated)", r"""
swapping A's indices and renaming dummies gives $P = \tau_{kl}S_{kl} - \tau_{kl}A_{kl}$; comparing with (2.28), $X =
\tau_{ij}A_{ij}$ satisfies $X = -X$, so $X = 0$. **Hence** `N61`: $\tau_{ij}B_{ij} = \tau_{ij}S_{ij} = \tfrac12
\tau_{ij}(B_{ij} + B_{ji})$ — like the integral of even × odd over a symmetric interval vanishing; Ch. 4 §4.5: the
dissipation $\tau_{ij}\,\partial u_i/\partial x_j$ sees only the strain rate.
""", equation=r"P = \tau_{kl}S_{kl} - \tau_{kl}A_{kl}", ref="2.29")
nb.worked_example("simple shear u₁ = Γx₂ with Γ = 2 s⁻¹", r"""
1. $\mathbf G = \begin{bmatrix}0 & 2\\ 0 & 0\end{bmatrix}$ s⁻¹ ($G_{12} = \partial u_1/\partial x_2$).
2. $\mathbf S = \tfrac12(\mathbf G + \mathbf G^{\rm T}) = \begin{bmatrix}0 & 1\\ 1 & 0\end{bmatrix}$, $\mathbf A =
   \tfrac12(\mathbf G - \mathbf G^{\rm T}) = \begin{bmatrix}0 & 1\\ -1 & 0\end{bmatrix}$; S + A = G ✓.
3. The vector of A (3-D, $A_{12} = -(\tfrac12\omega)_3$): $(\tfrac12\omega)_3 = -A_{12} = -1$ s⁻¹ — a material line turns
   clockwise at 1 rad/s. So the vorticity is $\omega_3 = -2$ s⁻¹ $= -\Gamma$ — exactly C11's curl of the shear flow ✓;
   the book's $\mathbf R = \mathbf G - \mathbf G^{\rm T} = 2\mathbf A$ packs this ω directly.
4. Check $\mathbf A\cdot\mathbf x = \tfrac12\boldsymbol\omega\times\mathbf x$ at x = (1, 0, 0): A·x = (0, −1, 0);
   $(0, 0, -1)\times(1, 0, 0) = (0\cdot0 - (-1)\cdot0,\ (-1)\cdot1 - 0\cdot0,\ 0) = (0, -1, 0)$ ✓.
5. S:S = 2, A:A (Frobenius) = 2, S:A = 0 — equal parts stretch and spin, orthogonal to each other.
""")
nb.code("""
G = ch02.velocity_gradient_preset("simple_shear", Gamma=2.0, dim=3)  # G[i, j] = ∂u_i/∂x_j of u_1 = Γ x_2: [[0, 2, 0], [0, 0, 0], [0, 0, 0]] s⁻¹
S, A_ = ch02.symmetric_part(G), ch02.antisymmetric_part(G)          # the two halves (D14): S = ½(G + Gᵀ), A = ½(G − Gᵀ)
print(S); print(A_)                                                 # the tiny example's matrices, padded to 3×3
assert np.allclose(S + A_, G)                                        # they add back to G
half_omega = ch02.vector_from_antisymmetric(A_)                     # (2.27) on A: the vector of A = ½ω = ½∇×u
print(half_omega)                                                   # [0 0 −1] s⁻¹: a material line turns clockwise at 1 rad/s
print(np.allclose(ch02.antisymmetric_from_vector(half_omega), A_))  # True: (2.26) rebuilds A from its vector
R_book = ch02.rotation_tensor(G)                                    # the book's rotation tensor R = G − Gᵀ = 2A (no ½)
print(ch02.vector_from_antisymmetric(R_book))                       # [0 0 −2] = ω = ∇×u: the vorticity, −Γ
x = np.array([1., 0., 0.])                                          # a test vector
print(A_ @ x, np.cross(half_omega, x))                              # [0 −1 0] twice: A·x = ½ω × x (D15 step 7)
print(ch02.independent_components(S), ch02.independent_components(A_), ch02.independent_components(G))   # 6 3 9 (N56)
print(ch02.symmetric_double_contraction(S, G))                      # (P, P_S, P_A) = (2, 2, 0): a symmetric tensor cannot see A (N59–N61)
Gnum = ch02.vector_gradient(u_sh, g3.h)[:, :, 5, 5, 5]              # the NUMERICAL velocity gradient of C11's shear flow (Γ = 1) at one node
print(ch02.vector_from_antisymmetric(ch02.antisymmetric_part(Gnum)), 0.5 * ch02.curl(u_sh, g3.h)[:, 5, 5, 5])   # [0 0 −0.5] twice: vector(A) = ½∇×u — the sign convention pinned
""", explain="""
1. The preset G of simple shear with $\\partial u_1/\\partial x_2 = 2$ s⁻¹.
2–3. The two parts and their sum.
4. The vector of A by (2.27) is ½ω = (0, 0, −1) s⁻¹ — never call it ω: the book's R = G − Gᵀ = 2A carries the vorticity
   ω = (0, 0, −2) s⁻¹ itself (`rotation_tensor`), and −2 = −Γ is C11's curl.
5. (2.26) inverts the map exactly.
6. A·x = ½ω × x — the antisymmetric part *acts* as a rotation at angular velocity ½ω.
7. Component counts 6, 3, 9.
8. (2.28)–(2.29): P = P_S and P_A = 0.
9. The same relation on a *numerical* velocity gradient from the grid of C11: the vector of A is half the curl — the
   convention G[i, j] = ∂u_i/∂x_j is what makes the sign come out right.
""")
nb.check_agree("""
S_mine, A_mine = 0.5 * (G + G.T), 0.5 * (G - G.T)                    # D14 step 1 by hand
half_om_mine = np.array([A_mine[2, 1], A_mine[0, 2], A_mine[1, 0]])  # reading (2.26): ω₁ = R₃₂, ω₂ = R₁₃, ω₃ = R₂₁ — applied to A, so this is ½ω
assert np.allclose(S_mine, S) and np.allclose(A_mine, A_)            # the library's split
assert np.allclose(half_om_mine, half_omega)                         # the library's (2.27)
print("½(G ± Gᵀ) and the three off-diagonal entries of A = the library")
""")
nb.animation("""
n_fr = 48 if not FAST else 24                                             # frames (FAST: 24)
ts = np.linspace(0, 0.8, n_fr)                                            # time [s]; Γ t up to 1.6
G2, S2, A2 = G[:2, :2], S[:2, :2], A_[:2, :2]                             # the 2×2 blocks (the flow is planar)
sq0 = ch02.deform_square(G2, 0.0, 10)                                     # a 10×10 lattice of material points at t = 0, shape (2, 100)
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.8))                        # G | S | A
fig.set_layout_engine("none")                                             # fixed layout: the automatic layout engine would re-run on every frame
titles =("G = S + A: the square leans", "S alone: stretch at +45°, squeeze at −45°", "A alone: rigid spin at ½ω = 1 rad/s")
dots = []                                                                 # the moving lattices
for ax, ttl in zip(axs, titles):
    ax.plot(*sq0, ".", color=COLORS["muted"], ms=2)                       # the t = 0 lattice (grey ghost)
    dots.append(ax.scatter(*sq0, s=8, color=COLORS["accent"]))            # the moving lattice
    ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.6, 2.6); ax.set_aspect("equal"); ax.set_title(ttl, fontsize=9)
    ax.set_xticks([-2, 0, 2]); ax.set_yticks([-2, 0, 2])                  # few tick labels: text is what makes video frames slow to draw
for sgn, col in ((1, COLORS["blue"]), (-1, COLORS["rose"])):              # the eigen-axes of S at ±45° (N62/N63 preview, C13)
    axs[1].plot([-2.2, 2.2], [-2.2 * sgn, 2.2 * sgn], "--", color=col, lw=1)
wheel = plt.Circle((1.6, 1.6), 0.35, fill=False, color=COLORS["rose"], lw=2); axs[2].add_patch(wheel)   # a paddle wheel on the A panel
(spoke,) = axs[2].plot([1.6, 1.95], [1.6, 1.6], color=COLORS["rose"], lw=2)   # its spoke, turning at ½ω
clock = fig.text(0.5, 0.01, "", ha="center")                              # the clock
def update(i):                                                            # frame i
    t = ts[i]
    for d, M in zip(dots, (G2, S2, A2)):                                  # each panel's tensor
        d.set_offsets(ch02.deform_square(M, t, 10).T)                     # x(t) = e^{Mt} x_0 for the whole lattice (P79)
    spoke.set_data([1.6, 1.6 + 0.35 * np.cos(-1.0 * t)], [1.6, 1.6 + 0.35 * np.sin(-1.0 * t)])   # the spoke turns at ½ω₃ = −1 rad/s
    clock.set_text(f"t = {t:.2f} s   (Γ = 2 s⁻¹)")
    return (*dots, spoke, clock)
show_animation(animate(update, frames=n_fr, fig=fig, interval=60), player="video")   # smooth MP4
""", explain="""
1. The same 10×10 square of material points is carried by three linear flows from the same start: the full velocity
   gradient G, its symmetric part S and its antisymmetric part A.
2. `deform_square(M, t)` moves every point exactly along $\\mathbf x(t) = e^{\\mathbf Mt}\\mathbf x_0$ (P79).
3. The dashed lines in the middle panel are the eigen-axes of S at ±45° (C13); the paddle wheel on the right turns at
   the vector of A, ½ω₃ = −1 rad/s (clockwise).
""")
see_read_change("Lean = stretch + spin: under G the square shears into a parallelogram; under S it stretches along +45° and squeezes along −45°, keeping those axes fixed; under A it turns rigidly clockwise at 1 rad/s.",
                "The S panel keeps the axes at 45° fixed and changes lengths (by $e^{\\pm t}$); the A panel keeps lengths and changes angles; G does both at once.",
                "…solid-body preset (`velocity_gradient_preset('solid_body_rotation')`): S = 0, the middle panel does nothing at all, and A turns the square at the full rate Γ.")
nb.explainer("strain_vs_rotation_split", heading="Is simple shear a rotation?",
             why="""The decomposition is a statement about motion: only watching the same square under G, S and A on one
             clock shows that shear is half stretch, half spin — and the term bars make τ:A = 0 a thing you see.""",
             tries=["Simple-shear preset: play; the stretch factor e^{Γt/2} and the turned angle Γt/2 appear on the end card.",
                    "Set G₁₂ = G₂₁: A vanishes, the wheel stops.",
                    "Solid-body preset: S = 0, the square keeps its shape.",
                    "Watch the S:A bar stay at zero while you drag any slider."])
nb.md(r"""
**What would change if…** you looked for the directions in which S *only* stretches, with no shear? Those are its
eigenvectors — C13.
""")

# =====================================================================================================================
# A.9 §2.11 Eigenvalues and Eigenvectors of a Symmetric Tensor — C13 (D17 ★★★)
# =====================================================================================================================
nb.section("2.11", "Eigenvalues and Eigenvectors of a Symmetric Tensor", intro="""
**What is this section about?** For a real symmetric tensor there are three perpendicular directions on which it acts
as a pure stretch; in those axes it is diagonal, and its three eigenvalues bound the normal stress (or strain rate) on
every plane. The book lists the facts; we prove them (D17) and then read Ex. 2.4 by hand.
""")
core("C13", "Principal axes of a symmetric tensor", r"""
Which planes through a point carry *no* shear — and why are the normal stresses on them the largest and smallest of
all?
""")
nb.md(r"""
#### The problem in plain words

A material fails on the plane where the shear is largest; a fluid element stretches fastest along one line and shrinks
along another. The C05 explainer let you *find* the shear-free planes by dragging n. The eigenvalue problem finds them
exactly, for any symmetric tensor, and tells you the extreme normal stresses without searching.
""")
nb.md(r"""
#### The idea

```
ask: is there a plane whose traction is purely normal?   τ·b = λ b   (f parallel to n = b)
→ (τ − λδ)·b = 0 has a nonzero b only if det(τ − λδ) = 0 → cubic in λ → three roots λ¹ λ² λ³ (real!)
→ three perpendicular b's ;  axes along them:  τ' = diag(λ¹, λ², λ³) ;  n·τ·n ∈ [λ_min, λ_max] for every n
```
""")
P("P80", "eigenvalues and eigenvectors", r"""
$\mathbf A\cdot\mathbf b = \lambda\mathbf b$: the matrix only *scales* the vector b by λ. Nonzero b exist exactly when
$\det(\mathbf A - \lambda\mathbf I) = 0$ (🔁 Ch. 1 P53: a zero determinant means dependent columns), a polynomial in λ
of degree n — the characteristic polynomial of D18. `np.linalg.eigh` solves the symmetric case (returns λ ascending and
unit eigenvectors as columns); `np.roots` finds polynomial roots.
""", code="""
M2 = np.array([[2., 1.], [1., 2.]])                        # a symmetric 2×2
lam, Bv = np.linalg.eigh(M2)                               # eigenvalues (ascending) and unit eigenvectors as columns (Ch. 1 P14: tuple unpacking)
print(lam, np.round(Bv, 4))                                # [1. 3.], columns (1,−1)/√2 and (1,1)/√2
print(np.allclose(M2 @ Bv[:, 1], lam[1] * Bv[:, 1]))       # True: A·b = λb
print(np.roots([1, -4, 3]))                                # [3. 1.]: the roots of the characteristic polynomial λ² − 4λ + 3
""")
P("P81", "complex conjugate", r"""
For $z = a + ib$ the conjugate is $\bar z = a - ib$; $z\bar z = a^2 + b^2 = |z|^2 \ge 0$, and $z = \bar z$ exactly when
z is real. Conjugating a product conjugates each factor; a real number is its own conjugate. (🔁 Ch. 1 P45 introduced
i² = −1.)
""", code="""
z = 3 + 4j                                                 # a complex number (Python writes i as j)
print(z.conjugate(), (z * z.conjugate()).real, abs(z)**2)  # (3−4j) 25.0 25.0
""")
P("P82", "quadratic form and the Rayleigh quotient", r"""
$n_i\tau_{ij}n_j = \mathbf n\cdot\boldsymbol\tau\cdot\mathbf n$ is a quadratic form — for a unit n it is the *normal*
component of the traction on the plane ⊥ n. Its values over all unit n range between the smallest and largest
eigenvalue (the Rayleigh quotient bound, proved in D17). *Gloss — repeated eigenvalues:* any orthonormal pair in the
eigenplane works (Gram–Schmidt picks one; `eigh` does it for you).
""", code="""
tau_ = np.array([[0., 1.], [1., 0.]])                      # pure shear (eigenvalues ±1)
for deg in (0, 30, 45, 90):                                # four unit normals
    n_ = np.array([np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))])
    print(deg, np.round(n_ @ tau_ @ n_, 3))                # 0.0 0.866 1.0 0.0 — never outside [−1, 1]
""")
nb.md(r"""
#### The maths, step by step

1. Eigen-equation $\tau_{ij}b_j = \lambda b_i$; nonzero b needs $\det|\tau_{ij} - \lambda\delta_{ij}| = 0$, the cubic
   $\lambda^3 - I_1\lambda^2 + I_2\lambda - I_3 = 0$ (C07, D18).
2. Facts (D17): (1) the three λ are real; (2) eigenvectors of distinct λ are orthogonal; (3) with $\mathbf C =
   [\mathbf b^1\,\mathbf b^2\,\mathbf b^3]$ (columns = new axes, exactly C02's C), $\boldsymbol\tau' = \mathbf C^{\rm T}
   \boldsymbol\tau\mathbf C = \mathrm{diag}(\lambda^1, \lambda^2, \lambda^3)$; (4) for every unit n, $\lambda_{\min}
   \le n_i\tau_{ij}n_j \le \lambda_{\max}$ and the shear on any plane is at most $(\lambda_{\max} - \lambda_{\min})/2$.
3. ⚠️ The book's fact (4) says the λ bound "the elements τ_ij"; the sharp and useful statement is about the *normal
   stress on every plane* and the shear bound (Mohr's circle, Ch. 4).
4. 2-D shortcut: principal angle $\tfrac12\arctan2(2S_{12}, S_{11} - S_{22})$, principal values $\tfrac12(\tau_{11} +
   \tau_{22}) \pm \sqrt{(\tfrac12(\tau_{11} - \tau_{22}))^2 + \tau_{12}^2}$.
""")
D("D17", "Real eigenvalues, orthogonal axes, diagonal form and the bounds — for a real symmetric tensor", ref="§2.11 facts (1)–(4)",
  goal="""Prove the four facts §2.11 states without proof: for a real symmetric τ (1) the eigenvalues are real, (2)
  eigenvectors of distinct eigenvalues are perpendicular, (3) in the axes of the eigenvectors τ is diagonal, (4) the
  normal stress on every plane lies between the smallest and largest eigenvalue (and the shear is bounded by half their
  spread).""",
  assumptions="""τ real (step 3) and symmetric (steps 4 and 7 — used twice, said both times) · the eigenvectors are
  normalised to unit length (steps 9, 12) · a right-handed choice of signs (step 9).""",
  start=(r"\tau_{ij}\,b_j = \lambda\,b_i \quad\text{with}\quad \tau_{ij} = \tau_{ji}\ \text{real}",
         "a direction b that τ only scales (no shear on the plane ⊥ b); λ and b may be complex until we prove otherwise."),
  plan=["Contract the eigen-equation with the conjugate eigenvector, conjugate, use symmetry: λ = λ̄.",
        "Do the same with two eigenpairs: (λ¹ − λ²) b¹·b² = 0.",
        "Put the unit eigenvectors as columns of C — exactly C02's C — and apply (2.12).",
        "Expand any unit n in the eigenbasis and read the normal stress as a weighted mean of the λ's."],
  uses=["eigenvalues and eigenvectors (P80)", "complex conjugate and |z|² = z z̄ (P81)", "symmetry τ_ij = τ_ji (N56)",
        "dummy renaming (N12)", "orthonormality and completeness (P65) and D02", "the transformation rule (2.12) (C06)",
        "quadratic form / Rayleigh quotient (P82)", "Gram–Schmidt for repeated λ (gloss in step 8)"],
  steps=[
      ("Contract the eigen-equation with the conjugate eigenvector",
       r"\bar b_i\,\tau_{ij}\,b_j = \lambda\,\bar b_i b_i",
       r"Multiply both sides by $\bar b_i$ and sum on i; $\bar b_ib_i = \sum|b_i|^2 > 0$ for a nonzero b. We conjugate one factor so that a real number appears on the right.",
       "A complex \"normal stress\" of τ along b equals λ times a positive number."),
      ("Take the complex conjugate of the whole line",
       r"b_i\,\bar\tau_{ij}\,\bar b_j = \bar\lambda\,b_i\bar b_i",
       "Conjugating a product conjugates each factor (P81); an equation stays true under conjugation. We do it to produce λ̄ next to λ.",
       "The mirror-image statement with λ̄."),
      ("Use that τ is real",
       r"b_i\,\tau_{ij}\,\bar b_j = \bar\lambda\,b_i\bar b_i",
       r"$\bar\tau_{ij} = \tau_{ij}$ because the components are real numbers. Assumption \"real\" enters here.",
       "τ is unchanged by the mirror."),
      ("Rename the dummies and use symmetry",
       r"b_i\,\tau_{ij}\,\bar b_j = \bar b_j\,\tau_{ji}\,b_i = \bar b_i\,\tau_{ij}\,b_j",
       r"Reorder the scalar factors, rename i ↔ j (both summed, N12), then $\tau_{ji} = \tau_{ij}$ (symmetry, first use). The left side of step 3 has become the left side of step 1.",
       "Thanks to symmetry the two \"normal stresses\" are the same number."),
      ("Subtract step 3 from step 1",
       r"(\lambda - \bar\lambda)\,\bar b_ib_i = 0\ \Rightarrow\ \lambda = \bar\lambda",
       r"Equal left sides give equal right sides; $\bar b_ib_i > 0$ can be divided out, leaving λ = λ̄ — λ is real (fact 1). A real λ makes $(\tau - \lambda\delta)b = 0$ a real linear system, so b can be chosen real.",
       "The eigenvalues of a real symmetric tensor are real numbers."),
      ("Contract two eigen-equations crosswise",
       r"b^2_i\,\tau_{ij}\,b^1_j = \lambda^1\,b^2_ib^1_i,\qquad b^1_i\,\tau_{ij}\,b^2_j = \lambda^2\,b^1_ib^2_i",
       r"Take eigenpairs (λ¹, b¹) and (λ², b²) (superscripts are labels, not powers); contract the first equation with $b^2_i$ and the second with $b^1_i$. Same trick as step 1 with two vectors.",
       "Two mixed \"normal stresses\", each equal to an eigenvalue times b¹·b²."),
      ("Show the two left sides are equal",
       r"b^1_i\,\tau_{ij}\,b^2_j = b^2_j\,\tau_{ji}\,b^1_i = b^2_i\,\tau_{ij}\,b^1_j",
       r"Rename i ↔ j, then $\tau_{ji} = \tau_{ij}$ (symmetry, second use). Without symmetry this step fails and eigenvectors need not be perpendicular.",
       "τ between b¹ and b² does not care about the order."),
      ("Subtract: orthogonality for distinct eigenvalues",
       r"(\lambda^1 - \lambda^2)\,\mathbf b^1\cdot\mathbf b^2 = 0\ \Rightarrow\ \mathbf b^1\cdot\mathbf b^2 = 0\ \ (\lambda^1 \neq \lambda^2)",
       "Equal left sides in step 6 force the right sides equal. If λ¹ ≠ λ² the dot product must vanish (fact 2). If λ¹ = λ², any vector in the eigenplane is an eigenvector and an orthogonal pair can be chosen (Gram–Schmidt — `eigh` does this).",
       "Principal axes of different eigenvalues are perpendicular; equal ones let us choose them so."),
      ("Build C from the unit eigenvectors",
       r"C_{ij} \equiv b^j_i = \mathbf e_i\cdot\mathbf b^j,\qquad \mathbf C^{\rm T}\mathbf C = \mathbf I",
       r"Take the new axes $\mathbf e'_j = \mathbf b^j$ (normalised, mutually perpendicular by step 8): then $C_{ij} = \mathbf e_i\cdot\mathbf e'_j$ is precisely C02's matrix, orthogonal by D02. Flip the sign of one b if det C = −1, to keep a rotation.",
       "The columns of C are the principal axes."),
      ("Apply the transformation rule (2.12)",
       r"\tau'_{mn} = C_{im}\,C_{jn}\,\tau_{ij} = b^m_i\,\tau_{ij}\,b^n_j",
       "(2.12) for a second-order tensor (D06) with step 9's C inserted. The right side is \"τ between b^m and b^n\" — the quantity of step 6.",
       "The new components are τ sandwiched between two principal axes."),
      ("Use the eigen-equation on b^n",
       r"\tau'_{mn} = \lambda^n\,b^m_ib^n_i = \lambda^n\,\delta_{mn}",
       r"$\tau_{ij}b^n_j = \lambda^nb^n_i$ (no sum on n), then orthonormality $\mathbf b^m\cdot\mathbf b^n = \delta_{mn}$ (steps 8–9). So τ' is diagonal with the eigenvalues on the diagonal (fact 3).",
       "In the principal axes there is no shear at all — only λ¹, λ², λ³."),
      ("Expand a unit normal in the eigenbasis",
       r"\mathbf n = c_k\,\mathbf b^k,\qquad c_k = \mathbf n\cdot\mathbf b^k,\qquad c_kc_k = 1",
       r"Completeness of the orthonormal set {b^k} (P65); $|\mathbf n|^2 = c_kc_k$ by orthonormality. We expand n because τ is simple on the b's.",
       "Any cut direction is a mix of the three principal axes, with weights whose squares add to 1."),
      ("Compute the normal stress on the plane ⊥ n",
       r"n_i\,\tau_{ij}\,n_j = c_kc_l\,\mathbf b^k\cdot\boldsymbol\tau\mathbf b^l = \lambda^k\,c_k^2",
       r"Insert step 12 twice; by steps 6–8 and 11, $\mathbf b^k\cdot\boldsymbol\tau\mathbf b^l = \lambda^l\delta_{kl}$, so only k = l survives. The normal stress is a weighted mean of the eigenvalues with weights $c_k^2$.",
       "The normal stress on any plane is an average of the principal stresses."),
      ("Bound the weighted mean",
       r"\lambda_{\min} \le n_i\tau_{ij}n_j \le \lambda_{\max}",
       r"Replace every λ^k by λ_min (or λ_max) in $\lambda^kc_k^2$: the sum can only decrease (increase), and $c_kc_k = 1$ gives the bound (Rayleigh quotient, P82). Equality when n is the corresponding b^k (fact 4, sharp form).",
       "No cut carries more pull than λ_max or more squeeze than λ_min."),
      ("Bound the shear too",
       r"\tau_s^2 = \lambda^k\lambda^k c_k^2 - (\lambda^kc_k^2)^2 \le \Big(\tfrac{\lambda_{\max} - \lambda_{\min}}{2}\Big)^2",
       r"$|\mathbf f|^2 = f_if_i = \lambda^k\lambda^kc_k^2$ by the same expansion, and $\tau_s^2 = |\mathbf f|^2 - \sigma_n^2$ is the *variance* of the values λ^k under the weights $c_k^2$; a variance of numbers confined to an interval of length L is at most (L/2)².",
       "The shear on any plane is at most half the spread of the eigenvalues (Mohr's circle radius)."),
  ],
  result=(r"\text{real } \lambda^k;\quad \text{orthonormal } \mathbf b^k;\quad \boldsymbol\tau' = \mathbf C^{\rm T}\boldsymbol\tau\mathbf C = \mathrm{diag}(\lambda^1, \lambda^2, \lambda^3),\ \mathbf C = [\mathbf b^1\ \mathbf b^2\ \mathbf b^3];\quad \lambda_{\min} \le \mathbf n\cdot\boldsymbol\tau\cdot\mathbf n \le \lambda_{\max},\ \ \tau_s \le \tfrac12(\lambda_{\max} - \lambda_{\min})",
          "a symmetric tensor is a pure stretch along three perpendicular axes, and those stretches bound what any plane can feel."),
  interpret="""Ch. 3: the strain rate has three perpendicular directions of pure stretching (no shearing), with rates
  $\\lambda^k$; Ch. 4: principal stresses and Mohr's circle; Ch. 11: normal modes of a symmetric operator. The book's loose
  "λ bound the elements τ_ij" is true but the useful statement is about the normal stress on planes and the shear bound.
  Every fact fails for a non-symmetric tensor (complex λ, skew eigenvectors) — `principal_axes` refuses one.""",
  check="""Units: those of τ ✓. Limit τ = −pδ: all λ = −p, every direction principal, τ_s = 0 ✓. Number (Ex. 2.4, Γ = 1):
  λ = ±1, b¹ = (1, 1)/√2, b² = (−1, 1)/√2, C = 45° rotation, S' = diag(1, −1); on the 30° plane n·S·n = sin 60° = 0.866
  ∈ [−1, 1], τ_s = cos 60° = 0.5 ≤ 1 ✓ (`principal_axes`, `diagonalize`, `normal_stress_bounds` below). The sympy cell
  below re-runs the construction for a general symmetric 2×2.""",
  traps="""symmetry is used twice (steps 4 and 7) — one use gives realness, the other orthogonality. Repeated eigenvalues:
  orthogonal eigenvectors exist but are not unique (test with projectors, not vectors). det C = −1 from an unlucky sign
  choice: flip one b. $\\lambda^k$ is a label, not a power.""",
  check_src="""
# D17 verified symbolically for a general real symmetric 2×2 — the construction, not only the result (Ch. 1 P40, P61)
p, q, r, th = sp.symbols('p q r theta', real=True)             # τ = [[p, q], [q, r]] real symmetric; θ = the cut angle
tau_s = sp.Matrix([[p, q], [q, r]])                             # the symmetric tensor
lam_s = sp.symbols('lambda')                                    # the eigenvalue variable
poly = sp.expand((tau_s - lam_s * sp.eye(2)).det())             # characteristic polynomial λ² − (p + r)λ + (pr − q²)  (D18 step 6 in 2-D)
disc = sp.expand(poly.coeff(lam_s, 1)**2 - 4 * poly.coeff(lam_s, 0))   # its discriminant b² − 4c
assert sp.simplify(disc - ((p - r)**2 + 4 * q**2)) == 0         # step 5: the discriminant is (p − r)² + 4q² ≥ 0 → both roots real
lam1, lam2 = [sp.simplify(s) for s in sp.solve(poly, lam_s)]    # the two eigenvalues
b1 = sp.Matrix([q, lam1 - p]); b2 = sp.Matrix([q, lam2 - p])    # eigenvectors from the first row of (τ − λI)·b = 0
assert sp.simplify(b1.dot(b2)) == 0                             # step 8: orthogonal (uses pr − q² = λ1 λ2, Vieta P71)
C_s = sp.Matrix.hstack(b1 / sp.sqrt(b1.dot(b1)), b2 / sp.sqrt(b2.dot(b2)))   # step 9: unit eigenvectors as columns (sqrt(b·b), not .norm(): .norm() adds Abs() and blocks simplify)
D_s = sp.simplify(C_s.T * tau_s * C_s)                          # step 11: τ' = CᵀτC
assert sp.simplify(D_s[0, 1]) == 0 and sp.simplify(D_s[1, 0]) == 0   # the off-diagonals vanish: no shear in the principal axes
assert sp.simplify(D_s[0, 0] - lam1) == 0 and sp.simplify(D_s[1, 1] - lam2) == 0   # the diagonal holds the eigenvalues
n_s = sp.Matrix([sp.cos(th), sp.sin(th)])                       # step 12: any unit normal
c1, c2 = n_s.dot(C_s[:, 0]), n_s.dot(C_s[:, 1])                 # its weights on the principal axes
assert sp.simplify((n_s.T * tau_s * n_s)[0] - (lam1 * c1**2 + lam2 * c2**2)) == 0   # step 13: σ_n = Σ λ_k c_k²
assert sp.simplify(c1**2 + c2**2 - 1) == 0                      # the weights add to 1 → the bounds of step 14 follow
print("D17 verified symbolically: real λ, orthogonal b, CᵀτC diagonal, σ_n = Σ λ_k c_k² with Σ c_k² = 1")
""")
note("N62", "Ex. 2.4, stated (traced in the tiny example below)", r"""
⚠️ The book writes "$2S_{12} = du_1/dx_2 = \Gamma$" *and* $\mathbf S = [[0, \Gamma], [\Gamma, 0]]$. Those disagree by a
factor 2 (for $u_1(x_2)$ alone, $S_{12} = \tfrac12\,du_1/dx_2$). Every later line of the example uses $S_{12} = \Gamma$,
so we take **Γ ≡ S₁₂** (`ch02.example_2_4(Gamma)`); with $du_1/dx_2 = 2$ s⁻¹ as in C12, Γ = 1 s⁻¹. Also
"$\det|E_{ij} - \lambda\delta_{ij}|$" means $S_{ij}$. Interpretation: stretching at rate Γ along $\mathbf b^1$ (45°),
compression at −Γ along $\mathbf b^2$ (135°) — the S panel of C12's animation; Ch. 3 §3.4 makes this the definition of
strain rate. **Fig. 2.8** `N63` (the 45° axes) is our slider figure below.
""")
nb.worked_example("Ex. 2.4 by hand with Γ = 1 s⁻¹", r"""
1. $\mathbf S = \begin{bmatrix}0 & 1\\ 1 & 0\end{bmatrix}$.
2. $\det\begin{bmatrix}-\lambda & 1\\ 1 & -\lambda\end{bmatrix} = \lambda^2 - 1 = 0 \Rightarrow \lambda^1 = 1,
   \lambda^2 = -1$ (Vieta: sum 0 = trace ✓, product −1 = det ✓).
3. Eigenvector for λ = 1: the first row of $(\mathbf S - \lambda\mathbf I)\cdot\mathbf b = 0$ gives $-b_1 + b_2 = 0
   \Rightarrow b_1 = b_2$; normalise $b_1^2 + b_2^2 = 1 \Rightarrow \mathbf b^1 = (1, 1)/\sqrt2$.
4. For λ = −1: $b_1 + b_2 = 0 \Rightarrow \mathbf b^2 = (-1, 1)/\sqrt2$ (sign chosen so that $\det[\mathbf b^1\,
   \mathbf b^2] = +1$: a 45° rotation, not a mirror).
5. $\mathbf C = \tfrac1{\sqrt2}\begin{bmatrix}1 & -1\\ 1 & 1\end{bmatrix}$ — C02's 2-D matrix with θ = 45° ✓.
6. $S'_{12} = C_{i1}C_{j2}S_{ij} = C_{11}C_{22}S_{12} + C_{21}C_{12}S_{21} = \tfrac12 - \tfrac12 = 0$; $S'_{11} =
   2C_{11}C_{21}S_{12} = 1$; $S'_{22} = 2C_{12}C_{22}S_{12} = -1$: $\mathbf S' = \mathrm{diag}(1, -1)$ ✓.
7. Bounds: on the 30° plane $\mathbf n\cdot\mathbf S\cdot\mathbf n = \sin60° = 0.866 \in [-1, 1]$ ✓; max shear
   $(1 - (-1))/2 = 1$, reached at φ = 0° and 90° (the original axes, where S has only shear) ✓.
""")
nb.code("""
res = ch02.example_2_4(1.0)                                          # Ex. 2.4 packaged: Γ = S₁₂ = 1 s⁻¹, the book's order (Γ, −Γ), b¹ first
print(res["lam"], np.round(res["B"], 4), np.round(np.rad2deg(res["angle_rad"]), 1))   # [1 −1], columns (1,1)/√2 and (−1,1)/√2, 45.0°
print(np.round(res["S_prime"], 12))                                  # diag(1, −1): no shear in the principal axes
lam, Bm = ch02.principal_axes(A)                                     # C07's 3×3 example: λ ascending, unit eigenvectors as columns, det +1
print(np.round(lam, 4), np.round(np.linalg.det(Bm), 12))             # [1.2679 3. 4.7321] 1.0
Cp, tp = ch02.diagonalize(A)                                         # (C, τ' = CᵀτC) — fact 3
print(np.round(tp, 10))                                              # diag(1.2679, 3, 4.7321)
print(np.round(np.roots(ch02.characteristic_polynomial(A)), 4))      # the same three λ from the cubic (D18): they are invariants
n_mc = 1000 if not FAST else 300                                     # random unit normals (FAST: 300)
print(ch02.normal_stress_bounds(A, np.random.default_rng(6), n=n_mc))   # (min σ_n, max σ_n, max τ_s) over random planes: inside [1.27, 4.73], shear ≤ 1.73
print(ch02.principal_angle_2d([[0., 1.], [1., 0.]]))                 # 0.7854 rad = 45°: the 2-D shortcut
try:                                                                 # a non-symmetric input …
    ch02.principal_axes([[0., 1.], [0., 0.]])
except ValueError as e:                                              # … is refused on purpose: facts (1)–(4) do not hold for it
    print("refused:", e)
""", explain="""
1. `example_2_4` returns the book's numbers: λ = (Γ, −Γ), the 45° axes, S' = diag(1, −1).
2. `principal_axes` on C07's 3×3 example: three real eigenvalues, an orthonormal C with det +1.
3. `diagonalize` shows fact 3: $\\mathbf C^{\\rm T}\\mathbf A\\mathbf C$ is diagonal.
4. The roots of D18's characteristic polynomial are the same λ — eigenvalues are invariants.
5. `normal_stress_bounds` tries 1000 random planes: no normal stress leaves [λ_min, λ_max] and no shear exceeds
   (λ_max − λ_min)/2 = 1.732 (fact 4).
6. The 2-D principal angle by the `arctan2` shortcut (P70).
7. Non-symmetric input is refused — the library raises on purpose (a `try/except ValueError` catches the message).
""")
nb.check_agree("""
S2 = np.array([[0., 1.], [1., 0.]])                                              # Ex. 2.4's S
tr, det = S2.trace(), np.linalg.det(S2)                                          # I₁ = 0 and I₂ = det = −1 (2-D invariants, D18)
lam_mine = np.array([(tr + np.sqrt(tr**2 - 4 * det)) / 2, (tr - np.sqrt(tr**2 - 4 * det)) / 2])   # the quadratic formula for λ² − I₁λ + I₂ = 0
b1 = np.array([S2[0, 1], lam_mine[0] - S2[0, 0]]); b1 /= np.linalg.norm(b1)       # eigenvector from the first row of (S − λI)b = 0, normalised (P67)
b2 = np.array([-b1[1], b1[0]])                                                   # the perpendicular unit vector (det [b¹ b²] = +1)
assert np.allclose(lam_mine, res["lam"])                                         # (1, −1)
assert np.allclose(abs(b1 @ res["B"][:, 0]), 1)                                  # the same first axis (up to sign)
assert np.allclose(np.c_[b1, b2].T @ S2 @ np.c_[b1, b2], np.diag(lam_mine))      # CᵀSC = diag(λ): fact 3 by hand (np.c_ stacks columns)
print("quadratic formula + first row of (S − λI)b = 0 reproduce example_2_4")
""")
nb.plotly("""
def strain_axes(Gamma):                                                       # one slider position: S = [[0, Γ], [Γ, 0]]
    S_ = np.array([[0., Gamma], [Gamma, 0.]])                                 # the strain rate [1/s]
    ex = ch02.example_2_4(Gamma)                                              # its principal axes and values
    t_ = 0.5                                                                  # a fixed deformation time [s]
    sq_0 = ch02.deform_square(S_, 0.0, 30, half_width=0.7, boundary_only=True)   # the square outline at t = 0
    sq_t = ch02.deform_square(S_, t_, 30, half_width=0.7, boundary_only=True)    # … after stretching along b¹ by e^{Γt}, squeezing along b² by e^{−Γt}
    b1_, b2_ = ex["B"][:, 0], ex["B"][:, 1]                                   # the eigen-axes
    return {"square at t = 0": (sq_0[0], sq_0[1]),
            "square under S at t = 0.5 s": (sq_t[0], sq_t[1]),
            "b¹ (stretch, λ¹ = +Γ)": ([0, ex["lam"][0] * b1_[0]], [0, ex["lam"][0] * b1_[1]]),    # arrow of length λ¹ along b¹
            "b² (squeeze, λ² = −Γ)": ([0, -ex["lam"][1] * b2_[0]], [0, -ex["lam"][1] * b2_[1]])}  # arrow of length |λ²| along b²
fig = slider_figure(strain_axes, "Γ", np.linspace(0.2, 3.0, 30 if not FAST else 15), unit="1/s", xlabel="x₁ [m]", ylabel="x₂ [m]",
                    title="Ex. 2.4: the principal axes stay at 45°; only the stretch (and λ = ±Γ) grows with Γ",
                    xrange=[-3.3, 3.3], yrange=[-3.3, 3.3], height=520, active=8)
fig.update_yaxes(scaleanchor="x", scaleratio=1)                               # equal axes
fig.show()
""", explain="""
1. For each Γ the strain rate S = [[0, Γ], [Γ, 0]] is diagonalised by `example_2_4`; the two arrows are the eigen-axes
   scaled by |λ| = Γ.
2. A square of material points is deformed by S alone for 0.5 s (`deform_square`, P79): it stretches along b¹ by
   $e^{\\Gamma t}$ and squeezes along b² by $e^{-\\Gamma t}$, never shearing those two axes.
""")
see_read_change("A square, its stretched-and-squeezed image, and two arrows at ±45° whose length grows with Γ; the arrows' directions never move.",
                "The principal *directions* depend on the shape of S (here always pure shear → 45°), the eigenvalues on its size (±Γ).",
                "…you added $S_{11} = 1$: the angle drops below 45° (`ch02.principal_angle_2d([[1, Γ], [Γ, 0]])`) and Mohr's centre moves right by ½.")
nb.md(r"""
🔁 *Reminder (Ch. 1 primer P47):* `live(fn, a=(lo, hi, step))` builds ipywidgets sliders that re-run `fn` while a kernel
runs; the published page freezes it, so the slider figure above carries the idea there.
""")
nb.live("""
from fluidpy.core.interact import live                                              # ipywidgets sliders (kernel only, Ch. 1 P47)
def principal_plot(S11=0.0, S22=0.0, S12=1.0):                                      # called again whenever a slider is released
    S_ = np.array([[S11, S12], [S12, S22]])                                         # a symmetric 2×2 [1/s or Pa]
    lam_, B_ = ch02.principal_axes(S_)                                              # eigenvalues (ascending) and unit eigenvectors as columns
    ang_ = np.rad2deg(ch02.principal_angle_2d(S_))                                  # the first principal angle [°]
    c_, r_ = ch02.mohr_circle_2d(S_)                                                # Mohr's circle
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4))                              # element + axes | Mohr circle
    sq_ = ch02.deform_square(S_, 0.4, 30, half_width=0.7, boundary_only=True)      # the element after 0.4 s under S
    a1.plot(*ch02.deform_square(S_, 0.0, 30, half_width=0.7, boundary_only=True), "--", color=COLORS["muted"]); a1.plot(*sq_, color=COLORS["accent"])
    for k, col in ((0, COLORS["rose"]), (1, COLORS["blue"])):                       # eigen-axes: smallest λ rose, largest blue
        a1.annotate("", 1.3 * B_[:, k], -1.3 * B_[:, k], arrowprops=dict(arrowstyle="<->", color=col, lw=2))
    a1.set_xlim(-2.5, 2.5); a1.set_ylim(-2.5, 2.5); a1.set_aspect("equal"); a1.set_title(f"λ = ({lam_[0]:.2f}, {lam_[1]:.2f}), first axis at {ang_:.1f}°", fontsize=9)
    th_ = np.linspace(0, 2 * np.pi, 181); a2.plot(c_ + r_ * np.cos(th_), r_ * np.sin(th_), color=COLORS["accent"])   # the circle
    a2.plot(lam_, [0, 0], "o", color=COLORS["ink"]); a2.axhline(0, color=COLORS["muted"], lw=0.8); a2.set_aspect("equal")
    a2.set_xlabel("σ_n"); a2.set_ylabel("τ_s"); a2.set_title(f"Mohr: centre {c_:.2f}, radius {r_:.2f} (max shear)", fontsize=9)
    plt.show()                                                                      # draw for these slider values
_ = live(principal_plot, S11=(-2.0, 2.0, 0.1), S22=(-2.0, 2.0, 0.1), S12=(-2.0, 2.0, 0.1))   # three sliders (_ hides the widget's return value)
""", explain="""
1. The callback diagonalises the symmetric 2×2 you set with the sliders, draws the deformed element with its eigen-axes,
   and Mohr's circle with the two principal values on the σ axis.
2. Set S₁₂ = 0: the axes align with x₁, x₂. Set S₁₁ = S₂₂ and S₁₂ = 0: the circle shrinks to a point — every plane is
   principal (the hydrostatic case).
""")
nb.md(r"""
**What would change if…** τ were not symmetric? Complex eigenvalues and skewed eigenvectors are possible and the bounds
fail — which is why `principal_axes` refuses, and why Ch. 4's proof that the stress *is* symmetric matters.
""")

# =====================================================================================================================
# A.10 §2.12 Gauss' Theorem — C14 (D25, E4), C15 (D21, D22)
# =====================================================================================================================
nb.section("2.12", "Gauss' Theorem", intro="""
**What is this section about?** A derivative integrated over a volume equals the field itself integrated over the
boundary with the outward normal — for a field of any order. Read backwards on a tiny box it *defines* the divergence
as outflux per unit volume, and the Cartesian formula (2.23) falls out face by face.
""")
core("C14", "Gauss' theorem: derivative inside = normal outside", r"""
Why does adding up $\partial Q/\partial x_i$ over a whole volume only depend on what Q does on the surface?
""")
nb.md(r"""
#### The problem in plain words

Ch. 4 will write conservation laws for a blob of fluid: the mass inside changes only by what flows through its skin;
the momentum inside changes by the forces on its skin (C05!). Turning "inside" into "on the skin" — and back — is
Gauss' theorem. Without it there is no differential form of any conservation law.
""")
nb.md(r"""
#### The idea

```
1-D:  ∫_a^b f'(x) dx = f(b) − f(a)                 (fundamental theorem: the inside adds up to the two ends)
3-D:  ∭_V ∂Q/∂x_i dV = ∯_A n_i Q dA                 (the ends become the boundary, with n_i = ±1 on two faces of a box)
any V:  tile it with boxes — every interior face is counted twice with opposite n and cancels
```
""")
P("P83", "volume and surface integrals as midpoint sums", r"""
$\iiint_V f\,dV$ adds f × (small volume) over a 3-D grid of cells; the midpoint rule uses f at each cell centre (error
∝ h²). A surface integral $\iint_A g\,dA$ does the same on a 2-D grid of the surface. Iterated integrals: integrate in
x first, then y, then z. (🔁 Ch. 1 P27 defined the definite integral, P37 the trapezoid rule; the midpoint rule is its
sibling.) For a *linear* integrand the midpoint rule is exact — the fact D22 leans on.
""", code="""
n = 20; xs = (np.arange(n) + 0.5) / n                      # cell centres on [0, 1]
X, Y, Z = np.meshgrid(xs, xs, xs, indexing='ij')           # a 20³ grid of cell centres
print(np.sum(X**2) / n**3)                                 # ∭ x² dV over the unit cube = 1/3 (0.3329 with n = 20)
print(np.sum(X[:, :, 0]**2) / n**2)                        # ∬ x² dA over the unit square = 1/3
""")
P("P84", "fundamental theorem of calculus", r"""
$\int_a^b f'(x)\,dx = f(b) - f(a)$: integrating a derivative gives the change of the function between the ends (🔁 Ch. 1
P27 defined the integral; this is its partner). It is the whole content of Gauss' theorem in one dimension.
""", code="""
f = lambda x: x**3                                         # a function
xs = np.linspace(0, 2, 2001)                               # a fine grid on [0, 2]
print(np.trapezoid(3 * xs**2, xs), f(2) - f(0))            # 8.000 8: ∫ f' dx = f(2) − f(0)  (trapezoid rule, Ch. 1 P37)
""")
nb.md(r"""
#### The maths, step by step

1. **Statement (2.30):** $\iiint_V \partial Q/\partial x_i\,dV = \iint_A n_i Q\,dA$ for Q of any order (scalar, vector,
   tensor — the other indices ride along), n outward.
2. For a vector, set Q → $Q_i$ and sum on i (allowed: (2.30) is linear in Q; N64): $\iiint_V \nabla\cdot\mathbf Q\,dV =
   \iint_A \mathbf n\cdot\mathbf Q\,dA$ — the **divergence theorem**: volume integral of the divergence = net outflux.
3. **Proof (D25):** on a box by the fundamental theorem of calculus along one axis; on any V by tiling and cancelling
   interior faces.
""")
D("D25", "Gauss' theorem, from the fundamental theorem of calculus", ref="2.30",
  goal="""Prove the theorem the book only states: the volume integral of a derivative equals the surface integral of the
  field with the outward normal — first for a box, then for any volume by tiling.""",
  assumptions="""Q continuously differentiable on V and its boundary (steps 2, 9) · the boundary A piecewise smooth, so
  that the tiling's outer faces approach it (step 9).""",
  start=(r"\int_a^b f'(x)\,dx = f(b) - f(a)", "the fundamental theorem of calculus (P84): integrating a derivative gives the values at the two ends."),
  plan=["On a box, write the volume integral of ∂Q/∂x₁ as an iterated integral and do the x₁ integral first.",
        "Read the right side of (2.30) face by face: n₁ = ±1 on two faces, 0 on four.", "Repeat for i = 2, 3.",
        "Tile a general volume with boxes; interior faces cancel."],
  uses=["fundamental theorem of calculus (P84)", "iterated integrals (P83)", "outward normals of a box (C04's faces)",
        "additivity of integrals over disjoint regions (step 7)", "vector area / opposite normals (P69; step 8)"],
  steps=[
      ("Write the i = 1 volume integral on a box as iterated integrals",
       r"\iiint_V \frac{\partial Q}{\partial x_1}\,dV = \int_{a_3}^{b_3}\!\int_{a_2}^{b_2}\Big[\int_{a_1}^{b_1}\frac{\partial Q}{\partial x_1}\,dx_1\Big]dx_2\,dx_3",
       r"On a box $[a_1, b_1]\times[a_2, b_2]\times[a_3, b_3]$ a volume integral is three nested one-dimensional integrals in any order; we do x₁ innermost because the derivative is in x₁.",
       "Integrate along x₁ first, on every line of the box."),
      ("Apply the fundamental theorem to the inner integral",
       r"\int_{a_1}^{b_1}\frac{\partial Q}{\partial x_1}\,dx_1 = Q(b_1, x_2, x_3) - Q(a_1, x_2, x_3)",
       "Along a line with x₂, x₃ fixed, ∂Q/∂x₁ is an ordinary derivative of a function of x₁; the fundamental theorem gives the two end values (Q ∈ C¹).",
       "Each line contributes Q at its far end minus Q at its near end."),
      ("Insert and split into two face integrals",
       r"\iiint_V \frac{\partial Q}{\partial x_1}\,dV = \iint Q(b_1, x_2, x_3)\,dx_2dx_3 - \iint Q(a_1, x_2, x_3)\,dx_2dx_3",
       r"Step 2 into step 1; the remaining double integrals run over the face $[a_2, b_2]\times[a_3, b_3]$ — the two faces perpendicular to x₁.",
       "The inside adds up to \"Q on the far face minus Q on the near face\"."),
      ("Read the right side of (2.30) face by face",
       r"n_1 = +1\ \text{on}\ x_1 = b_1,\quad n_1 = -1\ \text{on}\ x_1 = a_1,\quad n_1 = 0\ \text{on the other four faces}",
       "The outward normal of a box face is ±e_i; its first component is nonzero only on the two faces perpendicular to x₁. Four of the six faces contribute nothing to the i = 1 term.",
       "Only the two faces facing ±x₁ count for the x₁-derivative."),
      ("Conclude (2.30) for i = 1 on the box",
       r"\iint_A n_1Q\,dA = \iint Q(b_1,\cdot)\,dA - \iint Q(a_1,\cdot)\,dA = \iiint_V\frac{\partial Q}{\partial x_1}\,dV",
       "Step 4 makes the surface integral exactly the two face integrals of step 3 with the signs of n₁. Left and right sides of (2.30) coincide for i = 1.",
       "Gauss' theorem holds for the x₁-derivative on a box."),
      ("Repeat for i = 2 and i = 3",
       r"\iiint_V\frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA\quad\text{on any box, } i = 1, 2, 3",
       "Steps 1–5 with the roles of the axes relabelled (integrate along x_i first; the faces ⊥ x_i carry n_i = ±1). Nothing used which axis was first.",
       "The same for each of the three directions."),
      ("Tile a general volume with small boxes and add",
       r"\sum_k\iiint_{V_k}\frac{\partial Q}{\partial x_i}\,dV = \sum_k\iint_{A_k} n_iQ\,dA,\qquad \sum_k\iiint_{V_k} = \iiint_{V_{\rm tiles}}",
       "Step 6 holds on each box V_k; adding the equations adds both sides. Integrals over disjoint regions add up to the integral over their union.",
       "Gauss on every tile, summed."),
      ("Cancel the interior faces",
       r"\sum_k\iint_{A_k} n_iQ\,dA = \iint_{\text{outer faces}} n_iQ\,dA",
       "A face shared by two tiles appears twice with opposite outward normals (n and −n) and the same values of Q (continuity), so the two contributions sum to zero. Only faces on the outside of the tiling survive.",
       "Every inner wall is counted once from each side and cancels."),
      ("Refine the tiling",
       r"\iiint_V\frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA",
       "As the boxes shrink, the tiled region fills V and its staircase boundary tends to A with the outer normals tending to n (A piecewise smooth, Q continuous), so both sides converge to the integrals over V and A. Q may carry any other indices: nothing depended on its order.",
       "Gauss' theorem for any volume and any field."),
  ],
  result=(r"\iiint_V\frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA \quad\text{for Q of any order}",
          "a derivative integrated over the inside equals the field weighted by the outward normal on the boundary."),
  interpret="""Every conservation law of Ch. 4 passes through this line: the rate of change of what is inside a volume
  equals the flux through its boundary, and shrinking the volume (D21) turns the integral law into a partial
  differential equation. The theorem needs Q differentiable *inside* — a point source or a shock breaks it and adds a
  surface term.""",
  check="""Units: [Q]·m² both sides ✓. Limit Q constant: left 0; right ∮ n_i dA = 0 (vector area of a closed surface, P69)
  ✓. Number: Q = (x, y, z) on the unit cube: ∭ 3 dV = 3; faces +1 + 0 + 1 + 0 + 1 + 0 = 3 ✓ (`divergence_theorem_box`
  below); sphere benchmark F = (2x, y², z²): 8π/3 both sides (`divergence_theorem_sphere`); `divergence_theorem_tiled`
  reports interior faces summing to 0 to round-off.""",
  traps="""thinking all six faces contribute to each i. Cancelling interior faces without opposite normals and equal Q.
  Forgetting that the tiles' boundary is only an approximation of A until the limit.""")
note("N64", "The divergence theorem and the outflux reading", r"""
$\mathbf n\cdot\mathbf Q\,dA$ is how much Q leaves through dA. Benchmark (the standard divergence-theorem example):
$\mathbf F = (2x, y^2, z^2)$ through the unit sphere: $\nabla\cdot\mathbf F = 2 + 2y + 2z$, whose integral over the ball
is $2\cdot\tfrac43\pi + 0 + 0 = 8\pi/3 = 8.378$; `ch02.divergence_theorem_sphere` gives both sides. **Fig. 2.9** `N65`
is our 3-D box below, faces coloured by n·Q.
""", equation=r"\iiint_V \frac{\partial Q_i}{\partial x_i}\,dV = \iint_A n_i Q_i\,dA \quad\text{or}\quad \iiint_V \nabla\cdot\mathbf Q\,dV = \iint_A \mathbf n\cdot\mathbf Q\,dA", ref="§2.12")
nb.worked_example("Q = (x, y, z) on the unit cube [0, 1]³", r"""
1. Inside: $\nabla\cdot\mathbf Q = 1 + 1 + 1 = 3$; volume 1 → left side 3.
2. Outside, face by face: on x = 1, n = +e₁, n·Q = x = 1, area 1 → +1; on x = 0, n·Q = −x = 0 → 0; same for y and z →
   right side 1 + 0 + 1 + 0 + 1 + 0 = 3 ✓.
3. Q = (x², 0, 0): inside $\int_0^1 2x\,dx = 1$; outside +1 (x = 1) + 0 = 1 ✓ — only the two x-faces contribute, as D25
   says.
""")
nb.code("""
Q = lambda x, y, z: np.stack([x, y, z])                                  # the vector field Q = (x, y, z) as a function of grid arrays (components on axis 0)
lhs, rhs = ch02.divergence_theorem_box(Q, ((0, 1),) * 3, 16)             # both sides of (2.30) on the unit cube by midpoint sums (P83), 16 nodes per side
print(lhs, rhs)                                                          # 3.0 3.0
Q2 = lambda x, y, z: np.stack([x**2, 0 * y, 0 * z])                      # Q = (x², 0, 0)
print(ch02.divergence_theorem_box(Q2, ((0, 1),) * 3, 16))                # (1.0, 1.0): only the two x-faces leak
F = lambda x, y, z: np.stack([2 * x, y**2, z**2])                        # the sphere benchmark field
lhs_s, rhs_s = ch02.divergence_theorem_sphere(F, 1.0, 48 if not FAST else 24)   # ∭ ∇·F dV over the unit ball vs ∯ n·F dA over the sphere
print(lhs_s, rhs_s, 8 * np.pi / 3)                                       # 8.378 8.378 8.378 (quadrature error ~1e-4)
print(ch02.gauss_gradient_box(lambda x, y, z: x * y * z, ((0, 1),) * 3, 16))   # (2.30) on a SCALAR: both sides are the vector (¼, ¼, ¼)
Q2d = ch02.radial_field(1.0, dim=2)                                      # the plane field (x, y) for the 2-D pictures
print(ch02.flux_through_faces(Q2d, ((0, 1), (0, 1)), 16))                # {'+x': 1, '-x': 0, '+y': 1, '-y': 0}: the right and top faces leak 1 each
print(ch02.divergence_theorem_tiled(Q2d, ((0, 1), (0, 1)), 4, 16))       # 4×4 tiles: (sum over tiles, outer boundary, interior faces) = (2, 2, 0): D25 step 8 as a number
""", explain="""
1. `divergence_theorem_box` evaluates both sides of (2.30) for a vector Q on a box by midpoint sums: 3 = 3.
2. For Q = (x², 0, 0) both sides are 1: only the two x-faces contribute (D25 step 4).
3. The sphere benchmark: both sides equal 8π/3 to quadrature error.
4. `gauss_gradient_box` is the scalar form of (2.30): the volume integral of ∇(xyz) equals the surface integral of n·(xyz)
   — a vector on both sides.
5. `flux_through_faces` lists the four signed face fluxes of a plane field.
6. `divergence_theorem_tiled` splits the square into 4×4 tiles: the tile boundaries sum to the outer flux and the interior
   faces cancel to round-off — D25 step 8 made visible.
""")
nb.check_agree("""
n = 16; c = (np.arange(n) + 0.5) / n                                    # face-cell centres (midpoint rule, P83)
Ai, Bi = np.meshgrid(c, c, indexing='ij')                               # a 16×16 grid on one face
flux = 0.0                                                              # ∯ n·Q dA from scratch, face by face
for axis in range(3):                                                   # the pair of faces perpendicular to axis
    for side, sgn in ((1.0, +1), (0.0, -1)):                            # far face n = +e_axis, near face n = −e_axis
        coords = [Ai, Bi]; coords.insert(axis, np.full_like(Ai, side))  # the face's points: fixed coordinate = side
        flux += sgn * Q(*coords)[axis].sum() / n**2                     # n·Q = ±Q_axis, times the face area 1/n² per cell
Xc, Yc, Zc = np.meshgrid(c, c, c, indexing='ij')                        # cell centres inside
vol = (3 * np.ones_like(Xc)).sum() / n**3                               # ∭ ∇·Q dV with ∇·Q = 3
assert np.allclose((vol, flux), ch02.divergence_theorem_box(Q, ((0, 1),) * 3, n))   # the library does the same sums
print(vol, flux, "= 3: six face sums by hand")
""")
nb.plotly("""
from scripts.ch02_fig2_9_gauss import gauss_box_figure                   # drawing helper (P64): faces coloured by n·Q
fig = gauss_box_figure(F)                                                # F = (2x, y², z²) on the unit cube; the two totals in the title
fig.show()                                                               # rotate: blue = inflow, orange = outflow
""", explain="""
1. Each face of the unit cube is coloured by the local outflux n·F (blue in, orange out) and labelled with its total.
2. The title compares ∭ ∇·F dV with the sum of the six face fluxes — the two sides of (2.30).
""")
see_read_change("Orange faces where F leaves (x = 1, y = 1, z = 1), blue or pale where it enters or is tangent; the title's two totals agree.",
                "The sum of the face colours (weighted by area) equals the integral of the divergence inside — 4 = 4 for this F on the cube.",
                "…Q = b × x: every face is half blue, half orange, net zero — solenoidal (`gauss_box_figure(lambda x, y, z: np.stack([-y, x, 0*z]))`).")
nb.plotly("""
Fs = ch02.smooth_test_field(3)                                           # a sin/cos field with a known exact divergence
exact_vol = ch02.volume_integral_box(Fs.div_fn, ((0, 1),) * 3, 96)       # the exact ∭ ∇·Q dV (very fine grid)
ns_mesh = [4, 6, 8, 12, 16, 24, 32] if not FAST else [4, 8, 16, 32]      # mesh counts per side
def two_sides(n_):                                                       # one slider position: n nodes per direction
    l, r_ = ch02.divergence_theorem_box(Fs, ((0, 1),) * 3, int(n_))      # both sides by midpoint sums
    return {"∭ ∇·Q dV (inside)": ([0, 0], [0, l]), "∯ n·Q dA (boundary)": ([1, 1], [0, r_]),   # two bars
            "exact value": ([-0.5, 1.5], [exact_vol, exact_vol])}         # the reference line
fig = slider_figure(two_sides, "n", ns_mesh, unit="nodes/side", xlabel="", ylabel="value [Q·m²]",
                    title="Both sides are approximations of one number: the gap closes as the mesh refines (∝ 1/n²)",
                    xrange=[-0.5, 1.5], height=420, modes={"∭ ∇·Q dV (inside)": "lines+markers", "∯ n·Q dA (boundary)": "lines+markers"})
fig.update_xaxes(tickvals=[0, 1], ticktext=["inside", "boundary"]); fig.show()
gaps = [abs(np.subtract(*ch02.divergence_theorem_box(Fs, ((0, 1),) * 3, n_))) for n_ in ns_mesh]   # |lhs − rhs| per mesh
print({n_: f"{gp:.1e}" for n_, gp in zip(ns_mesh, gaps)})              # the gap falls by about 4 per doubling of n
""", explain="""
1. For the smooth sin/cos field the two sides of (2.30) are computed by midpoint sums with n nodes per side.
2. Two "bars" (vertical lines) and the exact value; the printed gaps shrink by ≈ 4 per doubling — both sides are
   second-order approximations of one number.
""")
see_read_change("Two bars of almost equal height and a dashed exact line; dragging n up brings the bars onto the line.",
                "The two bars meet as the mesh refines; the printed gap falls by about 4 per doubling of n.",
                "…Q were a polynomial of degree ≤ 2: the gap would be 1e-16 at every n (midpoint sums are exact for quadratics on each face).")
nb.explainer("gauss_flux_box", heading="What leaks out of a box?",
             why="""The theorem has a left side that lives inside and a right side on the boundary; moving and resizing a box
             over a field, watching the face fluxes re-sum to the integrated divergence, tiling it to see interior faces
             cancel, and shrinking it until (1/V)∮ settles on ∇·Q is the whole (2.30) → (2.32) chain as an experiment.""",
             tries=["Drag the box over the source: the total bar turns orange.",
                    "Tile 4×4: interior arrows appear in opposite pairs and cancel — D25 step 8.",
                    "Shrink h on the log slider and watch (1/V)∮ approach the dashed ∇·Q(x₀) with slope 2.",
                    "b × x preset: every box, anywhere, has zero net flux."])
nb.md(r"""
**What would change if…** the volume shrank to a point? Divide both sides by V and take the limit: the right side
becomes "outflux per unit volume" — a definition of the divergence that never mentions coordinates. C15.
""")

# ---------------------------------------------------------------------------------------------------------------------
core("C15", "Divergence as outflux per unit volume: the integral definitions", r"""
Can $\nabla\cdot\mathbf Q$ be defined without any coordinates at all — and how does the formula $\partial Q_i/\partial
x_i$ come back out?
""")
nb.md(r"""
#### The problem in plain words

The Cartesian formula (2.23) is a recipe, not a meaning. The meaning is: put a tiny closed surface around the point,
measure what leaks out, divide by the volume. That definition works in spherical coordinates on a planet as well as on
a Cartesian grid, and it is exactly the box argument Ch. 4 uses to derive the continuity equation.
""")
nb.md(r"""
#### The idea

```
(2.30) on a tiny box around x₀:  (1/V) ∯ n·Q dA  →  ∇·Q(x₀)  as V → 0          (2.32)
face by face (Ex. 2.5):  Q on the +e₁ face ≈ Q + (Δx₁/2) ∂Q/∂x₁ ;  on the −e₁ face ≈ Q − (Δx₁/2) ∂Q/∂x₁
                         opposite faces: the Q's cancel, the slopes add  →  ∂Q₁/∂x₁ Δx₁Δx₂Δx₃ ;  divide by V
```
""")
P("P85", "mean-value theorem for integrals", r"""
If f is continuous on a region V, then $\iiint_V f\,dV = f(\mathbf x^*)\,V$ for some point x* inside V: the integral
equals the volume times some intermediate value of f. As V shrinks to a point x₀, x* is squeezed onto x₀.
""", code="""
xs = np.linspace(0.9, 1.1, 2001)                          # a small interval around 1
f = xs**2                                                 # f = x²
mean = np.trapezoid(f, xs) / 0.2                          # the integral divided by the length = f(x*) for some x* (Ch. 1 P37)
print(mean, np.sqrt(mean))                                # 1.0033, and x* = 1.0017 lies inside [0.9, 1.1]
""")
nb.md(r"""
#### The maths, step by step

1. (2.31) $\mathcal DQ = \lim_{V\to0}(1/V)\iint_A n_iQ\,dA$ — the generalised derivative (gradient of any order), from
   (2.30) by the mean-value theorem (D21).
2. Contract: (2.32) $\nabla\cdot\mathbf Q = \lim (1/V)\iint \mathbf n\cdot\mathbf Q\,dA$; cross: (2.33) $\nabla\times
   \mathbf Q = \lim (1/V)\iint \mathbf n\times\mathbf Q\,dA$.
3. Ex. 2.5 (D22): evaluate (2.32) on a box — Taylor to the six face centres, opposite faces cancel the zeroth order,
   divide by the volume, take the limit: (2.23) returns.
""")
D("D21", "The integral definitions as small-volume limits of Gauss' theorem", ref="2.31–2.33",
  goal="""Turn Gauss' theorem, read on a shrinking volume, into coordinate-free definitions of the gradient, the
  divergence and the curl — the book calls (2.31) "the limiting form" and leaves the limit implicit.""",
  assumptions="""∂Q/∂x_i continuous near x₀ (steps 2 and 4) · V shrinks to x₀ in every direction (step 4; then the shape
  does not matter).""",
  start=(r"\iiint_V \frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA \quad\text{on a small V around } \mathbf x_0",
         "Gauss' theorem (2.30, D25) for any Q, applied to a small region."),
  plan=["Replace the volume integral by \"value at some interior point × volume\" (mean-value theorem).", "Divide by V.",
        "Shrink V to x₀.", "Contract with e_i for the divergence, cross with ε for the curl."],
  uses=["Gauss' theorem (2.30) (D25)", "mean-value theorem for integrals (P85)", "limits and continuity (P68)",
        "the index form of ∇·, (2.23) (C10)", "the index form of ∇×, (2.24) (C11)"],
  steps=[
      ("Apply Gauss' theorem to a small V around x₀",
       r"\iiint_V \frac{\partial Q}{\partial x_i}\,dV = \iint_A n_iQ\,dA",
       "(2.30) holds for any volume with a piecewise smooth boundary; we pick a small one around the point where we want a derivative.",
       "The inside sum of the derivative equals the boundary sum of Q."),
      ("Replace the volume integral by an interior value times V",
       r"\iiint_V \frac{\partial Q}{\partial x_i}\,dV = \frac{\partial Q}{\partial x_i}(\mathbf x^*)\;V",
       "Mean-value theorem for integrals: a continuous integrand over V equals its value at some point x* in V times the volume (P85). We do this to turn an integral into a single value.",
       "The integral is the volume times the derivative somewhere inside."),
      ("Divide by V",
       r"\frac{\partial Q}{\partial x_i}(\mathbf x^*) = \frac1V\iint_A n_iQ\,dA",
       "V > 0. The right side is now \"boundary flux of Q per unit volume\", a ratio that stays finite as V shrinks.",
       "The derivative at an interior point equals the flux per unit volume."),
      ("Shrink V to the point",
       r"\frac{\partial Q}{\partial x_i}(\mathbf x_0) = \lim_{V\to0}\frac1V\iint_A n_iQ\,dA \equiv \mathcal D Q",
       "x* is trapped in V, so x* → x₀; continuity of ∂Q/∂x_i lets the left side tend to its value at x₀. The left side never knew the shape of V, so neither does the limit. This is (2.31).",
       "The derivative of Q is the limit of the boundary flux per unit volume — no coordinates needed."),
      ("Specialise to a vector Q and contract",
       r"\frac{\partial Q_i}{\partial x_i} = \lim_{V\to0}\frac1V\iint_A n_iQ_i\,dA",
       "Put Q → Q_j in step 4 (Gauss holds component by component, the theorem being linear in Q), then set j = i and sum — contraction on both sides.",
       "Adding the three diagonal derivatives adds the three flux components."),
      ("Recognise the divergence and the outflux",
       r"\nabla\cdot\mathbf Q = \lim_{V\to0}\frac1V\iint_A \mathbf n\cdot\mathbf Q\,dA",
       r"$\partial Q_i/\partial x_i$ is (2.23) and $n_iQ_i = \mathbf n\cdot\mathbf Q$ is the outward flux density. This is (2.32).",
       "Divergence = what leaks out per unit volume."),
      ("Cross instead of dot",
       r"(\nabla\times\mathbf Q)_k = \varepsilon_{kij}\frac{\partial Q_j}{\partial x_i} = \lim_{V\to0}\frac1V\iint_A (\mathbf n\times\mathbf Q)_k\,dA",
       r"Multiply step 4 (with Q → Q_j) by $\varepsilon_{kij}$ and sum on i, j — a fixed linear combination passes through the limit; $\varepsilon_{kij}\partial_iQ_j$ is (2.24) and $\varepsilon_{kij}n_iQ_j = (\mathbf n\times\mathbf Q)_k$ by (2.21). This is (2.33).",
       "Curl = the boundary's \"n × Q\" per unit volume."),
  ],
  result=(r"\mathcal DQ = \lim_{V\to0}\frac1V\iint n_iQ\,dA\ (2.31),\quad \nabla\cdot\mathbf Q = \lim\frac1V\iint \mathbf n\cdot\mathbf Q\,dA\ (2.32),\quad \nabla\times\mathbf Q = \lim\frac1V\iint \mathbf n\times\mathbf Q\,dA\ (2.33)",
          "gradient, divergence and curl are boundary sums per unit volume in the limit of a vanishing volume."),
  interpret="""The operators of §2.9 are properties of the field, not of Cartesian axes: the same recipe in spherical or
  cylindrical coordinates gives the curvilinear formulas of Appendix B (N72), and Ch. 4 derives continuity from exactly
  (2.32) on a fluid box. The limit exists only where Q is smooth; at a point source (the explainer's preset) the flux is
  fixed while V → 0 and (1/V)∮ diverges — a delta function, not a derivative.""",
  check="""Units: (2.32) — [Q]·m²/m³ = [Q]/m on both sides ✓. Limit: Q constant: the boundary integral of n vanishes
  (vector area of a closed surface, P69), so all three derivatives are 0 ✓. Number: Q = (x², 0, 0), x₀ = (1, 0, 0), cube
  h = 0.2: (1/V)∮ = (1.21 − 0.81)(0.04)/0.008 = 2.000 = 2x₀ exactly ✓; Q = (x³, 0, 0): 3.01 vs 3 at h = 0.2, error h²/4,
  order 2 (`integral_divergence`, `integral_definition_convergence("divergence")` below).""",
  traps="""forgetting the 1/V. Assuming shape-independence without continuity. Reading (2.32) and (2.33) as new theorems —
  they are (2.31) contracted or crossed.""")
note("N66", "The generalised derivative", r"""
`ch02.integral_gradient(Q_fn, x0, h)` converges to `gradient` with order 2 in h (printed below).
""", equation=r"\mathcal D Q = \lim_{V\to 0}\frac{1}{V}\iint_A n_i Q\,dA", ref="2.31")
note("N67", "The curl as a boundary sum", r"""
`ch02.integral_curl` → `curl`; Ch. 5 uses this picture for vorticity. **Coordinate-free remark** `N72` (§2.13): the
integral definitions of ∇, ∇·, ∇× do not depend on the coordinate system (Exercises 2.16–2.18) — cylindrical and
spherical forms in Ch. 3 / Appendix B.
""", equation=r"\nabla\times\mathbf Q = \lim_{V\to 0}\frac{1}{V}\iint_A \mathbf n\times\mathbf Q\,dA", ref="2.33")
note("N68", "Ex. 2.5", r"""
is the derivation D22 below (its recipe animated after it); the same box argument gives the continuity equation in
Ch. 4 §4.2.
""")
D("D22", "Ex. 2.5: the Cartesian divergence recovered from the outflux definition", ref="2.32 → 2.23",
  goal="""Evaluate (2.32) on a small box and watch the formula $\\partial Q_i/\\partial x_i$ (2.23) come back — face by face
  — so that the coordinate formula is seen as a consequence of "outflux per volume".""",
  assumptions="""Q twice differentiable (Taylor remainder O(Δx²); steps 2–3) · the box shrinks in all three directions
  (step 8).""",
  start=(r"\nabla\cdot\mathbf Q = \lim_{V\to0}\frac1V\iint_A \mathbf n\cdot\mathbf Q\,dA \quad\text{on a box } \Delta x_1\Delta x_2\Delta x_3 \text{ centred at } \mathbf x",
         "the outflux definition (2.32), applied to a small rectangular box with faces perpendicular to the axes (the book's face letters: EADH/FBCG ⊥ x₁, ABCD/EFGH ⊥ x₂, ABFE/DCGH ⊥ x₃)."),
  plan=["Take the two faces ⊥ x₁; write their outward normals.", "Taylor-expand Q from the centre to each face centre.",
        "Integrate over each face (midpoint rule).", "Add the pair: constants cancel, slopes add.",
        "Same for the other two pairs; add; divide by the volume; take the limit."],
  uses=["first-order Taylor expansion (Ch. 1 primer P26)", "midpoint rule: a face integral of a linear function equals the centre value × area (P83)",
        "face normals ±e₁ (C04's faces)", "orders of smallness (P68)", r"$\mathbf e_1\cdot\partial\mathbf Q/\partial x_1 = \partial Q_1/\partial x_1$ (gloss in step 9)", "(2.23) (C10)"],
  steps=[
      ("Name the two faces perpendicular to x₁",
       r"\text{EADH}:\ \mathbf n = +\mathbf e_1\ \text{at}\ x_1 + \tfrac{\Delta x_1}2;\qquad \text{FBCG}:\ \mathbf n = -\mathbf e_1\ \text{at}\ x_1 - \tfrac{\Delta x_1}2",
       r"Outward normals of a box are ±e_i (C04's faces); each of these two faces has area $\Delta x_2\Delta x_3$. We start with one pair because the other two are copies.",
       "A \"front\" face facing +x₁ and a \"back\" face facing −x₁."),
      ("Taylor-expand Q to the two face centres",
       r"[\mathbf Q]_{\rm EADH} = \mathbf Q(\mathbf x) + \tfrac{\Delta x_1}2\,\frac{\partial\mathbf Q}{\partial x_1} + \dots,\qquad [\mathbf Q]_{\rm FBCG} = \mathbf Q(\mathbf x) - \tfrac{\Delta x_1}2\,\frac{\partial\mathbf Q}{\partial x_1} + \dots",
       "First-order Taylor expansion from the box centre along ±x₁ (P26); the dots are O(Δx₁²). We expand from the centre so the two faces get opposite first-order terms.",
       "The front face sees Q plus half a step of slope; the back face sees Q minus it."),
      ("Integrate n·Q over each face by the midpoint rule",
       r"\iint_{\rm EADH}\mathbf n\cdot\mathbf Q\,dA = \mathbf e_1\cdot[\mathbf Q]_{\rm EADH}\,\Delta x_2\Delta x_3 + O(\Delta x^4)",
       "Over a face, Q varies linearly to first order, and the integral of a linear function over a rectangle equals its centre value times the area (midpoint rule, exact for linear variation, P83); the neglected variation is O(Δx²) relative, i.e. O(Δx⁴) absolute here.",
       "Face integral ≈ centre value × face area."),
      ("Add the front and back faces",
       r"\big(\mathbf e_1\cdot[\mathbf Q]_{\rm EADH} - \mathbf e_1\cdot[\mathbf Q]_{\rm FBCG}\big)\Delta x_2\Delta x_3",
       "The back face has n = −e₁, so its n·Q is minus e₁·Q. Pairing opposite faces is the whole trick: their zeroth-order terms are about to cancel.",
       "Front minus back, times the face area."),
      ("Insert the Taylor values: the constants cancel, the slopes add",
       r"\Big(\mathbf e_1\cdot\frac{\partial\mathbf Q}{\partial x_1}\Big)\Delta x_1\Delta x_2\Delta x_3 + O(\Delta x^5)",
       r"$\mathbf Q(\mathbf x) - \mathbf Q(\mathbf x) = 0$; $\tfrac{\Delta x_1}2 + \tfrac{\Delta x_1}2 = \Delta x_1$. Only the first-order change of Q across the box survives, multiplied by the volume $V = \Delta x_1\Delta x_2\Delta x_3$.",
       "What leaks through this pair is the slope of Q₁ along x₁ times the volume."),
      ("Repeat for the pairs perpendicular to x₂ and x₃",
       r"\Big(\mathbf e_2\cdot\frac{\partial\mathbf Q}{\partial x_2}\Big)V + O(\Delta x^5),\qquad \Big(\mathbf e_3\cdot\frac{\partial\mathbf Q}{\partial x_3}\Big)V + O(\Delta x^5)",
       "Steps 1–5 with the axis relabelled: faces ABCD/EFGH (±e₂) and ABFE/DCGH (±e₃). Nothing in the argument used which axis was \"1\".",
       "Each pair of faces contributes its own slope times the volume."),
      ("Add the six faces",
       r"\iint_A \mathbf n\cdot\mathbf Q\,dA = \Big(\mathbf e_1\cdot\frac{\partial\mathbf Q}{\partial x_1} + \mathbf e_2\cdot\frac{\partial\mathbf Q}{\partial x_2} + \mathbf e_3\cdot\frac{\partial\mathbf Q}{\partial x_3}\Big)V + O(\Delta x^5)",
       "The closed surface of the box is exactly the six faces; the surface integral is the sum of the face integrals (steps 5–6).",
       "The total outflux is the three slopes times the volume, plus something much smaller."),
      ("Divide by V and take the limit",
       r"\lim_{V\to0}\frac1V\iint_A \mathbf n\cdot\mathbf Q\,dA = \mathbf e_i\cdot\frac{\partial\mathbf Q}{\partial x_i}",
       r"$O(\Delta x^5)/V = O(\Delta x^2) \to 0$ as all three sides shrink (orders of smallness, P68), while the slope terms do not depend on the box size. Summation convention on i.",
       "Outflux per unit volume tends to the sum of the three directional slopes."),
      ("Recognise the Cartesian divergence",
       r"\nabla\cdot\mathbf Q = \frac{\partial Q_1}{\partial x_1} + \frac{\partial Q_2}{\partial x_2} + \frac{\partial Q_3}{\partial x_3} = \frac{\partial Q_i}{\partial x_i}",
       r"$\mathbf e_i\cdot\partial\mathbf Q/\partial x_i = \partial(\mathbf e_i\cdot\mathbf Q)/\partial x_i = \partial Q_i/\partial x_i$ because the unit vectors are constant (Cartesian). This is (2.23): the recipe returns from the definition.",
       "Divergence = the three stretching rates added — now *because* that is what leaks out of a tiny box."),
  ],
  result=(r"\nabla\cdot\mathbf Q = \frac{\partial Q_i}{\partial x_i}\quad\text{(2.23) from (2.32)}",
          "the Cartesian formula is the outflux per unit volume of a shrinking box; opposite faces cancel the value of Q and keep its slope."),
  interpret="""Ch. 4 §4.2 runs this argument with Q = ρu on a fluid box to get the continuity equation; the same
  bookkeeping with n × Q gives the curl (Ex. 2.6, D26). In curvilinear coordinates the face areas themselves vary and
  extra terms appear — Appendix B.""",
  check="""Units: [Q]/m ✓. Limit: linear Q — the O(Δx⁵) terms vanish identically and the box formula is exact at any size
  ✓. Number: Q = (x², 0, 0), x₀ = 1, h = 0.2: front 1.21·0.04, back −0.81·0.04, net 0.016, /0.008 = 2.000 = 2x₀ ✓;
  Q = (x³, 0, 0): 3 + h²/4 → order 2 in h (`integral_divergence` below, and the animation).""",
  traps="""cancelling zeroth-order terms between *non*-opposite faces. Treating a face value as exact rather than a face
  average (the midpoint rule is what makes step 3 honest). Which face carries +e₁. Writing
  $\\mathbf e_1\\cdot\\partial\\mathbf Q/\\partial x_1 = |\\partial\\mathbf Q/\\partial x_1|$.""")
nb.worked_example("Q = (x², 0, 0) around x₀ = (1, 0, 0), box side h = 0.2", r"""
1. Faces ⊥ x: at x = 1.1, n·Q = 1.21, area h² = 0.04 → +0.0484; at x = 0.9, n·Q = −0.81 → −0.0324.
2. The other four faces: $Q_2 = Q_3 = 0$ → 0.
3. Net outflux 0.0160; volume h³ = 0.008; ratio 2.000.
4. Exact ∇·Q = 2x = 2 at x₀ ✓ — exact even at finite h, because $(x+h/2)^2 - (x-h/2)^2 = 2xh$ has no h³ term.
5. For Q = (x³, 0, 0): ratio $3x^2 + h^2/4 = 3.01$ vs exact 3 — the error is $h^2/4$: order 2.
""")
nb.code("""
Qf = lambda x, y, z: np.stack([x**2, 0 * y, 0 * z])                     # Q = (x², 0, 0) as a function of coordinate arrays
for h in (0.4, 0.2, 0.1):                                                # three box sizes
    print(h, ch02.integral_divergence(Qf, [1, 0, 0], h))                 # (2.32) by six face sums on a cube of side h: 2.0 exactly — a quadratic has no h² error
hs = (0.4, 0.2, 0.1) if FAST else (0.4, 0.2, 0.1, 0.05)                  # box sizes for the convergence check (FAST: three)
Qc = lambda x, y, z: np.stack([x**3, 0 * y, 0 * z])                     # Q = (x³, 0, 0): exact ∇·Q = 3 at x = 1
errs = [abs(ch02.integral_divergence(Qc, [1, 0, 0], h) - 3.0) for h in hs]   # the error for each h
print(np.round(errs, 6), round(observed_order(hs, errs), 2))             # [0.04 0.01 0.0025 …] = h²/4, order 2.0
print(ch02.integral_gradient(lambda x, y, z: x * y * z, [1, 2, 3], 0.1))    # (2.31) on a scalar: ∇(xyz) at (1, 2, 3) = (yz, xz, xy) = (6, 3, 2)
print(ch02.integral_curl(ch02.solid_body_rotation_field([0, 0, 1.0]), [0.3, 0.2, 0], 0.1))   # (2.33): the curl of b × x by boundary sums = (0, 0, 2)
st = {k: ch02.integral_definition_convergence(k, hs=hs) for k in ("divergence", "gradient", "curl", "curl_component")}   # cached studies on the smooth test field (reused in C16)
print({k: round(s["order"], 2) for k, s in st.items()})                  # ≈ 2.0 for all four
""", explain="""
1. `integral_divergence` evaluates (2.32) on a cube of side h by six face sums; for a quadratic it is exact at any h.
2. For a cubic the error is h²/4 — observed order 2 (`observed_order`, glossed in C09).
3. (2.31) on a scalar returns the gradient; (2.33) returns the curl of the rotation field.
4. `integral_definition_convergence` repeats the measurement for the smooth sin/cos field; all the integral definitions
   converge at second order (the dicts are reused in C16).
""")
nb.check_agree("""
def six_faces(Qfn, x0, h):                                               # (2.32) by hand: one midpoint per face
    tot = 0.0                                                            # the outflux
    for ax in range(3):                                                  # the pair of faces perpendicular to axis ax
        for sgn in (+1, -1):                                             # + face and − face
            p_ = np.array(x0, float); p_[ax] += sgn * h / 2              # the face centre
            tot += sgn * Qfn(*p_)[ax] * h**2                             # n·Q at the centre × face area (midpoint rule, exact for linear variation)
    return tot / h**3                                                    # divide by the volume
assert np.allclose(six_faces(Qf, [1, 0, 0], 0.2), ch02.integral_divergence(Qf, [1, 0, 0], 0.2, n_face=1))   # one-point face rule = the hand version: 2.0
errs_mine = [abs(six_faces(Qc, [1, 0, 0], h) - 3.0) for h in hs]        # the cubic's errors by hand
assert abs(observed_order(hs, errs_mine) - 2) < 0.15                     # order 2
print(six_faces(Qf, [1, 0, 0], 0.2), "= 2x₀ by six face values; cubic error order", round(observed_order(hs, errs_mine), 2))
""")
nb.animation("""
x0, h0 = 1.0, 0.6                                                        # the box centre (along x) and its initial side
Q1 = lambda x: x**2                                                      # Q₁(x) = x² (Q₂ = Q₃ = 0): only the x-faces leak
fig, ax = plt.subplots(figsize=(7.5, 4.2))                               # the box seen from above: x₁ horizontal, x₂ vertical
fig.set_layout_engine("none")                                            # fixed layout: the automatic layout engine would re-run on every frame
ax.set_xlim(0.2, 2.0); ax.set_ylim(-0.8, 0.9); ax.set_aspect("equal"); ax.set_xlabel("x₁ [m]"); ax.set_yticks([])   # fixed axes
(box,) = ax.plot([], [], color=COLORS["ink"], lw=2)                      # the box outline
(centre,) = ax.plot([x0], [0], "o", color=COLORS["accent"])              # x₀
tf = ax.text(0, 0, "", color=COLORS["orange"], ha="left", fontsize=9)    # front-face value (n = +e₁)
tb = ax.text(0, 0, "", color=COLORS["blue"], ha="right", fontsize=9)     # back-face value (n = −e₁)
ts_ = ax.text(0.25, 0.75, "", fontsize=9)                                # the running bookkeeping line
frames_txt = ["the box around x₀ = 1 with side Δx; Q = (x², 0, 0)",      # frame 1: the setting
              "front face (n = +e₁): n·Q = Q₁(x₀ + Δx/2) = Q₁ + (Δx/2)·Q₁′ + …",     # frame 2: D22 step 2, front
              "back face (n = −e₁): n·Q = −Q₁(x₀ − Δx/2) = −Q₁ + (Δx/2)·Q₁′ − …",    # frame 3: D22 step 2, back
              "the four other faces: Q₂ = Q₃ = 0 → they contribute nothing",           # frame 4: D22 step 6 (trivial here)
              "add front + back: the Q₁ constants cancel …",                            # frame 5: D22 step 4
              "… the slopes add: (Δx·Q₁′) × face area = Q₁′ · V",                       # frame 6: D22 step 5
              "divide by V = Δx³:  (1/V)∮ n·Q dA = Q₁′ + O(Δx²)",                        # frame 7: D22 step 8
              "shrink the box …", "… and again …", "… (1/V)∮ → ∂Q₁/∂x₁ = 2x₀ = 2.000"]   # frames 8–10: the limit
def update(i):                                                           # frame i of the bookkeeping
    h = h0 if i < 7 else h0 * (0.6 ** (i - 6))                           # the box shrinks in the last frames
    box.set_data([x0 - h/2, x0 + h/2, x0 + h/2, x0 - h/2, x0 - h/2], [-h/2, -h/2, h/2, h/2, -h/2])   # redraw the square
    front, back = Q1(x0 + h/2), -Q1(x0 - h/2)                            # the two face values of n·Q
    tf.set_position((x0 + h/2 + 0.03, 0)); tf.set_text(f"+{front:.3f}" if i >= 1 else "")   # front value appears from frame 2
    tb.set_position((x0 - h/2 - 0.03, 0)); tb.set_text(f"{back:.3f}" if i >= 2 else "")     # back value from frame 3
    est = (front + back) * h**2 / h**3                                   # (1/V)∮ n·Q dA for this box: (front + back) × face area / volume
    tail = f"   Δx = {h:.3f}: (1/V)∮ = {est:.3f}" if i >= 6 else ""      # the number, once the division is done
    ts_.set_text(frames_txt[i] + tail)                                   # the bookkeeping line
    return (box, tf, tb, ts_)                                            # the changed artists
show_animation(animate(update, frames=len(frames_txt), fig=fig, interval=900), player="frames")   # step through the bookkeeping
""", explain="""
1. Ex. 2.5's bookkeeping for Q = (x², 0, 0) around x₀ = 1, seen from above: the front face carries +Q₁(x₀ + Δx/2), the
   back face −Q₁(x₀ − Δx/2); the four other faces carry nothing.
2. Frames 5–7 add the pair: the constants die, the slopes survive, and dividing by V gives (1/V)∮ = 2.000 — exactly
   ∂Q₁/∂x₁, at *every* box size here because the field is quadratic (no h² error); the last frames shrink the box.
""")
see_read_change("A square box around x₀ with its two x-face values printed at the sides; a running line of bookkeeping that ends in (1/V)∮ = 2.000 at every size as the box shrinks.",
                "The divergence is the part of the outflux that opposite faces do *not* cancel — the first-order (slope) term; the constants cancel in pairs.",
                "…x₀ moved to (2, 0, 0): every face value doubles its slope term and (1/V)∮ → 4 (= 2x₀); for Q = (x³, 0, 0) the estimate would read 3 + Δx²/4 and only settle as the box shrinks.")
nb.plotly("""
xx = np.linspace(0.3, 1.7, 141)                                          # x-axis for the curve
def secant(h):                                                           # one slider position: box side h around x₀ = 1
    xa, xb = 1 - h / 2, 1 + h / 2                                        # the two x-faces
    est = (Q1c(xb) - Q1c(xa)) / h                                        # (1/V)∮ n·Q dA for Q = (x³, 0, 0): the secant slope of Q₁
    return {"Q₁(x) = x³": (xx, Q1c(xx)),                                 # the field along the axis (fixed)
            "the two face values": ([xa, xb], [Q1c(xa), Q1c(xb)]),       # what the box sees
            "secant = (1/V)∮ n·Q dA": ([0.3, 1.7], [Q1c(xa) + est * (0.3 - xa), Q1c(xa) + est * (1.7 - xa)]),   # the box estimate as a slope
            "tangent = ∂Q₁/∂x₁ = 3": ([0.3, 1.7], [1 + 3 * (0.3 - 1), 1 + 3 * (1.7 - 1)])}   # the exact derivative
Q1c = lambda x: x**3                                                     # Q₁ for the cubic field
fig = slider_figure(secant, "h", np.logspace(0, -2, 20), unit="m (box side)", xlabel="x₁ [m]", ylabel="Q₁ [–]",
                    title="(2.32) on a shrinking box is a secant slope: it glides onto the tangent ∂Q₁/∂x₁ with error h²/4",
                    xrange=[0.3, 1.7], yrange=[-1.5, 3.5], height=460, active=0,
                    modes={"the two face values": "markers"})
fig.show()
print({round(h, 3): round((Q1c(1 + h/2) - Q1c(1 - h/2)) / h - 3, 5) for h in (1.0, 0.5, 0.25, 0.1)})   # the error h²/4: 0.25, 0.0625, 0.0156, 0.0025
""", explain="""
1. For Q = (x³, 0, 0) only the two x-faces leak, so (1/V)∮ n·Q dA is the *secant slope* of Q₁ between the two face
   centres; the exact divergence is the tangent slope 3.
2. Dragging h down glides the secant onto the tangent; the printed error is h²/4 — (2.32) is a limit, approached at
   second order.
""")
see_read_change("A cubic curve, two dots on it (the face values), a secant line through them and the fixed tangent at x₀ = 1; as h shrinks the secant turns onto the tangent.",
                "(2.32) is a limit — the box estimate is the secant slope, the derivative is the tangent slope; the approximation improves as h².",
                "…Q were linear or quadratic: the secant would coincide with the tangent at every h (exact at any box size).")
nb.md(r"""
**What would change if…** you did the same with a *loop* instead of a closed surface? Circulation per unit area — the
curl, by Stokes' theorem. §2.13.
""")

# =====================================================================================================================
# A.11 §2.13 Stokes' Theorem — C16 (D26)
# =====================================================================================================================
nb.section("2.13", "Stokes' Theorem", intro="""
**What is this section about?** The circulation of a velocity field around a closed loop equals the flux of its curl
through any surface the loop bounds. Read on a tiny loop it defines the normal curl as circulation per unit area — and
Kelvin's theorem (Ch. 5), lift (Ch. 6, 14) and vortex tubes all live on it.
""")
core("C16", "Stokes' theorem and circulation", r"""
Why does the total spin inside a loop equal the flow's tendency to run around the loop — and when does that fail?
""")
nb.md(r"""
#### The problem in plain words

Put a closed loop into a flow and add up the velocity component along the loop all the way round: that is the
circulation. It is the quantity Kelvin's theorem conserves, the quantity that produces lift on a wing, the strength of
a tornado. Stokes' theorem says it equals the curl integrated over the loop's interior — so circulation is "total spin
inside", not "how curved the streamlines are".
""")
nb.md(r"""
#### The idea

```
one small rectangle:  ∮ u·t ds = (∂u₂/∂x₁ − ∂u₁/∂x₂) ΔA = (∇×u)·n ΔA        (Ex. 2.6 bookkeeping)
tile the surface with rectangles: every interior edge is walked twice, opposite ways → cancels
what is left:  ∮_C u·t ds = ∬_A (∇×u)·n dA        orientation: thumb along n, fingers along t (right hand)
```
""")
P("P86", "line integral of a vector field around a loop", r"""
Parametrise the loop by arc length s (or an angle); at each point take the velocity component along the unit tangent,
u·t, and add it up times the small length ds: $\oint_C \mathbf u\cdot\mathbf t\,ds$. (🔁 Ch. 1 P35 integrated p dv
along a path; this is the vector version.) A circle: x = c + R(cos θ, sin θ), t = (−sin θ, cos θ), ds = R dθ; a
rectangle: four straight sides, t = ±e₁ or ±e₂.
""", code="""
th = (np.arange(200) + 0.5) * 2 * np.pi / 200              # 200 midpoints round a unit circle
x, y = np.cos(th), np.sin(th)                              # the points
tx, ty = -np.sin(th), np.cos(th)                           # the unit tangent, counterclockwise
u1, u2 = -y, x                                             # u = b × x with b = e_3
Gamma_mine = np.sum((u1 * tx + u2 * ty) * 2 * np.pi / 200) # Σ u·t ds
print(Gamma_mine)                                          # 6.2832 = 2π
""")
nb.md(r"""
#### The maths, step by step

1. **Orientation** (N69): pick the outside of A → n; t runs counterclockwise about n (right-hand rule, P74);
   $\mathbf n_c = \mathbf n\times\mathbf t$ is the in-surface normal to C pointing **into** A, so $\mathbf t = \mathbf n_c
   \times\mathbf n$ and $(\mathbf n_c, \mathbf n, \mathbf t)$ is right-handed.
2. **Statement (2.34):** $\iint_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA = \oint_C \mathbf u\cdot\mathbf t\,ds$; the
   right side is the **circulation**.
3. **Small loop (2.35):** $\mathbf n\cdot(\nabla\times\mathbf u) = \lim_{A\to0}(1/A)\oint_C\mathbf u\cdot\mathbf t\,ds$
   — the normal curl is circulation per unit area (N71).
4. **Proof (D26):** one rectangle (Ex. 2.6 run backwards), tiling, cancellation; corollary $\oint\nabla\phi\cdot\mathbf t
   \,ds = 0 \Rightarrow \nabla\times\nabla\phi = 0$.
""")
note("N69", "⚠️ The orientation rule", r"""
the book writes $\mathbf t = \mathbf n_c\times\mathbf n$; this is right-handed only if $\mathbf n_c$ points along the
surface *into* A (up the cap in the book's Fig. 2.10). Our code fixes t by the right-hand rule about `normal` and
checks $\oint(\mathbf x - \mathbf c)\times\mathbf t\,ds = 2A\,\mathbf n$ (`ch02.planar_loop`, `ch02.boundary_tangent`).
**Fig. 2.10** `N70` is our 3-D cap below with n, n_c, t at a boundary point.
""")
note("N71", "The normal curl as circulation per unit area", r"""
`ch02.integral_curl_component(u, x0, n, h)` converges to n·curl with order 2; pointer: Ch. 3 vorticity = 2 × angular
velocity, Ch. 5.
""", equation=r"\mathbf n\cdot(\nabla\times\mathbf u) = \lim_{A\to 0}\frac{1}{A}\oint_C \mathbf u\cdot\mathbf t\,ds", ref="2.35")
note("N73", "Ex. 2.6, stated correctly", r"""
for the rectangle Δy × Δz in the plane x = const with $\mathbf n = \mathbf e_x$, the sides at y ± Δy/2 run along
$\pm\mathbf e_z$ (integrand $u_z$), the sides at z ∓ Δz/2 run along $\pm\mathbf e_y$ (integrand **$u_y$** — the book
prints $u_z$ there by mistake; its own limit $\partial u_z/\partial y - \partial u_y/\partial z$ confirms $u_y$);
midpoint values × side lengths, divide by ΔyΔz, take the limit: (2.25)'s first component. The y- and z-components follow
by cyclic relabelling. D26 runs the same bookkeeping in the x₃-plane.
""")
D("D26", "Stokes' theorem for a planar surface: one rectangle, then tiling", ref="2.34",
  goal="""Prove Stokes' theorem for a flat surface — the circulation round the boundary equals the curl integrated inside
  — by doing one tiny rectangle exactly (Ex. 2.6's bookkeeping run backwards) and tiling; then read off the corollary
  ∇×∇φ = 0.""",
  assumptions="""u continuously differentiable on A (steps 2–3, 8) · A a flat surface with piecewise smooth boundary C
  (steps 7–9; a curved A is done patch by patch with t = n_c × n — stated) · every tile oriented by the same n (step 8).""",
  start=(r"\text{a rectangle } \Delta x_1\times\Delta x_2 \text{ centred at } (x_1, x_2) \text{ in the plane } x_3 = \text{const},\quad \mathbf n = \mathbf e_3,\ \mathbf t \text{ counterclockwise}",
         "the smallest loop we can integrate around by hand."),
  plan=["Fix the orientation of the four sides.", "Midpoint values on each side.", "Add bottom + top and right + left: Taylor differences.",
        "Recognise the curl.", "Tile the surface: interior edges cancel; refine.", "Corollary for gradient fields."],
  uses=["line integral of a vector field round a loop (P86)", "orientation t counterclockwise about n (N69, P74)",
        "first-order Taylor expansion (Ch. 1 P26)", "midpoint values on a side (P83, as in D22)", "the curl's third component (2.25) (D12)",
        "Ex. 2.6's bookkeeping (N73)", "the chain rule for dφ along a curve (Ch. 1 P49; step 10)"],
  steps=[
      ("Fix the four sides and their tangents",
       r"\text{bottom}\ (x_2 - \tfrac{\Delta x_2}2):\ \mathbf t = +\mathbf e_1;\quad \text{right}:\ +\mathbf e_2;\quad \text{top}:\ -\mathbf e_1;\quad \text{left}:\ -\mathbf e_2",
       "Counterclockwise about n = e₃ (right-hand rule: thumb along n, fingers along t — N69, P74). The orientation decides every sign below, so we fix it first.",
       "Walk the rectangle anticlockwise seen from above."),
      ("Write each side's integral with its midpoint value",
       r"\int_{\rm side}\mathbf u\cdot\mathbf t\,ds = (\mathbf u\cdot\mathbf t)_{\rm mid}\times\text{length} + O(\Delta^3)",
       "Along a short side u varies linearly to first order and the integral of a linear function equals its midpoint value times the length (as in D22 step 3); the error is O(Δ³) per side.",
       "Each side contributes \"velocity along it at its middle times its length\"."),
      ("Add bottom and top",
       r"\big[u_1(x_1, x_2 - \tfrac{\Delta x_2}2) - u_1(x_1, x_2 + \tfrac{\Delta x_2}2)\big]\Delta x_1",
       r"On the bottom $\mathbf u\cdot\mathbf t = +u_1$, on the top $-u_1$ (step 1); both sides have length Δx₁ and midpoints directly below/above the centre.",
       "The horizontal sides see u₁ below minus u₁ above."),
      ("Taylor-expand the difference",
       r"-\,\frac{\partial u_1}{\partial x_2}\,\Delta x_1\Delta x_2 + O(\Delta^4)",
       r"$u_1(x_2 \mp \tfrac{\Delta x_2}2) = u_1 \mp \tfrac{\Delta x_2}2\,\partial_2u_1 + \dots$ (P26); the two u₁'s cancel and the halves add to Δx₂ with a minus sign.",
       "If u₁ grows upward, the top side walks against it: a clockwise contribution."),
      ("Add right and left the same way",
       r"\big[u_2(x_1 + \tfrac{\Delta x_1}2, x_2) - u_2(x_1 - \tfrac{\Delta x_1}2, x_2)\big]\Delta x_2 = +\frac{\partial u_2}{\partial x_1}\,\Delta x_1\Delta x_2 + O(\Delta^4)",
       r"On the right $\mathbf u\cdot\mathbf t = +u_2$, on the left $-u_2$; Taylor along x₁; this time the plus side is at $+\Delta x_1/2$ so the sign is positive.",
       "If u₂ grows to the right, the right side walks with it: an anticlockwise contribution."),
      ("Add the four sides and recognise the curl",
       r"\oint_{\rm rect}\mathbf u\cdot\mathbf t\,ds = \Big(\frac{\partial u_2}{\partial x_1} - \frac{\partial u_1}{\partial x_2}\Big)\Delta A + O(\Delta^4) = (\nabla\times\mathbf u)\cdot\mathbf n\,\Delta A",
       r"Steps 4 and 5 added, with $\Delta A = \Delta x_1\Delta x_2$; the bracket is $(\nabla\times\mathbf u)_3$ by (2.25) (D12), and n = e₃ picks that component. Dividing by ΔA and letting Δ → 0 gives (2.35) — Ex. 2.6 in the x₃-plane.",
       "One small loop's circulation is the normal curl times its area."),
      ("Tile the surface with small rectangles and add",
       r"\sum_k\oint_{C_k}\mathbf u\cdot\mathbf t\,ds = \sum_k(\nabla\times\mathbf u)\cdot\mathbf n\,\Delta A_k \to \iint_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA",
       "Step 6 holds for every tile, all oriented counterclockwise about the same n; the right side is a Riemann sum of a continuous function, which tends to the surface integral as the tiles shrink.",
       "The sum of all the little circulations is the total curl inside."),
      ("Cancel the interior edges",
       r"\sum_k\oint_{C_k}\mathbf u\cdot\mathbf t\,ds = \oint_{\text{outer edges}}\mathbf u\cdot\mathbf t\,ds",
       "An edge shared by two neighbouring tiles is traversed once in each direction (both tiles anticlockwise), so u·t flips sign while u is the same: the two contributions cancel. Only the edges on the outside survive.",
       "Every inner edge is walked twice, opposite ways, and drops out."),
      ("Refine: the outer edges become C",
       r"\oint_C\mathbf u\cdot\mathbf t\,ds = \iint_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA",
       "As the tiles shrink the staircase of outer edges tends to the boundary curve C with tangent t (u continuous, C piecewise smooth), so the left side of step 8 tends to the circulation round C. With step 7 this is (2.34). For a curved A: the same on each flat patch, with t = n_c × n at the boundary.",
       "Stokes' theorem: circulation round the edge equals curl flux through the surface."),
      ("Corollary for a gradient field",
       r"\oint_C\nabla\phi\cdot\mathbf t\,ds = \oint_C d\phi = 0\ \Rightarrow\ \nabla\times\nabla\phi = 0",
       r"Along the curve, $\nabla\phi\cdot\mathbf t\,ds = d\phi$ (chain rule, P49): the integral of an exact change round a closed loop is φ(end) − φ(start) = 0. Then (2.34) forces $\iint(\nabla\times\nabla\phi)\cdot\mathbf n\,dA = 0$ for every A, so the integrand vanishes (Exercise 2.20).",
       "A field that is the slope of something has no circulation and no curl."),
  ],
  result=(r"\iint_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA = \oint_C\mathbf u\cdot\mathbf t\,ds \quad(2.34),\qquad \nabla\times\nabla\phi = 0",
          "the circulation round a loop is the total curl inside it; gradient fields have none."),
  interpret="""Circulation is "total spin inside", not "curved streamlines": a straight shear flow has it, a whirlpool
  with a singular core has it only through the core. Kelvin's theorem (Ch. 5) follows the circulation of a material
  loop; the lift of a wing (Ch. 6, 14) is proportional to it; a vortex tube's strength is the same through every
  cross-section because ∇·(∇×u) = 0. This planar case is also known as Green's theorem. The theorem fails when A
  contains a point where u is not differentiable (the vortex core) — then the loop's circulation is fixed by the core,
  whatever the loop.""",
  check="""Units: (1/s)·m² = (m/s)·m = m²/s ✓. Limit: u uniform: both sides 0 ✓. Number: u = b × x, b = e₃, unit disc:
  circulation 2πR = 6.283, curl 2 × area π = 6.283 ✓ (`stokes_theorem_check` below); shear u₁ = Γx₂, unit square: −Γ
  both sides ✓; u = ∇(x² − y²): circulation ~1e-16 ✓; the irrotational vortex with the core inside: 2πK vs "curl flux 0"
  — the hypothesis (u differentiable on A) fails, and the code says so. On our grid ∇×∇φ = 0 and ∇·(∇×u) = 0 hold to
  round-off (< 1e-11), because the stencils commute exactly.""",
  traps="""the orientation of the horizontal sides (the top runs in −e₁). Cancelling interior edges when tiles are
  oriented differently. Applying (2.34) across a singular core. Confusing Γ (circulation) with Γ (shear rate) — the
  code names differ (`circulation` vs `Gamma`).""")
nb.worked_example("u = b × x with b = e₃ around the unit circle", r"""
1. $\mathbf u = (-x_2, x_1, 0)$; on the circle $\mathbf x = (\cos\theta, \sin\theta, 0)$, $\mathbf t = (-\sin\theta,
   \cos\theta, 0)$: $\mathbf u\cdot\mathbf t = \sin^2\theta + \cos^2\theta = 1$.
2. Circulation $= \oint 1\,ds = 2\pi R = 2\pi = 6.283$.
3. Curl (C11) = 2b = (0, 0, 2); flux through the disc $= 2\times\pi R^2 = 2\pi$ ✓.
4. Per unit area: $2\pi/\pi = 2 = (\nabla\times\mathbf u)_3$ ✓ (2.35).
5. Reverse n → −e₃: t reverses, circulation −2π, flux −2π: both sides flip together.
""")
nb.code("""
u_fn = ch02.solid_body_rotation_field([0, 0, 1.0])                        # u = b × x, b = e_3 (a callable VectorField)
loop = ch02.planar_loop([0, 0, 0], [0, 0, 1], 1.0, 200)                    # the unit circle, counterclockwise about e_3 (right-hand rule, N69)
disc = ch02.planar_disc([0, 0, 0], [0, 0, 1], 1.0, 20, 40)                 # the disc it bounds, normals e_3
print(ch02.circulation(u_fn, loop), 2 * np.pi)                             # ∮ u·t ds = 6.2832 = 2π
chk = ch02.stokes_theorem_check(u_fn, loop, disc)                          # both sides of (2.34): lhs = ∬ (∇×u)·n dA, rhs = ∮ u·t ds
print(chk)                                                                 # StokesCheck(lhs 6.283, rhs 6.283, hypothesis_ok=True, …)
print(ch02.integral_curl_component(ch02.shear_field(1.0), [0, 0], None, 0.1))   # (2.35) on the shear flow u_1 = Γx_2 (2-D field, square loop of side 0.1): −1.0 = −Γ
loop2 = ch02.planar_loop([0, 0], None, 1.0, 400)                           # a 2-D loop (plane z = 0) for the plane fields
print(ch02.circulation(ch02.potential_field(), loop2))                     # u = ∇(x₁² − x₂²): circulation 0.0 (D26 step 10, Exercise 2.20)
uv = ch02.irrotational_vortex_field(1.0)                                   # u_θ = K/r with K = 1: curl-free except at the core
off = ch02.planar_loop([2, 0], None, 0.5, 400)                             # a loop that does NOT enclose the core
print(ch02.circulation(uv, loop2), 2 * np.pi, ch02.circulation(uv, off))   # 6.2832 = 2πK around the core; ~1e-16 away from it
chk2 = ch02.stokes_theorem_check(uv, loop2, ch02.planar_disc([0, 0], None, 1.0, 20, 40))   # the core sits inside A …
print(chk2.hypothesis_ok, "—", chk2.note)                                  # … False: Stokes' hypothesis fails (u not differentiable on A); the code reports, never asserts
ok_lap = np.abs(ch02.curl(ch02.gradient(g3.X**2 * g3.Y - g3.Z**3, g3.h), g3.h)).max()   # ∇×∇φ on the C10 grid
ok_div = np.abs(ch02.divergence(w, g3.h)).max()                            # ∇·(∇×u) for the rotation field's curl
print(f"on the grid: max|∇×∇φ| = {ok_lap:.1e}, max|∇·(∇×u)| = {ok_div:.1e}   (round-off, < 1e-11: the stencils commute exactly)")
""", explain="""
1. The loop and the disc are built right-handed about $\\mathbf e_3$ (`planar_loop` checks ∮(x − c) × t ds = 2A n).
2. The circulation of the rotation field is 2π; `stokes_theorem_check` evaluates both sides of (2.34) and they agree.
3. Ex. 2.6's recipe on the shear flow gives −Γ per unit area — (2.35).
4. A gradient field has zero circulation (the D26 corollary).
5. The irrotational vortex has circulation 2πK around any loop enclosing the core and 0 around one that does not.
6. With the core inside A the check flags `hypothesis_ok=False`: the *hypothesis* fails, not the theorem.
7. On our grid the two identities ∇×∇φ = 0 and ∇·(∇×u) = 0 hold to round-off, not merely to truncation, because the
   central/one-sided stencils commute exactly.
""")
nb.check_agree("""
r_ = (np.arange(20) + 0.5) / 20; dth = 2 * np.pi / 40                              # a polar midpoint grid on the unit disc (20 radii × 40 angles)
flux_mine = np.sum(2.0 * r_[:, None] * (1 / 20) * dth * np.ones((20, 40)))         # ∬ (∇×u)_3 dA with (∇×u)_3 = 2: Σ 2 · r dr dθ
assert np.allclose(Gamma_mine, flux_mine, rtol=1e-3)                                # the primer's circulation sum = the curl flux (Stokes, by hand)
assert np.allclose(Gamma_mine, ch02.circulation(u_fn, loop), rtol=1e-6)             # = the library's circulation
print(Gamma_mine, flux_mine, "= 2π: circulation round the rim = curl × area")
""")
nb.plotly("""
from scripts.ch02_fig2_10_stokes import stokes_cap_figure                 # drawing helper (P64): Fig. 2.10 as ours
fig = stokes_cap_figure()                                                 # a hemispherical cap A, its rim C, and at one rim point the triad n (outward), n_c (up the cap), t (along the rim)
fig.show()                                                                # rotate
""", explain="""
1. The cap is A, the rim is C; small purple arrows along C show the sense of traversal (counterclockwise seen from the
   outside).
2. At the rim point P the three unit vectors are drawn: n (from the chosen outside), $\\mathbf n_c$ (perpendicular to C,
   tangent to A, pointing *into* A — up the cap) and t = $\\mathbf n_c\\times\\mathbf n$; the title confirms
   det[$\\mathbf n_c$, n, t] = +1 (`ch02.boundary_tangent`).
""")
see_read_change("A teal hemispherical cap, its purple rim with small arrows, and three perpendicular unit vectors at one rim point: n outward, $\\mathbf n_c$ up the cap, t along the rim.",
                "$\\mathbf n_c\\times\\mathbf n = \\mathbf t$ — check with the right hand (P74): fingers from $\\mathbf n_c$ to n, thumb along t.",
                "…you chose the inside as outside (`stokes_cap_figure(flip=True)`): n and t both reverse; $\\mathbf n_c$ stays. Both sides of (2.34) change sign together.")
nb.animation("""
u_sh2 = ch02.shear_field(1.0)                                              # u_1 = Γ x_2 with Γ = 1/s (a 2-D field)
hs_loop = np.geomspace(2.0, 0.1, 10 if not FAST else 6)                    # loop sides from 2 to 0.1 m
fig, ax = plt.subplots(figsize=(6.5, 4.6))                                 # one panel
fig.set_layout_engine("none")                                              # fixed layout: the automatic layout engine would re-run on every frame
ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4); ax.set_aspect("equal"); ax.set_xlabel("x₁ [m]"); ax.set_ylabel("x₂ [m]")   # fixed axes
yy = np.linspace(-1.3, 1.3, 9); ax.quiver(np.zeros(9) - 1.2, yy, yy, 0 * yy, color=COLORS["muted"], scale=12, width=0.004)   # the shear profile u_1 = x_2 (grey)
(sq_l,) = ax.plot([], [], color=COLORS["accent"], lw=2)                    # the square loop
wheel = plt.Circle((0.9, -0.9), 0.22, fill=False, color=COLORS["rose"], lw=2); ax.add_patch(wheel)   # a paddle wheel spinning at ½(∇×u)_3 = −Γ/2
(spoke,) = ax.plot([], [], color=COLORS["rose"], lw=2)
txt = ax.text(-1.35, 1.2, "", fontsize=9)
arrows = [ax.annotate("", (0, 0), (0, 0), arrowprops=dict(arrowstyle="->", color=COLORS["orange"], lw=2)) for _ in range(4)]   # the four side contributions
def update(i):                                                             # frame i: loop of side h
    h = hs_loop[i]
    lp = ch02.rectangle_loop([0, 0], None, h, h, 64)                       # square loop of side h, counterclockwise about e_3
    circ = ch02.circulation(u_sh2, lp)                                     # ∮ u·t ds
    sq_l.set_data([-h/2, h/2, h/2, -h/2, -h/2], [-h/2, -h/2, h/2, h/2, -h/2])
    # the four sides: bottom (t = +e_1, u·t = u_1 = −h/2 → walks against), top (t = −e_1, u·t = −u_1 = −h/2), right/left (u_2 = 0)
    for arr, (xa, ya, dx_, dy_) in zip(arrows, [(-h/4, -h/2, -0.25 * h, 0), (h/4, h/2, -0.25 * h, 0), (h/2, 0, 0, 0), (-h/2, 0, 0, 0)]):
        arr.set_position((xa, ya)); arr.xy = (xa + dx_, ya + dy_)          # the bottom and top sides both push clockwise; the vertical sides contribute nothing
    ang_ = -0.5 * 1.0 * i * 0.6                                            # the wheel turns clockwise at −Γ/2 (frame time 0.6 s)
    spoke.set_data([0.9, 0.9 + 0.22 * np.cos(ang_)], [-0.9, -0.9 + 0.22 * np.sin(ang_)])
    txt.set_text(f"side h = {h:.2f} m:  Γ_circ = ∮u·t ds = {circ:+.4f} m²/s,  A = {h*h:.4f} m²,  Γ_circ/A = {circ / h**2:+.3f} = (∇×u)₃ = −Γ")
    return (sq_l, spoke, txt, *arrows)
show_animation(animate(update, frames=len(hs_loop), fig=fig, interval=700), player="frames")   # step through the sizes
""", explain="""
1. A square loop of side h shrinks about the origin in the shear flow u₁ = Γx₂ (grey profile on the left); each frame
   prints the circulation, the area and their ratio.
2. The orange arrows show the two horizontal sides' contributions (both clockwise: the bottom side walks against the
   flow, the top side walks against it too because u₁ is larger there); the vertical sides contribute nothing (u₂ = 0).
3. Γ_circ/A is −Γ at *every* size — the field is linear, so Ex. 2.6's midpoint bookkeeping is exact; the paddle wheel
   turns clockwise at −Γ/2, half the curl.
""")
see_read_change("A shrinking square with two orange side-arrows both pointing clockwise; the text line shows Γ_circ shrinking with the area while the ratio stays at −1.000.",
                "(2.35) is a statement about the ratio circulation/area, not about the loop: here −Γ at every size because the field is linear.",
                "…the sin/cos smooth field: Γ_circ/A wanders at large h and settles at small h with error ∝ h² — the convergence plot below.")
nb.figure("""
fig, ax = plt.subplots(figsize=(6, 3.8))                                    # log–log
for k, col in (("divergence", COLORS["teal"]), ("curl_component", COLORS["accent"])):   # C15's cached studies
    s = st[k]                                                               # dict(h, err, order)
    lab = "(2.32) divergence: (1/V)∮ n·Q dA" if k == "divergence" else "(2.35) normal curl: (1/A)∮ u·t ds"   # legend text
    ax.loglog(s["h"], s["err"], "o-", color=col, label=f"{lab}, order {s['order']:.2f}")   # error vs size
hr = np.array([min(st["curl_component"]["h"]), max(st["curl_component"]["h"])])   # the h-range of the reference line
ax.loglog(hr, st["curl_component"]["err"].max() * (hr / hr.max())**2, "--", color=COLORS["muted"], label="slope 2 reference")   # ∝ h²
ax.set_xlabel("box / loop size h [m]"); ax.set_ylabel("|estimate − exact| [1/s]"); ax.legend(fontsize=8)   # labels
ax.set_title("Both integral definitions converge at second order in the size")   # the message
plt.show()                                                                  # display
""", see="Two straight lines of slope 2 on log–log axes, parallel to the dashed reference.",
   read="Halving the box or loop size quarters the error of the integral definitions — the O(Δ²) remainders of D22 step 8 and D26 step 6.",
   change="…a polynomial field of degree ≤ 2: both errors would sit at round-off for every h (the midpoint bookkeeping is exact).")
nb.md(r"""
🎮 *The explainer for this theorem* — `stokes_circulation_loop` — is embedded above in C11 (where the curl first bites):
use its rectangle mode for Ex. 2.6's four sides and the irrotational-vortex preset for the failing hypothesis.
""")
nb.md(r"""
**What would change if…** the loop enclosed a singular point (the vortex core)? The field is not differentiable inside,
Stokes' hypothesis fails, and the circulation is 2πK however small the loop — the "irrotational vortex with
circulation" of Ch. 5 and Ch. 6. Everywhere else, (2.34) holds to the last digit.
""")

# =====================================================================================================================
# A.12 §2.14 Comma Notation — N74 (tagged → C10)
# =====================================================================================================================
nb.section("2.14", "Comma Notation", intro="""
**What is this section about?** One more shorthand: a comma in the subscript means "partial derivative with respect to
the following index".
""")
note("N74", "Comma notation", r"""
the divergence and curl of C10/C11 in one line each. A comma index behaves like a tensor index (the chain rule, Ch. 1
P49, gives $\partial/\partial x'_j = C_{ij}\,\partial/\partial x_i$ — the vector rule (2.8) for ∂), but only in
Cartesian coordinates; the book adopts it in §5.6. Watch for the comma: $u_{i,j}$ is the velocity gradient of N52
(C10), $u_{ij}$ would be something else.
""", equation=r"A_{,i} \equiv \partial A/\partial x_i, \qquad \nabla\cdot\mathbf u = u_{i,i}, \qquad (\nabla\times\mathbf u)_i = \varepsilon_{ijk}u_{k,j}", ref="2.36")
nb.code("""
print(ch02.expand_indices_str("u_i,i"))                 # (2.36): the comma is parsed as ∂ → Derivative(u_1, x1) + Derivative(u_2, x2) + Derivative(u_3, x3) = (2.23)
print(ch02.comma_to_partial("u_i,j"))                   # the text form: ∂u_i/∂x_j — the velocity gradient of N52
print(ch02.expand_indices_str("eps_ijk u_k,j"))         # the curl in comma notation expands to the three lines of (2.25)
""", explain="""
1. `expand_indices_str` reads the comma as a partial derivative and expands the divergence.
2. `comma_to_partial` writes the derivative out in words.
3. The curl in comma notation expands to D12's three components.
""")

# =====================================================================================================================
# A.13 end matter — S01, S02, summary
# =====================================================================================================================
nb.pointer("""
Exercises 2.1–2.20 are not reproduced (the text is the book's). The identities we rely on from them — (2.7) and
orthogonality (Exercises 2.2, 2.8), the invariants (2.9), the tensor examples (2.10), isotropy (2.11), the dot and cross
products (2.12–2.14), ε–δ (2.5, 2.7), ∇×∇φ = 0 and ∇·∇×u = 0 (2.19, 2.20) — are derived above (D02, D03, D09, D18,
D26) and tested in `tests/test_ch02.py`. `S01`
""")
nb.pointer("""
Literature: Sommerfeld's tetrahedron argument for (2.12) is our D05 + D06; Aris and Prager are the classical tensor
references. Not needed to continue. `S02`
""")
nb.summary(
    clicked=[
        "A repeated index is a loop; the free indices are the shape of the answer (C01).",
        "C_ij = e_i·e'_j; the columns of C are the new axes; C Cᵀ = I, det C = +1 (C02).",
        "x' = Cᵀx defines a vector: any triple that fails the test is not one (C03).",
        "τ_ij: first index the face, second the force; positive out of a + face (C04).",
        "f = n·τ: the force per area on any plane, from a shrinking tetrahedron (C05).",
        "τ' = CᵀτC — one C per index — is what 'tensor' means (C06).",
        "Closed index chains (trace, I₂, det) are frame-independent (C07).",
        "ε_ijk encodes right-handed perpendicularity; ε ε = δδ − δδ closes vector identities (C08).",
        "∇φ is perpendicular to the level sets and points up the steepest slope (C09).",
        "∇·u is the total stretching rate — what leaks out per unit volume (C10).",
        "∇×u is twice the local spin; a straight shear flow has it (C11).",
        "G = S + A: stretch plus spin; the vector of A is ½∇×u (the book's R = 2A carries ∇×u); symmetric tensors cannot see A (C12).",
        "A symmetric tensor has three perpendicular shear-free axes; its eigenvalues bound every normal stress (C13).",
        "Gauss: derivative inside = normal outside; interior faces cancel (C14).",
        "Divergence = outflux per unit volume; the Cartesian formula falls out of a box (C15).",
        "Stokes: circulation round a loop = curl flux through it; curl = circulation per area (C16).",
    ],
    feeds_forward=[
        "Ch. 3: G = S + R from ∂u_i/∂x_j (C12), principal strain rates (C13), vorticity = 2 × angular velocity (C11, D15)",
        "Ch. 4: Cauchy's equation from ∮ n·τ dA and Gauss (C05, C14), continuity from the shrinking box (C15), dissipation τ:S (C12), the Newtonian stress from δ (C08)",
        "Ch. 5: vorticity identities with ε–δ (C08), Kelvin's circulation theorem (C16)",
        "Ch. 13: Coriolis 2Ω × u as ε_ijk (C08), rotating frames as C(t) (C02)",
    ],
    left_out=[
        "Exercises (S01) and the literature (S02)",
        "curvilinear coordinates: Ch. 3 and Appendix B (N72)",
    ],
)

path = nb.save()
print(f"wrote {path.relative_to(pathlib.Path(__file__).resolve().parents[1]).as_posix()} "
      f"({len(nb.cells)} cells, {len(nb.cores)} CORE blocks, {len(nb.derivations)} derivations, "
      f"{len(nb.primers)} primers, {len(nb.explainers)} explainers)")
