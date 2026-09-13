"""Dimensional analysis as linear algebra: dimension vectors, the dimensional matrix, its rank, and Buckingham Π groups.

Book: Kundu, Cohen & Dowling, *Fluid Mechanics* 5e, §1.11 — dimensional homogeneity, Buckingham's theorem (1.36)–(1.37),
the seven-step procedure on the pipe pressure drop (1.38)–(1.40), Examples 1.2–1.5. Reused by Ch. 4 §4.11 (Re, Fr, Ro,
...), Ch. 8, 9, 12, 13, 14, 15, 16.

The idea in one line: a monomial ``q1^a1 q2^a2 ... qn^an`` is dimensionless exactly when its exponent vector lies in the
null space of the dimensional matrix D (rows = base dimensions, columns = variables). Rank–nullity gives n − r
independent groups (Buckingham). All exponent arithmetic here is **exact** (``fractions.Fraction``), so the explainer's
JavaScript mirror (a small Fraction class) can reproduce every number.

Variables are given as an ordered mapping ``name -> dimension spec``; a spec is a pint unit string (``"Pa"``,
``"kg/m**3"``, ``"J/(kmol*K)"``), a pint Quantity/Unit, a mapping of base exponents (``{"M": 1, "L": -1, "T": -2}``)
or a sequence of exponents in the order of ``BASIS``.

Convention: the amount-of-substance dimension (kmol) is dropped with a warning, as the book does in Example 1.2
("a kmole is a pure number"). Because [substance] enters M_w and R_u only as kmol^-1, dropping it never changes the
groups of a problem that contains both.
"""
from __future__ import annotations

import itertools
import warnings
from fractions import Fraction
from typing import Any, Mapping, Sequence

import numpy as np

#: Base dimensions in the book's order (§1.11 Step 2): mass M, length L, time T, temperature Θ.
BASIS: tuple[str, ...] = ("M", "L", "T", "Θ")
_PINT_KEYS = {"M": "[mass]", "L": "[length]", "T": "[time]", "Θ": "[temperature]"}
_ALIASES = {"M": "M", "mass": "M", "L": "L", "length": "L", "T": "T", "time": "T", "Θ": "Θ", "θ": "Θ",
            "theta": "Θ", "Theta": "Θ", "temperature": "Θ"}

#: Default LaTeX symbols for the variable names used by the presets (group_latex).
LATEX_SYMBOLS: dict[str, str] = {
    "dp": r"\Delta p", "dx": r"\Delta x", "eps": r"\varepsilon", "rho": r"\rho", "mu": r"\mu", "nu": r"\nu",
    "lam": r"\lambda", "beta": r"\beta", "theta0": r"\theta_0", "tau": r"\tau", "T0": "T_0", "Mw": "M_w",
    "Ru": "R_u", "n_s": "n_s", "Lp": "L_p", "sigma": r"\sigma", "omega": r"\omega",
}


def _full_vector(q: Any, drop_substance: bool = True) -> np.ndarray:
    """Exponents over all four base dimensions (M, L, T, Θ) as float64."""
    if isinstance(q, Mapping):
        vec = np.zeros(4)
        for k, v in q.items():
            if k not in _ALIASES:
                raise ValueError(f"unknown base dimension {k!r}; use M, L, T, Θ")
            vec[BASIS.index(_ALIASES[k])] = float(v)
        return vec
    if isinstance(q, (int, float, np.integer, np.floating)):
        return np.zeros(4)
    if isinstance(q, (list, tuple, np.ndarray)):
        arr = np.asarray(q, dtype=float).ravel()
        if arr.size not in (3, 4):
            raise ValueError("a dimension sequence needs 3 (M, L, T) or 4 (M, L, T, Θ) exponents")
        return np.concatenate([arr, np.zeros(4 - arr.size)])
    from .units import ureg

    obj = ureg.parse_expression(q) if isinstance(q, str) else q
    dims = dict(getattr(obj, "dimensionality", {}))
    if "[substance]" in dims:
        if not drop_substance:
            raise ValueError(f"{q!r} has a [substance] dimension; pass drop_substance=True to treat kmol as a number")
        warnings.warn(f"{q!r}: the [substance] (mole) dimension is dropped — a kmol is treated as a pure number, as in "
                      "Example 1.2", UserWarning, stacklevel=3)
        dims.pop("[substance]")
    extra = set(dims) - set(_PINT_KEYS.values())
    if extra:
        raise ValueError(f"{q!r} involves dimensions outside M, L, T, Θ: {sorted(extra)}")
    return np.array([float(dims.get(_PINT_KEYS[b], 0.0)) for b in BASIS])


def _int_if_possible(arr: np.ndarray) -> np.ndarray:
    return arr.astype(int) if np.all(np.equal(np.mod(arr, 1), 0)) else arr


def dimension_vector(q: Any, basis: Sequence[str] = BASIS, drop_substance: bool = True) -> np.ndarray:
    """Exponents of the base dimensions in the dimension [q] of a quantity.

    Book: §1.11 Step 2 (notation [q]; base dimensions M, L, T, Θ; e.g. [p] = M/(L T^2), [C_p] = L^2/(Θ T^2)).

    Parameters
    ----------
    q : str, pint Quantity/Unit, mapping or sequence
        A pint unit expression (``"Pa"``), a pint object, a mapping ``{"M": 1, "L": -1, "T": -2}`` or exponents in
        ``BASIS`` order. Plain numbers, ``"dimensionless"`` and ``"rad"`` are dimensionless.
    basis : sequence of str, optional
        Which base dimensions to report and in which order (default ``BASIS`` = M, L, T, Θ).
    drop_substance : bool, optional
        Drop a [substance] exponent with a ``UserWarning`` (default, Example 1.2) or raise ValueError if False.

    Returns
    -------
    vec : ndarray of int (float if an exponent is fractional), shape (len(basis),)

    Raises
    ------
    ValueError
        For dimensions outside M, L, T, Θ, or a nonzero exponent of a base dimension missing from ``basis``.

    Validation (planned): V1 Pa → (1, −1, −2, 0); J/(kg K) → (0, 2, −2, −1); R_u in J/(kmol K) → (1, 2, −2, −1) as in
    the Example 1.2 matrix. Label: pending.
    """
    full = _full_vector(q, drop_substance)
    idx = [BASIS.index(_ALIASES.get(b, b)) for b in basis]
    missing = [BASIS[i] for i in range(4) if i not in idx and full[i] != 0]
    if missing:
        raise ValueError(f"{q!r} has dimensions {missing} that are not in basis {tuple(basis)}")
    return _int_if_possible(full[idx])


def dimensional_matrix(variables: Mapping[str, Any], basis: Sequence[str] = BASIS, drop_zero_rows: bool = True):
    """Dimensional matrix of a problem: one column of base-dimension exponents per variable.

    Book: §1.11 Step 2, Eq. (1.39) (pipe problem), and the matrices of Examples 1.2–1.5.

    Parameters
    ----------
    variables : mapping
        Ordered ``name -> dimension spec``, e.g. ``PIPE`` =
        ``{"dp": "Pa", "dx": "m", "d": "m", "eps": "m", "U": "m/s", "rho": "kg/m**3", "mu": "Pa*s"}``.
    basis : sequence of str, optional
        Row dimensions, default (M, L, T, Θ).
    drop_zero_rows : bool, optional
        Drop base dimensions that no variable involves (default True: the pipe problem gives the M, L, T rows of (1.39)).

    Returns
    -------
    (A, names, rows) : tuple
        ``A``: ndarray of int (float if exponents are fractional), shape (len(rows), n), column order = dict order;
        ``names``: list of variable names; ``rows``: list of row labels.

    Validation (planned): V1 pipe matrix equals (1.39); V6 Example 1.2, 1.4, 1.5 matrices (private JSON). Label: pending.
    """
    names = list(variables)
    if not names:
        raise ValueError("no variables")
    cols = np.column_stack([dimension_vector(variables[k], basis) for k in names]).astype(float)
    keep = [i for i in range(len(basis)) if (not drop_zero_rows) or np.any(cols[i] != 0)]
    A = _int_if_possible(cols[keep])
    return A, names, [basis[i] for i in keep]


# --------------------------------------------------------------------------------------------------------------------
# Exact linear algebra on small matrices (Fractions)
# --------------------------------------------------------------------------------------------------------------------
def _to_fraction_matrix(A: Any) -> list[list[Fraction]]:
    if hasattr(A, "tolist") and not isinstance(A, np.ndarray):  # sympy Matrix
        A = A.tolist()
    rows = [list(r) for r in (np.atleast_2d(np.asarray(A, dtype=object)))]
    return [[Fraction(x).limit_denominator(10**6) if isinstance(x, float) else Fraction(x) for x in r] for r in rows]


def _det_cofactor(M: list[list[Fraction]]) -> Fraction:
    """Determinant by cofactor (Laplace) expansion along the first row, skipping zero entries."""
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    total = Fraction(0)
    for j, a in enumerate(M[0]):
        if a != 0:
            sub = [row[:j] + row[j + 1:] for row in M[1:]]
            total += (-1) ** j * a * _det_cofactor(sub)
    return total


def _rank_exact(M: list[list[Fraction]]) -> int:
    """Rank by exact Gaussian elimination."""
    R = [row[:] for row in M]
    if not R:
        return 0
    m, n = len(R), len(R[0])
    rank, col = 0, 0
    while rank < m and col < n:
        piv = next((i for i in range(rank, m) if R[i][col] != 0), None)
        if piv is None:
            col += 1
            continue
        R[rank], R[piv] = R[piv], R[rank]
        for i in range(m):
            if i != rank and R[i][col] != 0:
                f = R[i][col] / R[rank][col]
                R[i] = [a - f * b for a, b in zip(R[i], R[rank])]
        rank += 1
        col += 1
    return rank


def _as_number(x: Fraction):
    return int(x) if x.denominator == 1 else x


def minor_determinant(A: Any, rows: Sequence[int], cols: Sequence[int]):
    """Exact determinant of the square sub-matrix of A formed by the given rows and columns (cofactor expansion).

    Book: §1.11 Step 3 (minors of (1.39): the first three columns give 0, the last three give −1).

    Parameters
    ----------
    A : array_like
        Matrix of integer (or rational) exponents.
    rows, cols : sequence of int
        0-based indices, same length.

    Returns
    -------
    det : int (Fraction if not integral)

    Validation (planned): V1 pipe minors 0 and −1; V2 agrees with numpy.linalg.det. Label: pending.
    """
    if len(rows) != len(cols):
        raise ValueError("a minor needs as many rows as columns")
    M = _to_fraction_matrix(A)
    sub = [[M[i][j] for j in cols] for i in rows]
    return _as_number(_det_cofactor(sub))


def rank_by_minors(A: Any):
    """Rank as the size of the largest square sub-matrix with a nonzero determinant, with a witness minor.

    Book: §1.11 Step 3 ("the rank r of any matrix is defined to be the size of the largest square submatrix that has a
    nonzero determinant").

    Parameters
    ----------
    A : array_like
        Dimensional matrix.

    Returns
    -------
    (r, rows, cols) : tuple
        Rank r and the 0-based row and column indices of the first nonzero r × r minor (sizes searched from min(m, n)
        downward, index combinations in lexicographic order); ``(0, (), ())`` for a zero matrix.

    Notes
    -----
    Method: exhaustive exact cofactor determinants over ``itertools.combinations`` — the transparent method the book
    teaches, fine for the n <= ~10 of dimensional analysis. Cross-check with ``numpy.linalg.matrix_rank``.

    Validation (planned): V1 ranks 3, 4, 1, 3, 2 for the pipe problem and Examples 1.2–1.5; the book's dependent-row
    example has rank 2; V2 agrees with numpy on 200 random integer matrices (seed 0). Label: pending.
    """
    M = _to_fraction_matrix(A)
    m, n = len(M), len(M[0])
    for k in range(min(m, n), 0, -1):
        for rows in itertools.combinations(range(m), k):
            for cols in itertools.combinations(range(n), k):
                if _det_cofactor([[M[i][j] for j in cols] for i in rows]) != 0:
                    return k, rows, cols
    return 0, (), ()


def solve_exponents(target: str, repeating: Sequence[str], variables: Mapping[str, Any]):
    """Exponents that make ``target · Π repeating_j^{a_j}`` dimensionless (exponent algebra).

    Book: §1.11 Step 5, exponent algebra: ``Π1 = Δp U^a d^b rho^c`` with
    ``M^{c+1} L^{a+b−3c−1} T^{−a−2} = M^0 L^0 T^0`` giving a = −2, b = 0, c = −1.

    Parameters
    ----------
    target : str
        Name of the non-repeating variable, raised to the first power.
    repeating : sequence of str
        The r repeating variables; their dimensional sub-matrix must have full column rank.
    variables : mapping
        All variables of the problem (``name -> dimension spec``).

    Returns
    -------
    group : dict[str, Fraction]
        ``{target: 1, repeating_1: a_1, ...}`` (zero exponents kept, so the book's b = 0 is visible).

    Raises
    ------
    ValueError
        ``"singular repeating set"`` if the repeating columns are dependent, or if they cannot cancel the target's
        dimensions (inconsistent system).

    Notes
    -----
    Method: exact Gauss–Jordan elimination on ``R a = −d_target`` (R = repeating columns), Fractions throughout.

    Validation (planned): V1 pipe (−2, 0, −1) with (U, d, rho); Example 1.2 (−1, 1, 1, −1); every result has a zero
    dimension vector; a singular repeating set raises. Label: pending.
    """
    repeating = list(repeating)
    sub = {k: variables[k] for k in [target] + repeating}
    A, names, _ = dimensional_matrix(sub, BASIS, drop_zero_rows=True)
    M = _to_fraction_matrix(A)
    m = len(M)
    r = len(repeating)
    R = [[M[i][names.index(k)] for k in repeating] for i in range(m)]
    if _rank_exact(R) != r:
        raise ValueError(f"singular repeating set {tuple(repeating)}: its dimensional sub-matrix has rank "
                         f"{_rank_exact(R)} < {r}")
    b = [-M[i][names.index(target)] for i in range(m)]
    aug = [R[i] + [b[i]] for i in range(m)]
    # Gauss–Jordan
    row = 0
    pivots = []
    for col in range(r):
        piv = next((i for i in range(row, m) if aug[i][col] != 0), None)
        if piv is None:
            continue
        aug[row], aug[piv] = aug[piv], aug[row]
        p = aug[row][col]
        aug[row] = [x / p for x in aug[row]]
        for i in range(m):
            if i != row and aug[i][col] != 0:
                f = aug[i][col]
                aug[i] = [x - f * y for x, y in zip(aug[i], aug[row])]
        pivots.append(col)
        row += 1
    if any(aug[i][r] != 0 for i in range(row, m)):
        raise ValueError(f"repeating set {tuple(repeating)} cannot cancel the dimensions of {target!r}")
    group: dict[str, Fraction] = {target: Fraction(1)}
    for i, col in enumerate(pivots):
        group[repeating[col]] = aug[i][r]  # R a = −d_target, solved exactly
    return group


def _auto_repeating(variables: Mapping[str, Any], exclude: Sequence[str], r: int) -> tuple[str, ...]:
    """Pick r repeating variables: prefer a velocity, then lengths, then mass-bearing properties (book, Step 5)."""
    vecs = {k: _full_vector(v) for k, v in variables.items()}
    cand = [k for k in variables if k not in exclude and np.any(vecs[k] != 0)]

    def pref(k):
        v = vecs[k]
        if np.allclose(v, [0, 1, -1, 0]):
            return 0  # velocity
        if np.allclose(v, [0, 1, 0, 0]):
            return 1  # length
        if v[0] != 0:
            return 2  # involves mass (density, viscosity, ...)
        return 3

    order = sorted(cand, key=lambda k: (pref(k), list(variables).index(k)))
    for combo in itertools.combinations(order, r):
        cols = np.column_stack([vecs[k] for k in combo])
        if _rank_exact(_to_fraction_matrix(cols)) == r:
            return combo
    raise ValueError("no repeating set of full rank exists among the non-solution variables")


def pi_groups(variables: Mapping[str, Any], solution: str | None = None, repeating: Sequence[str] | None = None):
    """A complete set of n − r independent dimensionless groups (Buckingham Π theorem).

    Book: §1.11, Eq. (1.37) ``φ(Π1, …, Π_{n−r}) = 0`` from Eq. (1.36) ``f(q1, …, qn) = 0``; Steps 1–6; (1.40) for the
    pipe with ``repeating=("U", "d", "rho")``.

    Parameters
    ----------
    variables : mapping
        Ordered ``name -> dimension spec`` of all n variables and parameters (Step 1).
    solution : str, optional
        The solution variable; it appears to the first power in the first group. Default: the first variable.
    repeating : sequence of str, optional
        r repeating variables (Step 5). Default: chosen automatically (prefer a velocity, then lengths, then a
        mass-bearing property; first full-rank combination). Pass the book's choice to reproduce its groups.

    Returns
    -------
    groups : list of dict[str, Fraction]
        n − r groups. The first contains ``solution``; each other non-repeating variable heads one group (in the
        order of ``variables``); a dimensionless variable is its own group.

    Notes
    -----
    The set is not unique (another repeating set gives another basis of the same null space). Assumptions: the
    variable list is complete and the relation is dimensionally homogeneous (D28).

    Validation (planned): V2 count n − r, every group dimensionless, groups independent; V4 group values invariant
    under a change of units (:func:`rescale_units`); V1 Example 1.3 gives a/C^2 and β. Label: pending.
    """
    names = list(variables)
    solution = names[0] if solution is None else solution
    if solution not in variables:
        raise ValueError(f"solution variable {solution!r} is not in the variable list")
    A, _, _ = dimensional_matrix(variables)
    r = _rank_exact(_to_fraction_matrix(A))
    if repeating is None:
        repeating = _auto_repeating(variables, exclude=[solution], r=r) if r > 0 else ()
    repeating = tuple(repeating)
    if len(repeating) != r:
        raise ValueError(f"need exactly r = {r} repeating variables, got {len(repeating)}")
    if solution in repeating:
        raise ValueError("the solution variable must not be a repeating variable")
    others = [k for k in names if k not in repeating and k != solution]
    groups = []
    for k in [solution] + others:
        if np.all(_full_vector(variables[k]) == 0):
            groups.append({k: Fraction(1)})  # already dimensionless
        else:
            groups.append(solve_exponents(k, repeating, variables))
    return groups


def group_dimension(group: Mapping[str, Any], variables: Mapping[str, Any]) -> np.ndarray:
    """Dimension vector (M, L, T, Θ) of a monomial group — zero for a true Π group.

    Book: §1.11 (dimensional homogeneity; the final check of Step 5).

    Parameters
    ----------
    group : mapping
        ``name -> exponent``.
    variables : mapping
        ``name -> dimension spec``.

    Returns
    -------
    vec : ndarray, shape (4,)

    Validation (planned): V1 zero for every output of :func:`pi_groups`. Label: pending.
    """
    return sum((float(e) * _full_vector(variables[k]) for k, e in group.items()), np.zeros(4))


def exponent_matrix(groups: Sequence[Mapping[str, Any]], names: Sequence[str]) -> np.ndarray:
    """Exponent matrix of a set of groups: rows = groups, columns = variable names.

    Book: §1.11 (combining groups; independence is the rank of this matrix).

    Parameters
    ----------
    groups : sequence of mappings
        ``name -> exponent``.
    names : sequence of str
        Column order.

    Returns
    -------
    E : ndarray of float, shape (len(groups), len(names))

    Validation (planned): V1 pipe groups give a 4 × 7 matrix of rank 4. Label: pending.
    """
    return np.array([[float(g.get(k, 0)) for k in names] for g in groups], dtype=float)


def groups_independent(groups: Sequence[Mapping[str, Any]], names: Sequence[str] | None = None) -> bool:
    """True if the groups are independent: their exponent vectors are linearly independent.

    Book: §1.11 Step 5 (combining groups: Δp d^2 rho/μ^2 = Π1/Π4^2 and ε/Δx = Π3/Π2 are groups, but only n − r are
    independent).

    Parameters
    ----------
    groups : sequence of mappings
        ``name -> exponent`` for each group.
    names : sequence of str, optional
        Column order of the exponent matrix; default the union of the groups' keys.

    Returns
    -------
    independent : bool
        ``rank(E) == len(groups)`` with exact arithmetic.

    Validation (planned): V1 {Π1..Π4} independent; {Π1, Π4, Π1/Π4^2} dependent. Label: pending.
    """
    if names is None:
        names = []
        for g in groups:
            for k in g:
                if k not in names:
                    names.append(k)
    E = [[Fraction(g.get(k, 0)) for k in names] for g in groups]
    return _rank_exact(E) == len(groups)


def _latex_power(sym: str, e: Fraction) -> str:
    if e == 1:
        return sym
    txt = str(e.numerator) if e.denominator == 1 else f"{e.numerator}/{e.denominator}"
    return f"{sym}^{{{txt}}}"


def group_latex(group: Mapping[str, Any], symbols: Mapping[str, str] | None = None) -> str:
    """LaTeX monomial of a group, e.g. ``{"dp": 1, "U": -2, "rho": -1}`` → ``\\frac{\\Delta p}{U^{2} \\rho}``.

    Book: §1.11 (Π groups as printed in (1.40) and Examples 1.2–1.5).

    Parameters
    ----------
    group : mapping
        ``name -> exponent`` (zero exponents omitted).
    symbols : mapping, optional
        LaTeX symbol per name; default ``LATEX_SYMBOLS`` then the name itself.

    Returns
    -------
    tex : str

    Validation (planned): V1 pipe Π1 string. Label: pending.
    """
    symbols = {**LATEX_SYMBOLS, **(symbols or {})}
    num, den = [], []
    for k, e in group.items():
        e = Fraction(e).limit_denominator(1000) if isinstance(e, float) else Fraction(e)
        if e > 0:
            num.append(_latex_power(symbols.get(k, k), e))
        elif e < 0:
            den.append(_latex_power(symbols.get(k, k), -e))
    top = " ".join(num) if num else "1"
    return rf"\frac{{{top}}}{{{' '.join(den)}}}" if den else top


def group_expression(group: Mapping[str, Any], labels: Mapping[str, str] | None = None):
    """Sympy expression of a group, e.g. ``dp/(U**2*rho)`` (sympy imported lazily).

    Parameters
    ----------
    group : mapping
        ``name -> exponent``.
    labels : mapping, optional
        Symbol names to use instead of the variable names.

    Returns
    -------
    expr : sympy.Expr

    Validation (planned): V1 pipe Π1. Label: pending.
    """
    import sympy as sp

    labels = labels or {}
    expr = sp.Integer(1)
    for k, e in group.items():
        e = sp.Rational(Fraction(e).numerator, Fraction(e).denominator)
        if e != 0:
            expr *= sp.Symbol(labels.get(k, k), positive=True) ** e
    return expr


def group_value(group: Mapping[str, Any], values: Mapping[str, float]) -> float:
    """Numerical value of a group: ``Π = Π_k values[k]^{exponent_k}``.

    Book: §1.11 (a Π group is a number, the same in every unit system).

    Parameters
    ----------
    group : mapping
        ``name -> exponent``.
    values : mapping
        ``name -> value`` in any **consistent** unit system (SI by default).

    Returns
    -------
    value : float

    Validation (planned): V4 unchanged by :func:`rescale_units`. Label: pending.
    """
    out = 1.0
    for k, e in group.items():
        out *= float(values[k]) ** float(e)
    return out


#: Base-unit scale factors: numeric value of one SI base unit (kg, m, s, K) in the target system's base units.
UNIT_SYSTEMS: dict[str, dict[str, float]] = {
    "SI": {"M": 1.0, "L": 1.0, "T": 1.0, "Θ": 1.0},
    "cgs": {"M": 1.0e3, "L": 1.0e2, "T": 1.0, "Θ": 1.0},                        # g, cm, s, K
    "imperial": {"M": 1.0 / 0.45359237, "L": 1.0 / 0.3048, "T": 1.0, "Θ": 1.8},  # lb, ft, s, °R (exact definitions)
}


def rescale_units(values: Mapping[str, float], variables: Mapping[str, Any], system: str | Mapping[str, float] = "cgs"):
    """Express the same physical values as numbers in another consistent system of base units.

    Book: §1.11 opening argument (natural laws do not depend on units: cgs, MKS or English), used to show that every
    Π group keeps its value when the units change.

    Parameters
    ----------
    values : mapping
        ``name -> SI magnitude``.
    variables : mapping
        ``name -> dimension spec``.
    system : {"SI", "cgs", "imperial"} or mapping, optional
        Target system (g-cm-s-K; lb-ft-s-°R) or scale factors ``{"M", "L", "T", "Θ"}`` = new units per SI unit.

    Returns
    -------
    new_values : dict[str, float]
        ``value · Π_i scale_i^{exponent_i}``.

    Notes
    -----
    Method: base-unit scale factors applied to the dimension vector (equivalent to pint conversion into coherent
    derived units of the target system; the pound here is the pound-mass).

    Validation (planned): V1 1 Pa = 10 dyn/cm^2 in cgs; V4 :func:`group_value` unchanged. Label: pending.
    """
    scales = UNIT_SYSTEMS[system] if isinstance(system, str) else {_ALIASES.get(k, k): v for k, v in system.items()}
    s = np.array([scales[b] for b in BASIS], dtype=float)
    return {k: float(values[k]) * float(np.prod(s ** _full_vector(variables[k]))) for k in values}


# --------------------------------------------------------------------------------------------------------------------
# Problem presets: variable name -> SI unit (keys are the code names used by every notebook cell and explainer)
# --------------------------------------------------------------------------------------------------------------------
PIPE = {"dp": "Pa", "dx": "m", "d": "m", "eps": "m", "U": "m/s", "rho": "kg/m**3", "mu": "Pa*s"}  # Eq. (1.38)
PENDULUM = {"tau": "s", "Lp": "m", "m": "kg", "g": "m/s**2", "theta0": "rad"}
SPHERE_DRAG = {"F": "N", "D": "m", "U": "m/s", "rho": "kg/m**3", "mu": "Pa*s"}
SCALE_HEIGHT = {"H": "m", "T0": "K", "Mw": "kg/kmol", "g": "m/s**2", "Ru": "J/(kmol*K)"}  # Example 1.2
PYTHAGORAS = {"a": "m**2", "beta": "rad", "C": "m"}  # Example 1.3
BLAST = {"E": "J", "D": "m", "rho": "kg/m**3", "t": "s"}  # Example 1.4
RAYLEIGH = {"S": "W/m**2", "I": "W/m**2", "lam": "m", "V": "m**3", "n_s": "dimensionless", "d": "m"}  # Example 1.5

#: Solution variable, the repeating set that reproduces the book's (or the classic) groups, and a title per preset.
PRESET_INFO: dict[str, dict[str, Any]] = {
    "pipe": {"variables": PIPE, "solution": "dp", "repeating": ("U", "d", "rho"),
             "title": "Pressure drop in a round pipe (§1.11, Eqs. 1.38–1.40)"},
    "pendulum": {"variables": PENDULUM, "solution": "tau", "repeating": ("Lp", "m", "g"),
                 "title": "Period of a simple pendulum"},
    "sphere_drag": {"variables": SPHERE_DRAG, "solution": "F", "repeating": ("D", "U", "rho"),
                    "title": "Drag force on a sphere"},
    "scale_height": {"variables": SCALE_HEIGHT, "solution": "H", "repeating": ("T0", "Mw", "g", "Ru"),
                     "title": "Scale height of an isothermal atmosphere (Example 1.2)"},
    "pythagoras": {"variables": PYTHAGORAS, "solution": "a", "repeating": ("C",),
                   "title": "Area of a right triangle (Example 1.3)"},
    "blast": {"variables": BLAST, "solution": "E", "repeating": ("D", "rho", "t"),
              "title": "Energy of a point blast (Example 1.4)"},
    "rayleigh": {"variables": RAYLEIGH, "solution": "S", "repeating": ("I", "lam"),
                 "title": "Light scattered by a small particle (Example 1.5)"},
}


__all__ = [
    "BASIS", "LATEX_SYMBOLS", "dimension_vector", "dimensional_matrix", "minor_determinant", "rank_by_minors",
    "solve_exponents", "pi_groups", "group_dimension", "exponent_matrix", "group_latex", "group_expression",
    "group_value", "groups_independent", "UNIT_SYSTEMS", "rescale_units",
    "PIPE", "PENDULUM", "SPHERE_DRAG", "SCALE_HEIGHT", "PYTHAGORAS", "BLAST", "RAYLEIGH", "PRESET_INFO",
]
