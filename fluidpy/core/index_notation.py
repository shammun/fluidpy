"""The summation convention as a machine: write an index expression as a string and get the implied sums back.

Book: Ch. 2 §2.1 (Einstein summation convention, Eq. (2.2)), §2.2 (free vs dummy indices, (2.5)–(2.8)), §2.3
((2.9)), §2.7 (δ_ij (2.16)–(2.17), ε_ijk (2.18)–(2.19)), §2.8 ((2.21)), §2.14 (comma notation (2.36)).

Grammar of a term (factors separated by spaces or ``*``; a leading ``-`` or a numeric factor is a coefficient)::

    a_i b_i              → a_1 b_1 + a_2 b_2 + a_3 b_3                        (2.2)
    A_ik B_kj            → Array P[i, j] = Σ_k A_ik B_kj                      (2.9)
    C_im C_jn tau_ij     → Array τ'[m, n] = Σ_ij C_im C_jn τ_ij               (2.12)
    delta_ij u_j         → Array u[i]                                          (2.17)
    eps_ijk u_i v_j      → Array (u × v)[k]                                    (2.21)
    eps_ijk eps_klm      → Array [i, j, l, m]  (compare with δ_il δ_jm − δ_im δ_jl)  (2.19)
    u_i,i                → Σ_i ∂u_i/∂x_i  (Derivative objects)                (2.36)
    eps_ijk u_k,j        → Array (∇×u)[i]                                      (2.36)
    -1/2 eps_ijk R_ij    → Array ω[k]                                          (2.27)
    A_ij + B_ij          → sums of terms with the same free indices

Rules enforced: an index letter may appear once (free) or twice (dummy) in a term — three times raises
``ValueError``; every term of a sum must carry the same free indices. Digits in a subscript are fixed index values
(``eps_ij1 u_i v_j`` is the k = 1 check after (2.21)). Recognised special names: ``delta`` (Kronecker), ``eps`` /
``epsilon`` (alternating tensor). Component symbols are ``Symbol("a_1")``, ``Symbol("tau_12")``; fields with a
comma derivative become ``Function("u_1")(x1, x2, x3)`` — build matching sympy objects with :func:`symbol_vector`,
:func:`symbol_matrix`, :func:`coordinates`, :func:`field_functions`.

Free indices of the result are ordered by **first appearance** in the term (see :func:`classify_indices`).
"""
from __future__ import annotations

import itertools
import re
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
import sympy as sp

__all__ = ["expand_indices", "expand_indices_str", "classify_indices", "tensor_order", "rename_dummy",
           "comma_to_partial", "coordinates", "symbol_vector", "symbol_matrix", "symbol_tensor", "field_functions",
           "SPECIAL_NAMES"]

SPECIAL_NAMES = {"delta": "kronecker", "eps": "levi_civita", "epsilon": "levi_civita"}

_FACTOR = re.compile(r"^(?P<name>[A-Za-z][A-Za-z0-9]*?)(?:_(?P<idx>[A-Za-z0-9]*))?(?:,(?P<didx>[A-Za-z0-9]+))?$")
_NUMBER = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)(/\d+)?$")


@dataclass
class _Factor:
    name: str
    idx: str  # subscript indices (letters or digits), '' for a scalar
    didx: str  # comma-derivative indices
    coeff: Fraction | float | None = None  # set for numeric factors


def _tokenize_term(term: str) -> list[str]:
    return [t for t in term.replace("*", " ").split() if t]


def _parse_term(term: str) -> tuple[list[_Factor], sp.Expr]:
    tokens = _tokenize_term(term)
    coeff = sp.Integer(1)
    factors: list[_Factor] = []
    for tok in tokens:
        if tok == "-":
            coeff = -coeff
            continue
        if _NUMBER.match(tok):
            coeff *= sp.Rational(tok) if "." not in tok else sp.Float(tok)
            continue
        neg = tok.startswith("-")
        if neg:
            tok = tok[1:]
            coeff = -coeff
        m = _FACTOR.match(tok)
        if not m:
            raise ValueError(f"cannot parse factor {tok!r} (expected name_indices[,derivindices])")
        factors.append(_Factor(m["name"], m["idx"] or "", m["didx"] or ""))
    if not factors:
        raise ValueError(f"term {term!r} has no tensor factors")
    return factors, coeff


def _split_terms(expr: str) -> list[tuple[int, str]]:
    parts = re.split(r"\s+([+-])\s+", expr.strip())
    out = [(1, parts[0])]
    for sign, term in zip(parts[1::2], parts[2::2]):
        out.append((1 if sign == "+" else -1, term))
    return out


def _letters(factors: list[_Factor]) -> list[str]:
    seq: list[str] = []
    for f in factors:
        seq += [c for c in f.idx + f.didx if c.isalpha()]
    return seq


def classify_indices(term: str) -> tuple[list[str], list[str]]:
    """Free and dummy (summed) indices of an index expression, each in order of first appearance.

    Book: §2.2 ("the free or not-summed-over index is j, while the repeated or summed-over index can be any letter
    other than j").

    Examples
    --------
    ``classify_indices("x_i C_ij") == (["j"], ["i"])``; ``classify_indices("A_ij B_kl") == (["i","j","k","l"], [])``.

    Raises ``ValueError`` if an index appears more than twice or the terms of a sum have different free indices.
    """
    result = None
    for _, t in _split_terms(term):
        factors, _ = _parse_term(t)
        seq = _letters(factors)
        counts = {c: seq.count(c) for c in seq}
        bad = [c for c, k in counts.items() if k > 2]
        if bad:
            raise ValueError(f"index {bad[0]!r} appears {counts[bad[0]]} times in {t!r}: the summation convention allows "
                             "one (free) or two (summed) occurrences")
        free = [c for c in dict.fromkeys(seq) if counts[c] == 1]
        dummy = [c for c in dict.fromkeys(seq) if counts[c] == 2]
        if result is None:
            result = (free, dummy)
        elif set(result[0]) != set(free):
            raise ValueError(f"terms of {term!r} carry different free indices {result[0]} vs {free}")
    return result  # type: ignore[return-value]


def tensor_order(term: str) -> int:
    """Order of the tensor an expression represents = number of free indices (§2.1: "the number of indices or
    subscripts clearly specifies the order of a tensor"). ``tensor_order("A_ij B_kl") == 4``, ``("a_i b_i") == 0``."""
    return len(classify_indices(term)[0])


def rename_dummy(term: str, old: str, new: str) -> str:
    """Rename a summed index letter throughout a term (allowed: x_i C_ij = x_k C_kj, §2.2 / Eq. (2.6)).

    Raises ``ValueError`` if ``old`` is not a dummy index or ``new`` already occurs.
    """
    free, dummy = classify_indices(term)
    if old not in dummy:
        raise ValueError(f"{old!r} is not a summed (dummy) index of {term!r}; free indices cannot be renamed alone")
    if new in free or new in dummy:
        raise ValueError(f"{new!r} already occurs in {term!r}")
    out = []
    for sign, t in _split_terms(term):
        toks = []
        for tok in _tokenize_term(t):
            if _NUMBER.match(tok) or tok == "-":
                toks.append(tok)
                continue
            m = _FACTOR.match(tok.lstrip("-"))
            if m is None:
                toks.append(tok)
                continue
            head = "-" if tok.startswith("-") else ""
            idx = (m["idx"] or "").replace(old, new)
            didx = (m["didx"] or "").replace(old, new)
            s = head + m["name"] + (f"_{idx}" if m["idx"] is not None else "") + (f",{didx}" if didx else "")
            toks.append(s)
        out.append((sign, " ".join(toks)))
    text = out[0][1]
    for sign, t in out[1:]:
        text += (" + " if sign > 0 else " - ") + t
    return text


def comma_to_partial(term: str, latex: bool = False) -> str:
    """Rewrite the comma-derivative factors of a term as explicit partial derivatives (2.36): ``u_i,i`` → ``∂u_i/∂x_i``.

    Book: §2.14, Eq. (2.36) A_,i ≡ ∂A/∂x_i. With ``latex=True`` returns ``\\frac{\\partial u_i}{\\partial x_i}``.
    """
    pieces = []
    for sign, t in _split_terms(term):
        factors, coeff = _parse_term(t)
        s = []
        if coeff != 1:
            s.append(str(coeff))
        for f in factors:
            base = f.name + (f"_{f.idx}" if f.idx else "")
            if f.didx:
                if latex:
                    den = " ".join(f"\\partial x_{c}" for c in f.didx)
                    s.append(f"\\frac{{\\partial{'^' + str(len(f.didx)) if len(f.didx) > 1 else ''} {base}}}{{{den}}}")
                else:
                    den = "".join(f"∂x_{c}" for c in f.didx)
                    s.append(f"∂{'^' + str(len(f.didx)) if len(f.didx) > 1 else ''}{base}/{den}")
            else:
                s.append(base)
        pieces.append((sign, " ".join(s)))
    text = pieces[0][1]
    for sign, t in pieces[1:]:
        text += (" + " if sign > 0 else " - ") + t
    return text


# ----------------------------------------------------------------------------------------------------------------------
# sympy building blocks shared with tests and notebooks
# ----------------------------------------------------------------------------------------------------------------------
def coordinates(dim: int = 3) -> tuple[sp.Symbol, ...]:
    """The coordinate symbols x1, x2[, x3] used for comma derivatives (``Symbol("x1")`` …)."""
    return sp.symbols(f"x1:{dim + 1}", real=True)


def _component_symbol(name: str, idx: str) -> sp.Symbol:
    return sp.Symbol(f"{name}_{idx}" if idx else name, real=True)


def symbol_vector(name: str, dim: int = 3) -> sp.Matrix:
    """Column matrix of component symbols ``name_1 … name_dim`` (the same names ``expand_indices`` produces)."""
    return sp.Matrix([_component_symbol(name, str(i)) for i in range(1, dim + 1)])


def symbol_matrix(name: str, dim: int = 3) -> sp.Matrix:
    """Matrix of component symbols ``name_ij`` (e.g. ``symbol_matrix("C")[0, 1] == Symbol("C_12")``)."""
    return sp.Matrix(dim, dim, lambda i, j: _component_symbol(name, f"{i + 1}{j + 1}"))


def symbol_tensor(name: str, order: int, dim: int = 3) -> sp.Array:
    """N-dimensional ``sympy.Array`` of component symbols ``name_ijk…`` for any order."""
    shape = (dim,) * order
    data = [_component_symbol(name, "".join(str(k + 1) for k in idx)) for idx in itertools.product(range(dim), repeat=order)]
    return sp.Array(data, shape) if order else sp.Array([_component_symbol(name, "")])[0]


def field_functions(name: str, order: int, dim: int = 3):
    """Undefined functions of position for a field: order 0 → ``name(x1, x2, x3)``; order 1 → list ``name_i(x…)``;
    order 2 → nested list ``name_ij(x…)`` — matching the objects ``expand_indices`` builds for comma derivatives."""
    X = coordinates(dim)
    if order == 0:
        return sp.Function(name)(*X)
    if order == 1:
        return [sp.Function(f"{name}_{i}")(*X) for i in range(1, dim + 1)]
    if order == 2:
        return [[sp.Function(f"{name}_{i}{j}")(*X) for j in range(1, dim + 1)] for i in range(1, dim + 1)]
    raise ValueError("field_functions supports orders 0, 1, 2")


# ----------------------------------------------------------------------------------------------------------------------
# the expander
# ----------------------------------------------------------------------------------------------------------------------
def _factor_value(f: _Factor, assign: dict[str, int], dim: int, values: dict | None):
    """Numeric/symbolic value of one factor under an index assignment (indices as 1-based ints)."""
    ints = [assign[c] if c.isalpha() else int(c) for c in f.idx]
    kind = SPECIAL_NAMES.get(f.name)
    if kind == "kronecker":
        if len(ints) != 2:
            raise ValueError("delta needs exactly two indices")
        return sp.Integer(1 if ints[0] == ints[1] else 0)  # Eq. (2.16)
    if kind == "levi_civita":
        if len(ints) != 3:
            raise ValueError("eps needs exactly three indices")
        return sp.LeviCivita(*ints)  # Eq. (2.18)
    if values is not None and f.name in values:
        arr = np.asarray(values[f.name])
        v = arr[tuple(i - 1 for i in ints)] if ints else arr[()]
        v = float(v)
        return sp.Integer(int(v)) if v.is_integer() else sp.Float(v)
    idx_str = "".join(str(i) for i in ints)
    if f.didx:
        X = coordinates(dim)
        base = sp.Function(f"{f.name}_{idx_str}" if idx_str else f.name)(*X)
        dints = [assign[c] if c.isalpha() else int(c) for c in f.didx]
        return sp.Derivative(base, *[X[d - 1] for d in dints])  # Eq. (2.36): A_,i = ∂A/∂x_i
    return _component_symbol(f.name, idx_str)


def expand_indices(term: str, dim: int = 3, values: dict | None = None, free_order: str | None = None):
    """Write out the sums the summation convention hides.

    Book: §2.1 ("Whenever an index is repeated in a term, a summation over this index is implied"); Eqs. (2.2), (2.5),
    (2.9), (2.12), (2.17), (2.19), (2.21), (2.27), (2.36) are all instances.

    Parameters
    ----------
    term : str
        Index expression (grammar in the module docstring), e.g. ``"a_i b_i"``, ``"C_im C_jn tau_ij"``, ``"u_i,i"``.
    dim : int
        Range of every index (3 in the book; 2 for plane examples).
    values : dict, optional
        ``{"a": [1, 2, 3], "tau": [[…]]}`` — numeric components substituted for the named factors (numpy indexing, the
        book's 1-based subscripts are shifted internally); factors not listed stay symbolic.
    free_order : str, optional
        The order of the free indices in the result, e.g. ``"ijlm"``; default = order of first appearance
        (``classify_indices(term)[0]``). Needed when comparing two expansions written in different orders, such as
        ``"eps_ijk eps_klm"`` with ``"delta_il delta_jm - delta_im delta_jl"``.

    Returns
    -------
    sympy.Expr if there are no free indices; otherwise ``sympy.Array`` of shape ``(dim,)*n_free`` indexed by the free
    indices in ``free_order``.

    Validation: V1 ``expand_indices("a_i b_i")`` == a_1 b_1 + a_2 b_2 + a_3 b_3; ``"delta_ij u_j"`` == u; V2
    ``"C_im C_jn tau_ij"`` == (Cᵀ τ C)_mn for all m, n; ``"eps_ijk u_i v_j"`` == sympy ``cross``; ``"eps_ijk eps_klm"``
    == δ_il δ_jm − δ_im δ_jl (81 cases); ``"u_i,i"`` == Σ ∂u_i/∂x_i; three repeats → ValueError. Label: symbolic.
    """
    free, dummy = classify_indices(term)
    if free_order is not None:
        if sorted(free_order) != sorted(free):
            raise ValueError(f"free_order {free_order!r} must be a permutation of the free indices {free}")
        free = list(free_order)
    terms = _split_terms(term)
    parsed = [(sign, *_parse_term(t)) for sign, t in terms]

    def value_at(free_assign: dict[str, int]):
        total = sp.Integer(0)
        for sign, factors, coeff in parsed:
            for dvals in itertools.product(range(1, dim + 1), repeat=len(dummy)):
                assign = dict(free_assign)
                assign.update(zip(dummy, dvals))
                prod = sign * coeff
                for f in factors:
                    prod = prod * _factor_value(f, assign, dim, values)
                total += prod
        return total

    if not free:
        return sp.expand(value_at({}))
    shape = (dim,) * len(free)
    data = [sp.expand(value_at(dict(zip(free, [i + 1 for i in idx])))) for idx in itertools.product(range(dim), repeat=len(free))]
    return sp.Array(data, shape)


def expand_indices_str(term: str, dim: int = 3, values: dict | None = None, lhs: str = "result",
                       free_order: str | None = None) -> str:
    """Plain-text form of :func:`expand_indices` for explainer parity rows and printing: one line per free-index
    combination (``lhs_12 = …``), or the single expanded sum when there are no free indices."""
    free, _ = classify_indices(term)
    free = list(free_order) if free_order else free
    res = expand_indices(term, dim, values, free_order)
    if not free:
        return str(res)
    lines = []
    for idx in itertools.product(range(dim), repeat=len(free)):
        lines.append(f"{lhs}_{''.join(str(i + 1) for i in idx)} = {res[idx]}")
    return "\n".join(lines)
