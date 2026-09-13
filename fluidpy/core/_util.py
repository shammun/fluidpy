"""Private helpers shared by the physics primitives (not part of the public API).

Physics functions in ``fluidpy`` are vectorised with numpy but must also be *scalar-callable*: an explainer's
``selftest()`` parity row calls e.g. ``ch01.capillary_rise(0.0728, 1.5708, 998.0, 1e-3)`` and expects a plain number.
``as_scalar_if_0d`` turns the 0-d arrays that numpy produces for scalar input back into Python floats.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def as_scalar_if_0d(x: Any) -> Any:
    """Return a Python float for a 0-d array / numpy scalar, otherwise the array unchanged."""
    arr = np.asarray(x)
    if arr.ndim == 0:
        if arr.dtype.kind in "biuf":
            return float(arr)
        return arr.item()
    return arr


def require_nonnegative(name: str, value: Any) -> None:
    """Raise ValueError if any element of ``value`` is negative (used for transport coefficients, second law iii)."""
    if np.any(np.asarray(value, dtype=float) < 0):
        raise ValueError(f"{name} must be >= 0 (second law: transport coefficients are positive); got {value!r}")


def require_positive(name: str, value: Any) -> None:
    """Raise ValueError if any element of ``value`` is not strictly positive."""
    if np.any(np.asarray(value, dtype=float) <= 0):
        raise ValueError(f"{name} must be > 0; got {value!r}")
