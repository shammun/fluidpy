"""Units and dimensional-homogeneity checking (seed file — extend as the book demands).

Dimensional analysis is the cheapest bug detector in this project: it catches nu-vs-mu, a missing rho, a missing 1/2 and
a gauge-vs-absolute pressure before any physics is compared. Use it in tests (V2 of the verify-implementation ladder)
and at script/notebook boundaries; keep the hot numerical paths as plain floats.

    from fluidpy.core.units import ureg, Q_, check_dimensions

    q = Q_(0.5, "") * Q_(998, "kg/m**3") * Q_(2, "m/s") ** 2
    assert q.check("[pressure]")
"""
from __future__ import annotations

from typing import Any, Callable

try:
    import pint
except ImportError as exc:  # pragma: no cover
    raise ImportError("pint is required: pip install pint") from exc

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# Common dimensionalities, so tests can say what they mean.
DIM = {
    "velocity": "[length] / [time]",
    "acceleration": "[length] / [time] ** 2",
    "density": "[mass] / [length] ** 3",
    "pressure": "[mass] / ([length] * [time] ** 2)",
    "dynamic_viscosity": "[mass] / ([length] * [time])",
    "kinematic_viscosity": "[length] ** 2 / [time]",
    "force": "[mass] * [length] / [time] ** 2",
    "energy": "[mass] * [length] ** 2 / [time] ** 2",
    "circulation": "[length] ** 2 / [time]",
    "vorticity": "1 / [time]",
    "mass_flow": "[mass] / [time]",
    "volume_flow": "[length] ** 3 / [time]",
}


def check_dimensions(value: Any, expected: str) -> bool:
    """True if ``value`` (a pint Quantity) has the dimensionality named in ``DIM`` or given as a pint expression."""
    dim = DIM.get(expected, expected)
    return bool(value.check(dim))


def assert_dimensions(value: Any, expected: str, name: str = "") -> None:
    """Raise AssertionError with a useful message if the dimensionality is wrong."""
    if not check_dimensions(value, expected):
        raise AssertionError(f"{name or 'value'} has dimensionality {value.dimensionality}, expected {expected}")


def dimensional_check(fn: Callable[..., Any], expected: str, **kwargs: Any) -> Any:
    """Call ``fn`` with pint Quantities and assert the result's dimensionality; returns the result."""
    out = fn(**kwargs)
    assert_dimensions(out, expected, getattr(fn, "__name__", "result"))
    return out


#: Offset between the Celsius and kelvin scales [K] (exact, SI definition of the degree Celsius).
KELVIN_OFFSET = 273.15


def celsius_to_kelvin(T_C: Any) -> Any:
    """Convert a temperature from degrees Celsius to kelvin.

    Book: §1.2 (the relation printed as "°C = K − 273.15", i.e. T[°C] = T[K] − 273.15).

    Parameters
    ----------
    T_C : float or array_like
        Temperature on the Celsius scale [°C].

    Returns
    -------
    T_K : float or ndarray
        Absolute temperature [K].

    Notes
    -----
    Assumptions: this converts a temperature *value*. A temperature *difference* has the same number in °C and K
    (pint: ``delta_degC``); never pass a difference through this function. Every ``fluidpy`` function takes kelvin.

    Validation: V1 round trip °C → K → °C (1e-12), 0 °C = 273.15 K returned as a float; agrees with pint
    ``Q_(T, "degC").to("K")`` (1e-12). Label: analytic.
    """
    import numpy as np

    out = np.asarray(T_C, dtype=float) + KELVIN_OFFSET  # §1.2: T[K] = T[°C] + 273.15
    return float(out) if out.ndim == 0 else out


def kelvin_to_celsius(T_K: Any) -> Any:
    """Convert an absolute temperature from kelvin to degrees Celsius.

    Book: §1.2 (T[°C] = T[K] − 273.15).

    Parameters
    ----------
    T_K : float or array_like
        Absolute temperature [K].

    Returns
    -------
    T_C : float or ndarray
        Temperature on the Celsius scale [°C].

    Notes
    -----
    Assumptions: a temperature value, not a difference (see :func:`celsius_to_kelvin`).

    Validation: V1 round trip with :func:`celsius_to_kelvin` (1e-12). Label: analytic.
    """
    import numpy as np

    out = np.asarray(T_K, dtype=float) - KELVIN_OFFSET  # §1.2: T[°C] = T[K] − 273.15
    return float(out) if out.ndim == 0 else out


def nondimensional(value: Any, name: str = "") -> float:
    """Assert a quantity is dimensionless (a Reynolds number, a similarity variable) and return it as a float."""
    if hasattr(value, "dimensionality") and not value.dimensionless:
        raise AssertionError(f"{name or 'value'} is not dimensionless: {value.dimensionality}")
    return float(value.magnitude) if hasattr(value, "magnitude") else float(value)
