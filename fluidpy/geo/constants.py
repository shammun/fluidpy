"""Planetary and gas constants for ``fluidpy.geo`` — the one home of every constant the Vallis modules use.

Vallis, *Atmospheric and Oceanic Fluid Dynamics* (2nd ed., 2017): §2.1 (Ω), §2.2–2.3 (a, f, β), §1.6 and §2.6 (R, c_p,
g). Where ``fluidpy.core`` already defines the value it is **re-exported**, never re-typed, so the Kundu and Vallis
modules can never disagree.

| name | value | unit | source |
|---|---|---|---|
| ``OMEGA`` | 7.292115e-5 | rad s⁻¹ | sidereal rotation rate, WGS-84 (``core.rotating.OMEGA_EARTH``) |
| ``EARTH_RADIUS`` (a) | 6.371e6 | m | mean radius of the spherical Earth (``core.gfd.EARTH_RADIUS_MEAN``) |
| ``SIDEREAL_DAY`` | 2π/Ω ≈ 86 164 | s | one turn relative to the stars (``core.gfd.SIDEREAL_DAY``) |
| ``G0`` | 9.80665 | m s⁻² | standard gravity, exact (``core.thermo.G0``); also the g₀ that defines geopotential height Z = Φ/g₀ |
| ``R_D`` | ≈ 287.05 | J kg⁻¹ K⁻¹ | gas constant of dry air R_u/M_air (``core.thermo.R_AIR``) |
| ``C_P`` | ≈ 1004.7 | J kg⁻¹ K⁻¹ | c_p of dry air, γR/(γ − 1) with γ = 1.4 (``core.thermo.CP_AIR``) |
| ``KAPPA`` | R_D/C_P = 2/7 | dimensionless | exponent of potential temperature |
| ``P_REF`` | 1.0e5 | Pa | reference pressure for potential temperature (``core.thermo.P_REF``) |
| ``RHO0_OCEAN`` | 1025.0 | kg m⁻³ | Boussinesq reference density of seawater (Vallis's typical ρ₀) |

Validation: V1 tests pin each value to its source (identity with the ``core`` constant; KAPPA = 2/7 exactly).
Note: ``SIDEREAL_DAY`` = 2π/Ω = 86 164.1006 s differs from the IERS sidereal day (≈ 86 164.09 s) by ~10 ms because Ω
is the rounded WGS-84 value (cosmetic; the table's "≈ 86 164" holds).
"""
from __future__ import annotations

from ..core.gfd import EARTH_RADIUS_MEAN, SIDEREAL_DAY
from ..core.rotating import OMEGA_EARTH
from ..core.thermo import CP_AIR, G0, P_REF, R_AIR

__all__ = ["OMEGA", "EARTH_RADIUS", "SIDEREAL_DAY", "G0", "R_D", "C_P", "KAPPA", "P_REF", "RHO0_OCEAN"]

OMEGA: float = OMEGA_EARTH            #: Earth's sidereal rotation rate Ω [rad/s] (WGS-84)
EARTH_RADIUS: float = EARTH_RADIUS_MEAN  #: mean Earth radius a [m]
R_D: float = R_AIR                    #: dry-air gas constant R [J/(kg K)]
C_P: float = CP_AIR                   #: dry-air specific heat at constant pressure c_p [J/(kg K)]
KAPPA: float = R_D / C_P              #: κ = R/c_p [dimensionless] (= 2/7 for γ = 1.4)
RHO0_OCEAN: float = 1025.0            #: reference seawater density ρ₀ [kg/m³]
