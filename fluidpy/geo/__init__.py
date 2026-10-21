"""fluidpy.geo — geophysical fluid dynamics for Vallis, *Atmospheric and Oceanic Fluid Dynamics* (2nd ed., 2017).

Reusable GFD physics used by the AOFD book (``shammun/aofd``): spherical geometry, rotating frames, balances, seawater,
vertical coordinates, vorticity and PV, QG, Ekman layers, waves, instabilities and reduced models. Where the physics
already exists in ``fluidpy.core`` (the Kundu modules ``gfd``, ``rotating``, ``shallow_water``, ``stratification``,
``vertical_modes``, ``thermo``, ``waves``) this package imports it rather than re-implementing it, and it never changes
a ``fluidpy.core`` signature (the Kundu notebooks import them).

Docstring contract (as fluidpy.core): Vallis section and equation, symbols with SI units, assumptions, validation label,
plus "Kundu: K§/eq" when the idea is shared with Kundu.
"""
