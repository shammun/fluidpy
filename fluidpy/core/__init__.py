"""Reusable pieces shared by two or more chapters.

Machinery (present from the start): ``project`` (paths, config, URLs), ``style`` (figure style, notebook setup),
``embed`` (show_viz), ``anim`` (animations), ``interact`` (plotly sliders, widgets), ``units`` (dimensional checks),
``refdata`` (reference data).

Physics primitives are added the moment a second chapter needs something (grids, operators, ode, potential, waves,
similarity, thermo, …). Anything used by one chapter only stays in that chapter's module.
"""
