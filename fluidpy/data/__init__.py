"""fluidpy.data — open-data access for the AOFD book: registry, fetchers, reductions and labelled pseudo-data.

Every download goes through one door, :func:`fluidpy.data.fetch` (added with the first fetcher), which checks the
disk guard before downloading. Raw files live under the data root (``AOFD_DATA_ROOT``; on Colab ``/content/aofd-data``),
never inside a repository. Remote OPeNDAP reads on Windows use ``engine="pydap"`` (netCDF4's bundled curl cannot find
TLS certificates there).
"""
