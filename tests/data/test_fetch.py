"""Smoke tests for fluidpy.data (no network). Derived-product checks skip when the AOFD repo is not beside fluidpy."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

import fluidpy.data as fd
from fluidpy.data import fronts, registry
from fluidpy.data.reduce import provenance, write_nc

AOFD = fd.aofd_repo()
DERIVED = AOFD / "data" / "derived" if AOFD else None


def _missing_drive() -> str:
    for d in "QRSTUVWXYZ":
        if not Path(f"{d}:/").exists():
            return f"{d}:/aofd-data"
    pytest.skip("no free drive letter to simulate a missing data drive")


def test_data_root_env(monkeypatch, tmp_path):
    monkeypatch.setenv("AOFD_DATA_ROOT", str(tmp_path))
    assert fd.data_root() == tmp_path


def test_data_root_colab(monkeypatch):
    monkeypatch.delenv("AOFD_DATA_ROOT", raising=False)
    monkeypatch.setenv("COLAB_RELEASE_TAG", "x")
    assert fd.data_root().as_posix() == "/content/aofd-data"


@pytest.mark.skipif(os.name != "nt", reason="drive letters are Windows-only")
@pytest.mark.parametrize("use_builtin", [True, False])
def test_missing_data_drive_never_falls_back(monkeypatch, use_builtin):
    monkeypatch.setenv("AOFD_DATA_ROOT", _missing_drive())
    if use_builtin:
        monkeypatch.setattr(fd, "_GUARD", False)
    elif fd._guard_module() is None:
        pytest.skip("AOFD tools/disk_guard.py not available")
    monkeypatch.setattr(fd, "_head_size", lambda url: 10)
    with pytest.raises(fd.DiskGuardError, match="not found"):
        fd.fetch("https://example.invalid/file.nc", quiet=True)


def test_builtin_guard_refuses_huge(monkeypatch, tmp_path):
    monkeypatch.setenv("AOFD_DATA_ROOT", str(tmp_path))
    monkeypatch.setattr(fd, "_GUARD", False)
    with pytest.raises(fd.DiskGuardError):
        fd._check(1e15, "absurd")


def test_unknown_size_is_refused(monkeypatch, tmp_path):
    monkeypatch.setenv("AOFD_DATA_ROOT", str(tmp_path))
    monkeypatch.setattr(fd, "_head_size", lambda url: None)
    with pytest.raises(fd.DiskGuardError, match="estimate"):
        fd.fetch("https://example.invalid/file.nc", quiet=True)


def test_open_remote_uses_pydap_for_opendap(monkeypatch):
    import xarray as xr

    seen = {}
    monkeypatch.setattr(xr, "open_dataset", lambda url, **kw: seen.update(url=url, **kw) or "ok")
    fd.open_remote("https://psl.noaa.gov/thredds/dodsC/Datasets/x.nc")
    if os.name == "nt":
        assert seen["engine"] == "pydap"
    seen.clear()
    fd.open_remote("C:/local.nc")
    assert "engine" not in seen


def test_registry_entries_have_required_fields():
    reg = registry.load()
    if not reg:
        pytest.skip("no AOFD registry beside fluidpy")
    for key, ent in reg.items():
        for f in ("url", "licence", "licence_url", "citation", "redistributable", "attribution", "fetch_date"):
            assert f in ent, (key, f)


def test_orsi_parser(tmp_path):
    p = tmp_path / "pf.txt"
    p.write_text("% seg 1\n 10.0 -50.0\n 20.0 -49.5\n%\n 350.0 -55.0\n 355.0 -56.0\n")
    lon, lat = fronts.read_orsi_txt(p)
    assert np.isnan(lon[2]) and lon.size == 5
    np.testing.assert_allclose(lon[[0, 1, 3, 4]], [10, 20, -10, -5])
    np.testing.assert_allclose(lat[[0, 1, 3, 4]], [-50, -49.5, -55, -56])


def test_write_nc_float32_and_provenance(tmp_path):
    import xarray as xr

    ds = xr.Dataset({"z": (("lat",), np.arange(3.0))}, coords={"lat": [-1.0, 0.0, 1.0]})
    p = write_nc(ds, tmp_path / "x.nc", provenance("K", "u", "PD", "cit", "none"))
    back = xr.open_dataset(p)
    assert back["z"].dtype == np.float32 and back.attrs["dataset_key"] == "K"


# ------------------------------------------------------------------ derived products (skip when absent)
def _need(rel):
    if DERIVED is None or not (DERIVED / rel).exists():
        pytest.skip(f"derived product {rel} not present")
    return DERIVED / rel


def test_z500_product():
    import xarray as xr

    ds = xr.open_dataset(_need("ch00/ncep_r1_z500_20140106.nc"))
    z = ds["z"]
    assert z.dtype == np.float32 and z.shape == (73, 144)
    assert np.all(np.diff(z.lat) > 0)                       # sorted S->N
    assert 4600 < float(z.min()) and float(z.max()) < 6000   # typical 500 hPa heights [m]
    # 2014-01-06: deep trough over the Great Lakes (~90 W) and ridge near the North American west coast (~120 W)
    r = z.sel(lat=50)
    assert 260 <= float(r.lon[int(r.argmin(dim="lon"))]) <= 290 and 225 <= float(r.lon[int(r.argmax(dim="lon"))]) <= 250
    assert os.path.getsize(_need("ch00/ncep_r1_z500_20140106.nc")) < 1_000_000


def test_fronts_product_latitudes_south_of_africa():
    f = np.load(_need("ch21/fig_21_23_fronts_orsi.npz"))
    near = {}
    for k in ("saf", "pf"):
        lo, la = f[f"{k}_lon"], f[f"{k}_lat"]
        m = np.isfinite(lo) & (np.abs(lo - 20) < 1.5)
        near[k] = float(np.mean(la[m]))
    assert -51.5 < near["pf"] < -48.5          # PF ~ 50 S at 20 E
    assert -49 < near["saf"] < -44             # SAF ~ 45-47 S at 20 E
    assert near["saf"] > near["pf"]
    assert "CC-BY" in json.loads(str(f["attrs"]))["licence"]


def test_gdp_product_drake_passage_eastward():
    import xarray as xr

    g = xr.open_dataset(_need("ch21/fig_21_23_gdp_mean_velocity.nc"))
    dp = g.sel(lon=slice(-66, -64), lat=slice(-60, -55))
    assert float(dp["U"].mean()) > 0.05                     # ACC flows east [m s-1]
    assert float(np.sqrt(g.U ** 2 + g.V ** 2).max()) < 3.0


def test_bathy_mask_classes():
    import xarray as xr

    m = xr.open_dataset(_need("ch21/fig_21_23_bathy3000_mask.nc"))["depth_class"]
    assert int(m.sel(lon=-62, lat=-58.5, method="nearest")) == 0   # Drake Passage deep
    assert int(m.sel(lon=70, lat=-49, method="nearest")) == 1      # Kerguelen Plateau shallow
    assert int(m.sel(lon=0, lat=-89, method="nearest")) == 2       # Antarctica land
