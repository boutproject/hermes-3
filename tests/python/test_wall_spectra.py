"""Absolute fixed-band RMS and physical wall-crossing conventions."""

import numpy as np
import pytest
from energy_analysis.dataset import AuditError
from energy_analysis.spectra import analyse_native_bands, band_rms, e_fold_width
from numpy.testing import assert_allclose


@pytest.mark.parametrize("nz", (128, 256, 512))
def test_fixed_bands_preserve_physical_amplitudes(nz):
    z = 2 * np.pi * np.arange(nz) / nz
    field = 2 * np.sin(20 * z) + 3 * np.cos(40 * z) + 7 * np.sin(64 * z)
    assert_allclose(band_rms(field, 16, 31), np.sqrt(2), rtol=2e-14)
    assert_allclose(band_rms(field, 32, 63), 3 / np.sqrt(2), rtol=2e-14)


def test_unresolved_band_is_rejected():
    with pytest.raises(AuditError, match="not resolved"):
        band_rms(np.ones(32), 32, 63)


def test_wall_width_records_physical_distance_and_dominance():
    dx, length = 0.1, 0.5
    distance = (np.arange(100) + 0.5) * dx
    profile = np.exp(-distance / length)
    width = e_fold_width(profile, dx, 0.01)
    assert width["wall_dominated"]
    assert_allclose(width["mm"], length + dx / 2)
    assert_allclose(width["cells"], length / dx + 0.5)
    assert not e_fold_width(profile, dx, 2)["wall_dominated"]
    assert e_fold_width(np.ones(100), dx, 0.01)["mm"] is None


@pytest.fixture
def native_bands(tmp_path):
    from netCDF4 import Dataset

    z = 2 * np.pi * np.arange(128) / 128
    for rank in range(2):
        with Dataset(tmp_path / f"BOUT.dmp.{rank}.nc", "w") as dump:
            for name, size in (("t", 1), ("x", 12), ("y", 1), ("z", 128)):
                dump.createDimension(name, size)
            dump.createVariable("t_array", "f8", ("t",))[:] = [2]
            for name, value in (
                ("evolve_vorticity_ion_energy_version", 1),
                ("NXPE", 2),
                ("NYPE", 1),
                ("MXG", 2),
            ):
                dump.createVariable(name, "i4")[...] = value
            for name, value in (
                ("Vort", 2 * np.sin(20 * z) + 3 * np.cos(40 * z)),
                ("Nh+", 1),
                ("Tion", 1),
            ):
                dump.createVariable(name, "f8", ("t", "x", "y", "z"))[:] = value
    return tmp_path


def test_native_band_reader_accepts_current_version_on_all_ranks(native_bands):
    report = analyse_native_bands(native_bands)
    assert report["n"] == 16
    assert report["endpoint"] == 2
    for band, amplitude in (("16-31", np.sqrt(2)), ("32-63", 3 / np.sqrt(2))):
        assert_allclose(report["vorticity"][band]["profile"], amplitude, rtol=2e-14)


@pytest.mark.parametrize("rank", (0, 1))
@pytest.mark.parametrize("version", (None, 2, 3, 4, 5, 6, 7, 8, 9, 99))
def test_native_band_reader_rejects_unsupported_version_on_each_rank(
    native_bands, rank, version
):
    from netCDF4 import Dataset

    with Dataset(native_bands / f"BOUT.dmp.{rank}.nc", "a") as dump:
        if version is None:
            dump.renameVariable("evolve_vorticity_ion_energy_version", "unused_tag")
        else:
            dump["evolve_vorticity_ion_energy_version"][...] = version
    with pytest.raises(
        AuditError, match="Unsupported evolve_vorticity_ion_energy version"
    ):
        analyse_native_bands(native_bands)
