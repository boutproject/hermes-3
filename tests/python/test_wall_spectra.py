"""Absolute fixed-band RMS and physical wall-crossing conventions."""

import numpy as np
import pytest
from energy_analysis.dataset import AuditError
from energy_analysis.spectra import band_rms, e_fold_width
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
