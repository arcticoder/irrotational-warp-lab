import numpy as np

from irrotational_warp.adm import compute_slice_z0


def test_v_zero_gives_zero_density():
    res = compute_slice_z0(rho=10.0, sigma=5.0, v=0.0, extent=10.0, n=81)
    assert np.allclose(res.rho_adm, 0.0, atol=1e-10)


def test_y_symmetry_on_z0_slice():
    res = compute_slice_z0(rho=10.0, sigma=5.0, v=1.2, extent=10.0, n=101)
    flipped = np.flip(res.rho_adm, axis=0)
    assert np.allclose(res.rho_adm, flipped, atol=1e-6)
