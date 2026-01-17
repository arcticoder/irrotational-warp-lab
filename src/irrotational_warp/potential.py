from __future__ import annotations

import numpy as np


def tanh_wall(r: np.ndarray, rho: float, sigma: float) -> np.ndarray:
    """Smooth top-hat-like wall function f(r/rho) in [0, 1].

    f(0) ~ 1 (inside), f(r>>rho) ~ 0 (outside)
    """
    xi = r / rho
    return 0.5 * (1.0 + np.tanh(sigma * (1.0 - xi)))


def phi_dipole_cartesian(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    rho: float,
    sigma: float,
    v: float,
    eps: float = 1e-12,
) -> np.ndarray:
    """Rodal-like axisymmetric dipole potential oriented along +x.

    Φ = v * rho * f(r/rho) * cosθ, with cosθ = x/r.

    Parameters are geometric; v is dimensionless (v/c).
    """
    r = np.sqrt(x * x + y * y + z * z)
    f = tanh_wall(r, rho=rho, sigma=sigma)
    costheta = x / (r + eps)
    return v * rho * f * costheta
