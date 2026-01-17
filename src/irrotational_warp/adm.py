from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .fd import central_diff_2d, central_diff_vec_2d
from .potential import phi_dipole_cartesian


@dataclass(frozen=True)
class SliceResult:
    x: np.ndarray  # 1D
    y: np.ndarray  # 1D
    rho_adm: np.ndarray  # 2D [y,x]
    phi: np.ndarray  # 2D [y,x]
    beta_x: np.ndarray  # 2D [y,x]
    beta_y: np.ndarray  # 2D [y,x]
    dx: float
    dy: float


def compute_slice_z0(
    *,
    rho: float,
    sigma: float,
    v: float,
    extent: float,
    n: int,
) -> SliceResult:
    """Compute z=0 slice of Φ, β, and ADM energy density ρ_ADM.

    Uses flat-slice, unit-lapse ADM Hamiltonian constraint:
      ρ_ADM = (K^2 - K_ij K^ij) / (16π)
    with K_ij = 1/2(∂i βj + ∂j βi).

    This is a fast diagnostic; it is not the Rodal/McMonigal invariant eigenvalue diagnostic.
    """
    x = np.linspace(-extent, extent, n)
    y = np.linspace(-extent, extent, n)
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    phi = phi_dipole_cartesian(X, Y, Z, rho=rho, sigma=sigma, v=v)
    dphi_dx, dphi_dy = central_diff_2d(phi, dx=dx, dy=dy)

    beta_x = dphi_dx
    beta_y = dphi_dy

    dbx_dx, dbx_dy, dby_dx, dby_dy = central_diff_vec_2d(beta_x, beta_y, dx=dx, dy=dy)

    # In 2D slice we approximate K_ij for i,j in {x,y} only.
    # This is an approximation to the 3D quantity; M1 will add a full 3D implementation.
    k_xx = dbx_dx
    k_yy = dby_dy
    k_xy = 0.5 * (dbx_dy + dby_dx)

    k_trace = k_xx + k_yy
    kij_kij = k_xx * k_xx + k_yy * k_yy + 2.0 * k_xy * k_xy

    rho_adm = (k_trace * k_trace - kij_kij) / (16.0 * math.pi)

    return SliceResult(
        x=x,
        y=y,
        rho_adm=rho_adm,
        phi=phi,
        beta_x=beta_x,
        beta_y=beta_y,
        dx=dx,
        dy=dy,
    )


def integrate_signed(field: np.ndarray, dx: float, dy: float) -> tuple[float, float, float]:
    """Integrate positive/negative parts of a 2D field over area."""
    dA = dx * dy
    pos = float(np.sum(field[field > 0.0]) * dA)
    neg = float(np.sum(-field[field < 0.0]) * dA)
    net = pos - neg
    return pos, neg, net
