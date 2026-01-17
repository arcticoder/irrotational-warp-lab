"""
Validation module: Lentz-style warp potentials.

Implements potentials from Celmaster & Rubin (2024) "Violations of the Weak 
Energy Condition for Lentz Warp Drives" for cross-validation.

References:
- Celmaster & Rubin Eq. (18): phi_L (Lentz's original, flawed)
- Celmaster & Rubin Eq. (phiMod): phi_rh (corrected rhomboidal potential)
- Celmaster & Rubin Eqs. (srcRhom)-(chgPflRhom): rhomboidal source
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass


@dataclass
class RhomboidalSource:
    """Rhomboidal source parameters from Celmaster & Rubin."""
    # Rhomboid geometry
    L_x: float = 1.0      # semi-length along x
    L_z: float = 0.125    # semi-length along z
    
    # Centroids (beta_i, xi_i)
    beta: np.ndarray = None  # [1, 2, 1, 0]
    xi: np.ndarray = None    # [-1, 0, 1, 2]
    alpha: np.ndarray = None # [25, -25, -25, 50]
    
    # Charge profile
    Q_max: float = 1.0
    Q_min: float = 0.25
    
    # Warp velocity
    v_h: float = 1.0
    
    def __post_init__(self):
        if self.beta is None:
            self.beta = np.array([1.0, 2.0, 1.0, 0.0])
        if self.xi is None:
            self.xi = np.array([-1.0, 0.0, 1.0, 2.0])
        if self.alpha is None:
            self.alpha = np.array([25.0, -25.0, -25.0, 50.0])
    
    def eval_rho(self, s: float, z: float) -> float:
        """
        Evaluate rhomboidal source at (s, z).
        
        Celmaster & Rubin Eq. (srcRhom):
          rho(s,z) = W * sum_i alpha_i * Theta(p_i(s,z)) * g_i(s)
        
        where:
          p_i(s,z) = L_x - (L_x/L_z)|z - xi_i| - ||s| - beta_i|
          g_i(s) = Q_min + (Q_max - Q_min)(|s| - beta_i)^2
        
        NOTE: For beta_i != 0, there are symmetric pairs at +/- beta_i.
        The source is defined with s = |x| + |y|, so it's already symmetric.
        We evaluate based on |s| - beta_i, which covers both sides.
        """
        rho = 0.0
        for i in range(len(self.alpha)):
            # For beta_i != 0, this automatically handles +/- beta_i due to |s|
            # For beta_i == 0, there's only one rhomboid at the origin
            
            # Shape function p_i (using |s| makes it symmetric)
            p = self.L_x - (self.L_x / self.L_z) * abs(z - self.xi[i]) - abs(abs(s) - self.beta[i])
            
            if p > 0:  # Heaviside step
                # Charge profile g_i
                g = self.Q_min + (self.Q_max - self.Q_min) * (abs(s) - self.beta[i])**2
                rho += self.alpha[i] * g
        
        return rho


def phi_L_lentz(x: float, y: float, z: float, source: RhomboidalSource) -> float:
    """
    Lentz's original potential (flawed).
    
    Celmaster & Rubin Eq. (phiLentz):
      phi_L(x,y,z) = (1/4v_h) * int dx' dz' 
                     * Theta(z - z' - (1/v_h)|x - x'|)
                     * rho(|x'| + |y|, z')
    
    NOTE: This has asymmetry between x and y and wrong sign (per Celmaster & Rubin).
    """
    v_h = source.v_h
    
    # Numerical quadrature: integrate over source domain
    # Source is bounded by |s| < some value, |z| < some value
    # For rhomboids: max s ~ 2*L_x + max(beta), max z ~ max(xi) + L_z
    
    s_max = 2.0 * source.L_x + np.max(source.beta)
    z_max = np.max(source.xi) + source.L_z
    z_min = np.min(source.xi) - source.L_z
    
    # Grid for integration
    n_s = 200
    n_z = 200
    x_vals = np.linspace(-s_max, s_max, n_s)
    z_vals = np.linspace(z_min, z_max, n_z)
    dx_val = x_vals[1] - x_vals[0]
    dz_val = z_vals[1] - z_vals[0]
    
    integral = 0.0
    for x_prime in x_vals:
        for z_prime in z_vals:
            # Heaviside condition: z - z' - (1/v_h)|x - x'| > 0
            if z - z_prime - (1.0 / v_h) * abs(x - x_prime) > 0:
                # Source argument: |x'| + |y| (asymmetric!)
                s_arg = abs(x_prime) + abs(y)
                rho_val = source.eval_rho(s_arg, z_prime)
                integral += rho_val * dx_val * dz_val
    
    return (1.0 / (4.0 * v_h)) * integral


def phi_rh_corrected(x: float, y: float, z: float, source: RhomboidalSource) -> float:
    """
    Corrected rhomboidal potential (what Lentz intended).
    
    Celmaster & Rubin Eq. (phiMod):
      phi_rh(x,y,z) = -(1/4) * v_h * int ds' dz'
                      * Theta(z - z' - (1/v_h)|s - s'|)
                      * rho(s', z')
    
    where s = |x| + |y|.
    
    Key differences from phi_L:
    1. Sign is negative (not positive)
    2. Factor is v_h (not 1/v_h)
    3. Uses s = |x| + |y| (symmetric)
    """
    v_h = source.v_h
    s = abs(x) + abs(y)
    
    # Integration bounds
    s_max = 2.0 * source.L_x + np.max(source.beta)
    z_max = np.max(source.xi) + source.L_z
    z_min = np.min(source.xi) - source.L_z
    
    # Grid
    n_s = 200
    n_z = 200
    s_vals = np.linspace(-s_max, s_max, n_s)
    z_vals = np.linspace(z_min, z_max, n_z)
    ds_val = s_vals[1] - s_vals[0]
    dz_val = z_vals[1] - z_vals[0]
    
    integral = 0.0
    for s_prime in s_vals:
        for z_prime in z_vals:
            # Heaviside: z - z' - (1/v_h)|s - s'| > 0
            if z - z_prime - (1.0 / v_h) * abs(s - s_prime) > 0:
                rho_val = source.eval_rho(s_prime, z_prime)
                integral += rho_val * ds_val * dz_val
    
    return -(1.0 / 4.0) * v_h * integral


def compute_shift_vector(phi_func, x: float, y: float, z: float, 
                         source: RhomboidalSource, 
                         h: float = 1e-5) -> Tuple[float, float, float]:
    """
    Compute shift vector N = grad(phi) via finite differences.
    
    Returns:
        (N_x, N_y, N_z)
    """
    # Central differences
    N_x = (phi_func(x + h, y, z, source) - phi_func(x - h, y, z, source)) / (2 * h)
    N_y = (phi_func(x, y + h, z, source) - phi_func(x, y - h, z, source)) / (2 * h)
    N_z = (phi_func(x, y, z + h, source) - phi_func(x, y, z - h, source)) / (2 * h)
    
    return N_x, N_y, N_z


def compute_extrinsic_curvature(N_x: float, N_y: float, N_z: float,
                                dN_x_dx: float, dN_x_dy: float, dN_x_dz: float,
                                dN_y_dx: float, dN_y_dy: float, dN_y_dz: float,
                                dN_z_dx: float, dN_z_dy: float, dN_z_dz: float) -> np.ndarray:
    """
    Compute extrinsic curvature K_ij = (1/2)(∂_i N_j + ∂_j N_i).
    
    Returns:
        3x3 symmetric matrix K
    """
    K = np.zeros((3, 3))
    
    K[0, 0] = dN_x_dx
    K[0, 1] = K[1, 0] = 0.5 * (dN_x_dy + dN_y_dx)
    K[0, 2] = K[2, 0] = 0.5 * (dN_x_dz + dN_z_dx)
    
    K[1, 1] = dN_y_dy
    K[1, 2] = K[2, 1] = 0.5 * (dN_y_dz + dN_z_dy)
    
    K[2, 2] = dN_z_dz
    
    return K


def compute_eulerian_energy(K: np.ndarray) -> float:
    """
    Compute Eulerian energy density from extrinsic curvature.
    
    E = (1/(8π)) * (1/2) * (-K^i_j K^j_i + K^2)
    
    where K = tr(K).
    """
    K_trace = np.trace(K)
    K_squared_trace = np.trace(K @ K)
    
    E = (1.0 / (8.0 * np.pi)) * 0.5 * (-K_squared_trace + K_trace**2)
    
    return E


def compute_energy_density_at_point(phi_func, x: float, y: float, z: float,
                                    source: RhomboidalSource, h: float = 1e-5) -> float:
    """
    Compute energy density at a point by computing shift vector derivatives
    and extrinsic curvature.
    
    Args:
        phi_func: Potential function (phi_L_lentz or phi_rh_corrected)
        x, y, z: Point coordinates
        source: Rhomboidal source configuration
        h: Finite difference step size
        
    Returns:
        E: Energy density at (x, y, z)
    """
    # Compute shift vector components
    N_x, N_y, N_z = compute_shift_vector(phi_func, x, y, z, source, h)
    
    # Compute derivatives of shift vector components via finite differences
    # dN_x/dx, dN_x/dy, dN_x/dz
    N_xp, _, _ = compute_shift_vector(phi_func, x + h, y, z, source, h)
    N_xm, _, _ = compute_shift_vector(phi_func, x - h, y, z, source, h)
    dN_x_dx = (N_xp - N_xm) / (2 * h)
    
    N_xp, _, _ = compute_shift_vector(phi_func, x, y + h, z, source, h)
    N_xm, _, _ = compute_shift_vector(phi_func, x, y - h, z, source, h)
    dN_x_dy = (N_xp - N_xm) / (2 * h)
    
    N_xp, _, _ = compute_shift_vector(phi_func, x, y, z + h, source, h)
    N_xm, _, _ = compute_shift_vector(phi_func, x, y, z - h, source, h)
    dN_x_dz = (N_xp - N_xm) / (2 * h)
    
    # dN_y/dx, dN_y/dy, dN_y/dz
    _, N_yp, _ = compute_shift_vector(phi_func, x + h, y, z, source, h)
    _, N_ym, _ = compute_shift_vector(phi_func, x - h, y, z, source, h)
    dN_y_dx = (N_yp - N_ym) / (2 * h)
    
    _, N_yp, _ = compute_shift_vector(phi_func, x, y + h, z, source, h)
    _, N_ym, _ = compute_shift_vector(phi_func, x, y - h, z, source, h)
    dN_y_dy = (N_yp - N_ym) / (2 * h)
    
    _, N_yp, _ = compute_shift_vector(phi_func, x, y, z + h, source, h)
    _, N_ym, _ = compute_shift_vector(phi_func, x, y, z - h, source, h)
    dN_y_dz = (N_yp - N_ym) / (2 * h)
    
    # dN_z/dx, dN_z/dy, dN_z/dz
    _, _, N_zp = compute_shift_vector(phi_func, x + h, y, z, source, h)
    _, _, N_zm = compute_shift_vector(phi_func, x - h, y, z, source, h)
    dN_z_dx = (N_zp - N_zm) / (2 * h)
    
    _, _, N_zp = compute_shift_vector(phi_func, x, y + h, z, source, h)
    _, _, N_zm = compute_shift_vector(phi_func, x, y - h, z, source, h)
    dN_z_dy = (N_zp - N_zm) / (2 * h)
    
    _, _, N_zp = compute_shift_vector(phi_func, x, y, z + h, source, h)
    _, _, N_zm = compute_shift_vector(phi_func, x, y, z - h, source, h)
    dN_z_dz = (N_zp - N_zm) / (2 * h)
    
    # Compute extrinsic curvature
    K = compute_extrinsic_curvature(
        N_x, N_y, N_z,
        dN_x_dx, dN_x_dy, dN_x_dz,
        dN_y_dx, dN_y_dy, dN_y_dz,
        dN_z_dx, dN_z_dy, dN_z_dz
    )
    
    # Compute energy density
    E = compute_eulerian_energy(K)
    
    return E
