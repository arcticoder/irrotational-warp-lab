from __future__ import annotations

import numpy as np


def central_diff_2d(a: np.ndarray, dx: float, dy: float) -> tuple[np.ndarray, np.ndarray]:
    """Central differences for a 2D scalar field a[y, x].

    Returns (da/dx, da/dy) on the same grid using edge-order=1 fallbacks.
    """
    dadx = np.gradient(a, dx, axis=1, edge_order=1)
    dady = np.gradient(a, dy, axis=0, edge_order=1)
    return dadx, dady


def central_diff_vec_2d(
    bx: np.ndarray, by: np.ndarray, dx: float, dy: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Central differences for 2D vector field b=(bx,by) with arrays [y,x].

    Returns: dbx/dx, dbx/dy, dby/dx, dby/dy.
    """
    dbx_dx, dbx_dy = central_diff_2d(bx, dx=dx, dy=dy)
    dby_dx, dby_dy = central_diff_2d(by, dx=dx, dy=dy)
    return dbx_dx, dbx_dy, dby_dx, dby_dy
