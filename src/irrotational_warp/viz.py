from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .adm import compute_slice_z0, integrate_signed
from .io import write_summary_json


def plot_slice(
    *,
    rho: float,
    sigma: float,
    v: float,
    extent: float,
    n: int,
    out_path: Path,
    json_out_path: Path,
) -> None:
    res = compute_slice_z0(rho=rho, sigma=sigma, v=v, extent=extent, n=n)

    epos, eneg, enet = integrate_signed(res.rho_adm, dx=res.dx, dy=res.dy)

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    vlim = float(np.nanmax(np.abs(res.rho_adm)))
    im = ax.imshow(
        res.rho_adm,
        origin="lower",
        extent=[res.x[0], res.x[-1], res.y[0], res.y[-1]],
        cmap="coolwarm",
        vmin=-vlim,
        vmax=vlim,
        interpolation="nearest",
    )
    ax.set_title(f"ADM energy density (z=0 slice) | v={v}, rho={rho}, sigma={sigma}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, label=r"$\rho_{\rm ADM}$ (geometric units)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    write_summary_json(
        json_out_path,
        {
            "params": {"rho": rho, "sigma": sigma, "v": v, "extent": extent, "n": n},
            "integrals_2d": {"E_pos": epos, "E_neg": eneg, "E_net": enet},
            "rho_adm": res.rho_adm,
        },
    )
