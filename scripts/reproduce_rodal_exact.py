from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np

from irrotational_warp.potential import phi_exact_rodal


def _resolve_backend(name: str):
    name = (name or "numpy").lower()
    if name == "numpy":
        return np, "numpy"
    if name == "cupy":
        try:
            import cupy as cp  # type: ignore

            return cp, "cupy"
        except Exception as exc:
            raise RuntimeError(
                "CuPy backend requested but CuPy could not be imported. "
                "Install cupy-cuda12x (or a matching build) and ensure WSL2 GPU passthrough works."
            ) from exc
    raise ValueError(f"Unknown backend: {name!r} (expected: numpy|cupy)")


def _git_info() -> dict:
    try:
        sha = (
            subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL)
            .decode("utf-8")
            .strip()
        )
        branch = (
            subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], stderr=subprocess.DEVNULL)
            .decode("utf-8")
            .strip()
        )
        dirty = subprocess.call(["git", "diff", "--quiet"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0
        return {"sha": sha, "branch": branch, "dirty": dirty}
    except Exception:
        return {"sha": None, "branch": None, "dirty": None}


def _two_point_tail_extrapolation(e_r1: float, e_r2: float, r1: float, r2: float) -> float:
    """Two-point 1/R extrapolation: E(∞) ≈ E(R2) + (R1/(R2-R1)) * (E(R2) - E(R1))."""
    factor = r1 / (r2 - r1)
    return float(e_r2 + factor * (e_r2 - e_r1))


def _split_pos_neg(xp, e_density) -> tuple[float, float, float]:
    pos = float(xp.sum(e_density * (e_density > 0.0)))
    neg = float(xp.sum(e_density * (e_density < 0.0)))
    net = float(pos + neg)
    return pos, neg, net


def _summary_metrics(e_pos: float, e_neg: float) -> dict:
    e_abs = e_pos + abs(e_neg)
    if e_abs == 0.0:
        return {"e_abs": 0.0, "neg_fraction": None, "net_over_abs": None, "ratio_e_pos_over_abs_e_neg": None}
    return {
        "e_abs": float(e_abs),
        "neg_fraction": float(abs(e_neg) / e_abs),
        "net_over_abs": float((e_pos + e_neg) / e_abs),
        "ratio_e_pos_over_abs_e_neg": float(e_pos / abs(e_neg)) if e_neg != 0.0 else None,
    }


def compute_energy_axisymmetric(
    *,
    rho: float,
    sigma: float,
    v: float,
    extent: float,
    nx: int,
    ny: int,
    backend: str,
    dtype: str,
) -> dict:
    """2.5D axisymmetric energy integral.

    This reconstructs a 3D volume integral assuming symmetry about the +x axis,
    using a meridional half-plane with cylindrical radius y and volume element:
      dV = 2π y dx dy.
    """
    xp, backend_name = _resolve_backend(backend)
    dtype_obj = xp.float32 if dtype == "float32" else xp.float64

    x = xp.linspace(-extent, extent, nx, dtype=dtype_obj)
    y = xp.linspace(0.0, extent, ny, dtype=dtype_obj)

    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])

    X, Y = xp.meshgrid(x, y, indexing="xy")
    Z = xp.zeros_like(X)

    t0 = time.perf_counter()
    phi = phi_exact_rodal(X, Y, Z, rho=rho, sigma=sigma, v=v, xp=xp)
    t_phi = time.perf_counter()

    grad_y, grad_x = xp.gradient(phi, dy, dx, edge_order=2)

    _, phi_xx = xp.gradient(grad_x, dy, dx, edge_order=2)
    phi_yy, _ = xp.gradient(grad_y, dy, dx, edge_order=2)
    _, phi_xy = xp.gradient(grad_y, dy, dx, edge_order=2)

    k_xx = -phi_xx
    k_yy = -phi_yy
    k_xy = -phi_xy

    try:
        k_zz = xp.zeros_like(grad_y)
        xp.divide(-grad_y, Y, out=k_zz, where=Y != 0.0)
        k_zz = xp.where(Y != 0.0, k_zz, k_yy)
    except TypeError:
        safe_y = xp.where(Y != 0.0, Y, 1.0)
        k_zz = xp.where(Y != 0.0, -grad_y / safe_y, k_yy)

    k_zz[0, :] = k_yy[0, :]

    k_trace = k_xx + k_yy + k_zz
    kij_kij = k_xx * k_xx + k_yy * k_yy + k_zz * k_zz + 2.0 * k_xy * k_xy
    rho_adm = (k_trace * k_trace - kij_kij) / (16.0 * xp.pi)

    dV = 2.0 * np.pi * Y * dx * dy
    e_density = rho_adm * dV
    e_pos, e_neg, e_net = _split_pos_neg(xp, e_density)

    r_grid = xp.sqrt(X * X + Y * Y)

    def e_at_r(r_lim: float) -> tuple[float, float]:
        mask = r_grid <= r_lim
        e_pos_sub = float(xp.sum(e_density * (e_density > 0.0) * mask))
        e_neg_sub = float(xp.sum(e_density * (e_density < 0.0) * mask))
        return e_pos_sub, e_neg_sub

    t1 = time.perf_counter()
    return {
        "mode": "axisym",
        "backend": backend_name,
        "dtype": dtype,
        "grid": {"nx": nx, "ny": ny, "extent": extent, "dx": dx, "dy": dy},
        "timing_s": {"phi": float(t_phi - t0), "total": float(t1 - t0)},
        "energies": {"e_pos": e_pos, "e_neg": e_neg, "e_net": e_net, **_summary_metrics(e_pos, e_neg)},
        "_e_at_r": e_at_r,
    }


def compute_energy_3d(
    *,
    rho: float,
    sigma: float,
    v: float,
    extent: float,
    n: int,
    backend: str,
    dtype: str,
) -> dict:
    """Full 3D Cartesian volume integration using dV = dx dy dz."""
    xp, backend_name = _resolve_backend(backend)
    dtype_obj = xp.float32 if dtype == "float32" else xp.float64

    x = xp.linspace(-extent, extent, n, dtype=dtype_obj)
    y = xp.linspace(-extent, extent, n, dtype=dtype_obj)
    z = xp.linspace(-extent, extent, n, dtype=dtype_obj)

    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    dz = float(z[1] - z[0])

    X, Y, Z = xp.meshgrid(x, y, z, indexing="ij")

    t0 = time.perf_counter()
    phi = phi_exact_rodal(X, Y, Z, rho=rho, sigma=sigma, v=v, xp=xp)
    t_phi = time.perf_counter()

    dphi_dx, dphi_dy, dphi_dz = xp.gradient(phi, dx, dy, dz, edge_order=2)

    g_dx = xp.gradient(dphi_dx, dx, dy, dz, edge_order=2)
    g_dy = xp.gradient(dphi_dy, dx, dy, dz, edge_order=2)
    g_dz = xp.gradient(dphi_dz, dx, dy, dz, edge_order=2)

    phi_xx = g_dx[0]
    phi_xy = g_dx[1]
    phi_xz = g_dx[2]
    phi_yy = g_dy[1]
    phi_yz = g_dy[2]
    phi_zz = g_dz[2]

    k_xx = -phi_xx
    k_yy = -phi_yy
    k_zz = -phi_zz
    k_xy = -phi_xy
    k_xz = -phi_xz
    k_yz = -phi_yz

    k_trace = k_xx + k_yy + k_zz
    kij_kij = (
        k_xx * k_xx
        + k_yy * k_yy
        + k_zz * k_zz
        + 2.0 * (k_xy * k_xy + k_xz * k_xz + k_yz * k_yz)
    )
    rho_adm = (k_trace * k_trace - kij_kij) / (16.0 * xp.pi)

    dV = dx * dy * dz
    e_density = rho_adm * dV
    e_pos, e_neg, e_net = _split_pos_neg(xp, e_density)

    r_grid = xp.sqrt(X * X + Y * Y + Z * Z)

    def e_at_r(r_lim: float) -> tuple[float, float]:
        mask = r_grid <= r_lim
        e_pos_sub = float(xp.sum(e_density * (e_density > 0.0) * mask))
        e_neg_sub = float(xp.sum(e_density * (e_density < 0.0) * mask))
        return e_pos_sub, e_neg_sub

    t1 = time.perf_counter()
    return {
        "mode": "3d",
        "backend": backend_name,
        "dtype": dtype,
        "grid": {"n": n, "extent": extent, "dx": dx, "dy": dy, "dz": dz},
        "timing_s": {"phi": float(t_phi - t0), "total": float(t1 - t0)},
        "energies": {"e_pos": e_pos, "e_neg": e_neg, "e_net": e_net, **_summary_metrics(e_pos, e_neg)},
        "_e_at_r": e_at_r,
    }


def _write_json(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Reproduce Rodal (2025) irrotational energetics.")
    parser.add_argument("--mode", choices=["axisym", "3d"], default="axisym")
    parser.add_argument("--backend", choices=["numpy", "cupy"], default="numpy")
    parser.add_argument("--dtype", choices=["float64", "float32"], default="float64")
    parser.add_argument("--rho", type=float, default=5.0)
    parser.add_argument("--sigma", type=float, default=4.0)
    parser.add_argument("--v", type=float, default=1.0)
    parser.add_argument("--extent-mult", type=float, default=12.0)
    parser.add_argument("--extent-pad", type=float, default=1.1)
    parser.add_argument("--nx", type=int, default=2400, help="axisym: x grid points")
    parser.add_argument("--ny", type=int, default=1200, help="axisym: y grid points")
    parser.add_argument("--n", type=int, default=100, help="3d: points per axis")
    parser.add_argument("--tail-r1-mult", type=float, default=8.0)
    parser.add_argument("--tail-r2-mult", type=float, default=12.0)
    parser.add_argument("--out", type=str, default=None, help="Write JSON summary to this path")

    args = parser.parse_args(argv)

    extent = args.extent_mult * args.rho * args.extent_pad
    r1 = args.tail_r1_mult * args.rho
    r2 = args.tail_r2_mult * args.rho
    if r2 <= r1:
        raise ValueError("tail-r2-mult must be > tail-r1-mult")

    started = time.perf_counter()
    if args.mode == "axisym":
        run = compute_energy_axisymmetric(
            rho=args.rho,
            sigma=args.sigma,
            v=args.v,
            extent=extent,
            nx=args.nx,
            ny=args.ny,
            backend=args.backend,
            dtype=args.dtype,
        )
    else:
        run = compute_energy_3d(
            rho=args.rho,
            sigma=args.sigma,
            v=args.v,
            extent=extent,
            n=args.n,
            backend=args.backend,
            dtype=args.dtype,
        )

    e_at_r = run.pop("_e_at_r")
    e_pos_r1, e_neg_r1 = e_at_r(r1)
    e_pos_r2, e_neg_r2 = e_at_r(r2)
    e_pos_inf = _two_point_tail_extrapolation(e_pos_r1, e_pos_r2, r1, r2)
    e_neg_inf = _two_point_tail_extrapolation(e_neg_r1, e_neg_r2, r1, r2)

    tail = {
        "r1": r1,
        "r2": r2,
        "at_r1": {"e_pos": e_pos_r1, "e_neg": e_neg_r1, "e_net": float(e_pos_r1 + e_neg_r1), **_summary_metrics(e_pos_r1, e_neg_r1)},
        "at_r2": {"e_pos": e_pos_r2, "e_neg": e_neg_r2, "e_net": float(e_pos_r2 + e_neg_r2), **_summary_metrics(e_pos_r2, e_neg_r2)},
        "at_inf": {
            "e_pos": e_pos_inf,
            "e_neg": e_neg_inf,
            "e_net": float(e_pos_inf + e_neg_inf),
            **_summary_metrics(e_pos_inf, e_neg_inf),
        },
    }

    payload = {
        "meta": {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "git": _git_info(),
            "python": sys.version,
        },
        "params": {"rho": args.rho, "sigma": args.sigma, "v": args.v},
        "run": run,
        "tail_2pt": tail,
        "wall_clock_s": float(time.perf_counter() - started),
    }

    print("-" * 60)
    print(f"Mode: {args.mode}")
    print(f"Params: rho={args.rho}, sigma={args.sigma}, v={args.v}")
    if args.mode == "axisym":
        print(f"Grid: nx={run['grid']['nx']}, ny={run['grid']['ny']}, extent={run['grid']['extent']}")
    else:
        print(f"Grid: n={run['grid']['n']}, extent={run['grid']['extent']}")
    r2_ratio = tail["at_r2"]["ratio_e_pos_over_abs_e_neg"]
    inf_abs = tail["at_inf"]["e_abs"]
    inf_net = tail["at_inf"]["e_net"]
    print(f"At R2={r2:.3f}: ratio E+/|E-| = {r2_ratio:.6g}")
    print(f"Tail net/abs (∞): {abs(inf_net)/inf_abs:.6%}")
    print("-" * 60)

    if args.out:
        _write_json(args.out, payload)
        print(f"Wrote JSON: {args.out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
