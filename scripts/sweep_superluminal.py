#!/usr/bin/env python3
"""
Superluminal parameter sweep (v > 1) for Rodal exact potential.

Sweeps over velocity values v ∈ [1, v_max] and monitors:
- Energy integrals (E+, E-, E_net)
- Tail-corrected totals
- Numerical stability indicators

Usage:
    python scripts/sweep_superluminal.py --rho 5 --sigma 4 --v-max 3.0 --v-steps 20 \
        --mode axisym --nx 1200 --ny 600 --out results/superluminal_sweep.json
"""

import argparse
import json
import sys
import time
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from reproduce_rodal_exact import (
    compute_energy_axisymmetric,
    compute_energy_3d,
    _two_point_tail_extrapolation,
    _git_info,
)


def sweep_velocity(
    *,
    mode: str,
    rho: float,
    sigma: float,
    v_values: list[float],
    extent: float,
    tail_r1_mult: float,
    tail_r2_mult: float,
    backend: str,
    dtype: str,
    nx: int = None,
    ny: int = None,
    n: int = None,
) -> dict:
    """
    Sweep over velocity values and compute energies.
    
    Returns dict with sweep results and metadata.
    """
    r1 = tail_r1_mult * rho
    r2 = tail_r2_mult * rho
    
    results = []
    
    for i, v in enumerate(v_values):
        print(f"[{i+1}/{len(v_values)}] v = {v:.3f} ...", end=" ", flush=True)
        t0 = time.perf_counter()
        
        if mode == "axisym":
            run = compute_energy_axisymmetric(
                rho=rho,
                sigma=sigma,
                v=v,
                extent=extent,
                nx=nx,
                ny=ny,
                backend=backend,
                dtype=dtype,
            )
        else:
            run = compute_energy_3d(
                rho=rho,
                sigma=sigma,
                v=v,
                extent=extent,
                n=n,
                backend=backend,
                dtype=dtype,
            )
        
        e_at_r = run.pop("_e_at_r")
        e_pos_r1, e_neg_r1 = e_at_r(r1)
        e_pos_r2, e_neg_r2 = e_at_r(r2)
        e_pos_inf = _two_point_tail_extrapolation(e_pos_r1, e_pos_r2, r1, r2)
        e_neg_inf = _two_point_tail_extrapolation(e_neg_r1, e_neg_r2, r1, r2)
        
        e_abs_inf = abs(e_pos_inf) + abs(e_neg_inf)
        e_net_inf = e_pos_inf + e_neg_inf
        tail_frac = abs(e_net_inf) / e_abs_inf if e_abs_inf > 0 else 0.0
        ratio_r2 = e_pos_r2 / abs(e_neg_r2) if e_neg_r2 != 0 else None
        
        elapsed = time.perf_counter() - t0
        
        results.append({
            "v": v,
            "e_pos_inf": e_pos_inf,
            "e_neg_inf": e_neg_inf,
            "e_net_inf": e_net_inf,
            "e_abs_inf": e_abs_inf,
            "tail_net_frac": tail_frac,
            "ratio_r2": ratio_r2,
            "timing_s": elapsed,
        })
        
        print(f"E_net/E_abs(∞)={tail_frac:.5%}, ratio(R2)={ratio_r2:.4f}, t={elapsed:.2f}s")
    
    return {
        "sweep_results": results,
        "params": {
            "mode": mode,
            "rho": rho,
            "sigma": sigma,
            "extent": extent,
            "r1": r1,
            "r2": r2,
            "backend": backend,
            "dtype": dtype,
        },
        "grid": run["grid"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Superluminal velocity sweep")
    parser.add_argument("--mode", choices=["axisym", "3d"], default="axisym")
    parser.add_argument("--backend", choices=["numpy", "cupy"], default="numpy")
    parser.add_argument("--dtype", choices=["float64", "float32"], default="float64")
    parser.add_argument("--rho", type=float, default=5.0)
    parser.add_argument("--sigma", type=float, default=4.0)
    parser.add_argument("--v-min", type=float, default=1.0)
    parser.add_argument("--v-max", type=float, default=3.0)
    parser.add_argument("--v-steps", type=int, default=20)
    parser.add_argument("--extent-mult", type=float, default=12.0)
    parser.add_argument("--extent-pad", type=float, default=1.1)
    parser.add_argument("--nx", type=int, default=1200)
    parser.add_argument("--ny", type=int, default=600)
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--tail-r1-mult", type=float, default=8.0)
    parser.add_argument("--tail-r2-mult", type=float, default=12.0)
    parser.add_argument("--out", type=str, required=True)
    
    args = parser.parse_args(argv)
    
    import numpy as np
    v_values = np.linspace(args.v_min, args.v_max, args.v_steps).tolist()
    extent = args.extent_mult * args.rho * args.extent_pad
    
    print("=" * 60)
    print("SUPERLUMINAL VELOCITY SWEEP")
    print("=" * 60)
    print(f"Mode: {args.mode}")
    print(f"Backend: {args.backend}")
    print(f"Velocity range: v ∈ [{args.v_min}, {args.v_max}] ({args.v_steps} steps)")
    print(f"Parameters: rho={args.rho}, sigma={args.sigma}")
    if args.mode == "axisym":
        print(f"Grid: {args.nx}×{args.ny}")
    else:
        print(f"Grid: {args.n}³")
    print("=" * 60)
    
    started = time.perf_counter()
    
    sweep_data = sweep_velocity(
        mode=args.mode,
        rho=args.rho,
        sigma=args.sigma,
        v_values=v_values,
        extent=extent,
        tail_r1_mult=args.tail_r1_mult,
        tail_r2_mult=args.tail_r2_mult,
        backend=args.backend,
        dtype=args.dtype,
        nx=args.nx,
        ny=args.ny,
        n=args.n,
    )
    
    total_time = time.perf_counter() - started
    
    payload = {
        "meta": {
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "git": _git_info(),
            "python": sys.version,
        },
        "sweep": sweep_data,
        "total_time_s": total_time,
    }
    
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    
    print("=" * 60)
    print(f"Sweep completed in {total_time:.1f}s")
    print(f"Results written to: {args.out}")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
