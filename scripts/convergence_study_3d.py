#!/usr/bin/env python3
"""
3D Convergence Study for Rodal Exact Potential

Runs the reproduction script at multiple grid resolutions to demonstrate
convergence behavior of tail-corrected energy integrals.

Usage:
    python scripts/convergence_study_3d.py --backend numpy --n-values 40,60,80,100 \
        --out results/experiments/convergence/study_3d.json
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


def run_convergence_study(
    *,
    mode: str,
    backend: str,
    dtype: str,
    n_values: list[int],
    rho: float,
    sigma: float,
    v: float,
    out_dir: Path,
):
    """Run reproduction at multiple resolutions."""
    results = []
    
    for i, n in enumerate(n_values):
        print(f"\n{'='*60}")
        print(f"[{i+1}/{len(n_values)}] Running n = {n}")
        print('='*60)
        
        out_path = out_dir / f"repro_n{n}.json"
        
        cmd = [
            sys.executable,
            "scripts/reproduce_rodal_exact.py",
            "--mode", mode,
            "--backend", backend,
            "--dtype", dtype,
            "--n", str(n),
            "--rho", str(rho),
            "--sigma", str(sigma),
            "--v", str(v),
            "--out", str(out_path),
        ]
        
        t0 = time.perf_counter()
        result = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.perf_counter() - t0
        
        if result.returncode != 0:
            print(f"❌ Failed with code {result.returncode}")
            print("STDERR:", result.stderr)
            continue
        
        # Load the result JSON
        with open(out_path) as f:
            data = json.load(f)
        
        tail = data["tail_2pt"]["at_inf"]
        
        results.append({
            "n": n,
            "n_points": n**3,
            "e_pos_inf": tail["e_pos"],
            "e_neg_inf": tail["e_neg"],
            "e_net_inf": tail["e_net"],
            "e_abs_inf": tail["e_abs"],
            "neg_fraction": tail["neg_fraction"],
            "net_over_abs_pct": tail["net_over_abs"] * 100,
            "ratio_r2": data["tail_2pt"]["at_r2"]["ratio_e_pos_over_abs_e_neg"],
            "wall_clock_s": elapsed,
            "file": str(out_path),
        })
        
        print(f"✅ Complete: net/abs={tail['net_over_abs']*100:.5f}%, "
              f"ratio(R2)={data['tail_2pt']['at_r2']['ratio_e_pos_over_abs_e_neg']:.4f}, "
              f"time={elapsed:.1f}s")
    
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description="3D convergence study")
    parser.add_argument("--mode", default="3d", choices=["3d"])
    parser.add_argument("--backend", default="numpy", choices=["numpy", "cupy"])
    parser.add_argument("--dtype", default="float64", choices=["float64", "float32"])
    parser.add_argument("--n-values", type=str, default="40,60,80,100",
                       help="Comma-separated list of grid resolutions")
    parser.add_argument("--rho", type=float, default=5.0)
    parser.add_argument("--sigma", type=float, default=4.0)
    parser.add_argument("--v", type=float, default=1.0)
    parser.add_argument("--out", type=str, required=True)
    
    args = parser.parse_args(argv)
    
    n_values = [int(x.strip()) for x in args.n_values.split(",")]
    out_path = Path(args.out)
    out_dir = out_path.parent / "convergence_runs"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("3D CONVERGENCE STUDY")
    print("=" * 60)
    print(f"Backend: {args.backend}")
    print(f"Grid resolutions: {n_values}")
    print(f"Parameters: rho={args.rho}, sigma={args.sigma}, v={args.v}")
    print("=" * 60)
    
    started = time.perf_counter()
    
    results = run_convergence_study(
        mode=args.mode,
        backend=args.backend,
        dtype=args.dtype,
        n_values=n_values,
        rho=args.rho,
        sigma=args.sigma,
        v=args.v,
        out_dir=out_dir,
    )
    
    total_time = time.perf_counter() - started
    
    # Analyze convergence
    print("\n" + "=" * 60)
    print("CONVERGENCE ANALYSIS")
    print("=" * 60)
    print(f"{'n':>6} {'Points':>10} {'|E_net/E_abs|%':>15} {'Ratio(R2)':>12} {'Time(s)':>10}")
    print("-" * 60)
    
    for r in results:
        print(f"{r['n']:>6} {r['n_points']:>10,} {abs(r['net_over_abs_pct']):>15.6f} "
              f"{r['ratio_r2']:>12.5f} {r['wall_clock_s']:>10.1f}")
    
    # Compute convergence rate estimate
    if len(results) >= 2:
        print("\nConvergence rate estimates (successive differences):")
        for i in range(1, len(results)):
            r1, r2 = results[i-1], results[i]
            d_net = abs(r2['net_over_abs_pct'] - r1['net_over_abs_pct'])
            d_ratio = abs(r2['ratio_r2'] - r1['ratio_r2'])
            print(f"  n={r1['n']}→{r2['n']}: Δ(net/abs%)={d_net:.6f}, Δ(ratio)={d_ratio:.6f}")
    
    print("=" * 60)
    
    payload = {
        "meta": {
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "python": sys.version,
        },
        "params": {
            "mode": args.mode,
            "backend": args.backend,
            "dtype": args.dtype,
            "rho": args.rho,
            "sigma": args.sigma,
            "v": args.v,
        },
        "results": results,
        "total_time_s": total_time,
    }
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    
    print(f"\nResults written to: {out_path}")
    print(f"Total time: {total_time:.1f}s")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
