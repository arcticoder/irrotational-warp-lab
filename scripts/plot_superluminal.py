#!/usr/bin/env python3
"""
Visualize superluminal velocity sweep results.

Creates plots showing:
- E+, E-, E_net vs v
- Tail imbalance percentage vs v
- Ratio E+/|E-| at R2 vs v

Usage:
    python scripts/plot_superluminal.py results/experiments/superluminal/sweep_v.json \
        --out results/experiments/superluminal/sweep_v_plot.png
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def plot_superluminal_sweep(data: dict, out_path: str):
    """Create multi-panel visualization of superluminal sweep."""
    results = data["sweep"]["sweep_results"]
    
    v_vals = [r["v"] for r in results]
    e_pos = np.array([r["e_pos_inf"] for r in results])
    e_neg = np.array([abs(r["e_neg_inf"]) for r in results])
    e_net = np.array([r["e_net_inf"] for r in results])
    tail_frac = np.array([abs(r["tail_net_frac"]) * 100 for r in results])  # to percent
    ratio_r2 = np.array([r["ratio_r2"] for r in results])
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(
        f"Superluminal Velocity Sweep (ρ={data['sweep']['params']['rho']}, "
        f"σ={data['sweep']['params']['sigma']}, {data['sweep']['params']['mode']} mode)",
        fontsize=14,
        fontweight="bold",
    )
    
    # Panel 1: Energy components
    ax = axes[0, 0]
    ax.plot(v_vals, e_pos, "b-o", label="$E^+$ (positive)", linewidth=2, markersize=4)
    ax.plot(v_vals, e_neg, "r-s", label="$|E^-|$ (negative)", linewidth=2, markersize=4)
    ax.plot(v_vals, e_net, "g-^", label="$E_{net}$ (net)", linewidth=2, markersize=4)
    ax.set_xlabel("Velocity $v/c$", fontsize=11)
    ax.set_ylabel("Energy (tail-corrected)", fontsize=11)
    ax.legend(frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3)
    ax.set_title("Energy Components vs Velocity")
    
    # Panel 2: Tail imbalance percentage
    ax = axes[0, 1]
    ax.plot(v_vals, tail_frac, "m-o", linewidth=2, markersize=6)
    ax.set_xlabel("Velocity $v/c$", fontsize=11)
    ax.set_ylabel("$|E_{net}| / E_{abs}$ (%)", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_title("Tail Imbalance Percentage")
    
    # Panel 3: Ratio at R2
    ax = axes[1, 0]
    ax.plot(v_vals, ratio_r2, "c-s", linewidth=2, markersize=6)
    ax.axhline(y=1.0, color="k", linestyle="--", alpha=0.5, label="Perfect balance")
    ax.set_xlabel("Velocity $v/c$", fontsize=11)
    ax.set_ylabel("$E^+ / |E^-|$ at $R_2$", fontsize=11)
    ax.legend(frameon=True)
    ax.grid(True, alpha=0.3)
    ax.set_title("Finite-Window Ratio at $R_2$")
    
    # Panel 4: Scaling analysis (log-log)
    ax = axes[1, 1]
    # Normalize energies to v=1 values
    e_pos_norm = e_pos / e_pos[0]
    e_neg_norm = e_neg / e_neg[0]
    
    ax.loglog(v_vals, e_pos_norm, "b-o", label="$E^+$ (norm)", linewidth=2, markersize=4)
    ax.loglog(v_vals, e_neg_norm, "r-s", label="$|E^-|$ (norm)", linewidth=2, markersize=4)
    
    # Add reference lines for v^2 and v^4 scaling
    v_ref = np.array(v_vals)
    ax.loglog(v_ref, v_ref**2 / v_ref[0]**2, "k--", alpha=0.5, label="$v^2$ scaling")
    ax.loglog(v_ref, v_ref**4 / v_ref[0]**4, "k:", alpha=0.5, label="$v^4$ scaling")
    
    ax.set_xlabel("Velocity $v/c$", fontsize=11)
    ax.set_ylabel("Normalized energy (v=1 baseline)", fontsize=11)
    ax.legend(frameon=True, fancybox=True, shadow=True, fontsize=9)
    ax.grid(True, alpha=0.3, which="both")
    ax.set_title("Energy Scaling (log-log)")
    
    plt.tight_layout()
    
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {out_path}")
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    print(f"Velocity range: [{min(v_vals):.2f}, {max(v_vals):.2f}]")
    print(f"Tail imbalance range: [{min(tail_frac):.4f}%, {max(tail_frac):.4f}%]")
    print(f"Ratio(R2) range: [{min(ratio_r2):.4f}, {max(ratio_r2):.4f}]")
    
    # Check if energies scale roughly as v^2
    v_doubled_idx = np.argmin(np.abs(np.array(v_vals) - 2 * v_vals[0]))
    if v_doubled_idx > 0:
        e_ratio = e_pos[v_doubled_idx] / e_pos[0]
        print(f"\nEnergy scaling check (v doubled):")
        print(f"  E+(v={v_vals[v_doubled_idx]:.2f}) / E+(v={v_vals[0]:.2f}) = {e_ratio:.2f}")
        print(f"  Expected for v^2 scaling: {(v_vals[v_doubled_idx]/v_vals[0])**2:.2f}")
    print("=" * 60)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot superluminal sweep results")
    parser.add_argument("input", type=str, help="JSON file from sweep_superluminal.py")
    parser.add_argument("--out", type=str, required=True, help="Output plot path")
    
    args = parser.parse_args(argv)
    
    with open(args.input) as f:
        data = json.load(f)
    
    plot_superluminal_sweep(data, args.out)
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
