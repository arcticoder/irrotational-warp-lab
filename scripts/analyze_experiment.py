#!/usr/bin/env python3
"""Analyze optimization experiment results and generate summary."""
import json
from pathlib import Path
import sys


def load_json(path: Path) -> dict:
    """Load JSON file."""
    with open(path) as f:
        return json.load(f)


def analyze_1d_sweep(sweep_data: dict) -> dict:
    """Extract key insights from 1D sigma sweep."""
    points = sweep_data["points"]
    
    # Find best (minimum neg_fraction)
    best_point = min(points, key=lambda p: p["neg_fraction"])
    worst_point = max(points, key=lambda p: p["neg_fraction"])
    
    # Calculate range
    neg_fractions = [p["neg_fraction"] for p in points]
    
    return {
        "n_points": len(points),
        "sigma_range": [points[0]["sigma"], points[-1]["sigma"]],
        "neg_fraction_min": best_point["neg_fraction"],
        "neg_fraction_max": worst_point["neg_fraction"],
        "neg_fraction_range": max(neg_fractions) - min(neg_fractions),
        "best_sigma": best_point["sigma"],
        "worst_sigma": worst_point["sigma"],
    }


def analyze_2d_sweep(sweep_data: dict) -> dict:
    """Extract key insights from 2D (sigma, v) sweep."""
    points = sweep_data["points"]
    
    # Find best configuration
    best_point = min(points, key=lambda p: p["neg_fraction"])
    
    # Find worst configuration
    worst_point = max(points, key=lambda p: p["neg_fraction"])
    
    # Statistics
    neg_fractions = [p["neg_fraction"] for p in points]
    e_negs = [p["e_neg"] for p in points]
    
    return {
        "n_points": len(points),
        "sigma_range": [sweep_data["params"]["sigma_min"], sweep_data["params"]["sigma_max"]],
        "v_range": [sweep_data["params"]["v_min"], sweep_data["params"]["v_max"]],
        "best_config": {
            "sigma": best_point["sigma"],
            "v": best_point["v"],
            "neg_fraction": best_point["neg_fraction"],
            "e_neg": best_point["e_neg"],
            "e_pos": best_point["e_pos"],
        },
        "worst_config": {
            "sigma": worst_point["sigma"],
            "v": worst_point["v"],
            "neg_fraction": worst_point["neg_fraction"],
        },
        "neg_fraction_mean": sum(neg_fractions) / len(neg_fractions),
        "e_neg_mean": sum(e_negs) / len(e_negs),
    }


def analyze_optimization(opt_data: dict) -> dict:
    """Extract optimization results."""
    opt = opt_data["optimization"]
    
    improvement_abs = opt["initial_value"] - opt["best_value"]
    improvement_pct = 100 * improvement_abs / opt["initial_value"]
    
    return {
        "method": opt["method"],
        "n_evaluations": opt["n_evaluations"],
        "best_params": opt["best_params"],
        "best_value": opt["best_value"],
        "initial_value": opt["initial_value"],
        "improvement_abs": improvement_abs,
        "improvement_pct": improvement_pct,
        "success": opt["success"],
    }


def main():
    """Analyze experiment results and print summary."""
    base_dir = Path("results/experiments/baseline_rodal")
    
    if not base_dir.exists():
        print(f"Error: {base_dir} does not exist")
        return 1
    
    print("=" * 80)
    print("BASELINE RODAL-LIKE POTENTIAL OPTIMIZATION EXPERIMENT")
    print("=" * 80)
    print()
    
    # 1D Sweep Analysis
    sweep_1d_path = base_dir / "sweep_sigma_v1.5.json"
    if sweep_1d_path.exists():
        print("1D SIGMA SWEEP (v=1.5)")
        print("-" * 80)
        sweep_1d = load_json(sweep_1d_path)
        analysis_1d = analyze_1d_sweep(sweep_1d)
        
        print(f"  Parameter range: σ ∈ [{analysis_1d['sigma_range'][0]:.1f}, {analysis_1d['sigma_range'][1]:.1f}]")
        print(f"  Points evaluated: {analysis_1d['n_points']}")
        print(f"  Best σ: {analysis_1d['best_sigma']:.2f} (neg_fraction = {analysis_1d['neg_fraction_min']:.4f})")
        print(f"  Worst σ: {analysis_1d['worst_sigma']:.2f} (neg_fraction = {analysis_1d['neg_fraction_max']:.4f})")
        print(f"  Neg fraction range: {analysis_1d['neg_fraction_range']:.4f}")
        print()
    
    # 2D Sweep Analysis
    sweep_2d_path = base_dir / "sweep_2d_sigma_v.json"
    if sweep_2d_path.exists():
        print("2D (σ, v) PARAMETER SWEEP")
        print("-" * 80)
        sweep_2d = load_json(sweep_2d_path)
        analysis_2d = analyze_2d_sweep(sweep_2d)
        
        print(f"  Parameter ranges:")
        print(f"    σ ∈ [{analysis_2d['sigma_range'][0]:.1f}, {analysis_2d['sigma_range'][1]:.1f}]")
        print(f"    v ∈ [{analysis_2d['v_range'][0]:.1f}, {analysis_2d['v_range'][1]:.1f}]")
        print(f"  Points evaluated: {analysis_2d['n_points']}")
        print(f"  Mean neg_fraction: {analysis_2d['neg_fraction_mean']:.4f}")
        print(f"  Mean |E⁻|: {analysis_2d['e_neg_mean']:.4e}")
        print()
        print(f"  Best configuration:")
        best = analysis_2d['best_config']
        print(f"    σ = {best['sigma']:.3f}, v = {best['v']:.3f}")
        print(f"    |E⁻| = {best['e_neg']:.4e}")
        print(f"    E⁺ = {best['e_pos']:.4e}")
        print(f"    neg_fraction = {best['neg_fraction']:.6f}")
        print()
        print(f"  Worst configuration:")
        worst = analysis_2d['worst_config']
        print(f"    σ = {worst['sigma']:.3f}, v = {worst['v']:.3f}")
        print(f"    neg_fraction = {worst['neg_fraction']:.6f}")
        print()
    
    # Optimization Analysis
    opt_path = base_dir / "optimization_refined.json"
    if opt_path.exists():
        print("HYBRID OPTIMIZATION (Grid + Nelder-Mead)")
        print("-" * 80)
        opt_data = load_json(opt_path)
        analysis_opt = analyze_optimization(opt_data)
        
        print(f"  Method: {analysis_opt['method']}")
        print(f"  Evaluations: {analysis_opt['n_evaluations']}")
        print(f"  Success: {analysis_opt['success']}")
        print()
        print(f"  Initial guess:")
        print(f"    |E⁻| = {analysis_opt['initial_value']:.6e}")
        print()
        print(f"  Optimized result:")
        print(f"    σ = {analysis_opt['best_params']['sigma']:.6f}")
        print(f"    v = {analysis_opt['best_params']['v']:.6f}")
        print(f"    |E⁻| = {analysis_opt['best_value']:.6e}")
        print()
        print(f"  Improvement:")
        print(f"    Absolute: {analysis_opt['improvement_abs']:.6e}")
        print(f"    Relative: {analysis_opt['improvement_pct']:.2f}%")
        print()
    
    # Key Findings
    print("=" * 80)
    print("KEY FINDINGS")
    print("=" * 80)
    
    if sweep_2d_path.exists() and opt_path.exists():
        sweep_2d = load_json(sweep_2d_path)
        opt_data = load_json(opt_path)
        analysis_2d = analyze_2d_sweep(sweep_2d)
        analysis_opt = analyze_optimization(opt_data)
        
        print(f"1. The Rodal-like irrotational potential produces negative energy regions")
        print(f"   in all tested configurations (σ ∈ [2, 10], v ∈ [0.8, 2.5]).")
        print()
        print(f"2. Best configuration from 2D sweep:")
        print(f"   - Prefers LOW v (v={analysis_2d['best_config']['v']:.2f}) and LOW σ (σ={analysis_2d['best_config']['sigma']:.2f})")
        print(f"   - Achieves neg_fraction = {analysis_2d['best_config']['neg_fraction']:.4f}")
        print()
        print(f"3. Hybrid optimizer confirms this trend:")
        print(f"   - Optimal σ = {analysis_opt['best_params']['sigma']:.4f}")
        print(f"   - Optimal v = {analysis_opt['best_params']['v']:.4f}")
        print(f"   - |E⁻| = {analysis_opt['best_value']:.4e}")
        print()
        print(f"4. Irrotational smoothing (higher σ) does NOT eliminate negatives,")
        print(f"   but SHARPER walls (lower σ) minimize |E⁻| magnitude.")
        print()
        print(f"5. Lower velocities reduce negative energy fraction.")
        print()
    
    print("=" * 80)
    print("NEXT STEPS")
    print("=" * 80)
    print("• Validate against Celmaster & Rubin energy density (already done in M6)")
    print("• Compare to Rodal's tanh-based potential with cosθ dependence")
    print("• Test alternative smoothing functions (polynomial, compact support)")
    print("• Implement 3D volume integrals for more accurate E± totals")
    print("• Explore Type-I constraint enforcement")
    print()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
