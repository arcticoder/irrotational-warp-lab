# Baseline Rodal-Like Potential Optimization Experiment

**Date:** January 16, 2026  
**Objective:** Establish baseline performance of Rodal-style irrotational warp potential and search parameter space for configurations minimizing negative energy.

## Methodology

### Potential Definition
Simple axisymmetric Rodal-like potential (oriented along +x):

$$\Phi(x,y,z) = v \cdot \rho \cdot f(r/\rho) \cdot \cos\theta$$

where:
- $f(\xi) = \frac{1}{2}(1 + \tanh(\sigma(1-\xi)))$ — smoothed wall function
- $r = \sqrt{x^2 + y^2 + z^2}$ — radial distance
- $\cos\theta = x/r$ — polar angle
- $\sigma$ — wall sharpness (inverse thickness)
- $v$ — dimensionless velocity ($v/c$)
- $\rho = 10$ — bubble radius (fixed)

### Diagnostics
Computed on z=0 2D slice (area integrals):
- $E^+ = \int_{\rho>0} \rho_{\rm ADM} \, dA$ — positive energy
- $E^- = \int_{\rho<0} |\rho_{\rm ADM}| \, dA$ — negative energy magnitude
- $\text{neg\_fraction} = E^- / (E^+ + E^-)$ — negative energy fraction

### Experiments

1. **1D Sigma Sweep** (`sweep_sigma_v1.5.*`)
   - Fixed $v = 1.5$, sweep $\sigma \in [1.0, 12.0]$ (25 points)
   - Grid: $n = 101$, extent = ±20 (geometric units)
   - Runtime: ~30 seconds

2. **2D Parameter Sweep** (`sweep_2d_*`)
   - $\sigma \in [2.0, 10.0]$ (15 points)
   - $v \in [0.8, 2.5]$ (15 points)
   - Grid: $n = 81$, extent = ±20
   - Total: 225 evaluations
   - Runtime: ~3 minutes

3. **Hybrid Optimization** (`optimization_refined.json`)
   - Method: Grid search (12×12) + Nelder-Mead local refinement
   - Objective: Minimize $|E^-|$
   - Grid: $n = 81$
   - Evaluations: 283 (144 grid + 139 refinement)
   - Runtime: ~3 minutes

## Key Results

### Finding 1: Balanced Energy Cancellation (Not Positive Dominance)
**neg\_fraction ≈ 0.50 across ALL tested parameters.**

This means the Rodal-like potential in this simple form produces **nearly equal positive and negative energy** — a cancellation regime, not the "predominantly positive" result claimed in some literature.

| Config | σ | v | neg\_fraction | |E⁻| |
|--------|---|---|---------------|------|
| 1D Best | 1.00 | 1.5 | 0.4999 | — |
| 2D Best | 2.00 | 1.65 | 0.4999995 | 2.04×10¹ |
| Optimizer | 2.00 | 0.80 | — | 4.80 |

### Finding 2: Lower σ and Lower v Minimize |E⁻|
- **Sharp walls** (low σ) produce less negative energy magnitude than smooth walls
- **Lower velocities** (v → 0.8) reduce |E⁻|
- But neg_fraction remains ~0.50 regardless

### Finding 3: No Parameter Regime Eliminates Negatives
Across the entire tested space:
- σ ∈ [1, 12]
- v ∈ [0.8, 2.5]

**All configurations show neg\_fraction ≈ 0.50.**

## Interpretation

The current implementation uses a **simplified Rodal ansatz** without:
1. Full 3D volume integration (current: 2D z=0 slice)
2. Tail corrections for far-field decay
3. Specific γ-correction or normalization from Rodal (2025) paper
4. Cosθ-dependent wall modulation (current: simple tanh wall)

**This explains the discrepancy with Rodal's "0.04% imbalance" claim.** The published result likely requires:
- Exact Rodal potential form (not just tanh smoothing)
- Full 3D integration with tail extrapolation
- Careful normalization choices

## Recommendations

### Immediate Next Steps
1. **Implement full Rodal (2025) potential** from `sn-article.tex`
   - Extract exact functional form with cosθ modulation
   - Verify normalization conventions
   - Compare side-by-side with this baseline

2. **Upgrade to 3D volume integrals** (M1 extension)
   - Replace 2D z=0 slice with full $dV$ integration
   - Add tail correction (M4 already implemented for 2D, extend to 3D)
   - Re-run optimization with 3D diagnostics

3. **Validate against Celmaster & Rubin energy density** (M6 complete)
   - Cross-check ADM diagnostics vs Einstein tensor eigenvalues
   - Already done for Lentz potentials; apply to Rodal-like case

### Scientific Implications
This baseline establishes that:
- Simple irrotational smoothing ≠ automatic positive energy dominance
- Specific functional forms and integration methods matter critically
- Claims of "predominantly positive" energy require careful reproduction

## Files

- `sweep_sigma_v1.5.json` — 1D sweep data
- `sweep_sigma_v1.5.png` — 1D sweep visualization (4 panels: E⁺, |E⁻|, E_net, neg_fraction)
- `sweep_2d_sigma_v.json` — 2D sweep data (225 points)
- `sweep_2d_heatmap.png` — Heatmap visualization (3 panels: |E⁻|, E⁺, neg_fraction)
- `optimization_refined.json` — Hybrid optimizer result with git provenance

## Reproducibility

```bash
# Run complete experiment suite
cd /path/to/irrotational-warp-lab

# 1D sweep
python -m irrotational_warp sweep --rho 10 --sigma-min 1.0 --sigma-max 12.0 \
  --sigma-steps 25 --v 1.5 --n 101 \
  --out results/experiments/baseline_rodal/sweep_sigma_v1.5.json \
  --out-plot results/experiments/baseline_rodal/sweep_sigma_v1.5.png

# 2D sweep
python -m irrotational_warp sweep-2d --rho 10 \
  --sigma-min 2.0 --sigma-max 10.0 --sigma-steps 15 \
  --v-min 0.8 --v-max 2.5 --v-steps 15 --n 81 \
  --out-json results/experiments/baseline_rodal/sweep_2d_sigma_v.json \
  --out-plot results/experiments/baseline_rodal/sweep_2d_heatmap.png

# Optimization
python -m irrotational_warp optimize --rho 10 \
  --sigma-min 2.0 --sigma-max 10.0 --sigma-steps 12 \
  --v-min 0.8 --v-max 2.5 --v-steps 12 --n 81 --refine \
  --out results/experiments/baseline_rodal/optimization_refined.json

# Analyze results
python scripts/analyze_experiment.py
```

## Git Provenance

All results include git metadata:
```json
{
  "git": {
    "sha": "...",
    "branch": "main",
    "dirty": "no"
  }
}
```

Run `git log` to see commit history at time of execution.
