# Irrotational Warp Lab

Research code for exploring irrotational (curl-free) shift-vector warp metrics and energy-condition diagnostics.

## Quickstart

### Installation
```bash
cd /home/echo_/Code/asciimath/irrotational-warp-lab

python -m venv .venv
. .venv/bin/activate
python -m pip install -e ".[dev]"

# Optional: Install GPU acceleration (requires CUDA 12.x)
pip install cupy-cuda12x

# Verify GPU availability (optional)
python scripts/check_gpu.py
```

### Quick Examples

**Basic 2D slice visualization:**
```bash
python -m irrotational_warp plot-slice --rho 10 --sigma 3 --v 1.5 --n 101 \
  --out results/slice.png --json-out results/summary.json
```

**Reproduce Rodal (2025) exact potential (CPU):**
```bash
# Axisymmetric (2.5D) - fast
python scripts/reproduce_rodal_exact.py --mode axisym --nx 2400 --ny 1200 \
  --rho 5 --sigma 4 --v 1 --out results/rodal_axisym.json

# Full 3D - slower but complete
python scripts/reproduce_rodal_exact.py --mode 3d --n 100 \
  --rho 5 --sigma 4 --v 1 --out results/rodal_3d.json
```

**GPU-accelerated 3D (requires CuPy):**
```bash
python scripts/reproduce_rodal_exact.py --mode 3d --backend cupy --dtype float32 --n 120 \
  --rho 5 --sigma 4 --v 1 --out results/rodal_3d_gpu.json
```

**Superluminal velocity sweep:**
```bash
python scripts/sweep_superluminal.py --mode axisym --nx 1200 --ny 600 \
  --v-min 1.0 --v-max 3.0 --v-steps 20 --rho 5 --sigma 4 \
  --out results/superluminal_sweep.json

python scripts/plot_superluminal.py results/superluminal_sweep.json \
  --out results/superluminal_plot.png
```

**Run tests:**
```bash
pytest -q
```

See [docs/TASKS.md](docs/TASKS.md) for the complete research plan and milestones.

## Latest Results

**Baseline Optimization Experiment** (Jan 16, 2026)  
[results/experiments/baseline_rodal/README.md](results/experiments/baseline_rodal/README.md)

Key finding: The simple Rodal-like potential (tanh wall, dipole) produces **balanced energy cancellation** (neg_fraction ≈ 0.50) across all tested parameters, NOT positive dominance. This establishes that functional form and integration details matter critically for reproducing literature claims.
