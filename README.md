# Irrotational Warp Lab

[![CI](https://github.com/DawsonInstitute/irrotational-warp-lab/actions/workflows/ci.yml/badge.svg)](https://github.com/DawsonInstitute/irrotational-warp-lab/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/github/v/release/DawsonInstitute/irrotational-warp-lab)](https://github.com/DawsonInstitute/irrotational-warp-lab/releases)

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

**Parameter optimization:**
```bash
# Grid search + Nelder-Mead refinement (default)
python -m irrotational_warp optimize --method hybrid --refine \
  --sigma-min 2 --sigma-max 8 --v-min 0.8 --v-max 2.0 \
  --sigma-steps 10 --v-steps 10 --n 71 \
  --out results/optimization.json

# Bayesian optimization with Gaussian Process (efficient)
python -m irrotational_warp optimize --method bayes \
  --sigma-min 2 --sigma-max 8 --v-min 0.8 --v-max 2.0 \
  --n-calls 50 --n-initial 10 --random-state 42 --n 71 \
  --out results/optimization_bayes.json

# Compare methods
python scripts/test_bayesian_optimization.py
```

**Run tests:**
```bash
pytest -q
```

**Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines.

See [docs/TASKS.md](docs/TASKS.md) for the complete research plan and milestones.

## Reproducing Paper Results

All figures and results in the paper are fully reproducible. See [results/README.md](results/README.md) for the complete results registry.

### Paper Build System

```bash
# Regenerate all figures and compile paper
make all

# Individual targets
make figures              # Regenerate all figures
make paper                # Compile LaTeX document
make test                 # Run test suite
make clean                # Remove build artifacts
```

### Key Paper Results

**Figure 1: 3D Convergence Study**
```bash
python scripts/convergence_study_3d.py \
  --n-values 40,60,80,100 \
  --out results/experiments/convergence/study_3d.json

python scripts/make_paper_figures.py --figure convergence_3d
```

**Figure 2: Superluminal Velocity Sweep**
```bash
python scripts/sweep_superluminal.py --mode axisym --nx 1200 --ny 600 \
  --v-min 1.0 --v-max 3.0 --v-steps 20 --rho 5 --sigma 4 \
  --out results/experiments/superluminal/sweep_v1_to_3.json

python scripts/make_paper_figures.py --figure superluminal
```

**Figure 3: Optimization Comparison**
```bash
python scripts/test_bayesian_optimization.py  # Generates comparison data
python scripts/make_paper_figures.py --figure optimization
```

### Computational Requirements

**Typical runtimes** (Intel i7, RTX 2060 Super):
- 2D slice (n=101): ~0.1s
- 3D volume (n=60, CPU): ~2s  
- 3D volume (n=100, GPU): ~0.5s
- Bayesian optimization (50 calls, n=71): ~30s
- Full test suite (39 tests): ~5s

---

## Latest Results

**Post-M6 Implementation** (Jan 17, 2026)

Completed milestones:
- ✅ GPU acceleration (CuPy backend, 5-10× speedup)
- ✅ 3D convergence validation (tail → 0.034%, approaching Rodal's ~0.04%)
- ✅ Superluminal studies (v ≤ 3, E ∝ v² confirmed)
- ✅ Bayesian optimization (5× fewer evaluations vs grid+NM)
- ✅ Physics invariant tests (11 regression tests)
- ✅ Paper infrastructure (LaTeX skeleton + figure pipeline)

**Test Coverage**: 39 tests (100% pass rate)

See [docs/session-progress-2026-01-17-continued.md](docs/session-progress-2026-01-17-continued.md) for details.

## Paper Status

Draft manuscript: [papers/irrotational_warp_metric.tex](papers/irrotational_warp_metric.tex)

To compile: `make paper` → Output: `papers/irrotational_warp_metric.pdf`

## Citation

If you use this code or results in your research, please cite:

```bibtex
@software{irrotational_warp_lab_2026,
  author = {Sherri and collaborators},
  title = {Irrotational Warp Lab: Validation and Optimization Framework for Curl-Free Shift-Vector Warp Metrics},
  year = {2026},
  url = {https://github.com/DawsonInstitute/irrotational-warp-lab},
  version = {0.1.0}
}
```


