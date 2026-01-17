# Irrotational Warp Lab

Research code for exploring irrotational (curl-free) shift-vector warp metrics and energy-condition diagnostics.

## Quickstart

```bash
cd /home/echo_/Code/asciimath/irrotational-warp-lab

python -m venv .venv
. .venv/bin/activate
python -m pip install -e ".[dev]"

python -m irrotational_warp plot-slice --out results/slice.png --json-out results/summary.json
pytest -q
```

See [docs/TASKS.md](docs/TASKS.md) for the research plan and milestones.

## Latest Results

**Baseline Optimization Experiment** (Jan 16, 2026)  
[results/experiments/baseline_rodal/README.md](results/experiments/baseline_rodal/README.md)

Key finding: The simple Rodal-like potential (tanh wall, dipole) produces **balanced energy cancellation** (neg_fraction ≈ 0.50) across all tested parameters, NOT positive dominance. This establishes that functional form and integration details matter critically for reproducing literature claims.
