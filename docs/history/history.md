- Created a new public repo: https://github.com/arcticoder/irrotational-warp-lab and cloned it to irrotational-warp-lab.
- Added it to the multi-root workspace in energy/energy.code-workspace.
- Wrote the handoff-grade backlog in irrotational-warp-lab/docs/TASKS.md, copied it into my internal task tracker, and implemented the first increment (M0 + minimal M1).

What you can run now:
- From the repo root: `.python -m irrotational_warp plot-slice --out slice.png --json-out results/summary.json`
- Tests: `.python -m pytest -q`

Notes:
- Current output is an ADM/Hamiltonian-constraint-based diagnostic on a 2D slice (fast, not Rodal’s invariant eigenvalue diagnostic yet); see irrotational-warp-lab/docs/NOTES.md.
<!-- ------ -->