# Irrotational Warp Lab — Task Backlog (handoff-ready)

## Mission
Build a reproducible, numerically robust research codebase to explore **irrotational (curl-free) shift-vector warp metrics** (Lentz-derivative / Rodal-style potentials) and quantify:

- local and global energy-condition diagnostics (WEC/NEC/DEC/SEC where applicable)
- the **size, shape, and scaling** of negative-energy regions (if any)
- parameter regimes that minimize negative contributions (e.g., wall thickness / smoothing, velocity, geometry)

Primary objective: produce **novel, defensible results** suitable for a paper.

## Research framing (working assumptions)
We start with ADM form with unit lapse and flat spatial slices:

$$ds^2 = -dt^2 + \delta_{ij}(dx^i + \beta^i dt)(dx^j + \beta^j dt)$$

Irrotational condition: $\nabla \times \beta = 0$ so we can set $\beta = \nabla \Phi$.

Rodal-style potentials are typically of the form:

$$\Phi \sim (v/c)\,\rho\, f(r/\rho)\cos\theta$$

with smoothed wall function $f(\xi)$ (e.g., tanh wall).

**Important:** many “positive energy” claims depend on *which observer / invariant is used.* We will explicitly compute multiple diagnostics:

- Eulerian energy density from 3+1 constraints (fast, but observer-dependent)
- invariant “proper energy density” from eigenvalues of $G^{\mu}{}_{\nu}$ (slower, but closer to Rodal/McMonigal framing)

## Non-goals (for now)
- Claims about realizable propulsion, engineering readiness, or experimental feasibility.
- Any “negative energy generator” design work.
- Quantum field theory backreaction or semiclassical stress tensors (unless later milestone).

## Repository layout (target)
- `docs/` — derivations, notes, validation reports
- `src/irrotational_warp/` — library code
- `scripts/` — reproducible runs (sweeps, plots, report generation)
- `notebooks/` — exploration notebooks (kept minimal; scripts are source of truth)
- `tests/` — unit/regression tests
- `data/` — small cached numerical outputs (optional; keep Git clean)

## Environment & tooling
- Python 3.10+ recommended
- Core deps: `numpy`, `scipy`, `matplotlib`
- Symbolics (optional): `sympy`
- Acceleration (optional): `numba`
- Formatting/lint: `ruff`, `black` (or just `ruff format`)

## Key inputs in the broader workspace
These are *not* copied into this repo; we reference them for reading/validation:

- Rodal (2025) LaTeX: `energy/papers/related/sn-article.tex`
- White (2025): `energy/papers/related/White_2025_Class._Quantum_Grav._42_235022.pdf`
- McMonigal et al. (2025) draft: `energy/papers/related/Comment_on_hyper-fast_solitons_draft_4.tex`
- Lentz (2021): `energy/papers/related/Lentz_2021.tex`
- Visser/Santiago/Schuster: `energy/papers/related/WarpPRD-2022-02-26.tex`
- Fuchs et al. (2024): `energy/papers/related/Fuchs_2024.tex`

## Milestones and tasks

### M0 — Project scaffold + reproducible entrypoints
**Goal:** anyone can run a baseline plot and a baseline integral in one command.

Tasks:
1. Create Python project scaffolding (`pyproject.toml`, package under `src/`).
2. Implement a CLI:
   - `plot-slice` (2D slice heatmap for chosen diagnostic)
   - `sweep` (parameter sweep grid output)
3. Add minimal tests:
   - shape/finite checks
   - regression snapshot (small grid) with tolerance

Acceptance:
- `python -m irrotational_warp plot-slice --diagnostic adm_rho --v 1.5 --rho 10 --sigma 5` writes a PNG.

### M1 — Fast 3+1 diagnostics (Eulerian observer)
**Goal:** compute fast, differentiable diagnostics suitable for optimization loops.

Core formulas (flat slices, $\alpha=1$):
- $K_{ij} = \tfrac12(\partial_i \beta_j + \partial_j \beta_i)$
- Hamiltonian constraint: $R + K^2 - K_{ij}K^{ij} = 16\pi\,\rho$; with $R=0$ gives
  $$\rho_{\rm ADM} = \frac{K^2 - K_{ij}K^{ij}}{16\pi}$$

Tasks:
1. Implement potential families:
   - Rodal-style dipole potential (axisymmetric)
   - optional: alternative smoothing families (polynomial, compact support)
2. Implement Cartesian evaluators for $\beta(x,y,z)$ and derivatives using finite differences (baseline) and optional analytic gradients.
3. Compute $K_{ij}$, $\rho_{\rm ADM}$, and track where it’s negative.
4. Global integrals:
   - $E^+ = \int_{\rho>0} \rho\, dV$
   - $E^- = \int_{\rho<0} |\rho|\, dV$
   - $E_{\rm net} = E^+ - E^-$

Acceptance:
- Produces stable values under grid refinement (convergence plot).

### M2 — Invariant diagnostics: eigenvalues of mixed Einstein tensor
**Goal:** approximate Rodal-style “proper energy density” $\rho_p$.

Two implementation tracks (do both; cross-check):

A) **Direct metric → Einstein tensor** (slower, clear):
- Construct 4D metric $g_{\mu\nu}$ from $\beta$ and compute $G^{\mu}{}_{\nu}$ numerically.
- Take eigenvalues; for Type I all eigenvalues real.

B) **3+1 stress-energy reconstruction** (faster, subtle):
- Use constraints + evolution terms with a chosen stationarity assumption to build an effective $T^{(a)}{}_{(b)}$ in an orthonormal frame.
- Then eigen-solve that 4×4 matrix.

Tasks:
1. Implement Track A for 2D axisymmetric case first (reduce cost).
2. Validate on known metrics / sanity checks:
   - Minkowski: all zeros
   - small-amplitude potential: scaling $\propto v^2$
3. Implement Type classification checks (Hawking–Ellis): ensure Type I region detection.

Acceptance:
- For at least one configuration, produces maps of eigenvalues with documented numerical error.

### M3 — Tail correction + finite-box error control
**Goal:** make global energies meaningful without absurd grid extents.

Tasks:
1. Compute radial shells of average density $\langle\rho\rangle(r)$.
2. Fit far-field decay (e.g., $\sim 1/r^4$) and extrapolate tail:
   $$E(\infty) \approx E(R) + \int_R^\infty 4\pi r^2\langle\rho\rangle(r)\,dr$$
3. Add uncertainty bars (fit residuals → tail uncertainty).

Acceptance:
- Report includes: chosen $R$, fit region, fitted exponent, tail fraction.

### M4 — Parameter sweeps + optimization
**Goal:** search for regimes minimizing negatives while keeping target kinematics.

Tasks:
1. Build sweep runner producing:
   - heatmaps of $E^-$ vs ($\sigma$, $v/c$)
   - Pareto fronts: minimize $E^-$ and peak negativity vs constraints
2. Add optimizer (start simple):
   - grid search + local Nelder–Mead refinement
   - optional Bayesian optimization later
3. Record full provenance to JSON/JSONL:
   - parameters
   - grid settings
   - diagnostics summary
   - code version (git SHA)

Acceptance:
- Reproduces same “best” config deterministically with fixed seed.

### M5 — Paper-grade validation against literature
**Goal:** “defensible claims” with cross-checks and reproduced plots.

Tasks:
1. Extract the exact potential definitions and parameter conventions from Rodal/McMonigal.
2. Ensure consistent units (geometric units vs SI; document conversions).
3. Replicate at least one key figure/metric trend:
   - sign/shape of negative regions
   - reported reduction factors (order-of-magnitude agreement)
4. Add a `docs/VALIDATION.md` with:
   - what matched
   - what didn’t
   - plausible reasons (grid, definition, invariants)

Acceptance:
- A single `scripts/reproduce_rodel.py` (name TBD) recreates a validation figure.

### M6 — Extensions (optional but high value)
Pick one after M5:

- **Modular sources / nacelle discretization:** test whether splitting sources changes $E^-$ scaling (White-style discretization but applied to irrotational flows).
- **Type-I enforcement:** explicitly constrain to regions where Type I holds globally.
- **Subluminal physical warp cross-check:** integrate ideas from Fuchs et al. (2024) to avoid energy-condition violations.

### M7 — Paper assembly pipeline
**Goal:** convert results into a paper quickly.

Tasks:
1. Create `docs/paper/` with a LaTeX skeleton.
2. Add a “results registry” in `results/` (small figures + JSON tables).
3. Script to regenerate all paper figures from raw results.

Acceptance:
- `make figures` or `python scripts/make_figures.py` regenerates the full figure set.

## Immediate next actions (this repo’s first increment)
1. Implement M0 + the minimal part of M1:
   - potential + shift
   - finite-difference derivatives
   - $K_{ij}$ + $\rho_{\rm ADM}$
   - a 2D slice plot
2. Add a tiny regression test.
3. Document numerical caveats in `docs/NOTES.md`.

## “Definition of done” for handoff
A less-capable model should be able to:
- run the CLI to generate a plot and a JSON summary
- extend the potential family by adding one function
- run tests and understand failures from the docs
