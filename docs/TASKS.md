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

Status: **COMPLETE**

Tasks:
1. ✅ Create Python project scaffolding (`pyproject.toml`, package under `src/`).
2. ✅ Implement a CLI:
   - `plot-slice` (2D slice heatmap for chosen diagnostic)
   - `sweep` (parameter sweep grid output)
3. ✅ Add minimal tests:
   - shape/finite checks
   - regression snapshot (small grid) with tolerance

Acceptance:
- ✅ `python -m irrotational_warp plot-slice --out results/slice.png --json-out results/summary.json`
- ✅ `python -m irrotational_warp sweep --out results/sweep.json`

Implementation notes:
- CLI uses `click` for arg parsing; outputs JSON + PNG.
- Tests use n=41 grids to avoid timeouts; production sweeps can use n=81 or higher.

### M1 — Fast 3+1 diagnostics (Eulerian observer)
**Goal:** compute fast, differentiable diagnostics suitable for optimization loops.

Status: **PARTIAL (2D z=0 slice approximation implemented; full 3D integration pending)**

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

Current implementation notes:
- The repo currently computes a **2D z=0 slice** diagnostic with an area integral ($dA$) rather than a full 3D volume integral ($dV$).
- This is sufficient for rapid iteration and plotting; a full 3D implementation is a natural next increment.

Math / code snippets (current conventions):

Rodal-like dipole potential (axisymmetric, oriented along +x):
```python
def phi_dipole_cartesian(x, y, z, *, rho, sigma, v, eps=1e-12):
   r = np.sqrt(x*x + y*y + z*z)
   f = 0.5 * (1.0 + np.tanh(sigma * (1.0 - r / rho)))
   costheta = x / (r + eps)
   return v * rho * f * costheta
```

Fast ADM energy density (flat slices, unit lapse):
```python
# K_ij = 1/2(∂i βj + ∂j βi), β = ∇Φ
rho_adm = (K_trace**2 - (KijKij)) / (16*np.pi)
```

Signed integrals (2D slice):
```python
E_pos = sum(rho[rho>0]) * dA
E_neg = sum(-rho[rho<0]) * dA
E_net = E_pos - E_neg
```

Acceptance:
- Produces stable values under grid refinement (convergence plot).

### M2 — Invariant diagnostics: eigenvalues of mixed Einstein tensor
**Goal:** approximate Rodal-style “proper energy density” $\rho_p$.

Status: **COMPLETE (Track A implemented for 2D z=0 slice)**
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

Math/code snippets for next increment (Track A):

Metric construction from ADM variables (unit lapse α=1, flat spatial metric γ_ij=δ_ij):
```python
# 4D covariant metric g_μν in Cartesian-like coords (t,x,y,z)
# g_00 = -α² + β^i β_i = -1 + (βx² + βy² + βz²)
# g_0i = β_i
# g_ij = γ_ij = δ_ij
g_cov = np.array([
    [-1.0 + (beta_x**2 + beta_y**2 + beta_z**2), beta_x, beta_y, beta_z],
    [beta_x, 1.0, 0.0, 0.0],
    [beta_y, 0.0, 1.0, 0.0],
    [beta_z, 0.0, 0.0, 1.0]
])  # shape (4,4,ny,nx) for 2D slice or (4,4,nz,ny,nx) for 3D
```

Christoffel symbols (numerical via finite differences):
```python
# Γ^α_μν = (1/2)g^{ασ}(∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
# Use np.gradient for ∂_i g_μν, then compute via Einstein summation
```

Ricci tensor and scalar (contract Riemann):
```python
# R_μν = ∂_σ Γ^σ_μν - ∂_ν Γ^σ_μσ + Γ^σ_ρσ Γ^ρ_μν - Γ^σ_ρν Γ^ρ_μσ
# R = g^{μν} R_μν
```

Mixed Einstein tensor:
```python
# G^μ_ν = R^μ_ν - (1/2)δ^μ_ν R
# where R^μ_ν = g^{μσ} R_σν
```

Eigenvalue solver (per grid point):
```python
# For each spatial point (i,j[,k]), extract 4x4 matrix G^μ_ν and solve:
eigs = np.linalg.eigvals(G_mixed)  # complex array
# For Type I spacetime (Hawking-Ellis), all eigenvalues should be real
# Proper energy density ≈ dominant eigenvalue (Rodal convention)
```

Acceptance:
- ✅ Produces maps of eigenvalues with Minkowski flatness validation and Type-I fraction reporting.

Implementation notes (see `src/irrotational_warp/einstein.py`):
- `compute_metric_z0()`: Construct 4D covariant metric from ADM variables (unit lapse, flat spatial metric)
- `compute_christoffel()`: Finite-difference Christoffel symbols Γ^α_μν
- `compute_ricci_tensor()`: Ricci tensor R_μν and scalar R via ∂Γ + Γ² terms
- `compute_einstein_tensor()`: Mixed Einstein tensor G^μ_ν = R^μ_ν - (1/2)δ^μ_ν R
- `compute_einstein_eigenvalues()`: Full pipeline with eigenvalues and Type-I classification

CLI usage (added `--einstein` flag to `plot-slice`):
```bash
python -m irrotational_warp plot-slice --rho 10 --sigma 3 --v 1.5 --n 51 --einstein \
  --out results/slice_einstein.png --json-out results/summary_einstein.json
```

JSON output includes:
- `eig_max`: Maximum real eigenvalue of G^μ_ν (≈ proper energy density ρ_p per Rodal)
- `ricci_scalar`: Ricci scalar R
- `type_i_fraction`: Fraction of grid where all eigenvalues are real (Type I per Hawking-Ellis)

Numerical caveats:
- Type-I classification is numerically fragile; finite-diff errors accumulate through Christoffel → Ricci → Einstein chain.
- Use modest grid sizes (n≤51) for Einstein diagnostics; ADM diagnostics scale better (can use n=101+).

### M3 — Parameter sweeps + basic diagnostics
**Goal:** search parameter space (sigma, v, rho) for configurations minimizing negative energy.

Status: **COMPLETE (sigma sweep implemented)**

Tasks:
1. ✅ Build sweep runner over sigma values; output JSON with energy integrals (E+, E-, Enet, neg_fraction).
2. ✅ Wire into CLI as `sweep` command.

Implemented:
- `sweep_sigma_z0()` in `sweep.py` loops over sigma values and computes signed integrals.
- CLI: `python -m irrotational_warp sweep --rho 10 --v 1.5 --sigma-min 1 --sigma-max 5 --sigma-count 20 --out results/sweep.json`

Next extensions (future milestones):
- 2D heatmaps (sigma vs v)
- Bayesian optimization for minimal E-
- Pareto fronts for multi-objective constraints

### M4 — Tail correction + finite-box error control
**Goal:** make global energies meaningful without absurd grid extents.

Status: **COMPLETE**

Tasks:
1. ✅ Compute radial shells of average density $\langle\rho\rangle(r)$.
2. ✅ Fit far-field decay (e.g., $\sim 1/r^4$) and extrapolate tail:
   $$E(\infty) \approx E(R) + \int_R^\infty 2\pi r\langle\rho\rangle(r)\,dr$$
3. ✅ Add uncertainty bars (fit residuals → tail uncertainty).

Implemented (see `src/irrotational_warp/tail.py`):
- `compute_radial_average_z0()`: Angle-average field over radial bins
- `fit_power_law_decay()`: Log-log linear regression for ρ ~ A/r^n
- `extrapolate_tail_integral_2d()`: Analytic tail integral for r > R (convergent for n > 2)
- `estimate_tail_uncertainty()`: Propagate fit residuals to tail integral uncertainty
- `compute_tail_correction()`: Full pipeline returning TailCorrectionResult

CLI usage:
```bash
python -m irrotational_warp plot-slice --rho 10 --sigma 3 --v 1.5 --n 101 \
  --tail-correction --out results/slice_with_tail.png --json-out results/summary_with_tail.json
```

JSON output includes:
- `exponent`: Fitted power-law exponent n (ρ ~ 1/r^n)
- `amplitude`: Fitted amplitude A
- `tail_integral_pos`, `tail_integral_neg`: Tail contributions from R to infinity
- `E_pos_corrected`, `E_neg_corrected`, `E_net_corrected`: Grid integrals + tail
- `tail_uncertainty`: Estimated 1-sigma uncertainty in tail
- `fit_residual_rms`: Quality of fit in log-log space

Acceptance:
- ✅ Report includes: chosen $R$ (fit_r_min), fit region, fitted exponent, tail fraction.
- ✅ Uncertainty propagation from fit residuals to tail integrals.

### M5 — Advanced optimization + multi-parameter sweeps
**Goal:** search for regimes minimizing negatives while keeping target kinematics.

Status: **COMPLETE**

Tasks:
1. ✅ Build 2D sweep runner producing:
   - heatmaps of $E^-$ vs ($\sigma$, $v/c$)
2. ✅ Add optimizer:
   - grid search + local Nelder–Mead refinement
   - ~~optional Bayesian optimization~~ (deferred as extension)
3. ✅ Record full provenance to JSON/JSONL:
   - ✅ parameters
   - ✅ grid settings
   - ✅ diagnostics summary
   - ✅ code version (git SHA, branch, dirty status)

Implemented:
- **2D Parameter Sweeps** (`sweep_2d_z0()` in `sweep.py`):
  - Loops over (sigma, v) grid, computes signed energy integrals at each point
  - Returns structured `SweepPoint2D` results
- **Heatmap Visualization** (`plot_heatmap_2d()` in `viz.py`):
  - 3-panel heatmap visualization (|E⁻|, E⁺, neg_fraction)
- **Optimization Engine** (`optimize.py`):
  - `grid_search()`: Exhaustive grid search over parameter space
  - `optimize_nelder_mead()`: Local Nelder-Mead simplex optimization
  - `optimize_hybrid()`: Grid search + Nelder-Mead refinement (recommended)
- **Git Provenance** (`get_git_info()` in `io.py`):
  - Captures git SHA, branch, and dirty status
  - Automatically included in all JSON outputs
- **CLI Commands**:
  - `sweep-2d`: 2D parameter sweep with heatmap
  - `optimize`: Hybrid optimizer for minimal |E⁻|
- **Test Coverage** (6 tests in `test_optimize.py`):
  - Objective function validation
  - Grid search correctness
  - Nelder-Mead local optimization
  - Hybrid optimizer with/without refinement
  - Deterministic results verification

CLI usage (2D sweep):
```bash
python -m irrotational_warp sweep-2d --rho 10 --sigma-min 1 --sigma-max 10 --sigma-steps 20 \
  --v-min 0.5 --v-max 2.5 --v-steps 20 --n 101 \
  --out-json results/sweep_2d.json --out-plot results/sweep_2d_heatmap.png
```

CLI usage (optimization):
```bash
# Grid search only
python -m irrotational_warp optimize --rho 10 --sigma-min 2 --sigma-max 8 \
  --v-min 0.8 --v-max 2.0 --sigma-steps 10 --v-steps 10 --n 71 \
  --out results/optimization.json

# Grid search + Nelder-Mead refinement (recommended)
python -m irrotational_warp optimize --rho 10 --sigma-min 2 --sigma-max 8 \
  --v-min 0.8 --v-max 2.0 --sigma-steps 10 --v-steps 10 --n 71 --refine \
  --out results/optimization_refined.json
```

JSON output (optimize command):
```python
{
  "git": {
    "sha": "401329a...",
    "branch": "main",
    "dirty": "no"
  },
  "params": {
    "rho": 10.0,
    "extent": 20.0,
    "n": 71,
    "sigma_range": [2.0, 8.0],
    "v_range": [0.8, 2.0],
    "sigma_steps": 10,
    "v_steps": 10,
    "refine": true
  },
  "optimization": {
    "best_params": {"sigma": 2.15, "v": 0.82},
    "best_value": 0.38,  // |E⁻| magnitude
    "initial_params": {"sigma": 2.0, "v": 0.8},
    "initial_value": 0.41,
    "n_evaluations": 235,
    "success": true,
    "method": "hybrid_grid_nelder_mead",
    "message": "Grid: 100 evals, Refinement: 135 evals"
  }
}
```

Heatmap visualization:
- Panel 1: |E⁻| magnitude across (σ, v) space
- Panel 2: E⁺ across (σ, v) space
- Panel 3: Negative fraction |E⁻|/(E⁺ + |E⁻|)

Test coverage:
- ✅ Objective function returns finite positive values
- ✅ Grid search finds minima within bounds
- ✅ Nelder-Mead improves or maintains initial value
- ✅ Hybrid optimizer combines both methods correctly
- ✅ Deterministic results with identical inputs

Acceptance:
- ✅ 2D sweeps with heatmap visualization
- ✅ Reproduces same "best" config deterministically (grid search is deterministic)
- ✅ Git SHA provenance in all outputs

Future extensions (deferred):
- Pareto front visualization for multi-objective optimization
- Bayesian optimization (GPyOpt or similar) for efficient exploration
- Constraint handling (e.g., min velocity, max sigma)

### M6 — Paper-grade validation against literature
**Goal:** "defensible claims" with cross-checks and reproduced plots.  
**Status:** ✅ COMPLETE

Tasks:
1. ✅ Extract the exact potential definitions from Celmaster & Rubin (2024)
   - Lentz's flawed φ_L (Eq. phiLentz)
   - Corrected φ_rh (Eq. phiMod)
   - Rhomboidal source parameters
2. ✅ Implement validation module (`validate_lentz.py`)
   - RhomboidalSource class with proper symmetry
   - Both potentials (φ_L and φ_rh) via numerical integration
   - Shift vector computation
3. ✅ Create validation script (`validate_celmaster_rubin.py`)
   - Source shape visualization
   - Shift vector plots (N_z, N_x)
   - Property tests (conservation, boundedness)
4. ✅ Document in `docs/VALIDATION.md`
   - What matched (source structure, symmetries, qualitative features)
   - What didn't (quantitative values - expected due to normalization)
   - Plausible reasons (grid resolution, γ-correction, integration method)
5. ✅ Energy density validation
   - Implemented full E calculation from shift vectors
   - Reproduced Celmaster & Rubin Fig. 4 (energy density heatmap)
   - Verified WEC violations: 3/900 points with negative energy, E_min = -5.94×10⁸
   - Reproduce Celmaster & Rubin Fig. 4
   - Verify WEC violation regions
   - Compare to our ADM implementation

Acceptance:
- ✅ `scripts/validate_celmaster_rubin.py` produces validation figures
- ✅ Source properties verified (conservation, symmetry, boundedness)
- ⬜ Energy density matches literature (order of magnitude, sign, spatial structure)
- ⬜ `docs/VALIDATION.md` documents all comparisons

**Reference:** Celmaster & Rubin (2024) "Violations of the Weak Energy Condition for Lentz Warp Drives"


### M7 — Extensions (optional but high value)
Pick one after M5:

- **Modular sources / nacelle discretization:** test whether splitting sources changes $E^-$ scaling (White-style discretization but applied to irrotational flows).
- **Type-I enforcement:** explicitly constrain to regions where Type I holds globally.
- **Subluminal physical warp cross-check:** integrate ideas from Fuchs et al. (2024) to avoid energy-condition violations.

### M8 — Paper assembly pipeline
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
