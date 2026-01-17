# Notes / Caveats

## Diagnostics we compute (first increment)

- `rho_adm` is computed from the 3+1 Hamiltonian constraint assuming **flat spatial slices** and **unit lapse**.
- This is a useful fast diagnostic but is **observer-dependent**; it is not the Rodal/McMonigal invariant “proper energy density”.

## Units

All quantities are currently in geometric units with $G=c=1$.

- `v` is dimensionless and represents $v/c$.
- Global integrals computed from `rho_adm` are in “geometric” energy units.

## Numerics

- Finite differences are used; grid resolution and domain size materially affect peak values.
- Use convergence checks before interpreting magnitude.

## Next

Milestone M2 will add mixed Einstein tensor eigenvalue diagnostics to match Rodal-style invariants.
