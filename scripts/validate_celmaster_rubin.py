#!/usr/bin/env python3
"""
Validation against Celmaster & Rubin (2024).

Reproduces key results from "Violations of the Weak Energy Condition 
for Lentz Warp Drives" to validate our implementation.

Cross-checks:
1. Rhomboidal source shape (Fig. 1 in Celmaster & Rubin)
2. Shift vector components N_z, N_x (Figs. 2 in C&R)
3. Energy density distributions (Figs. 3-4 in C&R)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from irrotational_warp.validate_lentz import (
    RhomboidalSource,
    phi_L_lentz,
    phi_rh_corrected,
    compute_shift_vector,
    compute_energy_density_at_point,
)


def plot_rhomboidal_source():
    """Reproduce Celmaster & Rubin Fig. 1: rhomboidal source."""
    source = RhomboidalSource()
    
    # Grid in (s, z) plane at y=0
    s_vals = np.linspace(-3, 5, 300)
    z_vals = np.linspace(-3, 4, 300)
    
    rho_grid = np.zeros((len(z_vals), len(s_vals)))
    
    for i, z in enumerate(z_vals):
        for j, s in enumerate(s_vals):
            rho_grid[i, j] = source.eval_rho(s, z)
    
    # Unscaled (divide by W factor would require N_z(0,0,0))
    # For now just plot raw source
    
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.contourf(s_vals, z_vals, rho_grid, levels=20, cmap='RdBu_r')
    ax.set_xlabel('s = |x| + |y|')
    ax.set_ylabel('z')
    ax.set_title('Rhomboidal Source ρ(s,z)\n(cf. Celmaster & Rubin Fig. 1)')
    plt.colorbar(im, ax=ax, label='ρ')
    
    out_dir = Path('validation_output')
    out_dir.mkdir(exist_ok=True)
    plt.savefig(out_dir / 'source_rhomboidal.png', dpi=150, bbox_inches='tight')
    print(f"✓ Saved source plot to {out_dir / 'source_rhomboidal.png'}")
    plt.close()


def plot_shift_vectors():
    """Reproduce Celmaster & Rubin Fig. 2: shift vector components."""
    source = RhomboidalSource()
    
    # Grid in (x, z) plane at y=0
    x_vals = np.linspace(-3, 5, 60)  # Coarser for speed
    z_vals = np.linspace(-3, 4, 60)
    
    print("Computing shift vectors (this may take a minute)...")
    
    N_z_grid = np.zeros((len(z_vals), len(x_vals)))
    N_x_grid = np.zeros((len(z_vals), len(x_vals)))
    
    for i, z in enumerate(z_vals):
        for j, x in enumerate(x_vals):
            # Use corrected potential phi_rh
            N_x, N_y, N_z = compute_shift_vector(phi_rh_corrected, x, 0.0, z, source)
            N_z_grid[i, j] = N_z
            N_x_grid[i, j] = N_x
        
        if (i + 1) % 10 == 0:
            print(f"  Row {i+1}/{len(z_vals)}")
    
    # Plot N_z
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    im1 = ax1.contourf(x_vals, z_vals, N_z_grid, levels=20, cmap='RdBu_r')
    ax1.set_xlabel('x')
    ax1.set_ylabel('z')
    ax1.set_title('N_z(x, y=0, z)\n(cf. Celmaster & Rubin Fig. 2 left)')
    plt.colorbar(im1, ax=ax1, label='N_z')
    
    # Plot N_x
    im2 = ax2.contourf(x_vals, z_vals, N_x_grid, levels=20, cmap='RdBu_r')
    ax2.set_xlabel('x')
    ax2.set_ylabel('z')
    ax2.set_title('N_x(x, y=0, z)\n(cf. Celmaster & Rubin Fig. 2 right)')
    plt.colorbar(im2, ax=ax2, label='N_x')
    
    out_dir = Path('validation_output')
    plt.savefig(out_dir / 'shift_vectors.png', dpi=150, bbox_inches='tight')
    print(f"✓ Saved shift vectors to {out_dir / 'shift_vectors.png'}")
    plt.close()


def test_source_properties():
    """
    Test source properties from Celmaster & Rubin:
    1. Linear cancellation condition
    2. Bounded domain
    3. Total integrated charge = 0
    """
    print("\n=== Testing Source Properties ===")
    
    source = RhomboidalSource()
    
    # Test 1: Bounded domain
    print("\n1. Bounded Domain Test:")
    print(f"   Source should be zero for |z| > {np.max(np.abs(source.xi)) + source.L_z}")
    print(f"   Source should be zero for |s| > {np.max(np.abs(source.beta)) + source.L_x}")
    
    # Sample outside bounds
    rho_far = source.eval_rho(10.0, 0.0)
    print(f"   ρ(10, 0) = {rho_far:.6f} (should be ~0)")
    
    # Test 2: Charge conservation
    print("\n2. Charge Conservation:")
    print(f"   Explicit sum of table α_i = {np.sum(source.alpha):.6f}")
    print(f"   Note: Symmetry creates 7 rhomboids from 4 table entries")
    print(f"   With symmetry: 2×25 + 2×(-25) + 2×(-25) + 1×50 = 0")
    
    # Verify by integrating total charge
    s_vals = np.linspace(-5, 5, 500)
    z_vals = np.linspace(-4, 4, 500)
    total_charge = 0.0
    for s in s_vals:
        for z in z_vals:
            total_charge += source.eval_rho(s, z) * (s_vals[1] - s_vals[0]) * (z_vals[1] - z_vals[0])
    print(f"   Integrated total charge: {total_charge:.6f} (should be ~0)")
    
    # Test 3: Linear cancellation (approximate numerical check)
    print("\n3. Linear Cancellation (diagonal integral):")
    s_vals = np.linspace(-5, 5, 500)
    alpha = 0.0  # offset
    integral_plus = sum(source.eval_rho(s, alpha + s) * (s_vals[1] - s_vals[0]) 
                        for s in s_vals)
    integral_minus = sum(source.eval_rho(s, alpha - s) * (s_vals[1] - s_vals[0]) 
                         for s in s_vals)
    print(f"   ∫ ρ(s', α+s') ds' = {integral_plus:.6f} (should be ~0)")
    print(f"   ∫ ρ(s', α-s') ds' = {integral_minus:.6f} (should be ~0)")


def plot_energy_density():
    """Reproduce Celmaster & Rubin Fig. 4: energy density for corrected potential."""
    source = RhomboidalSource()
    
    print("Computing energy density on coarse grid (this will take several minutes)...")
    print("Note: Energy computation requires 2nd derivatives, so it's computationally expensive")
    
    # Use a coarse grid for demonstration (Fig. 4 uses finer resolution)
    # Grid in (s, z) plane at y=0
    s_vals = np.linspace(-2, 4, 30)  # Coarse grid
    z_vals = np.linspace(-2, 3, 30)
    
    E_grid = np.zeros((len(z_vals), len(s_vals)))
    
    total_points = len(z_vals) * len(s_vals)
    count = 0
    
    for i, z in enumerate(z_vals):
        for j, s in enumerate(s_vals):
            # Convert (s, z) to (x, y, z) with y=0
            # s = x + |y| = x (at y=0)
            x = s
            y = 0.0
            
            # Compute energy density at this point
            E_grid[i, j] = compute_energy_density_at_point(
                phi_rh_corrected, x, y, z, source, h=1e-4
            )
            
            count += 1
            if count % 100 == 0:
                print(f"  {count}/{total_points} points computed...")
    
    print(f"  {total_points}/{total_points} points computed ✓")
    
    # Plot energy density
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Use symmetric colormap centered at zero to show negative energy regions
    vmax = np.abs(E_grid).max()
    im = ax.contourf(s_vals, z_vals, E_grid, levels=20, 
                     cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    
    ax.set_xlabel('s = x + |y| (at y=0)')
    ax.set_ylabel('z')
    ax.set_title('Energy Density E (Celmaster & Rubin corrected potential φ_rh)')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('E (energy density)')
    
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='k', linewidth=0.5, alpha=0.5)
    
    # Save
    output_dir = Path("validation_output")
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / "energy_density.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Saved energy density to {output_file}")
    
    # Report WEC violation
    negative_mask = E_grid < 0
    if np.any(negative_mask):
        num_negative = np.sum(negative_mask)
        total = E_grid.size
        print(f"\n✓ WEC violation detected:")
        print(f"  {num_negative}/{total} points have negative energy density")
        print(f"  E_min = {E_grid.min():.6e}")
    else:
        print(f"\n✗ No WEC violation detected (E_min = {E_grid.min():.6e})")


def main():
    """Run all validation checks."""
    print("=" * 60)
    print("Validation: Celmaster & Rubin (2024)")
    print("'Violations of the Weak Energy Condition for Lentz Warp Drives'")
    print("=" * 60)
    
    # Test source properties
    test_source_properties()
    
    print("\n=== Generating Validation Plots ===\n")
    
    # Plot source
    print("1. Plotting rhomboidal source...")
    plot_rhomboidal_source()
    
    # Plot shift vectors
    print("\n2. Computing and plotting shift vectors...")
    plot_shift_vectors()
    
    # Plot energy density
    print("\n3. Computing and plotting energy density...")
    plot_energy_density()
    
    print("\n" + "=" * 60)
    print("✓ Validation complete!")
    print("  See validation_output/ for plots")
    print("=" * 60)


if __name__ == '__main__':
    main()
