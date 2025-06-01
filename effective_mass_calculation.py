import numpy as np

# Constants
hbar_js = 1.05457e-34  # Reduced Planck constant in J·s
eV_to_J = 1.60218e-19  # Conversion factor from eV to Joules

def calculate_effective_mass(k_points, band_energies, lattice_parameter_angstrom, find_minimum=True):
    """
    Calculate effective mass considering periodic boundary conditions
    
    Parameters:
    k_points: array of k-points in Brillouin zone units (0 to 1)
    band_energies: array of energy values in eV
    lattice_parameter_angstrom: lattice parameter in Angstroms
    find_minimum: True for electrons (CBM), False for holes (VBM)
    """
    
    # Convert k-points from BZ units to m^-1
    k_points_m_inv = np.array(k_points) * 2 * np.pi / (lattice_parameter_angstrom * 1e-10)
    
    # Find extremum point
    extreme_index = np.argmin(band_energies) if find_minimum else np.argmax(band_energies)
    
    # Calculate second derivative considering periodicity
    if extreme_index == 0:  # At k = 0
        # Use last point before k = 1 and first point after k = 0
        d2E_dk2 = (band_energies[1] - 2 * band_energies[0] + band_energies[-2]) / (
            (k_points[1] - k_points[0]) ** 2
        )
        # Print points used for calculation
        print(f"\n{'Electron' if find_minimum else 'Hole'} calculation points:")
        print(f"k-points: {k_points[-2]:.6f}, {k_points[0]:.6f}, {k_points[1]:.6f}")
        print(f"energies: {band_energies[-2]:.6f}, {band_energies[0]:.6f}, {band_energies[1]:.6f}")
        
    elif extreme_index == len(band_energies) - 1:  # At k = 1
        # Use point before k = 1 and first point (k = 0) due to periodicity
        d2E_dk2 = (band_energies[0] - 2 * band_energies[-1] + band_energies[-2]) / (
            (k_points[1] - k_points[0]) ** 2
        )
        # Print points used for calculation
        print(f"\n{'Electron' if find_minimum else 'Hole'} calculation points:")
        print(f"k-points: {k_points[-2]:.6f}, {k_points[-1]:.6f}, {k_points[0]:.6f}")
        print(f"energies: {band_energies[-2]:.6f}, {band_energies[-1]:.6f}, {band_energies[0]:.6f}")
        
    else:  # Internal points
        d2E_dk2 = (band_energies[extreme_index + 1] - 2 * band_energies[extreme_index] + band_energies[extreme_index - 1]) / (
            (k_points[1] - k_points[0]) ** 2
        )
        # Print points used for calculation
        print(f"\n{'Electron' if find_minimum else 'Hole'} calculation points:")
        print(f"k-points: {k_points[extreme_index-1]:.6f}, {k_points[extreme_index]:.6f}, {k_points[extreme_index+1]:.6f}")
        print(f"energies: {band_energies[extreme_index-1]:.6f}, {band_energies[extreme_index]:.6f}, {band_energies[extreme_index+1]:.6f}")
    
    # Convert to proper units
    d2E_dk2_j_m2 = d2E_dk2 * eV_to_J / ((2 * np.pi / (lattice_parameter_angstrom * 1e-10)) ** 2)
    
    # Calculate effective mass
    effective_mass_kg = abs((hbar_js ** 2) / d2E_dk2_j_m2)
    return effective_mass_kg

def main():
    # Test data with custom values
    k_points = np.array([0.0000, 0.0833, 0.1667, 0.2500, 0.3333, 0.4167, 0.5000, 0.5833, 0.6667, 0.7500, 0.8333, 0.9167, 1.0000])
                         
    # Replace these arrays with your VBM and CBM energy values
    cbm_energies = np.array([3.4467, 3.4488, 3.4552, 3.4655, 3.4740, 3.4976, 3.5304, 3.4977, 3.4740, 3.4655, 3.4552, 3.4488, 3.4467])
    vbm_energies = np.array([0.0000, -0.0774, -0.2431, -0.3945, -0.4792, -0.6638, -0.8345, -0.6638, -0.4792, -0.3945, -0.2431, -0.0774, 0.000])

    # Lattice parameter - replace with your material's value
    lattice_parameter = 5.866  # Angstroms
    
    # Calculate effective masses
    electron_mass = calculate_effective_mass(k_points, cbm_energies, lattice_parameter, find_minimum=True)
    hole_mass = calculate_effective_mass(k_points, vbm_energies, lattice_parameter, find_minimum=False)
    
    # Print results
    print("\nTest Results:")
    print("-" * 50)
    print("K-points:")
    for k in k_points:
        print(f"{k:.12f}", end=", ")
    print("\n\nCBM energies:")
    for e in cbm_energies:
        print(f"{e:.12f}", end=", ")
    print("\n\nVBM energies:")
    for e in vbm_energies:
        print(f"{e:.12f}", end=", ")
    print("\n" + "-" * 50)
    print(f"Electron effective mass: {electron_mass:.6e} kg")
    print(f"Hole effective mass: {hole_mass:.6e} kg")
    
    # Plot the bands
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 6))
        plt.plot(k_points, cbm_energies, 'r-', label='CBM')
        plt.plot(k_points, vbm_energies, 'b-', label='VBM')
        plt.xlabel('k-point (Brillouin zone units)')
        plt.ylabel('Energy (eV)')
        plt.title('Test Band Structure')
        plt.grid(True)
        plt.legend()
        plt.show()
    except ImportError:
        print("Matplotlib not installed. Skipping plot.")

if __name__ == "__main__":
    main()