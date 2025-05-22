import numpy as np

# Constants
hbar_js = 1.05457e-34  # Reduced Planck constant in J·s
eV_to_J = 1.60218e-19  # Conversion factor from eV to Joules

def calculate_effective_mass(k_points, band_energies, lattice_parameter_angstrom, find_minimum=True):
    """
    Calculate effective mass from band structure data
    
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
    
    # Calculate second derivative based on position of extremum
    if extreme_index == 0:  # Left edge
        d2E_dk2 = (band_energies[2] - 2 * band_energies[1] + band_energies[0]) / (
            (k_points[1] - k_points[0]) ** 2
        )
    elif extreme_index == len(band_energies) - 1:  # Right edge
        d2E_dk2 = (band_energies[-3] - 2 * band_energies[-2] + band_energies[-1]) / (
            (k_points[-1] - k_points[-2]) ** 2
        )
    else:  # Internal points
        d2E_dk2 = (band_energies[extreme_index + 1] - 2 * band_energies[extreme_index] + band_energies[extreme_index - 1]) / (
            (k_points[extreme_index + 1] - k_points[extreme_index - 1]) ** 2
        )
    
    # Convert to proper units
    d2E_dk2_j_m2 = d2E_dk2 * eV_to_J / ((2 * np.pi / (lattice_parameter_angstrom * 1e-10)) ** 2)
    
    # Calculate effective mass
    effective_mass_kg = abs((hbar_js ** 2) / d2E_dk2_j_m2)
    return effective_mass_kg

def main():
    # Test data with custom values
    k_points = np.array([0, 0.08333334815419560, 0.1666666913083900, 0.250000039462586, 0.3333333613083900, 0.4166666781541960, 0.5000000000000000, 0.5833333481541960, 0.6666666913083910, 0.7500000394625860, 0.8333333613083900, 0.9166666781541960, 1])
    
    # Replace these arrays with your VBM and CBM energy values
    cbm_energies = np.array([0, -0.07740795032354430, -0.2431394253189060, -0.3945471001916490, -0.4792308302161480 ,-0.6638072766906750, -0.8344887744454070, -0.6638078209183430, -0.479228109077805, -0.3945487328746540, -0.2431396974327400, -0.07740822243737890, 0])
    vbm_energies = np.array([3.446652057339760, 3.448836587201350, 3.455227996940900, 3.465460293451290, 3.474003579392320, 3.497649999477340, 3.530395089839150, 3.497651360046510 ,3.474002763050820, 3.465461654020460, 3.455230718079240, 3.448834954518340, 3.446652057339760])
    
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