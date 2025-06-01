import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
from datetime import datetime

# Define the reduced Planck constant (in eV s for plotting)
hbar_evs = 6.582176807e-16  # eV s
# Define the reduced Planck constant in SI units for effective mass calculation
hbar_js = 1.05457e-34  # J s
# Conversion factor from eV to Joules
eV_to_J = 1.60218e-19
#Electron mass in kg
m_e_kg = 9.10938356e-31

class BandStructureGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Band Structure Analysis")
        self.root.geometry("800x600")
        
        # Variables to store data
        self.k_points = None
        self.bands = None
        self.total_bands = None
        self.is_metal = None
        self.conduction_bands_indices = None
        self.valence_bands_indices = None
        self.vbm_band_index = None
        self.cbm_band_index = None
        self.vbm_value = None
        self.cbm_value = None
        self.band_gap = None
        self.filename = None
        self.cbm_k_point = None
        self.vbm_k_point = None
        self.intermediate_bands_indices = None  # Add this line
        
        # Store session information
        self.session_time = "2025-05-17 18:16:07"
        self.current_user = "ShahnilZulkarnain"
        
        self.create_widgets()
        
    def create_widgets(self):
        # Create notebook for tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill='both', expand=True, padx=10, pady=5)
        
        # Main tab
        self.main_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.main_tab, text='Main')
        
        # Plot tab
        self.plot_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.plot_tab, text='Plot')
        
        # Main tab widgets
        self.create_main_tab_widgets()
        
    def create_main_tab_widgets(self):
        # File upload frame
        upload_frame = ttk.LabelFrame(self.main_tab, text="File Upload", padding=10)
        upload_frame.pack(fill='x', padx=10, pady=5)
        
        self.file_label = ttk.Label(upload_frame, text="No file selected")
        self.file_label.pack(side='left', padx=5)
        
        upload_btn = ttk.Button(upload_frame, text="Upload CSV File", command=self.upload_file)
        upload_btn.pack(side='right', padx=5)
        
        # Actions frame
        actions_frame = ttk.LabelFrame(self.main_tab, text="Actions", padding=10)
        actions_frame.pack(fill='x', padx=10, pady=5)
        
        # Create a frame for main action buttons
        main_buttons_frame = ttk.Frame(actions_frame)
        main_buttons_frame.pack(fill='x', pady=5)
        
        show_bands_btn = ttk.Button(main_buttons_frame, text="Show Band Structure", command=self.show_band_structure)
        show_bands_btn.pack(side='left', fill='x', expand=True, padx=2)
        
        calc_mass_btn = ttk.Button(main_buttons_frame, text="Calculate Effective Mass", command=self.calculate_effective_mass)
        calc_mass_btn.pack(side='left', fill='x', expand=True, padx=2)
        
        # Create a frame for clear buttons
        clear_buttons_frame = ttk.Frame(actions_frame)
        clear_buttons_frame.pack(fill='x', pady=5)
        
        clear_plot_btn = ttk.Button(clear_buttons_frame, text="Clear Plot", command=self.clear_plot)
        clear_plot_btn.pack(side='left', fill='x', expand=True, padx=2)
        
        clear_text_btn = ttk.Button(clear_buttons_frame, text="Clear Text", command=self.clear_text)
        clear_text_btn.pack(side='left', fill='x', expand=True, padx=2)
        
        clear_all_btn = ttk.Button(clear_buttons_frame, text="Clear All", command=self.clear_all)
        clear_all_btn.pack(side='left', fill='x', expand=True, padx=2)
        
        # Results frame
        self.results_frame = ttk.LabelFrame(self.main_tab, text="Results", padding=10)
        self.results_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        self.results_text = tk.Text(self.results_frame, height=10, wrap=tk.WORD)
        self.results_text.pack(fill='both', expand=True)
        
        # Add session information
        self.add_session_info()
    
    def add_session_info(self):
        """Add session information to the results text"""
        self.results_text.insert('end', f"Session started at: {self.session_time} UTC\n")
        self.results_text.insert('end', f"User: {self.current_user}\n\n")
    
    def upload_file(self):
        self.filename = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if self.filename:
            try:
                self.file_label.config(text=f"Selected: {self.filename.split('/')[-1]}")
                self.separate_bands_and_get_kpoints()
            except Exception as e:
                messagebox.showerror("Error", f"Error reading file: {str(e)}")
    
    def separate_bands_and_get_kpoints(self):
        try:
            # Read CSV without predefined column names first to check number of columns
            df = pd.read_csv(self.filename, header=None)
            
            if len(df.columns) == 2:
                # For 2-column files, proceed as before
                df.columns = ['k_point', 'energy']
                # Convert columns to float
                df = df.astype(float)
            elif len(df.columns) == 4:
                # For 4-column files, rename columns and concatenate the data
                df.columns = ['k_point1', 'energy1', 'k_point2', 'energy2']
                # Convert all columns to float
                df = df.astype(float)
                
                # Create new dataframe with concatenated data
                df_combined = pd.DataFrame({
                    'k_point': pd.concat([df['k_point1'], df['k_point2']], ignore_index=True),
                    'energy': pd.concat([df['energy1'], df['energy2']], ignore_index=True)
                })
                df = df_combined
            else:
                messagebox.showerror("Error", "CSV file must have either 2 or 4 columns")
                return
            
            # Find the index of the first occurrence of 0 and 1 in the 'k_point' column
            first_zero_index = df.loc[(df['k_point'] >= 0 - 1e-9) & (df['k_point'] <= 0 + 1e-9)].index.min()
            first_one_index = df.loc[(df['k_point'] >= 1 - 1e-9) & (df['k_point'] <= 1 + 1e-9)].index.min()
            
            if first_zero_index is None or first_one_index is None:
                messagebox.showerror("Error", "Could not find required k-points (0 and 1) in the data")
                return
            
            # The number of k-points is the difference between these indices
            num_k_points = first_one_index - first_zero_index + 1
            self.results_text.insert('end', f"Detected number of k-points per band: {num_k_points}\n")
            
            # Collect k-points for the first band
            self.k_points = df['k_point'][first_zero_index : first_one_index + 1].tolist()
            self.results_text.insert('end', "\nK-points for the bands (in Brillouin zone units):\n")
            self.results_text.insert('end', f"{self.k_points}\n")
            
            energy_values = df['energy'].tolist()
            self.bands = []
            start_index = 0
            while start_index < len(energy_values):
                end_index = start_index + num_k_points
                band = energy_values[start_index:end_index]
                # Remove the (number of k-points + 1)-th data point if it exists
                if len(band) > num_k_points:
                    band.pop()
                if band:  # Only add non-empty bands
                    self.bands.append(band)
                start_index += (num_k_points + 1)
            
            self.total_bands = len(self.bands)
            self.results_text.insert('end', f"\nTotal number of band sets created: {self.total_bands}\n")
            self.results_text.see('end')
            
        except Exception as e:
            messagebox.showerror("Error", f"Error processing file: {str(e)}\nPlease ensure all values in the CSV file are numbers.")
            self.clear_all()
            return
    
    def clear_plot(self):
        for widget in self.plot_tab.winfo_children():
            widget.destroy()
    
    def clear_text(self):
        self.results_text.delete(1.0, tk.END)
        self.add_session_info()
    
    def clear_all(self):
        self.clear_plot()
        self.clear_text()
        self.k_points = None
        self.bands = None
        self.total_bands = None
        self.is_metal = None
        self.conduction_bands_indices = None
        self.valence_bands_indices = None
        self.vbm_band_index = None
        self.cbm_band_index = None
        self.vbm_value = None
        self.cbm_value = None
        self.band_gap = None
        self.cbm_k_point = None
        self.vbm_k_point = None
        self.intermediate_bands_indices = None  # Add this line
        self.file_label.config(text="No file selected")
    
    def show_band_structure(self):
        if self.bands is None:
            messagebox.showerror("Error", "Please upload a file first")
            return
            
        self.classify_material()
        self.plot_bands()
        self.notebook.select(self.plot_tab)
    
    def show_extremum_info(self, k_point_value, carrier_type):
        """
        Display information about the extremum point and band gap type
        """
        if carrier_type == 'electron':
            self.cbm_k_point = k_point_value
            self.results_text.insert('end', f"CBM found at k = {k_point_value:.4f}\n")
        else:  # hole
            self.vbm_k_point = k_point_value
            self.results_text.insert('end', f"VBM found at k = {k_point_value:.4f}\n")
        
        # Check if we have both CBM and VBM k-points to determine band gap type
        if self.cbm_k_point is not None and self.vbm_k_point is not None:
            if abs(self.cbm_k_point - self.vbm_k_point) < 1e-6:  # Using small threshold for float comparison
                self.results_text.insert('end', "Direct Band Gap\n")
            else:
                self.results_text.insert('end', "Indirect Band Gap\n")
    
    def classify_material(self):
        self.is_metal = False
        self.conduction_bands_indices = []
        self.valence_bands_indices = []
        self.intermediate_bands_indices = [] 
        self.vbm_value = -float('inf')
        self.cbm_value = float('inf')
        self.vbm_band_index = None
        self.cbm_band_index = None
        
        for i, band in enumerate(self.bands):
            positive_energies = any(e > 1e-9 for e in band)
            negative_energies = any(e < -1e-9 for e in band)
            # Check for intermediate band (values between -0.1 and 0.1)
            if all(-0.1 <= e <= 0.1 for e in band):
                self.results_text.insert('end', f"\nBand {i+1} is an intermediate band (values between -0.1 and 0.1 eV).\n")
                self.intermediate_bands_indices.append(i) 
                continue
            if positive_energies and negative_energies:
                self.is_metal = True
                self.results_text.insert('end', f"\nBand {i+1} contains both positive and negative energy values.\n")
                break
            else:
                if not any(e > 1e-9 for e in band) and any(e < -1e-9 for e in band):  # Only negative or zero
                    self.valence_bands_indices.append(i)
                    max_energy = max(band)
                    if max_energy > self.vbm_value:
                        self.vbm_value = max_energy
                        self.vbm_band_index = i
                elif any(e > 1e-9 for e in band) and not any(e < -1e-9 for e in band):  # Only positive or zero
                    self.conduction_bands_indices.append(i)
                    min_energy = min(band)
                    if min_energy < self.cbm_value:
                        self.cbm_value = min_energy
                        self.cbm_band_index = i
        
        if not self.is_metal:
            self.band_gap = self.cbm_value - self.vbm_value if self.cbm_value != float('inf') and self.vbm_value != -float('inf') else None
            
            # Add this new section to print VBM and CBM information
            self.results_text.insert('end', "\n=== Band Structure Details ===\n")
            
            if self.vbm_band_index is not None:
                self.results_text.insert('end', f"\nVBM Band (Index {self.vbm_band_index}):\n")
                self.results_text.insert('end', "k-points: " + ", ".join(f"{k:.4f}" for k in self.k_points) + "\n")
                self.results_text.insert('end', "energies: " + ", ".join(f"{e:.4f}" for e in self.bands[self.vbm_band_index]) + "\n")
                self.results_text.insert('end', f"VBM Value: {self.vbm_value:.4f} eV\n")
                
            if self.cbm_band_index is not None:
                self.results_text.insert('end', f"\nCBM Band (Index {self.cbm_band_index}):\n")
                self.results_text.insert('end', "k-points: " + ", ".join(f"{k:.4f}" for k in self.k_points) + "\n")
                self.results_text.insert('end', "energies: " + ", ".join(f"{e:.4f}" for e in self.bands[self.cbm_band_index]) + "\n")
                self.results_text.insert('end', f"CBM Value: {self.cbm_value:.4f} eV\n")
    
    def plot_bands(self):
        # Clear previous plot
        self.clear_plot()
        
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111)
        
        if self.is_metal:
            for band_data in self.bands:
                ax.plot(self.k_points, band_data, color='blue')
            ax.legend(['Metallic Band'])
            
            plt.text(0.02, 0.98, "Material Type: METALLIC",
                    transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.8),
                    verticalalignment='top')
        else:
            # Plot conduction bands
            if self.conduction_bands_indices:
                for index in self.conduction_bands_indices:
                    color = 'yellow' if index == self.cbm_band_index else 'lightcoral'
                    label = 'CBM' if index == self.cbm_band_index else 'Conduction Band'
                    ax.plot(self.k_points, self.bands[index], color=color, 
                           label=label if index == self.cbm_band_index else None)
        
            # Plot valence bands
            if self.valence_bands_indices:
                for index in self.valence_bands_indices:
                    color = 'green' if index == self.vbm_band_index else 'cornflowerblue'
                    label = 'VBM' if index == self.vbm_band_index else 'Valence Band'
                    ax.plot(self.k_points, self.bands[index], color=color, 
                           label=label if index == self.vbm_band_index else None)
        
            # Plot intermediate bands
            if self.intermediate_bands_indices:
                for index in self.intermediate_bands_indices:
                    ax.plot(self.k_points, self.bands[index], color='gray', linestyle='--',
                           label='Intermediate Band' if index == self.intermediate_bands_indices[0] else None)
        
            # Create legend
            legend_handles = []
            legend_labels = []
            if self.conduction_bands_indices:
                legend_handles.append(Line2D([0], [0], color='yellow'))
                legend_labels.append('CBM')
                legend_handles.append(Line2D([0], [0], color='lightcoral'))
                legend_labels.append('Conduction Band')
            if self.valence_bands_indices:
                legend_handles.append(Line2D([0], [0], color='green'))
                legend_labels.append('VBM')
                legend_handles.append(Line2D([0], [0], color='cornflowerblue'))
                legend_labels.append('Valence Band')
            if self.intermediate_bands_indices:  # Add this block
                legend_handles.append(Line2D([0], [0], color='gray', linestyle='--'))
                legend_labels.append('Intermediate Band')
            
            # Remove duplicate handles and labels
            by_label = dict(zip(legend_labels, legend_handles))
            ax.legend(by_label.values(), by_label.keys())
            
            # Add band structure information
            info_text = f"Material Type: SEMICONDUCTOR\n"
            if self.vbm_value != -float('inf'):
                info_text += f"VBM: {self.vbm_value:.4f} eV\n"
            if self.cbm_value != float('inf'):
                info_text += f"CBM: {self.cbm_value:.4f} eV\n"
            if self.band_gap is not None:
                info_text += f"Band Gap: {self.band_gap:.4f} eV"
            
            plt.text(0.02, 0.98, info_text,
                    transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.8),
                    verticalalignment='top')
        
        ax.set_xlabel("K-point (Brillouin zone units)")
        ax.set_ylabel("Energy (eV)")
        ax.set_title("Energy Bands vs. K-point")
        ax.grid(True)
        
        canvas = FigureCanvasTkAgg(fig, master=self.plot_tab)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        toolbar = NavigationToolbar2Tk(canvas, self.plot_tab)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    
    def calculate_effective_mass(self):
        if self.bands is None:
            messagebox.showerror("Error", "Please upload a file first")
            return
        
        if self.is_metal:
            messagebox.showinfo("Info", "Cannot calculate effective mass for metallic materials")
            return
            
        # Create dialog for lattice parameter input
        dialog = tk.Toplevel(self.root)
        dialog.title("Enter Lattice Parameter")
        dialog.geometry("300x150")
        
        ttk.Label(dialog, text="Enter lattice parameter (Angstroms):").pack(pady=10)
        entry = ttk.Entry(dialog)
        entry.pack(pady=10)
        
        def calculate():
            try:
                lattice_parameter = float(entry.get())
                dialog.destroy()
                
                self.results_text.insert('end', "\n=== Effective Mass Calculation ===\n")
                
                if self.cbm_band_index is not None:
                    electron_mass = self.calculate_electron_effective_mass(lattice_parameter)
                    if electron_mass is not None:
                        self.results_text.insert('end', f"Electron Effective Mass ≈ {electron_mass:.2e} kg\n")
                
                if self.vbm_band_index is not None:
                    hole_mass = self.calculate_hole_effective_mass(lattice_parameter)
                    if hole_mass is not None:
                        self.results_text.insert('end', f"Hole Effective Mass ≈ {hole_mass:.2e} kg\n")
                
                self.results_text.see('end')
                
            except ValueError:
                messagebox.showerror("Error", "Please enter a valid number")
        
        ttk.Button(dialog, text="Calculate", command=calculate).pack(pady=10)
    
    def calculate_electron_effective_mass(self, lattice_parameter_angstrom):
        mass = self._calculate_effective_mass(
            self.k_points,
            self.bands[self.cbm_band_index],
            lattice_parameter_angstrom,
            find_minimum=True
        )
        if mass is not None:
            self.show_extremum_info(self.k_points[np.argmin(self.bands[self.cbm_band_index])], 'electron')
        return mass

    def calculate_hole_effective_mass(self, lattice_parameter_angstrom):
        mass = self._calculate_effective_mass(
            self.k_points,
            self.bands[self.vbm_band_index],
            lattice_parameter_angstrom,
            find_minimum=False
        )
        if mass is not None:
            self.show_extremum_info(self.k_points[np.argmax(self.bands[self.vbm_band_index])], 'hole')
        return mass
    
    def _calculate_effective_mass(self, k_points_bz, band, lattice_parameter_angstrom, find_minimum=True):
        """
        Calculate effective mass considering periodic boundary conditions
        """
        if len(k_points_bz) < 3 or lattice_parameter_angstrom <= 0:
            return None
        
        # Convert k-points from BZ units to m^-1
        k_points_m_inv = np.array(k_points_bz) * 2 * np.pi / (lattice_parameter_angstrom * 1e-10)
        
        # Find extremum point
        extreme_index = np.argmin(band) if find_minimum else np.argmax(band)
        
        # Calculate second derivative considering periodicity
        if extreme_index == 0:  # At k = 0
            # Use last point before k = 1 and first point after k = 0
            d2E_dk2 = (band[1] - 2 * band[0] + band[-2]) / (
                (k_points_bz[1] - k_points_bz[0]) ** 2
            )
            # Print points used for calculation
            self.results_text.insert('end', f"\n{'Electron' if find_minimum else 'Hole'} calculation points:\n")
            self.results_text.insert('end', f"k-points: {k_points_bz[-2]:.6f}, {k_points_bz[0]:.6f}, {k_points_bz[1]:.6f}\n")
            self.results_text.insert('end', f"energies: {band[-2]:.6f}, {band[0]:.6f}, {band[1]:.6f}\n")
            
        elif extreme_index == len(band) - 1:  # At k = 1
            # Use point before k = 1 and first point (k = 0) due to periodicity
            d2E_dk2 = (band[0] - 2 * band[-1] + band[-2]) / (
                (k_points_bz[1] - k_points_bz[0]) ** 2
            )
            # Print points used for calculation
            self.results_text.insert('end', f"\n{'Electron' if find_minimum else 'Hole'} calculation points:\n")
            self.results_text.insert('end', f"k-points: {k_points_bz[-2]:.6f}, {k_points_bz[-1]:.6f}, {k_points_bz[0]:.6f}\n")
            self.results_text.insert('end', f"energies: {band[-2]:.6f}, {band[-1]:.6f}, {band[0]:.6f}\n")
            
        else:  # Internal points
            d2E_dk2 = (band[extreme_index + 1] - 2 * band[extreme_index] + band[extreme_index - 1]) / (
                (k_points_bz[1] - k_points_bz[0]) ** 2
            )
            # Print points used for calculation
            self.results_text.insert('end', f"\n{'Electron' if find_minimum else 'Hole'} calculation points:\n")
            self.results_text.insert('end', f"k-points: {k_points_bz[extreme_index-1]:.6f}, {k_points_bz[extreme_index]:.6f}, {k_points_bz[extreme_index+1]:.6f}\n")
            self.results_text.insert('end', f"energies: {band[extreme_index-1]:.6f}, {band[extreme_index]:.6f}, {band[extreme_index+1]:.6f}\n")
        
        # Convert to proper units
        d2E_dk2_j_m2 = d2E_dk2 * eV_to_J / ((2 * np.pi / (lattice_parameter_angstrom * 1e-10)) ** 2)
        
        # Calculate effective mass
        effective_mass_kg = abs((hbar_js ** 2) / d2E_dk2_j_m2)
        return effective_mass_kg

if __name__ == "__main__":
    root = tk.Tk()
    app = BandStructureGUI(root)
    root.mainloop()