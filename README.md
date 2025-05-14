# Effective Mass Calculation from Band Structure

## Description
A Python-based GUI application for calculating effective mass from electronic band structures. This tool provides an intuitive interface for analyzing band structure data, determining band gaps, and calculating effective masses for both electrons and holes.

## Features
- Interactive GUI for data visualization and analysis
- Automatic band structure plotting
- Metal/Semiconductor classification
- VBM (Valence Band Maximum) and CBM (Conduction Band Minimum) identification
- Band gap calculation
- Effective mass calculation for both electrons and holes
- Clear and intuitive data presentation

## Requirements
See `requirements.txt` for detailed dependencies. Main requirements:
```bash
pip install -r requirements.txt
```

## Installation
1. Clone the repository:
```bash
git clone https://github.com/ShahnilZulkarnain/Effective-mass-Calculation-from-Band-structure.git
cd Effective-mass-Calculation-from-Band-structure
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

## Usage
1. Run the application:
```bash
python band_structure_gui.py
```

2. Use the GUI:
   - Click "Upload CSV File" to load your band structure data
   - Click "Show Band Structure" to visualize the bands
   - Use "Calculate Effective Mass" to compute effective masses
   - Use clear buttons to reset different parts of the interface

3. Input File Format:
   - CSV file with two columns (no headers)
   - First column: k-points in Brillouin zone units
   - Second column: energy values in eV

## Example Output
- Band structure plot with identified VBM and CBM
- Material classification (Metal/Semiconductor)
- Band gap value (for semiconductors)
- Effective mass values for electrons and holes

## Author
- **ShahnilZulkarnain**
- Created: 2025-05-14

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
- Built using Python's scientific computing stack (NumPy, Pandas, Matplotlib)
- Developed for electronic structure analysis in materials science
