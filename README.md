# Effective Mass Calculation from Band Structure

## Description
A tool for calculating effective mass from electronic band structures. It analyzes band
structure data, counts the individual bands, identifies the Conduction Band Minimum (CBM)
and Valence Band Maximum (VBM), reports the band gap (direct/indirect), classifies the
material as metallic or semiconducting, and calculates the electron and hole effective
masses.

There are two ways to use it:

- **`effective mass calculator.html` — the browser app (recommended).** A single self-contained file that runs
  in any modern web browser on Windows, macOS, or Linux. **No Python and no installation of
  any kind.** Just open the file.
- **`effective_mass_gui.py` — the original Python/tkinter app**, kept for reference and for
  users who prefer the desktop version. Requires Python and several libraries.

## Quick start (browser app — no install)
1. Download **`effective mass calculator.html`** (from this repository).
2. **Double-click it** — it opens in your default browser (Chrome, Edge, Safari, Firefox…).
3. Click **Load sample data** to try it immediately, or drag-and-drop your own CSV / click
   **Choose CSV…**.
4. Enter your material's **lattice parameter** (in Ångström).
5. Click **Analyze**.

It works completely offline and sends no data anywhere — everything runs locally in your
browser. You can also host this file on GitHub Pages to share a link (rename it to
`index.html` first if you want it to load automatically as the site's home page).

## Features
- Zero-install, cross-platform (runs in any browser)
- Drag-and-drop CSV upload + built-in sample data
- Automatic detection of the number of bands and k-points per band
- Metal / semiconductor classification and intermediate-band detection
- VBM and CBM identification with band-gap value and **direct/indirect** classification
- Effective mass for electrons and holes, shown in **kg and as m\*/mₑ**
- Interactive band-structure plot with highlighted VBM/CBM
- Export the report (`.txt`) and the parsed bands (`.csv`)

## Input file format
A CSV file with **no header row**:

- **2 columns:** `k, energy(eV)`
- **4 columns:** `k1, E1, k2, E2` (two spin components; the second pair is concatenated
  after the first)

k-points run in Brillouin-zone units from `0` to `1`. Energies should be shifted so the
Fermi level is at `0 eV` (positive → conduction, negative → valence). Bands may be separated
by a single delimiter row, matching the original tool's data layout.

## Method
The effective mass is obtained from a finite-difference second derivative of E(k) at the band
extremum, using periodic boundary conditions at the zone edges (k=0 and k=1):

```
m* = | ħ² / (d²E/dk²) |
```

with `d²E/dk²` evaluated on a 3-point stencil and converted to SI units using the lattice
parameter. This matches the physics of the original Python implementation, so results from
the browser app and the Python version agree.

## Using the original Python version (optional)
Requires Python 3 with the libraries in `requirements`:
```bash
pip install -r requirements
python effective_mass_gui.py
```
`effective_mass_calculation.py` is a standalone script with hard-coded example data that
also serves as the reference implementation of the physics.

## Author
- **ShahnilZulkarnain**

## License
This project is licensed under the MIT License.

## Acknowledgments
- Original desktop version built with Python's scientific stack (NumPy, Pandas, Matplotlib).
- The browser app reimplements the same analysis in dependency-free JavaScript.
