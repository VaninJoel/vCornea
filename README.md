📂 Data Availability and Reproducing Results

The raw simulation data used in the manuscript is not stored in this GitHub repository because it exceeds GitHub’s file size limits.
Instead, all datasets have been archived on Zenodo and are freely available at:

https://doi.org/10.5281/zenodo.16764319

📦 Contents of the Zenodo Archive

  Each experiment is stored in its own folder (e.g., (Mild) 1500_G/Comb0/...), with multiple replicates run using the same initial parameters.
  
  Within each replicate folder, you will find:
  
  Simulation/ — numerical outputs such as:
  
  cell_count_*.csv – cell count time series
  
  surfactant_*.csv – surfactant concentration over time
  
  thickness_rep_*.parquet – tissue thickness data
  
  LatticeData/ — CompuCell3D .vtk snapshots for visual inspection or rendering.

## 📊 Reproducing the Figures from the Paper

This repository allows you to either:

1. **Plot the results from the published raw data** (simplest option)  
2. **Run the simulation yourself** and plot your own outputs  

We recommend using a **separate conda environment** so that dependencies for plotting and simulation don’t interfere with other Python projects you may have.

---

### Option 1 — 📈 Plot from Published Data (No Simulation Required)

If you only want to generate the manuscript figures from the published data, follow these steps:

#### 1. Install Miniconda or Anaconda (if you don’t have it already)
[Miniconda download](https://docs.conda.io/en/latest/miniconda.html) (lightweight option)

#### 2. Create and activate a plotting environment
```bash
conda create -n vcornea_plot python=3.10
conda activate vcornea_plot
```
