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

📊 Plotting the Results from the Paper

Inside the data folder, you can find a ready-to-use plotting utility, plot_manager.py, that can reproduce several figures from the manuscript using the raw simulation data.

Run the Plotter

No path changes are needed — the script is already set up to find the data using relative paths.

From the Plot , simply run:

python plot_manager.py


This will:

Generate longitudinal cell count plots, histograms, and Q-Q plots for replicate variability.

Generate longitudinal and spatial thickness plots for the tissue.

Figures will display interactively; you can save them manually or modify the script to save automatically.


Running Your Own Simulations

You can also run the simulations on your local machine or HPC cluster, then replace the data in the Data folder with your own outputs to:

Compare your results with the manuscript data.

Perform your own statistical analyses.

Note: The simulations are stochastic, so results will vary slightly between runs.
