# V-Cornea: A Computational Model of Corneal Epithelium Homeostasis, Injury, and Recovery

This repository contains the source code and simulation files for **V-Cornea**, an agent-based computational model of the corneal epithelium, as described in the manuscript:

> Vanin J, Getz M, Mahony C, Knudsen TB, & Glazier JA. (2025). V-Cornea: A computational model of corneal epithelium homeostasis, injury, and recovery. *bioRxiv*. [https://doi.org/10.1101/2025.08.11.669602](https://doi.org/10.1101/2025.08.11.669602)

The model is implemented in CompuCell3D and simulates the dynamics of corneal epithelial tissue, including homeostasis, response to injury, and subsequent recovery processes.

---

## Model Description

V-Cornea is a 2D agent-based model that simulates a radial sagittal section of the corneal limbus and peripheral cornea. It captures the complex, emergent behaviors of the corneal epithelium by defining biologically-inspired rules for individual cells governing proliferation, differentiation, migration, and death.

The model explicitly represents key cellular and structural components:
* **Cell Types**: Limbal Epithelial Stem Cells (LESCs), Basal Cells, Wing Cells, and Superficial Cells.
* **Microenvironment**: The tear film (as a source of Epidermal Growth Factor - EGF), a simplified epithelial basement membrane (EpBM), and the stroma.
* **Signaling**: A reaction-diffusion system for EGF, a critical regulator of cell proliferation and wound healing.

By integrating these components, V-Cornea bridges the gap between cellular-level behaviors and tissue-level outcomes, providing a platform for in vitro to in vivo extrapolation (IVIVE).

---

## Features

* **Tissue Homeostasis**: Simulates the self-organization and long-term maintenance of a stable, stratified corneal epithelium, reproducing physiologically accurate cell turnover rates (7-14 days).
* **Injury and Recovery**: Models tissue response to different severities of injury (slight, mild, and moderate) based on the depth of cellular damage.
* **Predictive Healing**: Accurately predicts healing timeframes of 3-5 days for injuries confined to the epithelium.
* **Pathology Simulation**: Demonstrates emergent, dynamic structural disorganization that mimics recurrent corneal erosions (RCE) following moderate injuries that disrupt the basement membrane.
* **Extensible Framework**: Built with CompuCell3D, the model is modular and can be extended to include more complex biological mechanisms, such as immune response or basement membrane regeneration.

---

## Setup and Installation

### 1. Environment Setup with Conda

The most reliable way to set up the required environment is by using Conda. If you don't have it, please install [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Once Conda is installed, open the Anaconda Prompt (on Windows) or your terminal (on macOS/Linux) and run the following commands step-by-step:

```bash
# 1. Create a new conda environment named 'v_cornea' with Python 3.12
conda create -n v_cornea python=3.12

# 2. Activate the new environment
conda activate v_cornea

# 3. Install Mamba for faster package installation (recommended)
conda install -c conda-forge mamba

# 4. Install CompuCell3D using the specified channels
mamba install -c main -c conda-forge -c compucell3d compucell3d=4.7.0
```

### 2. Install Additional Libraries

After installing CompuCell3D, you'll need a few more libraries for data analysis and plotting. While still in your `v_cornea` environment, run this command:

```bash
# Install libraries for plotting and handling Parquet files
mamba install -c conda-forge pandas fastparquet matplotlib seaborn scipy
```

### 3. Clone the Repository

Finally, clone this repository to your local machine:

```bash
git clone [https://github.com/VaninJoel/vCornea.git](https://github.com/VaninJoel/vCornea.git)
cd vCornea
```
You are now ready to run the simulations and generate the plots.

---

## Usage

This project includes two primary simulation configurations located in `Local/Project/` that correspond to the experiments in the manuscript.

### Replicating the Manuscript Results

#### Scenario 1: Tissue Self-Organization (from Stem Cells)

This simulation demonstrates the self-organization of the stratified corneal tissue starting from only limbal epithelial stem cells (LESCs).

1.  Launch the CompuCell3D Player application.
2.  Go to `File` > `Open` and navigate to the `Local/Project/paper_version_STEM_ONLY/` directory.
3.  Select the `vCornea_v2.cc3d` project file to load the simulation.
4.  Run the simulation to observe the tissue developing to a stable, homeostatic state.

#### Scenario 2: Injury Experiments

This simulation starts with a pre-organized, stable tissue, ready for injury experiments. You can replicate the different injury scenarios from the paper by modifying one parameter.

1.  Launch the CompuCell3D Player application.
2.  Go to `File` > `Open` and navigate to the `Local/Project/paper_version/` directory.
3.  Select the `vCornea_v2.cc3d` project file.
4.  Before running, open the `Local/Project/paper_version/Simulation/Parameters.py` file.
5.  Locate the `SLS_Concentration` variable and set its value according to the desired injury level:
    * **Slight Injury**: `SLS_Concentration = 750.0`
    * **Mild Injury**: `SLS_Concentration = 1500.0`
    * **Moderate Injury**: `SLS_Concentration = 2500.0`
6.  Save the `Parameters.py` file and run the simulation in the CompuCell3D Player.

### Exploring New Scenarios

You are encouraged to design your own *in silico* experiments. The `Parameters.py` file is the central control panel for the simulation. You can adjust variables such as `InjuryTime`, `InjuryType` (True for chemical, False for ablation), `SLS_Gaussian_pulse` (True for localized - droplet-like or False for uniformly distributed coat-like initial chemical field distribution), and various cell properties to explore different biological questions.

### Running with CompuCell3D Player

1.  Launch the CompuCell3D Player application inside your active conda environment using:
```bash
python -m cc3d.player5
```     
2.  Go to `File` > `Open` and navigate to the `vCornea/Local/Project` directory.
3. And select one of the scenarios mentioned above, `paper_version_STEM_ONLY` for tissue self-organization or `paper_version` for injury and recovery.
4.  Select the `VCornea.cc3d` project file to load the simulation.
5.  Use the player controls to `Step`, `Run`, and `Pause` the simulation. The graphical window will display the real-time state of the corneal tissue.
**Note you might need to resize the window containing the cell_field, and or open new windows to display other fields like EGF and SLS.


## Generating Plots

The plots presented in the manuscript can be reproduced using the analysis scripts provided in the repository. The simulation output data is saved in the `Data/Plot Scripts` directory.

## Usage

You can run the simulations directly through the CompuCell3D Player or via the command line.



### Simulation Scenarios

The model is configured to run different experimental scenarios as described in the manuscript. You can switch between scenarios by editing the parameters in the Python scripts within the `Simulation` folder. Key scenarios include:

* **Homeostasis**: The default simulation demonstrates the development and long-term stability of the corneal epithelium.
* **Injury Simulation**: To simulate an injury, you can set the `IsInjury` parameter to `True` and specify the `InjuryType` (e.g., chemical or ablation) and severity.

---

## Data and Code Availability

All model code, simulation files, analysis scripts, and data used to generate the results presented in this study are openly available.

* **GitHub Repository (Development Version)**: [https://github.com/VaninJoel/vCornea](https://github.com/VaninJoel/vCornea)
* **Zenodo Archive (Permanent Version)**: [https://doi.org/10.5281/zenodo.16764319](https://doi.org/10.5281/zenodo.16764319)

---

## Citation

If you use this model or the associated code in your research, please cite both the manuscript and the Zenodo software archive.

### Manuscript

```bibtex
@article{Vanin2025VCornea,
  author = {Vanin, Joel and Getz, Michael and Mahony, Catherine and Knudsen, Thomas B. and Glazier, James A.},
  title = {{V-Cornea: A computational model of corneal epithelium homeostasis, injury, and recovery}},
  journal = {bioRxiv},
  year = {2025},
  doi = {10.1101/2025.08.11.669602},
  url = {[https://www.biorxiv.org/content/10.1101/2025.08.11.669602v1](https://www.biorxiv.org/content/10.1101/2025.08.11.669602v1)}
}
```

### Software

```bibtex
@software{Vanin2025VCorneaCode,
  author = {Vanin, Joel},
  title = {{vCornea: A computational model of corneal epithelium homeostasis, injury, and recovery}},
  month = jun,
  year = 2025,
  publisher = {Zenodo},
  version = {v1.0.0},
  doi = {10.5281/zenodo.16764319},
  url = {[https://doi.org/10.5281/zenodo.16764319](https://doi.org/10.5281/zenodo.16764319)}
}

