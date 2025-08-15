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

## Getting Started

### Prerequisites

To run the V-Cornea simulations, you will need to have **CompuCell3D** (version 4.6.0 or later) installed on your system. You can download it from the official website: [https://compucell3d.org](https://compucell3d.org).

### Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/VaninJoel/vCornea.git](https://github.com/VaninJoel/vCornea.git)
    ```
2.  **Navigate to the project directory:**
    ```bash
    cd vCornea
    ```

---

## Usage

You can run the simulations directly through the CompuCell3D Player or via the command line.

### Running with CompuCell3D Player

1.  Launch the CompuCell3D Player application.
2.  Go to `File` > `Open` and navigate to the `vCornea/Simulation` directory.
3.  Select the `VCornea.cc3d` project file to load the simulation.
4.  Use the player controls to `Step`, `Run`, and `Pause` the simulation. The graphical window will display the real-time state of the corneal tissue.

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

