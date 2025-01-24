# V-Cornea Simulation Model

## Overview
This repository hosts the CompuCell3D-based implementation of V-Cornea, an agent-based model of the corneal epithelium described in:
> Vanin J., Glazier J.A., Getz M., Knudsen T.B., Mahony C.  
> "Towards a computational model of corneal epithelium homeostasis, injury, and recovery" (*Paper details/link*)

The model simulates corneal epithelial homeostasis and injury response over time, capturing:
- Epithelial thickness maintenance (50–52 µm)
- 7–14 day turnover cycles
- Wound healing dynamics for superficial to moderate injuries
- Emergent recurrent corneal erosions with basement membrane disruption

## Installation

1. **Install CompuCell3D** (version >= 4.6.0) from [here](https://compucell3d.org/SrcBin). 
Optionally (create a env) 
before installing compucell and follow the process from there

2. **Clone this repo**:  
   ```bash
   git clone https://github.com/VaninJoel/vCornea.git 

3. **Install dependencies**

Create a conda environment:
bash
Copy
conda env create -f environment.yml
conda activate v-cornea-env
Usage
To run a basic simulation:

Launch CompuCell3D and open src/VCorneaModel.cc3d.
Press “Play” to start the simulation.
Monitor simulation progress in CompuCell3D Player.
Command-line approach (if supported):

bash
Copy
cc3drun --steppable=XMLScriptReader --parameter-file=src/VCorneaModel.cc3d
Parameter Adjustments:

EGF diffusion constants, cell volume constraints, etc., can be modified in src/params.py (see "Growth Dynamics" references in Section 2.5 of the paper).
Injury type (slight, mild, moderate) is set in simulation_config.json.
Reproducing Paper Figures
Figure 4: EGF Concentration Fields
After running VCorneaModel.cc3d, a results folder is generated in output/. Then run:

bash
Copy
python analysis/plot_egf_distribution.py --input output/run_000
This script regenerates the EGF heatmap shown in Figure 4 of the paper.

Injury Scenarios
Slight/mild/moderate injuries can be configured in src/injuries.py. The paper’s default scenario is mild injury radius = 30 µm, chemical threshold = 1.0.

Directory Structure
r
Copy
├── README.md               <- You are here
├── LICENSE
├── requirements.txt        <- Python dependencies
├── src/
│   ├── VCorneaModel.cc3d   <- Main CompuCell3D project
│   ├── params.py           <- Global parameters (EGF, growth, etc.)
│   ├── injuries.py         <- Config for ablation/chemical injuries
│   └── ...
├── analysis/
│   ├── plot_egf_distribution.py
│   ├── plot_tissue_thickness.py
│   └── ...
├── data/                   <- Sample data or smaller test files
└── output/                 <- Default output folder
Model Calibration & Validation
Calibration: Key parameters (EGF diffusion, contact energies) come from experimental data (Chan et al. 1991). Section 2.5 in the paper explains the scaling to simulation time steps.
Expected Results:
Epithelial thickness stabilizes ~50–52 µm after ~10 days (see paper Figure 5).
Slight injuries typically heal within ~3 days; moderate injuries produce recurrent erosions.
Contributing
Pull requests and issue reports are welcome. If you find a bug or want to extend the basement membrane regeneration model, open an issue or contact us at [email@example.com].

License
This project is licensed under the MIT License, unless stated otherwise.

Citation
If you use V-Cornea, please cite:

less
Copy
@article{Vanin2025VCornea,
  author  = {Vanin, Joel and Glazier, James A. and Getz, Michael and Knudsen, Thomas B. and Mahony, Catherine},
  title   = {Towards a computational model of corneal epithelium homeostasis, injury, and recovery},
  journal = {...},
  year    = {2025},
  doi     = {...}
}
vbnet
Copy

---

### Final Thoughts
- If anything is particularly unique or tricky in your build process (e.g., a custom version of CompuCell3D), that should be clearly stated.
- If there are multiple scenarios or advanced features (like coupling to other models or additional mechanical layers), point to the relevant scripts and/or sections in your paper.

**Answer to your direct question**:  
To craft the most effective README, it’s *extremely* helpful to have (1) **the paper** for accurate references to the biological motivations and parameter choices, (2) **the code** structure so we know how to instruct people on installation, usage, and how to tweak simulations, and (3) any **extra clarifications** about hidden complexities or known pitfalls. The more detail you can provide in all three areas, the stronger (and more useful) the README will be to your audience.





