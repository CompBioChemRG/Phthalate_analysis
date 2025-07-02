# Phthalate–Androgen Receptor Interaction Analysis

This repository contains the computational workflow and datasets for the study:  
**"In Silico Exploration of the Endocrine-Disrupting Potential of Phthalate Esters through Structural Interaction with the Human Androgen Receptor"**, which investigates the endocrine-disrupting potential of phthalate esters using molecular docking, molecular dynamics (MD), MM/PBSA free energy calculations, and non-covalent interaction analysis.

## 📋 Project Description

Phthalate esters (PAEs) are environmental pollutants that may disrupt androgen receptor (AR) signaling. This project explores how 88 phthalate compounds interact with the human androgen receptor ligand-binding domain (PDB ID: 2AM9), using:

- **Multi-platform virtual screening** (DockThor, AutoDock Vina, iGEMDOCK)
- **Molecular Dynamics Simulations** (GROMACS, 500 ns)
- **Binding Free Energy Calculations** (MM/PBSA)
- **Principal Component Analysis (PCA) and Free Energy Landscape (FEL)**
- **Non-covalent Interaction Profiling** (PLIP)

## 🧪 Methodology

1. **Ligand Preparation**
   - Source: [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
   - Standardized, optimized, and protonated at physiological pH.

2. **Receptor Preparation**
   - PDB structure: `2AM9`
   - Cleaned, hydrogenated, and minimized using CHARMM-GUI.

3. **Virtual Screening**
   - Docking tools: DockThor, AutoDock Vina (via PyRx), and iGEMDOCK.
   - Top ligands selected based on consensus scoring.

4. **Molecular Dynamics**
   - Software: GROMACS v2023 with CHARMM36 force field.
   - Solvation with TIP3P water, 500 ns production run per complex.

5. **Binding Energy**
   - Tool: `gmx_MMPBSA`
   - Free energy decomposition and entropy correction included.

6. **Structural Analysis**
   - RMSD, RMSF, Rg, SASA
   - PCA & FEL
   - B-factor and non-covalent interaction visualization (PLIP)

## 📌 Requirements

- Python ≥ 3.8  
- GROMACS ≥ 2023  
- RDKit, NumPy, Pandas, Matplotlib, SciPy  
- gmx_MMPBSA, AmberTools  
- UCSF ChimeraX 
- Schrödinger Release 2024-1

##🧑‍🔬 Authors

Haruna L. Barazorda-Ccahuana
hbarazorda@ucsm.edu.pe
Computational Biology and Chemistry Research Group, Catholic University of Santa María (UCSM), Peru

##📜 Acknowledgments

This work was supported by the internal research project of the Catholic University of Santa María (UCSM), which enabled the execution of the computational analyses and research presented in this study.
