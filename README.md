<p align="center">
  <img src="https://www.ukm.my/fst/wp-content/uploads/2022/01/logo_UKM_WATAN_FST_tulisanHitam_BI-01.png" alt="UKM Logo" width="200"/>
</p>

# CHARACTERIZATION OF OFF-TARGETED FENTANYL BINDING TO COMMON OPIOID RECEPTORS

**Final Year Project by Aina Rif'ah**  
Faculty of Science & Technology  
Universiti Kebangsaan Malaysia (UKM)  
Supervisor: [Prof. Firdaus Mohd Raih]  
Academic Year: 2024/2025

---

## 📌 Project Overview

This project investigates the potential off-target interactions of fentanyl with human proteins beyond the known μ-opioid receptor (μOR). By using structural motif search (ASSAM), filtering by low RMSD, and protein-ligand docking, we aim to identify unexpected protein targets that might explain the broad pharmacological effects of fentanyl.

---

## 🎯 Objectives

1. To identify known functional binding sites on fentanyl.
2. To characterize potential off-target protein interactions.
3. To determine if fentanyl may bind proteins contributing to additional clinical effects.

---

## 🧰 Tools & Methods

| Tool              | Purpose                                      |
|-------------------|----------------------------------------------|
| Python            | Data filtering, analysis                     |
| UCSF Chimera      | Protein structure extraction                 |
| ASSAM             | Structural motif search                      |
| AutoDock Vina     | Protein-ligand docking                       |
| GitHub            | Version control, result sharing              |

---

## 📁 Repository Structure

```plaintext
fentanyl_assam_fyp/
├── README.md                    <- Project overview and information
├── data/
│   ├── assam_csv_results/       <- Original ASSAM result CSVs
│   └── filtered_pdbs/           <- Filtered .pdb files (2 per protein)
├── scripts/
│   ├── filter_assam.py          <- Python script to filter CSVs
│   └── extract_chimera.py       <- Chimera script to extract structures
├── docking/
│   └── docking_results/         <- Folder for AutoDock or Vina results
└── images/
    └── ukm_logo.png             <- UKM logo used in README
