# PhD Research Repository – PK/PD Modeling and Clinical Applications

This repository contains the code, analyses, and applications developed as part of a **PhD thesis focused on pharmacokinetic/pharmacodynamic (PK/PD) modeling**, with a particular emphasis on **therapeutic drug monitoring (TDM)** and **model-informed precision dosing (MIPD)** in clinical settings.

The work includes PK/PD analyses, supporting utilities, protocols, and Shiny applications related to **dalbavancin** and **vancomycin**, developed in R.

---

## Repository Structure
```bash
├── Dalbavancin/          # PK/PD analyses and Shiny apps related to dalbavancin
├── Vancomycin/           # PK/PD analyses and tools related to vancomycin
├── Utils/                # Shared utility functions used across projects
├── Protocol/             # Study protocols and methodological documents
├── renv/                 # Project-specific R package library (managed by renv)
├── renv.lock             # Lockfile with exact package versions
├── .Rprofile             # Project configuration (renv activation)
├── .gitignore
└── PhD.Rproj             # RStudio project file
```
---

## System Requirements

To run this project successfully, you need:

- **R** (recommended ≥ 4.3)
- **RStudio** 
- **Rtools45** (required)

### Why Rtools45 is required

Several packages used in this project (e.g. `rxode2`, `lotri`, and other modeling or C++-backed packages) require compilation from source on Windows systems.  
**Rtools45 must be installed and properly configured** before restoring the environment.

👉 Download Rtools45 from:  
https://cran.r-project.org/bin/windows/Rtools/

---

## Dependency Management with `renv`

This project uses **`renv`** to ensure **full reproducibility** of the computational environment.

### What `renv` does

- Creates a **project-specific R library**
- Locks **exact package versions** in `renv.lock`
- Ensures analyses can be reproduced across machines and over time
- Avoids conflicts with your global R installation

This is particularly important for the long-term PhD research projects and the methodological reproducibility

---

## How to Set Up the Project

### 1. Clone the repository

```bash
git clone <repository-url>
```
### 2. Restore the R environment

Run the following command once:

```r
renv::restore()
```

This will:
- Install all required packages
- Use the exact versions specified in `renv.lock`
