# Covalent-Docking-Benchmark
Supporting data for benchmarking covalent docking tools and non-covalent docking-based optimization.

## 1. Overview

This repository contains structural data and analysis scripts associated with the research article:

**"Benchmarking Covalent Docking Tools: Enhancing Pose Prediction through Pre-reaction Non-covalent Docking"**

The repository provides docking results generated using different covalent docking programs, structural data from pre-reaction non-covalent docking and subsequent covalent docking, and Python scripts for ligand RMSD calculation.

This study evaluates the pose prediction performance of different covalent docking methods by comparing predicted ligand conformations with the corresponding reference crystal structures.

## 2. Data Description

### 2.1 docking/

This directory contains docking results and prepared structures generated using three covalent docking programs:

- **AutoDock/**: Docking results and related structural files generated using AutoDock4.
- **CovDock/**: Docking results and related structural files generated using Schrödinger CovDock.
- **MOE/**: Docking results and related structural files generated using MOE.

Each directory contains four datasets identified by the numbers `165`, `33`, `70`, and `88`.

### 2.2 noncov_opt/

This directory contains structural data associated with the pre-reaction non-covalent docking and subsequent covalent docking workflow.

- **noncovalent_docking/**: Input files and output structures from pre-reaction non-covalent docking.
- **covalent_docking/**: Structures obtained from subsequent covalent docking based on the non-covalent docking poses, as well as structures after covalent bond formation.

### 2.3 Scripts/

This directory contains computational analysis scripts and related files used in this study.

- **rmsd/**: Python scripts for calculating ligand RMSD from AutoDock4, MOE, and CovDock docking results.
- **qmmm_opt/**: Computational files related to QM/MM optimization and molecular dynamics simulations.


