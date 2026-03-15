# BARNO

**B**atch-**A**ware **R**egulatory **N**etwork **O**ptimization

BARNO is an R package for batch-aware optimization of transcription factor regulatory weights in single-cell transcriptomic analysis. The framework is designed to reduce sample-group-specific or batch-driven transcription factor signals while preserving biologically coherent regulatory structure.

<p align="center">
<img src="figures/pipeline.png" width="800">
</p>

<p align="center">
<img src="figures/calculation.png" width="800">
</p>

## Overview

BARNO builds on GENIE3-based regulatory network inference and refines transcription factor–target gene weights by integrating:

- regulatory importance from GENIE3
- transcription factor–target gene correlation
- temporal consistency inferred from transcriptional scoring
- batch-aware penalization based on spatial dispersion in PCA space

Using this framework, BARNO supports:

- construction of regulatory weight matrices
- batch-aware regulator prioritization
- survival-associated module selection
- core regulator identification
- layered gene regulatory network (GRN) reconstruction

## Repository structure

```text
BARNO/
├── R/            # Core package functions
├── Run/          # Example analysis scripts
├── man/          # Function documentation
├── figures/      # Figures used in the manuscript / README
├── data/         # Data notes and small example files
├── DESCRIPTION
├── NAMESPACE
└── README.md
