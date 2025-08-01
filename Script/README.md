# **Genome-scale metabolic modeling guided _Escherichia coli_ engineering for _de novo_ biosynthesis of chrysanthemic acid**

This project involves the analysis, modeling, and simulation of metabolic processes, focusing on data extraction, flux variability analysis, and model reconstruction. The project uses MATLAB scripts to manipulate metabolic models and evaluate performance under different conditions.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
    - [Main Scripts](#main-scripts)
    - [Additional Functions](#additional-functions)
- [License](#license)

## Overview

### This project aims to:
- Reconstruct and update metabolic models.
- Perform flux variability analysis (FVA) and sampling.
- Simulating the reduction of *ispA* expression to enhance the yield of downstream products
- Simulating the enhancement in  *ALDH1*, *ADH2*, or both genes expression and their effects.

### Key Scripts

- **Metabolic Model Reconstruction**:
	-Rebuild and refine metabolic models using `Reconstruction.m` and `UpdateModel.m`.
- **Flux Analysis**: 
	-Perform FVA and Jaccard index calculations using `FVAJaccardIndex.m` and `SamplingAndFVA.m`.
- **Memote Testing**:
	-Evaluate metabolic models for consistency and performance using `Memote_BL21DE3_CS3.m`, `Memote_BL21DE3_CS5.m`, and `Memote_BL21DE3_WT.m`.
- **Simulated gene knockdown**: 
	-The reduction of *ispA* expression was simulated via `SimulateReduceGeneExpression.m` to quantify its impact on downstream biosynthesis efficiency.
- **Simulated gene enhance**: 
    -The enhancement of *ALDH1*, *ADH2*, or both genes was simulated via `SimulateEnhanceGeneExpression.m` to quantify their impact on downstream product yield. 

### Key Functionalities

- **Data Extraction**: Extract protein-related information using `ExtractProteinInformation.m`.
- **Update the chemical formulas**: Modify and update the chemical formulas of metabolites in the model  using `ReviseMetsFormula.m`.
## Prerequisites

- MATLAB (version 2021a or higher)
- COBRA Toolbox 3.0
- solver : 'ibm_cplex'  

## Installation
1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/KeL1n1995/Chrysanthemic-GEMs.git
   ```

2. Set up required MATLAB toolboxes and dependencies
    
3. Add the necessary paths in MATLAB:
    `addpath(genpath('yourpath'));`
    
## Usage

### Main Scripts
- **`Main.m`**: The main entry point for running the core analyses.
- **`Reconstruction.m`**: Script for reconstructing the metabolic model.
- **`UpdateModel.m`**: Updates the model based on new constraints or conditions.
- **`SamplingIndividualData.m`**: Performs sampling for individual metabolic data.
- **`Memote_BL21DE3_CS3.m`**, **`Memote_BL21DE3_CS5.m`**, **`Memote_BL21DE3_WT.m`**: Perform `memote` tests on different versions of the BL21DE3 strain model.
- **`SimulateReduceGeneExpression.m`**: Evaluates the reduction in *ispA* gene expression and its effects.
- **`SimulateEnhanceGeneExpression.m`**:Evaluates the enhancement in  *ALDH1*, *ADH2*, or both genes expression and their effects.
- **`FVAJaccardIndex.m`**: Calculates the Jaccard Index based on FVA results.
- **`SamplingAndFVA.m`**: Combines flux variability analysis with sampling techniques.
- **`PredictedEnzymeReactionFlux.m`**:Predicted enzyme reaction flux as a function of substrate concentration and gene expression.
- **`VmaxRateRatiosCalculate`**:The catalytic rate ratios of different enzymes interacting with the DMAPP.

### Additional Functions
- **`ExtractProteinInformation.m`**: Extracts key protein information from the model.
- **`ReviseMetsFormula.m`** :Modify and update the chemical formulas of metabolites in the model.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.
