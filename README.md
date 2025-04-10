# PRE-ACT: Cost-Effectiveness Analysis of AI-Assisted Radiotherapy

This repository contains the source code and supporting scripts for Deliverable D7.2 of the PRE-ACT project. The aim of this deliverable is to develop and validate a probabilistic decision-analytic model that assesses the cost-effectiveness of integrating AI-assisted decision-making into current clinical practice for breast cancer patients. The model simulates patient progression using a Markov state-transition framework and estimates outcomes in terms of quality-adjusted life years (QALYs) and incremental cost-effectiveness ratios (ICERs).

The study focuses on women aged ≥18 who have undergone radical surgery for invasive breast cancer (stages pT1‐4, pN0–N3, M0). It evaluates the impact of adding personalized, data-driven AI insights, designed to predict and prevent arm lymphedema, against standard clinical protocols. The analysis spans a lifetime (a 40-year simulation with monthly cycles) from a UK healthcare perspective, while the model is readily adaptable to reflect Dutch (societal perspective) and French (collective perspective) settings.

## Table of Contents

- [Background](#background)
- [Installation](#installation)
- [Usage](#usage)
- [Input Data and Parameters](#input-data-and-parameters)
- [Model Description](#model-description)
- [Validation and Analyses](#validation-and-analyses)
- [Future Work](#future-work)
- [Contact](#contact)

## Background

The PRE-ACT project addresses the need for personalized radiotherapy enhancements by integrating advanced AI-assisted risk prediction with economic evaluation. This repository supports Deliverable D7.2, which presents a probabilistic decision-analytic model that leverages clinical evidence, toxicity data, and economic parameters to compare current practice with AI-assisted decision-making. The model captures the evolution of treatment toxicity over time and assesses the incremental cost per QALY gained.

## Installation

### Prerequisites

- **R** (version 4.0.0 or later)
- R packages:
  - dampack
  - docstring
  - DiagrammeR
  - dampack
  - BCEA
  - matrixStats
  - summarytools

### Steps

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/Bram-R/PRE-ACT.git](https://github.com/Bram-R/PRE-ACT.git)
    cd PRE-ACT
    ```
2.  **Install Required R Packages:**

    In R, run:
    ```r
    required_packages <- c("docstring", "DiagrammeR", "dampack", "BCEA", "matrixStats", "summarytools")
    install.packages(required_packages)
    ```
3.  **Set the Working Directory:**
    Ensure your R working directory is set to the repository folder.

## Usage

To replicate the model-based analysis:

1.  **Main results:**
    Run the script `Main.R` to generate the main results.

2.  **Sensitivity analyses:**
    Run the script `Sensitivity analyses.R` to run probabilistic and deterministic sensitivity analyses. This step provides outputs such as the incremental cost-effectiveness plane, cost-effectiveness acceptability curves, and tornado diagrams that help identify key drivers of cost-effectiveness.

3.  **Intermediate results:**
    Run the script `Intermediate results.R` to generate the intermediate results.

## Input Data and Parameters

All input parameters (including probabilities, costs, and utility values) are are based on published literature from key clinical trials and economic studies, and they have been validated through expert consultation and TECH-VER checklist procedures.

## Model Description

The model employs a probabilistic state-transition (Markov) cohort design to simulate patient journeys across defined health states: event-free, locoregional recurrence, distant metastasis, and death. Toxicity outcomes, particularly arm lymphedema, are dynamically incorporated over a 40-year lifetime horizon with monthly cycles. This design enables the assessment of incremental cost-effectiveness ratios and provides insights into the trade-offs between increased upfront costs and long-term quality-of-life benefits.

## Validation and Analyses

The model has been rigorously validated using the TECH-VER checklist to ensure compliance with best practices. Uncertainty is examined via Monte Carlo simulations (5,000 iterations), along with deterministic sensitivity analyses and various scenario analyses. Outputs include cost-effectiveness planes, cost-effectiveness acceptability curves, and tornado diagrams that identify key drivers of model outcomes.

## Future Work

This deliverable feeds directly into subsequent work in Work Package 7. In Task 7.4, the current model will be refined by incorporating updated input data from the PRE-ACT clinical trial and expanding analyses to include perspectives from the Dutch and French healthcare systems. These enhancements will improve the predictive accuracy and generalizability of the model, further informing treatment strategy adjustments and resource allocation across diverse clinical settings.

## Contact

For further information regarding this repository or the associated deliverable, please contact:

**Bram Ramaekers**
- GitHub: @Bram-R

**Willem Witlox**
- GitHub: @willemwitlox
