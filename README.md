# Delaying infections or reducing epidemic burden? The role of vaccination timing and non-pharmaceutical interventions during the Omicron wave

This repository contains the code and data required to reproduce all simulations and figures for the analysis of vaccination strategies and non-pharmaceutical interventions (NPIs) during the SARS-CoV-2 Omicron wave in Newfoundland and Labrador, Canada.

The analysis is implemented as a fully reproducible computational pipeline using GNU Make.

---

## Reproducibility

All results (including calibration, simulations, and figures) can be reproduced with a single command:

    make simulation_vax_scenarios.Rout

This command executes the full workflow end-to-end and regenerates all outputs from the raw inputs.

To force a full rebuild:

    make clean
    make simulation_vax_scenarios.Rout

---

## Computational workflow

The pipeline performs the following steps automatically:

1. **Model calibration**  
   Fits the epidemiological model to seroprevalence data and produces calibrated parameter objects.  

   Outputs:  
   - `calibrate.rds`  
   - `fitsero.rds`  

2. **Transmission reconstruction**  
   Estimates time-varying transmission and generates baseline and counterfactual transmission trajectories.  

   Outputs:  
   - `base_beta.rds`  
   - `beta_als4.rds`  
   - `policy_neutral.rds`  
   - `build_multipliers.rds`  

3. **Vaccination scenario construction**  
   Builds alternative vaccination trajectories, including observed rollout, earlier boosters, faster and slower uptake, and alternative initial coverage.  

   Output:  
   - `derive_vax_scenarios.rda`  

4. **Scenario simulation**  
   Simulates all combinations of vaccination and policy scenarios using the calibrated model.  

5. **Figure generation**  
   Produces all main and supplementary figures used in the analysis.  

---

## Pipeline specification

The workflow is defined by the following Make rule:

    simulation_vax_scenarios.Rout: simulation_vax_scenarios.R \
      calibrate.rds base_beta.rds beta_als4.rds policy_neutral.rds \
      build_multipliers.rds fitsero.rds seroprevdata.rds params.rda \
      derive_vax_scenarios.rda ../data/serop_avgcase_data.csv params.rda
    	$(pipeR)

All dependencies are resolved automatically. Intermediate objects are generated and reused when available.

---

## Repository structure

    .
    ├── scripts/              # Model specification, simulation, and plotting
    ├── data/                 # Input datasets
    ├── figures/              # Generated figures
    ├── *.rds / *.rda         # Intermediate model objects
    ├── Makefile              # Reproducible pipeline

---

## Software requirements

- R (≥ 4.2)  
- GNU Make  

---

## R package dependencies

Install required packages in R:

    install.packages(c(
      "tidyverse", "ggplot2", "dplyr", "tidyr", "purrr",
      "readr", "stringr", "lubridate", "scales",
      "patchwork", "cowplot", "zoo", "splines",
      "broom.mixed", "ggthemes", "readxl", "tibble"
    ))

    install.packages("remotes")
    remotes::install_github("canmod/macpan2")
    remotes::install_github("canmod/shellpipes")

---

## Inputs

The pipeline requires:

- seroprevalence data (`seroprevdata.rds`)  
- case data (`../data/serop_avgcase_data.csv`)  
- model parameter file (`params.rda`)  

All required inputs are included in the repository.

---

## Outputs

Running the pipeline produces:

### Figures (saved in `figures/`)

- vaccination coverage trajectories  
- daily incidence trajectories  
- cumulative infections  
- infections during vs post-policy periods  
- peak incidence comparisons  
- threshold-crossing analyses  
- reproduction number trajectories  
