# Delaying infections or reducing epidemic burden? The role of vaccination timing and non-pharmaceutical interventions during the Omicron wave

This repository contains the code and data required to reproduce all simulations and figures for the analysis of vaccination strategies and non-pharmaceutical interventions (NPIs) during the SARS-CoV-2 Omicron wave in Newfoundland and Labrador, Canada.

The analysis is implemented as a fully reproducible computational pipeline using GNU Make.

---

## Reproducibility

All results (including calibration, simulations, and figures) can be reproduced with a single command:

```bash
make simulation_vax_scenarios.Rout
