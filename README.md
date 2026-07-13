# The role of vaccination timing and non-pharmaceutical interventions during the Omicron wave

Reproducible analysis code for a **vaccination-stratified SEAIR** compartmental model of
COVID-19 transmission, dynamics, and severity under alternative vaccination and
transmission-control scenarios, applied to the SARS-CoV-2 Omicron wave in Newfoundland and
Labrador, Canada. The model is implemented in [`macpan2`](https://github.com/canmod/macpan2)
and the full analysis is orchestrated as a dependency-tracked `make` pipeline using the
[`makestuff`](https://github.com/dushoff/makestuff) framework.

> **Associated publication:** _Add the manuscript / preprint title, authors, and DOI here._
> **Archived release:** _Add the Zenodo/OSF DOI for the exact code version used in the paper._
> If you use this code, please cite both the publication and the archived release (see §13).

---

## Contents

1. [Model overview](#1-model-overview)
2. [Repository structure](#2-repository-structure)
3. [System requirements](#3-system-requirements)
4. [Installation and environment capture](#4-installation-and-environment-capture)
5. [Data availability](#5-data-availability)
6. [Reproducibility and determinism](#6-reproducibility-and-determinism)
7. [How the pipeline works](#7-how-the-pipeline-works)
8. [Reproducing the analysis, stage by stage](#8-reproducing-the-analysis-stage-by-stage)
9. [Pattern-rule (`%`) targets](#9-pattern-rule--targets)
10. [One-command reproduction and verification](#10-one-command-reproduction-and-verification)
11. [Runtime and hardware expectations](#11-runtime-and-hardware-expectations)
12. [Troubleshooting](#12-troubleshooting)
13. [Citation](#13-citation)
14. [License](#14-license)
15. [Funding and acknowledgements](#15-funding-and-acknowledgements)

---

## Model overview

The model tracks susceptible–exposed–infectious–recovered dynamics with the infectious
class split into asymptomatic and symptomatic compartments (**SEAIR**):

| Compartment | Meaning |
|-------------|---------|
| `S` | Susceptible |
| `E` | Exposed (latent, not yet infectious) |
| `A` | Asymptomatic infectious |
| `I` | Symptomatic infectious |
| `R` | Recovered / immune |

Each compartment is resolved across **three vaccination strata**, allowing dose-specific
susceptibility, severity, and force of infection:

| Stratum | Label | Interpretation |
|---------|-------|----------------|
| 1 | `S1` | Unvaccinated / partially vaccinated |
| 2 | `V2` | Double-vaccinated |
| 3 | `V3` | Boosted |

The force of infection is driven by symptomatic cases (`I`) plus down-weighted asymptomatic
cases (`A`, scaled by a relative-infectiousness parameter `zeta`). Vaccination moves
individuals between strata (`S1 → V2 → V3`), and the model accumulates incidence,
hospitalization, and death outputs separately for the ≤2-dose and booster groups. The
complete compartment and flow specification is in [`flows.R`](flows.R).

---

## Repository structure

```
.
├── Makefile                 # Pipeline definition (all build targets and dependencies)
├── flows.R                  # SEAIR compartment/flow model definition (macpan2)
├── constrain_params.R       # Parameter constraints
├── process_vac.R            # Vaccination-coverage data processing
├── params.R                 # Assembled model parameters
├── timevar_spec.R           # Time-varying parameter specification
├── seroprevdata.R           # Seroprevalence data preparation
├── fitsero.R                # Fit to seroprevalence data
├── calibrate.R              # Model calibration
├── calibrate_plots.R        # Calibration diagnostic plots
├── scenario_params.R        # Scenario parameter sets
├── base_beta.R              # Baseline transmission (beta) scenario
├── beta_als4.R              # Alternative transmission scenario
├── policy_neutral.R         # Policy-neutral transmission scenario
├── Rc.R, newRc.R            # Effective reproduction number calculations
├── betasim.R                # Transmission-scenario simulations
├── vax_*.R                  # Vaccination scenarios (see §8)
├── vacsim.R                 # Vaccination-scenario simulations
├── combo_*.R                # Roll-ups combining scenario outputs
├── incidence_plot.R         # Incidence figure
├── severity_plot.R          # Severity figure
├── makestuff/               # Build framework (cloned on first run; not version-controlled)
└── ../data/                 # Raw input data (sibling directory — see §5)
```

For each script `x.R`, running the pipeline produces a run log `x.Rout` and one or more
serialized outputs (`x.rds` / `x.rda`) that downstream steps consume.

---

## System requirements

**Software.** The pipeline has been tested on the versions listed below. _Replace the
placeholders with the exact versions used for the published results_ — pinning versions is
essential for reproducibility.

| Component | Tested version | Purpose |
|-----------|----------------|---------|
| Operating system | _e.g. Ubuntu 22.04 / macOS 14_ | — |
| **R** | _e.g. 4.3.2_ (≥ 4.2) | Runs every analysis script |
| **`macpan2`** | _commit / version_ | Compartmental modeling engine (`flows.R`) |
| **`shellpipes`** | _version_ | R-side helper for the makestuff pipeline |
| **GNU make** | _e.g. 4.3_ | Executes the pipeline and dependency graph |
| **git** | _e.g. 2.40_ | Fetches the `makestuff` build framework |
| **bash** | _e.g. 5.1_ | Shell used by the build rules |

**R package dependencies.** Install the analysis packages in R (verify against the
`library(...)` calls at the top of each `*.R` file — those are authoritative):

```r
install.packages(c(
  "tidyverse", "ggplot2", "dplyr", "tidyr", "purrr",
  "readr", "stringr", "lubridate", "scales",
  "patchwork", "cowplot", "zoo", "splines",
  "broom.mixed", "ggthemes", "readxl", "tibble"
))
install.packages("remotes")
remotes::install_github("canmod/macpan2")
remotes::install_github("canmod/shellpipes")
```

**Hardware.** _State the machine the published results were produced on (CPU, cores, RAM)
and any minimum requirements_, e.g. "Runs on a standard laptop with 8 GB RAM; no GPU
required." See §11 for approximate runtimes.

---

## Installation and environment capture

Clone the repository, then fetch the `makestuff` build framework, which is not committed to
the repository and is cloned automatically on first run:

```bash
git clone https://github.com/francisanokye/scenario-modelling-of-vaccination-and-NPIs.git>
cd scenario-modelling-of-vaccination-and-NPIs/scripts/

# Fetch the makestuff build framework
make makestuff/00.stamp
```

If that step does not trigger, clone `makestuff` manually into the project directory:

```bash
git clone --depth 1 https://github.com/dushoff/makestuff
```

You should now have a `makestuff/` folder beside the `Makefile`. Errors mentioning an
undefined `$(pipeR)` or missing includes mean `makestuff` is absent — repeat this step.

---

## Data availability

The pipeline reads raw inputs from a **sibling `../data/` directory** (one level above the
project folder):

| File | Consumed by | Description |
|------|-------------|-------------|
| `../data/vaccination-coverage-map.csv` | `process_vac`, `vax_national_rate` | _PHAC - https://health-infobase.canada.ca/covid-19/vaccination-coverage/_ |
| `../data/raw_citf_data.csv` | `seroprevdata`, `fitsero` | _COVID-19 Immunity Task Force (CITF) - https://www.covid19immunitytaskforce.ca/mandate-strategy/_ |


```bash
sha256sum ../data/vaccination-coverage-map.csv ../data/raw_citf_data.csv
```

---

## Reproducibility - How the pipeline works

Every analysis step is a `make` **target** named after its script, using the `.Rout`
extension. Requesting a target runs the corresponding `.R` script and produces its log and
outputs:

```bash
make process_vac.Rout      # runs process_vac.R
make calibrate.Rout        # runs calibrate.R and all of its upstream dependencies
```

---

### Parameters and vaccination processing

```bash
make constrain_params.Rout   # base parameter constraints
make process_vac.Rout        # process vaccination coverage (needs constrain_params + coverage CSV)
make params.Rout             # assemble model parameters (needs process_vac, constrain_params)
make flows.Rout              # SEAIR flow definition (needs params)
make timevar_spec.Rout       # time-varying parameter spec (needs flows, params, process_vac)
```

### Serology and calibration

```bash
make seroprevdata.Rout       # prepare seroprevalence data (needs constrain_params + CITF CSV)
make fitsero.Rout            # fit to seroprevalence (needs seroprevdata, CITF CSV, params)
make calibrate.Rout          # calibrate the model (needs timevar_spec, fitsero, params)
make calibrate_plots.Rout    # calibration diagnostics (needs calibrate, seroprevdata, fitsero, params)
```

### Transmission (beta) scenarios

```bash
make scenario_params.Rout    # scenario parameter sets (needs constrain_params)
make base_beta.Rout          # baseline transmission (needs calibrate, scenario_params)
make beta_als4.Rout          # alternative transmission (needs base_beta, scenario_params)
make policy_neutral.Rout     # policy-neutral transmission (needs base_beta, scenario_params)
```

### Reproduction number and transmission simulations

```bash
make Rc.Rout                 # effective reproduction number (needs calibrate, params)
make combo_newRc.Rout        # per-scenario Rc recomputation, combined
make combo_betasim.Rout      # per-scenario transmission simulations, combined
```

### Vaccination scenarios

```bash
make vax_actual_boosters.Rout   # observed booster rollout
make vax_early_boosters.Rout    # earlier booster rollout
make vax_faster_rate.Rout       # faster vaccination rate
make vax_slower_rate.Rout       # slower vaccination rate
make vax_national_rate.Rout     # national-rate rollout (needs coverage CSV)
```

Then the combined vaccination simulations:

```bash
make combo_vacsim.Rout
```

### Figures

```bash
make incidence_plot.Rout     # incidence panels
make severity_plot.Rout      # severity bar chart
```

---

## Pattern-rule (`%`) targets

Several rules use `%` as a wildcard **stem**. They are never run literally: substitute a
concrete name for `%`, and `make` supplies the corresponding dependencies. The stem must
correspond to an output (`.rds` / `.rda`) produced by an earlier step, because that file is
a dependency of the rule.

### `%.newRc.Rout` — reproduction number per transmission scenario

Rule: `%.newRc.Rout: newRc.R %.rds Rc.rda params.rda`. The stem is a transmission scenario.

```bash
make base_beta.newRc.Rout
make beta_als4.newRc.Rout
make policy_neutral.newRc.Rout
make combo_newRc.Rout        # builds the three above, then combines them
```

### `%.betasim.Rout` — transmission simulation per scenario

Rule: `%.betasim.Rout: betasim.R calibrate.rds %.rds constrain_params.rda`. Same three
scenario stems.

```bash
make base_beta.betasim.Rout
make beta_als4.betasim.Rout
make policy_neutral.betasim.Rout
make combo_betasim.Rout      # combines the three
```

### `%.vacsim.Rout` — vaccination simulation

Each vaccination simulation pairs a **vaccination scenario** (the stem `%`) with a
**transmission scenario** (selected by the middle keyword). The five vaccination-scenario
stems are `vax_actual_boosters`, `vax_early_boosters`, `vax_faster_rate`,
`vax_slower_rate`, and `vax_national_rate`. The middle keyword selects the transmission
simulation:

| Keyword | Rule | Pairs with |
|---------|------|------------|
| `.base` | `%.base.vacsim.Rout` | `base_beta.betasim.rda` |
| `.als4` | `%.als4.vacsim.Rout` | `beta_als4.betasim.rda` |
| `.pn`   | `%.pn.vacsim.Rout`   | `policy_neutral.betasim.rda` |

Each target therefore has the form `<vaccination_scenario>.<transmission>.vacsim.Rout`:

```bash
make vax_actual_boosters.base.vacsim.Rout   # actual boosters × baseline transmission
make vax_faster_rate.als4.vacsim.Rout       # faster rate     × alternative transmission
make vax_national_rate.pn.vacsim.Rout       # national rate   × policy-neutral transmission
```

The roll-up `combo_vacsim.Rout` requires **all 15 combinations** (5 vaccination scenarios ×
3 transmission scenarios). You do not need to enumerate them — request the combined target
and `make` builds any that are missing:

```bash
make combo_vacsim.Rout
```

---

## One-command reproduction and verification

For a guaranteed-clean rebuild, delete generated artifacts first and rerun the commands
above:

```bash
rm -f *.Rout *.rds *.rda *.args *.Rlog *.pdf    # removes build outputs only; scripts and data are untouched
```

Regenerate the full set of results:

```bash
make calibrate_plots.Rout
make combo_newRc.Rout
make combo_betasim.Rout
make combo_vacsim.Rout
make incidence_plot.Rout
make severity_plot.Rout
```

---

## Troubleshooting

- **A step fails.** Open the corresponding `*.Rout` log; the R error appears at the bottom.
- **Force a rebuild.** Use `make -B <target>.Rout` to rerun regardless of timestamps, or
  `touch <step>.R` to mark a step stale.
- **Inputs changed.** Editing a file in `../data/` causes the steps that read it, and
  everything downstream, to rebuild on the next `make`.
- **`$(pipeR)` undefined / missing includes.** The `makestuff` framework is not present;
  repeat the setup in §4.
- **`vim_session`** is a developer convenience target (`bash -cl "vmt"`) and is not part of
  the analysis.

---

## Citation

The citation for both the software and the associated publication:

```bibtex
@article{YOURKEY,
  title   = {The role of vaccinationtiming and non-pharmaceutical interventions during the Omicron wave},
  author  = {Author, A. and Author, B.},
  journal = {JOURNAL},
  year    = {YEAR},
  doi     = {DOI}
}

@software{YOURKEY_software,
  title   = {Vaccination-Stratified SEAIR Model of COVID-19},
  author  = {Author, A. and Author, B.},
  year    = {YEAR},
  version = {vX.Y.Z},
  doi     = {ZENODO_DOI}
}
```

---

