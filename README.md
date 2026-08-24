# High vaccine coverage in Newfoundland and Labrador, Canada decreased hospitalizations and deaths due to the spread of the Omicron variant

This repository contains the R code used to reproduce the analyses and figures for the manuscript:

**High vaccine coverage in Newfoundland and Labrador, Canada decreased hospitalizations and deaths due to the spread of the Omicron variant**

We used a vaccination-stratified epidemiological model calibrated to serological data to evaluate COVID-19 vaccination and non-pharmaceutical intervention (NPI) counterfactuals in Newfoundland and Labrador (NL), Canada.

## Analysis overview

The analysis compares the vaccination campaign implemented in NL with counterfactual scenarios in which booster doses were:

* administered earlier;
* administered at a faster rate;
* administered at a slower rate; or
* administered at the average Canadian rate.

Additional scenarios separate the effects of the vaccination rollout rate from differences in initial vaccination coverage between NL and Canada.

Each vaccination scenario is simulated under three NPI conditions:

* **Actual NPIs**
* **Strict NPIs**
* **No NPIs**

The workflow is managed with GNU Make. The `Makefile` tracks all dependencies from data processing and model calibration through counterfactual construction, simulation, and figure generation.

```text
Raw data
   |
   +-- Model and vaccination inputs
   |
   +-- Serological data
   |
   +-- Model calibration
           |
           +-- NPI counterfactuals
           |
           +-- Vaccination counterfactuals
                    |
                    +-- Vaccination × NPI simulations
                              |
                              +-- combo_vacsim.rds
                                        |
                                        +-- Figures
```

## Requirements

The analysis requires:

* R
* GNU Make
* `shellpipes`
* R packages loaded by the individual scripts

The Makefile uses `makestuff` from:

```text
https://github.com/dushoff/makestuff
```

Run all commands from the directory containing the `Makefile`.

The principal input data files are:

```text
../data/vaccination-coverage-map.csv
../data/raw_citf_data.csv
```

Intermediate `.rda`, `.rds`, and `.Rout` files are generated automatically.

---

## Model calibration

The calibrated model can be reproduced with:

```bash
make calibrate.Rout
```

Make automatically generates the required upstream objects from:

```text
constrain_params.R
process_vac.R
params.R
flows.R
timevar_spec.R
seroprevdata.R
fitsero.R
calibrate.R
```

The principal calibration output is:

```text
calibrate.rds
```

Calibration plots can be generated with:

```bash
make calibrate_plots.Rout
```

---

## NPI counterfactuals

Scenario parameters are generated with:

```bash
make scenario_params.Rout
```

The three NPI trajectories are then created using:

| Scenario    | Script          | Output            |
| ----------- | --------------- | ----------------- |
| Actual NPIs | `base.R`        | `base.rds`        |
| Strict NPIs | `strict_npis.R` | `strict_npis.rds` |
| No NPIs     | `no_npis.R`     | `no_npis.rds`     |

Generate all three with:

```bash
make base.Rout strict_npis.Rout no_npis.Rout
```

The exact transmission-rate assumptions defining the counterfactuals are implemented in `base.R`, `strict_npis.R`, and `no_npis.R`.

---

## Vaccination counterfactuals

The vaccination scenarios are defined in separate scripts:

| Scenario                                        | Script                       |
| ----------------------------------------------- | ---------------------------- |
| Actual NL                                       | `vax_actual_boosters.R`      |
| Earlier                                         | `vax_early_boosters.R`       |
| Faster rate                                     | `vax_faster_rate.R`          |
| Slower rate                                     | `vax_slower_rate.R`          |
| Canadian rate                                   | `vax_national_rate.R`        |
| Canadian rate with NL initial conditions        | `nl_init_national_rate.R`    |
| Actual NL rate with Canadian initial conditions | `vax_actual_rate_can_init.R` |

Generate all vaccination scenarios with:

```bash
make vax_actual_boosters.Rout \
     vax_early_boosters.Rout \
     vax_faster_rate.Rout \
     vax_slower_rate.Rout \
     vax_national_rate.Rout \
     nl_init_national_rate.Rout \
     vax_actual_rate_can_init.Rout
```

The exact timing shifts, rollout-rate modifications, coverage limits, and initial conditions are defined in the corresponding scenario scripts.

The principal vaccination scenarios can be checked with:

```bash
make audit_vaccination_scenarios.Rout
```

---

## Vaccination × NPI simulations

Each vaccination scenario is converted to the vaccination specification required by the model using `vacspec.R`.

For example:

```bash
make vax_actual_boosters.vacspec.Rout
make vax_early_boosters.vacspec.Rout
make vax_faster_rate.vacspec.Rout
```

Each vaccination scenario is then simulated under the three NPI conditions using `vacsim.R`.

The Makefile uses the following naming convention:

```text
*.base.vacsim.Rout          -> Actual NPIs
*.strict_npis.vacsim.Rout   -> Strict NPIs
*.no_npis.vacsim.Rout       -> No NPIs
```

For example, the actual NL vaccination campaign under each NPI condition is generated with:

```bash
make vax_actual_boosters.base.vacsim.Rout
make vax_actual_boosters.strict_npis.vacsim.Rout
make vax_actual_boosters.no_npis.vacsim.Rout
```

All vaccination × NPI simulations required for the manuscript are combined with:

```bash
make combo_vacsim.Rout
```

The resulting file is:

```text
combo_vacsim.rds
```

This is the principal simulation dataset used by the main and supplementary figure scripts.

Because Make resolves dependencies automatically, running:

```bash
make combo_vacsim.Rout
```

will generate any missing NPI scenarios, vaccination scenarios, vaccination specifications, and simulations required to create `combo_vacsim.rds`.

---

## Reproduction number

The reproduction number for the calibrated model is calculated with:

```bash
make Rc.Rout
```

Counterfactual reproduction-number trajectories for the actual, strict, and no-NPI scenarios are combined with:

```bash
make combo_newRc.Rout
```

---

## Reproducing the figures

### Main-text figures

```bash
make Figure_2.Rout
make Figure_3.Rout
make Figure_4.Rout
make Figure_5.Rout
```

Figures 2–4 use:

```text
combo_vacsim.rds
```

as their principal simulation input.

Figure 5 evaluates the duration of strict NPIs under alternative vaccination scenarios and uses the calibrated model and vaccination-scenario objects directly.

### Supplementary figures

```bash
make Fig_S2.Rout
make Fig_S3.Rout
make Fig_S4.Rout
make Fig_S5.Rout
```

These figures use `combo_vacsim.rds` as their principal simulation input.

All manuscript figures can be generated with:

```bash
make Figure_2.Rout Figure_3.Rout Figure_4.Rout Figure_5.Rout \
     Fig_S2.Rout Fig_S3.Rout Fig_S4.Rout Fig_S5.Rout
```

Make will automatically rebuild any missing or outdated upstream dependencies.

---

## Complete reproduction

From the required raw data, the principal counterfactual simulations and manuscript figures can be reproduced with:

```bash
make combo_vacsim.Rout

make Figure_2.Rout Figure_3.Rout Figure_4.Rout Figure_5.Rout \
     Fig_S2.Rout Fig_S3.Rout Fig_S4.Rout Fig_S5.Rout
```

To reproduce the NPI and vaccination scenarios explicitly before simulation:

```bash
# NPI scenarios
make base.Rout strict_npis.Rout no_npis.Rout

# Vaccination scenarios
make vax_actual_boosters.Rout \
     vax_early_boosters.Rout \
     vax_faster_rate.Rout \
     vax_slower_rate.Rout \
     vax_national_rate.Rout \
     nl_init_national_rate.Rout \
     vax_actual_rate_can_init.Rout

# Combined counterfactual simulations
make combo_vacsim.Rout

# Figures
make Figure_2.Rout Figure_3.Rout Figure_4.Rout Figure_5.Rout \
     Fig_S2.Rout Fig_S3.Rout Fig_S4.Rout Fig_S5.Rout
```

