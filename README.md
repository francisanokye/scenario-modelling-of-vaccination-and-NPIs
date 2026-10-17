# High vaccine coverage in Newfoundland and Labrador, Canada decreased hospitalizations and deaths due to the spread of the Omicron variant

This repository contains the R code and data used to reproduce the analyses and figures for the manuscript:

The workflow is managed with GNU Make. The `Makefile` tracks all dependencies from data processing and model calibration through counterfactual construction, simulation, and figure generation.

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

Calibration plots can be generated with:

```bash
make calibrate_plots.Rout
```

---

## NPI counterfactuals

The three NPI trajectories are created using:

| Scenario    | Script          | Output            |
| ----------- | --------------- | ----------------- |
| Actual NPIs | `base.R`        | `base.rds`        |
| Strict NPIs | `strict_npis.R` | `strict_npis.rds` |
| No NPIs     | `no_npis.R`     | `no_npis.rds`     |

Generate all three with:

```bash
make base.Rout strict_npis.Rout no_npis.Rout
```

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

## Vaccination × NPI simulations

Each vaccination scenario is converted to the vaccination specification required by the model using `vacspec.R`.

```bash
make vax_actual_boosters.vacspec.Rout
make vax_early_boosters.vacspec.Rout
make vax_faster_rate.vacspec.Rout
```

And then simulated under the three NPI conditions using `vacsim.R`.

```bash
make vax_actual_boosters.base.vacsim.Rout
make vax_actual_boosters.strict_npis.vacsim.Rout
make vax_actual_boosters.no_npis.vacsim.Rout
```

Running this below will generate all vaccination × NPI simulations required for the manuscript with:

```bash
make combo_vacsim.Rout
```

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

### Main-text figures

```bash
make Figure_4.Rout
make Figure_6.Rout
make Figure_7.Rout
make Figure_8.Rout
```

### Supplementary figures

```bash
make Fig_S2.Rout
make Fig_S3.Rout
make Fig_S4.Rout
make Fig_S5.Rout
```

## Complete reproduction

From the required raw data, the principal counterfactual simulations and manuscript figures can be reproduced with:

```bash
make combo_vacsim.Rout

make Figure_4.Rout Figure_6.Rout Figure_7.Rout Figure_8.Rout \
     Fig_S2.Rout Fig_S3.Rout Fig_S4.Rout Fig_S5.Rout Fig_S6.Rout
```