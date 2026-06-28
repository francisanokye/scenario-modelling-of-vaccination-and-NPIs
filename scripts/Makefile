## This is [project.Makefile] …

current: target
-include target.mk
Ignore = target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

### Makestuff

Sources += Makefile
Sources += $(wildcard *.R)

## Example of accumulator
sir.Rout: sir.R
	$(pipeR)

#####################################################################################################
## Processing vaccination data
#####################################################################################################

constrain_params.Rout: constrain_params.R 
	$(pipeR)

process_vac.Rout: process_vac.R ../data/vaccination-coverage-map.csv constrain_params.rda
	$(pipeR)

params.Rout: params.R process_vac.rda constrain_params.rda
	$(pipeR)

flows.Rout: flows.R params.rda 
	$(pipeR)

timevar_spec.Rout: timevar_spec.R flows.rda params.rda process_vac.rda # vacdat.rda combine_vax_scenarios.rda
	$(pipeR)
	
seroprevdata.Rout: seroprevdata.R constrain_params.rda ../data/raw_citf_data.csv    
	$(pipeR)

fitsero.Rout: fitsero.R seroprevdata.rds ../data/raw_citf_data.csv params.rda
	$(pipeR)

calibrate.Rout: calibrate.R timevar_spec.rds seroprevdata.rds fitsero.rds params.rda
	$(pipeR)

calibrate_plots.Rout: calibrate_plots.R calibrate.rds seroprevdata.rds fitsero.rds params.rda
	$(pipeR)

scenario_params.Rout: scenario_params.R constrain_params.rda
	$(pipeR)

base_beta.Rout: base_beta.R calibrate.rds scenario_params.rda
	$(pipeR)

beta_als4.Rout: beta_als4.R base_beta.rds scenario_params.rda
	$(pipeR)

policy_neutral.Rout: policy_neutral.R base_beta.rds scenario_params.rda
	$(pipeR)

## Everything above clean?

## Need to revisit later after figuring the vaccination stuff.
control_reprod_numb.Rout: control_reprod_numb.R base_beta.rds beta_als4.rds policy_neutral.rds seroprevdata.rds fitsero.rds calibrate.rds constrain_params.rda params.rda
	$(pipeR)

vax_early_boosters.Rout: vax_early_boosters.R process_vac.rda
	$(pipeR)

vax_faster_rate.Rout: vax_faster_rate.R process_vac.rda 
	$(pipeR)

vax_slower_rate.Rout: vax_slower_rate.R process_vac.rda
	$(pipeR)

vax_national_rate.Rout: vax_national_rate.R ../data/vaccination-coverage-map.csv
	$(pipeR)

## Need to redo this. 
scenarios.Rout: scenarios.R  process_vac.rda
	$(pipeR)

### Simulate trajectories for ALS policies and vaccination scenarios	
generate_figure1.Rout: generate_figure1.R calibrate.rds base_beta.rds beta_als4.rds policy_neutral.rds vax_slower_rate.rda vax_faster_rate.rda vax_early_boosters.rda \
	vax_national_rate.rda fitsero.rds seroprevdata.rds params.rda ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

# extract_booster_info.Rout: extract_booster_info.R calibrate.rds  
#	$(pipeR)

#extension_als4.Rout: extension_als4.R calibrate.rds fitsero.rds params.rda  ../data/serop_avgcase_data.csv 
#	$(pipeR)
	
Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/00.stamp
makestuff/%.stamp:
	- $(RM) makestuff/*.stamp
	(cd makestuff && $(MAKE) pull) || git clone --depth 1 $(msrepo)/makestuff
	touch $@

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
