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

# calibrate.Rout: calibrate.R constrain_params.R
calibrate.Rout: calibrate.R timevar_spec.rds fitsero.rds params.rda
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

## Need to revisit later after figuring the vaccination stuff.
# Rc.Rout: Rc.R constrain_params.R 
Rc.Rout: Rc.R calibrate.rds params.rda
	$(pipeR)

impmakeR += newRc
# base_beta.newRc.Rout: newRc.R 
# beta_als4.newRc.Rout: newRc.R 
# policy_neutral.newRc.Rout: newRc.R 
%.newRc.Rout: newRc.R %.rds Rc.rda params.rda
	$(pipeR)

combo_newRc.Rout: combo_newRc.R base_beta.newRc.rds beta_als4.newRc.rds policy_neutral.newRc.rds constrain_params.rda 
	$(pipeR)

impmakeR += betasim
# base_beta.betasim.Rout: betasim.R 
# beta_als4.betasim.Rout: betasim.R 
# policy_neutral.betasim.Rout: betasim.R 
%.betasim.Rout: betasim.R calibrate.rds %.rds constrain_params.rda
	$(pipeR)

combo_betasim.Rout: combo_betasim.R base_beta.betasim.rds beta_als4.betasim.rds policy_neutral.betasim.rds
	$(pipeR)

vax_actual_boosters.Rout: vax_actual_boosters.R process_vac.rda constrain_params.rda
	$(pipeR)

vax_early_boosters.Rout: vax_early_boosters.R process_vac.rda constrain_params.rda
	$(pipeR)

vax_faster_rate.Rout: vax_faster_rate.R process_vac.rda constrain_params.rda 
	$(pipeR)

vax_slower_rate.Rout: vax_slower_rate.R process_vac.rda constrain_params.rda
	$(pipeR)

vax_national_rate.Rout: vax_national_rate.R ../data/vaccination-coverage-map.csv constrain_params.rda
	$(pipeR)

## Need to clean up above, delete commented codes

impmakeR += vacsim
# vax_actual_boosters.vacsim.Rout: vacsim.R
# vax_early_boosters.vacsim.Rout: vacsim.R
# vax_faster_rate.vacsim.Rout: vacsim.R
# vax_slower_rate.vacsim.Rout: vacsim.R
# vax_national_rate.vacsim.Rout: vacsim.R
%.vacsim.Rout: vacsim.R calibrate.rds %.rda constrain_params.rda
	$(pipeR)

combo_vacsim.Rout: combo_vacsim.R vax_actual_boosters.vacsim.rds vax_early_boosters.vacsim.rds vax_faster_rate.vacsim.rds vax_slower_rate.vacsim.rds vax_national_rate.vacsim.rds
	$(pipeR)


#vax_national_rate.Rout: vax_national_rate.R ../data/vaccination-coverage-map.csv
#	$(pipeR)

# run simulations 
sim_results.Rout: sim_results.R sim_setup.rda params.rda
	$(pipeR)

## Need to redo this. 
scenarios.Rout: scenarios.R  process_vac.rda
	$(pipeR)
	
impmakeR += sim_setup

# simulation setup
sim_setup.Rout: sim_setup.R calibrate.rds policy_neutral.rds base_beta.rds beta_als4.rds vax_slower_rate.rda vax_faster_rate.rda vax_early_boosters.rda \
	vax_national_rate.rda constrain_params.rda params.rda
	$(pipeR)

## run simulations 
#sim_results.Rout: sim_results.R sim_setup.rda params.rda
#	$(pipeR)

# figure 2 (incidence panels)
figure2.Rout: figure2.R sim_results.rda sim_setup.rda params.rda
	$(pipeR)

# figure 3 (severity bar chart) 
figure3.Rout: figure3.R sim_results.rda sim_setup.rda params.rda
	$(pipeR)
	
# vaccination snapshot table 
vax_snapshot.Rout: vax_snapshot.R sim_setup.rda params.rda
	$(pipeR)


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
