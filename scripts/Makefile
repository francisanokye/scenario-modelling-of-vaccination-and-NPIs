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

vacdat.Rout: vacdat.R ../data/vaccination-coverage-map.csv 
	$(pipeR)

params.Rout: params.R vacdat.rda 
	$(pipeR)

flows.Rout: flows.R params.rda 
	$(pipeR)

timevar_spec.Rout: timevar_spec.R flows.rda vacdat.rda derive_vax_scenarios.rda params.rda 
	$(pipeR)
	
seroprevdata.Rout: seroprevdata.R params.rda ../data/raw_citf_data.csv    
	$(pipeR)

fitsero.Rout: fitsero.R seroprevdata.rds ../data/raw_citf_data.csv params.rda
	$(pipeR)

calibrate.Rout: calibrate.R timevar_spec.rds seroprevdata.rds fitsero.rds params.rda
	$(pipeR)
	
build_multipliers.Rout: build_multipliers.R seroprevdata.rds params.rda calibrate.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda
	$(pipeR)
	
%.simulate.Rout: simulate.R %.rds build_multipliers.rds calibrate.rds
	$(pipeRcall)
	
base_beta.Rout: base_beta.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
beta_als2.Rout: beta_als2.R build_multipliers.rds calibrate.rds
	$(pipeR)

beta_als3.Rout: beta_als3.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
beta_als4.Rout: beta_als4.R build_multipliers.rds calibrate.rds
	$(pipeR)

policy_neutral.Rout: policy_neutral.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
beta_ALS4_14days_only.Rout: beta_ALS4_14days_only.R build_multipliers.rds calibrate.rds
	$(pipeR)	
	
beta_ALS4_28days_only.Rout: beta_ALS4_28days_only.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
beta_ALS4_40days_only.Rout: beta_ALS4_40days_only.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
beta_ALS4_60days_only.Rout: beta_ALS4_60days_only.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
beta_ALS4_90days_only.Rout: beta_ALS4_90days_only.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
beta_ALS4_120days_only.Rout: beta_ALS4_120days_only.R build_multipliers.rds calibrate.rds
	$(pipeR)

beta_ALS4_150days_only.Rout: beta_ALS4_150days_only.R build_multipliers.rds calibrate.rds
	$(pipeR)

beta_ALS4_180days_only.Rout: beta_ALS4_180days_only.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
control_reprod_numb.Rout: control_reprod_numb.R beta_als4.rds policy_neutral.rds seroprevdata.rds fitsero.rds params.rda calibrate.rds
	$(pipeR)

initial_can_rate.Rout: initial_can_rate.R ../data/vaccination-coverage-map.csv 
	$(pipeR)

derive_vax_scenarios.Rout: derive_vax_scenarios.R initial_can_rate.rda ../data/vaccination-coverage-map.csv 
	$(pipeR)

simulation_vax_scenarios.Rout: simulation_vax_scenarios.R calibrate.rds base_beta.rds beta_als4.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda derive_vax_scenarios.rda  ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

extension_als4.Rout: extension_als4.R calibrate.rds base_beta.rds beta_als4.rds beta_ALS4_14days_only.rds beta_ALS4_28days_only.rds beta_ALS4_40days_only.rds beta_ALS4_60days_only.rds beta_ALS4_90days_only.rds beta_ALS4_120days_only.rds beta_ALS4_150days_only.rds beta_ALS4_180days_only.rds build_multipliers.rds fitsero.rds seroprevdata.rds derive_vax_scenarios.rda params.rda  ../data/serop_avgcase_data.csv params.rda 
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
