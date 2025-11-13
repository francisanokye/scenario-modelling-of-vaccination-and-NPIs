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

vacdat.Rout: vacdat.R ../data/vaccination-coverage-map.csv 
	$(pipeR)

params.Rout: params.R vacdat.rda
	$(pipeR)

flows.Rout: flows.R params.rda
	$(pipeR)

timevar_spec.Rout: timevar_spec.R flows.rda vacdat.rda params.rda
	$(pipeR)

seroprevdata.Rout: seroprevdata.R params.rda ../data/raw_citf_data.csv    
	$(pipeR)

fitsero.Rout: fitsero.R ../data/raw_citf_data.csv params.rda
	$(pipeR)

calibrate.Rout: calibrate.R timevar_spec.rds seroprevdata.rds fitsero.rds params.rda
	$(pipeRcall)
	
build_multipliers.Rout: build_multipliers.R seroprevdata.rds params.rda calibrate.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

base_beta.Rout: base_beta.R build_multipliers.rds calibrate.rds
	$(pipeR)

policy_neutral.Rout: policy_neutral.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
beta_als2.Rout: beta_als2.R build_multipliers.rds calibrate.rds
	$(pipeR)

beta_als4.Rout: beta_als4.R build_multipliers.rds calibrate.rds
	$(pipeR)
	
switch_als4_als2.Rout: switch_als4_als2.R build_multipliers.rds calibrate.rds
	$(pipeR)

%.simulate.Rout: simulate.R %.rds build_multipliers.rds calibrate.rds
	$(pipeRcall)

vax_scenarios.Rout: vax_scenarios.R calibrate.rds seroprevdata.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda 
	$(pipeR)

#  Strictest ALS before relaxation at varying vaccination rates
switch_als42_counterfactual.Rout: switch_als42_counterfactual.R calibrate.rds switch_als4_als2.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

# Earlier booster rollout 
# early_boosters.Rout: early_boosters.R calibrate.rds seroprevdata.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda 
#	$(pipeR)

# Combined scenarios (varying vaccination rates and alternate over ALS)
als_counterfactual.Rout: als_counterfactual.R calibrate.rds beta_als2.rds beta_als4.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

# 100% Vaccination with Different ALS (compared to baseline)
als_vax100.Rout: als_vax100.R calibrate.rds base_beta.rds beta_als2.rds beta_als4.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

# Earlier booster rollout 
early_booster_als.Rout: early_booster_als.R calibrate.rds base_beta.rds beta_als2.rds beta_als4.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

# Early booster vaccination, with/without restrictions & varying vaccination rates 
early_vax_als.Rout: early_vax_als.R calibrate.rds base_beta.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda
	$(pipeR)
	
# init_cond_scenario.Rout: init_cond_scenario.R calibrate.rds seroprevdata.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda 
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
