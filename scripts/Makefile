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

can_vac.Rout: can_vac.R ../data/vaccination-coverage-map.csv 
	$(pipeR)
	
params.Rout: params.R vacdat.rda
	$(pipeR)

flows.Rout: flows.R params.rda
	$(pipeR)

timevar_spec.Rout: timevar_spec.R flows.rda vacdat.rda params.rda
	$(pipeR)

seroprevdata.Rout: seroprevdata.R params.rda ../data/raw_citf_data.csv    
	$(pipeR)

fitsero.Rout: fitsero.R seroprevdata.rds ../data/raw_citf_data.csv params.rda
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
# alternate over ALS with no vaccination 	
als_novax_scenario.Rout: als_novax_scenario.R calibrate.rds base_beta.rds beta_als2.rds beta_als4.rds policy_neutral.rds build_multipliers.rds seroprevdata.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda 
	$(pipeR)

# alternate over varying vaccination rates with no ALS 	
vax_scenarios.Rout: vax_scenarios.R calibrate.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda 
	$(pipeR)

# Combined scenarios (varying vaccination rates and alternate over ALS)
comb_vax_als_scen.Rout: comb_vax_als_scen.R calibrate.rds base_beta.rds beta_als2.rds beta_als4.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

# Early booster vaccination, with/without restrictions & varying vaccination rates 
early_vax_als.Rout: early_vax_als.R calibrate.rds base_beta.rds beta_als2.rds beta_als4.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

# Combined scenarios (varying vaccination rates and alternate over ALS)
canada_vax_scen.Rout: canada_vax_scen.R calibrate.rds base_beta.rds beta_als2.rds beta_als4.rds policy_neutral.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda can_vac.rda ../data/serop_avgcase_data.csv params.rda
	$(pipeR)
	
# compute baseline vs early-booster differences
compare_bars.Rout: compare_bars.R ../output/als_vax_bar.csv ../output/early_booster_bar.csv ../output/natl_vax_rate_bar.csv 
	$(pipeR)

change_baseline_early.Rout: change_baseline_early.R ../output/comb_vax_als_bar.csv ../output/early_boosters_bar.csv 
	$(pipeR)

change_baseline_canadavax.Rout: change_baseline_canadavax.R ../output/comb_vax_als_bars.csv ../output/canada_vax_bars.csv
	$(pipeR)

#  Strictest ALS before relaxation at varying vaccination rates
#switch_als42_counterfactual.Rout: switch_als42_counterfactual.R calibrate.rds switch_als4_als2.rds build_multipliers.rds fitsero.rds seroprevdata.rds params.rda  ../data/serop_avgcase_data.csv params.rda
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
