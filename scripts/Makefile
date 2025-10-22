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
	
build_multipliers.Rout: build_multipliers.R seroprevdata.rds params.rda calibrator_opt.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda
	$(pipeR)
	
scenario_vax_scaled.Rout: scenario_vax_scaled.R calibrator_opt.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda 
	$(pipeR)

counterfinal.Rout: counterfinal.R calibrator_opt.rds fitsero.rds seroprevdata.rds params.rda ../data/serop_avgcase_data.csv params.rda
	$(pipeR)

counterfactuals.Rout: counterfactuals.R seroprevdata.rds params.rda calibrator_vax_scale_100.rds fitsero.rds ../data/serop_avgcase_data.csv params.rda
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
