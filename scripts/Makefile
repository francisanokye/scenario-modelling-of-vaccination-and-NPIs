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

constrain_params.Rout: constrain_params.R 
	$(pipeR)

process_vac.Rout: process_vac.R ../data/vaccination-coverage-map.csv constrain_params.rda
	$(pipeR)

params.Rout: params.R process_vac.rda constrain_params.rda
	$(pipeR)

flows.Rout: flows.R params.rda 
	$(pipeR)

timevar_spec.Rout: timevar_spec.R flows.rda params.rda process_vac.rda 
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

 
Rc.Rout: Rc.R calibrate.rds params.rda
	$(pipeR)

impmakeR += newRc

# prefix with beta type and run eg. policy_neutral.newRc.Rout 
%.newRc.Rout: newRc.R %.rds Rc.rda params.rda
	$(pipeR)

combo_newRc.Rout: combo_newRc.R base_beta.newRc.rds beta_als4.newRc.rds policy_neutral.newRc.rds constrain_params.rda 
	$(pipeR)

impmakeR += betasim

# prefix with beta type and run eg. policy_neutral.betasim.Rout 
%.betasim.Rout: betasim.R calibrate.rds %.rds constrain_params.rda
	$(pipeR)

combo_betasim.Rout: combo_betasim.R base_beta.betasim.rds beta_als4.betasim.rds policy_neutral.betasim.rds
	$(pipeR)

# vaccination scenarios
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

impmakeR += vacsim

%.vacsim.Rout: vacsim.R policy_neutral.betasim.rda %.rda constrain_params.rda
	$(pipeR)

%.base.vacsim.Rout: vacsim.R base_beta.betasim.rda %.rda constrain_params.rda
	$(pipeR)

%.als4.vacsim.Rout: vacsim.R beta_als4.betasim.rda %.rda constrain_params.rda
	$(pipeR)

%.pn.vacsim.Rout: vacsim.R policy_neutral.betasim.rda %.rda constrain_params.rda
	$(pipeR)

combo_vacsim.Rout: combo_vacsim.R vax_actual_boosters.base.vacsim.rds vax_early_boosters.base.vacsim.rds \
	vax_faster_rate.base.vacsim.rds vax_slower_rate.base.vacsim.rds vax_national_rate.base.vacsim.rds \
	vax_actual_boosters.als4.vacsim.rds vax_early_boosters.als4.vacsim.rds vax_faster_rate.als4.vacsim.rds \
	vax_slower_rate.als4.vacsim.rds vax_national_rate.als4.vacsim.rds vax_actual_boosters.pn.vacsim.rds \
	vax_early_boosters.pn.vacsim.rds vax_faster_rate.pn.vacsim.rds vax_slower_rate.pn.vacsim.rds vax_national_rate.pn.vacsim.rds
	$(pipeR)

# incidence panels
incidence_plot.Rout: incidence_plot.R params.rda
	$(pipeR)

# severity bar chart 
severity_plot.Rout: severity_plot.R params.rda
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
