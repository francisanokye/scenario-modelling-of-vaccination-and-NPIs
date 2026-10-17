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

calibrate_plots.Rout: calibrate_plots.R calibrate.rds seroprevdata.rds fitsero.rds params.rda constrain_params.rda
	$(pipeR)

scenario_params.Rout: scenario_params.R constrain_params.rda
	$(pipeR)

base.Rout: base.R calibrate.rds scenario_params.rda
	$(pipeR)

strict_npis.Rout: strict_npis.R base.rds scenario_params.rda
	$(pipeR)

no_npis.Rout: no_npis.R base.rds scenario_params.rda
	$(pipeR)

Rc.Rout: Rc.R calibrate.rds params.rda
	$(pipeR)

impmakeR += newRc

# prefix with beta type and run eg. no_npis.newRc.Rout 
%.newRc.Rout: newRc.R %.rds Rc.rda params.rda
	$(pipeR)

combo_newRc.Rout: combo_newRc.R base.newRc.rds strict_npis.newRc.rds no_npis.newRc.rds constrain_params.rda 
	$(pipeR)

impmakeR += vacspec

%.vacspec.Rout: vacspec.R calibrate.rds %.rda constrain_params.rda
	$(pipeR)

impmakeR += betasim

# prefix with beta type and run eg. no_npis.betasim.Rout 
%.betasim.Rout: betasim.R calibrate.rds %.rds constrain_params.rda
	$(pipeR)

combo_betasim.Rout: combo_betasim.R base.betasim.rds strict_npis.betasim.rds no_npis.betasim.rds
	$(pipeR)

# vaccination scenarios
vax_actual_boosters.Rout: vax_actual_boosters.R process_vac.rda constrain_params.rda
	$(pipeR)

vax_actual_rate_can_init.Rout: vax_actual_rate_can_init.R process_vac.rda constrain_params.rda
	$(pipeR)
	
vax_early_boosters.Rout: vax_early_boosters.R process_vac.rda constrain_params.rda
	$(pipeR)

vax_faster_rate.Rout: vax_faster_rate.R process_vac.rda constrain_params.rda 
	$(pipeR)

vax_slower_rate.Rout: vax_slower_rate.R process_vac.rda constrain_params.rda vax_actual_boosters.rda
	$(pipeR)

vax_national_rate.Rout: vax_national_rate.R process_vac.rda ../data/vaccination-coverage-map.csv constrain_params.rda
	$(pipeR)

nl_init_national_rate.Rout: nl_init_national_rate.R process_vac.rda ../data/vaccination-coverage-map.csv constrain_params.rda
	$(pipeR)
	
audit_vaccination_scenarios.Rout: audit_vaccination_scenarios.R  vax_actual_boosters.rda vax_faster_rate.rda vax_slower_rate.rda \
	vax_early_boosters.rda vax_national_rate.rda
	$(pipeR)

impmakeR += vacsim

%.vacsim.Rout: vacsim.R no_npis.betasim.rda %.rda constrain_params.rda
	$(pipeR)

%.base.vacsim.Rout: vacsim.R %.vacspec.rda base.rds constrain_params.rda
	$(pipeR)

%.strict_npis.vacsim.Rout: vacsim.R %.vacspec.rda strict_npis.rds constrain_params.rda
	$(pipeR)

%.no_npis.vacsim.Rout: vacsim.R %.vacspec.rda no_npis.rds constrain_params.rda
	$(pipeR)

combo_vacsim.Rout: combo_vacsim.R vax_actual_boosters.base.vacsim.rds vax_early_boosters.base.vacsim.rds \
	vax_faster_rate.base.vacsim.rds vax_slower_rate.base.vacsim.rds vax_national_rate.base.vacsim.rds \
	vax_actual_boosters.strict_npis.vacsim.rds vax_early_boosters.strict_npis.vacsim.rds vax_faster_rate.strict_npis.vacsim.rds \
	vax_slower_rate.strict_npis.vacsim.rds vax_national_rate.strict_npis.vacsim.rds vax_actual_boosters.no_npis.vacsim.rds \
	vax_early_boosters.no_npis.vacsim.rds vax_faster_rate.no_npis.vacsim.rds vax_slower_rate.no_npis.vacsim.rds vax_national_rate.no_npis.vacsim.rds \
	nl_init_national_rate.no_npis.vacsim.rds nl_init_national_rate.base.vacsim.rds nl_init_national_rate.strict_npis.vacsim.rds \
	vax_actual_rate_can_init.no_npis.vacsim.rds vax_actual_rate_can_init.base.vacsim.rds vax_actual_rate_can_init.strict_npis.vacsim.rds
	$(pipeR)

###################### Main Text Figures ################### 

# incidence plot
Figure_4.Rout: Figure_4.R combo_vacsim.rds constrain_params.rda
	$(pipeR)
	
Figure_6.Rout: Figure_6.R combo_vacsim.rds constrain_params.rda
	$(pipeR)

Figure_7.Rout: Figure_7.R combo_vacsim.rds constrain_params.rda
	$(pipeR)

# strict npi duration
Figure_8.Rout: Figure_8.R calibrate.rds base.rds no_npis.rds scenario_params.rda vax_actual_boosters.rda vax_early_boosters.rda vax_faster_rate.rda \
	vax_slower_rate.rda vax_national_rate.rda constrain_params.rda fitsero.rds
	$(pipeR)

##### Supplementary Material Figures ################### 
Fig_S2.Rout: Fig_S2.R combo_vacsim.rds constrain_params.rda
	$(pipeR)
	
Fig_S3.Rout: Fig_S3.R combo_vacsim.rds constrain_params.rda
	$(pipeR)

Fig_S4.Rout: Fig_S4.R combo_vacsim.rds constrain_params.rda
	$(pipeR)

Fig_S5.Rout: Fig_S5.R combo_vacsim.rds constrain_params.rda
	$(pipeR)

Fig_S6.Rout: Fig_S6.R combo_vacsim.rds constrain_params.rda
	$(pipeR)
	
checks.Rout: checks.R combo_vacsim.rds constrain_params.rda
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
