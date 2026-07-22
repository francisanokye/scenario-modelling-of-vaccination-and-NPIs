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

impmakeR += betasim

# prefix with beta type and run eg. no_npis.betasim.Rout 
%.betasim.Rout: betasim.R calibrate.rds %.rds constrain_params.rda
	$(pipeR)

combo_betasim.Rout: combo_betasim.R base.betasim.rds strict_npis.betasim.rds no_npis.betasim.rds
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

%.vacsim.Rout: vacsim.R no_npis.betasim.rda %.rda constrain_params.rda
	$(pipeR)

%.base.vacsim.Rout: vacsim.R base.betasim.rda %.rda constrain_params.rda
	$(pipeR)

%.strict_npis.vacsim.Rout: vacsim.R strict_npis.betasim.rda %.rda constrain_params.rda
	$(pipeR)

%.no_npis.vacsim.Rout: vacsim.R no_npis.betasim.rda %.rda constrain_params.rda
	$(pipeR)

combo_vacsim.Rout: combo_vacsim.R vax_actual_boosters.base.vacsim.rds vax_early_boosters.base.vacsim.rds \
	vax_faster_rate.base.vacsim.rds vax_slower_rate.base.vacsim.rds vax_national_rate.base.vacsim.rds \
	vax_actual_boosters.strict_npis.vacsim.rds vax_early_boosters.strict_npis.vacsim.rds vax_faster_rate.strict_npis.vacsim.rds \
	vax_slower_rate.strict_npis.vacsim.rds vax_national_rate.strict_npis.vacsim.rds vax_actual_boosters.no_npis.vacsim.rds \
	vax_early_boosters.no_npis.vacsim.rds vax_faster_rate.no_npis.vacsim.rds vax_slower_rate.no_npis.vacsim.rds vax_national_rate.no_npis.vacsim.rds
	$(pipeR)


# incidence plot
incidence_plot.Rout: incidence_plot.R combo_vacsim.rds constrain_params.rda
	$(pipeR)
	
incidence_plot_v2.Rout: incidence_plot_v2.R combo_vacsim.rds constrain_params.rda
	$(pipeR)

# severity bar chart 
severity_plot.Rout: severity_plot.R combo_vacsim.rds constrain_params.rda
	$(pipeR)

# strict npi duration
npi_duration.Rout: npi_duration.R calibrate.rds base.rds no_npis.rds scenario_params.rda vax_actual_boosters.rda vax_early_boosters.rda vax_faster_rate.rda \
	vax_slower_rate.rda vax_national_rate.rda constrain_params.rda fitsero.rds
	$(pipeR)

# percentage of the population vaccinated at start
vax_snapshot.Rout: vax_snapshot.R combo_vacsim.rds constrain_params.rda
	$(pipeR)

diagnose_early_actual.Rout: diagnose_early_actual.R combo_vacsim.rds constrain_params.rda
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
