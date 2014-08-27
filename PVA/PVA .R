s.dir<-"~/GitHub/gtq/ABC/" #source directory

source(paste(s.dir, "Quoll_functions.R", sep=""))


###################################################################################
## Model Runs ##

# fitted values
load(paste(s.dir, "Kept_sims_poponlyhalfpcED.RData", sep=""))
pars<-apply(gold[,1:6], 2, mean)


system.time(test<-mother(dem.pars=pars, gens=100, sel.time=101))