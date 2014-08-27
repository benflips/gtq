s.dir<-"~/GitHub/gtq/ABC/" #source directory

source(paste(s.dir, "Quoll_functions.R", sep=""))


###################################################################################
## Model Runs ##

# fitted values
load(paste(s.dir, "Kept_sims_poponlyhalfpcED.RData", sep=""))
pars<-apply(gold[,1:6], 2, mean)

s.dir<-"~/GitHub/gtq/Priors/" #source directory
load(paste(s.dir, "init.density", sep=""))
init.pop<-sum(test3[,3])

test<-mother(dem.pars=pars, gens=10, sel.time=101, h=0.3, plot=T)