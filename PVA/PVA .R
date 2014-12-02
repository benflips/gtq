
##fitted values FOR MY LAPTOP
s.dir<-"/Users/ellakelly/Documents/GitHub/gtq/" #source directory

source(paste(s.dir, "Quoll_functions_outbreeding.R", sep=""))

#### FOR BEN's COMP 

#s.dir<-"/Users/ellakelly/Documents/GitHub/gtq/ABC/" #source directory

#source(paste(s.dir, "Quoll_functions_outbreeding.R", sep=""))

###################################################################################
## Model Runs ##

# fitted values
load(paste(s.dir, "ABC/Kept_sims_poponlyhalfpcED.RData", sep=""))

#load(paste(s.dir, "Kept_sims_poponlyhalfpcED.RData", sep=""))
pars<-apply(gold[,1:6], 2, mean)

#equilibrium population density 
s.dir<-"/Users/ellakelly/GitHub/gtq/Priors/" #source directory
load(paste(s.dir, "init.density", sep=""))
init.pop<-sum(test3[,3])

## PVA ##
###GET B VALUE ###
#loop through values of init.b "n" times

#input
n <-200 #no. of iterations 	
b<- c(-12,-11.5,-11,-10.5,-10) #values of b

#output
e <-vector() # i loop output 
pr<- vector() # j loop output

for(j in b) {
		for (i in 1:n){
			run<-mother(dem.pars=pars, init.b=-10.75, gens=50, sel.time=20, h=0.3, plot=F)
			e[i]<-c(run)
		}
	prob.e<- length(e[e==TRUE])/n
	pr<- c(pr, prob.e)
}
output<- cbind(pr, -12)
output
#write(output, "pva results")

####GET SMART TOADS#####

### turn on write.table in Quoll_functions
n<-2
for (i in 1:n){
  mother(dem.pars=pars, init.b=-10.75, gens=5, sel.time=20, h=0.3, plot=F)
}
