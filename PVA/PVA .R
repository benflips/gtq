s.dir<-"/Users/ellakelly/Documents/GitHub/gtq/ABC/" #source directory

source(paste(s.dir, "Quoll_functions.R", sep=""))


###################################################################################
## Model Runs ##

# fitted values
load(paste(s.dir, "Kept_sims_poponlyhalfpcED.RData", sep=""))
pars<-apply(gold[,1:6], 2, mean)

#equilibrium population density 
s.dir<-"/Users/ellakelly/Documents/GitHub/gtq/Priors/" #source directory
load(paste(s.dir, "init.density", sep=""))
init.pop<-sum(test3[,3])

## PVA ##
#loop through values of init.b "n" times

#input
n <-10 #no. of iterations 	
b<- c(-8.5) #values of b
#output
e <-vector() # i loop output 
pr<- vector() # j loop output

for(j in b) {
		for (i in 1:n){
			run<-mother(dem.pars=pars, init.b=j, gens=50, sel.time=20, h=0.3)	
			e[i]<-c(run)
		}
	prob.e<- length(e[e==TRUE])/n
	pr<- c(pr, prob.e)
}
output<- cbind(pr, b)
output
write(output, "pva results")