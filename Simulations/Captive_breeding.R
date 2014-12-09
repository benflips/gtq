### Functions and model for doing PVA of quolls under evolution ###

s.dir<-"~/GitHub/gtq/" #source directory FOR MY LAPTOP

#s.dir<-"/Users/ellakelly/GitHub/gtq/" #Ben's comp source directory

#########################################
########### CAPTIVE BREEDING ############
#########################################

#NAIVE POP
init.inds<-function(n, spX, spY, init.b, init.p.var, h){
	Y<-rep(1,times=n) #X and Y from equilibrium pop
	X<-rep(1,times=n)
	S<-rbinom(n, 1, 0.5) #sex
	A<-rep(1, n) #age
	b.var<-h*init.p.var
	B<-rnorm(n, init.b, b.var^0.5) # breeding values
	P<-rnorm(n, B, (init.p.var-b.var)^0.5)
	HI<-rep(0, n) #hybrid index, all start at 0
	cbind(X, Y, S, A, B, P, HI)	
}

naive<-init.inds(n=72, spX=1, spY=1, init.b=-10.75, init.p.var=10, h=0.3)

#SMART POP
setwd("~/GitHub/gtq/Simulations/")
smart<-read.csv("ToadSmart.csv")
randomRows = function(dframe,n){
	   dframe[sample(nrow(dframe),n,replace=FALSE),]
	}
smart<-randomRows(smart, 72)

# split into breeding groups 
m.naive<-subset(naive, naive[,"S"]==1)
f.naive<-subset(naive, naive[,"S"]==0)

m.smart<-subset(smart, smart[,"S"]==1)
f.smart<-subset(smart, smart[,"S"]==0)

# CAPTIVE BREEDING repro
captive.repro<-function(m, f, n, spX, spY, init.p.var, h){
b.var<-h*init.p.var
f<- randomRows(f,n)
m<- randomRows(m,n)
B.off<-rowMeans((cbind(m[,"B"],f[,"B"])), na.rm = FALSE, dims = 1)
HI.off<-rowMeans((cbind(m[,"HI"],f[,"HI"])), na.rm = FALSE, dims = 1)
mated<-cbind(B.off, HI.off)
off<-mated[rep(1:nrow(mated), 8),]
offspring<- matrix(nrow=length(off), ncol=ncol(f), dimnames=list(NULL, colnames(f)))
offspring[,"B"]<-off[,"B.off"]
offspring[,"B"]<-rnorm(nrow(offspring), offspring[,"B"], (b.var/2)^0.5) #mid parent breeding values with variance of surviving offspring
offspring[,"S"]<-rbinom(nrow(offspring), 1, 0.5) # random sex allocation
offspring[,"A"]<-rep(0, nrow(offspring))	# age=0
offspring[,"X"]<-rep(spX, nrow(offspring))
offspring[,"Y"]<-rep(spX, nrow(offspring))
offspring[,"P"]<-rnorm(nrow(offspring), offspring[,"B"], (init.p.var-b.var)^0.5) #offspring phenotypes
offspring[,"HI"]<-off[,"HI.off"]
offspring
}
}

#SMARTxSMART
sxs<- captive.repro(m.smart, f.smart, 12, spX=5, spY=5, init.p.var=10, h=0.3)

#NAIVExNAIVE
nxn<- captive.repro(m.naive, f.naive, 12, spX=5, spY=5, init.p.var=10, h=0.3)

#SMARTxNAIVE
nxs<- captive.repro(m.naive, f.smart, 6, spX=5, spY=5, init.p.var=10, h=0.3)
sxn<- captive.repro(m.smart, f.naive, 6, spX=5, spY=5, init.p.var=10, h=0.3)
nxs<- rbind(nxs,sxn)

#########################################
### TRANSLOCATION TO VANDERLIN ISLAND ###
#########################################

#initial population
cap.pop<- rbind(nxn,sxs,nxs)

#load functions
source(paste(s.dir, "Simulations/Quoll_functions.R", sep=""))

# fitted values
load(paste(s.dir, "ABC/Kept_sims_poponlyhalfpcED.RData", sep=""))
pars<-apply(gold[,1:6], 2, mean)

#equilibrium population density 
s.dir<-"/Users/ellakelly/GitHub/gtq/Priors/" #source directory
load(paste(s.dir, "init.density", sep=""))
init.pop<-sum(test3[,3])

#run model
mother<-function(spX=10, spY=10, dem.pars, gens=50, sel.time=20, plot=FALSE, init.inds=cap.pop){
	#browser()
	alpha<-dem.pars[1]
	fsurv1<-dem.pars[2]
	fsurv2<-dem.pars[3]
	msurv<-dem.pars[4]
	beta<-dem.pars[5]
	prob.d<-dem.pars[6] #assign dem.pars
	pop<-init.inds # create a population
	n.list<-neighbours.init(spX, spY) # create a list of neighbours for each cell (individual?)
	popsize<-n
	sel<-FALSE
	for (g in 2:gens){
		#if (g==trans.time) pop<-rbind(pop, trans.pop)
		pop<-repro(popmat=pop, spX=spX, spY=spY, init.p.var=10, h=0.3, alpha=alpha, beta=beta) # females reproduce (density dep)
		pop<-mdieoff(pop, msurv) # males die off
		if (g>sel.time) sel<-TRUE # have toads arrived?
		pop<-age(pop, selection=sel, alpha=alpha, fsurv1=fsurv1, fsurv2=fsurv2) # everyone gets a little older
		pop<-disperse(pop, n.list, prob.d, spX) # and the Juveniles trudge off
		if (length(pop[,1])==0) { # did we go extinct?
			print(paste("The population went extinct at generation ", g, ":-("))
			break	
			} 
		popsize<-c(popsize, length(pop[,1])) # new population size appended to popsize vector
		if (plot==T) plotter(pop, popsize, spX, spY, sel.time, gens, fid=g)		
	}
print(mean(pop[,"B"]))
print(mean(pop[,"HI"]))
}

mother(dem.pars=pars, gens=20, sel.time=2, plot=T)