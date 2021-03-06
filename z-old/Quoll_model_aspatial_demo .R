## Load functions ##
source(file="~/evo-dispersal/Quoll_model/Quoll_model_functions.R", echo=F)
###################################################################################
## Model Runs ##

#set the working directory
setwd("~/Documents/PhD/Figures")

mother<-function(n=100, fsurv1=0.23, fsurv2=0.02, msurv=0.042, init.b=0, init.p.var=10, h=0.3, gens=50, fec=6.725, sel.time=20, plot=FALSE){
	pop<-init.inds(n, init.b, init.p.var, h) # create a population
	popsize<-n
	sel<-FALSE
	for (g in 2:gens){
		pop<-repro(popmat=pop, fec=fec, init.p.var=init.p.var, h=h) # females reproduce (density dep)
		pop<-mdieoff(pop, msurv) # males die off
		if (g>sel.time) sel<-TRUE # have toads arrived?
		pop<-age(pop, selection=sel, alpha=alpha, fsurv1=fsurv1, fsurv2=fsurv2, jsurv=jsurv, msurv=msurv) # everyone gets a little older
		if (length(pop[,1])==0) { # did we go extinct?
			extinct<-TRUE
			print(paste("The population went extinct at generation ", g, ":-("))
			break	
		} 
		else{
		extinct<-FALSE
		}
		popsize<-c(popsize, length(pop[,1])) # new population size appended to popsize vector
		if (plot==T) plotter(pop, popsize, sel.time, gens, fid=g)
	}
	list(pop, popsize)
	print(extinct)	
}

plotter<-function(popmatrix, popsize, sel.time, gens, fid){
	fname<-paste("simgen", fid, ".png", sep="")
	#png(fname, width=12, height=12, units="cm", res=300)	#p.summ<-tapply(popmatrix[,"P"], list(factor(popmatrix[,"X"], levels=1:spX), factor(popmatrix[,"Y"], levels=1:spY)), mean)
	par(mfrow=c(2, 2))
	hist(popmatrix[,"P"], xlim=c(-30, 30), main="Phenotype")
	plot(popsize, type="l", ylab="Population Size", xlab="Generations", bty="l", xlim=c(0, gens))
	arrows(sel.time,0,sel.time,100000, length=0, col="red")
	hist(popmatrix[,"B"], xlim=c(-30, 30), main="Genotype")
	#dev.off()
}



mother(h=0.3, init.b=0, plot=T, gens=20)

#loop? 
#values<-seq(from = 0, to = 1, length.out = 10)
#rep(values, times=100)
#for(i in values){
 # mother(h=0, plot=T, init.b= i, gens=50)
#}
