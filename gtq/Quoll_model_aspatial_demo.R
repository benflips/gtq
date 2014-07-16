## Load functions ##
source(file="~/evo-dispersal/Quoll_model/Quoll_model_functions.R", echo=F)
###################################################################################
## Model Runs ##

#set the working directory
setwd("~/Documents/PhD/Figures")

mother<-function(n=2000, spX=10, spY=10, alpha=-1, fsurv1=0.4, fsurv2=0.1, msurv=0.1, init.b=0, init.p.var=10, h=0.3, gens=10, K=20, fec=4, beta=100, prob.d=0.1, sel.time=10, plot=FALSE){
	browser()
	pop<-init.inds(n, spX, spY, init.b, init.p.var, h) # create a population
	n.list<-neighbours.init(spX, spY) # create a list of neighbours for each cell (individual?)
	popsize<-n
	sel<-FALSE
	for (g in 2:gens){
		pop<-repro(popmat=pop, spX=spX, spY=spY, fec=fec, init.p.var=init.p.var, h=h) # females reproduce (density dep)
		pop<-mdieoff(pop, msurv) # males die off
		if (g>sel.time) sel<-TRUE # have toads arrived?
		pop<-age(pop, selection=sel, alpha=alpha, fsurv1=fsurv1, fsurv2=fsurv2, msurv=msurv,spX=spX, spY=spY, beta=beta, K=K) # everyone gets a little older
		pop<-disperse(pop, n.list, prob.d, spX) # and the Juveniles trudge off
		if (length(pop[,1])==0) { # did we go extinct?
			popsize<-c(popsize, 0)
			print(paste("The population went extinct at generation ", g, ":-("))
			break	
		}
		popsize<-c(popsize, length(pop[,1])) # new population size appended to popsize vector
		if (plot==T) plotter(pop, popsize, spX, spY, sel.time, gens, fid=g)
	}
	list(pop, popsize)	
}

plotter<-function(popmatrix, popsize, spX, spY, sel.time, gens, fid){
	fname<-paste("simgen", fid, ".png", sep="")
	#png(fname, width=12, height=12, units="cm", res=300)
	p.summ<-table(factor(popmatrix[,"X"], levels=1:spX), factor(popmatrix[,"Y"], levels=1:spY))
	#p.summ<-tapply(popmatrix[,"P"], list(factor(popmatrix[,"X"], levels=1:spX), factor(popmatrix[,"Y"], levels=1:spY)), mean)
	par(mfrow=c(2, 2))
	image(p.summ, col=rev(heat.colors(10)), breaks=seq(0, 1000, length=11), main="Population density")
	hist(popmatrix[,"P"], xlim=c(-30, 30), main="Phenotype")
	plot(popsize, type="l", ylab="Population Size", xlab="Generations", bty="l", xlim=c(0, gens))
	arrows(sel.time,0,sel.time,100000, length=0, col="red")
	hist(popmatrix[,"B"], xlim=c(-30, 30), main="Genotype")
	#dev.off()
}



mother(h=0.3, init.b=-5, plot=T, gens=20)
	

