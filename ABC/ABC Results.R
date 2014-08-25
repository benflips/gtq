d.dir<-"~/gtq/"


load(paste(d.dir, "ABC/Kept_sims_poponlyhalfpcED.RData", sep=""))

load(paste(d.dir, "Priors/SurvivalPriors.rdata", sep=""))

#generalised density dependent decay (where alpha < 0)
ddep1<- function(x, beta, alpha) {
	1/(1+exp(-alpha*(x-beta)))
}

pdf(paste(d.dir, "ABC/Posteriors.pdf", sep=""))
	#Alpha
	plot(density(gold[,"alpha"]), 
		xlim=c(-2, 0), 
		type="l", 
		main="Alpha",
		xlab="Alpha",
		bty="l")
	X<-seq(-2, 0, 0.1)
	lines(X, dunif(X, -2, 0),
		col="grey70")
		
	#Beta
	plot(density(gold[,"beta"]), 
		xlim=c(15, 60), 
		type="l", 
		main="Beta",
		xlab="Beta",
		bty="l")
	X<-seq(15, 60, 0.1)
	lines(X, dunif(X, 15, 60),
		col="grey70")
	
	#Prob.d
	X<-seq(0, 1, 0.01)
	plot(density(gold[,"prob.d"]),
		xlim=c(0,1),
		type="l",
		main="Dispersal probability",
		xlab="Dispersal probability",
		bty="l")
	lines(X, dunif(X, 0, 1),
		col="grey70")
	
	#Fsurv1
	yl<-range(c(dbeta(X, output[2,1], output[2,2])),
		gold[,"fsurv1"])
	plot(density(gold[,"fsurv1"]),
		xlim=c(0,1),
		ylim=yl,
		type="l",
		main="Female survival, year 1",
		xlab="Female survival, year 1",
		bty="l")
	lines(X, dbeta(X, output[2,1], output[2,2]),
		col="grey70")
		
	#Fsurv2
	yl<-range(c(dbeta(X, output[3,1], output[3,2])),
		gold[,"fsurv2"])
	plot(density(gold[,"fsurv2"]),
		xlim=c(0,1),
		ylim=yl,
		type="l",
		main="Female survival, year 2",
		xlab="Female survival, year 2",
		bty="l")
	lines(X, dbeta(X, output[3,1], output[3,2]),
		col="grey70")
	
	#Msurv
	yl<-range(c(dbeta(X, output[1,1], output[1,2])),
		gold[,"msurv"])
	plot(density(gold[,"msurv"]),
		xlim=c(0,1),
		ylim=yl,
		type="l",
		main="Male survival",
		xlab="Male survival",
		bty="l")
	lines(X, dbeta(X, output[1,1], output[1,2]),
		col="grey70")
		
	X<-seq(0, 100, 0.2)
	Y<-ddep1(X, mean(gold[,"beta"]), mean(gold[,"alpha"]))
	plot(Y~X, type="l")	
	
dev.off()