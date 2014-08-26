d.dir<-"~/GitHub/gtq/"


load(paste(d.dir, "ABC/Kept_sims_poponlyhalfpcED.RData", sep=""))

load(paste(d.dir, "Priors/SurvivalPriors.rdata", sep=""))

#generalised density dependent decay (where alpha < 0)
ddep1<- function(x, beta, alpha) {
	1/(1+exp(-alpha*(x-beta)))
}

obs<-cbind(Year=c(2003, 2006:2009), Pob=c(19, 766, 818, 470, 300), Ast=c(45, 3228, 4820, 2590, 2890))
pred<-apply(gold[,7:14], 2, mean)
pred.sd<-apply(gold[,7:14], 2, sd)
pred.x<-rep(2006:2009, times=2)


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
		xlim=c(20, 180), 
		type="l", 
		main="Beta",
		xlab="Beta",
		bty="l")
	X<-seq(20, 180, 0.1)
	lines(X, dunif(X, 20, 180),
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
	yl<-range(c(dbeta(X, output[2,1], output[2,2]),
		density(gold[,"fsurv1"])$y))
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
	yl<-range(c(dbeta(X, output[3,1], output[3,2]),
		density(gold[,"fsurv2"])$y))
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
		
	X<-seq(0, 180, 0.2)
	Y<-ddep1(X, mean(gold[,"beta"]), mean(gold[,"alpha"]))
	plot(Y~X, type="l",
		xlab="Local density",
		ylab=expression(italic(p)),
		bty="l",
		main="Density dependence curve")
		
	yl<-c(0, max(pred+2*pred.sd))	
	matplot(x=obs[,1], y=obs[,2:3],
		pch=1:2,
		col=1,
		ylim=yl,
		bty="l",
		xlab="Year",
		ylab="Population size")
	points(pred.x, pred, 
		pch=rep(1:2, each=4),
		col="grey70")
	arrows(pred.x, pred-2*pred.sd, y1=pred+2*pred.sd,
		length=0,
		col="grey70")
	legend("topleft", 
		legend=c("Pobasoo", "Astell", "Observed", "Simulated"), 
		pch=c(1:2, NA, NA),
		fill=c(NA, NA, "black", "Grey70"),
		border=NA)
dev.off()