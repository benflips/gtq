### Functions and model for doing PVA of quolls under evolution ###

###Note: execution of quantitative genetics is (badly) incorrect and needs to be fixed before evolutionary scenarios can be implemented

# Initialises a matrix with n individuals
# individuals have a location (spX, spY) and are generated from
# initial mean breeding value, initial variance in breeding value, and heritability
init.inds<-function(n, init.b, init.p.var, h){
	S<-rbinom(n, 1, 0.5) #sex
	A<-rep(1, n) #age
	b.var<-h*init.p.var
	B<-rnorm(n, init.b, b.var^0.5) # breeding values
	P<-rnorm(n, B, (init.p.var-b.var)^0.5)
	cbind(S, A, B, P)		
}

# Function that defines relationship between Phenotype and Fitness (survival)
# threshold: P>0 survives, P<0 dies
fit.func<-function(P){
	P>0
}

#generalised density dependent decay of offspring
ddep<-function(x){
	1/(1+exp(3*(x-7)))	
} 

#function that reproduces individuals

# takes population matrix
repro<-function(popmat, fec, init.p.var, h, alpha){
	b.var<-h*init.p.var
	if (length(popmat[,1])==0) return(popmat) #if all dead then don't bother
# collect breeding females and their mates' B values
	male<-subset(popmat, popmat[,"S"]==1 & popmat[,"A"]>0) #matrix of all breeding males
	if (length(male)==0) return(popmat) #if are none then don't bother
	female<-subset(popmat, popmat[,"S"]==0 & popmat[,"A"]>0) #matrix of all adult females
	femB<- female[,c("B")] #extract fem breeding values
	malB<- male[,c("B")] #extract male breeding values
	mat<-list(femB, malB) #bind
	mat<-sapply(mat, '[', seq(min(sapply(mat,length)))) # crop parent breeding values to min length
# calc no. of offspring per couple and their survival (density dependant)
	B<-(mat[,1]+mat[,2])/2
	off.no<- rpois(length(mat), fec) # no. offspring per couple (stochastic around expected fecundity)
	osurv<- lapply(off.no, ddep) #density dep survival (off.no determines off survival prob)
	len<- sum(off.no) #sum of total offspring
	B<-rep(B, times=off.no, length.out=len) #replicate B for no. of offspring per parent (length is sum of off.no)
	osurv<- rep(osurv, times=off.no, length.out=len) #replicate osurv for no. of offspring per parent
	B<-subset(B, rbinom(length(B[,"1"]), 1, osurv)==1) #surviving juveniles
# collect surviving offspring 
	B<-rnorm(length(Juv), Juv, (b.var/2)^0.5) #mid parent breeding values with variance of surviving offspring
	S<-rbinom(length(B), 1, 0.5) # random sex allocation
	A<-rep(0, length(B))	# age=0
	P<-rnorm(length(B), B, (init.p.var-b.var)^0.5) #phenotypic variation
	off<-cbind(S, A, B, P)
	popmat<-rbind(popmat, off)
	popmat
}
	 

# runs the male die off
mdieoff<-function(popmatrix, msurv){
	psurv<-rep(1, length(popmatrix[,1]))
	psurv[which(popmatrix[,"S"]==1 & popmatrix[,"A"]==1)]<-msurv #one year old males
	psurv[which(popmatrix[,"S"]==1 & popmatrix[,"A"]>1)]<-0 #two year old males
	popmatrix<-subset(popmatrix, rbinom(length(popmatrix[,1]), 1, psurv)==1) # probabilistic survival
	popmatrix
}

#ages or kills individuals based on age-specific survival plus fitness
age<-function(popmatrix, selection=F, alpha, fsurv1, fsurv2, jsurv, msurv, beta, K){
	if (length(popmatrix[,1])==0) return(popmatrix)
	#Survival of adults
	Ad<-subset(popmatrix, popmatrix[,"A"]>0) #Grabs all adults
	psurv<-rep(0, length(Ad[,1])) 
	psurv[which(Ad[,"S"]==0 & Ad[,"A"]==1)]<-fsurv1 #one year old females
	psurv[which(Ad[,"S"]==0 & Ad[,"A"]==2)]<-fsurv2 #two year old females
	psurv[which(Ad[,"S"]==1 & Ad[,"A"]==1)]<-1 #one year old males (already gone through die off)
	if (selection==T) psurv<-psurv*fit.func(Ad[,"P"]) #toad relative fitness
	Ad<-subset(Ad, rbinom(length(Ad[,1]), 1, psurv)==1) # probabilistic survival
#Survival of juvs (not density dependent)
	Juv<-subset(popmatrix, popmatrix[,"A"]==0)
	psurv<- rep(0, length(Juv[,1]))
	psurv[which(Juv[,"A"]==0)]<-jsurv #juvenile survival probability (supposedly quite high)
	if (selection==T) psurv<-psurv*fit.func(Juv[,"P"]) #toad relative fitness
	Juv<-subset(Juv, rbinom(length(Juv[,1]), 1, psurv)==1) # probabilistic survival
# gather survivors, age them and return
	popmatrix<-rbind(Ad, Juv)
	popmatrix[,"A"]<-popmatrix[,"A"]+1
	popmatrix
}


# samples a matrix and returns one row of the matrix
matrix.sample<-function(mtrx, size=1){
	mtrx[sample(seq(length(mtrx[,1])), size=size),]	
}


#plotter function
plotter<-function(popmatrix, popsize, sel.time){
	#p.summ<-table(factor(popmatrix[,"X"], levels=1:spX), factor(popmatrix[,"Y"], levels=1:spY))
	#p.summ<-tapply(popmatrix[,"P"], list(factor(popmatrix[,"X"], levels=1:spX), factor(popmatrix[,"Y"], levels=1:spY)), mean)
	par(mfrow=c(2, 2))
	image(p.summ, col=rev(heat.colors(10)), breaks=seq(0, 100, length=11), main="Population density")
	hist(popmatrix[,"P"], xlim=c(-30, 30), main="Phenotype")
	plot(popsize, type="l", ylab="Population Size", xlab="Generations", bty="l")
	arrows(sel.time,0,sel.time,100000, length=0, col="red")
	hist(popmatrix[,"B"], xlim=c(-30, 30), main="Genotype")
}

