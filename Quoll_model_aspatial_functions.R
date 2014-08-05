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
#ddep<-function(x){
#	(1-0.1)/((1)+exp(1*(x-2000)))+0.01
	#} 
	
bev.holt<-function(N, R, K){
 a<-(R-1)/K
 R/(1+a*N)
}  


#function to sample across a list
sample2<-function(v){
	if (is.null(v)==T) return(NULL)
	else sample(v, 1, T)	
}

#function that reproduces individuals

# takes population matrix
repro<-function(popmat, fec, init.p.var, h, alpha, beta){
	b.var<-h*init.p.var
	if (length(popmat[,1])==0) return(popmat) #if all dead then don't bother
# collect breeding females and their mates' B values
	male<-subset(popmat, popmat[,"S"]==1 & popmat[,"A"]>0) #matrix of all breeding males
	if (length(male)==0) return(popmat) #if are none then don't bother
	female<-subset(popmat, popmat[,"S"]==0 & popmat[,"A"]>0) #matrix of all adult females
	B.off<-sample(male[,"B"], length(female), replace=T) #sample a male for each female (replace=T)
	#if (length(B.off)!=length(female[,"B"])) print("error")#browser() #error catcher!
	B<-(B.off+female[,"B"])/2 #calculate offspring breeding values and cbind (midparent value)
	off.no<- bev.holt(length(female), fec, 2000) 
	off.no<- rpois(length(off.no), off.no) 
	print(off.no)
	B<-rep(B, off.no) #replicate midparent breeding values by the number of offspring each had
	B<-rnorm(length(B), B, (b.var/2)^0.5) #mid parent breeding values with variance of surviving offspring
	A<-rep(0, length(B)) # age=0
	S<-rbinom(length(A), 1, prob=0.5) # random sex allocation
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
	Juv<-subset(popmatrix, popmatrix[,"A"]==0) #Grabs little ones
	jsurv<- 1 #all survive
	if (selection==T) jsurv<-jsurv*fit.func(Juv[,"P"]) #toad relative fitness
	Juv<- subset(Juv, rbinom(length(Juv[,1]), 1, prob=jsurv)==1) #selects surviving offspring based on prob 
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

