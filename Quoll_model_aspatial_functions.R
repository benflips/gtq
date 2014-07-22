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

####DONT NEED ####
#function that creates a list of vectors of male breeding values for each grid cell
#b.list.maker<-function(popmat, spX, spY){
#	male<-subset(popmat, popmat[,"S"]==1 & popmat[,"A"]>0)
#	out<-vector("list", spX*spY)
#	for (i in 1:length(male[,1])){
#		out[[(male[i,"Y"]-1)*spX+male[i,"X"]]]<-c(out[[(male[i,"Y"]-1)*spX+male[i,"X"]]], male[i,"B"])
#	}	
#	out
#}


#generalised density dependent decay (where alpha < 0) 
ddep<-function(beta, Nt){
	1/(1+exp(3*(Nt-beta)))	
}

#function to sample across a list
sample2<-function(v){
	if (is.null(v)==T) return(NULL)
	else sample(v, 1, T)	
}


#function that reproduces individuals

#install.packages("survival")
#library(survival)

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
# calc no. of offspring per couple
	B<-(mat[,1]+mat[,2])/2
	off.no<- rpois(length(mat), fec) # no. offspring per couple (stochastic around expected fecundity)
	#if off.no>6 return off.no=6
	#parents<- cbind(B, off.no)
	
	len<- sum(off.no) #sum of total offspring
	#index<- rep(0, each=sum)
	#index<- cumsum(index, rep(1:end-1)
	
	#index <- zeros(1,sum(off.no))
	#index(cumsum([1 rep(1:end-1)])) = 1
	#index <- cumsum(index)
	#off <- B(index)
	
	
	off<-rep(B, times=off.no, length.out=len) #replicate B for no. of offspring per parent (length is sum of off.no)
	#N<- rep(off[off.no], length(B))
	
	#sur<- cbind(off, N)
	#osurv<-sapply(sur, ddep(alpha, beta, N))	
	#off<-subset(off, rbinom(length(off[,1]), 1, osurv)==1) # probabilistic survival

	
	B<-rnorm(length(off), B, (b.var/2)^0.5) #mid parent breeding values with vairance
	S<-rbinom(length(off), 1, 0.5) # random sex allocation
	A<-rep(0, length(off))	# age=0
	P<-rnorm(length(B), B, (init.p.var-b.var)^0.5)
	off<-cbind(S, A, B, P)
	popmat<-rbind(popmat, off)
	popmat
}
	 
#discarded script
	#osurv<-1/(1+exp(3*(off.no[,"1"]-7)	
	#parents <- cbind(B, off.no) #parent breeding values, offspring no., mid parent values, offspring survival prob
	#off<-subset(parents, rbinom(length(parents[,1]), 1, osurv)==1) #surviving juveniles
	
	#off.no.num<- data.matrix(off.no, rownames.force = NA)	
	#osurv<- ddep(beta, off.no.num) # offspring survival prob
	#if (off.no > 8){ sur = 1/off.no} else if (off.no<8) {sur=1} #sur = survival probability 
#create offspring
	#off<-rep(parents, each=off.no) #replicate B.off and osurv for no of offspring per parent
	#off<- subset(off, osurv>0) # get rid of dead ones
	
	
	##PROBLEM: non-numeric argument for ddep formula
	#don't need this line? B.off<- rnorm(length(B.off), B.off, (b.var/2)^0.5) #variance in offspring b values ## mistake? might need to be fixed. 
	
	#E<- subset(E<8) #carrying capacity=8
	#off<- cbind(B.off, E)
	#osurv<- survfit(B.off~1/(1+exp(3*(E-7))), data=off)
	#off<- subset(off[,E>0])
	#osurv<- ddep(alpha, beta, E)
	#off<- subset(osurv>1)
	#if E>6 osurv<- ddep(alpha, beta, E)
	#B1<-rep(B.off, each=E)
	#osurv<- ddep(alpha, beta, E) 
	#osurv<- ddep(alpha, beta, E) #density dependant offspring survival depending on fec
	#off<-subset(osurv, rbinom(length(E[,1]), 1, osurv)==1) # gather surviving offspring
	#S<-rbinom(length(B1), 1, 0.5) # random sex allocation
	#A<-rep(0, length(B1))	# age=0
	#B<-rnorm(length(B1), B1, (b.var/2)^0.5) #offspring breeding values = female B.off* b.var # Needs to have variation within familes
	#P<-rnorm(length(B), B, (init.p.var-b.var)^0.5)
	#off<-cbind(S, A, B, P)
	#popmat<-rbind(popmat, off)
	#popmat
#}


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
	#Survival of juvs (density dependent)
	Juv<-subset(popmatrix, popmatrix[,"A"]==0)
	psurv<- rep(0, length(Juv[,1]))
	psurv[which(Juv[,"A"]==0)]<-jsurv
	Juv<-subset(Juv, rbinom(length(Juv[,1]), 1, psurv)==1) # probabilistic survival
	#no Juv survival added - die off according to density in repro 
	#density<- nrow(Juv)
	#turn off density dep
	#density<-table(factor(Ad[,"X"], levels=1:spX), factor(Ad[,"Y"], levels=1:spY)) # adult density through space
	#density<-density+table(factor(Juv[,"X"], levels=1:spX), factor(Juv[,"Y"], levels=1:spY)) #plus juvenile density
	#density<-density/K # Density is measured relative to carrying capacity/habitat size
	#temp<-(Juv[,"Y"]-1)*spX+Juv[,"X"] #matrix indexes for density matrix
	#psurv<-ddep(alpha, beta, density[,1]) #juvenile survival
	#if (selection==T) psurv<-psurv*fit.func(Juv[,"P"]) #toad relative fitness
	#Juv<-subset(Juv, rbinom(length(Juv[,1]), 1, psurv)==1) #surviving juveniles
	# gather survivors, age them and return
	popmatrix<-rbind(Ad, Juv)
	popmatrix[,"A"]<-popmatrix[,"A"]+1
	popmatrix
}


# samples a matrix and returns one row of the matrix
matrix.sample<-function(mtrx, size=1){
	mtrx[sample(seq(length(mtrx[,1])), size=size),]	
}

#dont need dispersal
# executes nearest neighbour dispersal where individuals disperse to neighbouring cell with probability=prob.d 
# is called after age!
#disperse<-function(popmatrix, n.list, prob.d, spX){
	#if (length(popmatrix[,1])==0) return(popmatrix)
	#disp<-rbinom(length(popmatrix[,1]), 1, prob.d) # calculate a probability of dispersal
	#disp<-which(disp==1 & popmatrix[,"A"]==1) #get indices of surviving dispersing juveniles
	#sub.list<-n.list[(popmatrix[disp,"Y"]-1)*spX+popmatrix[disp,"X"]] #generate list of neighbour matrices
	#sub.list<-lapply(sub.list, matrix.sample) #samples one row from each neighbour matrix
	#sub.list<-do.call("rbind", sub.list) #places result into a matrix
	#popmatrix[disp, 1:2]<-sub.list #assigns new locations back to popmatrix
	#popmatrix
#}

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




######################

#runs model over generations
mother<-function(n=2000, spX=10, spY=10, alpha=-1, fsurv1=0.4, fsurv2=0.1, msurv=0.1, init.b=0, init.p.var=10, h=0.3, gens=40, beta=100, K=400, fec=4, prob.d=0.1, sel.time=10, plot=FALSE){
	pop<-init.inds(n, spX, spY, init.b, init.p.var, h)
	n.list<-neighbours.init(spX, spY)
	popsize<-n
	sel<-FALSE
	for (g in 2:gens){
		pop<-repro(pop, spX, spY, fec, init.p.var, h)
		pop<-mdieoff(pop, msurv)
		if (g>sel.time) sel<-TRUE
		pop<-age(pop, selection=sel, alpha=alpha, fsurv1=fsurv1, fsurv2=fsurv2, msurv=msurv,spX=spX, spY=spY, beta=beta, K=K)
		if (spX!=1 & spY!=1) pop<-disperse(pop, n.list, prob.d, spX)
		if (length(pop[,1])==0) {
			popsize<-c(popsize, 0)
			print(paste("The population went extinct at generation ", g, ":-("))
			break	
		}
		popsize<-c(popsize, length(pop[,1]))
		if (plot==T) plotter(pop, popsize, spX, spY, sel.time)
	}
	list(pop, popsize)	
}

small.mother<-function(spX=1, spY=1, alpha, fsurv1, fsurv2, msurv, beta, K, fec, init.b=0, init.p.var=10, h=0.3, gens=7, prob.d=0.5, sel=FALSE, plot=FALSE){
	#Pobassoo
	pop.pob<-init.inds(19, spX, spY, init.b, init.p.var, h) #founders
	pop.pob[1:11,"S"]<-0
	pop.pob[12:19, "S"]<-1
	popsize.pob<-11 #to collect outputs
	sex.rat.pob<-11/8
	K.pob<-392
	#Astell
	spX.as<-3 #treat astell as 1x3 space
	pop.as<-init.inds(45, spX, spY, init.b, init.p.var, h) #founders
	pop.as[,"X"]<-2 #start all founders at (2,1)
	pop.as[1:34,"S"]<-0
	pop.as[35:45, "S"]<-1
	popsize.as<-45 #to collect ouztputs
	sex.rat.as<-34/11
	K.as<-1265/spX.as
	n.list.as<-neighbours.init(spX=spX.as, spY)
	for (g in 2:gens){
		pop.pob<-repro(pop.pob, spX, spY, fec, init.p.var, h)
		pop.pob<-mdieoff(pop.pob, msurv)
		popsize.pob<-c(popsize.pob, length(which(pop.pob[,"S"]==0 & pop.pob[,"A"]>0)))
		sex.rat.pob<-c(sex.rat.pob, length(which(pop.pob[,"S"]==0 & pop.pob[,"A"]>0))/length(which(pop.pob[,"S"]==1 & pop.pob[,"A"]>0)))
		pop.pob<-age(pop.pob, sel, alpha, fsurv1, fsurv2, msurv,spX=spX, spY, beta, K.pob)
		pop.as<-repro(pop.as, spX=spX.as, spY, fec, init.p.var, h)
		pop.as<-mdieoff(pop.as, msurv)
		popsize.as<-c(popsize.as, length(which(pop.as[,"S"]==0 & pop.as[,"A"]>0)))
		sex.rat.as<-c(sex.rat.as, length(which(pop.as[,"S"]==0 & pop.as[,"A"]>0))/length(which(pop.as[,"S"]==1 & pop.as[,"A"]>0)))
		pop.as<-age(pop.as, sel, alpha, fsurv1, fsurv2, msurv,spX=spX.as, spY, beta, K.as)
		pop.as<-disperse(pop.as, n.list.as, prob.d, spX.as)
		if (plot==T) plotter(pop.as, popsize.as, spX, spY, sel.time=999)
	}	
	keep<-4:7
	popsize.pob<-popsize.pob[keep]
	sex.rat.pob<-sex.rat.pob[keep]
	popsize.as<-popsize.as[keep]
	sex.rat.as<-sex.rat.as[keep]
	c(popsize.pob, sex.rat.pob, popsize.as, sex.rat.as)
}



