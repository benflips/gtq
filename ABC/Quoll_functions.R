### Functions and model for doing PVA of quolls under evolution ###

###Note: execution of quantitative genetics is (badly) incorrect and needs to be fixed before evolutionary scenarios can be implemented

#### Ella edited 12/8/14

# Initialises a matrix with n individuals
# individuals have a location (spX, spY) and are generated from
# initial mean breeding value, initial variance in breeding value, and heritability
init.inds<-function(n, spX, spY, init.b, init.p.var, h){
	X<-runif(n, 1, spX+1)%/%1 #grid ref
	Y<-runif(n, 1, spY+1)%/%1
	S<-rbinom(n, 1, 0.5) #sex
	A<-rep(1, n) #age
	b.var<-h*init.p.var
	B<-rnorm(n, init.b, b.var^0.5) # breeding values
	P<-rnorm(n, B, (init.p.var-b.var)^0.5)
	cbind(X, Y, S, A, B, P)		
}

# Function that defines relationship between Phenotype and Fitness (survival)
# threshold: P>0 survives, P<0 dies
fit.func<-function(P){
	P>0
}

#function that creates a list of vectors of male breeding values for each grid cell
b.list.maker<-function(popmat, spX, spY){
	male<-subset(popmat, popmat[,"S"]==1 & popmat[,"A"]>0)
	out<-vector("list", spX*spY)
	for (i in 1:length(male[,1])){
		out[[(male[i,"Y"]-1)*spX+male[i,"X"]]]<-c(out[[(male[i,"Y"]-1)*spX+male[i,"X"]]], male[i,"B"])
	}	
	out
}

#generalised density dependent decay (where alpha < 0)
ddep1<- function(x, beta, alpha) {
	1/(1+exp(-alpha*(x-beta)))
}

#function to sample across a list
sample2<-function(v){
	if (is.null(v)==T) return(NULL)
	else sample(v, 1, T)	
}

#function that reproduces individuals
#density dependance acting on fecundity according to grid cell K (beta)
# takes population matrix
repro<-function(popmat, spX, spY, init.p.var, h, alpha, beta){
	#browser()
	b.var<-h*init.p.var
	if (length(popmat[,1])==0) return(popmat)
	# collect breeding females and their mates' B values
	male<-subset(popmat, popmat[,"S"]==1 & popmat[,"A"]>0)
	if (length(male)==0) return(popmat)
	male<-table(factor(male[,"X"], levels=1:spX), factor(male[,"Y"], levels=1:spY)) #number of males present in each cell
	b.list<-b.list.maker(popmat, spX, spY) #list male Bs by space
	female<-subset(popmat, popmat[,"S"]==0 & popmat[,"A"]>0) #matrix of all adult females
	dens<-table(factor(female[,"X"], levels=1:spX), factor(female[,"Y"], levels=1:spY)) #work out the density of females in each grid cell	
	temp<-(female[,"Y"]-1)*spX+female[,"X"] #matrix indexes for male matrix and blist
	female<-subset(female, male[temp]>0) # matrix of females with available males
	if (nrow(female)==0) return(NULL)
	temp<-temp[which(male[temp]>0)] #update temp
	B.off<-unlist(lapply(b.list[temp], sample2)) #sample a male for each female (replace=T)
	if (length(B.off)==0) return(NULL)
	B.off<-(B.off+female[,"B"])/2 #calculate midparent value	
#matrix of females ordered by grid number
	female<- cbind(B.off, female) # bind offspring B values to their mums
	fdens<-dens[(female[,"Y"]-1)*spX+female[,"X"]] # Density each mum experiences
	dd.off<- ddep1(x=fdens, beta, alpha)
	off.no<-rbinom(nrow(female), 10, dd.off) # max offspring = 10

#get info for next generation and combine
	offspring<-matrix(female[rep(1:nrow(female), off.no),], ncol=7, dimnames=list(NULL, colnames(female)))
	if (is.null(nrow(offspring)) | is.na(nrow(offspring))) return(popmat)
	offspring[,"B"]<-offspring[,"B.off"]
	offspring[,"B"]<-rnorm(nrow(offspring), offspring[,"B"], (b.var/2)^0.5) #mid parent breeding values with variance of surviving offspring
	offspring[,"S"]<-rbinom(nrow(offspring), 1, 0.5) # random sex allocation
	offspring[,"A"]<-rep(0, nrow(offspring))	# age=0
	offspring[,"P"]<-rnorm(nrow(offspring), offspring[,"B"], (init.p.var-b.var)^0.5) #offspring phenotypes
	offspring<-offspring[,-1]
	popmat<-rbind(popmat, offspring)
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
age<-function(popmatrix, selection=F, alpha, fsurv1, fsurv2){
	if (length(popmatrix[,1])==0) return(popmatrix)
	#Survival of adults
	Ad<-subset(popmatrix, popmatrix[,"A"]>0) #Grabs all adults
	psurv<-rep(0, length(Ad[,1])) 
	psurv[which(Ad[,"S"]==0 & Ad[,"A"]==1)]<-fsurv1 #one year old females
	psurv[which(Ad[,"S"]==0 & Ad[,"A"]==2)]<-fsurv2 #two year old females
	psurv[which(Ad[,"S"]==1 & Ad[,"A"]==1)]<-1 #one year old males (already gone through die off)
	if (selection==T) psurv<-psurv*fit.func(Ad[,"P"]) #toad relative fitness
	Ad<-subset(Ad, rbinom(length(Ad[,1]), 1, psurv)==1) # probabilistic survival
	#Survival of juvs only if toads are present (otherwise survival=1)
	Juv<-subset(popmatrix, popmatrix[,"A"]==0)
	psurv<-1
	if (selection==T) psurv<-psurv*fit.func(Juv[,"P"]) #toad relative fitness
	Juv<-subset(Juv, rbinom(length(Juv[,1]), 1, psurv)==1) #surviving juveniles
	# gather survivors, age them and return
	popmatrix<-rbind(Ad, Juv)
	popmatrix[,"A"]<-popmatrix[,"A"]+1
	popmatrix
}


#calculates neighbours on a grid exclusing nonsense neighbours and returns a list of matrices of grid refs
#doesn't include source cell
neighbours.init<-function(spX, spY){
	neigh<-vector("list", length=spX*spY)
	for (x in 1:spX){
		for (y in 1:spY){
			out<-matrix(nrow=8, ncol=2)
			ifelse((x-1)==0, X<-NA, X<-x-1) #correct for indices going to zero
			ifelse((y-1)==0, Y<-NA, Y<-y-1)
			ifelse((y+1)>spY, Yp<-NA, Yp<-y+1) # correct for indices going to > space
			ifelse((x+1)>spX, Xp<-NA, Xp<-x+1)
			# Assign neighbours
			out[1,]<-c(Xp, y)
			out[2,]<-c(X, y)
			out[3,]<-c(x, Y)
			out[4,]<-c(Xp, Y)	
			out[5,]<-c(X, Y)
			out[6,]<-c(x, Yp)
			out[7,]<-c(Xp, Yp)
			out[8,]<-c(X, Yp)
			out<-subset(out, !is.na(apply(out, 1, sum)))
			neigh[[(y-1)*spX+x]]<-out
		}	
	}
	return(neigh)	#return the list
}

# samples a matrix and returns one row of the matrix
matrix.sample<-function(mtrx, size=1){
	mtrx[sample(seq(length(mtrx[,1])), size=size),]	
}

# executes nearest neighbour dispersal where individuals disperse to neighbouring cell with probability=prob.d 
# is called after age!
disperse<-function(popmatrix, n.list, prob.d, spX){
	if (length(popmatrix[,1])==0) return(popmatrix)
	disp<-rbinom(length(popmatrix[,1]), 1, prob.d) # calculate a probability of dispersal
	disp<-which(disp==1 & popmatrix[,"A"]==1)# & popmatrix[,"S"]==1) #get indices of surviving dispersing MALE juveniles
	sub.list<-n.list[(popmatrix[disp,"Y"]-1)*spX+popmatrix[disp,"X"]] #generate list of neighbour matrices
	sub.list<-lapply(sub.list, matrix.sample) #samples one row from each neighbour matrix
	sub.list<-do.call("rbind", sub.list) #places result into a matrix
	popmatrix[disp, 1:2]<-sub.list #assigns new locations back to popmatrix
	popmatrix
}

# runs the model
# dem.pars is a vector with alpha, fs1, fs2, msurv, beta, prob.d
mother<-function(n=2000, spX=10, spY=10, dem.pars, init.b=0, init.p.var=10, h=0.3, gens=50, sel.time=20, plot=FALSE){
	alpha<-dem.pars[1]
	fsurv1<-dem.pars[2]
	fsurv2<-dem.pars[3]
	msurv<-dem.pars[4]
	beta<-dem.pars[5]
	prob.d<-dem.pars[6]
	pop<-init.inds(n, spX, spY, init.b, init.p.var, h) # create a population
	n.list<-neighbours.init(spX, spY) # create a list of neighbours for each cell (individual?)
	popsize<-n
	sel<-FALSE
	for (g in 2:gens){
		pop<-repro(popmat=pop, spX=spX, spY=spY, init.p.var=init.p.var, h=h, alpha=alpha, beta=beta) # females reproduce (density dep)
		pop<-mdieoff(pop, msurv) # males die off
		if (g>sel.time) sel<-TRUE # have toads arrived?
		pop<-age(pop, selection=sel, alpha=alpha, fsurv1=fsurv1, fsurv2=fsurv2) # everyone gets a little older
		pop<-disperse(pop, n.list, prob.d, spX) # and the Juveniles trudge off
		if (length(pop[,1])==0) { # did we go extinct?
			pe<-TRUE
			print(paste("The population went extinct at generation ", g, ":-("))
			break	
		} else {pe<-FALSE }
		popsize<-c(popsize, length(pop[,1])) # new population size appended to popsize vector
		if (plot==T) plotter(pop, popsize, spX, spY, sel.time, gens, fid=g)
	}
	list(pop, popsize)	
	#return(pe)
}

# Runs the model wrt Pobassoo and Astell islands
small.mother<-function(alpha, fsurv1, fsurv2, msurv, beta, init.b=0, init.p.var=10, h=0.3, gens=7, prob.d, sel=FALSE, plot=FALSE){
	#Pobassoo
	# Treat Pobassoo 392 Ha as 5x5 space = 15.68 ha per cell
	spX.pob<-5
	pop.pob<-init.inds(19, spX=spX.pob, spY=spX.pob, init.b, init.p.var, h) #founders
	pop.pob[1:11,"S"]<-0
	pop.pob[12:19, "S"]<-1
	#pop.pob[,"X"]<-sample(1:4, nrow(pop.pob), replace=TRUE) #introduce into 4x4 space
	#pop.pob[,"Y"]<-sample(1:4, nrow(pop.pob), replace=TRUE)
	popsize.pob<-11 #to collect outputs
	sex.rat.pob<-11/8
	
	#Astell
	# treat Astell as 9x9 space = 15.61 ha per cell
	spX.as<-9 
	pop.as<-init.inds(45, spX=spX.as, spY=spX.as, init.b, init.p.var, h) #founders
	pop.as[,"X"]<-2 #start all founders at (2,1)
	pop.as[1:34,"S"]<-0
	pop.as[35:45, "S"]<-1
	#pop.as[,"X"]<-sample(1:7, nrow(pop.as), replace=TRUE) #introduce into 7x7 space
	#pop.as[,"Y"]<-sample(1:7, nrow(pop.as), replace=TRUE)
	popsize.as<-45 #to collect outputs
	sex.rat.as<-34/11
	
	n.list.pob<-neighbours.init(spX=spX.pob, spY=spX.pob)
	n.list.as<-neighbours.init(spX=spX.as, spY=spX.as)
	for (g in 2:gens){
		pop.pob<-repro(pop.pob, spX.pob, spY=spX.pob, init.p.var, h, alpha=alpha, beta=beta)
		pop.pob<-mdieoff(pop.pob, msurv)
		popsize.pob<-c(popsize.pob, length(which(pop.pob[,"S"]==0 & pop.pob[,"A"]>0)))
		sex.rat.pob<-c(sex.rat.pob, length(which(pop.pob[,"S"]==0 & pop.pob[,"A"]>0))/length(which(pop.pob[,"S"]==1 & pop.pob[,"A"]>0)))
		pop.pob<-age(pop.pob, sel, alpha, fsurv1, fsurv2)
		pop.as<-repro(pop.as, spX=spX.as, spY=spX.as, init.p.var, h, alpha=alpha, beta=beta)
		pop.as<-mdieoff(pop.as, msurv)
		popsize.as<-c(popsize.as, length(which(pop.as[,"S"]==0 & pop.as[,"A"]>0)))
		sex.rat.as<-c(sex.rat.as, length(which(pop.as[,"S"]==0 & pop.as[,"A"]>0))/length(which(pop.as[,"S"]==1 & pop.as[,"A"]>0)))
		pop.as<-age(pop.as, sel, alpha, fsurv1, fsurv2)
		pop.as<-disperse(pop.as, n.list.as, prob.d, spX.as)
		if (plot==T) plotter(pop.as, popsize.as, spX=spX.as, spY=spX.as, sel.time=999, gens, 22)
	}	
	keep<-4:7
	popsize.pob<-popsize.pob[keep]
	sex.rat.pob<-sex.rat.pob[keep]
	popsize.as<-popsize.as[keep]
	sex.rat.as<-sex.rat.as[keep]
	c(popsize.pob, sex.rat.pob, popsize.as, sex.rat.as)
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

