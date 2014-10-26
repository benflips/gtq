### Functions and model for doing PVA of quolls under evolution ###

###Note: execution of quantitative genetics is (badly) incorrect and needs to be fixed before evolutionary scenarios can be implemented

##fitted values FOR MY LAPTOP
s.dir<-"/Users/ellakelly/GitHub/gtq/ABC/" #source directory

#### FOR BEN's COMP 

#s.dir<-"/Users/ellakelly/Documents/GitHub/gtq/ABC/" #source directory

###################################################################################
## Model Runs ##

# fitted values for survival and other things
load(paste(s.dir, "Kept_sims_poponlyhalfpcED.RData", sep=""))
pars<-apply(gold[,1:6], 2, mean)

# equilibrium population density 
s.dir<-"/Users/ellakelly/GitHub/gtq/Priors/" #source directory
load(paste(s.dir, "init.density", sep=""))
init.pop<-sum(test3[,3])

# Initialises a matrix with n individuals
# individuals have a location (spX, spY) and are generated from
# initial mean breeding value, initial variance in breeding value, and heritability
init.inds<-function(n, spX, spY, init.b, init.p.var, h){
	random<-data.frame(Y=test3[,1], X=test3[,2], dens=sample(test3[,3]))
	Y<-rep(random[,1],times=random[,"dens"]) #X and Y from equilibrium pop
	X<-rep(random[,2],times=random[,"dens"])
	S<-rbinom(n, 1, 0.5) #sex
	A<-rep(1, n) #age
	b.var<-h*init.p.var
	B<-rnorm(n, init.b, b.var^0.5) # breeding values
	P<-rnorm(n, B, (init.p.var-b.var)^0.5)
	HI<-rep(0, n) #hybrid index, all start at 0
	cbind(X, Y, S, A, B, P, HI)	
}

###translocated population: ALL PUT IN 5x5 ON GRID
trans.inds<-function(n=nT, spX, spY, init.b, init.p.var, h){
	Y<-rep(5, n) 
	X<-rep(5, n) 
	S<-rbinom(n, 1, 0.5) #sex
	A<-rep(1, n) #age
	b.var<-h*init.p.var
	B<-rnorm(n, init.b, b.var^0.5) # breeding values
	P<-rnorm(n, B, (init.p.var-b.var)^0.5)
	HI<-rep(1, n) #hybrid index, all start at 0
	cbind(X, Y, S, A, B, P, HI)	
}

# Function that defines relationship between Phenotype and Fitness (survival)
# threshold: P>0 survives, P<0 dies
fit.func<-function(P){
	P>0
}

#generalised density dependent decay (where alpha < 0)
ddep1<- function(x, beta, alpha) {
	1/(1+exp(-alpha*(x-beta)))
}

#sampling males
randomRows = function(df,n){
   return(df[sample(nrow(df),n,replace=TRUE),])
}

# subset males based on available females (replace=true)
sample3<-function(list, density, spX, spY){
	out<-vector("list", length=length(density))
	for (i in 1:length(list)){
	out[[i]]<-
	if ((density[i,])=="0") return(NULL)
	else randomRows(list[[i]], density[i,])
	}
	out
}
	
#subset females if there are no males in a cell 	
densequal<-function(mdens,fdens){
	out<-matrix(nrow=length(fdens), ncol=1)
for(i in 1:(length(fdens))){
	if((mdens[i,])=="0") {out[i,]<-0}
	else {out[i,]<- fdens[i,]}
	}
	out
}
	
#create matrix with midparent value and hybrid index	
mating<-function(m,f,x){
	out<-matrix("numeric", ncol=2, nrow=length(x))
	mB<- unlist(sapply(m, `[[`, 1))
	fB<- unlist(sapply(f, `[[`, 1))
	mHI<- unlist(sapply(m, `[[`, 2))
	fHI<- unlist(sapply(f, `[[`, 2))
	fX<- unlist(sapply(f, `[[`, 4))
	fY<- unlist(sapply(f, `[[`, 5))	
	B.off<-rowMeans((cbind(mB,fB)), na.rm = FALSE, dims = 1)
	HI.off<-rowMeans((cbind(mHI,fHI)), na.rm = FALSE, dims = 1)
	out<-cbind(B.off, HI.off,fX,fY)
	out
}

#function that reproduces individuals
#density dependance acting on fecundity according to grid cell K (beta)
# takes population matrix
repro<-function(popmat, spX, spY, init.p.var, h, alpha, beta){
	#browser()
	b.var<-h*init.p.var
	if (length(popmat[,1])==0) return(popmat)
# collect breeding females and their mates' B and HI values
	#males
	male<-subset(popmat, popmat[,"S"]==1 & popmat[,"A"]>0) 
	if (length(male)==0) return(popmat)
	mdens<-as.matrix(as.vector(table(factor(male[,"X"], levels=1:spX), factor(male[,"Y"], levels=1:spY)))) #number of males present in each cell
	B<-male[,"B"]
	HI<-male[,"HI"]
	XY<-(male[,"Y"]-1)*spX+male[,"X"]
	df<-data.frame(B,HI,XY)
	M.list <- split(df , f = df$XY)	#list males HI and Bd by space
	#females
	female<-subset(popmat, popmat[,"S"]==0 & popmat[,"A"]>0) #matrix of all adult females
	if (length(female)==0) return(popmat)
	fdens<-as.matrix(as.vector(table(factor(female[,"X"], levels=1:spX), factor(female[,"Y"], levels=1:spY)))) #work out the density of females in each grid cell	
	fdens<- densequal(mdens, fdens) #subset females if males in cell = zero
	if (sum(fdens)==0) return(popmat)
# who's available for breeding?	
	male.list<-sample3(M.list, fdens, spX, spY) #sample males depending on available females
	fXY<-(female[,"X"]-1)*spX+female[,"Y"]
	df<-data.frame(female[,"B"],female[,"HI"],fXY, female[,"X"], female[,"Y"])
	female.list <- split(df , f= df$fXY) #list female B and HI values by space 
	female.list <- sample3(female.list, fdens, spX, spY) #sample females depending on fdens
#mating (mean of B and HI values) 
	mated<-mating(male.list, female.list, female)
#density dependance
	dd.off<- ddep1(x=fdens, beta, alpha)
	off.no<-rbinom(nrow(mated), 10, dd.off) # max offspring = 10
	off<-mated[rep(1:nrow(mated), off.no),]
#get info for next generation and combine
	offspring<- matrix(nrow=length(off), ncol=ncol(female), dimnames=list(NULL, colnames(female)))
	if (is.null(nrow(offspring)) | is.na(nrow(offspring))) return(popmat)
	offspring[,"B"]<-off[,"B.off"]
	offspring[,"B"]<-rnorm(nrow(offspring), offspring[,"B"], (b.var/2)^0.5) #mid parent breeding values with variance of surviving offspring
	offspring[,"S"]<-rbinom(nrow(offspring), 1, 0.5) # random sex allocation
	offspring[,"A"]<-rep(0, nrow(offspring))	# age=0
	offspring[,"X"]<-off[,"fX"]
	offspring[,"Y"]<-off[,"fY"]
	offspring[,"P"]<-rnorm(nrow(offspring), offspring[,"B"], (init.p.var-b.var)^0.5) #offspring phenotypes
	offspring[,"HI"]<-off[,"HI.off"]
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

#Hybrid index fitness (F1 lowered fitness)
HI_fitness<-function(HI, s, beta){
  W<-1-s*(4*HI*(1-HI))^beta
  W
}
#Or other function (where F2 have lowered fitness)
#HI_fitness<-function(HI, z){
# W<-(1-z)+zCos(x*4pi)
 #  W
#}

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
#Juvenile survival
	#hybrid index
	#Juv<-subset(popmatrix, popmatrix[,"A"]==0)
	psurv<- HI_fitness(popmatrix[,"HI"], 1, 10)
	popmatrix<-subset(popmatrix, rbinom(length(popmatrix[,1]), 1, psurv)==1) #surviving juveniles
	#if toads are present	
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
mother<-function(n=init.pop, spX=10, spY=10, dem.pars, init.b=-5, init.p.var=10, h=0.3, gens=50, sel.time=20, plot=FALSE, hybrid=0, nT=100, trans.time=14){
	alpha<-dem.pars[1]
	fsurv1<-dem.pars[2]
	fsurv2<-dem.pars[3]
	msurv<-dem.pars[4]
	beta<-dem.pars[5]
	prob.d<-dem.pars[6] #assign dem.pars
	pop<-init.inds(n, spX, spY, init.b, init.p.var, h) # create a population
	n.list<-neighbours.init(spX, spY) # create a list of neighbours for each cell (individual?)
	popsize<-n
	sel<-FALSE
	trans.pop<- trans.inds(nT, spX, spY, init.b, init.p.var, h)# init.b will be something specific to "QLD" quolls, and n will be variable
	for (g in 2:gens){
		if (g==trans.time) pop<-rbind(pop, trans.pop)
		pop<-repro(popmat=pop, spX=spX, spY=spY, init.p.var=init.p.var, h=h, alpha=alpha, beta=beta) # females reproduce (density dep)
		pop<-mdieoff(pop, msurv) # males die off
		if (g>sel.time) sel<-TRUE # have toads arrived?
		pop<-age(pop, selection=sel, alpha=alpha, fsurv1=fsurv1, fsurv2=fsurv2) # everyone gets a little older
		pop<-disperse(pop, n.list, prob.d, spX) # and the Juveniles trudge off
		if (length(pop[,1])==0) { # did we go extinct?
			pe<-TRUE
			print(paste("The population went extinct at generation ", g, ":-("))
			break	
		} else {pe<-FALSE			
 			}
		popsize<-c(popsize, length(pop[,1])) # new population size appended to popsize vector
		if (plot==T) plotter(pop, popsize, spX, spY, sel.time, gens, fid=g)
	}
if (length(pop[,1])>0) {		
		print(mean(pop[,"B"])) ## print final B values for surviving populations
 		}
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

mother(dem.pars=pars, init.b=-8.5, gens=50, sel.time=20, h=0.3, plot=T)