

###################################################################################
## Model Runs ##
## Model Runs ##
args=(commandArgs(TRUE))

#evaluate the arguments
# input arguments will be number of sample runs (sr) and file.ID a number to differentiate files
for(i in 1:length(args)) {
	 eval(parse(text=args[[i]]))
}

#set the working directory
setwd("/home1/30/jc227089/Quoll_model/Random_outputs")
print(file.ID)
print(sr)

# fitted values
load("/home1/30/jc227089/Quoll_model/Posteriors halfpc ED Summary.RData")
ps<-out
rm(out)

size<-1:5 #size of space, in grid cells where each grid cell = 400 ha
extime<-c()
for (s in 1:sr){
	#Sample from posteriors
	alpha.s<-runif(1, ps["tf.alpha",5], min(ps["tf.alpha",6], 0)) #flat posterior
	fsurv1.s<-rnorm(1, ps["tf.fsurv1",1], ps["tf.fsurv1", 2]) #others
	while (fsurv1.s<0 | fsurv1.s>1) fsurv1.s<-rnorm(1, ps["tf.fsurv1",1], ps["tf.fsurv1", 2]) #catch stupid values
	fsurv2.s<-rnorm(1, ps["tf.fsurv2",1], ps["tf.fsurv2", 2]) 
	while (fsurv2.s<0 | fsurv2.s>1) fsurv2.s<-rnorm(1, ps["tf.fsurv2",1], ps["tf.fsurv2", 2]) 
	msurv.s<-rnorm(1, ps["tf.msurv",1], ps["tf.msurv", 2]) 
	while (msurv.s<0 | msurv.s>1) msurv.s<-rnorm(1, ps["tf.msurv",1], ps["tf.msurv", 2]) 
	beta.s<-rnorm(1, ps["tf.beta",1], ps["tf.beta", 2])
	fec.s<-rnorm(1, ps["tf.fec",1], ps["tf.fec", 2])
	dprob.s<-rnorm(1, ps["tf.probd",1], ps["tf.probd", 2])
	while (dprob.s<0 | dprob.s>1) dprob.s<-rnorm(1, ps["tf.probd",1], ps["tf.probd", 2]) 
	rep<-c() #to catch each size rep
	for (ss in 1:length(size)){
		temp<-mother(n=50*size[ss]^2, spX=size[ss], spY=size[ss], gens=100, alpha=alpha.s, fsurv1=fsurv1.s, fsurv2=fsurv2.s, msurv=msurv.s, beta=beta.s, fec=fec.s, prob.d=dprob.s, K=400)[[2]]
		temp<-c(length(temp), size[ss])
		rep<-rbind(rep, temp)
	}
	extime<-rbind(extime, rep)
}
	
save(extime, file=paste("extime", file.ID, ".RData", sep=""))