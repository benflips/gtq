## Load functions ##
source(file="/home/jc227089/Quoll_model/ABC/Quoll_functions.R", echo=F)

load("/home/jc227089/Quoll_model/Priors/SurvivalPriors.rdata")
###################################################################################
## Model Runs ##
args=(commandArgs(TRUE))

#evaluate the arguments
# input arguments will be number of sample runs (sr) and file.ID a number to differentiate files
for(i in 1:length(args)) {
	 eval(parse(text=args[[i]]))
}

#set the working directory
setwd("/home/jc227089/Quoll_model/Random_outputs")

out<-c() #variable to take output
for (s in 1:sr){
	prob.d.s<-runif(1, 0, 1) #flat prior
	alpha.s<-runif(1, -2, 0) #flat prior
	fsurv1.s<-rbeta(1, output[2,1], output[2,2]) #other priors
	fsurv2.s<-rbeta(1, output[3,1], output[3,2]) 
	msurv.s<-rbeta(1, output[1,1], output[1,2]) 
	beta.s<-runif(1, 20,100)
	temp<-small.mother(alpha=alpha.s, fsurv1=fsurv1.s, fsurv2=fsurv2.s, msurv=msurv.s, beta=beta.s, plot=F, prob.d=prob.d.s)
	temp<-c(alpha.s, fsurv1.s, fsurv2.s, msurv.s, beta.s, prob.d.s, temp)
	out<-rbind(out, temp)
}
colnames(out)<-c("alpha", "fsurv1", "fsurv2", "msurv", "beta", "prob.d", "pob4", "pob5", "pob6", "pob7", 
		"pob.sr4", "pob.sr5", "pob.sr6", "pob.sr7", "as4", "as5", "as6", "as7", "as.sr4", "as.sr5", "as.sr6", "as.sr7")
save(out, file=paste("ABCSample", file.ID, ".RData", sep=""))

