## Load functions ##
source(file="/home/jc227089/Quoll_model/Quoll_model_functions.R", echo=F)
###################################################################################
## Model Runs ##
args=(commandArgs(TRUE))

#evaluate the arguments
# input arguments will be number of sample runs (sr) and file.ID a number to differentiate files
for(i in 1:length(args)) {
	 eval(parse(text=args[[i]]))
}

#set the working directory
setwd("/home1/30/jc227089/Quoll_model/Random_outputs")

out<-c() #variable to take output
for (s in 1:sr){
	prob.d.s<-runif(1, 0, 1) #flat prior
	alpha.s<-runif(1, -2, 0) #flat prior
	fsurv1.s<-rnorm(1, 0.2338186, 0.1235535) #other priors
	while (fsurv1.s<0 | fsurv1.s>1) fsurv1.s<-rnorm(1, 0.2338186, 0.1235535) #catch stupid values
	fsurv2.s<-rnorm(1, 0.02020202, 0.03499093)
	while (fsurv2.s<0 | fsurv2.s>1) fsurv2.s<-rnorm(1, 0.02020202, 0.03499093)
	msurv.s<-rnorm(1, 0.04211957, 0.05893253)
	while (msurv.s<0 | msurv.s>1) msurv.s<-rnorm(1, 0.04211957, 0.05893253)
	beta.s<-runif(1, 5,40)
	fec.s<-rnorm(1, 6.725, 1.206140)
	temp<-small.mother(alpha=alpha.s, fsurv1=fsurv1.s, fsurv2=fsurv2.s, msurv=msurv.s, beta=beta.s, fec=fec.s, plot=F, prob.d=prob.d.s)
	temp<-c(alpha.s, fsurv1.s, fsurv2.s, msurv.s, beta.s, fec.s, prob.d.s, temp)
	out<-rbind(out, temp)
}
colnames(out)<-c("alpha", "fsurv1", "fsurv2", "msurv", "beta", "fec", "prob.d", "pob4", "pob5", "pob6", "pob7", 
		"pob.sr4", "pob.sr5", "pob.sr6", "pob.sr7", "as4", "as5", "as6", "as7", "as.sr4", "as.sr5", "as.sr6", "as.sr7")
save(out, file=paste("ABCSample", file.ID, ".RData", sep=""))

