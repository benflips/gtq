# Here we estimate survival probability from the assembled data
# survival is treated as Beta distributed
# Use JAGS to sample the posterior of the beta and then fitdistr to generate parameters that can be fed to the ABC process

library(MASS)
library(rjags)

#setwd("/Users/ellakelly/Documents/PhD/gtq/Priors")
setwd("/home/jc227089/Quoll_model/Priors")


#litter size
d<-read.table("Litter size.txt", header=T, sep="\t")

# These data aren't easily turned into useful numbers for the Binomial
# reading through Braithwaite and Griffiths, max litter size seems to be 10.


#Male survival
d<-read.table("Male survival.txt", header=T, sep="\t")
nrw<-nrow(d)
s<-d$surv
n<-d$n
a<-jags.model("BetaEst.txt", data=list(nrw=nrw, s=s, n=n), n.adapt=10000)
a2<-coda.samples(a, c("p"), n.iter=50000)
msFit<-fitdistr(a2[[1]][,1], densfun="beta", start=list(shape1=1, shape2=1))
# check fit
pdf()
	hist(a2[[1]][,1], freq=F)
	lines(x=seq(0, 1, 0.01), y=dbeta(seq(0, 1, 0.01), msFit$estimate[1], msFit$estimate[2]))
dev.off()


#Female survival year 1
d<-read.table("female survival yr1.txt", header=T, sep="\t")
nrw<-nrow(d)
s<-d$surv
n<-d$n
a<-jags.model("BetaEst.txt", data=list(nrw=nrw, s=s, n=n), n.adapt=10000)
a2<-coda.samples(a, c("p"), n.iter=50000)
fs1Fit<-fitdistr(a2[[1]][,1], densfun="beta", start=list(shape1=1, shape2=1))
pdf()
	hist(a2[[1]][,1], freq=F)
	lines(x=seq(0, 1, 0.01), y=dbeta(seq(0, 1, 0.01), fs1Fit$estimate[1], fs1Fit$estimate[2]))
dev.off()


#Female survival year 2
d<-read.table("female survival yr2.txt", header=T, sep="\t")
nrw<-nrow(d)
s<-d$surv
n<-d$n
a<-jags.model("BetaEst.txt", data=list(nrw=nrw, s=s, n=n), n.adapt=10000)
a2<-coda.samples(a, c("p"), n.iter=50000)
fs2Fit<-fitdistr(a2[[1]][,1], densfun="beta", start=list(shape1=1, shape2=1))
pdf()
	hist(a2[[1]][,1], freq=F)
	lines(x=seq(0, 1, 0.01), y=dbeta(seq(0, 1, 0.01), fs2Fit$estimate[1], fs2Fit$estimate[2]))
dev.off()


output<-msFit$estimate
output<-rbind(output, fs1Fit$estimate)
output<-rbind(output, fs2Fit$estimate)
row.names(output)<-c("ms", "fs1", "fs2")

save(output, file="SurvivalPriors.rdata")