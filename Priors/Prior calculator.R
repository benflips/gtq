setwd("/Users/ellakelly/Documents/PhD/gtq/Priors")


#litter size
d<-read.table("Litter size.txt", header=T, sep="\t")
mean(d$mean)

# mean=6.725
# sd = 1.206140

#Male survival
d<-read.table("Male survival.txt", header=T, sep="\t")
p.ms<-d$surv/d$n
mean(p.ms)
sd(p.ms)
# mean = 0.04211957
# sd = 0.05893253

#Female survival year 1
d<-read.table("Female survival yr1.txt", header=T, sep="\t")
p.fs1<-d$surv/d$n
mean(p.fs1)
sd(p.fs1)
#mean= 0.2338186
#sd = 0.1235535

#Female survival year 2
d<-read.table("Female survival yr2.txt", header=T, sep="\t")
p.fs2<-d$surv/d$n
mean(p.fs2)
sd(p.fs2)
#mean= 0.02020202
#sd = 0.03499093