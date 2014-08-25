setwd("/home/jc227089/Quoll_model/Random_outputs")

# vector of observations: popsize at Pob and Astell islands
obs<-c(766, 818, 470, 300, 3228, 4820, 2590, 2890)
lobs<-log(obs)
# function to standardize a vector (mean=0, sd=1)
stdz<-function(vec){
	(vec-mean(vec))/sd(vec)	
}
								
gold<-c() # place to put the best models
for (i in 1:200){
	fname=paste("ABCSample", i, ".RData", sep="")
	load(fname)
	out<-out[,!grepl("sr", colnames(out))]
	
	test<-apply(out[,7:14]==0, 1, sum)==0
	out<-out[test,]
	
	qtile<-0.005*5000/nrow(out) #work out quantile to sample
	if (is.na(qtile)==T) next
	
	temp<-(out[,7:14])
	
	temp<-sweep(temp, 2, obs) # calculate sum of squares
	temp<-temp^2 
	
	ss.pop.pb<-apply(temp[,1:4], 1, sum) 
	ss.pop.as<-apply(temp[,1:4], 1, sum) 
	ss.pop<-stdz(ss.pop.pb)+stdz(ss.pop.as) #sum standardised deviations
	
	quant.pop<-quantile(ss.pop, qtile) #find quantiles
	keepers<-ss.pop<=quant.pop
	gold<-rbind(gold, out[keepers,]) #keep these rows and rbind to gold
	rm(out)	
}

save(gold, file="../ABC/Kept_sims_poponlyhalfpcED.RData")