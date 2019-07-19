rm(list=ls())
install.packages("plyr"); install.packages("dplyr"); install.packages("foreign")
library(plyr); library(dplyr); library(foreign)

setwd(getwd())
# INPUT DATA 

numHousesPerVillage<-30 # number of households per village
numDates<-72 # number of experimental dates

# Sort data according to the experimental date
exptDayList<-seq(1:numDates)

# make shell for iterations below with house cordinates
# RANDOM HOUSE LOCATIONS with a range matching original study for calculation of distance differences
mosqDataShell<-array(-9,dim=c(30, 6))
mosqDataShell[,1]<-seq(1:numHousesPerVillage)
mosqDataShell[,2]<-runif(numHousesPerVillage,0.00,0.05) 
mosqDataShell[,3]<-runif(numHousesPerVillage,0.00,0.05) 
mosqDataShell[,4]<--9
mosqDataShell[,5]<--9
mosqDataShell[,6]<--9

mosqDataShell<-data.frame(mosqDataShell)
colnames(mosqDataShell)<-c("house","south","east","combo","exptDay","obsmosq")

# simulate data for all houses
mosqDataShell$notmissing<-1
mosqDataShell<-data.frame(mosqDataShell)

# sort data by house number and experimental day
mosqDataShell<-arrange(mosqDataShell, house, exptDay)

# specify inputs 
# SIMULATION FUNCTION FOR MOSQUITO DENSITIES AND COVERAGE SCENARIOS
mosqfunction = function(bbeta, lambda,  coverage, hseFactor){
# bbeta - proportion of mosquitoes repelled
# lambda - mean distance moved between households by the repelled mosquitoes
# coverage - number of households using repellents per day
# hseFactor - proportion of mosquitoes diverted elsewhere other than to a house
  
  loglik<-0 			
  NmosqPerDay<-1000 # total number of mosquitoes per day
  p<-rep(1/30, 30) # baseline distributions of mosquitoes in houses
  numDays<-72
  numHouses<-length(mosqDataShell[,1])
  
  # create an array for combo (presence of spatial repellent) per day
  comboArray<-array(0,dim=c(numHouses,numDays))

  # first 18 days have 0 coverage, next 54 days have some houses per day with combo==1
  for (t in 19:numDays) {
    selectedHouse<-sample(seq(1:numHouses),coverage,replace=FALSE)
    comboArray[selectedHouse,t]<-1
  }
  
  # for each experimental day, get predicted numbers of mosquitoes  
  for (t in 1:numDays) {
    
    mosqData1<-mosqDataShell
    mosqData1$combo<-comboArray[,t]
    mosqData1$exptDay<-t
    
    # set up a dataframe with every pair of houses 
    pairs<-array(-9,dim=c((30*30),(4*2)))
    pairs<-data.frame(pairs)
    #house1 = blocks of house1..., house2...
    for (h1 in 1:30) {
      for (h2 in 1:30) {
        pairs[(((h1-1)*30)+h2),1:4]<-mosqData1[h1,1:4]
      }
    }
    #house2 = house 1,2,3.., house 1,2,3... 
    for (h2 in 1:30) { 
      pairs[(1+(h2-1)*30):(30+(h2-1)*30),5:8]<-mosqData1[1:30,1:4]
    }
    colnames(pairs)<-c("house1","south1","east1","combo1","house2_", "south2_", "east2_", "combo2")
    # delete unnecessary pairs
    pairs<-pairs[pairs$house1!=pairs$house2,]
    
    # calculate harvesian distances
    x<-0.0174532925199433
    dlong<-x*(pairs$south2 - pairs$south1)
    dlat<-x*(pairs$east2 - pairs$east1)
    a<-((sin(dlat/2))^2) + cos(pairs$east1*x) * cos(pairs$east2*x) * ((sin(dlong/2))^2)
    pairs$distance<-6367 * (2 * atan2(sqrt(a), sqrt(1-a)))
    
    # get overall outgoing from house h1 and incoming to house h2 for calculations
    numPairs<-length(pairs[,1])
    for (d in 1:numPairs) {
      pairs$baseh1[d] <- p[pairs$house1[d]]
      pairs$outgoingh1[d] <- pairs$baseh1[d] * (bbeta*pairs$combo1[d]) 
      
      pairs$probGoToh2[d]<-exp(-(pairs$distance[d]^2)/(2*(lambda^2))) * (1-(pairs$combo2[d]))
        if (pairs$house2[d] == pairs$house1[d]) {pairs$probGoToh2[d]=0.0}
    }
    
    for (d in 1:numPairs) {
      # scale to one
      pairs$sumProbGoToh2[d]<-sum(pairs$probGoToh2[pairs$house1==pairs$house1[d]])
      pairs$incomingh2[d]<-0
      pairs$sumProbGoToh2[d][is.na(pairs$sumProbGoToh2[d])]<-0.00001 #to avoid missing values in denominator
      
      if (pairs$sumProbGoToh2[d]>0) {
        pairs$incomingh2[d] <- (pairs$probGoToh2[d]/pairs$sumProbGoToh2[d]) * (pairs$outgoingh1[d])
      }
    }
    
    # derive the total outgoing and incoming mosquitoes per house
    for (h in 1:30) {
      mosqData1$outgoing[h]<-mean(pairs$outgoingh1[pairs$house1==mosqData1$house[h]])
      mosqData1$incoming[h]<-sum(pairs$incomingh2[pairs$house2==mosqData1$house[h]])
      
    }

    # obtain the predicted mosquitoes in each house 
	mosqData1$predp<- (p - mosqData1$outgoing + (mosqData1$incoming * hseFactor)) * mosqData1$notmissing
    mosqData1$predp[is.na(mosqData1$predp)]<-0.00001 #to avoid missing values in denominator

	# obtain the adjusted predicted mosquitoes  
	mosqData1$adjpredp <- (mosqData1$predp/sum(mosqData1$predp))

	# simulate observed mosquitoes in each house from a multinomial distribution
    observedMosq<-rmultinom(1,NmosqPerDay,mosqData1$adjpredp)
    mosqData1$obsmosq<-observedMosq
	
	# save data by experimental day
    if (t==1) {
      mosqDataAll<-mosqData1
    }
    else if (t>1) {
      mosqDataAll<-rbind(mosqDataAll, mosqData1)
    }
    for (h in 1:30) {
      # ignore if adjpredprob is zero (missings)
     if (mosqData1$adjpredp[h]>0) {
      loglik <- loglik + (mosqData1$obsmosq[h]*log(mosqData1$adjpredp[h]))
    }
	 }
    
	# simulated data is mosqDataAll
    mosq<-mosqDataAll[,c("house", "east", "south", "exptDay", "obsmosq", "combo", "notmissing")]
  }
  return(mosq)
}

# small test run (since 100 takes a little time to run)
head(mosqfunction(0.2,0.5,24,0.8))
# example simulation for one combination of parameter values with 100 replicates
sim1<-data.frame(replicate(100,  mosqfunction(0.2, 0.5,  24, 0.8)))
# bbeta, lambda, coverage(as number of houses), hseFactor

for(i in 1:length(sim1)) {
  sim <- data.frame(sim1[,i])
  write.table(sim, file = paste0("agg1.1.05.24.2_", i,".txt"), row.names = F, col.names=F)
  sim0<-read.table(file = paste0("agg1.1.05.24.2_",  i,".txt"))
  colnames(sim0)<-c("house", "east", "south", "exptDay", "obsmosq", "combo", "notmissing")
  sim0<-arrange(sim0, exptDay,  house)
  write.table(sim0, file = paste0("agg1.1.05.24.2_",  i,".txt"), row.names = F, col.names=F)
}

