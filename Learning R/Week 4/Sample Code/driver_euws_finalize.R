#----------COMMENTS
#---R program to generate plots of OEP/AEP for poisson and clustered models and compare to historical data
#---code is specific to one country, plan is to launch multiple jobs in parallel using a script, inputing arguments
#---note that for efficiency, we may eventually want to split up some of the simulation computation across many processors
#---which can be done by bringing in the seed initialization by using an argument
#---CLEANED UP VERSION OF ORIGINAL driver.R FOR CREATION OF FINALIZED EWS model
#---NOTE THAT UNLIKE THE GENERATION OF THE SIMULATION FILES --- I AM NOT (PAINSTAKINGLY) TAKING OUT OPTIONS IN THE CODE - THESE ARE HARD CODED AND OK
#---THIS CODE IS FOR MODEL DEVELOPMENT ONLY and NOT INTENDED to be USED FOR FINALIZATION OF SIMULATION FILES
#----------












#----------CLEAR + SOURCE
rm(list=ls())

#---analytic solutions obtained from wang paper
source('wang_utilitiesK.R')

#---routines consistent with RL implementations
source('white_paper_routines.R')

#---plot routines
source('plots.R')
#----------











#----------PARAMETERS TO SPECIFY, INCLUDING ON/OFF SWITCHES

#---various on/off switches
LOAD.IN.HFREQ=FALSE
RUN.SIMS=TRUE
RUN.WANG=FALSE
optimal.beta.config=1 #dummy optimal configuration for eu windstorm
RUN.RL.CLUSTERS=FALSE
RUN.RL=FALSE
INTENSITY.FILTER=TRUE
NAHU=FALSE
EUWS=TRUE

#---following is a condition put in specifically to allow manipulation of mid-atlantic losses in NAHU clustering, and other versions of sub-clustering
MA.SPECIAL.CLUSTER=FALSE
INTENSITY.SPECIAL.CLUSTER=FALSE

data.table = 'ClusteredData_IED_STOCH_EUWS2011_HIST_ERA40INT_winter_test5v_3_split_w4_9'

#---file for the new EWS2011 model
eventid.rates = 'EUWS11.RATES.FINAL'


#---csv file where we load in a regions elt which contains the high freq events not in the clustered data
elt.with.high.freq = 'hfelt_'

#---number of grid points for integral/distribution approximations --- for analytic calculations only
ngrid = 2^12

#---new EWS2011 model --- noting the change of formatting
#---we note that these data were extracted from the database
#---note that this file does not have a header as we are specifying explicitly
ssi.filename='/home/skhare/clustering/get_euws2011_data/ssiews11.csv'; ssi=read.csv(ssi.filename,header=F)
ssi.eu=ssi[,2]
#----------









#----------READ COMMAND LINE ARGUMENTS, LOAD IN HISTORICAL LOSS DATA, CLUSTER DEFS, EVENT EXPOSURES ETC.
#---region gotten via command line argument  (commented out when testing)
args=commandArgs()
args = args[length(args)]
args<-strsplit(args,",")

#---set the variable which indicates the geographic region
#region <- as.character(args[[1]][1]); region; print('the region')
region='eu'

#---set the seed for the random number generation
#seed <- as.numeric(args[[1]][2]); seed; print('the seed'); set.seed(seed); seed.number = as.character(seed)
seed=7

#---set a file description of which cluster we are using
#description <- as.character(args[[1]][3])
description='test'

#---read ELT table with cluster definitions
filename = paste('data/',data.table,'.csv',sep='')
print(filename)
data = read.csv(filename,header=T)

#---get the pressures as an indicator of intensity --- vorticity in future models
pressure.header=paste('pressure')
print(pressure.header)
pressure=data[,pressure.header]


#---get event ids
eventids = data[,'event_id']

#---get losses
loss.header=paste('mean.loss.',region,sep='')
print(loss.header)
loss=data[,loss.header]

#---get the mid-atlantic losses, only if we are creating a special mid-atlantic cluster
if(MA.SPECIAL.CLUSTER){
loss.header=paste('mean.loss.','ma',sep='')
lossma=data[,loss.header]
}

#---if for european windstorm, get losses on target countries, france, germany, and uk
#if(EUWS){
#  loss.header=paste('mean.loss.','de',sep='')
#  lossde=data[,loss.header]
#  loss.header=paste('mean.loss.','fr',sep='')
#  lossfr=data[,loss.header]
#  loss.header=paste('mean.loss.','gb',sep='')
#  lossgb=data[,loss.header]
#}

#---get stdi
stdi.header=paste('stdi.loss.',region,sep='')
std.ind=data[,stdi.header]

#---get stdc
stdc.header=paste('stdc.loss.',region,sep='')
std.corr=data[,stdc.header]

#---get exposure
exposure.header=paste('exposure.',region,sep='')
expsr=data[,exposure.header]

#---get cluster designations
cluster.num=data[,'Cluster']
#---record the orignal cluster definitions in anticipation of upcoming manipulations
cluster.num.original=data[,'Cluster']





#---some filtering based on pressure/vorticity (need to ensure consistency when we define the cluster search space)
if(INTENSITY.FILTER){
  if(NAHU){
    print('intensity filtering on hurricane')
  }
  if(EUWS){

    
    #REMEBER CRUCIALLY THAT THE TAKE CONDITION IS SPECIFYING MAKING THIS POISSON
    #below is a statement about only clustering storms with appreciable intensity over countries of interest
    #dumping storms out of clusters which do not cause any appreciable loss in france, germany and the uk
    #take = which((lossde < 5e+8 ) & (lossgb < 1e+9 ) & (lossfr < 5e+8))
    #take = which((data[,"max.int.de"] < 30 | data[,"max.int.fr"] < 30 | data[,"max.int.gb"] < 30) )
    #take = which(data[,"max.int.wind"] < 37.5)
    #as of august 24th, 2010 --- we have a cutoff of 1.0 on the ssi for the old/current model and i believe that this was the
    #threshold used in the proposal results
    #take = which(ssi.eu < 4.5)
    #the following ssi cutoff was used for the old model
    ####take = which(ssi.eu < 1.0)
    #take = which((lossde < 5e+8 ) | (lossgb < 7e+8 ) | (lossfr < 2e+8))
    #take = which((lossde < 5e+8 | lossde > 1.5e+9) & (lossgb < 1e+9 | lossgb > 2e+9) & (lossfr < 5e+8 | lossfr > 1.5e+9))
    #take = which(pressure > 950)
    #take = which(data[,"max.int.country"]=='' & data[,"pressure"] > 950 )
    #take=which(data[,"pressure"] > 960)
    #take=which(data[,"mean.loss.eu"] < 5e+8)


    
    #LOGIC OF FOLLOWING CODE
    #FOR EACH CLUSTER, WE TAKE THE 80TH PERCENTILE IN SSI. EVERYTHING UPTO THE 80TH PERCENTILE IN SSI IS TREATED AS POISSON
    #SO EACH ORIGINAL CLUSTER, IN EFFECT, HAS A SUB-CLUSTER OF STRONG STORMS, WHICH WE MODULATE USING THE GAMMA DISTRIBUTION
    #THE IDEA IS THAT WHEN VIEWED AS A WHOLE, THE ENTIRE CLUSTER, WHICH IS AN ADMIXTURE OF POISSON AND NB, CAN BE VIEWED TOGETHER, AND WE CAN LOOK AT ITS OVERDISPERSION
    #AS A WHOLE --- THE TRICK OF ALL THIS IS THAT THE SUB-CLUSTERS ARE STRONG STORMS, WHICH WE CAN MANIPULATE QUITE STRONGLY, AND THUS SEE APPRECIABLE IMPACTS
    #ON THE TAIL STATISTICS
    
    num.base.clusters = length(unique(cluster.num))
    #dump out 3/4's of the clusters based entirely on ssi (percentage is tunable)
    alltake=c()
    for (i in 1:num.base.clusters){
      #convoluted logic to take 
      take = which(cluster.num == i)
      nnn = length(take)
      #note that for the final EWS implementation we have nnn*0.8 and if (threshold < 3.0){threshold = 3.0}
      threshold = sort(ssi.eu[take])[ceiling(nnn*0.80)]
      if (threshold < 3.0){threshold = 3.0}
      indices = which(cluster.num == i)
      ssis = ssi.eu[indices]
      takessi = which(ssis < threshold)
      take = indices[takessi]
      alltake=c(alltake,take)
    }
    take=alltake
    #dumping the unwanted clusters the 'cluster that we will treat as solely as poisson
    cluster.num[take] = num.base.clusters+1
    
  }
  
}


browser()





#---filter out zero loss events --- actually filtering out more severely to avoid getting eroneous beta distribution parameters
take = which(loss>10000)
loss=loss[take]; std.ind=std.ind[take]; std.corr=std.corr[take]; expsr=expsr[take]; cluster.num=cluster.num[take]; eventids=eventids[take]; pressure=pressure[take]; cluster.num.original=cluster.num.original[take];
#---also get the mid-atlantic losses
if(MA.SPECIAL.CLUSTER){lossma=lossma[take]}


#---get rates --- first we read in all possible event ids and rates, which crucially, ***are sorted in ascending eventid***
filename = paste('data/',eventid.rates,'.csv',sep='')
print(filename)
classes = c(rep('numeric',2))
eventids.and.rates = read.csv(filename,colClasses=classes,skip=0,header=F)

#---take the appropriate set of event ids and rates from the full set of event ids, with an error check
take=is.element(eventids.and.rates[,1],eventids)
take=which(take)
eventids.subset=eventids.and.rates[take,1]
rates.subset=eventids.and.rates[take,2]
if (length(eventids) == length(eventids.subset)){
  eventids=eventids.subset
  rate=rates.subset
  print('finished getting associated rates')
} else {
  print('!!!program failure, event ids for which we do not have a rate!!!')
  browser()
}


#---load in the high frequency event losses, and append to clustered data (no filtering necessary) - optional, only for EU data
if(LOAD.IN.HFREQ){
filename = paste('data/',elt.with.high.freq,region,'.csv',sep='')
print(filename)
classes = c(rep('numeric',2))
high.freq.elt = read.csv(filename,colClasses=classes,skip=0,header=F)
hfeventids = high.freq.elt[,1]
hflosses = high.freq.elt[,2]
hfstd.ind = high.freq.elt[,3]
hfstd.corr = high.freq.elt[,4]
hfexpsr = high.freq.elt[,5]
hfrate = high.freq.elt[,6]
n.hf.events = length(hfeventids)
eventids=c(eventids, hfeventids)
loss = c(loss, hflosses)
std.ind = c(std.ind, hfstd.ind)
std.corr = c(std.corr, hfstd.corr)
expsr = c(expsr, hfexpsr)
rate = c(rate, hfrate)
#---augment the cluster.num array to represent the hf events in EUWS, if we are doing intensity filtering, we throw them in with the previously defined poisson cluster
if(INTENSITY.FILTER){
  num.clusters = length(unique(cluster.num))
} else {
  num.clusters = length(unique(cluster.num))+1
}
add.high.freq.cluster = array(num.clusters,n.hf.events)
cluster.num = c(cluster.num, add.high.freq.cluster)
cluster.num.original = c(cluster.num.original,add.high.freq.cluster)
}
num.clusters = length(unique(cluster.num))


#---load in the historical data for this region



#---get mean loss ratio and standard deviations of the event loss, and also maximum loss
mlr = loss / expsr
std = (std.corr + std.ind) / expsr
#---get the maximum loss over the entire event set for grid discretization
MAXLOSS = max(loss)

#----------





#browser()



































#----------DEFINE WHICH CLUSTERS ARE POISSON or NB (HF EVENTS ARE TREATED AS POISSON), AND SET THE VARIANCES

if(FALSE){
#---code for our favourite cluster on the current european windstorm model
cluster.num = array(0,c(length(loss))); 
take = !loss > 5e8; cluster.num[take] = 1
take = loss > 5e8; cluster.num[take] = 2
num.beta.combinations = 1
num.clusters = length(unique(cluster.num))
poisson.or.nb <- array('dummy', c(num.clusters))
betas <- array(0,c(num.clusters,num.beta.combinations))
betas[] = 6
poisson.or.nb[1] = 'poisson' 
poisson.or.nb[2] = 'negative binomial'}


if(FALSE){
#---code for search space! naive european cluster (1,2,3) are cluster, the 4th weak storm cluster is treated as poisson
#---load in the beta search space
#---assign the poisson cluster a dummy variance of one
num.beta.combinations = 1000
betas = array(0,c(num.clusters,num.beta.combinations))  
betas[num.clusters,] = 1
ijk = 1
for (i in seq(10,200,by=20)){
  for (j in seq(10,200,by=20)){
    for (k in seq(10,200,by=20)){
      betas[1,ijk] = i; betas[2,ijk] = j; betas[3,ijk] = k
      ijk = ijk + 1
    }
  }
}
poisson.or.nb = array(0,c(num.clusters,num.beta.combinations))
poisson.or.nb[1] = 'negative binomial'
poisson.or.nb[2] = 'negative binomial'
poisson.or.nb[3] = 'negative binomial'
poisson.or.nb[4] = 'poisson'
}







#---next case is where we might consider doing some kind of sub-clustering to make MA come out correctly
#---note that this is a beta search space definition for US hurricane
if(FALSE){

#SUCCESSIVE CLUSTER DEFINITIONS WHICH ARE GETTING OVER-WRITTEN

  
#---code for search space for NAHU cluster defined in R:\Clustering.Simulations\NAHU2011\CLUSTERS\ClusteredData_IED_STOCH_HIST_test1v_4.csv
#---load in the beta search space
#---first cluster 
num.beta.combinations = 10^4
betas = array(0,c(num.clusters,num.beta.combinations))  
ijk = 1
for (i in seq(0.5,5,by=1)){
  for (j in seq(0.5,5,by=1)){
    for (k in seq(0.5,5,by=1)){
      for (n in seq(0.5,5,by=1)){
        betas[1,ijk] = i; betas[2,ijk] = j; betas[3,ijk] = k; betas[4,ijk] = n
        ijk = ijk + 1
      }
    }
  }
}
poisson.or.nb = array(0,c(num.clusters,num.beta.combinations))
poisson.or.nb[1] = 'negative binomial'
poisson.or.nb[2] = 'negative binomial'
poisson.or.nb[3] = 'negative binomial'
poisson.or.nb[4] = 'negative binomial'



#MA SPECIAL CLUSTER
#temp code to put in special clustering for the MA - essentially redefining what is above, and adding in a special MA cluster
if(MA.SPECIAL.CLUSTER){
  num.clusters = num.clusters+1
  num.beta.combinations=5^5
  betas = array(0,c(num.clusters, num.beta.combinations))
  #have a special subcluster of storms causing at least X=50000000 loss in MA
  take = which(lossma > 5e+6)
  cluster.num[take] = 5
  #recompute the beta parameter space
  ijk = 1
  for (i in seq(0.5,5,by=1)){
    for (j in seq(0.5,5,by=1)){
      for (k in seq(0.5,5,by=1)){
        for (n in seq(0.5,5,by=1)){
          for (m in seq(0.5,100,by=20)){
            betas[1,ijk] = i; betas[2,ijk] = j; betas[3,ijk] = k; betas[4,ijk] = n; betas[5,ijk] = m
              ijk = ijk + 1
          }
        }
      }
    }
  }
  poisson.or.nb = array(0,c(num.clusters,num.beta.combinations))
  poisson.or.nb[1] = 'negative binomial'
  poisson.or.nb[2] = 'negative binomial'
  poisson.or.nb[3] = 'negative binomial'
  poisson.or.nb[4] = 'negative binomial'
  poisson.or.nb[5] = 'negative binomial'
}

#NAHU --- note that this is the definition of the clusters i used in generating the 'proposal' results for NAHU
#the basic idea for the results i have sent is to start with the kossin-like definitions of the superclusters
#do the sub-clustering on intensity (which is done very simply with the landfalling pressure which correlates very well with loss)
#well, this isn't really sub-clustering, cause we could treat the low intensity storms using a negative binomial, but i didn't feel this
#was necessary. the idea was then to play around with the beta factors on the clusters 1 and 4, so that we end up matching the overall
#over-dispersion
if(INTENSITY.SPECIAL.CLUSTER){
  #---in each cluster, only accept storms beyond a certain intensity                                        
  num.clusters = num.clusters+1
  take = which(pressure > 950) 
  cluster.num[take] = 5
#  num.beta.combinations=10^4
  num.beta.combinations=1
  betas = array(0,c(num.clusters, num.beta.combinations))
  #---set the poisson cluster gamma distribution variance to dummy constant value
  betas[5,] = 1
  #---recompute the beta parameter space
  ijk = 1
  
#  for (i in seq(0.5,10,by=1)){
#    for (j in seq(0.5,10,by=1)){
#      for (k in seq(0.5,10,by=1)){
#        for (n in seq(0.5,10,by=1)){
#          betas[1,ijk] = i; betas[2,ijk] = j; betas[3,ijk] = k; betas[4,ijk] = n
#          ijk = ijk + 1
#        }
#      }
#    }
#  }
  betas[1,] = 1; betas[2,] = 1; betas[3,] = 1; betas[4,] = 1;
  
  poisson.or.nb = array(0,c(num.clusters,num.beta.combinations))
  poisson.or.nb[1] = 'negative binomial'
  poisson.or.nb[2] = 'poisson'
  poisson.or.nb[3] = 'poisson'
  poisson.or.nb[4] = 'negative binomial'
  poisson.or.nb[5] = 'poisson'
}

}









#4-cluster case EUROWIND
if(FALSE){
#---search space definition for ERA 40 mailier clusters
num.beta.combinations = 1
betas = array(0,c(num.clusters,num.beta.combinations))
#---set the betas of the poisson clusters to 1
betas[num.clusters,] = 1
betas[1,] = 1
betas[2,] = 1
betas[3,] = 1
#ijk = 1
#for (i in seq(0.5,60,by=2)){
#  for (j in seq(0.5,60,by=2)){
#      betas[2,ijk] = i; betas[3,ijk] = j
#      ijk = ijk + 1
#    }
#  }
poisson.or.nb = array(0,c(num.clusters,num.beta.combinations))
poisson.or.nb[1] = 'poisson'
poisson.or.nb[2] = 'negative binomial'
poisson.or.nb[3] = 'poisson'
poisson.or.nb[4] = 'poisson'
}







#THE FOLLOWING IS THE PARAMETERIZATION OF THE EWS CLUSTERS


if(TRUE){
#---search space definition for ERA40 (mailier type) clusters
num.beta.combinations = 1
betas = array(0,c(num.clusters,num.beta.combinations))
#---set the betas of the poisson clusters to 1
betas[num.clusters,] = 1
betas[1,] = 1
betas[2,] = 1
betas[3,] = 1
betas[4,] = 1
betas[5,] = 1
betas[6,] = 1
betas[7,] = 1
betas[8,] = 1
betas[9,] = 1
#ijk = 1
#for (i in seq(0.5,60,by=2)){
#  for (j in seq(0.5,60,by=2)){
#      betas[2,ijk] = i; betas[3,ijk] = j
#      ijk = ijk + 1
#    }
#  }
poisson.or.nb = array(0,c(num.clusters,num.beta.combinations))
poisson.or.nb[1] = 'negative binomial'
poisson.or.nb[2] = 'negative binomial'
poisson.or.nb[3] = 'poisson'
poisson.or.nb[4] = 'negative binomial'
poisson.or.nb[5] = 'negative binomial'
poisson.or.nb[6] = 'negative binomial'
poisson.or.nb[7] = 'negative binomial'
poisson.or.nb[8] = 'negative binomial'
poisson.or.nb[9] = 'negative binomial'
poisson.or.nb[10]= 'poisson'

}

#----------




































#----------COMMENTS
#---we now have:
#--->>> mlr, std, exposure, cluster.num, rate + cluster definitions <<<
#---which are all the arrays we require to compute the oep and aep
#----------









































#----------WANG ANALYTIC THEORY (!UPDATE FOR ZERO LOSS THRESHOLD!)
#---original.poisson.mixture = poisson.mixture( mlr, std, rate, expsr, loss, 5e8 , 6 )
if(RUN.WANG){
print('before analytic solution for K clusters')
print(date())
WANG = poisson.mixtureK( mlr, std, rate, expr, loss, num.clusters, cluster.num, betas, poisson.or.nb, ngrid )
print(date())
print('after analytic solution for K clusters')
#list k.clusters.analytic contains:
#xx(ngrid) - losses
#EEF(ngrid)
#OEPpoi(ngrid)
#OEP2poi(ngrid)
#OEP3poi(ngrid)
#AEPpoi(ngrid)
#OEP(num.beta.combinations, ngrid)
#AEP(num.beta.combinations, ngrid)
#---code to dump OEP curves to file for the different beta parameter combinations (and associated losses)
filename = paste('OEPclustered.', description, '.', region, '.csv', sep='')
write.csv(WANG$OEP, file=filename)
filename = paste('Losses.', description, '.', region, '.csv', sep='')
write.csv(WANG$xx, file=filename)
filename = paste('OEPpoisson.', description, '.', region, '.csv', sep='')
write.csv(WANG$OEPpoi, file=filename)
filename = paste('Betas.', description,'.', region,'.csv', sep='')
write.csv(betas, filename)
filename = paste('OverD.', description,'.', region,'.csv', sep='')
write.csv(WANG$overdispersion, file=filename)
filename = paste('RatesInClusters.', description,'.', region,'.csv', sep='')
write.csv(WANG$rates.in.clusters, file=filename)
filename = paste('AEPclustered.', description, '.', region, '.csv', sep='')
write.csv(WANG$AEP, file=filename)
filename = paste('AEPpoisson.', description, '.', region, '.csv', sep='')
write.csv(WANG$AEPpoi, file=filename)
#---record the aal due for each cluster, over this entire region
aal.bycluster = array(0,c(num.clusters))
for (i in 1:num.clusters){take = which(cluster.num == i); aal.bycluster[i] = sum( rate[take] * loss[take])}
filename = paste('AALbycluster.', description, '.', region, '.csv', sep='')
write.csv(aal.bycluster, file=filename)
}
#----------



















#for physical sub-clustering on the original euws case
#betas[1,]=16
#betas[2,]=1
#betas[3,]=1
#betas[4,]=4
#betas[5,]=4
#betas[6,]=4
#betas[7,]=2
#betas[8,]=1
#betas[9,]=4
#betas[10,]=1
#temp code for automated choice of over-dispersions for the 10 cluster case
#for EWS, the target over-dispersions come from analysis of the ERA 40 data
target.od=c(6.47/6.32,2.26/2.22,1,9.59/8.86,13.44/9.56,10.33/5.86,4.35/3.88,11.92/9.58,12.44/10.64,1)
if(TRUE){
  for(i in 1:num.clusters){
    take=which(cluster.num==i)
    rate.clusi = sum(rate[take])
    if (rate.clusi > 0){betas[i,]= (target.od[i] - 1)/sum(rate.clusi)}
    #note, for empty clusters, we set the beta to one --- instead of using the above formula which would generate an infinity
    if (rate.clusi == 0){betas[i,] = 1} #note that it does then matter what we set it to if the rate is truly zero
  }
}

#betas=betas*15 #####*3
#betas[5,]=betas[5,]*4
#betas[6,]=betas[6,]*.7
#betas[7,]=betas[7,]*1

#multiplication for the weak case with the new model
betas=betas*18
betas[4,] = betas[4,] * 0.4
betas[5,] = betas[5,] * 0.5
betas[7,] = betas[7] * 1.1
betas[8,] = betas[8,] * 0.20
betas[9,] = betas[9,]*1.4
betas[6,] = betas[6,] * 0.60


#---permanently fix the EU tuned betas (these yield acceptable results as of feb 14th, 2011)
#---FOLLOWING GAMMA DISTRIBUTIONS GENERATE SETS OF RESULTS 
#betas[1]  = 153.396622
#betas[2,] = 11.366536
#betas[3,] =  0.000000
#betas[4,] = 2.858247
#betas[5,] = 41.169611
#betas[6,] = 573.982675
#betas[7,] = 313.943917
#betas[8,] = 10.196645
#betas[9,] = 102.287620
#betas[10,] = 0.000000



#weak case #1 --- as of february 21st, 2011
#betas[1,] = 47.936
#betas[2,] = 3.552
#betas[3,] = 0.00
#betas[4,] = 0.714
#betas[5,] = 12.8655
#betas[6,] = 107.62
#betas[7,] = 98.107
#betas[8,] = 2.124
#betas[9,] = 31.96
#betas[10,] = 0.00




#weak case #2 --- as of february 21st, 2011
#THIS HAS BEEN SIGNED OFF BY ROBERT-MUIR-WOOD
betas[1,] = 86.3
betas[2,] = 6.393
betas[3,] = 0.00
betas[4,] = 1.286
betas[5,] = 23.157
betas[6,] = 193.719
betas[7,] = 176.59
betas[8,] = 3.823442
betas[9,] = 57.536786
betas[10,] = 0.00

















#---JUNE 1ST, 2010 --- MAKING MAJOR CODE CHANGE TO RUN SIMULATIONS IN MEAN MODE






#----------CODE TO GENERATE VARIOUS SIMULATIONS---***some early versions of simulation codes copy and pasted to bottom of this file*** 
#---this version is optimized for speed, *except* we still need to loop over the number of years, and sample from the events, as it does not appear
#---that we can insert arrays of probabilities into the R sample function
if(RUN.SIMS){

  #---set the number of years in this simulation
  nyears=10^4
  #---set parameter which selects out the configuration of beta parameters we'll want to look at
  beta.config=optimal.beta.config
  EPs = array(0, c(nyears)); for(i in 1:nyears){EPs[i] = 1 - i/nyears}

  #---for the K cluster model
  gammas = array(0,c(num.clusters,nyears))
  for (i in 1:num.clusters){
    if (poisson.or.nb[i] == 'poisson')          { gammas[i,] = 1 }
    if (poisson.or.nb[i] == 'negative binomial'){ gammas[i,]  = rgamma(nyears, shape=1/betas[i,beta.config], scale=betas[i,beta.config])}
  }
  rates.of.clusters = array(0,c(num.clusters,nyears))
  for (i in 1:num.clusters){ rates.of.clusters[i,] = sum( rate[which(cluster.num==i)] ) }
  for (i in 1:num.clusters){ rates.of.clusters[i,] = rates.of.clusters[i,] * gammas[i,]  }
  num.events.in.clusters = array(0,c(num.clusters,nyears))
  for (i in 1:num.clusters){ num.events.in.clusters[i,] = rpois(nyears,rates.of.clusters[i,]) }



  
  #---compute the time series of events, corresponding to the orignal cluster defs, but maintaining sub-clustering
  rates.of.clusters.original = array(0,c(num.clusters,nyears))
  for (i in 1:num.clusters){
    #---for the case where we have the original cluster defined to be poisson
    if (poisson.or.nb[i] == 'poisson'){rates.of.clusters.original[i,] = sum(rate[which(cluster.num.original==i)])}
    #---for the case where we take the original cluster, but do some negative binomial sub-clustering
    if (poisson.or.nb[i] == 'negative binomial'){
      #get the positions associated with the original cluster
      take=which(cluster.num.original==i)
      all.rates.i=rate[take]
      all.sub.cluster.num = cluster.num[take]
      take=which(all.sub.cluster.num==i)
      #we treat the ones that we are sub-clustering as negative binomial, drawn from the prescribed gamma
      rates.of.clusters.original[i,] = sum(all.rates.i[take]) * gammas[i,]
      take=which(all.sub.cluster.num != i)
      #the events do not remain in this cluster i after sub-clustering are treated as poisson
      #my rule is that if we have an original cluster, some of which gets sub-clustered, the remaining events are poisson
      rates.of.clusters.original[i,] = rates.of.clusters.original[i,] + sum(all.rates.i[take])
    }
  }
  num.events.in.clusters.original = array(0, c(num.clusters,nyears))
  for (i in 1:num.clusters){ num.events.in.clusters.original[i,] = rpois(nyears,rates.of.clusters.original[i,]) }

  
  all.losses.cluster = NULL
  for (i in 1:num.clusters){

    take = which(cluster.num == i)
    ids = eventids[take]
    mlrs = mlr[take]
    stds = std[take]
    rate.base = rate[take]
    exps = expsr[take]

    all.ids.cluster.i = NULL
    #can we somehow pass arrays of prob's to the sample function so we need not loop over the number of years
    for (j in 1:nyears){
      if(num.events.in.clusters[i,j] > 0){all.ids.cluster.i = c(all.ids.cluster.i, sample(ids, num.events.in.clusters[i,j], replace = TRUE, prob = rate.base * gammas[i,j])) }
    }

    take = match(all.ids.cluster.i,ids)
    means = mlrs[take]
    sigmas = stds[take]
    exposures = exps[take]
    beta.coeffs = get.alpha.beta.vec(means, sigmas)
    
    all.non.zero.losses = rbeta( sum(num.events.in.clusters[i,]), shape1 = beta.coeffs[,1], shape2 = beta.coeffs[,2] ) * exposures

    
    all.losses.cluster.i = array(0, c(nyears, max(num.events.in.clusters[i,])) )
    endj = 0
    for (j in 1:nyears){
      if(num.events.in.clusters[i,j] > 0){
        endj = endj + num.events.in.clusters[i,j]
        startj = endj - num.events.in.clusters[i,j] + 1
        all.losses.cluster.i[j,1:num.events.in.clusters[i,j]] =  all.non.zero.losses[startj:endj]
      }
    }
    all.losses.cluster = cbind(all.losses.cluster, all.losses.cluster.i)
  }








   

  #---ensure the losses in each year are sorted in descending order
  all.losses.cluster = t(apply(all.losses.cluster,1,sort))
  num.cols = dim(all.losses.cluster)[2]
  all.losses.cluster = all.losses.cluster[,num.cols:1]
  AGG = apply(all.losses.cluster,1,sum)
  OMX = sort(all.losses.cluster[,1])
  OE2 = sort(all.losses.cluster[,2])
  OE3 = sort(all.losses.cluster[,3])
  OE4 = sort(all.losses.cluster[,4])
  
  #---for the usual Poisson model (agrees with analytic solution)
  num.poisson.events = rpois(nyears,sum(rate))
  total.num.poisson.events = sum(num.poisson.events)
  all.ids = sample( eventids, total.num.poisson.events, replace = TRUE, prob = rate )
  take = match(all.ids,eventids)
  means = mlr[take]
  sigmas = std[take]
  exposures = expsr[take]
  beta.coeffs = get.alpha.beta.vec(means, sigmas)
  all.non.zero.losses = rbeta( total.num.poisson.events, shape1 = beta.coeffs[,1], shape2 = beta.coeffs[,2] ) * exposures

  
  all.losses.poisson = array(0,c(nyears,max(num.poisson.events)))
  endj = 0
  for (j in 1:nyears){
    if(num.poisson.events[j] > 0){
      endj = endj + num.poisson.events[j]
      startj = endj - num.poisson.events[j] + 1
      all.losses.poisson[j,1:num.poisson.events[j]] = sort( all.non.zero.losses[startj:endj], decreasing = TRUE)
    }
  }
  AGGP = apply(all.losses.poisson,1,sum)
  OMXP = sort(all.losses.poisson[,1])
  OE2P = sort(all.losses.poisson[,2])
  OE3P = sort(all.losses.poisson[,3])
  OE4P = sort(all.losses.poisson[,4])

  #save the clustered and poisson assumption results
  filename = paste('ALossesclustered.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(AGG, file=filename)
  filename = paste('OLossesclustered.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(OMX, file=filename)
  filename = paste('OE22222clustered.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(OE2, file=filename)
  filename = paste('OE33333clustered.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(OE3, file=filename)
  filename = paste('OE44444clustered.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(OE4, file=filename)
  filename = paste('ALossespoissonnn.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(AGGP, file=filename)
  filename = paste('OLossespoissonnn.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(OMXP, file=filename)
  filename = paste('OE22222poissonnn.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(OE2P, file=filename)
  filename = paste('OE33333poissonnn.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(OE3P, file=filename)
  filename = paste('OE44444poissonnn.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(OE4P, file=filename)
  filename = paste('numofstormssssss.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(t(num.events.in.clusters), file=filename)
  filename = paste('numofstormspoiss.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(num.poisson.events, file=filename)
  filename = paste('allcluslossessss.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(all.losses.cluster, file=filename)
  filename = paste('allpoislossessss.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(all.losses.poisson, file=filename)
  filename = paste('numofstormsorign.', description, '.', region, '.seed', seed, '.csv', sep='')
  write.csv(t(num.events.in.clusters.original), file=filename)
}
#----------







































































#----------GET EQUIVALENT CALCS TO WHAT IS IN RL
if(RUN.RL){
print('before RL poisson calculations')
print(date())
RL8.0.Poisson = wp_oep_aep  ( rate, loss, expsr, std, mlr, 'poisson', NA, MAXLOSS, ngrid )
print(date())
print('after the RL poisson calculations')
}
#list RL8.0.Poisson contains:
#xx(ngrid) - losses, and does not include the zero loss level
#OEP(ngrid + 1) - includes zero loss level
#AEP(ngrid + 1) - includes zero loss level
#----------
#----------GET RESULTS USING CONVOLUTIONS AS IN RL FOR K CLUSTERS (!UPDATE FOR ZERO LOSS THRESHOLD!)
#---notes about the RL implementation: note that we get exactly the same OEP, but different
#---AEPs when we go to the wang implementation, especially for high rate events.
#---is the difference in AEP due to some truncation in the convolution? or the fact that
#---we're not summing to an infinite number events. does the RL implementation allow for
#---expanded convolutions? if we get time, we should understand a bit more about what may
#---be wrong with the implementation of the AEP in RL
if(RUN.RL.CLUSTERS){
print('before the RL implementation for K clusters')
print(date())
#---first loop over clusters to get cluster CEPs and convolutions, which need not be repeated
CEP.clusters = array(0,c(num.clusters,ngrid))
MAXEVENTS.clusters = array(0,c(num.clusters))
rates.clusters = array(0,c(num.clusters))
#initialize collection of aggregate density matrix with null
dAEP.clusters = c()
for (j in 1:num.clusters){
  take = which(cluster.num==j)
  cep.daep = get.cep.daep(rate[take], loss[take], expsr[take], std[take], mlr[take], MAXLOSS, ngrid, betas[j,i])
  rates.clusters[j] = sum(rate[take])
  #list cep.daep contains:
  #xx(ngrid) - losses (!same for each cluster, since input MAXLOSS is the same for each cluster!)
  #CEP(ngrid) 
  #dAEP(MAXEVENTS,ngrid)
  #MAXEVENTS
  CEP.clusters[j,] = cep.daep$CEP
  dAEP.clusters=rbind(dAEP.clusters, cep.daep$dAEP)
  MAXEVENTS.clusters[j] = cep.daep$MAXEVENTS
}
RL.clusters = get.RL.clusters( CEP.clusters, dAEP.clusters, MAXEVENTS.clusters, num.beta.combinations, num.clusters, poisson.or.nb, betas, ngrid, rates.clusters, cep.daep$xx)
#list RL.clusters contains:
#xx(ngrid) - losses and does not include the zero loss level
#OEP(ngrid + 1) - includes the zero loss level
#AEP(ngrid + 1) - includes the zero loss level
print(date())
print('after the RL implemenation for K clusters')
}
#----------






















































































































#---some junk/commented out code below
#---simulation gotten from a pure poisson process, and outputs a year event table
#if (FALSE){
#nyears=10^5
#nevents.per.year = rpois(nyears,sum(rate))
#table.events = array(0,c(nyears,max(nevents.per.year)))
#table.losses = array(0,c(nyears,max(nevents.per.year)))
#for (i in 1:nyears){
#   npicks = nevents.per.year[i]
#   the.events = sample(eventids, npicks, replace=TRUE, prob = rate)
#   table.events[i,1:npicks] = the.events
#   take = which(is.element(eventids,the.events))
#   table.losses[i,1:length(take)] = loss[take]
#}
#write.csv(table.events, file = "euws.ied.with.pla.csv")
#aal = sum(table.losses/nyears)
#}

#---naive simulation code to verify analytic theories and compute OEPN curves - for cluster and poisson
#if (FALSE){
#  nyears=10^4
#  beta.config=1
#  gammas=array(1,c(num.clusters))
#  AGG <- AGGP <- OMX <- OMXP <- OE2 <- OE2P <- OE3 <- OE3P <- array(0,c(nyears))
#  gamma.check = array(0,c(num.clusters,nyears))
#  EPs = array(0, c(nyears))
#  for (i in 1:nyears){
#    #---get random losses for the conditionally Poisson cluster, and the Poisson assumption as well
#    ratei = rate
#    for (j in 1:num.clusters){
#      if (poisson.or.nb[j] == 'poisson') {gammas[j] = 1}
#      if (poisson.or.nb[j] == 'negative binomial') {gammas[j] = rgamma(1,shape=1/betas[j,beta.config],scale=betas[j,beta.config]); take = which(cluster.num == j); ratei[take] = gammas[j] * ratei[take]; gamma.check[j,i] = gammas[j] }
#    }
#    nevents.yeari  = rpois(1,sum(ratei))
#    nevents.yeariP = rpois(1,sum(rate))
#    the.eventsi  = sample(eventids, nevents.yeari, replace = TRUE, prob = ratei)
#    the.eventsiP = sample(eventids, nevents.yeariP, replace = TRUE, prob = rate)
#    take  = match(the.eventsi ,eventids)
#    takeP = match(the.eventsiP,eventids)
#    std.losses = std[take]
#    std.lossesP = std[takeP]  
#    mlrs = mlr[take]
#    mlrsP = mlr[takeP]
#    exposures.i = expsr[take]
#    exposures.iP = expsr[takeP]
#    the.random.losses.i = array(0, c(length(mlrs)) )
#    the.random.losses.iP = array(0, c(length(mlrsP)) )
#    #---cover the zero event case
#    if (nevents.yeari  == 0){AGG[i]  <- OMX[i]  <- OE2[i]  <- OE3[i]  <- 0}
#    if (nevents.yeariP == 0){AGGP[i] <- OMXP[i] <- OE2P[i] <- OE3P[i] <- 0}
#    #---cover the case where we have at least one event - clusters
#    if (nevents.yeari > 0){
#      for (ijk in 1:nevents.yeari) { beta.coeffs=get.alpha.beta(mlrs[ijk],std.losses[ijk]); the.random.losses.i[ijk] = rbeta(1, shape1 = beta.coeffs[1], shape2 = beta.coeffs[2]) * exposures.i[ijk] }
#      RLC = sort(the.random.losses.i, decreasing=TRUE)
#      if (nevents.yeari == 1){ AGG[i] <- OMX[i] <- RLC[1]; OE2[i] <- OE3[i] <- 0 }
#      if (nevents.yeari == 2){ AGG[i] = sum(RLC); OMX[i] = RLC[1]; OE2[i] = RLC[2]; OE3[i] = 0 }
#      if (nevents.yeari >  2){ AGG[i] = sum(RLC); OMX[i] = RLC[1]; OE2[i] = RLC[2]; OE3[i] = RLC[3] }
#    }
#    #---cover the case where we have at least one event - Poisson
#    if (nevents.yeariP > 0){
#      for (ijk in 1:nevents.yeariP) { beta.coeffs=get.alpha.beta(mlrsP[ijk],std.lossesP[ijk]); the.random.losses.iP[ijk] = rbeta(1, shape1 = beta.coeffs[1], shape2 = beta.coeffs[2]) * exposures.iP[ijk] }
#      RLP = sort(the.random.losses.iP, decreasing=TRUE)
#      if (nevents.yeariP == 1) { AGGP[i] <- OMXP[i] <- RLP[1]; OE2P[i] <- OE3P[i] <- 0 }
#      if (nevents.yeariP == 2) { AGGP[i] = sum(RLP); OMXP[i] = RLP[1]; OE2P[i] = RLP[2]; OE3P[i] = 0 }
#      if (nevents.yeariP >  2) { AGGP[i] = sum(RLP); OMXP[i] = RLP[1]; OE2P[i] = RLP[2]; OE3P[i] = RLP[3] }
#    }
#    EPs[i] = 1 - i/nyears
#  }
#  #---sort all the simulated losses, in ascending order
#  #---clustered simulations
#  AGG  = sort(AGG)
#  OMX  = sort(OMX)
#  OE2  = sort(OE2)
#  OE3  = sort(OE3)
#  #---poisson assumption simulations
#  AGGP = sort(AGGP)
#  OMXP = sort(OMXP)
#  OE2P = sort(OE2P)
#  OE3P = sort(OE3P)
#}
