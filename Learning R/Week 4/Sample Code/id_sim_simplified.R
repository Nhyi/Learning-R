







set.seed(1)


##########
# simplified simulation routine for the high frequency poisson case
##########
get.id.sim.HF <- function( ratek, idsk, Nyears ){


  ##########
  # determine number of years with nevents given our choice of frequency distribution
  ##########
  #---get average annual rate applicable across frequency distributions
  lamda = sum(as.numeric(ratek))
  
  #---lamda above specifies the one parameter poisson distribution
  nevents.in.sim = NULL
  ni.array = NULL
  ni = 0
  while (ni >= 0){


      #---get the probability of ni events from the poisson distribution and multiply by Nyears
      #---good aspect of R is that there is nothing wrong with analytic evaluation numerically
      #nyears.with.ni = Nyears * exp(-lamda) * lamda^ni / factorial( ni )
      nyears.with.ni = Nyears * dpois(ni, lambda = lamda)
      
      rounded.nyears.with.ni = round(nyears.with.ni)
      
      if (rounded.nyears.with.ni >= 1){nevents.in.sim = rbind(nevents.in.sim, rounded.nyears.with.ni);
                                       ni.array = c(ni.array,ni);
                                       ni = ni + 1}



      print(ni); print(rounded.nyears.with.ni)
      #browser()
      # if we are on the left side of the distribution, we want to keep going, make sure ni>5 so we keep going
      # these manipulations are inexact under some circumstances and used with caution in testing phase
      if (rounded.nyears.with.ni < 1) {
        if (ni > 66){
          ni = -99
        }else{
          ni = ni + 1
        }
      }

      
  }
  #--if we have not filled up all the years then we increase the number of years with zero losses as an approximation
  if (sum(nevents.in.sim) < Nyears){nevents.in.sim[1] = nevents.in.sim[1] + (Nyears - sum(nevents.in.sim))}


  

  ##########
  # stratified year, event id and genesis date simulation
  ##########
  

  #---get the number of events per year over the Nyears
  num.per.year = NULL
  for (i in 1:length(ni.array)){
#    print(i)
    num.per.year = c( num.per.year, array(ni.array[i],c(nevents.in.sim[i])) )
  }



#  browser()

  #---randomly shuffle the time series over the Nyears
  random.years = sample(num.per.year,Nyears,replace=FALSE,prob=NULL)

  #---get the list of years with non-zero events, and repeat year number M times if there are M events in a particular year
  allyears = seq(1,Nyears,by=1)
  take = which(random.years > 0)
  allyears = allyears[take]
  num.events.non.zero.years = random.years[take]
  #---create an array holding the list of years
  yearssim = array( 0, c(sum(num.events.non.zero.years)) )
  #---get the mother list of years
  dumi = 1
  for (i in 1:length(allyears)){
    yearssim[dumi:(dumi+num.events.non.zero.years[i]-1)] = allyears[i]
    dumi = dumi + num.events.non.zero.years[i]
  }
  
  #---get representative sample of the conditional event distribution (not sorting in terms of loss as unnecessary)
  num.samples.ced = sum(num.events.non.zero.years)
  #---partition 0-->>1 with midpoints
  right.partition  = seq(1,  num.samples.ced,    by=1) / num.samples.ced
  left.partition = seq(0, (num.samples.ced-1), by=1) / num.samples.ced
  midpoints.partition = (left.partition + right.partition)/2.0
  
  #---now getting the cummulative event distribution
  event.probs = ratek / sum( ratek )
  cumsum.events = cumsum(event.probs)

  
  #---sweep over all possible events (picking events consistent with our partitioning of the uniform distribution)
  #---adding in simulation of base genesis dates
  idsimulation = array(0,c(num.samples.ced))
#  genesis.simulation = array("0",c(num.samples.ced))
  dumi = 1
  for (i in 1:length(ratek)){
#    print(i)
    take = which( midpoints.partition <= cumsum.events[i] )
    #---if there are any events to take, add to the idsimulation list, and then ignore by setting to value greater than 1
    if(length(take) > 0){idsimulation[dumi:(dumi + length(take) - 1)] = idsk[i]
                      #   genesis.simulation[dumi:(dumi + length(take) - 1)] = as.character(as.Date(datesk[i],"%d/%m/%y"))
                         midpoints.partition[take] = 99.0}
    dumi = dumi + length(take)
  }
  
  #---spread the sample ids across the years randomly (keeping track of the associated genesis dates by using the match function)
  random.ids = sample(idsimulation, num.samples.ced, replace=FALSE, prob=NULL)
  take = match(random.ids, idsimulation)
 # associated.genesis.dates = genesis.simulation[take]

  print('completed event id simulation')


  ##########
  # stratified sampling of associated loss quantiles 
  ##########
  
  #---get the simulation of the loss quantiles (obtained from random shuffling of the midpoints.partition used for the conditional event distribution)
  #---using left/right partitions stored in memory
  allsortedquantiles = (right.partition + left.partition)/2.0
  randomquantiles = sample(allsortedquantiles, num.samples.ced, replace=FALSE, prob=NULL)


  print('completed stratified sampling of the loss quantiles')

  


  ##########
  # return stratified years simulation, event/source ids, random loss quantiles, genesis dates and randomized dates
  ##########

  #---return the list of simulation year, ids, and loss quantiles
  return(list( yearssim=yearssim, random.ids=random.ids, randomquantiles=randomquantiles ))
  
  
  
}




