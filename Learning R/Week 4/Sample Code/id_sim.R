##########
# routine to generate Nyears simulations of event ids based on a poisson or negative binomial distribution
# used for any sub-event set whether negative binomial (for overall rate) or poisson
# note that for years without events, nothing is output
##########
get.id.sim.SCk <- function( ratek, idsk, datesk, Nyears, poisson.or.nb, var.gamma ){


  ##########
  # determine number of years with nevents given our choice of frequency distribution
  ##########

  #---get average annual rate applicable across frequency distributions
  lamda = sum(as.numeric(ratek))

  if ( poisson.or.nb == "poisson" ){

    #---lamda above specifies the one parameter poisson distribution
    nevents.in.sim = NULL
    ni.array = NULL
    ni = 0
    while (ni >= 0){

      #---get the probability of ni events from the poisson distribution and multiply by Nyears
      #---good aspect of R is that there is nothing wrong with analytic evaluation numerically
      nyears.with.ni = Nyears * exp(-lamda) * lamda^ni / factorial( ni )
      # nyears.with.ni = Nyears * dpois(ni, lambda = lamda)
      
      rounded.nyears.with.ni = round(nyears.with.ni)
      
      if (rounded.nyears.with.ni >= 1){nevents.in.sim = rbind(nevents.in.sim, rounded.nyears.with.ni);
                                       ni.array = c(ni.array,ni);
                                       ni = ni + 1}

      # if we are on the left side of the distribution, we want to keep going, make sure ni>5 so we keep going
      # these manipulations are inexact under some circumstances and used with caution in testing phase
      # the final test of whether or not this is technically correct will be the sim platform results
      if (rounded.nyears.with.ni < 1) {
        if (ni > 5){
          ni = -99
        }else{
          ni = ni + 1
        }
      }

      
    }
    if (sum(nevents.in.sim) < Nyears){nevents.in.sim[1] = nevents.in.sim[1] + (Nyears - sum(nevents.in.sim))}

    
  }




  
  
  if ( poisson.or.nb == "negative binomial" ){
    
    #---frequency distribution is negative binomial with mean = lamda and variance = lamda + var.gamm*lamda^2
    #---using definition 6.3 in the Klugman loss model book ... using r and B (beta) definitions
    r = 1.0/var.gamma
    B = var.gamma * lamda
    
    #---compute num of years with 0,1,2, ... storms until rounded number of years to nearest integer is 0
    #---init an nevents array which stores number of years with 0,1,2,3,...,MAX number of storms
    #---start by computing number of years with zero events
    #---recovering missing years by adding in a few zero years
    #---i don't believe that our use of the analytic formula is leading to any numerical difficulty due to op order
    #---note that for the ews clustered cases, i have not noticed any problem of this code 'stopping' on the left hand side of the distribution
    nevents.in.sim = NULL
    ni.array = NULL 
    ni = 0
    while (ni >= 0){
      
      nyears.with.ni = Nyears * ( choose( (ni+r-1), ni )  * (1/(1+B))^r * (B/(1+B))^ni )

      rounded.nyears.with.ni = round(nyears.with.ni)
      
      if (rounded.nyears.with.ni >= 1){nevents.in.sim = rbind(nevents.in.sim, rounded.nyears.with.ni);
                                       ni.array = c(ni.array,ni);
                                       ni = ni + 1}

      if (rounded.nyears.with.ni < 1) {ni = -99}

      print(rounded.nyears.with.ni)

    }
    if (sum(nevents.in.sim) < Nyears){nevents.in.sim[1] = nevents.in.sim[1] + (Nyears - sum(nevents.in.sim))}
  }


  

  

  ##########
  # stratified year, event id and genesis date simulation
  ##########
  

  #---get the number of events per year over the Nyears
  num.per.year = NULL
  for (i in 1:length(ni.array)){
#    print(i)
    num.per.year = c( num.per.year, array(ni.array[i],c(nevents.in.sim[i])) )
  }



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
  genesis.simulation = array("0",c(num.samples.ced))
  dumi = 1
  for (i in 1:length(ratek)){
#    print(i)
    take = which( midpoints.partition <= cumsum.events[i] )
    #---if there are any events to take, add to the idsimulation list, and then ignore by setting to value greater than 1
    if(length(take) > 0){idsimulation[dumi:(dumi + length(take) - 1)] = idsk[i]
                         genesis.simulation[dumi:(dumi + length(take) - 1)] = as.character(as.Date(datesk[i],"%d/%m/%y"))
                         midpoints.partition[take] = 99.0}
    dumi = dumi + length(take)
  }
  
  #---spread the sample ids across the years randomly (keeping track of the associated genesis dates by using the match function)
  random.ids = sample(idsimulation, num.samples.ced, replace=FALSE, prob=NULL)
  take = match(random.ids, idsimulation)
  associated.genesis.dates = genesis.simulation[take]

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
  # randomization of genesis dates by monthly bin
  ##########
  

  #---randomize the associated genesis dates within monthly bins (adding +- 15 days) by doing a stratified sampling
  #---filling remainder days with zero increments. note that our base year is 2011 and some increments will push us into 2010 and 2012
  #---changing back to year 2011 will be done at a script level
  min.day = -15
  max.day = 15
  num.days = 1 - min.day + max.day
  num.each.increment = floor(num.samples.ced/num.days)
  num.missing = num.samples.ced - num.each.increment*num.days
  plus.minus.days = array(-99,c(num.samples.ced))
  plus.minus.days[1:num.missing] = 0
  dumi = num.missing+1
  for (i in seq(-15,15,by=1)){
    plus.minus.days[dumi:(dumi+num.each.increment-1)] = i
    dumi = dumi + num.each.increment
  }
  random.day.increments = sample(plus.minus.days)
  ran.dates = as.Date(associated.genesis.dates) + random.day.increments

  print('completed date randomization')



  ##########
  # return stratified years simulation, event/source ids, random loss quantiles, genesis dates and randomized dates
  ##########

  #---return the list of simulation year, ids, and loss quantiles
  return(list( yearssim=yearssim, random.ids=random.ids, randomquantiles=randomquantiles, associated.genesis.dates=associated.genesis.dates, ran.dates=ran.dates ))
  
  
  
}



















##########
# routine to generate ep statistics from simulation array (which includes randomized losses)
# ***assumes that simyear array is in ascending order***
# ***need to be very careful when filling in the agg.loss array to fill in the appropriate years***
# ***above will not matter for the overall statistics, but is crucial when we are considering ordering***
##########
get.AEPOEP.from.simarray <- function ( simarray, Nyears ){

  #---get simyear and lossses
  simyear    = simarray[,1]
  losses     = as.numeric(simarray[,5])

  #---define aggregate and occurrence loss arrays over all years and initializing to zero
  agg.loss = array(0,c(Nyears))
  occ.loss = array(0,c(Nyears))

  #---get starting indeces of unique years ... this will result in skipping 'missing' years
  duplicates    = duplicated(simyear)
  start.indeces = which(duplicates == FALSE)
  unique.years = unique(simyear)
  num.unique.years = length(unique.years)

  #---loop over all years except the last (filling in correct years!)
  for (i in 1:(num.unique.years-1)){
    #---note again that these year.to.fill's will on occassion skip years when there are not events, but those years are zeroed at the start so ok
    year.to.fill = unique.years[i]
    #---in what follows we are taking the data that belongs to the year of interest
    agg.loss[year.to.fill] = sum(losses[start.indeces[i]:(start.indeces[i+1]-1)])
    occ.loss[year.to.fill] = max(losses[start.indeces[i]:(start.indeces[i+1]-1)])
  }

  #---get the last year inefficiently but only once
  #---we are getting the last year for which there are simulation data available
  final.year = unique.years[ num.unique.years ]
  take = which( simyear == final.year )
  agg.loss[ final.year ] = sum( losses[take] )
  occ.loss[ final.year ] = max( losses[take] )

  return(list( occ.loss=occ.loss, agg.loss=agg.loss ))
  
}






get.AEPOEP.HF <- function ( year, sim.losses, Nyears ){

  #---get simyear and lossses
  simyear    = year
  losses     = sim.losses

  #---define aggregate and occurrence loss arrays over all years and initializing to zero
  agg.loss = array(0,c(Nyears))
  occ.loss = array(0,c(Nyears))

  #---get starting indeces of unique years ... this will result in skipping 'missing' years
  duplicates    = duplicated(simyear)
  start.indeces = which(duplicates == FALSE)
  unique.years = unique(simyear)
  num.unique.years = length(unique.years)

  #---loop over all years except the last (filling in correct years!)
  for (i in 1:(num.unique.years-1)){
    #---note again that these year.to.fill's will on occassion skip years when there are not events, but those years are zeroed at the start so ok
    year.to.fill = unique.years[i]
    #---in what follows we are taking the data that belongs to the year of interest
    agg.loss[year.to.fill] = sum(losses[start.indeces[i]:(start.indeces[i+1]-1)])
    occ.loss[year.to.fill] = max(losses[start.indeces[i]:(start.indeces[i+1]-1)])
    print(i)
  }

  #---get the last year inefficiently but only once
  #---we are getting the last year for which there are simulation data available
  final.year = unique.years[ num.unique.years ]
  take = which( simyear == final.year )
  agg.loss[ final.year ] = sum( losses[take] )
  occ.loss[ final.year ] = max( losses[take] )

  return(list( occ.loss=occ.loss, agg.loss=agg.loss ))
  
}




























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
      browser()
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



  browser()

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




