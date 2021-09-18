
# SPREAD COEFFICIENTS
my_acoeff <- function ( n,M ) {
   ac = array(0,dim=c(M,n))
   for ( i in 1:M ) {
      ac[i,1]= 1/factorial(i)
   }
   for ( j in 2:n ) {
      ac[1,j]=0
   }
   for ( i in 2:M ) {
      for ( j in 2:n ) {
         ac[i,j] = 1/i * ( ( i - ( j-1) ) * ac[(i-1),(j-1)] +
                           ( (j-1) + 1  ) * ac[(i-1),j] )
      }
   }
   return(ac)
}







#small coding changes made by Shree Khare to correct a few small bugs
wp_oep_aep <- function ( rate, loss, expsr, std, mlr, DIST, alpha, maxloss , ngrid) {


   N = length(rate)
   nthrs = ngrid 

#-----------------------------------------
# FREQUENCY DISTRIBUTION GENERATION
#-----------------------------------------
   LAMBDA = sum(rate)

   if ( DIST == 'poisson' ) {
      MAXEVENTS = pmax(20,qpois((1-1/20000),LAMBDA))
      feq.dist = vector()
      tmp = vector()
      for ( k in 1:(MAXEVENTS+1) ){
         #feq.dist[k] = exp(-LAMBDA) * LAMBDA^(k-1) / factorial(k-1)
         feq.dist[k] = dpois((k-1),LAMBDA)
      }
   } else {

     #crucially we are picking from a negative binomial
     
      mu = LAMBDA
      va = LAMBDA + LAMBDA^2 * (alpha)
      r  = mu^2/(va-mu)
      p  = r/(r+mu)
      MAXEVENTS = pmax(20,qnbinom((1-1/20000),prob=p,size=r))
      feq.dist = vector()
      for ( k in 1:(MAXEVENTS+1) ) {
         feq.dist[k] = dnbinom((k-1),prob=p,size=r)
      }
   }
   print(paste('MAXEVENTS', MAXEVENTS))
#-----------------------------------------
# LOSS THRESHOLD TABLE (LTT)
#-----------------------------------------
   dx = maxloss*2 / nthrs
   xx = c((1:nthrs)*dx)

   LTT = array(0,dim=c(nthrs,2))
   LTT[,1] = xx


#-----------------------------------------
# BETA COEFFICIENTS from ELT
#-----------------------------------------

   beta.coeff = get.alpha.beta(mlr,std)

#-----------------------------------------------------------
# CONDITIONAL EXCEEDANCE PROBABILITY (CEP)
#----------------------------------------------------------
   FX <- CEP <- array(0,nthrs)
   for ( i in 1:N ) {
     if ( i %% 100 == 0 ) { print(i) }
      this.aa = beta.coeff[i,1]
      this.bb = beta.coeff[i,2]
      this.lambda = rate[i]
      this.xx = LTT[,1]/expsr[i]
      this.FX = pbeta(this.xx , shape1 = this.aa, shape2 = this.bb )
      FX = FX + (this.lambda/LAMBDA) * this.FX
      CEP = CEP + (this.lambda/LAMBDA) * (1-this.FX)
   }

#-----------------------------------------------------------
# OCCURRENCE EXCEEDANCE PROBABILITY (OEP)
#----------------------------------------------------------

   OEP = vector('numeric',nthrs+1)
   if ( DIST == 'poisson' ) {
      OEP[1] = 1 - exp(-LAMBDA)
   
      OEP[2:(nthrs+1)] = 1 - exp( - LAMBDA * CEP )
   } else {

     OEP[1] = 1 - (1 + LAMBDA * alpha * (CEP) )^(-1/alpha)
     #why is the index here starting at 2?
     OEP[2:(nthrs+1)] =  1 - ( 1 + LAMBDA * alpha *  ( CEP )) ^ (-(1/alpha))
   }

#-----------------------------------------------------------
# SEVERITY DENSITY FUNCTION (SDF)
#----------------------------------------------------------

   nthrs = length(CEP)
   SDF = vector('numeric',nthrs)
   hold = c(1,CEP)
   for ( i in 1:nthrs ) {
     SDF[i] = hold[i] - hold[i+1]
   }


#-----------------------------------------------------------
# AGGREGATE EXCEEDANCE PROBABILITY (AEP)
#----------------------------------------------------------

   # zero filling for no-wrap convolution
   tmp = vector('numeric',nthrs)
   tmp[1:length(SDF)]=SDF
   SDF = tmp

   AC = my_acoeff(nthrs,MAXEVENTS)
   Sconv = SDF
   dAEP = array(0,dim=c(MAXEVENTS,nthrs))
   for ( k in 1:(MAXEVENTS-1) ) {
      print(paste('EVENT ; ',k))
      hold = convolve(AC[k,],rev(Sconv),type='o')[1:nthrs]
      dAEP[k,] = feq.dist[k+1]*hold[1:nthrs]
      Sconv = convolve(SDF,rev(Sconv),type='o')[1:nthrs]
   }
   #Shree Khare has made change to the code below that that we pick MAXEVENTS
   hold = convolve(AC[k+1,],rev(Sconv),type='o')[1:nthrs]
   dAEP[MAXEVENTS,] = feq.dist[MAXEVENTS+1]*hold[1:nthrs]
   cAEP = colSums(dAEP,na.rm=TRUE)
   
   AEP= vector('numeric',nthrs+1)
   AEP[1] = 1-feq.dist[1]
   for ( i in 2:(nthrs+1) ) {
       AEP[i] = AEP[i-1]-cAEP[i-1]
   }

   #note that the AEP is of ngrid + 1 dimension, and includes the zero loss level, also the OEP
   return(list(xx=LTT[,1],OEP=OEP,AEP=AEP))

}
























#RL.clusters = get.RL.clusters( CEP.clusters, dAEP.clusters, MAXEVENTS.clusters, num.beta.combinations, num.clusters, poisson.or.nb, betas, ngrid, cep.daep$xx, rates.clusters)
#function written by Shree Khare to get OEP/AEP for K clusters consistent with convolutions in RL
#passing in CEP.clusters(num.clusters,ngrid), dAEP.clusters(MAXEVENTS.clusters*num.clusters,ngrid)
#note that we have written the code to include the zero loss level
#----------
get.RL.clusters <- function( CEP.clusters, dAEP.clusters, MAXEVENTS.clusters, num.beta.combinations, num.clusters, poisson.or.nb, betas, ngrid, rates.clusters, xx) {



  OEP.clusters <- AEP.clusters <- array(0,c(num.beta.combinations,(ngrid+1)))

  
  for (i in 1:num.beta.combinations){
    
    #--->>>
    #OEP (coded to include the zero loss level)
    cdf.mix = array(1,c(ngrid+1))
    for (cn in 1:num.clusters ){
      CEP.cn = c(1,CEP.clusters[cn,])
      if (poisson.or.nb[cn] == 'poisson') {cdf.mix = cdf.mix * exp(-rates.clusters[cn] * CEP.cn) }
      if (poisson.or.nb[cn] == 'negative binomial') {cdf.mix = cdf.mix * (1 + rates.clusters[cn] * betas[cn,i] * CEP.cn)^( -1/betas[cn,i] )}
    }
    OEP.clusters[i,] = 1. - cdf.mix
    #--->>>

    #--->>>
    #AEP (coded to include the zero loss level)
    dum.int = 0
    charac.func = array(1,c(ngrid+1))
    for (cn in 1:num.clusters){
 
      max.events.cn = MAXEVENTS.clusters[cn]
      daep.cn = array(0,c(max.events.cn, ngrid)); daep.cn = dAEP.clusters[(dum.int + 1) : (dum.int + max.events.cn),];  dum.int = dum.int + max.events.cn
      
      if (poisson.or.nb[cn] == 'poisson') {
        freq.dist = array( dpois(c(1:max.events.cn), rates.clusters[cn]), c(max.events.cn, ngrid) )
        freq.dist.0 = dpois(0,rates.clusters[cn])
        daep.cn = freq.dist* daep.cn
      }

      
      if (poisson.or.nb[cn] == 'negative binomial') {
        mu = rates.clusters[cn]; va = rates.clusters[cn] + rates.clusters[cn]^2 * (betas[cn,i]); r = mu*mu/(va-mu); p = r/(r+mu)
        freq.dist = array( dnbinom( c(1:max.events.cn), prob=p, size=r), c(max.events.cn, ngrid) )
        freq.dist.0 = dnbinom(0, prob=p, size=r)
        daep.cn = freq.dist *  daep.cn
      }

      caep.cn = colSums(daep.cn, na.rm=TRUE)
      
      #get the aep for this particular cluster
      aep = vector('numeric', (ngrid + 1))
      aep[1] = 1 - freq.dist.0
      for (ijk in 2:(ngrid+1)){ aep[ijk] = aep[ijk-1] - caep.cn[ijk-1] }
      
      cdf.nc = 1 - aep
      #note that in the original code, the charac.func is only on the ngrid (not ngrid + 1)
      #should be 100% that this is not having an impact
      daep = cdf.nc - c(0,cdf.nc[1:ngrid])
      charac.func = charac.func * fft(daep)
    }
    AEP.clusters[i,] = 1 - cumsum(Re(fft(charac.func,inverse=TRUE))/(ngrid+1))  # finally compute the AEP
    #--->>>
    
  }


   return(list(xx=xx,OEP=OEP.clusters,AEP=AEP.clusters))
  

  
}
#----------

















#routine to get the cep and daep (with no multiplication by the frequencies)
#idea is to, for one particular definition of cluster groupings, strip out the calculations that need not
#be repeated, so we can try many different configurations of the variances
#note that we are always getting the MAXEVENTS using a poisson distribution with twice the input rate
get.cep.daep <- function ( rate, loss, expsr, std, mlr, maxloss , ngrid, alpha) {
  
   N = length(rate)
   nthrs = ngrid
   LAMBDA = sum(rate)
   
   #determine the MAXEVENTS, determined from poisson with TWICE the input LAMBDA, which should 'cover' NB cases
   MAXEVENTS = pmax(20,qpois((1-1/20000),2*LAMBDA))
   
   dx = maxloss*2 / nthrs
   xx = c((1:nthrs)*dx)
   LTT = array(0,dim=c(nthrs,2))
   LTT[,1] = xx
   beta.coeff = get.alpha.beta(mlr,std)

   #get the CEP array CEP(ngrid)
   FX <- CEP <- array(0,nthrs)
   for ( i in 1:N ) {
#     if ( i %% 100 == 0 ) { print(i) }
      this.aa = beta.coeff[i,1]
      this.bb = beta.coeff[i,2]
      this.lambda = rate[i]
      this.xx = LTT[,1]/expsr[i]
      this.FX = pbeta(this.xx , shape1 = this.aa, shape2 = this.bb )
      FX = FX + (this.lambda/LAMBDA) * this.FX
      CEP = CEP + (this.lambda/LAMBDA) * (1-this.FX)
   }

   nthrs = length(CEP)
   SDF = vector('numeric',nthrs)
   hold = c(1,CEP)
   #now difference the cdf (1-ep) to get the value of the density on the grid
   for ( i in 1:nthrs ) { SDF[i] = hold[i] - hold[i+1] }

   #now convolve the severity distribution n=1,n=2,n=3 etc. times to take care of all possible aggregations, upto MAXEVENTS
   AC = my_acoeff(nthrs,MAXEVENTS)
   Sconv = SDF
   dAEP = array(0,dim=c(MAXEVENTS,nthrs))
   for ( k in 1:(MAXEVENTS-1) ) {
#      print(paste('EVENT ; ',k))
      hold = convolve(AC[k,],rev(Sconv),type='o')[1:nthrs]
      dAEP[k,] = hold[1:nthrs]
      #dAEP[k,] = feq.dist[k+1]*hold[1:nthrs]
      Sconv = convolve(SDF,rev(Sconv),type='o')[1:nthrs]
   }	
   #Shree Khare has made change to the code below so that we pick MAXEVENTS
   hold = convolve(AC[k+1,],rev(Sconv),type='o')[1:nthrs]
   dAEP[MAXEVENTS,] = hold[1:nthrs]

   return(list(xx=LTT[,1],CEP=CEP,dAEP=dAEP, MAXEVENTS = MAXEVENTS))

}

   

