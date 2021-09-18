



















#---function written by Shree Khare which generalizes wang's analytical theory to K clusters
#---cluster number is found in the incoming array clutser.num
#---array poisson.or.nb == 1 if we have a cluster, and 0 is subset of events is treated as poisson
#---we return both pure poisson and clustered aep and oeps (as a function of the loss discretization)
#---this function generalizes the poisson.mixture function
#---we will later create a vectorized version to avoid looping over all the events
#---note, i have tested replacement of the loop using vectorized code and it was slower
#---hypothesis for why vectorized code slower is the increase in the amount of matrix multiplies
#---note that for the case of the favourite german cluster, we get exact agreement between this generalized
#---code and the original poisson.mixture routine
#---should update for zero loss thresholds
poisson.mixtureK <- function( mlr, std, rate, expsr, loss, nclusters, cluster.num, betas, poisson.or.nb, ngrid ) {
  
  #---some parameters defining computational grid (for computation of subset CEPs)
  nn = ngrid
  max.ll = 2.0*max(loss)
  dx = max.ll / nn
  xx = c((1:nn)*dx)

  #---get the overall rate and the beta distribution parameters, rates for the cluster subsets
  alpha.beta = get.alpha.beta(mlr,std)
  lambda=sum(rate)
  LC = array(0,c(nclusters))
  for (k in 1:nclusters){ take = which(cluster.num == k); LC[k] = sum(rate[take]) }
  
  #---init entire event set weighted cummulative loss distribution - 1 - FX is the CEP
  FX <- aep.fx <- array(0,nn)
  #---init weighted cummulative loss distributions per cluster
  FC <- aep.fc <- array(0,c(nclusters,nn))
  
  #---loop event by event, first get overall number of events
  fn.nevents = length(mlr)
  for ( ievent in 1:fn.nevents ) {

    #---get the ceps by integration
    this.lambda = rate[ievent]
    this.aa = alpha.beta[ievent,1]
    this.bb = alpha.beta[ievent,2]
    this.mlr = mlr[ievent]
    this.std = std[ievent]
    this.xx = xx/expsr[ievent]
    this.dx = dx/expsr[ievent]
    this.FX = pbeta(this.xx[1:nn] , shape1 = this.aa, shape2 = this.bb )
    #---FX for poisson model
    FX = FX + (this.lambda/lambda) * this.FX
    #---FC for clustered model
    FC[cluster.num[ievent],] = FC[cluster.num[ievent],] + (this.lambda/LC[cluster.num[ievent]]) * this.FX 
    
    #---do the differencing to get approximation weighted event loss pdf (note that we want to discretize for the fft)
    FX2 = pbeta(this.xx[2:nn-1] + this.dx/2, shape1= this.aa, shape2 = this.bb)
    FX1 = pbeta(this.xx[2:nn-1] - this.dx/2, shape1= this.aa, shape2 = this.bb)
    this.fx = array(0.,nn)
    this.fx[1] = pbeta(this.xx[1] - this.dx/2, shape1= this.aa, shape2 = this.bb)
    this.fx[2:nn] = FX2 - FX1
    #---aep.fx for poisson model
    aep.fx = aep.fx + (this.lambda/lambda) * this.fx
    #---aep.fc for clustered model
    aep.fc[cluster.num[ievent],] = aep.fc[cluster.num[ievent],] + (this.lambda/LC[cluster.num[ievent]]) * this.fx 
  }


  #---perform poisson calculations (have added in the OEP2 and OEP3 for the poisson assumption)
  OEP = 1 - exp(-lambda*(1-FX))
  OEP2 = 1 - exp(-lambda*(1-FX))*(1 + lambda*(1-FX))
  OEP3 = 1 - exp(-lambda*(1-FX))*(1 + lambda*(1-FX)) - (1/2)*exp(-lambda*(1-FX))*(lambda*(1-FX))^2
  EEF = (1-FX)*lambda
  fx.hat.poisson = exp( -lambda * ( 1-fft(aep.fx) ) )
  AEP = 1 - cumsum(Re(fft(fx.hat.poisson,inverse=TRUE)/nn))
  
  
  #---loop over all different beta combinations
  num.beta.combinations = length(betas[1,])
  OEPcluster <- AEPcluster <- array(0,c(num.beta.combinations,nn))
  for (ijk in 1:num.beta.combinations){
    
     #---perform calculations for the clusters - oep
     oep.factors = array(0,c(nclusters,nn))
     for (k in 1:nclusters){
       #---NB
       if (poisson.or.nb[k] == 'negative binomial') oep.factors[k,] = ( 1 - LC[k]*betas[k,ijk]*(FC[k,]-1) ) ^ (-(1/betas[k,ijk]))
       #---Poisson
       if (poisson.or.nb[k] == 'poisson') oep.factors[k,] = exp(-LC[k]*(1 - FC[k,]))
     }
     oep.product.pgfs = array(1.,nn); for (k in 1:nclusters){ oep.product.pgfs = oep.product.pgfs * oep.factors[k,] }
     OEPcluster[ijk,] = 1 - oep.product.pgfs

  
     #---perform calculations for the clusters - aep
     #---here we are using the probability generating function of the negative binomial distribution (derived from mixing gamma+poisson)
     aep.factors = array(0,c(nclusters,nn))
     for (k in 1:nclusters){
       #---NB
       if (poisson.or.nb[k] == 'negative binomial') aep.factors[k,] = ( 1 - LC[k]*betas[k,ijk]*(fft(aep.fc[k,]) - 1) ) ^ (-(1/betas[k,ijk]))
       #---Poisson
       if (poisson.or.nb[k] == 'poisson') aep.factors[k,] = exp(-LC[k] * (1 - fft(aep.fc[k,])) ) 
     }
     aep.product.pgfs = array(1.,nn); for (k in 1:nclusters){ aep.product.pgfs = aep.product.pgfs * aep.factors[k,] }
     AEPcluster[ijk,] = 1 - cumsum(Re(fft(aep.product.pgfs,inverse=TRUE)/nn))
}



  #---compute the implied overdispersion parameters for each cluster (gotten from the analytic formula)
  overdispersion = betas; for (i in 1:nclusters){ overdispersion[i,] = LC[i] * overdispersion[i,] }; overdispersion = overdispersion + 1
  
  
  #---return the clustered, and the unclustered EPs and also the EEFs
  return(list(xx=xx,EEF=EEF,OEPpoi=OEP,AEPpoi=AEP,OEP=OEPcluster,AEP=AEPcluster,OEP2poi=OEP2,OEP3poi=OEP3,overdispersion=overdispersion,rates.in.clusters=LC))
  
}






























#original version of poisson mixture code
poisson.mixture <- function ( mlr, std, rate, expsr, loss, thrs, alpha ) {
    # alpha is 1/alpha in 7.1 wang
    li = sum(rate[loss<thrs])
    lc = sum(rate[!loss<thrs])
    lambda = sum(rate)

    #assuming some grid to get the weighted event loss distributions, could expand this grid
    #note that we are essentially differencing the cdf
    alpha.beta = get.alpha.beta(mlr,std)
    ## computation of grid of losses
    nn = 2^12
    max.ll = 2.0*max(loss)
    dx = max.ll / nn
    xx = c((1:nn)*dx)
    fn.nevents = length(mlr)
    FI <- FC <- FX <- aep.fx <- aep.fi <- aep.fc <- array(0,nn)
    for ( ievent in 1:fn.nevents ) {

        this.lambda = rate[ievent]
        this.aa = alpha.beta[ievent,1]
        this.bb = alpha.beta[ievent,2]
        this.mlr = mlr[ievent]
        this.std = std[ievent]
        this.xx = xx/expsr[ievent]
        this.dx = dx/expsr[ievent]
        this.FX = pbeta(this.xx[1:nn] , shape1 = this.aa, shape2 = this.bb )
        FX = FX + (this.lambda/lambda) * this.FX
        if ( loss[ievent] < thrs ) {
          FI = FI + (this.lambda/li) * this.FX
        } else {
          FC = FC + (this.lambda/lc) * this.FX
        }

        FX2 <- pbeta(this.xx[2:nn-1] + this.dx/2, shape1= this.aa, shape2 = this.bb)
        FX1 <- pbeta(this.xx[2:nn-1] - this.dx/2, shape1= this.aa, shape2 = this.bb)

        this.fx = array(0.,nn)
        #it works perfectly fine, but why do we integrate from 0 to x(1) - dx/2 and not x(1) + dx/2 ???
        this.fx[1] = pbeta(this.xx[1] - this.dx/2, shape1= this.aa, shape2 = this.bb)

        this.fx[2:nn] <- FX2 - FX1
        aep.fx <- aep.fx + (this.lambda/lambda) * this.fx # pdf weighted average
        if ( loss[ievent] < thrs ) {
           aep.fi = aep.fi + ( this.lambda/li ) * this.fx
        } else {
           aep.fc = aep.fc + ( this.lambda/lc ) * this.fx
        }

    }

    OEP = 1 - exp(-lambda*(1-FX))
    fl_c = ( 1 - lc*alpha * ( FC - 1)) ^ (-(1/alpha))
    fl_i = exp(-li*(1-FI))
    OEPmix = 1- fl_c*fl_i

    fx.hat.poisson = exp( -lambda * ( 1-fft(aep.fx) ) )
    AEP = 1 - cumsum(Re(fft(fx.hat.poisson,inverse=TRUE)/nn))

    fx.nc = ( 1 - lc*alpha * ( fft(aep.fc) - 1)) ^ (-(1/alpha))
    fx.ni = exp( -li * ( 1-fft(aep.fi) ) )
    fx.nb = array(0.,nn)
    for ( i in 1:nn) { fx.nb[i] = fx.ni[i]*fx.nc[i] }
    AEPmix = 1 - cumsum(Re(fft(fx.nb,inverse=TRUE))/nn) #inverse fourier transform of the product of pgs of characteristics
    # key thing for a correlated cluster, is that you can't multiply the pgs for each event together, you need to use the joint expression
    AEPpoi =  1 - cumsum(Re(fft(fx.ni,inverse=TRUE))/nn)
    AEPneg =  1 - cumsum(Re(fft(fx.nc,inverse=TRUE))/nn)
    EEF = (1-FX)*lambda

    return(list(xx=xx,EEF=EEF,OEP=OEP,OEPmix=OEPmix,AEP=AEP,AEPmix=AEPmix,
                fx.ni=fx.ni,fx.nb=fx.nb,AEPpoi=AEPpoi,AEPneg=AEPneg))
}

#routine to get the alpha and beta coefficients of the BETA distribution
get.alpha.beta <- function ( mu, sig ) {
   N = length(mu)
   beta.coeff = array(0,dim=c(N,2))
   for ( i in 1:N ) {
       beta.coeff[i,1] = (mu[i]^2 * ( 1-mu[i] ))/sig[i]^2 - mu[i]
       beta.coeff[i,2] = ( beta.coeff[i,1] *( 1- mu[i]))/mu[i]
   }
   return(beta.coeff)
}

#vectorized version of get.alpha.beta routine
get.alpha.beta.vec <- function ( mu, sig ) {
  beta.coeff.mat = array(0,c(length(mu),2))
  beta.coeff.mat[,1] = (mu^2*(1-mu))/sig^2 - mu
  beta.coeff.mat[,2] = beta.coeff.mat[,1]*(1-mu)/mu
  return(beta.coeff.mat)
}



