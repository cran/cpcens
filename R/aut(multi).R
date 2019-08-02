#' Most recent changepoints from mrc method using censored AR/MA timeseries.

#' @description Detecting most recent changepoints from mrc method consisting
#' of many related univariate timeseries (Bardwell, Eckley, Fearnhead, and Smith, 2016)
 #'after generating censored data from AR/MA model
#'and pools information across the time-series by
#'solving the K-median problem using tb.raw (Teitz and Bart, 1968).
#' @param mrc data obtained from mrc.mean
#'@param pmax  Maximum number of most recent changepoints
#' to search for. Default value pmax=10.
#'@param alpha The variable specific penalty used to penalise
#' the addition of a given changepoint into a given variable. Default
#' value alpha = 2.
#' @param n length of series
#'@param elbow.thresh default 0.5.
#'@return  indicates the most recent changepoint in each series .
#' @importFrom Rdpack reprompt
#' @importFrom tbart tb.raw
#'@export

#'@seealso mrc
#'@references Teitz, M. B. and Bart, P. (1968). Heuristic methods for estimating the generalized vertex
#'median of a weighted graph. Operations Research, 16(5):955–961.
#'
#'Bardwell, L., Fearnhead, P., Eckley, I. A., Smith, S., and Spott, M. (2019). Most recent
#'changepoint detection in panel data. Technometrics, 61(1):88–98.


#'@examples

#'#'#example(left censoring)
#'library(cpcens)
#'n=300
#'N=100
#'# Generate censored data using MA model
#' sim=MA1.data(n = 300, N = 100, K = 5, eps = 1,
#'rho = 0.4, mu = 0, siga = 1, rates = c(0.6,NA), Mrate = 0)
#'data=sim$data
#'mrc = mrc.mean( data , beta = 1.5*log(n) )
#'c = multiple.mrc(  mrc , pmax=10, alpha = 2 , elbow.thresh = 0.5, n=500 )
#'p.hat = c$MDL
#'mrc.chpts = c$locs[[p.hat]][ c$affected[[p.hat]] ]
#'mrc.chpts


multiple.mrc = function( mrc , pmax=10  , alpha = 2 , elbow.thresh = 0.5, n=500){

  # no of series
  N = dim(mrc)[1]
  location.vec = 0:(n-1)

  cost = numeric( pmax )
  mmrc = vector( "list" , pmax )
  affected = vector( "list" , pmax )
  # p=1 separate as simpler
  mmrc[[1]] = tb.raw( mrc , c(1) )
  # say all series are affected by this 1 change
  index = rep(1,times=N)
  # find which affected (more evidence above threshold)
  #index[ data[ , mmrc[[1]] ] > - alpha ] <- 0
  affected[[1]] = index
  # objective cost
  cost[1] =  sum( mrc[ , mmrc[[1]] ] )

  if (pmax>1){

    for (p in 2:pmax){

      # mmrc[[i]] gives locations of the i best locations
      mmrc[[p]] = tb.raw( mrc , c(1:p) )
      # affected[[i]], gives each dimension a label from 1:i
      # depending on which change it is associated with
      affected[[p]] = apply( mrc[ , mmrc[[p]] ] , 1 , which.min )
      # cost[[i]] gives the objective cost for solving with i different changes/sets
      csum = 0
      for (i in 1:p){
        wa = which( affected[[p]] == i)
        csum = csum + sum( mrc[ wa , mmrc[[p]][i] ] )
      }
      cost[p] = csum

    }

  }

  locations = vector("list",pmax)
  for ( i in 1:pmax ){
    locations[[i]] = location.vec[ mmrc[[i]] ]
  }

  BIC = which.min( cost + ( (1:pmax) )*log(N)*log(N*n) )  ###
  MDL = which.min(cost + (N*log(1:pmax)+ (1:pmax)*log(n))/log(2)) ##MDL

  # selecting best p - Lavielle
  J = numeric(pmax)
  for (k in 1:pmax){
    J[k] = ( ( cost[pmax] - cost[k] )/( cost[pmax] - cost[1] ) ) * ( pmax - 1 ) + 1
  }
  D = numeric(pmax-1)
  D[1] = Inf
  for ( k in 2:(pmax-1) ){
    D[k] = J[k-1] - 2*J[k] + J[k+1]
  }
  elbow = max( which( D > elbow.thresh ) )

  newlist = list( "locs" = locations , "affected" = affected ,
                  "cost" = cost , "bic" = BIC, "elbow" = elbow, "MDL"=MDL )
  return(newlist)

}
