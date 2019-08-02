#' Most recent changepoints from MV method using censored AR timeseries.

#' @description Detecting most recent changepoints from MV methd (Lavielle and Teyssiere, 2006) deal with
#' multivariate data which is modeling the data within each segment as a
#'  multivariate (MV) Gaussian having a given covariance
#'   after generating censored data from AR model.
#' @param beta   default 101*log(dim(data)[2])). Here dim(data)[2] means consider
#' size(length) of series (n).
#' @param data a censored data matrix obtained from AR1.data .
#'@return  indicates the most recent changepoint in each series .
#'@export
#'@seealso AR1.data
#'@references Lavielle, M. and Teyssiere, G. (2006). Detection of multiple changepoints in multivariate
#'time series.Lithuanian Mathematical Journal, 46(3):287-306

#' @examples
#' # example (Right censoring)
#'library(cpcens)
#'# The size of series(n) should be greater than 200.
#'sim=AR1.data(n = 500, N = 100, K = 5, eps = 1,
#' rho = 0.6, mu = 0, siga = 1, rates = c(NA, 0.2), Mrate = 0)
#'data=sim$data
#'N=100
#'pmv = PELT.MVar( data , 101*log(dim(data)[2]) )
#'mv.chpts =  rep( rev( pmv$cpts )[1] , N )

# PELT-MVar

PELT.MVar = function( data , beta =  101*log(dim(data)[2]) ){

  N = dim(data)[1]
  n = dim(data)[2]

  # F[t] = optimal value of segmentation upto time t
  F = numeric(n)
  F[1] = -beta
  pen = beta
  # chpts[[t]] = a vector of changepoints upto time t (optimal)
  chpts = vector("list",n)
  chpts[[1]] = c()

  R = 1
  # useful for calculating seg costs
  cd = t( apply(data,1,cumsum) )
  cd_2 = t( apply(data^2,1,cumsum) )

  for (t in 2:n){

    cpt_cands = R
    seg_costs = numeric(length(cpt_cands))
    for (i in 1:N){
      seg_costs = seg_costs + ( cd_2[i,t] - cd_2[i,cpt_cands] ) - (cd[i,t] - cd[i,cpt_cands])^2/(t-cpt_cands)
    }

    f = F[cpt_cands] + seg_costs + pen
    F[t] = min(f)
    tau = cpt_cands[ which.min(f) ]

    chpts[[t]] = c( chpts[[ tau ]] , tau )

    # pruning step
    ineq_prune = F[cpt_cands] + seg_costs <= F[t]
    R = c( cpt_cands[ineq_prune] , t )

  }

  newList <- list("cpts" = chpts[[n]] , "F" = F )
  return(newList)

}

