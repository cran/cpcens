#' Most recent changepoints from MV method.

#' @description Detecting most recent changepoints from MV methd (Lavielle and Teyssiere, 2006) deal with
#' multivariate data which is modeling the data within each segment as a
#'  multivariate (MV) Gaussian having a given covariance.
#' @param beta   default 101*log(dim(data)[2])). Here dim(data)[2] means consider
#' size(length) of series (n).
#' @param data a censored data matrix.
#'@return  indicates the most recent changepoint in each series .
#'@export

#' @examples
#' # example
#'library(cpcens)
#'data("censoredex")
#'data=censoredex
#'N=100
#'n=144
#'pmv = PELT.MV( data , 101*log(dim(data)[2]) )
#'mv.chpts =  rep( rev( pmv$cpts )[1] , N )

# PELT-MV

PELT.MV = function( data , beta =  101*log(dim(data)[2]) ){

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

