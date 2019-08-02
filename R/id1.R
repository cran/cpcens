#' Most recent changepoints from ind method using censored AR timeseries.
#'
#' @description Detecting most recent changepoints from ind method
#'  (analyzes all the series independently in the panel data and in
#'  each given series estimate the most recent changepoint) after
#' generating censored data from AR model. We use PELT for segmenting a time series
#'  into changing means, assumes normally distributed observations
#' with changing mean but constant variance.
#' @param pen default 2*log(n). If pen is equal to zero, penalty term will be equal
#' to 2*log(n).
#' @param data a censored data matrix obtained from AR1.data .
#'@return  indicates the most recent changepoint in each series .
#'@export
#'@seealso AR1.data
#' @examples
#' #Default example
#' library(cpcens)
#' sim=AR1.data()
#'data=sim$data
#' ans = indAR(data,pen)
#' #example(right censoring)
#' #The length of series(n) should be greater than 200.
#' sim=AR1.data(n = 500, N = 100, K = 5, eps = 1,
#'rho = 0.6, mu = 0, siga = 1, rates = c(NA, 0.2), Mrate = 0)
#'data=sim$data
#' ans = indAR(data,pen=0)

indAR = function( data , pen=0){

  n = dim(data)[2]
  N = dim(data)[1]

  icp = numeric(N)
  for (j in 1:N){
    cpts = PELT(data[j,],2*log(n))$cpts
    icp[j] = cpts[length(cpts)]
  }

  return( icp )

}


PELT <- function(data,pen){

  n = length(data)
  # if unspecified penalty make it BIC
  if (pen==0){
    pen = 2*log(n)
  }

  # F[t] = optimal value of segmentation upto time t
  F = numeric(n+1)
  F[1:2] = c(-pen,0)

  # chpts[[t]] = a vector of changepoints upto time t (optimal)
  chpts = vector("list",n+1)
  chpts[[1]] = NULL
  chpts[[2]]= c(0)

  R = c(0,1)

  # useful for calculating seg costs
  cd = cumsum(c(0,data))
  cd_2 = cumsum(c(0,data^2) )

  for (t in 2:n){

    cpt_cands = R
    seg_costs =  cd_2[t+1] - cd_2[cpt_cands+1] - ((cd[t+1] - cd[cpt_cands+1])^2/(t-cpt_cands))

    f = F[cpt_cands+1] + seg_costs + pen
    F[t+1] = min(f)
    tau = cpt_cands[ which.min(f) ]

    chpts[[t+1]] = c( chpts[[ tau+1 ]] , tau )

    # pruning step
    ineq_prune = F[cpt_cands+1] + seg_costs <= F[t+1]
    R = c( cpt_cands[ineq_prune] , t )

  }

  # if only a single chpt detected at 1 then no changepoint in series, i.e., 0.
  cpts = chpts[[n+1]]
  newList <- list("cpts" = cpts , "F" = F )
  return(newList)

}


