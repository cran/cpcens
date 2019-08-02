#' Most recent changepoints from AGG method.

#' @description Detecting most recent changepoints from AGG method (detect
#'  changepoint in univariate time series). We use PELT
#'  for segmenting a time series into changing means, assuming normally
#'  distributed observations with changing mean but constant variance.
#' @param pen (penalty term)  default 200*log(dim(data)[2]. Here dim(data)[2] means
#' length of series (n).
#' @param data a censored data matrix. And then we add this data matrix column wise and resulting
#'  data use as first argument in PELT function.
#'@return  indicates the most recent changepoint in each series .
#'@export

#' @examples
#' #example
#' library(cpcens)
#' data("censoredex")
#' data=censoredex
#' n=144
#' N=100
#' agg = apply( data , 2 , sum )
#'pagg = PELT( agg , 200*log(dim(data)[2]) )
#'agg.chpts = rep( rev( pagg$cpts )[1] , N )

PELT <- function (data, pen=200*log(dim(data)[2])) {
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


