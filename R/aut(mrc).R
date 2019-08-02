#' Segmenting data generated from AR1.data/MA1.data using PELT funtion.

#' @description To find changepoints using mrc method, segmenting
#'the data (obtained from AR1.data/MA1.data) using PELT (Killick, Fearnhead
#' and Eckley 2012) function in such a way
#' that cost is minimum .
#'@param data a censored data matrix obtained from AR1.data/
#' MA1.data .
#'@param beta  default 1.5*log(n).

#'@return  data

#'@export
#'@seealso AR1.data, MA1.data
#'@references Killick, R., Fearnhead, P., and Eckley, I. A. (2012). Optimal detection of changepoints
#'with a linear computational cost. Journal of the American Statistical Association,
#'107(500):1590–1598.
#' @examples

#' #example(right censoring)
#'library(cpcens)
#'n=500
#'N=100
#'# Generate censored data using AR model
#'# The size of series(n) should be greater than 200.
#' sim=AR1.data(n = 500, N = 100, K = 5, eps = 1,
#'rho = 0.4, mu = 0, siga = 1, rates = c(NA, 0.4), Mrate = 0)
#'data=sim$data
#'mrc = mrc.mean( data , beta = 1.5*log(n) )
#'mrc

#'#example(left censoring)
#'library(cpcens)
#'n=500
#'N=100
#'# Generate censored data using MA model
#'# The size of series(n) should be greater than 200.
#' sim=MA1.data(n = 500, N = 100, K = 5, eps = 1,
#'rho = 0.4, mu = 0, siga = 1, rates = c(0.6,NA), Mrate = 0)
#'data=sim$data
#'mrc = mrc.mean( data , beta = 1.5*log(n) )
#'mrc

mrc.mean = function( data , beta =  1.5*log(n) ){

  n = dim(data)[2]

  # first find the dimensions which arent all zero
  non_null_dims = which( apply( data , 1 , sum ) != 0 )

  # mean change in individual series
  PELT_data = matrix( nrow = length(non_null_dims) , ncol = n )
  for ( i in 1:length(non_null_dims) ){
    d = non_null_dims[i]
    PELT_data[i,] = PELT( data[d,] , beta )$F[1:(n)] + vec_C( data[d,] )
  }

  return(PELT_data)

}

# vector of costs for (t+1):n FOR MEAN
vec_C = function(data){

  n = length(data)
  cd = cumsum( c(0,data) )
  cd_2 = cumsum( c(0,data^2) )
  t = 0:(n-1)
  seg_cost = cd_2[(n+1)] - cd_2[(t+1)] - ( ( cd[(n+1)] - cd[(t+1)] )^2 )/(n-t)
  return(seg_cost)

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
