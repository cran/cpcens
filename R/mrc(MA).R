#'Most recent changepoints
#'
#'@description Detecting most recent changepoints using censored data generated from MA model.
#'@param n	length of series, default 500. The size of series(n) should
#' be greater than 200.
#'@param N number of series, default 100.
#'@param K number of most recent changepoints, default 5.
#'@param eps  size of the mean change at the most recent changepoint.
#'@param rho ma coefficients
#'@param mu  mean
#'@param siga  standard deviation of innovations
#'@param rates 	either a vector of length 2 or a matrix with n rows and 2 columns.
#'In the vector case, the first element indicates the left-censor rate and
#'the second element indicates the right-censor rate. Set to NA if there
#'is no censoring. Interval censored data corresponds to setting both a
#'left-censor rate and right-censor rate. The default setting indicates
#'a right-censor rate 0.2 with no left censoring. The vector case handles
#'single censoring and the matrix case is for multiple censor points.
#'In this case each column indicates the corresponding censoring for
#'each observation.
#'@param Mrate   fraction of missing values. Default is 0
#'@return an object of class 'censored' which is a list with four elements.
#' First element, 'data', is the censored time series. Second element,
#' 'mrc', indicates the most recent changepoints. Third element, 'series.mrc'
#' indicates which series is affecting from the most recent changepoint.
#'Fourth element, 'series.chpts' indicates the changepoints in each series.
#'@export
#'@examples
#' #Default example
#' library(cpcens)
#' ans<-MA1.data()
#' #example (right censoring)
#' # The size of series(n) should be greater than 200.
#' ans<-MA1.data ( n=500 , N = 100 , K = 5 , eps = 1 , rho=0.2,
#' mu = 0,  siga = 1, rates = c(NA,0.4), Mrate=0 )
## function to sim a single series with given chpts and eps
## new chpt sim function
# # length of time series
# n = 500
# # dimension
# N = 100
# # number of MRC's
# K = 5
# # mu + eps - mean of last seg
# eps = 10


MA1.data = function( n=500 , N = 100 , K = 5 , eps = 1 ,
                              rho=0.6, mu = 0,  siga = 1, rates = c(NA,0.2), Mrate=0 ) {

  ### alternative K<=10##
  true.mrc.chpts = n-sample(20*(1:10) , K , replace = FALSE)

  # which series carry MRC's
  f = floor( N/K )
  # reorder series
  tsr = sample(1:N,N)

  # locations of ordinary chpts
  chpt.pot.locs = rbinom( min(true.mrc.chpts) , 1 , prob = 0.02)
  chpt.locs = which( chpt.pot.locs == 1 )
  # prop of series each chpt affects
  alpha = runif(length(chpt.locs))

  chpts.each.series = vector("list",N)
  series.which.mrc = numeric(N)
  data = matrix(nrow=N,ncol=n)
  for (i in 1:N){

    # which of the chpts are in this series
    probs = runif(length(chpt.locs))
    wc = which( probs < alpha )

    # which most recent chpt is series affected by
    w = which(tsr == i)
    m = ceiling(w/f)
    if (m >K){
      m <- K
    }
    # which MRC affects ith series
    series.which.mrc[i] = m
    # changepoints in each series
    chpts.each.series[[i]] = c( chpt.locs[wc] , true.mrc.chpts[m] )

    sim_MA1_series_chpts = function( n , chpts , eps , rho ){


      ##########Randomly Generate Censored AR

      rcma<-function (n , ma, mu , siga, rates, Mrate )
      {
        Rates <- rates
        if (is.vector(rates))
          Rates <- matrix(rep(rates, n), byrow = TRUE, ncol = 2)
        y <- z <- mu + siga * as.vector(arima.sim(model = list(ma = ma
        ), n = n))
        iy <- yL <- yR <- rep(NA, n)
        cL <- quantile(z, Rates[, 1])
        indL0 <- z > cL
        indL <- !ifelse(is.na(indL0), TRUE, indL0)
        y <- ifelse(indL, cL, z)
        cR <- quantile(z, 1 - Rates[, 2])
        indR0 <- z < cR
        indR <- !ifelse(is.na(indR0), TRUE, indR0)
        y <- ifelse(indR, cR, y)
        indMissing <- is.element(1:n, sample(1:n, size = floor(Mrate *
                                                                 n)))
        y[indMissing] <- yL[indMissing] <- yR[indMissing] <- NA
        indL <- indL & !indMissing
        indR <- indR & !indMissing
        indo <- !(indMissing | indL | indR)
        iy <- rep("na", n)
        iy[indo] <- "o"
        iy[indL] <- "L"
        iy[indR] <- "R"
        ans <- list(y = y, iy = iy, censorPts = matrix(c(cL, cR),

                                                       ncol = 2), z = z)
        class(ans) <- "censored"
        ans
      }

      out <-rcma(n , ma=rho, mu , siga, rates, Mrate )

      data<-out$y



      mu = rnorm(1,0,2)

      data[1:chpts[1]] = data[1:chpts[1]] + mu
      if ( length(chpts) > 2 ){
        for (i in 2:( length(chpts) - 1) ){
          mu = rnorm(1,0,2)
          data[ (chpts[i] + 1):( chpts[(i+1)] ) ] =  data[ (chpts[i] + 1):( chpts[(i+1)] ) ] + mu
        }
      }
      data[ ( tail(chpts,1) + 1 ):n] = data[ ( tail(chpts,1) + 1 ):n] + mu + eps
      return(data)
    }
    data[i,] = sim_MA1_series_chpts( n , chpts.each.series[[i]] , eps , rho )

  }

  newlist = list("data" = data , "mrc" = true.mrc.chpts , "series.mrc" = series.which.mrc , "series.chpts" =  chpts.each.series )
  return(newlist)

}






