#'Simulate censored AR time series
#'
#'@description Randomly Generate Censored AR
#'@param n	length of series, default 100.
#'@param ar  ar coefficients
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
#'@return an object of class 'censored' which is a list with three elements.
#' First element, 'y', is the censored time series. Second element, 'iy',
#'indicates for each observed valued "o", "L", "R", NA according to
#'whether the value is fully observed, left-censored, right-censored,
#' or missing. Third element, 'censorPts', is a matrix with 2 columns
#' indicating the censor point or NA if no censoring is applicable.
#'Note that censorPts does not indicate if the observation was actually
#' censored since this depends on the unknown latent variable.
#'  An observation is censored if and only if the corresponding entry in
#' iy is either "L" or "R".  See example below
#'
#'@export
#'@examples
#' #Default example
#' library(cpcens)
#' ans<-rcar()
#' #example (right censoring)
#' ans = rcar (n=100 , ar=0.2, mu=0 , siga=1, rates=c(NA,0.7), Mrate=0 )
#' #example (left censoring)
#' ans = rcar (n=100 , ar=0.7, mu=0 , siga=1, rates=c(0.3,NA), Mrate=0 )
#' #example (interval censoring)
#' ans = rcar (n=100 , ar=0.7, mu=0 , siga=1, rates=c(0.1,0.1), Mrate=0 )

rcar<-function (n=100 , ar=0.6, mu=0 , siga=1, rates=c(NA,0.2), Mrate=0 )
{
  Rates <- rates
  if (is.vector(rates))
    Rates <- matrix(rep(rates, n), byrow = TRUE, ncol = 2)
  y <- z <- mu + siga * as.vector(arima.sim(model = list(ar = ar
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


