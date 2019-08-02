#'Mean Squared Error using censored data generated from AR model

#' @description  (Accuracy of prediction)
#' One can find the mean squared error (MSE) to check how these
#' different methods (ind, dcbs, mrc, agg, mv) perform if the objective is to make prediction.
#'@param chpts changepoints that are obtained using the discussed method (ind, dcbs, mrc, agg, mv)
#' @param data.train divide generated censored data from AR model into data.train and
#' data.test. Here, we consider n=500 (size of each series) and N=100 (number of series)
#' so we have a matrix of N*n. In data.train we leave out the five data points at the end
#' of each series.
#'@param data.test Remaining dataset (five time points at the end of each series) will
#'be considered as data.test.

#'@return  return mean squared error (MSE)
#'@export
#'@seealso AR1.data, indAR,  Bin_segAR, PELT.MVar
#' @examples
#' # example
#' #mean squared error to check the accuracy of ind method using
#' #censored data generated from AR model.
#' # data generated through AR model considering 60% censoring rate
#' #(Left censoring) and missing rate is equal to zero
#'library(cpcens)
#'sim = AR1.data ( n=500 , N = 100 , K = 5 , eps = 1 , rho=0.6,
#'mu = 0,  siga = 1, rates = c(0.6,NA), Mrate=0 )
#'data=sim$data
#'n=500
#'N=100
#'# training and test
#'data.train = sim$data[,1:(n-5)]
#'data.test = sim$data[,(n-4):n]
#' ##If pen is equal to zero, penalty term will be equal to 2*log(n)
#'indar.chpts=indAR(data.train, pen=0)
#'indar.mse = predar.mse( indar.chpts , data.train , data.test )
#'indar.mse
#'#example
#'#mean squared error to check the accuracy of dcbs method using
#'#censored data generated from AR model.
#'library(cpcens)
#'# data generated through AR model considering 20% censoring rate
#' #(Right censoring) and missing rate is equal to zero
#'sim = AR1.data ( n=500 , N = 100 , K = 5 , eps = 1 , rho=0.4,
#'mu = 0,  siga = 1, rates = c(NA,0.2), Mrate=0 )
#'data=sim$data
#'n=500
#'N=100
#'# training and test
#'data.train = sim$data[,1:(n-5)]
#'data.test = sim$data[,(n-4):n]

#'dcbsar.chpts= Bin_segAR(data.train, 10)
#'dcbsar.mse = predar.mse( dcbsar.chpts , data.train , data.test )
#'dcbsar.mse




predar.mse = function( chpts , data.train , data.test ){

  N = dim(data.train)[1]
  train.n = dim(data.train)[2]

  # if MRC does not retur full vector
  if (length(chpts)<N){
    chpts = rep( unique(chpts)[1] , 100 )
  }

  ms.errs = numeric(N)
  for (i in 1:N){
    # mean of last segment
    mu = mean( data.train[ i , chpts[i]:train.n ]  )
    ms.errs[i] = mean( (data.test[i,] - mu)^2 )
  }

  return( mean(ms.errs) )
}




