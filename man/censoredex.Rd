\name{censoredex}
\alias{censoredex}
\docType{data}
\title{Successive readings of a toxic substance in the Niagara River near Fort Erie, Ontario and then converting it into panel data.}
\description{ Niagara River at Fort Erie, successive readings of 12-Dichloro in units of ng/L measured approximately biweekly.
}
\usage{data("censoredex")}
\format{
  A data matrix with 100 rows annd 144 columns.
}
\details{Dr. Abdel El-Shaarwai provided through Environment Canada some a special water quality time series that is of great practical interest. The time series is from Station ON02HA0019 (Fort Erie) on the water quality of the Niagara River. There are more than 500 water quality parameters or variables of interest in this river. The water quality in this river is montiored by a joint U.S./Canada committee. The time series data are about water quality of the Niagara River where toxicity is of great interest,i.e., the chemical known as 12-Dichloro which is measured in units of ng/L (nanogram/liter) when dissolved in water. The data are measured approximately every two weeks over the period from March 1, 2001 to March 22, 2007, and in-total we have 144 left censored values. The observed censoring rate is CR = (21*100)/144 = 14\%. The detection level for 12 Dichloro after March 24, 2005 dropped from 0.214 to 0.0878. There is only one censored value at Julian day number 1807, 75 complete observations and 20 censored ones before the change in censoring while there were 48 complete observations and only one censored observation from March 24, 2005 to the last observation on March 22, 2007. Although the data are time series, not panel data and to convert it into panel we simulated the panel data.}
\source{ Dr. Abdel El-Shaarwai, Environment Canada}
\references{N. M. Mohammad (2014). Censored time series analysis. Ph.D. Thesis, Western University
}
\examples{
library(cpcens)
data("censoredex")
data=censoredex
dim(censoredex)
}
\keyword{datasets}
