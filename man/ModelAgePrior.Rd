% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model_AgeOnly.R
\docType{data}
\name{ModelAgePrior}
\alias{ModelAgePrior}
\title{JAGS Models for OSL Age Estimation in \code{\link{Compute_AgeS_D}}}
\format{
BayLum's Old Age Model (wrong vector law)

Jeffrey's Age Model with log-uniform order settings

Nicholls' Age Model applied on ages directly

Jeffrey's Age Model with conditional setting

Independance Age Model
}
\usage{
ModelAgePrior
}
\description{
JAGS models used to estimate true OSL ages based on data obtained from the Bayesian OSL analysis performed
by the function  \code{\link{Palaeodose_Computation}}.
}
\details{
These models take as input the estimated dose response ($D$) from  \code{\link{Palaeodose_Computation}}
along with the structured data matrix computed by  \code{\link{create_MeasuresDataFrame}}.
The models are designed to refine age estimation by integrating these measurements into a Bayesian framework.
}
\references{
To cite this package, please use: citation("BayLum")
}
\keyword{datasets}
