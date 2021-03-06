% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/normalizeSensitivity.R
\name{normalizeSensitivity}
\alias{normalizeSensitivity}
\title{Normalize the drug sensitivity data}
\usage{
normalizeSensitivity(IC50, method = "minMax")
}
\arguments{
\item{IC50}{a vector contains the drug sensitivity in the form of IC50.}

\item{method}{a string to specify the method used to normalize the sensitivity data. If it is
"minMax", the sensitivity is scaled by (Max_IC50-IC50)/(Max_IC50-Min_IC50). If it is "logistic",
it is scaled by 1/(1+exp(-1/IC50)). If it is "hyperbolic", it is scaled by tanh(1/IC50).}
}
\value{
A vector contains the normalized drug sensitivity data.
}
\description{
A function to normalize the drug sensitivity data to [0,1]
}
\examples{
data(tyner_sensitivity)
normalizedSensitivity<-normalizeSensitivity(tyner_sensitivity[,1])
}
\author{
Jing Tang \email{jing.tang@helsinki.fi}
}

