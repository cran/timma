% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/timmaSearchBinary.R
\name{timmaSearchBinary}
\alias{timmaSearchBinary}
\title{Prediction in the search space with one.sided TIMMA model}
\usage{
timmaSearchBinary(profile_k, space, sens, loo = TRUE)
}
\arguments{
\item{profile_k}{current selected drug-target interaction data}

\item{space}{the search space returned by \code{\link{searchSpace}} function}

\item{sens}{drug sensitivity data}

\item{loo}{a logical value indicating whether to use the leave-one-out cross-validation in the model
selection process. By default, loo = TRUE.}
}
\value{
the prediction error
}
\description{
A function to return the prediction error in the search space for sffs
}
\examples{
data(tyner_interaction_binary)
data(tyner_sensitivity)
num<-length(tyner_sensitivity[,1])
k_set<-rep(0, dim(tyner_interaction_binary)[2])
k_set[c(1,2,3)]<-1
space<-searchSpace(num, k_set, tyner_interaction_binary, tyner_sensitivity[,1])
profile_k<-tyner_interaction_binary[, which(k_set==1)]
error<-timmaSearchBinary(profile_k, space, tyner_sensitivity[,1])
}
\author{
Liye He \email{liye.he@helsinki.fi}
}

