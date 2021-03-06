% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/timmaCategory.R
\name{timmaCategory}
\alias{timmaCategory}
\title{Predicting drug sensitivity with multi-class drug-target interaction data using one.sided TIMMA model}
\usage{
timmaCategory(drug_target_profile, sens, loo = TRUE, class)
}
\arguments{
\item{drug_target_profile}{the drug-target interaction data. See \code{\link{timma}}.}

\item{sens}{a drug sensitivity vector.}

\item{loo}{a logical value indicating whether to use the leave-one-out cross-validation in the model
selection process. By default, loo = TRUE.}

\item{class}{the number of classes in the drug-target interaction data}
}
\value{
A list containing the following components:
\item{dummy}{the predicted efficacy for target combinations that can be found from the training data}
\item{error}{the prediction errors}
\item{prediction}{predicted drug sensitivity}
}
\description{
A function to predict the drug sensitivity with multi-class drug-target interaction data using the
one.sided TIMMA model
}
\examples{
data(tyner_interaction_multiclass)
data(tyner_sensitivity)
results<-timmaCategory(tyner_interaction_multiclass[, 1:6], tyner_sensitivity[,1], class = 6)
}
\author{
Liye He \email{liye.he@helsinki.fi}
}

