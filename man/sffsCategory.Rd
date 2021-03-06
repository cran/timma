% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/sffsCategory.R
\name{sffsCategory}
\alias{sffsCategory}
\title{Model selection with sffs for the multi-class drug-target interaction data using one.sided TIMMA model}
\usage{
sffsCategory(profile_data, sens, sp = 1, max_k = 2, loo = TRUE, class,
  verbosity = FALSE)
}
\arguments{
\item{profile_data}{drug-target interaction data which is a matrix with drugs as row indexes and targets
as column indexes.}

\item{sens}{a drug sensitivity vector.}

\item{sp}{an integer to specify the starting point for sequential forward floating search (sffs) search
algorithm to navigate the target set space. By default, sp = 1.}

\item{max_k}{an integer to sepcify the maximum number of targets that can be selected by the sffs
algorithm. By default, max_k = 2. In practice it should not be over than 10 as the number of target combinations will increase exponentially.}

\item{loo}{a logical value indicating whether to use the leave-one-out cross-validation in the model
selection process. By default, loo = TRUE.}

\item{class}{an integer to specify the number of classes in the drug-target interaction data. For a binary drug-target
interaction data, class = 2. For a multi-class drug-target interaction data, class should be the number of classes.}

\item{verbosity}{a boolean value to decide if the information should be displayed. If it is TRUE, the information
will be displayed while the model is running. Otherwise, the information will not be displayed. By default, it is
FALSE.}
}
\value{
A list containing the following components:
\item{timma}{a list contains: the predicted efficacy matrix, prediction error and predicted drug sensitivity}
\item{k_sel}{the indexes for selected targets}
}
\description{
A function to select the most predictive targets with sffs for the multi-class drug-target interaction data using
the one.sided TIMMA model
}
\examples{
\dontrun{
data(tyner_interaction_multiclass)
data(tyner_sensitivity)
results<-sffsCategory(tyner_interaction_multiclass, tyner_sensitivity[, 1], max_k = 2, class = 6)
}
}
\author{
Jing Tang \email{jing.tang@helsinki.fi}
}

