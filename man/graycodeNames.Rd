% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/graycodeNames.R
\name{graycodeNames}
\alias{graycodeNames}
\title{Names for the predicted sensitivity matrix}
\usage{
graycodeNames(m, names, gc_row, gc_col)
}
\arguments{
\item{m}{an integer to specify the number of targets}

\item{names}{a vector of the names of the targets}

\item{gc_row}{the gray code as row indexes. It can be returned by \code{\link{graycode3}}.}

\item{gc_col}{the gray code as column indexes. It can be returned by \code{\link{graycode3}}.}
}
\value{
a list of the following components:
\item{nr}{the gray code format target names as row names.}
\item{nc}{the gray code format target names as row names.}
}
\description{
A function to make the target names in the format of gray code for the predected sensitivity matrix
}
\examples{
\dontrun{
data(tyner_interaction_binary)
data(tyner_sensitivity)
k_select<-sffsBinary(tyner_interaction_binary, tyner_sensitivity[, 1])$k_sel
gc_timma<-graycode3(length(k_select))
select_kinase_names<-dimnames(tyner_interaction_binary)[[2]][k_select]
gc_names<-graycodeNames(length(k_select), select_kinase_names, gc_timma$gc_row, gc_timma$gc_col)
}
}
\author{
Liye He \email{liye.he@helsinki.fi}
}

