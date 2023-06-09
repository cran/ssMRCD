% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions_application.R
\name{restructure_as_list}
\alias{restructure_as_list}
\title{Restructure Data Matrix as List}
\usage{
restructure_as_list(data, neighborhood_vec)
}
\arguments{
\item{data}{data matrix with all observations.}

\item{neighborhood_vec}{numeric neighborhood assignment vector.
Should contain numbers from \code{1} to \code{N} and not leave integers out.}
}
\value{
Returns a list containing the observations per neighborhood assignment.
}
\description{
This function restructures neighborhood information given by a data matrix
containing all information and one neighborhood assignment vector. It returns a list
of data matrices used in \code{\link[ssMRCD]{ssMRCD}}.
}
\examples{

# data matrix
data = matrix(rnorm(n = 3000), ncol = 3)
N_assign = sample(x = 1:10, size = 1000, replace = TRUE)

restructure_as_list(data, N_assign)

}
