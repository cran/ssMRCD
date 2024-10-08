% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparsePCA_helpers.R
\name{eval_objective}
\alias{eval_objective}
\title{Objective function value for local sparse PCA}
\usage{
eval_objective(PC, eta, gamma, COVS)
}
\arguments{
\item{PC}{vectorised component to evaluate.}

\item{eta}{degree of sparsity.}

\item{gamma}{distribution of sparsity between groupwise (\eqn{\gamma = 1}) and entrywise (\eqn{\gamma = 0}) sparsity.}

\item{COVS}{list of covariance matrices used for PCA}
}
\value{
Returns value of the objective function for given \code{v}.
}
\description{
Objective function value for local sparse PCA
}
\examples{
S1 = matrix(c(1, 0.9, 0.8, 0.5,
              0.9, 1.1, 0.7, 0.4,
              0.8, 0.7, 1.5, 0.2,
              0.5, 0.4, 0.2, 1), ncol = 4)
S2 = t(S1)\%*\% S1
S2 = S2/2

eval_objective(PC = c(1,0,0,0,sqrt(2),0,0,-sqrt(2)),
               eta = 1, gamma = 0.5,
               COVS = list(S1, S2))
}
