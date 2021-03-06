% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateNu.R
\name{estimateNu}
\alias{estimateNu}
\title{Estimate noise intensity (nv) for non-DMCs,
using maximum likelihood estmiation.}
\usage{
estimateNu(z, phi, maxit = 50, beginP = 20)
}
\arguments{
\item{z}{A matrix. Observated mixed turmor samples.}

\item{phi}{mode of beta-values of each row in pure nomral samples y.}

\item{maxit}{A postive integer. The iteration number used in maximum likelihood.}

\item{beginP}{A number, where the method start to search from for root.

return  estimated nv (noise intensity)}
}
\description{
Estimate noise intensity (nv) for non-DMCs,
using maximum likelihood estmiation.
}
