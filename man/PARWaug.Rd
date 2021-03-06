% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abrw.R
\name{PARWaug}
\alias{PARWaug}
\title{Graph Laplacian Augmented screening method}
\usage{
PARWaug(amatrix, utility, alpha, max_it = 100, threshold = 0)
}
\arguments{
\item{amatrix}{the adjacency matrix with diagonal elements equal 0.}

\item{utility}{the original screening measure.}

\item{alpha}{the paramter balances the original ranking measure and the network information}

\item{max_it}{the maximum of iteration}

\item{threshold}{make all marignal screening measure that are smaller than the threshold equals 0 to denoise.}
}
\value{
A vector of size p contains the Graph Laplacian Augmented screening measure
}
\description{
The algorithm based on the effective power iteration method for the partial absorbing random walk algorithm
}
