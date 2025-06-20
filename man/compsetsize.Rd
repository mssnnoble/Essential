% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compsetsize.R
\name{compsetsize}
\alias{compsetsize}
\title{computes number of complete sets size k for each k}
\usage{
compsetsize(CG)
}
\arguments{
\item{CG}{- an essential graph object}
}
\value{
L is a vector where L[k] gives the number of complete sets of size k
in layer 1.
}
\description{
This simply computes for each k, the number of complete sets of size k in
layer 1 and puts in a vector L, where L[k] is number of complete sets
of size k in layer 1.
}
\details{
This is used to compute probability of a chain graph; if we have N vertices
in layer 1, and a CC in layer 2 has m parents, we are not allowed a complete
set; there are (N choose m) - L[m] possible ways to take a valid subset of m
vertices
}
