% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/removeedgemat.R
\name{removeedgemat}
\alias{removeedgemat}
\title{Edges which can be removed from a decomposable graph to give a decomposable graph}
\usage{
removeedgemat(M)
}
\arguments{
\item{M}{- lower triangular adjacency matrix for decomposable graph}
}
\value{
RUE - lower triangular matrix with 1 if there exists an edge
 which is removable and 0 otherwise.
 @export
}
\description{
Given a lower triangular matrix M giving directed edges of a decomposable
graph, this routine computes which are removable so that their removal
results in a decomposable graph. This uses the theorem from Green and Thomas
pointing out that an edge is removable if and only if it does not appear in
any separator.
}
