% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centrescale.R
\name{centrescale}
\alias{centrescale}
\title{Centres and scales a data matrix}
\usage{
centrescale(new)
}
\arguments{
\item{new:}{n by p data set each row is an instantiation of a p-variate
observtion.}
}
\value{
out - n by p matrix where out[i,j] obtained from new[i,j] by
subtracting column average and dividing through by
column standard deviation
}
\description{
Centres and scales a data matrix
}
