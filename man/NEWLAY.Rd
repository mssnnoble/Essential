% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NEWLAY.R
\name{NEWLAY}
\alias{NEWLAY}
\title{Computes layering of a graph when CCs and their parents are established}
\usage{
NEWLAY(CG)
}
\arguments{
\item{CG}{- essential graph object where CCs, CCueds, TOTPARs are correct}
}
\value{
L - the correct essential graph object
}
\description{
Input is output in the style of an essential graph object,
where the  CCs, CCueds, TOTPAR are correct and  are the CCs, CCueds, TOTPARs
of an essential graph. This routine computes
the layering, the i2ip1 parents, the DEex parents, the VIL and prVIL to give
a full essential graph object
}
