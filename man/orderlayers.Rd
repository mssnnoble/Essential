% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orderlayers.R
\name{orderlayers}
\alias{orderlayers}
\title{orders layers at random ensuring a suitable ground layer is chosen}
\usage{
orderlayers(z, q)
}
\arguments{
\item{z}{- this is the list of CCs, output of CCsets}

\item{q}{- adjusted, unordered layering of the CCs (adjusted to ensure
there is a suitable layer for ground layer)}
}
\value{
b - the ordered layers, b[[1]] the ground layer. b[[j]] lists labels
of CCs in layer j.
}
\description{
The function orderlayers orders the layers at random, ensuring that a suitable
ground layer is chosen (i.e one that can admit immoralities for layer 2 CCs).
This means taking the specially created layer if this was necessary,
otherwise choosing a ground layer at random from among
suitable layers (i.e. not a layer which is a single CC which is either a
singleton or double). Other layers selected at random from the unselected layers
}
