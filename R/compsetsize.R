#' computes number of complete sets size k for each k
#'
#'This simply computes for each k, the number of complete sets of size k in
#'layer 1 and puts in a vector L, where L[k] is number of complete sets
#'of size k in layer 1.
#'
#'This is used to compute probability of a chain graph; if we have N vertices
#'in layer 1, and a CC in layer 2 has m parents, we are not allowed a complete
#'set; there are (N choose m) - L[m] possible ways to take a valid subset of m
#'vertices

#'
#' @param CG - an essential graph object

#'
#'
#' @return L is a vector where L[k] gives the number of complete sets of size k
#' in layer 1.
#'
#'
#'@export

compsetsize <- function(CG)
  {
  #L a list, L[k] gives number of complete sets of size k in layer 1

  CCs = CG$CC
  CCued = CG$CCued
  ordlay = CG$ordlay
  VIL = CG$VIL

g <- function(j){NumberComplete(CCued[[ordlay[[1]][j]]])}
g1 <- lapply(1:length(ordlay[[1]]),g)
len = function(j){length(g1[[j]])}
len1 = sapply(1:length(g1),len)

dim = length(VIL[[1]])
Hmat = matrix(0,length(ordlay[[1]]),dim)

for(i in 1:length(ordlay[[1]]))
{
  for(j in 1:length(g1[[i]]))
  {
    Hmat[i,j] = g1[[i]][j]
  }
}
L = colSums(Hmat)
  return(L)
}
