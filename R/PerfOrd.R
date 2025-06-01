#' Computes a Perfect Order of the MCs
#'
#' @param J junction tree is output of Kr(C) (Kruskal algorithm)
#' this routine computes a perfect order of the maximal cliques
#' for a CONNECTED decomposable graph.
#'
#' @return PO
#' PO[1,j] gives index of jth MC
#' PO[2,j] gives index of MC that separates PO[1,j] from earlier MCs
#' @export

perford <- function(J)
{
  a = sum(J)
  #number of edges

  Mat = J + t(J)

  PO <- matrix(0,nrow = 2,ncol=ncol(J))

  #initialising the PO. PO[1,j] gives index of jth
  #MC in the PO; PO[2,j] gives index of the MC
  #that separates the jth MC from 1,...,j-1
  #a value 0 denotes empty set.

  WHOLE = seq(1,nrow(J))
  count <- 1

  PO[,count] = c(1,0)

  set = c(1)
  SEQSET = setdiff(WHOLE,set)
  repeat
  {
    if(count == length(WHOLE))
    {
      break
    }
    count <- count+1
    f <- function(j){sum(Mat[j,SEQSET])}
    f1 <- sapply(set,f)
    ind1 = min(which(f1 >=1))
    g1 <- Mat[set[ind1],SEQSET]
    ind0 <- min(which(g1>=1))
    PO[,count] = c(SEQSET[ind0],set[ind1])
    set <- c(set,SEQSET[ind0])
    SEQSET <- SEQSET[-ind0]
}
return(PO)
}


