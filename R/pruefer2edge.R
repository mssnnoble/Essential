#'Computes the lower triangular adjacency matrix for tree defined by Pruefer code.
#'
#' @param p2e: a pruefer code of length n-2
#' @param CComp: corresponding chain component of length n
#'
#'
#' @return mat: an n by n lower triangular adjacency matrix where
#' 1's correspond to edges given by the pruefer code.
#' @export



pruefer2edge <- function(p2e,CComp){
  #p2e is a pruefer code of length n-2
  #CC is the chain component
  #pruefer code is integer(0) if CC has length 1 or 2
  #this function turns it into the spanning tree
  #and expresses the spanning tree as a lower triangular matrix.

  n = length(CComp)
  if(n==1)
  {
    mat <- matrix(0,1,1)
    return(mat)
  }
  else
  {
    if(n==2)
    {
      mat <- matrix(0,2,2)
      mat[2,1] <- 1
      return(mat)
    }
    else
    {
      if(n >= 3)
      {
        parents = integer(n - 1)
        kids = integer(n - 1)
        PC = setdiff(seq_len(n), p2e)
        P = p2e
        for (k in 1:(n - 2)) {
          j = which.min(PC)
          j1 = PC[j]
          kids[k] = j1
          parents[k] = P[k]
          PC = PC[-j]
          if (k == (n - 2)) {
            PC = c(PC, P[k])
            break
          }
          if (!any(P[k] == P[(k + 1):(n - 2)])) {
            PC = c(PC, P[k])
          }
        }
        kids[n - 1L] = PC[1L]
        parents[n - 1L] = PC[2L]
        edge = cbind(kids, parents)
        edges = t(apply(edge,1,sort))

        mat <- matrix(0,n,n)
        mat[edges]<-1
        mat <-t(mat)
        return(mat)
      }
    }
  }
}
