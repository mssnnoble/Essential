#' Computes number of complete sets size k for decomposable graph
#'
#' input is a matrix M adjacency for decomposable graph
#' on d vertices. Output is number of complete
#' sets of size k for k = 1,....,d
#'
#' @param M adjacency matrix for a decomposable graph
#'
#' @return b a vector, where b[j] gives number of complete sets size j.
#' @export

NumberComplete <- function(M)
{


  M1 = MC(M)
  S1 = Separators(M)

  #now create the mcv

  if(length(S1)==0)
  {
    a = vector(mode="integer",length=nrow(M))
    a[nrow(M)]=1
    b = vector(mode="integer",length=nrow(M))
    d = matrix(0,nrow(M),nrow(M))
    for(j in 1:nrow(M))
    {
      for(k in j:nrow(M))
      {
        d[j,k] = choose(k,j)*a[k]
      }
    }
    for(j in 1:nrow(M))
    {
      b[j]=sum(d[j,])
    }
    return(b)
  }

  f <- function(j){length(M1[[j]])}
  f1 = sapply(1:length(M1),f)
  #lengths of the MCs
  g <- function(j){length(S1[[j]])}
  g1 = sapply(1:length(S1),g)
  #lengths of the separators
  a = vector(mode="integer",length = max(f1))
  #a will become the mcv
  b = vector(mode="integer",length = ncol(M))
  #b[j] will give number of complete sets size j
  d = matrix(0,max(f1),max(f1))
  fn <- function(j){length(which(f1==j)) - length(which(g1==j))}
  a = sapply(1:(max(f1)),fn)
  for(j in 1:(max(f1)))
  {
  for(k in j:max(f1))
  {
    d[j,k] = choose(k,j)*a[k]
  }
  }
  for(j in 1:max(f1))
  {
    b[j]=sum(d[j,])
  }
  return(b)
  }
