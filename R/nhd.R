#' Computes the neighbour sets of the vertices
#'
#' Input is a lower triangular adjacency matrix for a decomposable graph.
#' Output gives neighbours.
#'
#' @param M - lower triangular adjacency matrix for a decomposable graph
#'
#'
#' @return O - list
#' \itemize{
#' \item O$L where O$L[[x]] gives neighbours of x, not including x
#' \item O$C where O$C[[x]] gives neighbours of x, including x
#' }
#' @export


nhd <- function(M){
  #M is a lower triangular matrix giving edges of a CC
  #output is a list of neighbours for each vertex
  O = list()
  N = M + t(M)
  L<- c()
  C<-c()
  f <-function(x)
    {
    if(sum(N[,x])==0)
      {
      vector(mode="integer",length=0L)
      }
    else
    {
        which(N[,x]==1)
      }
  }
  #for vertex x, f lists its neighbours
  L = lapply(1:nrow(M), f)

  g = function(x)
    {
    if(sum(N[,x])==0)
    {
      vector(mode="integer",length=0L)
      }
    else
    {
      sort(c(f(x),x))
        }
  }
  #for vertex x, g(x) lists the vertex
  #together with its neighbours
  C = lapply(1:nrow(M),g)
  #C[[i]] gives vertex i together with its neighbours
  #(which is a maximal clique if and only if vertex is simplicial)

  #L[[i]] gives only neighbours of i


  O$L = L
  O$C = C
  return(O)
}
