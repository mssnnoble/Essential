#' Edges which can be removed from a decomposable graph to give a decomposable graph
#'
#'Given a lower triangular matrix M giving directed edges of a decomposable
#'graph, this routine computes which are removable so that their removal
#'results in a decomposable graph. This uses the theorem from Green and Thomas
#'pointing out that an edge is removable if and only if it does not appear in
#'any separator.
#'
#' @param M - lower triangular adjacency matrix for decomposable graph
#' @return RUE - lower triangular matrix with 1 if there exists an edge
#'  which is removable and 0 otherwise.
#'  @export

removeedgemat <- function(M)
{
  #this function computes the undirected edges
  #that can be removed from an undirected decomposable graph
  #so that the result is an undirected decomposable graph.
  #edge can be deleted if and only if it appears in exactly
  #one max clique.

  Cl <- MC(M)
  Sep = Separators(M)
  Jtr <- Kr(Cl)

  RUE <- M
  for(x in 1:nrow(M))
  {
    for(y in 1:nrow(M))
    {
      if(M[x,y] == 1)
      {
        for(k in 1:length(Sep))
        {
          if(all(c(x,y) %in% Sep[[k]]))
          {
            RUE[x,y]=0
          }
        }
      }
    }
  }
  return(RUE)
}

