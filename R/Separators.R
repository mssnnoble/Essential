#' Computes list of separators
#'
#' @param M is the lower triangular adjacency matrix for a connected
#' decomposable graph
#'
#' @return S the list of separators (according to multiplicity)
#' @export


Separators <- function(M)
{
  C <- MC(M)
  JTr = kruskal(C)

  #JTr gives a lower triangular matrix which defines a junction tree

  if(nrow(JTr)==1)
  {
    S = vector(mode="integer",length=0L)
    return(S)
  }
  S <-  vector(mode="list",length=nrow(JTr)-1)
  S <-lapply(S,as.integer)
  #S will list the separators according to multiplicity

  count = 1
  for(x in 2:length(C))
  {
    for(y in 1:(x-1))
    {
      if(JTr[x,y]==1)
      {
        S[[count]] = intersect(C[[x]],C[[y]])
        count <- count+1
      }
    }
  }

  #S now lists all the  separators according to multiplicity
  return(S)
}
