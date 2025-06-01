#' Generates parents from one layer back
#'
#' This routine generates the directed edges to each CC in layer i from
#' vertices in layer i-1, for each i.
#'
#' Routine is called from essgraph (generating essential graph)
#'
#' @param CC - list where CC[[k]] denotes vertices in the kth CC
#' @param ordlay where ordlay [[j]] gives the labels of the CCs in layer j
#' @param VIL where VIL[[k]] denotes the collection of all vertices in
#' layer k
#' @param CCued - list of matrices CCued[[k]] denotes undirected edges
#' for CC[[k]]
#' @param p - the 'extra edge' probability
#'
#' @return L - list of vectors,
#' where, for CC[[k]] in layer i, M[[k]] gives the parents
#' from layer i-1
#'
#' use *rejection method* for layer 1 parents of layer 2 CC (reject if set
#' of l1 vertices is complete)
#' @export



DEi2iplus1 <- function(CC,ordlay,VIL,CCued,p){


L <- vector(mode="list",length=length(CC))
L <-lapply(L,as.integer)

if(length(ordlay)==1){
    return(L)
    #this if deals with the situation where everything
    #is in a single layer - i.e. no directed edges
  }
  if(length(ordlay) >= 2)
    {
    x = 2 + rbinom((length(ordlay[[2]])),(length(VIL[[1]])-2),p)
    #we're dealing with CCs in layer 2.
    #x lists the number of layer 1 vertices to be chosen as parents
    #for each CC in layer 2.
    set = 1:length(VIL[[1]])
    for(i in 1:length(ordlay[[2]]))
      {
      repeat{
      h = sort(sample(set,x[i],replace = FALSE))
      z = VIL[[1]][h]
      if (checkcomplete(z,CC,ordlay,CCued) == 0)
        {
          L[[ordlay[[2]][i]]] = z
          #need to check that collection of nodes is not complete
          break
          }
      }
    }
    if(length(ordlay)==2)
      {
      return(L)
    }
    for(j in 3:length(ordlay))
      {
      a = 1 + rbinom(length(ordlay[[j]]),length(VIL[[j-1]])-1,p)
      #a lists the number of layer j-1 vertices to be chosen as parents
      #for each CC in layer j.
      for(i in 1:length(ordlay[[j]]))
        {
        set = 1:length(VIL[[j-1]])
        h = sort(sample(set,a[i],replace=FALSE))
        L[[ordlay[[j]][i]]] = VIL[[j-1]][h]
      }
      }
    return(L)
  }
}
