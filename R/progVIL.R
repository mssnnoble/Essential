#' computes union of layer together with all previous layers
#'
#' @param VIL - list where VIL[[j]] gives vertices in each layer
#'
#'
#' @return g3 where g3[[j]] is union of VIL[[i]] for i=1,...,j
#'
#'@export

progVIL <- function(VIL){
  #This gives the union of the vertices of layers 1,...,i
  #for i = 1,...,length(ordlay)

  g3 <- vector(mode="list",length=length(VIL))
  g3 <-lapply(g3,as.integer)

  g3[[1]] = VIL[[1]]
  if(length(VIL)==1)
  {return(g3)}
  else{
    for(j in 2:length(VIL)){
      g3[[j]] <- sort(union(g3[[j-1]],VIL[[j]]))}
    return(g3)
  }
}
