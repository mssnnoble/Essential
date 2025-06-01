#' Creates a random tree using Pruefer codes
#'
#'This algorithm creates a random tree, each tree equally likely for
#'each CC. This uses pruefer codes.
#'
#' @param CC this is a list of all the chain components
#'
#'
#' @return CCedge list of lower triangular matrices one for each CC
#' where 1 denotes edge, 0 denotes no edge
#' @export

CCedge <- function(CC){
  if(is.list(CC)==FALSE)
  {
   a = list()
   a[[1]] = CC
   CC = list()
   CC = a
  }

  #CC is the list of CCs
  #if CC has 3 or more nodes, use matrix from pruefer2edge

  pruf <-  createpruf(CC)
  #list of Pruefer codes

  a <- vector(mode="list",length=length(CC))
  g <- function(i){matrix(0,length(CC[[i]]),length(CC[[i]]))}
  a <- lapply(1:length(CC),g)

  f <- function(i){pruefer2edge(pruf[[i]],CC[[i]])}
  a <- lapply(1:length(CC),f)
  CCedge <- a
  return(CCedge)
}
