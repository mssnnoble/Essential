#' Checks if a set in the ground layer is complete
#'
#' checks if a set x in layer 1 is complete. This is because, when generating
#' an essential graph, we use the rejection method (for each vertex in layer 2
#' we choose a set from layer 1 at random and reject if the set is complete)
#'
#' @param x - set to be checked
#' @param CC - chain components
#' @param ordlay - layering of the CCs
#' @param CCued - collection of undirected edge matrices for the CCs
#'
#' @return a, which is 1 if complete and 0 otherwise
#'
#'@export


checkcomplete <- function(x,CC,ordlay,CCued){
  #check if x is a complete set in layer 1.
  #this is to ensure that the 1 to 2 directed edges
  #are permissible.
  #the size of the set is determined
  #each of that size which is not complete
  #equally likely
  #rejection method used.
  #first if it is complete then all vertices must
  #belong to the same CC
  #if they do, we call the function 'complete'
  a = 0
  f <- function(i){1*(all(x %in% CC[[ordlay[[1]][i]]]))}
  f1 = sapply(1:length(ordlay[[1]]),f)
  if(sum(f1)>=1)
    {
    g1 = which(f1==1)
      z = match(x,CC[[ordlay[[1]][g1[1]]]])
      t <- length(z)
      y <-CCued[[ordlay[[1]][g1[1]]]][z,z]
      if(sum(y) == 0.5*t*(t-1))
      {
        a = 1
      }
    }
  return(a)
}
