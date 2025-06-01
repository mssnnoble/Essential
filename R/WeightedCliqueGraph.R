#' Computes the weighted maximal clique graph
#'
#' This function takes a collection of maximal cliques
#' and returns a lower triangular matrix W
#' where W[i,j] for i > j is the number of vertices in
#' C[[i]] intersect C[[j]]
#'
#' @param C collection of maximal cliques of a decomposable graph
#'
#' @return W lower triangular graph W[i,j] = |C[[i]] intersect C[[j]]|
#' for i > j
#' @export

WeightedCliqueGraph <- function(C){
  W <- matrix(0,length(C),length(C))
  f <- function(i,j){length(intersect(C[[i]],C[[j]]))}
  for(i in 2:length(C))
  {
    g <- function(j){f(i,j)}
    W[i,1:(i-1)] = sapply(1:(i-1),g)
  }
  return(W)
}
