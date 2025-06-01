 #' Runs a Kruskal algorithm
 #'
 #' runs a Kruskal algorithm on a set of maximal cliques of a decomposable graph
#' to get a junction tree
#'
#'@param C collection of maximal cliques of a decomposable graph
 #'
#' @return J lower triangular matrix 1 if edge in junction tree 0 otherwise.
#' @export


kruskal <- function(C){
  #Kruskal for getting a junction tree
  #C collection of maximal cliques
  if(length(C)==1)
  {
    J <- matrix(0,1,1)
    return(J)
  }
  W <- matrix(0,length(C),length(C))
  for(i in 1:length(C)){
    for(j in 1:length(C)){
      W[j,i] <- length(intersect(C[[i]],C[[j]]))
    }
  }
  W[upper.tri(W,diag=TRUE)] <- 0
  #W is weighted max clique graph
  #(only lower triangular entries non-zero)
  J <- matrix(0,length(C),length(C))

   #J returns 1 for used edges and 0 otherwise
  #Kruskal is a greedy algorithm: choose largest unused
  #edge provided it does not form a cycle.

  #ConComp lists connected components. Each new edge
  #must join a connected component

  ConComp <-  vector(mode="list",length=length(C))
  ConComp <-lapply(ConComp,as.integer)
  f <- function(j){j}
  ConComp <- lapply(1:length(C),f)
  #(initially, each vertex is in a separate CC)

  repeat{
    a = which(W==max(W),arr.ind=TRUE)
    #a gives list of maximal elements

    f <- function(j){1*(all(a[1,] %in% ConComp[[j]]))}
    f1 <-sapply(1:length(ConComp),f)
    if(sum(f1)==1){W[a[1,1],a[1,2]] <- 0}
    else{
      J[a[1,1],a[1,2]] <- 1
      W[a[1,1],a[1,2]] <- 0
      if(sum(J)==length(C)-1)
      {
        return(J)
      }
      gone <- function(j){1*(a[1,1] %in% ConComp[[j]])}
      gtwo <- function(j)(1*(a[1,2] %in% ConComp[[j]]))
      g1 <- sapply(1:length(ConComp),gone)
      g2 <- sapply(1:length(ConComp),gtwo)
      h1 = which(g1 == 1)[1]
      h2 = which(g2 == 1)[1]
      ConComp[[h1]] <- sort(c(ConComp[[h1]],ConComp[[h2]]))
      ConComp <- ConComp[-h2]}
  }
}
