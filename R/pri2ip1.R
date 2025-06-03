#' Probability of one-layer-back parent configurations
#'
#'computes the probability of the i to i+1 directed edges
#'
#' @param CG an essential graph object
#' @param p extra edge parameter
#' @return L list of probabilities of one-layer-back parent
#' configurations for each CC
#' @export

pri2ip1 <- function(CG,p){
i2ip1 = CG$DEi2ip1
ordlay = as.list(CG$ordlay)
CCs = CG$CC
CCued = CG$CCued
VIL = CG$VIL

  L <-rep(0,length(CCs))

  #eventually, L returns the list of probabilities
  #of the arrows into each CC conditioned on layering
  #and undirected edges.
  #the ground layer CCs initiated as 1 will remain as 1
  #

  CP <- compsetsize(CG)
  f <- function(j){1}
  L[ordlay[[1]]] = sapply(ordlay[[1]],f)

    #CCs in layer 1 have no parents with probability 1

  n <- length(VIL[[1]])

  if(length(ordlay)>=2)
  {
  f <- function(j){length(i2ip1[[ordlay[[2]][j]]])}
  k = sapply(1:length(ordlay[[2]]),f)
  g <- function(j){((choose(n-2,k[j]-2))*(p^(k[j]-2))*(((1-p)^(n-k[j]))))/((choose(n,k[j]))-(CP[k[j]]))}
  L[ordlay[[2]]] = sapply(1:length(ordlay[[2]]),g)
}

  if(length(ordlay)>=3)
  {
  for(j in 3:length(ordlay)){
    n <- length(VIL[[j-1]])
    h <- function(a){length(i2ip1[[ordlay[[j]][a]]])}
    k <- sapply(1:length(ordlay[[j]]),h)
    fn <- function(j){(k[j]/n)*(p^(k[j]-1))*((1-p)^(n-k[j]))}
    L[ordlay[[j]]] <- sapply(1:length(ordlay[[j]]),fn)

  }
    }
  return(L)
}
