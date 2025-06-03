#'Computes the probability of an essential graph
#'
#'
#'@param GRAPH an essential graph object
#' @param theta - Chinese Restaurant concentration parameter for CCs
#' @param alpha - Chinese Restaurant concentration parameter for layers
#' @param p - extra edge parameter
#'
#' @return L
#' L$PCCs probability of assignment of vertices to CCs
#' L$PCCued vector, length is number of CCs, each entry is probability
#' of edge configuration for the CC
#' L$Pordlay probability of assignment of CCs to layers (conditioned on CCs)
#' L$PDEi2ip1 vector, each entry is probability of one-layer-back parent sets
#' L$PDEex vector, each entry is probability of two-or-more layer back
#' parent sets (conditioned on preceding probabilities)
#'
#' All the numbers, when multiplied together, give the probability
#' of the essential graph.
#' @export

GraphPR <- function(GRAPH,theta,alpha,p){
L <- list()
CCs <- as.list(GRAPH$CC)
s <- sapply(CCs,length)
PCCs <- prs(s,theta)
L$PCCs <- PCCs

ordlay <- GRAPH$ordlay
Pordlay <- prlay(GRAPH,alpha)
L$Pordlay <- Pordlay

PCCued <- prCC(GRAPH,p)
L$PCCued <- PCCued

PDEi2iP1 <- pri2ip1(GRAPH,p)
L$PDEi2iP1 <- PDEi2iP1
PDEex <- DEexPr(GRAPH,p)
L$PDEex <- PDEex
return(L)
}
