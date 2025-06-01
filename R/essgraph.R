#' Generates random essential graph
#'
#' generates a random essential graph - the undirected edges of each CC
#' form a spanning tree.
#
#' @param theta - Chinese Restaurant concentration parameter for CCs
#' @param alpha - Chinese Restaurant concentration parameter for layers
#' @param p - extra edge parameter
#' @param d - number of vertices
#'
#' @return M
#' M$CC assignment of vertices to chain components
#' M$CCued lower triangular matrices, one for each CC
#'  where 1 denotes undirected edge
#' M$ordlay list where M$ordlay[[j]] gives the CC labels for layer j
#' M$VIL where M$VIL[[j]] gives list of vertices in layer j
#' M$prVIL where M$prVIL[[j]] gives list of vertices from layers 1,...,j
#' M$DEi2ip1 where  M$DEi2ip1[[j]] gives parents in layer i-1 if CC[[j]]
#' is in layer i.
#' M$DEex where M$DEex[[j]] gives parents in layers 1,...,i-2 if CC[[j]]
#' is in layer i
#' M$TOTPAR where M$TOTPAR[[j]] gives all parents of CC[[j]]
#'
#'@export


essgraph <- function(theta,alpha,p,d){

  M <- list()

  tables <- CRP(theta,d)
  CC <- CCsets(tables)

  l <- CRP(alpha,length(CC))
  lay <- CCsets(l)
  lad <- CCsetsADJ(CC,lay)
  ordlay <- orderlayers(CC,lad)

   CCued <- CCedge(CC)

  VIL <- vertinlayer(ordlay, CC)

  DEi2ip1 <- DEi2iplus1(CC,ordlay,VIL,CCued,p)

  prVIL <- progVIL(VIL)

  DEex <- DEextra(DEi2ip1,CC,CCued,prVIL,ordlay,p)

  TOTPAR <- vector(mode="list",length=length(CC))
  TOTPAR <- lapply(TOTPAR,as.integer)
  f <- function(j){sort(c(DEi2ip1[[j]],DEex[[j]]))}
  TOTPAR <- lapply(1:length(CC),f)

  M$CC <- CC
  M$CCued <- CCued
  M$ordlay <- ordlay
  M$VIL <- VIL
  M$prVIL <- prVIL
  M$DEi2ip1 <- DEi2ip1
  M$DEex <-DEex
  M$TOTPAR <- TOTPAR


  return(M)
}
