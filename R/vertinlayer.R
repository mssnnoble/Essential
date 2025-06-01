#' Union of all CCs in each layer
#'
#'This simply takes the union of all CCs in a layer, for each layer.
 #' @param CC - collection of chain components
 #' @param ordlay - which CCs are in each layer
#'
#' @return
#' g2 = list g2[[i]] all vertices in layer i
#' @export

vertinlayer <- function(ordlay,CC){

  g2 <- vector(mode="list",length=length(ordlay))
  g2 <-lapply(g2,as.integer)

  for(j in 1:length(ordlay))
    {
    L = vector(mode="list",length = length(ordlay[[j]]))
    L = lapply(L,as.integer)
    f <- function(k){CC[[ordlay[[j]][k]]]}
    L = lapply(1:length(ordlay[[j]]),f)
    a <- sort(unlist(L))
    g2[[j]] <- a
  }
  return(g2)
}
